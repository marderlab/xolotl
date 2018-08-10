// _  _ ____ _    ____ ___ _
//  \/  |  | |    |  |  |  |
// _/\_ |__| |___ |__|  |  |___
//
// class that defines a network
// a network can either be a network of
// single AbstractCompartment neurons, or a
// multi-compartment neuron

#ifndef NETWORK
#define NETWORK
#include <cmath>
#include <vector>
#include "AbstractCompartment.hpp"
#include "mex.h"

using namespace std;


class network
{
protected:
    // housekeeping
    AbstractCompartment* temp_comp;
    AbstractCompartment* last_valid_comp;

    // store pointers to all soma
    vector<AbstractCompartment*> soma_comp;

    // keeps track of the # of terminals in terminal_comp
    int n_terminals;
public:
     // pointers to all compartments in network
    vector<AbstractCompartment*> comp;

    // temperature
    double temperature;
    double temperature_ref;

    // housekeeping
    int n_comp = 0;
    int n_soma = 0; // will be used in the Crank-Nicholson scheme

    double verbosity;

    // constructor
    network() {}

    // function declarations
    void integrate(double,double *, double);
    void integrateClamp(double, double *, double);
    void addCompartment(AbstractCompartment*);
    void resolveTree(void);

};

void network::resolveTree(void)
{
    AbstractCompartment * connected_comp = NULL;

    if (verbosity > 0)
    {
        mexPrintf("[C++] network::resolveTree() called\n");
    }


    // ttl =  this_tree_level
    for (int ttl = 0; ttl < n_comp; ttl++)
    {
        // find all AbstractCompartments with this_tree_level
        for (int i = 0; i < n_comp; i++)
        {
            if (isnan(comp[i]->tree_idx)) {continue;}
            if ((comp[i]->tree_idx) != ttl) {continue;}

            // OK, this AbstractCompartment has the tree level we
            // are currently interested in

            if (comp[i]->tree_idx == 0)
            {
                // mexPrintf("found a soma, calling it = %i\n",n_soma);
                comp[i]->neuron_idx = n_soma;
                n_soma++;
                soma_comp.push_back(comp[i]);

            }

            // now go over every synapse that impinges on
            // this AbstractCompartment and check if they're Axial
            // and if so get pointers to those presyn AbstractCompartments

            for (int j = 0; j < comp[i]->n_axial_syn; j ++ )
            {
                connected_comp = comp[i]->getConnectedCompartment(j);
                if (!connected_comp){continue;}
                if (isnan(connected_comp->tree_idx))
                {
                    double child_tree_idx = ttl+1;

                    // set it
                    (connected_comp->tree_idx) = child_tree_idx;

                    // wire up stream pointers
                    (comp[i]->downstream) = connected_comp;
                    (connected_comp->upstream) = comp[i];

                    connected_comp->neuron_idx = comp[i]->neuron_idx;
                }
                else if ((connected_comp->tree_idx) == (ttl+1)) {
                    // connected_comp already has a tree_idx
                    // possibly manually entered, or from a previous
                    // integrate. if compatible, wire up stream
                    // pointers
                    (comp[i]->downstream) = connected_comp;
                    (connected_comp->upstream) = comp[i];

                    connected_comp->neuron_idx = comp[i]->neuron_idx;

                }
            }
        }
    }



    // OK, now we have resolved the tree.
    // now, we need to mark the downstream_g and
    // upstream_g for every AbstractCompartment

    for (int i = 0; i < n_comp; i ++)
    {
        comp[i]->resolveAxialConductances();
    }



    // go over every AbstractCompartment, and check that stream
    // pointers and gs match up

    if (verbosity > 0)
    {

        for (int i = 0; i < n_comp; i++)
        {
            mexPrintf("---------------\n");
            mexPrintf("this comp tree_idx = %f\n",comp[i]->tree_idx);
            if (comp[i]->downstream)
            {
                mexPrintf("downstream pointer exists\n");

            } else {
                mexPrintf("NO downstream pointer\n");
            }
            mexPrintf("downstream_g =  %f\n", comp[i]->downstream_g);
            if (comp[i]->upstream)
            {
                mexPrintf("upstream pointer exists\n");

            } else {
                mexPrintf("No upstream pointer\n");
            }
            mexPrintf("upstream_g =  %f\n", comp[i]->upstream_g);

        }
    }

}


// add a AbstractCompartment to the network -- network contains
// a vector of pointers to AbstractCompartments
void network::addCompartment(AbstractCompartment *comp_)
{
    comp.push_back(comp_);
    n_comp++;
    comp_->verbosity = verbosity;
    comp_->RT_by_nF = (0.0431)*(temperature + 273.15);

    if (verbosity > 0)
    {
        mexPrintf("[C++] adding AbstractCompartment to network. \n");
    }
}

// this integrate method works for networks
// of single AbstractCompartments, or cells with
// multiple AbstractCompartments under normal
// conditions. Don't use if something is
// being voltage clamped!
void network::integrate(double dt, double * I_ext_now, double delta_temperature)
{

    // we will use Exponential Euler for single-AbstractCompartment
    // models and networks, and Crank-Nicholson for
    // multi-AbstractCompartment models. How do we know which is which?
    // all AbstractCompartments with a neuron_idx will be assumed to
    // be part of a multi-AbstractCompartment model
    // all AbstractCompartments where neuron_idx is NaN will be treated
    // as single AbstractCompartments and integrated normally

    // integrate all channels in all AbstractCompartments
    for (int i = 0; i < n_comp; i++)
    {

        // move current values to previous values
        comp[i]->V_prev = comp[i]->V;
        comp[i]->Ca_prev = comp[i]->Ca;
        comp[i]->i_Ca_prev = comp[i]->i_Ca;

        comp[i]->i_Ca = 0;
        comp[i]->I_ext = I_ext_now[i];

        // integrate controllers
        comp[i]->integrateMechanisms(dt);

        comp[i]->integrateChannels(dt, delta_temperature);



        // integrate synapses
        if (isnan(comp[i]->neuron_idx))
        {
            comp[i]->integrateSynapses(dt, delta_temperature);
        }

    }

    // integrate voltages in all single AbstractCompartments
    // and calcium in all AbstractCompartments
    for (int i = 0; i < n_comp; i++)
    {
        if (isnan(comp[i]->neuron_idx)) {
            comp[i]->integrateVoltage(dt, delta_temperature);
        } else {
            // this is a multi-AbstractCompartment model,
            // so just integrate the calcium
            // this is assumed to be already done in the
            // mechanism integration, so do nothing here
        }

    }




    // OK, now we have integrated all single AbstractCompartments,
    // but in multi AbstractCompartment models, the
    // voltages are un-integrated. since they are coupled to each
    // other, we use the Crank Nicholson scheme to integrate them
    // I will follow the scheme described in Dayan and Abbott
    // (Chapter 6 Appendix B)
    // Convention will be the same, primes in the text will
    // be replaced with _
    // integration starts at the soma (mu = 0), and proceeds
    // down the cable

    // for each cable
    for (int nidx = 0; nidx < n_soma; nidx++)
    {


        temp_comp = soma_comp[nidx];


        // go down the cable -- away from soma
        // to terminal, increasing in tree_idx
        while (temp_comp)
        {
            temp_comp->integrateCNFirstPass(dt);
            last_valid_comp = temp_comp;
            temp_comp = temp_comp->downstream;
            // nothing after this--temp_comp may be NULL
        }


        temp_comp = last_valid_comp;


        // go up the cable -- towards soma
        // away from terminal, decreasing in tree_idx
        while (temp_comp)
        {
            temp_comp->integrateCNSecondPass(dt);
            temp_comp = temp_comp->upstream;
            // nothing after this, becaue temp_comp may be NULL
        }
    }
}

// integrate while voltage clamping some AbstractCompartments
void network::integrateClamp(double dt, double *V_clamp, double delta_temperature)
{

    // integrate all channels in all AbstractCompartments
    for (int i = 0; i < n_comp; i++)
    {

        // move current values to previous values
        comp[i]->V_prev = comp[i]->V;
        comp[i]->Ca_prev = comp[i]->Ca;
        comp[i]->i_Ca_prev = comp[i]->i_Ca;


        comp[i]->i_Ca = 0;
        comp[i]->I_ext = 0;
        comp[i]->I_clamp = 0;

        comp[i]->integrateChannels(dt, delta_temperature);

        // integrate synapses
        comp[i]->integrateSynapses(dt, delta_temperature);
    }

    // integrate all voltages and Ca in all AbstractCompartments
    for (int i = 0; i < n_comp; i++)
    {

        if (isnan(V_clamp[i])){

            comp[i]->integrateVoltage(dt, delta_temperature);
        } else {
            comp[i]->integrateV_clamp(V_clamp[i], dt, delta_temperature);
        }

    }

}


#endif
