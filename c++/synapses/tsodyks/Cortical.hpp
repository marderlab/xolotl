// Cholingeric Synapse
#ifndef CORTICAL
#define CORTICAL
#include "synapse.hpp"

class Cortical: public synapse {

public:

    // specify parameters + initial conditions
    Cortical(double g_, double s_)
    {
        gbar = g_;
        E = -80.0;
        Delta = 5.0;
        Vth = -35.0;
        k_ = 0.01;
        s = s_;

        // defaults
        if (isnan (s)) { s = 0; }
        if (isnan (gbar)) { gbar = 0; }
        is_electrical = false;
    }

    void integrate(double dt);
    int getFullStateSize(void);
    void connect(compartment *pcomp1_, compartment *pcomp2_);
    double getCurrent(double V_post);
    int getFullState(double*, int);
};

int Cortical::getFullStateSize()
{
    return 2;
}

void Cortical::integrate(double dt)
{
    // figure out the voltage of the pre-synaptic neuron
    double V_pre = pre_syn->V;

    // find s_inf
    double s_inf = 1.0/(1.0+exp((Vth - V_pre)/Delta));

    // integrate using exponential Euler
    double tau_s = (1 - s_inf)/k_;

    s = s_inf + (s - s_inf)*exp(-dt/tau_s);

}


int Cortical::getFullState(double *syn_state, int idx)
{
    // give it the current synapse variable
    syn_state[idx] = s;

    idx++;

    // also return the current from this synapse
    syn_state[idx] = gbar*s*(post_syn->V - E);
    idx++;
    return idx;
}

void Cortical::connect(compartment *pcomp1_, compartment *pcomp2_)
{
    pre_syn = pcomp1_;
    post_syn = pcomp2_;

    // tell the post-synaptic cell that we're connecting to it
    post_syn->addSynapse(this);
}




#endif


tinact = 3 msec,
trec = 800 msec, USE = 0.67, ASE = 250 pA, tmem = 50 msec
