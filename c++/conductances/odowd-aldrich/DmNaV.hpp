// _  _ ____ _    ____ ___ _    
//  \/  |  | |    |  |  |  |    
// _/\_ |__| |___ |__|  |  |___ 
//
// Drosophila NaV
// Voltage-clamp analysis of sodium channels in wild-type and mutant Drosophila neurons (1988)
// https://www.ncbi.nlm.nih.gov/pubmed/2848103
// assumes a reversal potential of +47mV

#ifndef DMNAV
#define DMNAV
#include "../../conductance.hpp"

//inherit conductance class spec
class DmNaV: public conductance {

public:

    // specify parameters + initial conditions 
    DmNaV(double g_, double E_, double m_, double h_)
    {
        gbar = g_;
        E = E_;
        m = m_;
        h = h_;
    }
    
    void integrate(double V, double Ca, double dt);
    void connect(compartment *pcomp_);
    double m_inf(double V);
    double h_inf(double V);
    double tau_m(double V);
    double tau_h(double V); 
};

void DmNaV::connect(compartment *pcomp_) {container = pcomp_; }

void DmNaV::integrate(double V, double Ca, double dt)
{
    m = m_inf(V) + (m - m_inf(V))*exp(-dt/tau_m(V));
    h = h_inf(V) + (h - h_inf(V))*exp(-dt/tau_h(V));
    g = gbar*m*m*m*m*m*h;
}

double DmNaV::m_inf(double V) {return 1.0/(1.0+exp((V+34.3)/-8.79));}
double DmNaV::h_inf(double V) {return 1.0/(1.0+exp((V+48)/6));}
double DmNaV::tau_m(double V) {return 3 - 2.4/(1+exp((V+33.6)/-9.0));}
double DmNaV::tau_h(double V) {return 3 - 2.53/(1+exp((V+22.8)/-3.5));}

#endif