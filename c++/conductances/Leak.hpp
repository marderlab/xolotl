// _  _ ____ _    ____ ___ _
//  \/  |  | |    |  |  |  |
// _/\_ |__| |___ |__|  |  |___
//
// Slow Calcium conductance
#ifndef LEAK
#define LEAK
#include "conductance.hpp"

//inherit conductance class spec
class Leak: public conductance {

public:

    // specify parameters + initial conditions
    Leak(double gbar_, double E_)
    {
        gbar = gbar_;
        g = gbar; // this is important as integrate doesn't do anything in the leak channels
        E = E_;

        if (isnan (E)) { E = -55; }

    }


    void integrate(double, double);
    void integrateMS(int, double, double);
    void integrateLangevin(double, double);

    string getClass(void);

};

string Leak::getClass(){return "Leak";}

void Leak::integrate(double V, double Ca) {
    // do nothing
    // because there is nothing to integrate
    g = gbar;
}

// Runge-Kutta 4 integrator
void Leak::integrateMS(int k, double V, double Ca) {

    // do nothing
    // because there is nothing to integrate
    g = gbar;
}

void Leak::integrateLangevin(double V, double Ca) {
    // do nothing 
    g = gbar;
}



#endif
