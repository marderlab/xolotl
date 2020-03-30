// _  _ ____ _    ____ ___ _
//  \/  |  | |    |  |  |  |
// _/\_ |__| |___ |__|  |  |___
//
// Cell Int1 Potassium Conductance
// http://www.jneurosci.org/content/18/13/5053
// Nadim 1998
#ifndef KD
#define KD
#include "conductance.hpp"

//inherit conductance class spec
class Kd: public conductance {

public:


    //specify both gbar and erev and initial conditions
    Kd(double gbar_, double E_, double m_)
    {
        gbar = gbar_;
        E = E_;
        m = m_;

         // defaults
        if (isnan(gbar)) { gbar = 0; }
        
        if (isnan (E)) { E = 0; }

        p = 4;

        // allow this channel to be approximated
        approx_m = 1;

    }

    double m_inf(double, double);
    double tau_m(double, double);
    string getClass(void);

};

string Kd::getClass(){return "Kd";}

double Kd::m_inf(double V, double Ca) {return (1.0/(1.0+exp(-0.045*((V)+25.0))));}
double Kd::tau_m(double V, double Ca) {return (154.0/(1.0+exp(0.065*((V)+30.0))));}


#endif
