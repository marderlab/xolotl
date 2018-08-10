// _  _ ____ _    ____ ___ _
//  \/  |  | |    |  |  |  |
// _/\_ |__| |___ |__|  |  |___
//
// Abstract class for defining compartments
#ifndef CONTAINER
#define CONTAINER
#include <cmath>
#include <vector>
// #include "conductance.hpp"
// #include "synapse.hpp"
// #include "mechanism.hpp"

using namespace std;

class container {
protected:


public:

  // general parameters
  double verbosity = 0;
  double RT_by_nF;
  double V;
  double Ca;
  double i_Ca;
  double I_ext;

  // spatial variables
  container* upstream;
  container* downstream;

  // this int stores an integer that indicates
  // the hierarchy of this container in a multi-comp
  // neuron tree. tree_idx of 0 means it is a soma
  // container
  double tree_idx;

  // this number stores a numeric value
  // that corresponds to the neuron # of this compartment
  // in the whole network. automatically assigned
  double neuron_idx;

  // housekeeping variables
  int n_axial_syn;
  double V_prev;
  double Ca_prev;
  double i_Ca_prev;


  // constructor
  container()
  {
    container = 0; // null pointer for safety
  }

  ~container() {}

  // DECLARE FUNCTIONS HERE

  // integration methods
  void integrateMechanisms(double);
  void integrateChannels(double, double);
  void integrateSynapses(double, double);
  void integrateVoltage(double, double);
  void integrateV_clamp(double, double, double);

  // methods for integrating using Crank-Nicholson
  // and methods for multi-compartment models
  double getBCDF(int);
  void integrateCNFirstPass(double);
  void integrateCNSecondPass(double);
  void resolveAxialConductances(void);

};

// DETAIL CONTAINER:: FUNCTIONS HERE



#endif
