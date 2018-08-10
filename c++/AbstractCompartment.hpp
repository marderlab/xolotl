// _  _ ____ _    ____ ___ _
//  \/  |  | |    |  |  |  |
// _/\_ |__| |___ |__|  |  |___
//
// Abstract class for defining compartments
#ifndef ABSTRACTCOMPARTMENT
#define ABSTRACTCOMPARTMENT
#include <cmath>
#include <vector>
// #include "conductance.hpp"
// #include "synapse.hpp"
// #include "mechanism.hpp"

using namespace std;

class AbstractCompartment {
protected:


public:

  // general parameters
  double verbosity = 0;
  double RT_by_nF;
  double V;
  double Ca;
  double i_Ca;
  double I_ext;
  double I_clamp; // this is the current required to clamp it
  double E_Ca;
  double Ca_average;

  // spatial variables
  AbstractCompartment* upstream;
  AbstractCompartment* downstream;

  // this int stores an integer that indicates
  // the hierarchy of this AbstractCompartment in a multi-comp
  // neuron tree. tree_idx of 0 means it is a soma
  // AbstractCompartment
  double tree_idx;

  // this number stores a numeric value
  // that corresponds to the neuron # of this compartment
  // in the whole network. automatically assigned
  double neuron_idx;

  // housekeeping variables
  int n_cond; // this keep tracks of the # channels
  int n_cont; // # of mechanisms
  int n_syn; // # of synapses
  int n_axial_syn;
  double V_prev;
  double Ca_prev;
  double i_Ca_prev;

  // conductances downstream and upstream
  // will be generated on initialization
  double downstream_g;
  double upstream_g;


  // constructor
  AbstractCompartment() {}

  ~AbstractCompartment() {}

  // DECLARE FUNCTIONS HERE
  virtual AbstractCompartment* getConnectedCompartment(int) = 0;

  // methods to retrieve information from compartment
  virtual int getFullMechanismState(double*, int) = 0;
  virtual int getFullCurrentState(double*, int) = 0;
  virtual int getFullSynapseState(double*, int) = 0;
  virtual int getFullMechanismSize(void) = 0;
  virtual int getFullSynapseSize(void) = 0;

  // integration methods
  virtual void integrateMechanisms(double) = 0;
  virtual void integrateChannels(double, double) = 0;
  virtual void integrateSynapses(double, double) = 0;
  virtual void integrateVoltage(double, double) = 0;
  virtual void integrateV_clamp(double, double, double) = 0;

  // methods for integrating using Crank-Nicholson
  // and methods for multi-compartment models
  virtual double getBCDF(int) = 0;
  virtual void integrateCNFirstPass(double) = 0;
  virtual void integrateCNSecondPass(double) = 0;
  virtual void resolveAxialConductances(void) = 0;

};

// DETAIL CONTAINER:: FUNCTIONS HERE



#endif
