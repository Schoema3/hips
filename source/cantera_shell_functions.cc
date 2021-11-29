/**
 * @file cantera_shell_functions.cc
 * Source file for class cantera_shell_functions
 */

#ifndef DOCANTERA

// This file is used to store dummy function definitions when cantera is not available
// Just for compilation

#include "cantera_shell_functions.h"
#include <cstdlib>
#include <iostream>

using namespace std;
using namespace Cantera;

//-------------------------------------------------------------------------------

void errorMessage(string name) {
    cout << "\nERROR, FUNCTIONS NOT IMPLEMENTED IN CANTERA_SHELL_FUNCTIONS " 
         << name << endl;
    exit(0);
}
void warningMessage(string name) {
    cout << "\nWarning (" << name <<"), code was built without cantera."
         << "\n  Using shell functions to facilitate compilation," 
         << "\n  but these functions are not implemented so don't run cases that need them. " 
         << endl;
}

//-------------------------------------------------------------------------------

void           IdealGasPhase::setState_PY(double p, double* y)             { errorMessage("setState_PY"); }
void           IdealGasPhase::setState_HP(double h, double p, double tol)  { errorMessage("setState_HP"); }
void           IdealGasPhase::setState_TPY( double T, double p, double *y) { errorMessage("setState_TPY"); }
double         IdealGasPhase::temperature()                                { errorMessage("temperature"); return 0; }
double         IdealGasPhase::density()                                    { errorMessage("density"); return 0;}
double         IdealGasPhase::cp_mass()                                    { errorMessage("cp_mass"); return 0;}
double         IdealGasPhase::enthalpy_mass()                              { errorMessage("enthalpy_mass"); return 0;}
void           IdealGasPhase::getEnthalpy_RT(double* hrt)                  { errorMessage("getEnthalpy_RT"); }
void           IdealGasPhase::getNetProductionRates(double* rr)            { errorMessage("getNetProductionRates"); }
double         IdealGasPhase::atomicWeight(int k)                          { errorMessage("atomicWeight"); return 0; }
double         IdealGasPhase::molecularWeight(int k)                       { errorMessage("molecularWeight"); return 0; }
vector<double> IdealGasPhase::molecularWeights()                           { errorMessage("molecularWeights"); return vector<double>(); }
double         IdealGasPhase::meanMolecularWeight()                        { errorMessage("meanMolecularWeight"); return 0; }
int            IdealGasPhase::speciesIndex(string name)                    { errorMessage("speciesIndex"); return 0; }
int            IdealGasPhase::elementIndex(string name)                    { errorMessage("elementIndex"); return 0; }
int            IdealGasPhase::nSpecies()                                   { errorMessage("nSpecies"); return 0; }
int            IdealGasPhase::nElements()                                  { errorMessage("nElements"); return 0; }
void           IdealGasPhase::setMoleFractions(const double *x)            { errorMessage("setMoleFractions"); }
void           IdealGasPhase::setMassFractions(const double *y)            { errorMessage("setMassFractions"); }
void           IdealGasPhase::getMassFractions(double *y)                  { errorMessage("getMassFractions"); }
void           IdealGasPhase::getMoleFractions(double *x)                  { errorMessage("getMoleFractions"); }
void           IdealGasPhase::getAtoms(int k, double *Z)                   { errorMessage("getAtoms"); }
double         IdealGasPhase::moleFraction(int k)                          { errorMessage("moleFraction"); return 0; }
string         IdealGasPhase::speciesName(int k)                           { errorMessage("speciesName"); return ""; }
int            IdealGasPhase::nAtoms(int j, int k)                         { errorMessage("nAtoms"); return 0; }
void           IdealGasPhase::equilibrate(const string& XY,
                   const string& solver, double rtol, int max_steps, 
                   int max_iter, int estimate_equil, int log_level)        { errorMessage("equilibrate"); return; }   // inherited memeber from ThermoPhase



IdealGasPhase::IdealGasPhase(string mechfile)                              { warningMessage("IdealGasPhase"); }


//-------------------------------------------------------------------------------

double Transport::viscosity()                   { errorMessage("viscosity"); return 0; }
double Transport::thermalConductivity()         { errorMessage("thermalCond"); return 0; }
void   Transport::getMixDiffCoeffs(double *D)   { errorMessage("getMixDiffC"); }

//-------------------------------------------------------------------------------

Transport* Cantera::newTransportMgr(string s1, IdealGasPhase *gas) { warningMessage("newTransportMgr"); return 0;}

#endif

