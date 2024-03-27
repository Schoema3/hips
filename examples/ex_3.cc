#include "hips.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo.h"

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;


///////////////////////////////////////////////////////////////////////////////

/**
 * @file main.cpp
 * @brief Main program for HiPS simulation with chemical reactions.
 */

/**
 * @brief Main function for HiPS simulation with chemical reactions.
 * @details This program sets up and runs a HiPS simulation with chemical reactions using Cantera.
 * 
 * @return 0 on successful execution.
 */
int main() {
  
    //---------- setting  HiPS tree and constructor

    int            nLevels      = 9;           // # HiPS tree levels
    double         domainLength = 1.0;         // HiPS tree lengthscale
    double         tau0         = 1.0;         // HiPS tree timescale
    double         C_param      = 0.5;         // HiPS eddy rate multiplier
    double         tRun         = 4.0;       // tree evolution time
    int            forceTurb    = 0;           // Force turbulent profile
    vector<double> ScHips(54,1);

    //---------- setting gas solution

    //auto   cantSol = Cantera::newSolution("2S_CH4_BFER.yaml", "2S_CH4_BFER", "None");
    auto   cantSol = Cantera::newSolution("gri30.yaml");
    auto   gas = cantSol->thermo();
    size_t nsp = gas->nSpecies();
    int    nVar = nsp+1;           // 1 + nsp: h, y 

    hips HiPS(nLevels,
          domainLength,
          tau0, 
          C_param,
          forceTurb,
          nVar,
          ScHips, 
          true
#ifdef REACTIONS_ENABLED
          , cantSol
#endif
        );

    int nparcels = HiPS.nparcels;

    //------------- declaring state vectors

    vector<vector <double > > ysp(nsp, vector<double> (nparcels, 0));
    vector<double>            y0(nsp);
    vector<double>            y1(nsp);
    double                    h0;
    
    vector<double>            T(nparcels); 
    double                    T0;
    double                    T1;
    double                    h1;

    vector<string>            variableNames(nsp+1);

    //-------------- set state vectors in each parcel

    T0 = 700.0;
    gas->setState_TPX(T0, Cantera::OneAtm, "CH4:1");
    gas->getMassFractions(&y0[0]);
    h0 = gas->enthalpy_mass();


    T1 = 1200;
    gas->setState_TPX(T1, Cantera::OneAtm, "O2:2, N2:7.52");
    gas->getMassFractions(&y1[0]);
    h1 = gas->enthalpy_mass();

    for (int i=0; i<nparcels/2; i++) {       // left parcels are stream 0 (air)
        T[i] = T0;
        for (int k=0; k<nsp; k++)
            ysp[k][i] = y0[k];
    }

    for (int i=nparcels/2; i<nparcels; i++) { // right parcels are stream 1 (fuel)
        T[i] = T1;
        for (int k=0; k<nsp; k++)
            ysp[k][i] = y1[k];
    }

    //----------------- set array of hips variables from state variables
     
   variableNames[0] = "Temperature";
   for(int i=0; i<ysp.size(); i++)
       variableNames[i+1] = gas->speciesName(i);

   vector<double>           weight(nparcels, 1.0/nparcels);

    //----------------- set array of hips variables from state variables

    HiPS.set_varData(T,weight, variableNames[0],  0);
    for (int k=0; k<ysp.size(); k++) 
        HiPS.set_varData(ysp[k],weight, variableNames[k+1]);

    //----------------- advancing HiPS to do rxn and mixing 
    HiPS.calculateSolution(tRun);

    return 0;

}  
