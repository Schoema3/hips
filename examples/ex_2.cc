/**
 * @file main.cpp
 * @brief Main program for HiPS simulation with chemical reactions.
 */

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
 * @brief Main function for HiPS simulation with chemical reactions.
 * @details This program sets up and runs a HiPS simulation with chemical reactions using Cantera.
 * 
 * @return 0 on successful execution.
 */
int main() {
  
    //---------- setting  HiPS tree and constructor

    int            nLevels      = 6;           // # HiPS tree levels
    double         domainLength = 0.01;        // HiPS tree lengthscale
    double         tau0         = 0.0005;    // HiPS tree timescale
    double         C_param      = 0.5;         // HiPS eddy rate multiplier
    double         tRun         = 0.05;      // tree evolution time
    int            forceTurb    = 0;           // Force turbulent profile
    vector<double> ScHips(54,1);               // Schmidt number is unity for all species

    //---------- setting gas solution

    auto cantSol = Cantera::newSolution("gri30.yaml");
    auto gas = cantSol->thermo();
    size_t nsp = gas->nSpecies();
    int nVar = nsp + 1; // 1 (enthalpy) + number of species 

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

    vector<vector<double>> ysp(53, vector<double>(nparcels, 0)); // Vector to store species mass fractions for each parcel
    vector<double> h(nparcels);                                   // Vector to store enthalpy for each parcel
    vector<double> T(nparcels);                                   // Vector to store temperature for each parcel

    vector<double> y0(nsp);                                       // Initial state vector for fresh reactants
    double T0;
    double h0;

    vector<double> y1(nsp);                                       // Initial state vector for pre-combusted parcels
    double T1;
    double h1;
    vector<string> variableNames(nsp + 1);                        // Names of variables (enthalpy and species)

    //-------------- set state vectors in each parcel
    double fracBurn = 0.25; // 25% of parcels are pre-combusted

    T0 = 300.0;
    gas->setState_TPX(T0, Cantera::OneAtm, "C2H4:1, O2:3, N2:11.25");
    h0 = gas->enthalpy_mass();

    gas->getMassFractions(y0.data());

    for (int i = 0; i < nsp; i++) {
        for (int j = 0; j < 25; j++) { // Set state for fresh reactants
            h[j] = h0;
            T[j] = gas->temperature();
            ysp[i][j] = y0[i];
            //cout << " i " << i << " j " << j << " ysp " << gas->temperature() << endl;
        }
    }
 
    T1 = 300.0;
    gas->setState_TPX(T0, Cantera::OneAtm, "C2H4:1, O2:3, N2:11.25");
    h1 = gas->enthalpy_mass();

    gas->setState_HP(h1, gas->pressure());
    gas->equilibrate("HP");
    gas->getMassFractions(y1.data());

    for (int i = 0; i < nsp; i++) {
        for (int j = 25; j < nparcels; j++) { // Set state for pre-combusted parcels
            h[j] = h1;
            T[j] = gas->temperature();
            ysp[i][j] = y1[i];
            //cout << " i " << i << " j " << j << "  " << gas->speciesName(i) << " ysp " << gas->temperature() << endl;
        }
    }

    HiPS.Temp = T;

    variableNames[0] = "enthalpy";
    for(int i = 0; i < ysp.size(); i++)
        variableNames[i + 1] = gas->speciesName(i);

    vector<double> weight(nparcels, 1.0 / nparcels); // Set weight for each parcel

    HiPS.set_varData(h, weight, variableNames[0]);
    for (int k = 0; k < ysp.size(); k++) 
        HiPS.set_varData(ysp[k], weight, variableNames[k + 1]);

    //  //----------------- advancing HiPS to do rxn and mixing
    HiPS.calculateSolution(tRun, true);

    return 0;
}

