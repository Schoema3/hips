
/// \file ex_3.cc
/// \brief Example demonstrating non-premixed combustion using HiPS and Cantera.
/// 
/// This example simulates non-premixed turbulent flames using a HiPS domain with six levels. 
/// It initializes separate parcels for fuel and oxidizer with a temperature gradient. 
/// The simulation uses Cantera to handle chemical kinetics and calculates the resulting 
/// temperature and species mass fractions.
/// 
/// ## Compilation:
/// ```bash
/// cd build
/// cmake .. -DREACTIONS_ENABLED=ON
/// make
/// sudo make install
/// ```
/// 
/// ## Execution:
/// ```bash
/// cd ../run
/// ./ex_3.x
/// ```
/// 
/// ## Output:
/// The simulation generates a file containing temperature and species concentration profiles.
/// 
/// ## Dependencies:
/// - HiPS library (compiled with reaction support)
/// - Cantera library (for chemical kinetics)

///////////////////////////////////////////////////////////////////////////////////////
#include "hips.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo.h"

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>  // For setting precision in output
#include <cmath>    // For isnan and isfinite

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////
int main() {
    // HiPS simulation parameters for non-premixed combustion
    int nLevels = 6;             // Number of hierarchical levels
    double domainLength = 0.01;  // HiPS domain length scale
    double tau0 = 0.0005;        // Mixing timescale
    double C_param = 0.5;        // Eddy rate multiplier
    double tRun = 0.05;          // Simulation runtime
    int forceTurb = 0;           // No forced turbulence
    vector<double> ScHips(54, 1); // Unity Schmidt number for all species

    // Initialize Cantera for combustion chemistry
    auto cantSol = Cantera::newSolution("gri30.yaml");
    auto gas = cantSol->thermo();
    size_t nsp = gas->nSpecies();
    int nVar = nsp + 1; // Number of variables (species + enthalpy)

    // Create HiPS instance with reaction support
    hips HiPS(nLevels, domainLength, tau0, C_param, forceTurb, nVar, ScHips, true
    #ifdef REACTIONS_ENABLED
        , cantSol
    #endif
    );

    int nparcels = HiPS.nparcels;

    // Species mass fractions, temperature, and enthalpy for each parcel
    vector<vector<double>> ysp(nsp, vector<double>(nparcels, 0));
    vector<double> h(nparcels, 0.0), T(nparcels, 0.0), Z(nparcels, 0.0); // Enthalpy, temperature, mixture fraction
    vector<string> variableNames(nsp + 1); // Variable names (enthalpy + species)

    // Debugging: print the number of species
    cout << "Number of species (nsp) = " << nsp << endl;

    // Define temperature profile parameters
    double T_ox = 300.0;  // Oxidizer temperature
    double T_fuel = 1500.0; // Fuel temperature
    double n_temp = 1.0;   // Exponent for interpolation

    // Initialize mixture fraction (Z) from 0 (oxidizer) to 1 (fuel)
    for (int j = 0; j < nparcels; j++) {
        Z[j] = static_cast<double>(j) / (nparcels - 1);
    }

    double min_oxidizer = 0.05; // Minimum oxidizer fraction for stability

    for (int j = 0; j < nparcels; j++) {
        double Zj = Z[j];
        double fuel_fraction = Zj;
        double oxidizer_fraction = max(1.0 - Zj, min_oxidizer); // Prevent zero oxidizer

        // Adjusted composition for fuel and oxidizer
        string composition = "C2H4:" + to_string(fuel_fraction) + 
                             ", O2:" + to_string(oxidizer_fraction * 3.5) +
                             ", N2:" + to_string(oxidizer_fraction * 13.2); // Realistic air composition

        // Calculate initial temperature based on mixture fraction
        double T_init = T_ox + (T_fuel - T_ox) * pow(Zj, n_temp);

        // Set gas state and calculate enthalpy
        gas->setState_TPX(T_init, Cantera::OneAtm, composition);
        h[j] = gas->enthalpy_mass();
        gas->setState_HP(h[j], gas->pressure());
        gas->equilibrate("HP", "auto", 50); // Equilibrate the gas state

        // Store calculated temperature
        T[j] = gas->temperature();

        // Retrieve species mass fractions
        for (int i = 0; i < nsp; i++) {
            ysp[i][j] = gas->massFraction(i);
        }
    }

    // Set temperature and variable names in HiPS
    HiPS.Temp.assign(T.begin(), T.end());
    variableNames[0] = "enthalpy";
    for (int i = 0; i < nsp; i++) {
        variableNames[i + 1] = gas->speciesName(i);
    }

    // Set uniform weight for each parcel
    vector<double> weight(nparcels, 1.0 / nparcels);

    // Assign variables to HiPS
    HiPS.set_varData(h, weight, variableNames[0]);
    for (int k = 0; k < nsp; k++) 
        HiPS.set_varData(ysp[k], weight, variableNames[k + 1]);

    // Run the combustion simulation
    HiPS.calculateSolution(tRun, true);

    return 0;
}

