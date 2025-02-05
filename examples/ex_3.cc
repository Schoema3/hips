////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file ex_3.cc
/// \brief Example 3: HiPS simulation with non-premixed combustion.
///
/// Simulates a non-premixed ethylene/air mixture in a HiPS domain. The mixture fraction varies across 
/// the HiPS tree, with fuel and oxidizer fractions based on this profile. Schmidt numbers are unity for 
/// all species, and the simulation runs for a specified time to observe mixing and combustion phenomena.
////////////////////////////////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Main function for HiPS simulation with non-premixed combustion.
///
/// Sets up and runs the HiPS simulation with non-premixed combustion using Cantera. Initializes the state 
/// of each parcel based on its mixture fraction profile, computes species mass fractions, enthalpy, and 
/// temperature for each parcel, and runs the simulation.
////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
    int            nLevels      = 6;           // Number of HiPS tree levels
    double         domainLength = 0.01;        // HiPS tree length scale
    double         tau0         = 0.0005;      // HiPS tree timescale
    double         C_param      = 0.5;         // HiPS eddy rate multiplier
    double         tRun         = 0.05;        // Simulation runtime
    int            forceTurb    = 0;           // Force turbulent profile (0 = no)
    vector<double> ScHips(54, 1);              // Schmidt number (unity for all species)

    auto cantSol = Cantera::newSolution("gri30.yaml");
    auto gas = cantSol->thermo();
    size_t nsp = gas->nSpecies();
    int nVar = nsp + 1; // Number of variables (species + enthalpy)

    hips HiPS(nLevels, domainLength, tau0, C_param, forceTurb, nVar, ScHips, true
#ifdef REACTIONS_ENABLED
              , cantSol
#endif
             );

    int nparcels = HiPS.nparcels;

    vector<vector<double>> ysp(nsp, vector<double>(nparcels, 0)); // Species mass fractions
    vector<double> h(nparcels, 0.0), T(nparcels, 0.0), Z(nparcels, 0.0); // Enthalpy, temperature, mixture fraction
    vector<string> variableNames(nsp + 1); // Variable names (enthalpy + species)

    cout << "Number of species (nsp) = " << nsp << endl; // Debugging species count

    // Temperature profile settings
    double T_ox = 300.0;  // Oxidizer temperature
    double T_fuel = 1500.0; // Fuel temperature
    double n_temp = 1.0;   // Exponent for interpolation

    for (int j = 0; j < nparcels; j++) {
        Z[j] = static_cast<double>(j) / (nparcels - 1); // Z ranges from 0 (oxidizer) to 1 (fuel)
    }

    double min_oxidizer = 0.05; // Ensure at least 5% oxidizer for stability

    for (int j = 0; j < nparcels; j++) {
        double Zj = Z[j];
        double fuel_fraction = Zj;
        double oxidizer_fraction = max(1.0 - Zj, min_oxidizer); // Prevent zero oxidizer

        // Adjusted stoichiometric oxidizer and nitrogen amounts
        string composition = "C2H4:" + to_string(fuel_fraction) + 
                             ", O2:" + to_string(oxidizer_fraction * 3.5) +
                             ", N2:" + to_string(oxidizer_fraction * 13.2); // Maintain realistic air composition

        // Compute initial temperature based on mixture fraction
        double T_init = T_ox + (T_fuel - T_ox) * pow(Zj, n_temp);

        gas->setState_TPX(T_init, Cantera::OneAtm, composition);

        // Get enthalpy after setting temperature
        h[j] = gas->enthalpy_mass();
        gas->setState_HP(h[j], gas->pressure());
        gas->equilibrate("HP", "auto", 50);

        T[j] = gas->temperature();

        // Properly assign species fractions in a nested loop
        for (int i = 0; i < nsp; i++) {
            ysp[i][j] = gas->massFraction(i);
        }
    }

    HiPS.Temp.assign(T.begin(), T.end());

    variableNames[0] = "enthalpy";
    for (int i = 0; i < nsp; i++) {
        variableNames[i + 1] = gas->speciesName(i);
    }

    vector<double> weight(nparcels, 1.0 / nparcels); // Set weight for each parcel

    HiPS.set_varData(h, weight, variableNames[0]);
    for (int k = 0; k < nsp; k++) 
        HiPS.set_varData(ysp[k], weight, variableNames[k + 1]);

    HiPS.calculateSolution(tRun, true);

    return 0;
}

