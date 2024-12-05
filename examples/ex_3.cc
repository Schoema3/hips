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

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Main function for HiPS simulation with non-premixed combustion.
///
/// Sets up and runs the HiPS simulation with non-premixed combustion using Cantera. Initializes the state 
/// of each parcel based on the mixture fraction profile, computes species mass fractions, enthalpy, and 
/// temperature for each parcel, and runs the simulation.
////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
    int            nLevels      = 6;           // Number of HiPS tree levels
    double         domainLength = 0.01;        // HiPS tree length scale
    double         tau0         = 0.0005;      // HiPS tree timescale
    double         C_param      = 0.5;         // HiPS eddy rate multiplier
    double         tRun         = 0.05;        // Simulation runtime
    int            forceTurb    = 0;           // Force turbulent profile (0 = no)
    vector<double> ScHips(54, 1);               // Schmidt number (unity for all species)

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
    vector<double> h(nparcels), T(nparcels), Z(nparcels);         // Enthalpy, temperature, and mixture fraction
    vector<string> variableNames(nsp + 1);                         // Variable names (enthalpy + species)

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Define mixture fraction profile.
    ///
    /// Defines a linear mixture fraction profile for the parcels, ranging from 0 (oxidizer) to 1 (fuel).
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    for (int j = 0; j < nparcels; j++) {
        Z[j] = static_cast<double>(j) / (nparcels - 1); // Z ranges from 0 (oxidizer) to 1 (fuel)
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Initialize parcels based on mixture fraction.
    ///
    /// Initializes the state of each parcel based on its mixture fraction. Fuel and oxidizer fractions 
    /// are computed from the mixture fraction, and the gas state is set accordingly for each parcel.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    for (int j = 0; j < nparcels; j++) {
        double Zj = Z[j];
        double fuel_fraction = Zj;
        double oxidizer_fraction = 1.0 - Zj;

        string composition = "C2H4:" + to_string(fuel_fraction) + ", O2:" + to_string(oxidizer_fraction * 3.0) +
                             ", N2:" + to_string(oxidizer_fraction * 11.52); // Air composition
        gas->setState_TPX(400.0, Cantera::OneAtm, composition);
        gas->equilibrate("HP");

        // Store parcel properties
        h[j] = gas->enthalpy_mass();
        T[j] = gas->temperature();
        gas->getMassFractions(ysp[j].data());
    }

    HiPS.Temp = T;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Set variable names and HiPS data.
    ///
    /// Assigns variable names and sets the HiPS data (enthalpy and species mass fractions) for each parcel.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    variableNames[0] = "enthalpy";
    for (int i = 0; i < ysp.size(); i++) {
        variableNames[i + 1] = gas->speciesName(i);
    }

    vector<double> weight(nparcels, 1.0 / nparcels); // Set weight for each parcel

    HiPS.set_varData(h, weight, variableNames[0]);
    for (int k = 0; k < ysp.size(); k++) 
        HiPS.set_varData(ysp[k], weight, variableNames[k + 1]);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Advance the HiPS tree.
    ///
    /// Advances the HiPS tree to simulate chemical reactions and mixing for the specified runtime.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    HiPS.calculateSolution(tRun, true);

    return 0;
}

