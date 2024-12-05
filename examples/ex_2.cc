////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file ex_2.cc
/// \brief Example 2: HiPS simulation of ethylene/air combustion.
///
/// Simulates a premixed stoichiometric ethylene/air mixture in a HiPS domain with six levels.
/// 25% of parcels are pre-combusted and 75% contain fresh reactants. Schmidt number is set to unity
/// for all species, with a domain length scale of 0.01. The simulation runs for various mixing timescales
/// to observe the effects on ignition and extinction.
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
/// \brief Main function for HiPS simulation with chemical reactions.
///
/// Sets up and runs a HiPS simulation using Cantera, initializes the HiPS tree, sets the gas state,
/// and performs the simulation. Results include species mass fractions, enthalpy, and temperature for each parcel.
///
/// \return 0 on success.
////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Initialize HiPS tree parameters.
    ///
    /// Configures the HiPS tree with parameters such as levels, domain length, mixing timescale,
    /// Schmidt numbers, and runtime. The HiPS tree models turbulent mixing with chemical reactions.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    int            nLevels      = 6;           // Number of HiPS tree levels
    double         domainLength = 0.01;        // HiPS domain length scale
    double         tau0         = 0.0005;      // HiPS timescale
    double         C_param      = 0.5;         // Eddy rate multiplier
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

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Set up Cantera gas solution.
    ///
    /// Initializes Cantera gas solution using "gri30.yaml", sets thermodynamic state, and retrieves species count.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    vector<vector<double>> ysp(53, vector<double>(nparcels, 0)); // Species mass fractions for each parcel
    vector<double> h(nparcels), T(nparcels);                     // Enthalpy and temperature for each parcel

    vector<double> y0(nsp), y1(nsp);                              // Initial species mass fractions for reactants
    double T0 = 300.0, T1 = 300.0;                                // Initial temperature
    double h0, h1;                                               // Enthalpy for fresh reactants and pre-combusted parcels

    vector<string> variableNames(nsp + 1);                        // Variable names (enthalpy and species)

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Initialize parcel states.
    ///
    /// Sets the state for fresh reactants and pre-combusted parcels, using the fraction of pre-combusted parcels.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    double fracBurn = 0.25; // Fraction of parcels pre-combusted
    gas->setState_TPX(T0, Cantera::OneAtm, "C2H4:1, O2:3, N2:11.25");
    h0 = gas->enthalpy_mass();
    gas->getMassFractions(y0.data());

    // Initialize fresh reactants
    for (int j = 0; j < fracBurn * nparcels; j++) {
        h[j] = h0;
        T[j] = gas->temperature();
        ysp[j] = y0;
    }

    // Initialize pre-combusted parcels
    gas->setState_TPX(T1, Cantera::OneAtm, "C2H4:1, O2:3, N2:11.25");
    h1 = gas->enthalpy_mass();
    gas->setState_HP(h1, gas->pressure());
    gas->equilibrate("HP");
    gas->getMassFractions(y1.data());

    for (int j = fracBurn * nparcels; j < nparcels; j++) {
        h[j] = h1;
        T[j] = gas->temperature();
        ysp[j] = y1;
    }

    HiPS.Temp = T;
    variableNames[0] = "enthalpy";
    for (int i = 0; i < ysp.size(); i++) {
        variableNames[i + 1] = gas->speciesName(i);
    }

    vector<double> weight(nparcels, 1.0 / nparcels); // Parcel weight

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Set variable data in HiPS tree.
    ///
    /// Sets the variable data (enthalpy and species mass fractions) for each parcel in the HiPS tree.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

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

