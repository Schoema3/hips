////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file ex_2.cc
/// \brief Example 2: HiPS simulation with chemical reactions focusing on combustion processes.
///
/// This example simulates a premixed stoichiometric ethylene/air mixture in a HiPS domain.
/// The HiPS tree is initialized with six levels, where 25% of parcels are pre-combusted, 
/// and 75% contain fresh reactants. The Schmidt number is set to unity for all species, 
/// and the domain length scale is 0.01. The simulation is run for different values of 
/// mixing timescales to observe the impact on ignition and extinction phenomena.
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
/// This function sets up and runs a HiPS simulation with chemical reactions using Cantera.
/// It initializes the HiPS tree, sets the gas state, and runs the simulation for a specified
/// time period. The results include species mass fractions, enthalpy, and temperature for each parcel.
///
/// \return 0 on successful execution.
////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Set up the HiPS tree parameters and constructor.
    ///
    /// Initializes the HiPS tree with parameters such as the number of levels, domain length scale,
    /// mixing timescale, eddy rate multiplier, and Schmidt numbers. The HiPS tree is used to simulate
    /// turbulent mixing coupled with chemical reactions.
    ///
    /// \param nLevels        Number of HiPS tree levels (hierarchical levels).
    /// \param domainLength   Length scale of the HiPS domain.
    /// \param tau0           Timescale for the HiPS tree.
    /// \param C_param        Multiplier for the eddy rate in the HiPS tree.
    /// \param forceTurb      Flag to force a turbulent profile (set to 0 for no forcing).
    /// \param ScHips         Vector of Schmidt numbers for each species (unity for all species).
    /// \param tRun           Simulation runtime.
    /// \param nparcels       Number of parcels in the HiPS tree.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    int            nLevels      = 6;           // Number of HiPS tree levels
    double         domainLength = 0.01;        // HiPS tree lengthscale
    double         tau0         = 0.0005;      // HiPS tree timescale
    double         C_param      = 0.5;         // HiPS eddy rate multiplier
    double         tRun         = 0.05;        // Tree evolution time
    int            forceTurb    = 0;           // Force turbulent profile
    vector<double> ScHips(54,1);               // Schmidt number is unity for all species

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Set up the gas solution using Cantera.
    ///
    /// Initializes the Cantera gas solution using the "gri30.yaml" mechanism. This involves setting
    /// up the thermodynamic state of the gas and retrieving the number of species. The number of 
    /// variables (nVar) is determined by the number of species plus one for enthalpy.
    ///
    /// \param cantSol  Pointer to the Cantera solution object.
    /// \param nsp      Number of species in the Cantera gas solution.
    /// \param nVar     Number of variables in the simulation (enthalpy + species).
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

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

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Declare state vectors for species mass fractions, enthalpy, and temperature.
    ///
    /// This section sets up vectors to store the species mass fractions, enthalpy, and temperature 
    /// for each parcel in the HiPS tree. Additionally, initial state vectors are prepared for fresh 
    /// reactants and pre-combusted parcels.
    ///
    /// \param ysp   2D vector to store species mass fractions for each parcel.
    /// \param h     Vector to store enthalpy for each parcel.
    /// \param T     Vector to store temperature for each parcel.
    /// \param y0    Initial state vector for fresh reactants.
    /// \param h0    Initial enthalpy for fresh reactants.
    /// \param T0    Initial temperature for fresh reactants.
    /// \param y1    Initial state vector for pre-combusted parcels.
    /// \param h1    Initial enthalpy for pre-combusted parcels.
    /// \param T1    Initial temperature for pre-combusted parcels.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    
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

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Set state vectors for each parcel in the HiPS tree.
    ///
    /// This section initializes the state of each parcel based on whether it is pre-combusted or contains 
    /// fresh reactants. The state includes the temperature, enthalpy, and species mass fractions.
    ///
    /// \param fracBurn  Fraction of parcels that are pre-combusted.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

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
        }
    }

    HiPS.Temp = T;

    variableNames[0] = "enthalpy";
    for(int i = 0; i < ysp.size(); i++)
        variableNames[i + 1] = gas->speciesName(i);

    vector<double> weight(nparcels, 1.0 / nparcels); // Set weight for each parcel

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Set variable data in the HiPS tree.
    ///
    /// This section sets the variable data (enthalpy and species mass fractions) in the HiPS tree 
    /// for each parcel, using the previously initialized state vectors.
    ///
    /// \param weight          Vector of weights for each parcel.
    /// \param variableNames   Vector of variable names (enthalpy and species).
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    HiPS.set_varData(h, weight, variableNames[0]);
    for (int k = 0; k < ysp.size(); k++) 
        HiPS.set_varData(ysp[k], weight, variableNames[k + 1]);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Advance the HiPS tree to simulate reaction and mixing.
    ///
    /// This final section advances the HiPS tree to perform the chemical reactions and mixing 
    /// for the specified runtime.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    HiPS.calculateSolution(tRun, true);

    return 0;
}

