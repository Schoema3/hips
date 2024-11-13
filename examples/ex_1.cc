
////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file ex_1.cc
/// \brief HiPS simulation example for scalar mixing in turbulent flow.
/// Initializes fluid parcels with distinct scalar values and observes mixing influenced by Schmidt numbers.
////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include "hips.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Initializes fluid parcel scalar values (mixing fractions) as 0 and 1.
/// \param numParcels Total number of fluid parcels.
/// \return Vector of initial scalar values for each parcel.
////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<double> initializeMixingFractions(int numParcels) {
    vector<double> mixingFractions(numParcels);
    
    for (int i = 0; i < numParcels; i++) {
        // Assign 0.0 to the first half and 1.0 to the second half
        mixingFractions[i] = (i < numParcels / 2) ? 0.0 : 1.0;
    }
    return mixingFractions;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Runs the HiPS simulation for scalar mixing dynamics.
/// Sets up HiPS tree parameters, initializes scalar values, and runs the simulation.
////////////////////////////////////////////////////////////////////////////////////////////////////////

int main() {
    // HiPS tree and simulation parameters
    int nLevels = 9;                           // Hierarchical levels
    double domainLength = 1.0;                 // Simulation domain length
    double tau0 = 1.0;                         // Initial timescale for largest eddies
    double C_param = 0.5;                      // Eddy turnover rate multiplier
    double tRun = 300.0;                       // Total simulation time
    int forceTurb = 2;                         // Forcing turbulent profile
    vector<double> ScHips = {0.0625, 16.0};    // Schmidt numbers for scalars
    int numVariables = 2;                      // Number of scalar fields

    // Set up HiPS tree and parcel count
    hips HiPS(nLevels, domainLength, tau0, C_param, forceTurb, numVariables, ScHips, false);
    int numParcels = HiPS.nparcels; 

    // Initialize mixing fractions for each scalar
    vector<vector<double>> mixingFractions(numVariables);
    vector<double> weights(numParcels, 1.0 / numParcels);

    for (int i = 0; i < numVariables; ++i) {
        mixingFractions[i] = initializeMixingFractions(numParcels);
        HiPS.set_varData(mixingFractions[i], weights, "mixf_0" + to_string(i));
    }

    // Run simulation to observe mixing dynamics
    HiPS.calculateSolution(tRun, true);

    return 0;
}

