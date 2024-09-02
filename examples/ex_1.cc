
////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file ex_1.cc
/// \brief Example 1: HiPS simulation of scalar mixing dynamics.
///
/// This example demonstrates scalar mixing dynamics in a turbulent flow using HiPS.
/// The simulation initializes fluid parcels with two distinct scalar values (0 and 1),
/// representing an initially unmixed state. The scalar values mix over time, showcasing
/// the effects of different Schmidt numbers on the mixing process.
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include "hips.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Initializes mixing fractions for fluid parcels.
///
/// This function divides the fluid parcels into two groups:
/// - The first half is assigned a scalar value of 0.
/// - The second half is assigned a scalar value of 1.
///
/// The function returns a vector containing the initial scalar values (mixing fractions)
/// for each fluid parcel in the simulation.
///
/// \param numParcels The total number of fluid parcels.
/// \return A vector of mixing fractions corresponding to each parcel.
/////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<double> initializeMixingFractions(int numParcels) {
    vector<double> mixingFractions(numParcels);
    
    for (int i = 0; i < numParcels; i++) {
        // Assign 0.0 to the first half and 1.0 to the second half
        mixingFractions[i] = (i < numParcels / 2) ? 0.0 : 1.0;
    }
    return mixingFractions;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Main function for the HiPS simulation of scalar mixing.
///
/// This function sets up and runs a HiPS simulation to observe scalar mixing dynamics.
/// It initializes the HiPS tree, assigns initial scalar values (mixing fractions) to the
/// fluid parcels, and runs the simulation over a specified period.
///
/// The simulation is designed to highlight the impact of different Schmidt numbers on the
/// scalar mixing process.
///
/// \return 0 on successful execution.
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Set up the HiPS tree parameters and constructor.
    ///
    /// Initializes the HiPS tree with parameters such as the number of hierarchical levels,
    /// domain length, initial timescale, eddy rate multiplier, Schmidt numbers, and the
    /// number of variables being tracked. The HiPS tree is used to simulate scalar mixing
    /// in a turbulent flow environment.
    ///
    /// \param nLevels        Number of hierarchical levels in the HiPS tree.
    /// \param domainLength   Physical length scale of the simulation domain.
    /// \param tau0           Initial timescale associated with the largest eddies.
    /// \param C_param        Multiplier for the eddy turnover rate.
    /// \param forceTurb      Flag to force a turbulent profile (set to 2 for forcing).
    /// \param ScHips         Vector of Schmidt numbers for each scalar being tracked.
    /// \param numVariables   Number of variables (e.g., scalar fields) being simulated.
    /// \param tRun           Total runtime for the simulation.
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    int nLevels = 9;                           // Number of levels in the HiPS tree
    double domainLength = 1.0;                 // Length scale of the HiPS domain
    double tau0 = 1.0;                         // Initial timescale
    double C_param = 0.5;                      // Eddy turnover rate multiplier
    double tRun = 300.0;                       // Simulation runtime
    int forceTurb = 2;                         // Force turbulent profile
    vector<double> ScHips = {0.0625, 16.0};    // Schmidt numbers for each scalar
    int numVariables = 2;                      // Number of variables (scalars)

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Initialize the HiPS tree.
    ///
    /// The HiPS object is created with the specified parameters, representing the physical
    /// and turbulent characteristics of the simulation domain. The number of parcels in the
    /// simulation is determined based on the HiPS tree configuration.
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    hips HiPS(nLevels, domainLength, tau0, C_param, forceTurb, numVariables, ScHips, false);
    int numParcels = HiPS.nparcels; // Number of fluid parcels in the HiPS tree

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Initialize and assign mixing fractions to fluid parcels.
    ///
    /// The scalar fields (mixing fractions) are initialized for each parcel, and the HiPS
    /// tree is updated with these values. Two scalar fields are simulated, corresponding to
    /// two different Schmidt numbers.
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    // Initialize vectors to store mixing fractions for each scalar variable
    vector<vector<double>> mixingFractions(numVariables);

    // Initialize vector for parcel weights (assumed uniform)
    vector<double> weights(numParcels, 1.0 / numParcels);

    // Initialize and set mixing fractions for each scalar field
    for (int i = 0; i < numVariables; ++i) {
        mixingFractions[i] = initializeMixingFractions(numParcels);
        HiPS.set_varData(mixingFractions[i], weights, "mixf_0" + to_string(i));
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Run the HiPS simulation to observe scalar mixing.
    ///
    /// The HiPS tree is advanced in time to simulate the mixing process, where the scalar fields
    /// evolve according to the specified Schmidt numbers. The results demonstrate how scalars
    /// with different Schmidt numbers mix at different rates in a turbulent flow.
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    HiPS.calculateSolution(tRun, true);

    return 0;
}

