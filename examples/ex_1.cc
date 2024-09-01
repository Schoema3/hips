
////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file ex_1.cc
/// \brief Example simulation of scalar mixing dynamics using HiPS.
///
/// This example demonstrates scalar mixing dynamics, highlighting the effects of 
/// different Schmidt numbers on the mixing process. Initially, fluid parcels are 
/// divided into two groups with scalar values of 0 and 1. These values mix over time 
/// to reach an average value, showcasing the impact of high and low Schmidt numbers.
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include "hips.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Initializes mixing fractions for fluid parcels.
///
/// Splits fluid parcels into two groups:
/// - First half assigned a scalar value of 0.
/// - Second half assigned a scalar value of 1.
///
/// \param numParcels Number of fluid parcels.
/// \return A vector of mixing fractions for each parcel.
/////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<double> initializeMixingFractions(int numParcels) {
    vector<double> mixingFractions;
    for (int i = 0; i < numParcels; i++) {
        mixingFractions.push_back(i < numParcels / 2 ? 0.0 : 1.0);
    }
    cout << "\n\n";
    return mixingFractions;
}

std::vector<std::string> variableNames = {"mixf_00", "mixf_01"};

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Runs the HiPS simulation for scalar mixing dynamics.
///
/// Initializes the HiPS tree, assigns mixing fractions, and runs the simulation 
/// over a specified time to observe the mixing process.
///
/// \return 0 on successful execution.
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
    // HiPS tree setup parameters
    int nLevels = 9;
    double domainLength = 1.0;
    double tau0 = 1.0;
    double C_param = 0.5;
    double tRun = 300.0;
    int forceTurb = 2;
    vector<double> ScHips = {0.0625, 16.0};
    int numVariables = 2;

    // Initialize HiPS tree
    hips HiPS(nLevels, domainLength, tau0, C_param, forceTurb, numVariables, ScHips, false);
    int numParcels = HiPS.nparcels;

    // Initialize and set mixing fractions for each parcel
    vector<vector<double>> tot_vec;
    vector<double> weights(numParcels, 1.0 / numParcels);

    for (int i = 0; i < variableNames.size(); ++i) 
        tot_vec.push_back(initializeMixingFractions(numParcels));
    
    for (int i = 0; i < tot_vec.size(); i++) 
        HiPS.set_varData(tot_vec[i], weights, variableNames[i]);

    // Run the simulation
    HiPS.calculateSolution(tRun, true);

    return 0;
}

