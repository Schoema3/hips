#include <iostream>
#include <vector>
#include "hips.h"

using namespace std;

// Function to initialize mixing fractions
vector<double> initializeMixingFractions(int numParcels) {
    vector<double> mixingFractions;
    for (int i = 0; i < numParcels; i++) {
        mixingFractions.push_back(i < numParcels / 2 ? 0.0 : 1.0);
    }
    cout << "\n\n";
    return mixingFractions;
}

std::vector<std::string> variableNames = {"mixf_00", "mixf_01", "mixf_02"};

int main() {
    // HiPS tree and constructor parameters
    int nLevels = 6;
    double domainLength = 1.0;
    double tau0 = 1.0;
    double C_param = 0.5;
    double tRun = 300.0;
    int forceTurb = 2;
    vector<double> ScHips;
    
    ScHips.push_back(0.25);
    ScHips.push_back(1.0);
    ScHips.push_back(4.0);

    int numVariables = 3;

    hips HiPS(nLevels, domainLength, tau0, C_param, forceTurb, numVariables, ScHips, false);
    int numParcels = HiPS.nparcels;

    // Initialize mixing fractions for each parcel
    vector<vector<double>> tot_vec;
    vector<double> weights(numParcels, 1.0 / numParcels);  // Initialize weights once

    for (int i = 0; i < variableNames.size(); ++i) 
        tot_vec.push_back(initializeMixingFractions(numParcels));
    

    // Set state vectors in each parcel
    for (int i = 0; i < tot_vec.size(); i++) 
        HiPS.set_varData(tot_vec[i], weights, variableNames[i]);

    HiPS.calculateSolution(tRun, true);

    return 0;
}

