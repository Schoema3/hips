
#ifdef REACTIONS_ENABLED
#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#endif
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

    // for (int i =0; i<7;i++)
    //    cout<<"first vat "<<mixingFractions[i] <<endl;

    return mixingFractions;
}
    
vector<double> initializeMixingFractions_2(int numParcels) {
    vector<double> mixingFractions_2;
    for (int i = 0; i < numParcels; i++) {
        mixingFractions_2.push_back(i < numParcels / 2 ? -1000.0 : -500.0);
     //cout<<"i is "<<i<<"     "<<   mixingFractions_2[i]<<endl;

    }
   cout<<"\n\n";
   return mixingFractions_2;
}
   
std::vector<std::string> variableNames={"mixf_00","mixf_01"};

int main() {
    // HiPS tree and constructor parameters
    //hips Hi(200.08);
    int numLevels = 4;
    double domainLength = 1.0;
    double tau0 = 1.0;
    double C_param = 0.5;
    double tRun = 30.0;
    int forceTurb = 0;
    vector<double> ScHips(2, 1 );

    // Gas solution setup

#ifdef REACTIONS_ENABLED
    auto cantSol = Cantera::newSolution("gri30.yaml");
    auto gas = cantSol->thermo();
    size_t numSpecies = gas->nSpecies();
 #endif
  

    int numVariables = 2;

    // HiPS tree creation
    hips HiPS(numLevels, domainLength, tau0, C_param, forceTurb, numVariables, ScHips,
              #ifdef REACTIONS_ENABLED
                 cantSol,
              #endif
              false);
    int numParcels = HiPS.nparcels;


    cout<<"number of parcels --> "<<numParcels<<endl;
    // Initialize mixing fractions
    vector<double> mixingFractions = initializeMixingFractions(8);
    vector<double> mixingFractions_2 = initializeMixingFractions_2(8);
    vector<double> weight(8, 1.0/8);


   // for (int i =0; i<weight.size();i++)
   //     cout<<"weight "<<weight[i]<<endl;


    // Create a vector to hold all mixing fractions
    vector<vector<double>> tot_vec;
    tot_vec.push_back(mixingFractions);
    tot_vec.push_back(mixingFractions_2);

    // Set state vectors in each parcel
    for (int i = 0; i < tot_vec.size(); i++) 
        HiPS.set_varData(tot_vec[i],weight, variableNames[i], i);

    HiPS.calculateSolution(tRun);

    for (int i = 0; i < tot_vec.size(); i++) 
        HiPS.get_varData();

    return 0;
}

