#pragma once

#ifdef  REACTIONS_ENABLED
    #include "cantera/base/Solution.h"
    #include "cantera/thermo.h"
    #include "batchReactor.h"
#endif

#include <iostream>
#include <vector>
#include <string>
#include <vector>
#include "randomGenerator.h"

class hips {

    ////////////////////////////// DATA MEMBERS: PUBLIC /////////////////////////////

public:
  
    int nparcels;                                     ///< number of parcels
    std::vector<std::vector<double>*> varData;        ///< vector of pointers to vector
 
    #ifdef REACTIONS_ENABLED
        std::shared_ptr<Cantera::ThermoPhase> gas;    ///< Shared pointer to a Cantera thermochemistry object
        std::unique_ptr<batchReactor> bRxr;           ///< Unique pointer to the integrator object
    #endif

    double domainLength;                              ///< length of domain (m)
    double tau0;                                      ///< integral timescale
    double C_param;                                   ///< Eddy frequency parameter
    std::vector<double> Temp;                         ///< Vector containg temperature in each parcel;
    
    /////////////////////////// DATA MEMBERS: PRIVATE //////////////////////////////////

private:
    
    int currentIndex = 0;                              ///< member variable to keep track of current index of variables
    int nLevels;                                       ///< number of tree levels
    int nLevels_;                                      ///< number of tree levels?
    int forceTurb;                                     ///< forcing function for statistically stationary: -1 = none, 1 = source term, 2 = dir
    int nVar;                                          ///< number of parcel variables (e.g., h, ysp)
    int nsp;                                           ///< number of species
    int Nm1;                                           ///< nLevels - 1 
    int Nm2;                                           ///< nLevels - 2                     
    int Nm3;                                           ///< nLevels - 3 
    int iEta;                                          ///< Kolmogorov level (needed for variable Sc scalars)
    
    bool LScHips;                                      ///< hips schmidt number
    bool performReaction;                              ///< flag indicating whether chemical reactions are performed in the simulation 
    bool LrandSet;                                     ///< flag indicating new randomGen  --> allow deletion
    
    double time;                                       ///< current simulation time
    double eddyRate_total;                             ///< total rate of all eddies 0 through nLevels-3
    double eddyRate_inertial;                          ///< total rate of all eddies 0 through iEta (= eddyRate_total if Sc=1) 
    double Afac = 0.5;                                 ///< level lengthscale reduction factor (0.5)
    double Re;                                         ///< Reynolds number
    double dtEE;                                       ///< time increment to next eddy event 
    
    randomGenerator rand;                                   
    
    std::vector<int> i_plus;                            ///< ceil(i_batchelor)
    std::vector<int> pLoc;                              ///< parcel index array for fast implementation of swaps
    std::vector<double> varRho;                           
    std::vector<double> ScHips;                         ///< vector containing Schmidt numbers related to each variable
    std::vector<std::string> varName;                   ///< vector containing the names of parcel variables
    std::vector<double> parcelTimes;                    ///< current times corresponding to the parcel states
    std::vector<double> levelRates;                     ///< list of eddy event rates at each level
    std::vector<double> i_batchelor;                    ///< Batchelor level for variable Sc scalars; NOTE: double, as in, between levels
    std::vector<double> xc;                             ///< vector containing physical domain of flow particles
    std::vector<double> xh;                             ///< vector containing physical domain of HiPS parcels
    
      /////////////////////////  STATIC MEMBERS  ///////////////////////// 

    static int nL;                                      ///< adjusted number of levels based on the Reynolds number
    static double Prob;                                 ///< probability value for probability-based solution
    static double lStar;                                ///< length of the level associated with the Reynolds number 
    static double Anew;                                 ///< adjusted level lengthscale reduction factor for dynamic adjustment of reduction factor

      ////////////////////////////// MEMBER FUNCTIONS: PUBLIC /////////////////////////////

public:
    
         /////////////////////////  COSTRUCTORS  ///////////////////////// 
    hips(double C_param_,
         int forceTurb_,
         int nVar_,
         bool performReaction,
         #ifdef REACTIONS_ENABLED
             std::shared_ptr<Cantera::Solution> cantSol = nullptr,
         #endif
         int seed = 10);


    hips(int nLevels_,
         double domainLength_,
         double tau0_,
         double C_param_,
         int forceTurb_,
         int nVar_,
         std::vector<double> &ScHips_,
         bool performReaction,
         #ifdef REACTIONS_ENABLED
             std::shared_ptr<Cantera::Solution> cantSol = nullptr,
         #endif
         int seed = 10);
 
    virtual ~hips() {
        for(auto& data : varData)
            delete data;
    }                                                                                               // Destructor for the hips class

    void set_tree(int nLevels_, double domainLength_, double tau0_, std::vector<double> &ScHips_);  // setting the HiPS tree based on the number of levels
    void  set_tree(double Re_, double domainLength_, double tau0_, std::vector<double> &ScHips_, std::string approach);
   
    void set_varData(std::vector<double> &v, std::vector<double> &w, const std::string &varN);      // passing all variables to vector of pointer 
    void set_varData(std::vector<double> &v, std::vector<double> &w, const std::string &varN,
                     const std::vector<double> &rho);                                               // passing all variables to vector of pointer

    std::vector<std::vector<double>> get_varData();                                                 // Retrieves modified data from the HiPS library and stores it in the provided vector. 
    
    void calculateSolution(const double tRun, bool shouldWriteData =false);                         // Running simulations 

private:

    ////////////////////////////// MEMBER FUNCTIONS: PRIVATE /////////////////////////////


    std::vector<double> projection(std::vector<double> &vcfd, std::vector<double> &weight);           //Perform vector projection of flow particles onto hips parcels operation without density 
    
    std::pair<std::vector<double>, std::vector<double>>  projection(std::vector<double> &vcfd, std::vector<double> &weight,                     
                                   const std::vector<double> &density);                             // Perform vector projection flow particles onto hips parcels operation with density 
    
    std::vector<double> setGridHips(int N);                                                          // Set Hips grid with a specified number of grid points equal to number of parcels   
    std::vector<double> setGridCfd(std::vector<double> &w);                                          // Set CFD grid using provided weight vector
    std::vector<double>  projection_back(std::vector<double> &vb);                                   // Perform vector projection of hips parcels onto flow particles operation without density
    
    void sample_hips_eddy(double &dt, int &iLevel);                                                  // Sample hips eddy with specified time step and level                                                
    void selectAndSwapTwoSubtrees(const int iLevel, int &iTree);                                     // Select and swap two subtrees in the level tree
    void advanceHips(const int iLevel, const int iTree);                                             // Advancing simulations to do mixing and reaction
   
    void reactParcels_LevelTree(const int iLevel, const int iTree);                                  // Reacting parcels involved in micro-mixing
    void mixAcrossLevelTree(int kVar, const int iMixLevel, const int iTree);                         // Mixing paecels involved in micr0-mixing.
   
    void forceProfile();
    
    void writeData(const int ifile, const double outputTime);                                        // Writing the results for a user-defined number of eddies in the data folder.
    void writeInputParameters();
};



















