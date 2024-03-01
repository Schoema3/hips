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
public:
    // DATA MEMBERS
    int nparcels;                                                       ///< number of parcels
    std::vector<std::vector<double>*> varData;                          ///< vector of pointers to vector
    #ifdef REACTIONS_ENABLED
        std::shared_ptr<Cantera::ThermoPhase> gas;
        std::unique_ptr<batchReactor> bRxr;
    #endif

    // MEMBER VARIABLES
    int nL;                                                              ///< adjusted number of levels based on Reynolds.
    double domainLength;                                                 ///< length of domain (m)
    double tau0;                                                         ///< integral timescale
    double C_param;                                                      ///< Eddy frequency parameter
    double time;                                                         ///< current simulation time
    double eddyRate_total;                                               ///< total rate of all eddies 0 through nLevels-3
    double eddyRate_inertial;                                            ///< ?????
    double Afac;                                                         ///< level lengthscale reduction factor (0.5)

private:
    // Private data members
    int nLevels;                                                         ///< number of tree levels
    int nLevels_;                                                        ///< number of tree levels?
    int forceTurb;                                                       ///< forcing function for statistically stationary: -1 = none, 1 = source term, 2 = dir
    int nVar;                                                            ///< number of parcel variables (e.g., h, ysp)
    int nsp;                                                             ///< number of species
    int Nm1;                                                             ///< nLevels - 1 
    int Nm2;                                                             ///< nLevels - 2                     
    int Nm3;                                                             ///< nLevels - 3 
    int iEta;                                                            ///< Kolmogorov level (needed for variable Sc scalars)
    bool LScHips;                                                        ///< hips schmidt number
    bool performReaction;                                                ///< Flag indicating whether chemical reactions are performed in the simulation 
    bool LrandSet;                                                       ///< flag indicating new randomGen  --> allow deletion
    randomGenerator rand;
    std::vector<int> i_plus;                                             ///< ceil(i_batchelor)
    std::vector<int> pLoc;                                               ///< parcel index array for fast implementation of swaps
    std::vector<double> varRho;
    std::vector<double> ScHips;                                          ///< Vector containing Schmidt numbers related to each variable
    std::vector<std::string> varName;                                    ///< Vector containing the names of parcel variables
    std::vector<double> parcelTimes;                                     ///< current times corresponding to the parcel states
    std::vector<double> levelRates;                                      ///< list of eddy event rates at each level
    std::vector<double> i_batchelor;                                     ///< Batchelor level for variable Sc scalars; NOTE: double, as in, between levels

public:
    // MEMBER FUNCTIONS

    /**
    * @brief Constructor for the hips class.
    * Initializes a hips object with specified parameters.
    */            
    hips(int nLevels_,
         double domainLength_,
         double tau0_,
         double C_param_,
         int forceTurb_,
         int nVar_,
         std::vector<double> &ScHips_,
         #ifdef REACTIONS_ENABLED
            std::shared_ptr<Cantera::Solution> cantSol,
         #endif
         bool performReaction,
         int seed = 10);

    hips(double Re);                                                       ///< Constructor for the hips class with a Reynolds number parameter

    virtual ~hips() {
        for(auto& data : varData)
            delete data;
    }                                                                      ///< Destructor for the hips class

    void set_varData(std::vector<double> &v, std::vector<double> &w, const std::string &varN, int i);    ///< passing all variables to vector of pointer 
    void set_varData(std::vector<double> &v, std::vector<double> &w, const std::string &varN,
                     const std::vector<double> &rho, int i);                                             ///< passing all variables to vector of pointer                                            
    std::vector<double> projection(std::vector<double> &vcfd, std::vector<double> &weight);              ///< Perform vector projection operation without density 

    std::pair<std::vector<double>, std::vector<double>>  projection(std::vector<double> &vcfd, std::vector<double> &weight,                     
                                   const std::vector<double> &density);                                   ///< Perform vector projection operation with density 

    std::vector<double> setGridHips(int N);                                  ///< Set Hips grid with a specified number of grid points equal to number of parcels   

    std::vector<double> setGridCfd(std::vector<double> &w);                  ///< Set CFD grid using provided weight vector

    std::vector<std::vector<double>> get_varData();                          ///< Retrieves modified data from the HiPS library and stores it in the provided vector. 

    std::vector<double>  projection_back(std::vector<double> &vb);

    void calculateSolution(const double tRun, bool shouldWriteData =false);                                ///< Running simulations 

private:
    void sample_hips_eddy(double &dt, int &iLevel);                           ///< Sample hips eddy with specified time step and level                                                
    void selectAndSwapTwoSubtrees(const int iLevel, int &iTree);              ///< Select and swap two subtrees in the level tree
    void advanceHips(const int iLevel, const int iTree);                      ///< Advancing simulations to do mixing and reaction
    void reactParcels_LevelTree(const int iLevel, const int iTree);           ///< Reacting parcels involved in micro-mixing
    void mixAcrossLevelTree(int kVar, const int iMixLevel, const int iTree);  ///< Mixing paecels involved in micr0-mixing.
    void forceProfile();
    void writeData(const int ifile, const double outputTime);                 ///< Writing the results for a user-defined number of eddies in the data folder.
};
























