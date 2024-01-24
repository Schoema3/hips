#pragma once

#include "cantera/base/Solution.h"
#include "cantera/thermo.h"

#include <iostream>
#include <vector>
#include <string>
#include <vector>
#include "randomGenerator.h"
#include "batchReactor.h"

/////////////////////////////////////////////////////////////////////////////


class hips {


////////////////////////////// DATA MEMBERS /////////////////////////////

    public:

        int                                    nparcels; ///< number of parcels
        std::vector<std::vector<double > * >   varData;  ///< vector of pointers to vector
        std::shared_ptr<Cantera::ThermoPhase>  gas;      ///< cantera thermo object
        std::unique_ptr<batchReactor>          bRxr;     ///< chemistry integrator (constant pressure)

                        

    private:

        int                 nLevels;                     ///< number of tree levels
        int                 nLevels_;                    ///< number of tree levels?
        double              domainLength;                ///< length of domain (m)
        double              tau0;                        ///< integral timescale
        double              C_param;                     ///< Eddy frequency parameter
        int                 forceTurb;                   ///< forcing function for statistically stationary: -1 = none, 1 = source term, 2 = dir
        int                 nVar;                        ///< number of parcel variables (e.g., h, ysp)
        int                 nsp;                         ///< number of species

        double              time;          
        int                 Nm1;                         ///< nLevels - 1 
        int                 Nm2;                         ///< nLevels - 2                     
        int                 Nm3;                         ///< nLevels - 3 
        int                 iEta;                        ///< Kolmogorov level (needed for variable Sc scalars)
        double              eddyRate_total;              ///< total rate of all eddies 0 through nLevels-3
        double              eddyRate_inertial;           ///< ?????
        bool                LScHips;                     ///< hips schmidt number
        std::vector<double> i_batchelor;                 ///< Batchelor level for variable Sc scalars; NOTE: double, as in, between levels
        std::vector<int>    i_plus;                      ///< ceil(i_batchelor)
        std::vector<double> parcelTimes;                 ///< current times corresponding to the parcel states
      
        std::vector<double> levelRates;                  ///< list of eddy event rates at each level
        std::vector<int>    pLoc;                        ///< parcel index array for fast implementation of swaps
        std::vector<double> ScHips;                      ///< Vector containing Schmidt numbers related to each variable


        bool                performReaction;             ///< Flag indicating whether chemical reactions are performed in the simulation 
        bool                LrandSet;                    ///< flag indicating new randomGen  --> allow deletion
        randomGenerator     rand;

        const double        Afac = 0.5;                 ///< level lengthscale reduction factor (0.5) 
        std::vector<std::string>  varName;              ///< Vector containing the names of parcel variables

  /////////////////////////////// MEMBER FUNCTIONS//////////////////////////////////////////

    public:
       
        /**
        * @brief Constructor for the hips class.
        * Initializes a hips object with specified parameters.
        */            

        hips(int             nLevels_,
             double          domainLength_,
             double          tau0_,
             double          C_param_,
             int             forceTurb_,
             int             nVar_,
             std::vector<double> &ScHips_,
             std::shared_ptr<Cantera::Solution> cantSol,
             bool            performReaction,
             int             seed = 10);

        hips(double          Re);                              ///< Constructor for the hips class with a Reynolds number parameter
        virtual ~hips() { for( auto& data : varData)
                              delete data;
        }                                                                                         ///< Destructor for the hips class

        void set_varData(std::vector<double> &v, std::vector<double> &w, const std::string &varN, int i);          ///< passing all variables to vector of pointer 
        std::vector<double> projection(std::vector<double> &vcfd, std::vector<double> &weight);                    ///< Perform vector projection operation
        std::vector<double> setGridHips(int N);                                                                    ///< Set Hips grid with a specified number of  \
                                                                                                                   ///< grid points equal to number of parcels   
        std::vector<double> setGridCfd(std::vector<double> &w);                                                    ///< Set CFD grid using provided weight vector
        std::vector<std::vector<double>> get_varData();                                                //Retrieves modified data from the HiPS library and stores it in the provided vector. 
       std::vector<double>  projection_back(std::vector<double> &vb);
        void calculateSolution(const double tRun);                                                                 ///< Running simulations 


    private:

        void sample_hips_eddy(double &dt, int &iLevel);                                                            ///< Sample hips eddy with specified time step and level                                                
        void selectAndSwapTwoSubtrees(const int iLevel, int &iTree);                                               ///< Select and swap two subtrees in the level tree
        void advanceHips(const int iLevel, const int iTree);
        void reactParcels_LevelTree(const int iLevel, const int iTree);
        void mixAcrossLevelTree( int kVar, const int iMixLevel, const int iTree);                       // const int iTree=-1);
        void forceProfile();
        void writeData(const int ifile, const double outputTime);
        
            
};


