using namespace std;
/**
 * @file solver.h
 * Header file for class solver
 */

#pragma once

#include <vector>


class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing solver object
 *
 *  @author David O. Lignell
 */

class solver {


    //////////////////// DATA MEMBERS //////////////////////

    public:

        domain         *domn;          ///< pointer to domain object

        double         time;           ///< odt time (during sampling)
        double         t0;             ///< time of last eddy event; diffusion left off here.
        double         dtSmean;        ///< initial mean eddy sample time
        double         dtCUmax;        ///< max time before catch up diff/eddy

        bool           LeddyAccepted;  ///< flag for accepted eddy
        int            iEtrials;       ///< number of eddy trials

        double         PaSum;          ///< sum of Pa of eddies
        int            nPaSum;         ///< number going into PaSum
        int            neddies;        ///< number of eddies accepted
        double         PaSumC;         ///< sum of Pa of eddies
        int            nPaSumC;        ///< number going into PaSum

        //---------- for hips interface (inherited)

        double         tMix;           ///< parcel mixing timescale
        vector<int>    pLoc;           ///< parcel index array for fast implementation of swaps
        int iSc;                       ///< for Sc < 1, tree level index at/below which all parcels are mixed

        int Nm1;                       ///< nLevels - 1
        int Nm2;                       ///< nLevels - 2
        int Nm3;                       ///< nLevels - 3
        int iTree;                     ///< base node (at iLevel) for the swap


    //////////////////// MEMBER FUNCTIONS /////////////////

    public:

 virtual void init(domain *p_domn);
        virtual ~solver();    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////
        virtual void calculateSolution();

    private:

        bool   sampleEddyAndImplementIfAccepted();
        bool   sampleAndImplementLEMeddy();
        void   computeDtSmean();
        void   computeDtCUmax();
        double sampleDt();
        void   diffusionCatchUpIfNeeded(bool Ldoit=false);
        void   raiseDtSmean();
        void   lowerDtSmean();
        bool   testLES_elapsedTime(const double time, const double tauEddy);
        bool   testLES_fracDomain( const double eSize);
        bool   testLES_integralLength(const double time, const double eSize);
        bool   testLES_thirds();

//-----------------------------------------------------------------------------------
        vector<double> levelRates;     ///< list of eddy event rates at each level
        vector<pair<double,int> > eTL; ///< list of eddy times and levels
        int iEta;                      ///< Kolmogorov level (needed for variable Sc scalars)
        double eddyRate_total;         ///< total rate of all eddies 0 through nLevels-3
        double eddyRate_inertial;      ///< tota    public:

        void setEddyEventTimes();
        void selectAndSwapTwoSubtrees(const int iLevel, int &Qstart, int &Rstart, int &nPswap);
        void reset_rates_for_Sc(const vector<double> &levelTaus);
        void sample_hips_eddy(double &dt, double &iLevel);
//---------------------------------------------------------------------------------------
    protected:
};


////////////////////////////////////////////////////////////////////////////////


