/**
 * @file solver.h
 * Header file for class \ref solver
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
 };


////////////////////////////////////////////////////////////////////////////////


