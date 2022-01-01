/**
 * @file solver.h
 * Header file for class solver
 */

#pragma once"

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
        double         dtCUmax;        ///< max time before catch up diff/eddy

        int            nPaSum;         ///< number going into PaSum
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

        virtual void calculateSolution();

    private:

        void   computeDtSmean();
        void   computeDtCUmax();
        double sampleDt();
        void   diffusionCatchUpIfNeeded(bool Ldoit=false);
        void   raiseDtSmean();
        void   lowerDtSmean();
       

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        virtual void init(domain *p_domn);
        virtual ~solver();

};


////////////////////////////////////////////////////////////////////////////////


