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
       

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        solver(){ ed3=0; eddl3=0; }
        virtual void init(domain *p_domn);
        virtual ~solver();

};


////////////////////////////////////////////////////////////////////////////////


