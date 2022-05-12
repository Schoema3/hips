/**
 * @file param.h
 * Header file for class \ref param
 */

#pragma once

#include "inputoutput.h"
#include <string>
#include <cstdlib>

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing inputoutput object
 *
 *  @author David O. Lignell
 */

class param {

    public:

    //////////////////// DATA MEMBERS //////////////////////

        domain  *domn;           ///< pointer to domain object
        inputoutput *io;         ///< pointer to io object (has the input file)

        int     seed;            ///<  random number generator seed (negative to randomize it)
        double  tEnd;            ///<  ending time of realization
        double  domainLength;    ///<  length of domain (m)
        int     ngrd0;           ///<  initial grid points
        double  rho0;            ///<  initial uniform density (kg/m^3)
        double  kvisc0;          ///<  initial uniform kinematic viscosity (m^2/s)
        double  pres;            ///<  initial pressure (Pa)
        string  chemMechFile;    ///<  name of chemical mechanism file
        string  probType;        ///<  problem type: CHANNEL, CHANNEL_SCALAR, JETMIXL_RXN, COUETTE

        double  C_param;         ///<  Eddy frequency parameter
        double  diffCFL;         ///<  multiplies min diffusion timestep
        double  cvode_atol;      ///<  absolute tolerace atol for cvode
        double  cvode_rtol;      ///<  relative tolerace rtol for cvode

        bool    LdoDL;           ///<  flag to do the DL energy from the DL instability
        double  g;               ///<  gravity (default -9.81)
        string  Lsolver;         ///<  EXPLICIT, SEMI-IMPLICIT, or STRANG
        bool    Lperiodic;       ///<  periodic if true
        bool    LisFlmlt;        ///<  true if solving an unsteady flamelet
        bool    LisFlmltX;       ///<  true if solving an unsteady flamelet in the physical domain (Pierce 2004)
        int     modDump;         ///<  accepted eddies before output file
        bool    Ltecplot;        ///<  set TRUE for tecplot friendly output
        bool    Lrestart;        ///<  true to restart from file, else false
        string  rstType;         ///<  "single" or "multiple"
        double  trst;            ///<  restart time (from restart file), default is 0.0;

        //----------------- Soot variables

        bool           Lsoot;               ///< true for soot, false for no soot
        int            nsvar;               ///< number of soot variables transported (# soot moments)
        string         PSD_method;          ///< method name for soot PSD: MONO, QMOM, MOMIC

        //----------------- HIPS quantities

        bool    LisHips;         ///<  true if solving hips
        int     nLevels;         ///< number of levels in the tree: 0, 1, 2, ... N-1
        double  Afac;            ///< level lengthscale reduction factor (0.5)
        double  L0;              ///< tree lengthscale
        double  tau0;            ///< integral timescale
        double  fmix;            ///< timescale factor for micromixing
        double  LScHips;         ///< hips schmidt number
        bool    LsimpleMix;      ///< true for simple instantaneous mixing of parcel pairs
        int     forceHips;       ///< forcing function for statistically stationary: -1 = none, 1 = source term, 2 = direct profile

    //////////////////// MEMBER FUNCTIONS /////////////////

    private:

        template <class T>
        T errMsg(const string param) {
            *io->ostrm << endl << "ERROR: missing parameter: " + param << endl;
            exit(0);
            T dummy = static_cast<T> (0);
            return dummy;
        }

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        param(inputoutput *p_io);
        void init(domain *p_domn);
        ~param(){}

};


////////////////////////////////////////////////////////////////////////////////


