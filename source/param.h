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
        string  chemMechFile;    ///<  name of chemical mechanism file
        string  probType;        ///<  problem type: CHANNEL, CHANNEL_SCALAR, JETMIXL_RXN, COUETTE

<<<<<<< HEAD
        double  Z_param;         ///<  Viscous penalty parameter
        double  A_param;         ///<  Energy Distribution parameter alpha
        string  LES_type;        ///<  NONE, THIRDS, ELAPSEDTIME, FRACDOMAIN, INTEGRALSCALE
        double  diffCFL;         ///<  multiplies min diffusion timestep
        double  cvode_atol;      ///<  absolute tolerace atol for cvode
        double  cvode_rtol;      ///<  relative tolerace rtol for cvode
        double  x0virtual;       ///<  LES virtual origin

        bool    LdoDL;           ///<  flag to do the DL energy from the DL instability
        bool    Lrad;            ///<  radiation flag
        bool    Lbuoyant;        ///<  flag to turn on bouyancy (horizontal domain)
        bool    LPeEddy;         ///<  flag to turn on potential energy for eddies (vertical domain)
        bool    LplanarExpCent0; ///<  flag: for planar cases (C=1) set the expansion center at 0 for outflow cases (normally expand about the expansion center.
        double  g;               ///<  gravity (default -9.81)
        string  Lsolver;         ///<  EXPLICIT, SEMI-IMPLICIT, or STRANG
        bool    Lperiodic;       ///<  periodic if true
        bool    Lspatial;        ///<  spatial formulation if true
        bool    LTMA;            ///<  true for the triplet map TMA: 3 = vol segments; false for TMB: 3 equal length segments
        bool    LplanarTau;      ///<  true for computing cylindrical/spherical tau_eddy using a planar formulation. If accepted, a cylindrical eddy is implemented
        bool    Lignition;        ///<  true if starting with unreacted mixing profile to allow ignition

        string  bcType;          ///<  OUTFLOW, PERIODIC, WALL, WALL_OUT
        int     cCoord;          ///<  1 = planar, 2 = cylindrical, 3 = spherical
        double  xDomainCenter;   ///<  position of the center of the domain

        double  gDens;           ///<  grid density for mesher
        double  dxmin;           ///<  min grid spacing: = dxmin / domain length
        double  dxmax;           ///<  max grid spacing = dxmax / domain length

        double  Pmax;            ///<  maximum eddy acceptance probability
        double  Pav;             ///<  Average acceptance probability
        double  dtfac;           ///<  maximum factor to increase dtSmean
        int     nDtSmeanWait;    ///<  number of eddy samples before increase dtSmean
        int     eddyMinCells;    ///<  eddy must overlap at least this many cells
        double  DAtimeFac;       ///<  time until catch-up adaption is DAtimeFac * dtCUmax
        double  tdfac;           ///<  factor between dtCUmax and dtCFL for temporal flows; DEFAULT = 1.0
        int     sLastDA;         ///<  size of the lastDA vector for timing adaptmesh after diff
        double  Lp;              ///<  Most probable eddy size frac of domainLength
        double  Lmax;            ///<  Max eddy size frac of domainLength
        double  Lmin;            ///<  Min eddy size frac of domainLength
=======
        double  C_param;         ///<  Eddy frequency parameter
        double  diffCFL;         ///<  multiplies min diffusion timestep
        double  cvode_atol;      ///<  absolute tolerace atol for cvode
        double  cvode_rtol;      ///<  relative tolerace rtol for cvode
        string  Lsolver;         ///<  EXPLICIT, SEMI-IMPLICIT, or STRANG 
>>>>>>> Edit_hips

        int     modDump;         ///<  accepted eddies before output file
        bool    Lrestart;        ///<  true to restart from file, else false
        string  rstType;         ///<  "single" or "multiple"

        double  trst;            ///<  restart time (from restart file), default is 0.0;
        double  pres;            ///<  initial pressure (Pa)
        bool    Ltecplot;        ///<  set TRUE for tecplot friendly output


       
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

<<<<<<< HEAD

=======
>>>>>>> Edit_hips
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


