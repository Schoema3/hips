/**
 * @file domain.h
 * Header file for class \ref domain
 */

#pragma once

#include "dv.h"
#include "domaincase.h"
#include "inputoutput.h"
#include "param.h"
#include "streams.h"
#include "micromixer.h"
#include "solver.h"
#include "randomGenerator.h"

#ifdef DOCANTERA
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/transport.h"
#else
#include "cantera_shell_functions.h"
#endif

#include <vector>
#include <string>
#include <map>

using namespace std;
using namespace Cantera;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing domain object
 *
 *  @author David O. Lignell
 */

class domain {

    public:

    //////////////////// DATA MEMBERS //////////////////////

        domain                  *domn;     ///< (for one domain to point to another (eddl))
        int                     ngrd;      ///< number of parcels
        vector<dv*>             v;         ///< All domain variables are stored in here.
        dv*                     rho;
        dv*                     enth;
        dv*                     temp;
        dv*                     mixf;
        dv*                     chi;
        vector<dv*>::iterator   ysp;       ///< access as: ysp=v.begin(), (*ysp)->d[i] or (*(ysp+k))->d[i], or ysp[k]->d[i].
        vector<dv*>::iterator   svar;      ///< iterator for increment to go through moments (*(ysp+k))->d[i];)
        
        map<string,dv*>         varMap;

        IdealGasPhase           *gas;        ///< pointer to cantera thermochemistry object (reaction rates, Cp, etc.)
        Transport               *tran;       ///< pointer to cantera transport object (viscosity, diffusivity, etc.)
        streams                 *strm;       ///< pointer to gas stream properties
        inputoutput             *io;         ///< pointer to input/output object
        param                   *pram;       ///< pointer to the parameters object
        micromixer              *mimx;       ///< pointer to micromixer for diffusion, reaction, domain evolution.
        solver                  *solv;       ///< pointer to solver object
       
        randomGenerator         *rand;

        int                     nTrans;      ///< number of transported variables on the domain.

        domaincase              *domc;       ///< domaincase class: set specific vars...






        bool                    LdomcSet;    ///< flag indicating new domainCase --> allow deletion
       // bool                    LstrmSet;    ///< flag indicating new streams    --> allow deletion
        bool                    LmimxSet;    ///< flag indicating new micromixer --> allow deletion
        bool                    LsolvSet;    ///< flag indicating new solver     --> allow deletion
        bool                    LrandSet;    ///< flag indicating new randomGen  --> allow deletion
        bool                    LioSet;      ///< flag indicating new inputoutput--> allow deletion
        bool                    LpramSet;    ///< flag indicating new param      --> allow deletion

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

     
              void init(
                int   nShiftFileNumbers,
                  string caseName); 
        domain(domain *p_domn);
        virtual ~domain();
        //{
        //    for(int k=0; k<v.size(); k++)
        //        delete v.at(k);
        //    delete domc;     
        //}

};


////////////////////////////////////////////////////////////////////////////////


