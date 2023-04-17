/**
 * @file domain.cc
 * Source file for class \ref domain
 */

#include "domain.h"
#include "dv.h"
#include "dv_rho.h"
<<<<<<< HEAD
#include "dv_dvisc.h"
#include "solver.h"
=======
>>>>>>> Edit_hips
#include "domaincase_hips.h"
#include "domaincase_hips_comb.h"
#include "domaincase_hips_simpleRxn.h"
#include <cmath>
#include <iomanip>
#include "processor.h"
#include "domaincase.h"
#include "inputoutput.h"
#include "param.h"
#include "streams.h"
#include "micromixer.h"
#include "solver.h"
#include "randomGenerator.h"


extern processor proc;
/////////////////////////////////////////////////////////////////////
/** Constructor
 */

domain::domain(domain *p_domn,   int   nShiftFileNumbers, string caseName){

    domn = p_domn;
    LdomcSet = false;
    LmimxSet = false;
    LsolvSet = false;
    LrandSet = false;
 
    LioSet   = false;
    LpramSet = false;


<<<<<<< HEAD
/////////////////////////////////////////////////////////////////////
/** Initializer
 */

void domain::init(inputoutput     *p_io,
                  streams         *p_strm,
                  IdealGasPhase   *p_gas,
                  Transport       *p_tran,
                  micromixer      *p_mimx,
                  domain          *p_eddl,
                  solver          *p_solv,
                  randomGenerator *p_rand,
                  bool             LisEddyDomain) {

    //----------------------

    io     = p_io;
    mesher = p_mesher;
    gas    = p_gas;
    tran   = p_tran;
    strm   = p_strm;
    mimx   = p_mimx;
    solv   = p_solv;
    rand   = p_rand;
    prb    = p_prb;

    //----------------------
=======
  //----------------------
>>>>>>> Edit_hips

    io   = new inputoutput(this, caseName, nShiftFileNumbers);
     
    pram = new param(io);
    ngrd    = pram->ngrd0;
 
    rand = new randomGenerator(pram->seed + nShiftFileNumbers);
    LrandSet = true;  


    solv= new solver(this, pram);

    LsolvSet = true;


<<<<<<< HEAD
    //----------------------
    io->init(this);
    pram->init(this);
    prb->init(this);
    solv->init(this);
    // mesher is init below in caseinit for phi
    // strm is init below in caseinit  (domc), (if needed)
    // mimx is init below since it needs v[] set for cvode
=======


>>>>>>> Edit_hips

    //---------------------- Continue setting up the case using the case_somecase class.
    // Adds to the above variable list, and initializes solution for the run

<<<<<<< HEAD
     else if(pram->probType == "HIPS")
         domc = new domaincase_hips(); // hips

     else if(pram->probType == "HIPS_COMB")
=======

    if(pram->probType == "HIPS")
        domc = new domaincase_hips(); // hips

      else if(pram->probType == "HIPS_COMB")
>>>>>>> Edit_hips
         domc = new domaincase_hips_comb(); // hips_comb

     else if(pram->probType == "HIPS_SIMPLERXN")
         domc = new domaincase_hips_simpleRxn(); // hips_simpleRxn

     else {
         cout << endl << "ERROR, probType UNKNOWN" << endl;
         exit(0);
     }

   LdomcSet   = true;
    domc->init(this);

    //----------------------

    for(int k=0; k<v.size(); k++)
        varMap[v.at(k)->var_name] = v.at(k);
    nTrans = 0;
    for(int k=0; k<v.size(); k++)
        if(v.at(k)->L_transported)
            nTrans++;

    //----------------------
    mimx = new micromixer();
    LmimxSet = true;
    mimx->init(this);
    //----------------------

    if(pram->Lrestart) {
        io->loadVarsFromRestartFile();
        io->set_iNextDumpTime(pram->trst);
    }

}

/////////////////////////////////////////////////////////////////////
/** Destructor
 */
<<<<<<< HEAD

double domain::Ldomain() {
     return posf->d.at(ngrd) - posf->d.at(0);
}

/////////////////////////////////////////////////////////////////////
/** Initialize data members of the eddy domain.
 *  Note, none of the other members of this domain should be used (like random).
 *  Note, all variables here should have corresponding variables (by var_name) in the
 *     main domn. This is needed for using the eddl object.
 */

void domain::initEddyDomain() {

    v.push_back(new dv_pos(  this, "pos",   false, true));
    v.push_back(new dv_posf( this, "posf",  false, true));
    v.push_back(new dv_rho(  this, "rho",   false, false));
    v.push_back(new dv_dvisc(this, "dvisc", false, false));
    if(domn->pram->LdoDL)
       v.push_back(new dv_aDL(this, "aDL",   false, false));

    int k = 0;
    rho   = v.at(k++);
    dvisc = v.at(k++);
    if(domn->pram->LdoDL)
        aDL   = v.at(k++);

}

/////////////////////////////////////////////////////////////////////
/** Set the domain from a region of the domn.  Normally called by eddy domain.
 *  @param i1 \input index of starting cell of domn to build from
 *  @param i2 \input index of ending cell of domn to build from
 *  If i2 < i1, we have a periodic region (wrap around the domain).
 *     This only happens in planar cases, not cylindrical or sphericial.
 *  nonwrap: |   | * | * | * | * | * | * |   |   |
 *                i1                  i2
 *  new domain consists of *'d cells
 *
 *  Wrap: | 4 | 5 |   |   |   |   | 1 | 2 | 3 |
 *             i2                  i1
 *  New domain consists of #'d cells: 1 2 3 4 5}
 */

void domain::setDomainFromRegion(const int i1, const int i2) {

    ngrd  = i2-i1+1;
    ngrdf = ngrd+1;

=======
domain::~domain(){
>>>>>>> Edit_hips
    for(int k=0; k<v.size(); k++)
        delete v.at(k);
    if(LdomcSet) delete domc; 
    if(LmimxSet) delete mimx;
    if(LsolvSet) delete solv;
    if(LrandSet) delete rand;

    if(LioSet)   delete io;
    if(LpramSet) delete pram;





}

void domain::hips_advance(){

    solv->calculateSolution();

}





