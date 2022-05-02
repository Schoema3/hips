/**
 * @file domain.cc
 * Source file for class \ref domain
 */

#include "domain.h"
#include "dv.h"
#include "dv_rho.h"
#include "dv_dvisc.h"
#include "domaincase_hips.h"
#include "domaincase_hips_comb.h"
#include "domaincase_hips_simpleRxn.h"
#include <cmath>
#include <iomanip>
#include "processor.h"

extern processor proc;
/////////////////////////////////////////////////////////////////////
/** Constructor
 */

domain::domain(domain *p_domn, param *p_pram) {

    domn = p_domn;
    pram = p_pram;
    domc = 0;               // initialize for destruction of eddy domains

 }

/////////////////////////////////////////////////////////////////////
/** Initializer
 */

void domain::init(inputoutput     *p_io,
                  streams         *p_strm,
                  IdealGasPhase   *p_gas,
                  Transport       *p_tran,
                  micromixer      *p_mimx,
                  solver          *p_solv,
                  randomGenerator *p_rand) {

    //----------------------

    io     = p_io; 
    gas    = p_gas;
    tran   = p_tran;
    strm   = p_strm;
    mimx   = p_mimx;
    solv   = p_solv;
    rand   = p_rand;
 

    //----------------------

    ngrd    = pram->ngrd0;
    ngrdf   = ngrd + 1;

    //----------------------

    io->init(this);
    pram->init(this);
    solv->init(this);



    //---------------------- Continue setting up the case using the case_somecase class.
    // Adds to the above variable list, and initializes solution for the run


    if(pram->probType == "HIPS")
         domc = new domaincase_hips(); // hips

     else if(pram->probType == "HIPS_COMB")
         domc = new domaincase_hips_comb(); // hips_comb

     else if(pram->probType == "HIPS_SIMPLERXN")
         domc = new domaincase_hips_simpleRxn(); // hips_simpleRxn

     else {
         cout << endl << "ERROR, probType UNKNOWN" << endl;
         exit(0);
     }

    domc->init(this);

    //----------------------

    for(int k=0; k<v.size(); k++)
        varMap[v.at(k)->var_name] = v.at(k);

    nTrans = 0;
    for(int k=0; k<v.size(); k++)
        if(v.at(k)->L_transported)
            nTrans++;

    //----------------------

    mimx->init(this);

    //----------------------

    if(pram->Lrestart) {
        io->loadVarsFromRestartFile();
        io->set_iNextDumpTime(pram->trst);
    }

}

/////////////////////////////////////////////////////////////////////
/** Compute size of domain based on faces.
 */

double domain::Ldomain() {
    return posf->d.at(ngrd) - posf->d.at(0);
}




