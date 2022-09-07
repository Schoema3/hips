/**
 * @file domain.cc
 * Source file for class \ref domain
 */

#include "domain.h"
#include "dv.h"
#include "dv_rho.h"
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


  //----------------------

    io   = new inputoutput(this, caseName, nShiftFileNumbers);
     
    pram = new param(io);
    ngrd    = pram->ngrd0;
 
    rand = new randomGenerator(pram->seed + nShiftFileNumbers);
    LrandSet = true;  


    solv= new solver(this, pram);

    LsolvSet = true;





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
domain::~domain(){
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





