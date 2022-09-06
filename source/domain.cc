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

domain::domain(domain *p_domn){
    domn = p_domn;
    LdomcSet = false;
    //LstrmSet = false;
    LmimxSet = false;
    LsolvSet = false;
    LrandSet = false;
 
    LioSet   = false;
    LpramSet = false;
}

/////////////////////////////////////////////////////////////////////
/** Destructor
 */
domain::~domain(){
    for(int k=0; k<v.size(); k++)
        delete v.at(k);
    if(LdomcSet) delete domc;
    //if(LstrmSet) delete strm;
    if(LmimxSet) delete mimx;
    if(LsolvSet) delete solv;
    if(LrandSet) delete rand;

    if(LioSet)   delete io;
    if(LpramSet) delete pram;
}


///////////////////////////////////////////////////////////////////:q//
/** Initializer
 *  @param \input nShiftFileNumbers increment the seed if starting MPI multiple times (standalone ODT realizations)
 */
void domain::init(int   nShiftFileNumbers,
                         string caseName
                        ) {

  


    //----------------------

   // gas    = p_gas;
   // tran   = p_tran;
    //strm   = p_strm;
 
    //----------------------

     io   = new inputoutput(this, caseName, nShiftFileNumbers);
     
    pram = new param(io);

   cout << endl << "here a" << endl; //doldb
    ngrd    = pram->ngrd0;
   // ngrdf   = ngrd + 1;

   cout << endl << "here b" << endl; //doldb
    //----------------------
 rand = new randomGenerator(pram->seed + nShiftFileNumbers);
    LrandSet = true;  

cout << endl << "here cm" << endl; //doldb

solv= new solver(this, pram);

 LsolvSet = true;

cout << endl << "here ebbb" << endl; //doldb




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
    cout << endl << "here ekkk" << endl; //doldb

   LdomcSet   = true;
    domc->init(this);

cout << endl << "here ejj" << endl; //doldb
    //----------------------

    for(int k=0; k<v.size(); k++)
        varMap[v.at(k)->var_name] = v.at(k);
    nTrans = 0;
    for(int k=0; k<v.size(); k++)
        if(v.at(k)->L_transported)
            nTrans++;

    //----------------------
cout << endl << "here m " << endl; //doldb
    mimx = new micromixer();
    LmimxSet = true;
    mimx->init(this);
 cout << endl << "here n" << endl; //doldb
    //----------------------

    if(pram->Lrestart) {
        io->loadVarsFromRestartFile();
        io->set_iNextDumpTime(pram->trst);
    }

}




