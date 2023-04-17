/**
 * @file main.cc
 * @brief main driver for HiPS 
 */

#include "domain.h"
#include "yaml-cpp/yaml.h"
#include <iostream>
#include <string>
#include <ctime>
#include <sstream>

using namespace std;


processor proc;
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////

int main(int argc, char*argv[]) {

 if(argc<3) {
        cout << endl << "ERROR: code needs caseName and shift arguments" << endl;
        return 1;
    }
    string caseName= argv[1];            // example: temporalJet (../input/temporalJet, without the ../input/)

<<<<<<< HEAD
    int nShiftFileNumbers = 0;
    stringstream ss1;
    string       s1;
    ss1.clear(); ss1 << argv[2];
    ss1 >> nShiftFileNumbers;

    inputoutput   io(caseName, nShiftFileNumbers);
    param         pram(&io);
    probes        prb;
    streams       strm;
    IdealGasPhase gas("../input/gas_mechanisms/"+pram.chemMechFile);
    Transport    *tran = newTransportMgr("Mix", &gas);
    solver       *solv;
    micromixer   *mimx;

    
    else if(pram.LisHips) {
        solv = new solver_hips();
        mimx = new micromixer_hips();
    }
    else {
        solv = new solver();
        mimx = new micromixer();
    }

    domain domn(NULL,  &pram);
    domain eddl(&domn, &pram);

    // we should increment the seed if we are starting MPI multiple times
    if ( pram.seed >= 0 ) pram.seed += nShiftFileNumbers;
    randomGenerator rand(pram.seed);

    domn.init(&io, &strm, &gas, tran, mimx, solv, &rand, &prb);
    eddl.init(NULL, NULL, NULL,  NULL, NULL, NULL, NULL, NULL,  NULL, NULL,  NULL, true);
    //
=======
    int nShiftFileNumbers = 0;          // we should increment the seed if we are starting MPI multiple times
    stringstream ss1; ss1.clear(); ss1 << argv[2]; ss1 >> nShiftFileNumbers;
   
    domain domn(NULL,nShiftFileNumbers, caseName);
  
>>>>>>> Edit_hips
    //-------------------

  time_t mytimeStart, mytimeEnd;
   mytimeStart = time(0);
  *domn.io->ostrm << endl << "#################################################################";
   *domn.io->ostrm << endl << "#  Start Time = " << ctime(&mytimeStart);
   *domn.io->ostrm << endl << "#################################################################";

   cout << endl << "made it c" << endl; //doldb

    //-------------------

   cout << endl << "made it d" << endl; //doldb

   domn.hips_advance();
   cout << endl << "made it e" << endl; //doldb

  mytimeEnd = time(0);
  *domn.io->ostrm << endl << "#################################################################";
  *domn.io->ostrm << endl << "#  Start Time = " << ctime(&mytimeStart);
  *domn.io->ostrm << endl << "#  End Time   = " << ctime(&mytimeEnd);
  *domn.io->ostrm << endl << "#################################################################";
  *domn.io->ostrm << endl;

    return 0;


}
