/**
 * @file main.cc
 * @brief main driver for HiPS 
 */

#include "domain.h"
#include "streams.h"
#include "param.h"
#include "micromixer.h"
#include "processor.h"
#include "solver.h"
#include "randomGenerator.h"
#ifdef DOCANTERA
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/transport.h"
#else
#include "cantera_shell_functions.h"
#endif
#include <iostream>
#include <string>
#include <ctime>
#include <sstream>

using namespace std;
using namespace Cantera;


processor proc;
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////

int main(int argc, char*argv[]) {


    if(argc<3) {
        cout << endl << "ERROR: code needs caseName and shift arguments" << endl;
        return 1;
    }
    string caseName= argv[1];            // example: temporalJet (../input/temporalJet, without the ../input/)

    int nShiftFileNumbers = 0;
    stringstream ss1;
    string       s1;
    ss1.clear(); ss1 << argv[2];
    ss1 >> nShiftFileNumbers;


    inputoutput   io(caseName, nShiftFileNumbers);
    param         pram(&io);
    streams       strm;
    IdealGasPhase gas("../input/gas_mechanisms/"+pram.chemMechFile);
    Transport    *tran = newTransportMgr("Mix", &gas);
   
    solver       *solv;
    micromixer   *mimx;


    solv = new solver();
    mimx = new micromixer();
    domain domn(NULL,  &pram);

  

    // we should increment the seed if we are starting MPI multiple times
    if ( pram.seed >= 0 ) pram.seed += nShiftFileNumbers;
    randomGenerator rand(pram.seed);

    domn.init(&io, &strm, &gas, tran, mimx, solv, &rand);
        //
    //-------------------

    time_t mytimeStart, mytimeEnd;
    mytimeStart = time(0);
    *io.ostrm << endl << "#################################################################";
    *io.ostrm << endl << "#  Start Time = " << ctime(&mytimeStart);
    *io.ostrm << endl << "#################################################################";


    //-------------------

    domn.solv->calculateSolution();

    mytimeEnd = time(0);
    *io.ostrm << endl << "#################################################################";
    *io.ostrm << endl << "#  Start Time = " << ctime(&mytimeStart);
    *io.ostrm << endl << "#  End Time   = " << ctime(&mytimeEnd);
    *io.ostrm << endl << "#################################################################";
    *io.ostrm << endl;

    return 0;


}
