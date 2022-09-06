/**
 * @file main.cc
 * @brief main driver for HiPS 
 */

#include "domain.h"
#include "yaml-cpp/yaml.h"
//#include "streams.h"
//#include "param.h"
//#include "micromixer.h"
//#include "processor.h"
//#include "solver.h"
//#include "randomGenerator.h"
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
//using namespace Cantera;


processor proc;
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////

int main(int argc, char*argv[]) {

 if(argc<3) {
        cout << endl << "ERROR: code needs caseName and shift arguments" << endl;
        return 1;
    }
    string caseName= argv[1];            // example: temporalJet (../input/temporalJet, without the ../input/)

    int nShiftFileNumbers = 0;          // we should increment the seed if we are starting MPI multiple times
    stringstream ss1; ss1.clear(); ss1 << argv[2]; ss1 >> nShiftFileNumbers;

    string inputFileDir  = "../data/"+caseName+"/input/";
    YAML::Node inputFile = YAML::LoadFile(inputFileDir+"input.yaml");
    YAML::Node params    = inputFile["params"];
    string chemMechFile  = params["chemMechFile"] ? params["chemMechFile"].as<string>() : "";

   // auto sol  = newSolution("../input/gas_mechanisms/"+chemMechFile, "", "mixture-averaged");
   // auto gas  = sol->thermo();
   // auto kin  = sol->kinetics();
   // auto tran = sol->transport();

   cout << endl << "made it a" << endl; //doldb

    domain domn(NULL);
    //domn domn(NULL, gas, tran, nShiftFileNumbers, caseName);
   cout << endl << "made it b" << endl; //doldb
   domn.init(nShiftFileNumbers, caseName);
   cout << endl << "made it bb" << endl; //doldb




    //-------------------

    time_t mytimeStart, mytimeEnd;
   mytimeStart = time(0);
  *domn.io->ostrm << endl << "#################################################################";
   *domn.io->ostrm << endl << "#  Start Time = " << ctime(&mytimeStart);
   *domn.io->ostrm << endl << "#################################################################";

   cout << endl << "made it c" << endl; //doldb

    //-------------------

   cout << endl << "made it d" << endl; //doldb

    domn.solv->calculateSolution();

   cout << endl << "made it e" << endl; //doldb

  mytimeEnd = time(0);
  *domn.io->ostrm << endl << "#################################################################";
  *domn.io->ostrm << endl << "#  Start Time = " << ctime(&mytimeStart);
  *domn.io->ostrm << endl << "#  End Time   = " << ctime(&mytimeEnd);
  *domn.io->ostrm << endl << "#################################################################";
  *domn.io->ostrm << endl;

    return 0;


}
