
#include "ame.h"
#include "yaml-cpp/yaml.h"

#ifdef DOCANTERA
#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "cantera/transport.h"
#include "cantera/kinetics.h"
#else
#include "cantera_shell_functions.h"
#endif

#include <iostream>
#include <string>
#include <ctime>
#include <sstream>

using namespace std;
using namespace Cantera;

//////////////////////////////////////////////////////////////

int driver_hips(int argc, char*argv[]) {

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

    auto sol  = newSolution("../input/gas_mechanisms/"+chemMechFile, "", "mixture-averaged");
    auto gas  = sol->thermo();
    auto tran = sol->transport();
    auto kin  = sol->kinetics();

    ame AME(gas, tran, kin, caseName);
    //-------------------

    for(auto d : AME.ogrd) //we need this
        d->solv->calculateSolution();

    //-------------------

  
    return 0;
}
