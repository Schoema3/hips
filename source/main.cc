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

    int nShiftFileNumbers = 0;          // we should increment the seed if we are starting MPI multiple times
    stringstream ss1; ss1.clear(); ss1 << argv[2]; ss1 >> nShiftFileNumbers;
   
    domain domn(NULL,nShiftFileNumbers, caseName);
  
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
