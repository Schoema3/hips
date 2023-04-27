#include "hips.h"

#include "cantera/thermo.h"
#include "cantera/base/Solution.h"
#include "cantera/numerics/Integrator.h"

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

///////////////////////////////////////////////////////////////////////////////

int main() {
  
    //---------- setting  HiPS tree and constructor

    int            nLevels      = 9;           // # HiPS tree levels
    double         domainLength = 1.0;         // HiPS tree lengthscale
    double         tau0         = 1.0;         // HiPS tree timescale
    double         C_param      = 0.5;         // HiPS eddy rate multiplier
    double         tRun         = 400.0;       // tree evolution time
    int            forceTurb    = 0;           // Force turbulent profile
    int            nVar         = 7;           // 1 + nsp: T, y 
    vector<double> ScHips(7,1);

    //---------- setting gas solution

    //auto   cantSol = Cantera::newSolution("2S_CH4_BFER.yaml", "2S_CH4_BFER", "None");
    auto   cantSol = Cantera::newSolution("gri30.yaml");
    auto   gas = cantSol->thermo();
    size_t nsp = gas->nSpecies();

    hips HiPS(nLevels,
              domainLength,
              tau0, 
              C_param,
              forceTurb,
              nVar,
              ScHips,
              cantSol);

    int nparcels = HiPS.nparcels;

    //------------- declaring state vectors

    vector<vector <double > > ysp(6, vector<double> (nparcels, 0));
    vector<double>            h(nparcels); 

    vector<double>            y0(nsp);
    double                    T0;
    double                    h0;

    vector<double>            y1(nsp);
    double                    T1;
    double                    h1;

    //-------------- set state vectors in each parcel

    T0 = 700.0;
    gas->setState_TPX(T0, Cantera::OneAtm, "CH4:1");
    gas->getMassFractions(&y0[0]);
    h0 = gas->enthalpy_mass();

    T1 = 1200;
    gas->setState_TPX(T1, Cantera::OneAtm, "O2:2, N2:7.52");
    gas->getMassFractions(&y1[0]);
    h1 = gas->enthalpy_mass();


    for (int i=0; i<nparcels/2; i++) {       // left parcels are stream 0 (air)
        h[i] = h0;
        for (int k=0; k<nsp; k++)
            ysp[k][i] = y0[k];
    }

    for (int i=nparcels/2; i<nparcels; i++) { // right parcels are stream 1 (fuel)
        h[i] = h1;
        for (int k=0; k<nsp; k++)
            ysp[k][i] = y1[k];
    }

    //----------------- set array of hips variables from state variables

    HiPS.set_varData(h, 0);
    for (int k=0; k<ysp.size(); k++) 
        HiPS.set_varData(ysp[k],k+1);

    //----------------- advancing HiPS to do rxn and mixing

    HiPS.calculateSolution(tRun);

    return 0;

}  
