#include <catch2/catch_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "hips.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo.h"

#include <memory>
#include <cmath>
#include <format>
using Catch::Approx;

using namespace std;
using Catch::Matchers::WithinRel;

//////////////////////////////////////////////////////////////////////////

TEST_CASE( "Test HiPS library" ) {

    //int            nLevels = 6;          // Number of hierarchical levels; nLevels = 6 gives 32 parcels: 2^{nLevels-1}
    int            nLevels = 8;          // Number of hierarchical levels; nLevels = 6 gives 32 parcels: 2^{nLevels-1}
    double         domainLength = 0.01;  // Domain length scale
    //double         tau0 = 0.0018;        // Mixing timescale
    double         tau0 = 0.000005;        // Mixing timescale
    double         C_param = 0.5;        // Eddy rate multiplier
    double         tRun = 0.00012;   // Total simulation runtime
    //double         tRun = 0.000000018;   // Total simulation runtime
    bool           forceTurb = false;        // No forced turbulence
    vector<double> ScHips(54, 1);        // Schmidt number (unity for all species)

    auto   cantSol = Cantera::newSolution("gri30.yaml");
    auto   gas     = cantSol->thermo();  // cantera thermo object
    size_t nsp     = gas->nSpecies();    // # of gas species
    int    nVar    = nsp + 1;            // Number of variables (species + enthalpy)

    //////////////////////////////////////////////////////////////////////////

    SECTION("Test variable setup: projection and back projection") {

        vector<double> ScHips = {1.0}; // Schmidt numbers for low and high diffusivity
        int            nvars = 1;      // Number of scalar fields

        hips H(nLevels, domainLength, tau0, C_param, forceTurb, nvars, ScHips, false);

        vector<double> var = {0.0, 0.1, 0.2, 0.8, 0.9, 1.0};               // mixture fraction
        vector<double> w   = {0.5, 0.6, 0.7, 0.7, 0.6, 0.2};               // test: not a power of 2, not uniform, don't sum to 1
        vector<double> rho = {1.300, 0.474, 0.303, 0.171, 0.335, 1.057};   // equilibrium for O2, CH4/5N2 mix

        //--------- compute raw input sum

        double sumw = 0.0;                  // normalize w
        for(int i=0; i<var.size(); i++)
            sumw += w[i];
        for(int i=0; i<var.size(); i++)
            w[i] /= sumw;
        double sum1 = 0.0;                  // compute raw input sum
        for(int i=0; i<var.size(); i++)
            sum1 += var[i]*w[i]*rho[i];

        //--------- compute sum of HiPS projected variable

        H.set_varData(var, w, "test", rho);
        const auto& pLoc    = H.get_pLoc();
        const auto& varData = H.get_HipsVarData_ptr();

        double sum2 = 0.0;                  // compute HiPS projected sum
        for(int i=0; i<H.get_nparcels(); i++)
            sum2 += (*varData[0])[pLoc[i]] * H.varRho[pLoc[i]] * H.wPar[pLoc[i]];

        //------------ compute back_projected sum

        auto pairs = H.get_varData_with_density();
        auto var2 = pairs.first[0];      // First variable (back-projected values)
        auto rho2 = pairs.second;        // Back-projected density (shared for all variables)

        double sum3 = 0.0;
        for (int i = 0; i < var2.size(); i++) {
            sum3 += var2[i] * rho2[i] * w[i];
        }

        //------------ test sum1=sum2 and sum2=sum3

        REQUIRE_THAT( sum1, WithinRel(sum2, 1E-14) );
        REQUIRE_THAT( sum2, WithinRel(sum3, 1E-14) );

    }

    /////////////////////////////////////////////////////////////////////////

    SECTION("Test overall solve") {

        vector<double> Sc(54, 1.0);
    
        hips H(nLevels, domainLength, tau0, C_param, forceTurb, nVar, Sc, true, cantSol, 11);

        H.setOutputIntervalTime(tRun/12);  // Save results every 100 seconds

        //--------- initialize vars

        int npts = 40;           // not a power of 2, not equal to nparcels
    
        vector<vector<double>> ysp(nsp, vector<double>(npts, 0.0));
        vector<double>         h(npts);
        vector<double>         rho(npts);
        vector<double>         weight(npts, 1.0/npts);     // use uniform weights

        vector<double> y0(nsp);
        vector<double> y1(nsp);
        double T0 = 300.0; 
        double T1 = 300.0; 

        double fracBurn = 0.25;
    
        gas->setState_TPX(T0, Cantera::OneAtm, "C2H4:1, O2:3, N2:11.25");
        gas->getMassFractions(y0.data());
        double rho0 = gas->density();
        for(int j=0; j<=(1-fracBurn)*npts; ++j) { 
            h[j] = gas->enthalpy_mass();
            rho[j] = gas->density();
            for(int k=0;k<nsp;++k) 
                ysp[k][j] = y0[k]; 
        }
    
        gas->setState_TPX(T1, Cantera::OneAtm, "C2H4:1, O2:3, N2:11.25");
        gas->equilibrate("HP");
        gas->getMassFractions(y1.data());
        for(int j=((1 - fracBurn)*npts + 1); j<npts; ++j) { 
            h[j] = gas->enthalpy_mass();
            rho[j] = gas->density();
            for(int k=0; k<nsp; ++k)
                ysp[k][j] = y1[k];
        }

        //---- set conservation quantities 

        double M1 = 0.0;        // mass, raw input
        double M2 = 0.0;        // mass, hips tree before run
        double M3 = 0.0;        // mass, hips tree after run
        double M4 = 0.0;        // mass, recovered variable

        double H1 = 0.0;        // enthalpy, raw input
        double H2 = 0.0;        // enthalpy, hips tree before run
        double H3 = 0.0;        // enthalpy, hips treee after run
        double H4 = 0.0;        // enthalpy, recovered variable

        vector<double> Y1(nsp, 0.0);
        vector<double> Y2(nsp, 0.0);
        vector<double> Y3(nsp, 0.0);
        vector<double> Y4(nsp, 0.0);

        //---- conservation 1, raw input
        
        for(int j=0; j<npts; j++) {
            M1 += rho[j]*weight[j];
            H1 += rho[j]*weight[j]*h[j];
            for(int k=0; k<nsp; k++)
                Y1[k] += rho[j]*weight[j]*ysp[k][j];
        }
    
        //---- load into HiPS (enthalpy + species)

        vector<string> names = gas->speciesNames();
        names.insert(names.begin(), "enthalpy");

        H.set_varData(h, weight, names[0], rho);
        for(int k=0; k<nsp; ++k) 
            H.set_varData(ysp[k], weight, names[k+1], rho);

        //---- conservation 2, hips tree before run
        
        auto varData = H.get_HipsVarData_ptr();
        auto pLoc    = H.get_pLoc();
        for(int j=0; j<H.get_nparcels(); j++) {
            M2 += H.varRho[pLoc[j]]*H.wPar[pLoc[j]];
            H2 += H.varRho[pLoc[j]]*H.wPar[pLoc[j]]*(*varData[0])[pLoc[j]];
            for(int k=0; k<nsp; k++)
                Y2[k] += H.varRho[pLoc[j]]*H.wPar[pLoc[j]]*(*varData[k+1])[pLoc[j]];
        }

        //--------- run
        H.writeData(1, 0, 0.0 );

        H.calculateSolution(tRun, true);

        //---- conservation 3, hips tree after run
        
        varData = H.get_HipsVarData_ptr();
        pLoc    = H.get_pLoc();
        for(int j=0; j<H.get_nparcels(); j++) {
            M3 += H.varRho[pLoc[j]]*H.wPar[pLoc[j]];
            H3 += H.varRho[pLoc[j]]*H.wPar[pLoc[j]]*(*varData[0])[pLoc[j]];
            for(int k=0; k<nsp; k++)
                Y3[k] += H.varRho[pLoc[j]]*H.wPar[pLoc[j]]*(*varData[k+1])[pLoc[j]];
        }
    
        //---- conservation 4, recovered (back_projected) variables

        auto snap = H.get_varData_with_density();
        auto var2 = snap.first;
        auto rho2 = snap.second;
        for(int j=0; j<npts; ++j) { 
            M4 += rho2[j]*weight[j];
            H4 += rho2[j]*weight[j]*var2[0][j];
            for(int k=0; k<nsp; k++)
                Y4[k] += rho2[j]*weight[j]*var2[k+1][j];
        }

        //--------- test conservation

        REQUIRE_THAT( M1, WithinRel(M2, 1E-14));
        REQUIRE_THAT( M2, WithinRel(M3, 1E-14));
        REQUIRE_THAT( M3, WithinRel(M4, 1E-14));

        REQUIRE_THAT( H1, WithinRel(H2, 1E-14));
        REQUIRE_THAT( H2, WithinRel(H3, 1E-14));
        REQUIRE_THAT( H3, WithinRel(H4, 1E-14));

        for(int k=0; k<nsp; k++) {
            REQUIRE_THAT( Y1[k], WithinRel(Y2[k], 1E-14));     // species are not conserved by reaction
            REQUIRE_THAT( Y3[k], WithinRel(Y4[k], 1E-14));     // so don't compare Y2 and Y3
        }

    }
}
