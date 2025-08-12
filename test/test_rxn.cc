#include <catch2/catch_all.hpp>

#include "hips.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo.h"

#include <memory>
#include <cmath>
#include <format>
using Catch::Approx;

using namespace std;

//////////////////////////////////////////////////////////////////////////

TEST_CASE( "Test HiPS library" ) {

    int            nLevels = 6;          // Number of hierarchical levels
    double         domainLength = 0.01;  // Domain length scale
    double         tau0 = 0.0018;      // Mixing timescale (fast mixing)
    double         C_param = 0.5;        // Eddy rate multiplier
    double         tRun = 0.000000018;       // Total simulation runtime
    int            forceTurb = 0;        // No forced turbulence
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
        //H.setDensityWeightedMixing(true);


        vector<double> var = {0.0, 0.1, 0.2, 0.8, 0.9, 1.0};               // mixture fraction
        //vector<double> w   = {0.5, 0.6, 0.7, 0.7, 0.6, 0.2};               // test: not a power of 2, not uniform, don't sum to 1
        vector<double> w(6,1.0/6); //   = {0.5, 0.6, 0.7, 0.7, 0.6, 0.2};               // test: not a power of 2, not uniform, don't sum to 1
        vector<double> rho = {1.300, 0.474, 0.303, 0.171, 0.335, 1.057};   // equilibrium for O2, CH4/5N2 mix
        int nparcels = H.get_nparcels();


        double sumw = 0.0;                  // normalize w
        for(int i=0; i<var.size(); i++)
            sumw += w[i];
        for(int i=0; i<var.size(); i++)
            w[i] /= sumw;
        double sum1 = 0.0;                  // compute raw input sum
        for(int i=0; i<var.size(); i++)
            sum1 += var[i]*w[i]*rho[i];

        H.set_varData(var, w, "test", rho);
        const auto& pLoc    = H.get_pLoc();
        const auto& varData = H.get_varData_ptr();

        double sum2 = 0.0;                  // compute HiPS projected sum
        for(int i=0; i<nparcels; i++)
            sum2 += (*varData[0])[pLoc[i]] * H.varRho[pLoc[i]] * H.wPar[pLoc[i]];



        REQUIRE( abs((sum1 - sum2)/sum1) < 1E-8 );   // results equal to within roundoff error

        //------------

        auto pairs = H.get_varData_with_density();
        auto var2 = pairs.first[0];      // First variable (back-projected values)
        auto rho2 = pairs.second;        // Back-projected density (shared for all variables)

        double sum3 = 0.0;
        for (int i = 0; i < var2.size(); i++) {
            sum3 += var2[i] * rho2[i] * w[i];
        }
 

        REQUIRE( abs((sum3 - sum2)/sum3) < 1E-8 );   // results equal to within roundoff error
    }

/////////////////////////////////////////////////////////////////////////
    SECTION("Test overall solve") {

        vector<double> Sc(54, 1.0);
        int nVar = nsp + 1;
    
        hips H(nLevels, domainLength, tau0, C_param, forceTurb, nVar, Sc, true, cantSol, 10);
        H.setDensityWeightedMixing(true);                                // Perform mass_weighted_mixing
    
        int npar = H.get_nparcels();
        vector<double> weight(npar, 1.0 / npar);
    
        // ---- init fresh / burnt

        vector<vector<double>> ysp(nsp, vector<double>(npar, 0.0));
        vector<double> h(npar), rho(npar), y0(nsp), y1(nsp);
        double T0 = 300.0, T1 = 300.0, h0, h1;
        double fracBurn = 0.25;
    
        gas->setState_TPX(T0, Cantera::OneAtm, "C2H4:1, O2:3, N2:11.25");
        h0 = gas->enthalpy_mass();  gas->getMassFractions(y0.data());
        double rho0 = gas->density();
        for (int j = 0; j <= (1 - fracBurn) * npar; ++j) { h[j] = h0; for (int k=0;k<nsp;++k) ysp[k][j] = y0[k]; }
    
        gas->setState_TPX(T1, Cantera::OneAtm, "C2H4:1, O2:3, N2:11.25");
        gas->equilibrate("HP");     gas->getMassFractions(y1.data());
        double rho1b = gas->density(); h1 = gas->enthalpy_mass();
        for (int j = ((1 - fracBurn) * npar + 1); j < npar; ++j) { h[j] = h1; for (int k=0;k<nsp;++k) ysp[k][j] = y1[k]; }
    
        int fresh_end = (int)((1 - fracBurn) * npar);
        for (int j=0; j<fresh_end; ++j) rho[j] = rho0;
        for (int j=fresh_end; j<npar;   ++j) rho[j] = rho1b;
    
        // ---- load into HiPS (enthalpy + species)
        vector<string> names(nsp + 1);
        names[0] = "enthalpy";
        for (int k=0; k<nsp; ++k) names[k+1] = gas->speciesName(k);
    
        H.set_varData(h, weight, names[0], rho);
        for (int k=0; k<nsp; ++k) H.set_varData(ysp[k], weight, names[k+1], rho);
    
        // ---- Before snapshot (CFD space)

        auto snap1 = H.get_varData_with_density();
        const auto& var1 = snap1.first;
        const auto& rho1 = snap1.second;
    
        // run
        H.calculateSolution(tRun, true);
    
        // ---- AFTER snapshot (CFD space)
        auto snap2 = H.get_varData_with_density();
        const auto& var2 = snap2.first;
        const auto& rho2 = snap2.second;
    
        // mass conservation in CFD space
        double M1 = 0.0, M2 = 0.0;

        for (int i=0;i<npar;++i) { M1 += rho1[i]*weight[i]; M2 += rho2[i]*weight[i]; }

        REQUIRE( M2 == Approx(M1).epsilon(1e-9).margin(1e-12) );
    
        // Sum of Y = 1 per parcel (parcel space)
        const auto& pLoc    = H.get_pLoc();
        const auto& varData = H.get_varData_ptr();

        for (int j=0; j<npar; ++j) {
            int p = pLoc[j];
            long double sumY = 0.0L;

            for (int k=1; k<nVar; ++k) sumY += (long double)(*varData[k])[p];

            REQUIRE( (double)sumY == Approx(1.0).epsilon(1e-8).margin(1e-12) );
        }
    
    // CO2 behavior

    int idxCO2 = gas->speciesIndex("CO2") + 1;
    double CO2_before = 0.0, CO2_after = 0.0;

    for (int i=0;i<npar;++i) {
        CO2_before += var1[idxCO2][i]*rho1[i]*weight[i];
        CO2_after  += var2[idxCO2][i]*rho2[i]*weight[i];
    }
    REQUIRE( CO2_after >= CO2_before - 1e-12 );

    // some CO2 must exist somewhere
    double yco2_max = 0.0;
    for (int j=0;j<npar;++j) yco2_max = std::max(yco2_max, (*varData[idxCO2])[pLoc[j]]);
    REQUIRE( yco2_max > 1e-6 );

    // make sure number of parcels matches structure
    REQUIRE( npar == (1 << (nLevels-1)) );
    }

    
}
