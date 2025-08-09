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
    double         tau0 = 0.0005;      // Mixing timescale (fast mixing)
    double         C_param = 0.5;        // Eddy rate multiplier
    double         tRun = 0.0013;       // Total simulation runtime
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

        vector<double> var = {0.0, 0.1, 0.2, 0.8, 0.9, 1.0};               // mixture fraction
        //vector<double> w   = {0.5, 0.6, 0.7, 0.7, 0.6, 0.2};               // test: not a power of 2, not uniform, don't sum to 1
        vector<double> w(6,1.0/6); //   = {0.5, 0.6, 0.7, 0.7, 0.6, 0.2};               // test: not a power of 2, not uniform, don't sum to 1
        vector<double> rho = {1.300, 0.474, 0.303, 0.171, 0.335, 1.057};   // equilibrium for O2, CH4/5N2 mix

        double sumw = 0.0;                  // normalize w
        for(int i=0; i<var.size(); i++)
            sumw += w[i];
        for(int i=0; i<var.size(); i++)
            w[i] /= sumw;
        double sum1 = 0.0;                  // compute raw input sum
        for(int i=0; i<var.size(); i++)
            sum1 += var[i]*w[i]*rho[i];

        H.set_varData(var, w, "test", rho);
        double sum2 = 0.0;                  // compute HiPS projected sum
        for(int i=0; i<H.nparcels; i++)
            sum2 += (*H.varData[0])[H.pLoc[i]] * H.varRho[H.pLoc[i]] * H.wPar[H.pLoc[i]];



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

    //////////////////////////////////////////////////////////////////////////

    SECTION("Test overall solve") {
    
        vector<double> ScHips(54, 1);     // Schmidt number (unity for all species)
        int nVar = nsp + 1;               // Number of variables (species + enthalpy)

        //---------
    
        hips H(nLevels, domainLength, tau0, C_param, forceTurb, nVar, ScHips, true, cantSol, 10);
    
        //--------- initialize vars

        vector<vector<double>> ysp(nsp, vector<double>(H.nparcels, 0));     // species
        vector<double> h(H.nparcels);                                       // enthalpy
        vector<double> T(H.nparcels);                                       // temperature

        vector<double> y0(nsp);    // fresh reactant mass fractions                              
        vector<double> y1(nsp);    // burnt reactant mass fractions;
        vector<double> rho(H.nparcels);  // density array

        double T0 = 300.0;         // fresh temperature
        double T1 = 300.0;         // burnt temperature
        double h0;                 // fresh enthalpy
        double h1;                 // burnt enthalpy

        double fracBurn = 0.25;                // fraction of parcels pre-combusted

        // Initialize fresh reactants
        gas->setState_TPX(T0, Cantera::OneAtm, "C2H4:1, O2:3, N2:11.25");
        h0 = gas->enthalpy_mass();
        gas->getMassFractions(y0.data());
        double rho0 = gas->density();

        // Assign fresh reactant properties to parcels
        for (int i = 0; i < nsp; i++) {
            for (int j = 0; j <= (1 - fracBurn) * H.nparcels; j++) {
                h[j] = h0;
                T[j] = gas->temperature();
                ysp[i][j] = y0[i];
            }
        }

        // Initialize pre-combusted parcels
        gas->setState_TPX(T1, Cantera::OneAtm, "C2H4:1, O2:3, N2:11.25");
        gas->equilibrate("HP");
        gas->getMassFractions(y1.data());
        double rho_burnt = gas->density();  // burnt gas density
        h1 = gas->enthalpy_mass();  //  assign h1 for burnt gas



        // Assign pre-combusted properties to parcels
        for (int i = 0; i < nsp; i++) {
            for (int j = ((1 - fracBurn) * H.nparcels + 1); j < H.nparcels; j++) {
                h[j] = h1;
                T[j] = gas->temperature();
                ysp[i][j] = y1[i];
            }
        }
        int fresh_end = static_cast<int>((1 - fracBurn) * H.nparcels);  // upper index for fresh

        for (int j = 0; j < fresh_end; j++) {
            rho[j] = rho0;
        }
        for (int j = fresh_end; j < H.nparcels; j++) {
            rho[j] = rho_burnt;
        }
        
        // Assign variables to HiPS
        vector<string> variableNames(nsp + 1); // Variable names: enthalpy and species
        variableNames[0] = "enthalpy";
        for (int k=0; k<ysp.size(); k++)
            variableNames[k+1] = gas->speciesName(k);

        vector<double> weight(H.nparcels, 1.0 / H.nparcels); // Uniform weights

        H.set_varData(h, weight, variableNames[0], rho);  //

          //todo: compute vector of sum1 for all of the variables
        for (int k = 0; k < ysp.size(); k++) 
            H.set_varData(ysp[k], weight, variableNames[k + 1], rho);  // includes density

        //--------- solve
       
        auto initial = H.get_varData_with_density();
        const std::vector<std::vector<double>>& var1 = initial.first;
        const std::vector<double>& rho1 = initial.second;

        std::vector<double> sum1(nVar, 0.0);
        for (int k = 0; k < nVar; ++k) {
            for (int i = 0; i < H.nparcels; ++i) {
                sum1[k] += var1[k][i] * rho1[i] * H.wPar[H.pLoc[i]];

            } 
        }

        for (int k = 0; k < nVar; ++k) {
            std::cout << "sum1[" << k << "] = " << sum1[k] << std::endl;
        }

        H.calculateSolution(tRun, true);
        
        auto result = H.get_varData_with_density();
        const std::vector<std::vector<double>>& var2 = result.first;
        const std::vector<double>& rho2 = result.second;
        
        std::vector<double> sum2(nVar, 0.0);
        for (int k = 0; k < nVar; ++k) {
            for (int i = 0; i < H.nparcels; ++i) {
                sum2[k] += var2[k][i] * rho2[i] * H.wPar[H.pLoc[i]];

            }
        }
        
        for (int k = 0; k < nVar; ++k) {
            REQUIRE(std::abs(sum2[k] - sum1[k]) < 1e-12);


           // REQUIRE(sum2[k] == Approx(sum1[k]).epsilon(1e-8));
        }

        //--------- test

        REQUIRE( H.nparcels == (1 << (nLevels-1)) );                     // based on ScHips set above
        REQUIRE( abs((*H.varData[16])[H.pLoc[11]] - 0.11312593835096947) < 1E-5);   // CO2 mass fraction at parcel index 11
        REQUIRE( abs(H.Temp[H.pLoc[11]] - 1887.8571573335937) < 1E-0 );             // Temperature at parcel index 11    

       
        }
}
