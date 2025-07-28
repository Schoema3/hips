#include <catch2/catch_test_macros.hpp>
#include "hips.h"

#include <memory>

using namespace std;

//////////////////////////////////////////////////////////////////////////

TEST_CASE( "Test HiPS library" ) {

    int            nLevels = 9;                  // Number of hierarchical levels
    double         domainLength = 1.0;           // Simulation domain length
    double         tau0 = 1.0;                   // Initial time scale for the largest eddies
    double         C_param = 0.5;                // Eddy turnover rate multiplier
    double         tRun = 300.0;                 // Total simulation time
    int            forceTurb = 2;                // Forcing parameter to impose turbulent profile
    vector<double> ScHips = {0.0625, 1.0, 16.0}; // Schmidt numbers for low and high diffusivity
    int            nvars = 3;                    // Number of scalar fields

    //--------- create HiPS object

    hips H(nLevels, domainLength, tau0, C_param, forceTurb, nvars, ScHips, false);

    //--------- initialize vars

    vector<vector<double>> vars(nvars, vector<double>(H.nparcels));
    vector<double> weights(H.nparcels, 1.0/H.nparcels);  // Uniform weights

    for (int k=0; k<nvars; k++) {
        for (int i=0; i<H.nparcels; i++)
            vars[k][i] = (i<H.nparcels) ? 0 : 1;
        H.set_varData(vars[k], weights, "mixf_0" + to_string(k));
    }
    //--------- solve

    H.calculateSolution(tRun, false);

    auto v = H.varData[0];
    for(int i=0; i<H.nparcels; i++)
        cout << endl << "v = " << (*v)[i];
    cout << endl;

    //--------- test

    REQUIRE( H.nparcels == (1 << (nLevels+1)) );
}
