
#include <iostream>
#include "cvodeDriver.h"

#include <vector>

using namespace Cantera;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

cvodeDriver::cvodeDriver(shared_ptr<Solution> sol) {
    gas = sol->thermo();
    kin = sol->kinetics();
    Neq = gas->nSpecies();
    rr.resize(gas->nSpecies());

    integrator = unique_ptr<Integrator>(newIntegrator("CVODE"));
    integrator->setTolerances(1E-4, 1E-10);       // rtol, atol
    integrator->setMaxSteps(2000);
    integrator->initialize(0.0, *this);
}

////////////////////////////////////////////////////////////////////////////////

void cvodeDriver::integrate(double dt) {

    h_fixed = gas->enthalpy_mass();
    P_fixed = gas->pressure();

    integrator->reinitialize(0.0, *this);

    integrator->integrate(dt);

    double *solution = integrator->solution();

    gas->setMassFractions(solution);
    gas->setState_HP(h_fixed, P_fixed);

}

////////////////////////////////////////////////////////////////////////////////

void cvodeDriver::eval(double t, double *y, double *dydt, double *not_used) {

    double *massFracs = &y[0];
    double *dYdt      = &dydt[0];
   
    gas->setMassFractions_NoNorm(massFracs);
    gas->setState_HP(h_fixed, P_fixed);

    double rho = gas->density();
    kin->getNetProductionRates(&rr[0]);
    for (size_t k=0; k < gas->nSpecies(); k++)
        dYdt[k] = rr[k] * gas->molecularWeight(k) / rho;
}
