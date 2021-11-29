/**
 * @file cantera_shell_functions.h
 * Classes for Cantera shell functions
 */

#ifndef DOCANTERA

#pragma once

// This file is used to store dummy function definitions when cantera is not available

#include <string>
#include <vector>

using std::string;
using std::vector;

namespace Cantera {

const double GasConstant = 8314.46261815324;

//-------------------------------------------------------------------------------

/** Cantera shell functions when not compiling with Cantera.  
 *  When cantera is not needed, compile with these
 *  routines that have the same headers as those of Cantera.  Add new routines
 *  when new cantera functions are called.  This also provides a template for
 *  implementing user-defined functions instead of calling Cantera.  However,
 *  the routines are not meant to be called, and return errors if they are.
 *  
 *  @author David O. Lignell
 */

class IdealGasPhase {

    public :

        void            setState_PY(double p, double* y);
        void            setState_HP(double h, double p, double tol = 1.e-8);
        void            setState_TPY( double T, double p, double *y);
        double          temperature();
        double          density();
        double          cp_mass();
        double          enthalpy_mass();
        void            getEnthalpy_RT(double* hrt);
        void            getNetProductionRates(double* rr);
        double          atomicWeight(int k);
        double          molecularWeight(int k);
        double          meanMolecularWeight();
        int             speciesIndex(string name);
        int             elementIndex(string name);
        int             nSpecies();
        int             nElements();
        void            setMoleFractions(const double *x);
        void            setMassFractions(const double *y);
        void            getMassFractions(double *y);
        void            getMoleFractions(double *x);
        void            getAtoms(int k, double *Z);
        double          moleFraction(int k);
        string          speciesName(int k);
        vector<double>  molecularWeights();
        int             nAtoms(int j, int k);
        void            equilibrate(const std::string& XY, const std::string& solver="auto",
                           double rtol=1e-9, int max_steps=50000, int max_iter=100,
                           int estimate_equil=0, int log_level=0);  // inherited member from ThermoPhase

        IdealGasPhase(string mechfile);

};

//-------------------------------------------------------------------------------


/** Cantera shell functions when not compiling with Cantera.  
 *  When cantera is not needed, compile with these
 *  routines that have the same headers as those of Cantera.  Add new routines
 *  when new cantera functions are called.  This also provides a template for
 *  implementing user-defined functions instead of calling Cantera.  However,
 *  the routines are not meant to be called, and return errors if they are.
 *  
 *  @author David O. Lignell
 */

class Transport {

    public :

        double viscosity();
        double thermalConductivity();
        void getMixDiffCoeffs(double *D);

};

//-------------------------------------------------------------------------------

Transport* newTransportMgr(string s1, IdealGasPhase *gas);

} // end namespace Cantera

#endif

