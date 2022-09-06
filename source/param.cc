/**
 * @file param.cc
 * @brief Source file for class \ref param
 */

#include "param.h"
#include "domain.h"

///////////////////////////////////////////////////////////////////////////////
/** param initialization function
 * @param p_domn  \input set domain pointer with.
 */

void param::init(domain *p_domn) {
    domn = p_domn;
}

///////////////////////////////////////////////////////////////////////////////
/** param constructor function
 *
 * @param p_io  \input set inputoutput pointer with.
 */

param::param(inputoutput *p_io) {

    io = p_io;

    seed           = io->params["seed"]           ? io->params["seed"].as<int>()             : -1;
    tEnd           = io->params["tEnd"]           ? io->params["tEnd"].as<double>()          : errMsg<double>("tEnd");
    domainLength   = io->params["domainLength"]   ? io->params["domainLength"].as<double>()  : errMsg<double>("domainLength");
    ngrd0          = io->params["ngrd0"]          ? io->params["ngrd0"].as<int>()            : 1000;     //errMsg<int>("ngrd0");
    rho0           = io->params["rho0"]           ? io->params["rho0"].as<double>()          : 1.0;      //errMsg<double>("rho0");
    kvisc0         = io->params["kvisc0"]         ? io->params["kvisc0"].as<double>()        : 0.001694; //errMsg<double>("kvisc0");
    chemMechFile   = io->params["chemMechFile"]   ? io->params["chemMechFile"].as<string>()  : errMsg<string>("chemMechFile");
    probType       = io->params["probType"]       ? io->params["probType"].as<string>()      : errMsg<string>("probType");

    C_param        = io->params["C_param"]        ? io->params["C_param"].as<double>()       : 5.0;      //errMsg<double>("C_param");
    diffCFL        = io->params["diffCFL"]        ? io->params["diffCFL"].as<double>()       : errMsg<double>("diffCFL");
    cvode_atol     = io->params["cvode_atol"]     ? io->params["cvode_atol"].as<double>()    : 1.0E-10;
    cvode_rtol     = io->params["cvode_rtol"]     ? io->params["cvode_rtol"].as<double>()    : 1.0E-4;
    Lsolver        = io->params["Lsolver"]        ? io->params["Lsolver"].as<string>()       : errMsg<string>("Lsolver");
    modDump        = io->params["modDump"]        ? io->params["modDump"].as<int>()          : 1000000; //errMsg<int>("modDump");
   
    pres           = io->params["pres"]           ? io->params["pres"].as<double>()          : 101325.0;

    Ltecplot       = io->params["Ltecplot"]       ? io->params["Ltecplot"].as<bool>()        : false; //errMsg<int>("modDump");


    Lrestart       = io->params["Lrestart"]       ? io->params["Lrestart"].as<bool>()        : false;
    rstType        = io->params["rstType"]        ? io->params["rstType"].as<string>()       : "single";    // "single" or "multiple"
    trst = 0.0; // (dont read this in, it comes from the restart file




    // HIPS variables ---------------------

    LisHips        = io->params["LisHips"]        ? io->params["LisHips"].as<bool>()         : false;
    nLevels        = io->params["nLevels"]        ? io->params["nLevels"].as<int>()          : 0;
    Afac           = io->params["Afac"]           ? io->params["Afac"].as<double>()          : 0.5;
    tau0           = io->params["tau0"]           ? io->params["tau0"].as<double>()          : 0.0;
    fmix           = io->params["fmix"]           ? io->params["fmix"].as<double>()          : 0.0;
    LScHips        = io->params["LScHips"]        ? io->params["LScHips"].as<bool>()         : false;
    LsimpleMix     = io->params["LsimpleMix"]     ? io->params["LsimpleMix"].as<bool>()      : false;
    forceHips      = io->params["forceHips"]      ? io->params["forceHips"].as<int>()        : -1;

}

