/**
 * @file dv_chi.cc
 * Source file for class dv_chi
 */

#include "dv_chi_hips.h"
#include "domain.h"
#include <cstdlib>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
/** dv_chi  constructor function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

dv_chi_hips::dv_chi_hips(domain     *line,
                         const      string s,
                         const bool Lt,
                         const bool Lo) {

    domn          = line;
    var_name      = s;
    L_transported = Lt;
    L_output      = Lo;

    d             = vector<double>(domn->ngrd, 0.0);

    tLC           = vector<double>(domn->ngrd, 0.0);

}

////////////////////////////////////////////////////////////////////////////////
/*! Set chi (scalar dissipation rate)
 *  @param vtilde \input variable state before mixing
 *  @param vhat   \input variable state after mixing
 *  @param ipt    \input optional point to compute at given index
 */

void dv_chi_hips::setVar(const double vtilde, const double vhat, const int ipt){

    if(ipt == -1){
        cout << endl << "ERROR in dv_chi_hips::setVar; only call ipt != -1" << endl;
        exit(0);
    }
    
    double delta_t_chi = domn->solv->time - tLC[ipt];
    d[ipt] = 2*(vtilde - vhat)*(vtilde-vhat)/delta_t_chi;
    cout<<"d-->  "<<d[ipt]<<endl;
    tLC[ipt] = domn->solv->time;

}


