
#include "micromixer.h"
#include "domain.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <numeric> //accumulate

#include "interp_linear.h"

///////////////////////////////////////////////////////////////////////////////
/** micromixer constructor function
 */

void micromixer::setNominalStepSize() {
    dtStepNominal = domn->pram->diffCFL * domn->solv->tMix;
}


micromixer::micromixer() {
    cvode = new cvodeDriver();
    nsteps = 0;
}

///////////////////////////////////////////////////////////////////////////////
/** micromixer initialization function
 *
 * @param p_domn  \input set domain pointer with.
 */

void micromixer::init(domain *p_domn) {
    domn = p_domn;

    bool LincludeRhsMix = (domn->pram->Lsolver=="SEMI-IMPLICIT") ? true : false;
    cvode->init(domn, LincludeRhsMix);
}

///////////////////////////////////////////////////////////////////////////////

void micromixer::mixAcrossLevelTree(const int kVar, const int iLevel, const int iTree) {

    int istart;
    int iend;

    int nPmix = 1<<(domn->pram->nLevels - iLevel - 2);     // number of parcels mixed together
                                                           // also = # to skip in loop
                                  
    int ime;            // index of first cell (mapped to istart)
    double s;           // initialize sum to 0

    bool LdoChi = (domn->v[kVar]->var_name == "mixf" || 
                   domn->v[kVar]->var_name == "mixf_000") ? true : false;

    // TODO: generalize this to compute chi for any scalar, not just mixf, mixf_0

    //--------- Mix left branch of iTree

    istart = iTree << (domn->solv->Nm1-iLevel);  // same thing
    iend   = istart + nPmix;

    s = 0;            // initialize sum to 0
    for(int i=istart; i<iend; i++) {
        ime = domn->solv->pLoc[i];
        s += domn->v[kVar]->d[ime];
    }
    for(int i=istart; i<iend; i++) {
        ime = domn->solv->pLoc[i];

        if(LdoChi)
            domn->chi->setVar(domn->v[kVar]->d[ime], s/nPmix, ime);

        domn->v[kVar]->d[ime] = s / nPmix;
       

    }

    //--------- Mix right branch of iTree

    istart = iend;        // same thing
    iend   = istart + nPmix;

    s = 0;            // initialize sum to 0
    for(int i=istart; i<iend; i++) {
        ime = domn->solv->pLoc[i];
        s += domn->v[kVar]->d[ime];
    }
    for(int i=istart; i<iend; i++) {
        ime = domn->solv->pLoc[i];

        if(LdoChi)
            domn->chi->setVar(domn->v[kVar]->d[ime], s/nPmix, ime);

        domn->v[kVar]->d[ime] = s / nPmix;
    }
}

//--------------------------------------------------------------------
void micromixer::advanceOdt(const double p_tstart, const double p_tend, const int iLevel) { // iLevel is for hips

 tstart = p_tstart;
    tend   = p_tend;

    cout << endl << "here 1" << endl; //doldb 
    if(domn->pram->forceHips==2 && iLevel==0)   // forcing for statistically stationary
        forceProfile();
    cout << endl << "here 2" << endl; //doldb 


    if(domn->pram->LsimpleMix){
        for(int k=0; k<domn->v.size(); k++){
    cout << endl << "here 3 " << domn->v[k]->var_name << " " << domn->v.size() << " " << endl; //doldb 
            if(!domn->v[k]->L_transported)
                continue;
            if(iLevel == domn->v[k]->i_plus-1 &&                                       // mix scalar across level i_minus
               domn->rand->getRand() <= domn->v[k]->i_plus-domn->v[k]->i_batchelor)    // with probability i_plus - i_batchelor
               // for(int iTree=0; iTree<(1<<iLevel); iTree++) 
                    mixAcrossLevelTree(k, iLevel,domn->solv-> iTree);  
            else if(iLevel >= domn->v[k]->i_plus ){                                    // scalars at or below i_plus are fully mixed
               // for(int iTree=0; iTree<(1<<iLevel); iTree++)                           // note: the >= could be more efficient as == if 
                    mixAcrossLevelTree(k, iLevel,domn->solv-> iTree);                              //    the initial condition has uniform mixing below Batchelor
            }
        }
    }

    //if(domn->pram->LsimpleMix){{{
    //    for(int k=0; k<domn->v.size(); k++)
    //        if(domn->v.at(k)->L_transported){
    //            for(int i=0; i<domn->ngrd; i+=2){
    //                int ime = domn->solv->pLoc[i];
    //                int inb = domn->solv->pLoc[i+1];
    //                double val = (domn->v[k]->d[ime] + domn->v[k]->d[inb])/2.0;
    //                domn->v[k]->d[ime] = val;
    //                domn->v[k]->d[inb] = val;
    //            }
    //        }}}}

    setNominalStepSize();

    for(time=tstart; time<tend; time+=dt) {

        setStepSize();
        if(domn->pram->Lsolver=="EXPLICIT")
            advanceOdtSingleStep_Explicit();
        else if(domn->pram->Lsolver=="SEMI-IMPLICIT")
            advanceOdtSingleStep_SemiImplicit();
        else if(domn->pram->Lsolver=="STRANG")
            advanceOdtSingleStep_StrangSplit();

        domn->io->dumpDomainIfNeeded();

    }
   
}


///////////////////////////////////////////////////////////////////////////////

void micromixer::forceProfile() {

    for(int k=0; k<domn->v.size(); k++){
        if(!domn->v[k]->L_transported)
            continue;

        double s;
        //---------- force the left half of parcels to average 0

        s=0.0;
        for(int i=0; i<domn->ngrd>>1; i++)
            s += domn->v[k]->d[domn->solv->pLoc[i]];
        s /= (domn->ngrd>>1);
        for(int i=0; i<domn->ngrd>>1; i++)
            domn->v[k]->d[domn->solv->pLoc[i]] += (-s - 0.0);

        //---------- force the right half of parcels to average 1

        s=0.0;
        for(int i=domn->ngrd>>1; i<domn->ngrd; i++)
            s += domn->v[k]->d[domn->solv->pLoc[i]];
        s /= (domn->ngrd>>1);
        for(int i=domn->ngrd>>1; i<domn->ngrd; i++)
            domn->v[k]->d[domn->solv->pLoc[i]] += (-s + 1.0);
    }
}


////////////////////////////////////////////////////////////////////////////////
/** Set time step size.
 *  This is based on a diffusive (or other) timescale.
 *  This is a uniform step size for the given integration period.
 *  The actual step size will be rest based on a data dump or tend.
 */

void micromixer::setStepSize() {

    dt = dtStepNominal;

    if (time+dt > domn->io->dumpTimes.at(domn->io->iNextDumpTime)) {
        dt = domn->io->dumpTimes.at(domn->io->iNextDumpTime) - time;
        domn->io->LdoDump = true;
    }

    if (time + dt > tend) {
        dt = (tend - time)*(1.0+1.0E-8);
        domn->io->LdoDump = false;
    }
}

///////////////////////////////////////////////////////////////////////////////
/** Advance ODT solution: diffusion and reaction
 */

void micromixer::advanceOdtSingleStep_Explicit(){

 


    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported) {
            domn->v.at(k)->getRhsMix(gf, dxc);
            domn->v.at(k)->getRhsSrc();
        }

    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported)
            for(int i=0; i < domn->ngrd; i++)
                domn->v.at(k)->d.at(i) = domn->v.at(k)->d.at(i) + dt*( domn->v.at(k)->rhsMix.at(i) + domn->v.at(k)->rhsSrc.at(i));

  


}

///////////////////////////////////////////////////////////////////////////////
/** Advance ODT solution: diffusion and reaction; Some terms are implicit, others explicit.
 *  Nominally the mixing terms are explicit. Calling the cvode driver.
 *  First order.
 *  dphi/dt =  D(phi_0) + S(phi) : solving from t0 to t1.
 *  Here, D() is the diffusive term, and S() is the (stiff) source term.
 *  We solve the whole RHS implicitly, but the D(phi_0) is fixed at time 0.
 */

void micromixer::advanceOdtSingleStep_SemiImplicit() {

    if(domn->pram->Lsolver!="SEMI-IMPLICIT")
        return;

   // setGridDxcDx();
   // setGf();
    domn->domc->setCaseSpecificVars();
    //set_oldrho_or_rhov();
    if(domn->pram->Lspatial) transform(oldrho_or_rhov.begin(), oldrho_or_rhov.end(), domn->uvel->d.begin(), oldrho_or_rhov.begin(), multiplies<double>());
    //if(domn->pram->LdoDL) do_DL("set DL_1");

    //--------------- Set the explicit (mixing) terms

    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported)
            domn->v.at(k)->getRhsMix(gf, dxc);

    //--------------- Perform the implicit integration on each cell

    for(int i=0; i<domn->ngrd; i++)
        cvode->integrateCell(i, dt);

    //---------------

  }

///////////////////////////////////////////////////////////////////////////////
/** Advance ODT solution: diffusion and reaction; Some terms are implicit, others explicit.
 *  Nominally the mixing terms are explicit. Calling the cvode driver.
 *  First order.
 *  dphi/dt =  D(phi_0) + S(phi) : solving from t0 to t1.
 *  Here, D() is the diffusive term, and S() is the (stiff) source term.
 *  We solve the whole RHS implicitly, but the D(phi_0) is fixed at time 0.
 */

void micromixer::advanceOdtSingleStep_StrangSplit() {

    if(domn->pram->Lsolver!="STRANG")
        return;

    //--------------- First step: phi_1 = phi_0 + 0.5*dt*D(phi_0)

    domn->domc->setCaseSpecificVars();
 
    if(domn->pram->Lspatial) transform(oldrho_or_rhov.begin(), oldrho_or_rhov.end(), domn->uvel->d.begin(), oldrho_or_rhov.begin(), multiplies<double>());

    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported)
            domn->v.at(k)->getRhsMix(gf, dxc);

    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported)
            for(int i=0; i < domn->ngrd; i++)
                domn->v.at(k)->d.at(i) = domn->v.at(k)->d.at(i) + 0.5*dt*domn->v.at(k)->rhsMix.at(i);

 

  

    domn->domc->setCaseSpecificVars();

    if(domn->pram->Lspatial) transform(oldrho_or_rhov.begin(), oldrho_or_rhov.end(), domn->uvel->d.begin(), oldrho_or_rhov.begin(), multiplies<double>());

    domn->v[0]->resetSourceFlags();             // sets L_source_done = false for all transported vars
    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported)
            domn->v.at(k)->getRhsSrc();

    for(int i=0; i<domn->ngrd; i++)
        cvode->integrateCell(i, dt);

 
    domn->domc->setCaseSpecificVars();

    if(domn->pram->Lspatial) transform(oldrho_or_rhov.begin(), oldrho_or_rhov.end(), domn->uvel->d.begin(), oldrho_or_rhov.begin(), multiplies<double>());

    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported)
            domn->v.at(k)->getRhsMix(gf, dxc);

    for(int k=0; k<domn->v.size(); k++)
        if(domn->v.at(k)->L_transported)
            for(int i=0; i < domn->ngrd; i++)
                domn->v.at(k)->d.at(i) = domn->v.at(k)->d.at(i) + 0.5*dt*domn->v.at(k)->rhsMix.at(i);


    //-------------------------

}
///////////////////////////////////////////////////////////////////////////////////



