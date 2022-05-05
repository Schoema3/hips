/**
 * @file domaincase_hips_soot.cc
 * Source file for class domaincase_hips_soot
 */

#include "domaincase_hips_soot.h"
#include "domain.h"
#include "dv.h"
#include "dv_rho.h"
#include "dv_dvisc.h"
#include "dv_temp.h"
#include "dv_mixf_hips.h"
#include "dv_enth_hips.h"
#include "dv_ygas_hips.h"


#include <string>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
/** domaincase_hips_soot initialization function
 *
 * @param p_domn  \input set domain pointer with.
 * @param p_phi   \input set vector pointer with.
 */

void domaincase_hips_soot::init(domain *p_domn) {

    domn = p_domn;

    vector<double> gammas(4,0.0);
    gammas[0] = 2.0/domn->gas->atomicWeight(domn->gas->elementIndex("C"));
    gammas[1] = 0.5/domn->gas->atomicWeight(domn->gas->elementIndex("H"));
    gammas[2] = -1.0/domn->gas->atomicWeight(domn->gas->elementIndex("O"));
    gammas[3] = 0.0;
    domn->strm->init(domn, gammas);

    domn->v.push_back(new dv_mixf_hips(    domn, "mixf", true,  true ));    // use same order as for Sc below
    domn->v.push_back(new dv_rho(          domn, "rho",  false, true ));
    domn->v.push_back(new dv_dvisc(        domn, "dvisc",false, true ));
    domn->v.push_back(new dv_temp(         domn, "temp", false, true ));
    domn->v.push_back(new dv_enth_hips(    domn, "enth", true,  true ));
    for(int k=0; k<domn->gas->nSpecies(); k++)
        domn->v.push_back(new dv_ygas_hips(domn, "y_"+domn->gas->speciesName(k), true, true ));

    // Add soot moments to variable list
    if (domn->pram->Lsoot) {

        string PSD_method = domn->io->sootParams["PSD_method"].as<string>();
        stringstream ss;

     
      
    }

    int ii = 0;
    domn->mixf   = domn->v.at(ii++);
    domn->rho    = domn->v.at(ii++);
    domn->dvisc  = domn->v.at(ii++);
    domn->temp   = domn->v.at(ii++);
    domn->enth   = domn->v.at(ii++);
    domn->ysp    = domn->v.begin()+ii;          // access as domn->ysp[k]->d[i], etc. where k is the species starting from 0.
    ii += domn->gas->nSpecies();
    if (domn->pram->Lsoot) {
        domn->svar  = domn->v.begin()+ii;      // access as domn->svar[k]->d[i], etc. where k is the species starting from 0
        ii += domn->pram->nsvar;
    }
                

    //------------- set Sc of transported scalars: input file needs to match this order

    if(domn->pram->LScHips){
        domn->mixf->ScHips = domn->io->scalarSc[0].as<double>();
        domn->enth->ScHips = domn->io->scalarSc[1].as<double>();
        for(int k=0; k<domn->gas->nSpecies(); k++)
            domn->ysp[k]->ScHips = domn->io->scalarSc[2].as<double>();
        if(domn->pram->Lsoot)
            for(int k=0; k<domn->pram->nsvar; k++) 
                domn->svar[k]->ScHips = domn->io->scalarSc[3].as<double>();
    }

    //-------------------- initialize profiles

    double premixed_mixf      = domn->io->initParams["premixed_mixf"].as<double>();
    double frac_burnt         = domn->io->initParams["frac_burnt"].as<double>();

    for(int i=0; i<domn->ngrd; i++)
        //domn->mixf->d.at(i) = premixed_mixf;
        //domn->mixf->d.at(i) = i<domn->ngrd/2 ? 0.055 : 0.055;
        domn->mixf->d.at(i) = i<domn->ngrd/2 ? 0.0: 1.0;

    int nsp = domn->gas->nSpecies();
    vector<double> ysp(nsp);               // dummy storage

    for(int i=0; i<domn->ngrd; i++) {
        if(i < domn->ngrd*(1.0-frac_burnt))
            domn->strm->getMixingState(premixed_mixf, ysp, domn->enth->d.at(i), domn->temp->d.at(i));
            //domn->strm->getMixingState(domn->mixf->d.at(i), ysp, domn->enth->d.at(i), domn->temp->d.at(i));
        else
            domn->strm->getProdOfCompleteComb(premixed_mixf, ysp, domn->enth->d.at(i), domn->temp->d.at(i));
            //domn->strm->getProdOfCompleteComb(domn->mixf->d.at(i), ysp, domn->enth->d.at(i), domn->temp->d.at(i));
        for(int k=0; k<nsp; k++)
            domn->ysp[k]->d.at(i) = ysp.at(k);
    }

    //-------------------- set initial soot profile

    if (domn->pram->Lsoot) {
        if (domn->pram->PSD_method == "MOMIC") { // for other types: || domn->pram->PSD_method == "QMOM" etc.
            double M0 = 1.0E0;
            double sigL = 3.0;
            double mavg = 1.0E-21;
            for (int k=0; k<domn->pram->nsvar; k++) {
                for (int j=0; j<domn->ngrd; j++) {
                    domn->svar[k]->d[j] = M0 * pow(mavg, k) * exp(0.5 * pow(k,2) * pow(sigL,2));
                }
            }
        }
    }




    //for(int i=0; i<domn->ngrd; i++)
    //    //domn->mixf->d.at(i) = double(i)/(domn->ngrd-1);
    //    if(i<256-23)
    //        domn->mixf->d.at(i) = 0.0;
    //    else if(i < 256-10)
    //        domn->mixf->d.at(i) = double(i-256+23)/12;
    //    else
    //        domn->mixf->d.at(i) = 1.0;
    //
    //int nsp = domn->gas->nSpecies();
    //vector<double> ysp(nsp);               // dummy storage
    //
    //for(int i=0; i<domn->ngrd; i++) {
    //    domn->strm->getProdOfCompleteComb(domn->mixf->d.at(i), ysp, domn->enth->d.at(i), domn->temp->d.at(i));
    //    for(int k=0; k<nsp; k++)
    //        domn->ysp[k]->d.at(i) = ysp.at(k);
    //}



    enforceMassFractions();

    domn->rho->setVar();
    domn->dvisc->setVar();

    //------------------- set minimial mesher

    vector<dv*> phi;
    

    //------------------- for variable Sc, set the Batchelor level and i_plus level for each scalar

    double iEta = domn->io->params["nLevels"].as<int>() - 3;       // nLevels may have been increased in solver_hips initializer; want the original
    for(int k=0; k<domn->v.size(); k++) {
        if(domn->v[k]->L_transported){
            if(domn->v[k]->ScHips < 1.0)
                domn->v[k]->i_batchelor = iEta + 1.5*log(domn->v[k]->ScHips)/log(4);
            else if(domn->v[k]->ScHips > 1.0)
                domn->v[k]->i_batchelor = iEta +     log(domn->v[k]->ScHips)/log(4);
            else
                domn->v[k]->i_batchelor = iEta;
            domn->v[k]->i_plus = ceil(domn->v[k]->i_batchelor);
            cout << endl << "Scalar index, Sc, i_batchelor, i_plus: " << k << " " << domn->v[k]->ScHips << " " << domn->v[k]->i_batchelor << " " << domn->v[k]->i_plus;
        }
    }
    cout << endl;

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the soltuion.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */

void domaincase_hips_soot::setGasStateAtPt(const int &ipt) {

    int nsp = domn->gas->nSpecies();
    vector<double> yi(nsp);
    for(int k=0; k<nsp; k++)
        yi.at(k) = domn->ysp[k]->d.at(ipt);

    domn->gas->setState_PY(domn->pram->pres, &yi.at(0));
    domn->gas->setState_HP(domn->enth->d.at(ipt), domn->pram->pres, 1.E-10);

}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion.
 *  These should not be transported. Those are already set.
 */
void domaincase_hips_soot::setCaseSpecificVars() {

    enforceSootMom();
    enforceMassFractions();
    domn->rho->setVar();
    domn->dvisc->setVar();
    domn->temp->setVar();
}

////////////////////////////////////////////////////////////////////////////////
/** Update/set variables that are needed in the solution.
 *  Especially for diffusion. See cvodeDriver.cc
 *  These should not be transported. Those are already set.
 */
void domaincase_hips_soot::setCaseSpecificVars_cvode(const int &ipt) {

    domn->rho->setVar(ipt);
    domn->temp->setVar(ipt);
}
