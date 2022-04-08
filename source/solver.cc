
#include "solver.h"

#include "domain.h"
#include <cmath>     // pow, log, ceil
#include <algorithm> // copy, sort
#include <sstream>
#include <string>
#include <iomanip>   // precision///////////////////////////////////////////////////////////////////////////////
/** solver initialization function
 *
 * @param p_domn  \input set domain pointer with.
 */


///////////////////////////////////////////////////////////////////////////////

/** destructor
 */

bool sortFunc(pair<double,int> &a, pair<double,int> &b){
    return a.first < b.first;
}


void solver::init(domain *p_domn) {

    domn    = p_domn;

//------------------- Set number of parcels, and level lengthscales, timescales, and rates

    iEta = domn->pram->nLevels - 3;
    if(domn->pram->LScHips){
        int maxSc=1.0;
        for(int i=0; i<domn->io->scalarSc.size(); i++)
            maxSc = domn->io->scalarSc[i].as<double>() > maxSc ? domn->io->scalarSc[i].as<double>() : maxSc;
        if(maxSc > 1.0) {
            domn->pram->nLevels += ceil(log(maxSc)/log(4));
            cout << endl << endl << "maxSc > 1.0, increasing nLevels to " << domn->pram->nLevels << endl << endl;
        }
    }

    Nm1 = domn->pram->nLevels - 1;
    Nm2 = domn->pram->nLevels - 2;
    Nm3 = domn->pram->nLevels - 3;

    domn->ngrd = static_cast<int>(pow(2, Nm1));

    vector<double> levelLengths(domn->pram->nLevels);    // including all levels, but last 2 don't count:
    vector<double> levelTaus(domn->pram->nLevels);       //     smallest scale is 2 levels up from bottom
    levelRates   = vector<double>(domn->pram->nLevels);

    for(int i=0; i<domn->pram->nLevels; i++){
        levelLengths[i] = domn->pram->domainLength * pow(domn->pram->Afac,i);
        levelTaus[i]    = domn->pram->tau0 *
                          pow(levelLengths[i]/domn->pram->domainLength, 2.0/3.0) /
                          domn->pram->C_param;
        levelRates[i]   = 1.0/levelTaus[i] * pow(2.0,i);
    }
    if(domn->pram->LScHips){     // correct levels for high Sc (levels above Kolmogorov)
        for(int i=iEta+1; i<domn->pram->nLevels; i++) {
            levelTaus[i]    = domn->pram->tau0 *
                              pow(levelLengths[iEta]/domn->pram->domainLength, 2.0/3.0) /
                              domn->pram->C_param;
            levelRates[i]   = 1.0/levelTaus[i] * pow(2.0,i);
            cout<<"just print<"<<i<<endl;
                    }
    }

    //-------------------

    eddyRate_total = 0.0;
    for(int i=0; i<=Nm3; i++)
        eddyRate_total += levelRates[i];

    eddyRate_inertial = 0.0;
    for(int i=0; i<=iEta; i++)
        eddyRate_inertial += levelRates[i];

    //-------------------

    tMix = domn->pram->fmix * levelTaus[domn->pram->nLevels-2];

    //------------------- Set the parcel addresses (index array)

    pLoc.resize(domn->ngrd);
    for(int i=0; i<domn->ngrd; i++)
        pLoc[i] = i;

    //------------------- set the eddy event times at all levels

    // setEddyEventTimes();          // precomputed list: not needed when sampling eddies on the fly

    //------------------- verify settings

    if(domn->pram->LScHips && !domn->pram->LsimpleMix){
        cout << endl << "\nERROR: LScHips requires LsimpleMix\n" << endl;
        exit(0);
    }}

solver::~solver(){};

void solver::setEddyEventTimes(){

    double dmb;

    for(int i=0; i<levelRates.size()-2; i++) {
        double rate = levelRates[i];
        double t = 0.0;
        for(;;) {
            double r = domn->rand->getRand();
            t += -log(r)/rate;
            if(t > domn->pram->tEnd)
                break;
            eTL.push_back(make_pair(t, i));
        }
    }

    //--------- now sort the list of times and the associated levels

    sort(eTL.begin(), eTL.end(), sortFunc);

}  


void solver::sample_hips_eddy(double &dt, double &iLevel) {

    static double c1 = 1.0 - pow(2.0, 5.0/3.0*(iEta+1));
    static double c2 = pow(2.0, Nm2) - pow(2.0, iEta+1);
    static double c3 = pow(2.0, iEta+1);

    //--------------- time to next eddy

    double r = domn->rand->getRand();
    dt = -log(r)/eddyRate_total;

    //--------------- get eddy level

    r = domn->rand->getRand();

    if(r <= eddyRate_inertial/eddyRate_total) {   // inertial region
        r = domn->rand->getRand();
        iLevel = ceil(3.0/5.0*log2(1.0-r*c1) - 1.0);
        if(iLevel < 0)    iLevel = 0;
        if(iLevel > iEta) iLevel = iEta;
    }
    else {                                        // "Batchelor" region
        r = domn->rand->getRand();
        iLevel = ceil(log2(r*c2 + c3) - 1.0);
        if(iLevel < iEta+1) iLevel = iEta+1;
        if(iLevel > Nm3) iLevel = Nm3;
    }

    return;
}

void solver::selectAndSwapTwoSubtrees(const int iLevel, int &Qstart, int &Rstart, int &nPswap){

    iTree = domn->rand->getRandInt((1 << iLevel)-1);          
    int zero_q = domn->rand->getRandInt(1);                       // 0q where q is 0 or 1
    int one_r  = 2+domn->rand->getRandInt(1);                     // 1r where r is 0 or 1

    Qstart = (zero_q << (Nm3-iLevel)) + (iTree << (Nm1-iLevel));  // starting index of Q parcels
    Rstart = (one_r  << (Nm3-iLevel)) + (iTree << (Nm1-iLevel));  // starting index of R parcels
    nPswap = 1 << (Nm3-iLevel);                                   // number of parcels that will be swapped

    int Qend = Qstart + nPswap;                      // inclusive indices are Qstart to Qend-1
    int Rend = Rstart + nPswap;                      // inclusive indices are Rstart to Rend-1

    vector<int> aa(pLoc.begin()+Qstart, pLoc.begin()+Qend);
    copy(pLoc.begin()+Rstart, pLoc.begin()+Rend, pLoc.begin()+Qstart);  // python: pLoc[Qstart:Qend]=pLoc[Rstart:Rend]
    copy(aa.begin(), aa.end(), pLoc.begin()+Rstart);                    // python: pLoc[Rstart:Rend]=aa


}

///////////////////////////////////////////////////////////////////////////////
/** The actual solver
 *
 * Advance in time only sampling eddies till get one (EE).  Then diffuse the system
 * until you catch up to the eddy time:
 * <code><pre>
 *     last EE    this EE       diffuse to catch up
 * .....|..........|           ---------------------->  ................|
 *      t0         time                                                t0,time
 *                          (last EE)   this EE         next EE
 * Then advance some more ......;.........|................|         etc.
 *                                        t0              time
 * </pre></code>
 * Eddy timesteps are smaller than required diffusive steps, so if we were to
 * lock-step the two processes we would do too much diffusive work (i.e., take
 * more diffusive steps than needed)
 *
 */

void solver::calculateSolution() {

    domn->io->writeDataFile("hips_init.dat", 0.0);

    //--------------------

    unsigned long  iEE = 0;                // eddy counter
    int    QS, RS, nPs;
    //double time;                   // current simulation time
    double dt;                     // time to next eddy event (EE)
    double iLevel;                 // current level of EE
    double iLevel_p;               // previous level of EE

    time = 0.0;                          // init time
    iLevel_p  = -1;                      // init init prev level
    sample_hips_eddy(dt, iLevel);        // init next EE

    while(time+dt <= domn->pram->tEnd) {
        domn->mimx->advanceOdt(time, time+dt, iLevel_p);       //--- ADVANCE to sampled EE ---
        selectAndSwapTwoSubtrees(iLevel, QS, RS, nPs);         //--- IMPLEMENT  sampled EE ---

        if(++iEE % 1000 == 0)
            *domn->io->ostrm << endl << "eddy #, time, level: " << iEE << "  " << time+dt << "  " << iLevel;
        if(iEE % domn->pram->modDump == 0){
            stringstream ss1; string s1;
            ss1.clear();  ss1 << setfill('0') << setw(4) << iEE; ss1 >> s1;
            domn->io->writeDataFile("hips_eddy_"+s1+".dat", time+dt);
        }

        time += dt;                      // init time
        iLevel_p = iLevel;               // init prev level
        sample_hips_eddy(dt, iLevel);    // init next EE

    }
    domn->mimx->advanceOdt(time, domn->pram->tEnd, iLevel_p);    //--- ADVANCE ---

    //--------------------

    *domn->io->ostrm << endl;
    domn->io->writeDataFile("hips_final.dat", 0.0);

}

///////////////////////////////////////////////////////////////////////////////
/**dtSmean is computed, which is the mean eddy sample time.  The Poisson
 * process draws dt's with this mean.
 */

void solver::computeDtSmean() {

    if(domn->pram->Llem) ;
    //    dtSmean = 1.0/(domn->pram->lemRateParam * domn->Ldomain());

    else if(!domn->pram->Lspatial)
        dtSmean = 0.1*domn->pram->Pav * domn->Ldomain() * domn->Ldomain() /
                  domn->pram->kvisc0 / domn->ngrd / domn->ngrd / domn->ngrd;
    else
        dtSmean = 10.*domn->pram->Pav * domn->Ldomain() / domn->ngrd / domn->ngrd;

}

///////////////////////////////////////////////////////////////////////////////
/**Sample the eddy trial time step with mean dtSmean.
 * @return Poisson sampled timestep.
 */

double solver::sampleDt() {

    return -dtSmean*log( max(1.0e-14, domn->rand->getRand()) );

}

///////////////////////////////////////////////////////////////////////////////
/**Computes dtCUmax (as the name suggests).  This variable is the time
 * increment of eddy trial time advancement before we diffuse to catch
 * up to that time in the event of no eddy before that time.
 */

void solver::computeDtCUmax() {

    double dxmin = 1.0E10;
    double d1;
    for(int i=1; i<domn->ngrdf; i++) {
        d1 = domn->posf->d.at(i) - domn->posf->d.at(i-1);
        if(d1 < dxmin)
            dxmin = d1;
    }
    if(!domn->pram->Lspatial)
        dtCUmax = domn->pram->tdfac * dxmin * dxmin / domn->pram->kvisc0;
    else
        dtCUmax = 50.0 * domn->pram->tdfac * dxmin;

}

///////////////////////////////////////////////////////////////////////////////
/** Diffuse the domain to catch up t0 to time if we have not had eddies for a while.
 *  @param Ldoit  \input Flag with default false
 */

void solver::diffusionCatchUpIfNeeded(bool Ldoit) {

    if(!Ldoit && time-t0 < dtCUmax)
        return;

#ifndef SILENT
    *domn->io->ostrm << endl << "# Catching up diffusion to eddies: time = " << time;
#endif

    domn->mimx->advanceOdt(t0, time);

    t0 = time;

 

}

///////////////////////////////////////////////////////////////////////////////
/**Every once in a while (nDtSmeanWait) test the mean acceptance probability.
 * Increase dtSmean if its too small.
 */

void solver::raiseDtSmean() {

    if(domn->pram->Llem)
        return;

    if(iEtrials % domn->pram->nDtSmeanWait != 0)
        return;                              // only do this once in a while

    if(nPaSum > 0)
        PaSum /= nPaSum;                     // PaSum is now average Pa of eddies

    if(PaSum < domn->pram->Pav) {
       if(PaSum < domn->pram->Pav/domn->pram->dtfac)
           dtSmean *= domn->pram->dtfac;     // increase by at most dtfac
       else
           dtSmean *= domn->pram->Pav/PaSum; // increase dtSmean to target value

       *domn->io->ostrm << endl << "#  increasing dtSmean to " << dtSmean
                        << " (neddies = " << neddies << "  eddy trails = " << iEtrials << ")";
    }

    PaSum  = 0.0;                          // reset values
    nPaSum = 0;
    if (iEtrials > 10000*domn->pram->nDtSmeanWait) {
        *domn->io->ostrm << endl << "#  reset iEtrials, PaSumC, nPaSumC after "
                         << 10000*domn->pram->nDtSmeanWait << " eddy trials.";
        iEtrials = 0;
        PaSumC   = 0.0;
        nPaSumC  = 0;
    }

}

///////////////////////////////////////////////////////////////////////////////
/**Sample an eddy size and position. Fill the eddl from domn.
 * Then triplet map the eddy, compute the eddy timescale, then the
 * acceptance probability.  Roll the dice and if you win (rand # < prob)
 * then accept the eddy.  This means, apply velocity kernels, then
 * insert the eddl into the domn.
 * Note, this function may be better as a member of eddy
 *
 * @return true if the sampled eddy was implemented.
 */

bool solver::sampleEddyAndImplementIfAccepted() {



   
}

///////////////////////////////////////////////////////////////////////////////
/** Sample and implement an LEM eddy
 * @return true if the sampled eddy was implemented.
 */

bool solver::sampleAndImplementLEMeddy() {

   
}

///////////////////////////////////////////////////////////////////////////////
/** Apply a large eddy suppression test based on elapsed time (Echekki 2001)
 *  @param time    \input current time.
 *  @param tauEddy \input eddy timescale, or in spatial cases, the eddy size
 *  Note, this function may be better as a member of eddy
 */

bool solver::testLES_elapsedTime(const double time, const double tauEddy) {

    if(domn->pram->LES_type != "ELAPSEDTIME")
        return true;

    if( (time-domn->pram->x0virtual) < domn->pram->Z_LES * tauEddy ) //doldb
        return false;

    return true;
}

///////////////////////////////////////////////////////////////////////////////
/** Large eddy suppression test on fraction of domain size
 *  @param eSize \input eddy size
 *  Note, this function may be better as a member of eddy
 */
bool solver::testLES_fracDomain(const double eSize) {
    if(domn->pram->LES_type != "FRACDOMAIN")
        return true;
    if(eSize/domn->Ldomain() > domn->pram->Z_LES)
        return false;
    return true;
}

///////////////////////////////////////////////////////////////////////////////
/** Apply a large eddy suppression test based on integral length scale
 *  integral length scale L = L0 * (t/t0)^(1-n/2)
 *  t is elapsed time;
 *  n is between 1.15 and 1.45, usually 1.3
 *  @param time  \input current time.
 *  @param eSize \input eddy size
 *  Guangyuan Sun 06/2013
 *  Note, this function may be better as a member of eddy
 */

bool solver::testLES_integralLength(const double time, const double eSize) {

    if(domn->pram->LES_type != "INTEGRALSCALE")
        return true;

    double n = 1.1;
    double t0 = 0.15899;
    double L0 = 0.028323;
    double integralLength = domn->pram->Z_LES * L0 * pow(time/t0, 1-0.5*n);

    if(eSize > integralLength)
        return false;
    return true;
}

///////////////////////////////////////////////////////////////////////////////
/** Apply the large-eddy suppression test.
 *  Note, this function may be better as a member of eddy
 */

bool solver::testLES_thirds() {

    

    
}

///////////////////////////////////////////////////////////////////////////////
/** Reduce dtSmean if it is resulting in too high an acceptance probability.
 */

void solver::lowerDtSmean() {

  

}

