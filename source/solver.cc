
#include "solver.h"
#include "domain.h"
#include <sstream>
#include <string>
#include <iomanip>
#include <numeric> //accumulate

///////////////////////////////////////////////////////////////////////////////
/** solver initialization function
 *
 * @param p_domn  \input set domain pointer with.
 */


///////////////////////////////////////////////////////////////////////////////

/** destructor
 */


void solver::init(domain *p_domn) {

    domn    = p_domn;}

solver::~solver(){};
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


    //-------------------------------------------------------------------------

    double tLastDA         = 0.0;            ///< time of last diffusive mesh adaption
    int    cLastDA         = 0;              ///< for adaption
    bool   LeddyAccepted   = false;

    stringstream ss1;
    string       s1;

    time     = domn->pram->trst;    // 0.0 unless you are restarting
    t0       = domn->pram->trst;    // 0.0 unless you are restarting

    PaSum    = 0.0;
    nPaSum   = 0;
    neddies  = 0;
    PaSumC   = 0.0;
    nPaSumC  = 0;

    iEtrials = 0;

    //-------------------------------------------------------------------------

    computeDtSmean();
    computeDtCUmax();


    domn->io->writeDataFile("odt_init.dat", time);



    domn->io->writeDataFile("odt_init_adpt.dat", time);

    domn->io->outputHeader();

    //-------------------------------------------------------------------------

    while(time <= domn->pram->tEnd) {

        diffusionCatchUpIfNeeded();


        computeDtCUmax();

        if(domn->pram->Llem)
            LeddyAccepted = sampleAndImplementLEMeddy();
        else
            LeddyAccepted = sampleEddyAndImplementIfAccepted();  // ODT may reduce dtSmean; sets Pa

        iEtrials++;

        //----------------------

        if(LeddyAccepted) {

            if(++neddies % (domn->pram->modDisp*50) == 0)
                domn->io->outputHeader();
            if(neddies   % domn->pram->modDisp == 0)
                domn->io->outputProgress();

            //if (neddies % domn->pram->modDump == 0) {
            //    ss1.clear();  ss1 << setfill('0') << setw(4) << neddies; ss1 >> s1;
            //    domn->io->writeDataFile("odt_"+s1+"_eddy.dat", time);
            //}

            

            //if (neddies % domn->pram->modDump == 0) {
            //    ss1.clear();  ss1 << setfill('0') << setw(4) << neddies; ss1 >> s1;
            //    domn->io->writeDataFile("odt_"+s1+"_adptEd.dat", time);
            //}

            diffusionCatchUpIfNeeded(true);

            //if (neddies % domn->pram->modDump == 0) {
            //    ss1.clear();  ss1 << setfill('0') << setw(4) << neddies; ss1 >> s1;
            //    domn->io->writeDataFile("odt_"+s1+"_diffuse.dat", time);
            //}



            if (neddies % domn->pram->modDump == 0) {
                ss1.clear();  ss1 << setfill('0') << setw(4) << neddies; ss1 >> s1;
                domn->io->writeDataFile("odt_"+s1+"_adptDif.dat", time);
            }

        }

        //----------------------

        time += sampleDt();             // advance the time

        raiseDtSmean();                 // may reset PaSum, nPaSum, dtSmean
    }

    time = domn->pram->tEnd;
    if(t0 < time)
        diffusionCatchUpIfNeeded(true);

    //-------------------------------------------------------------------------


    domn->io->writeDataFile("odt_end.dat", time);
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

