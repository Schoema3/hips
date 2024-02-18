#include "hips.h"
#include "batchReactor_cvode.h"
#include "batchReactor_cantera.h"
#include "randomGenerator.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
/**
 * \brief Constructor for initializing 'nLevels_' based on Reynolds number.
 * \param Re                Reynolds number for turbulence simulation.
 */
hips::hips(double Re) {

    double original_N = (3.0/4) * log(1/Re) / log(Afac);                            // Calculate the original N

    double closest_N = round(original_N);                                           // Round the original N to the closest integer

    double desired_exponent = ceil(log(1/Re) / log(Afac));                          // Calculate the desired exponent that makes log_A(1/Re) an integer

    Afac = pow(1/Re, 1/desired_exponent);                                           // Calculate the adjusted A

    nL = (3.0/4) * log(1/Re) / log(Afac) + 3;                                       // Recalculate N

    if (abs(nL - original_N) > abs(closest_N - original_N)) {                        // Check if adjusting A to a lower exponent brings N closer to the original N
        desired_exponent -= 1;
        Afac = pow(1/Re, 1/desired_exponent);
        nL = (3.0/4) * log(1/Re) / log(Afac);
        nL = nL + 3;
    }

      
}

///////////////////////////////////////////////////////////////////////////////
/**
 * \Constructor for initializing required parameters to create the HiPS tree. 
 * \param nLevels_                     Input number of tree levels.
 * \param domainLength_                Input length scale of the domain.
 * \param tau0_                        Input time scale of the domain.
 * \param C_param_                     Input c parameter to control eddy rate
 * \param forceTurb_                   Input flag for forcing turbulence.
 * \param nVar_                        Input number of variables.
 * \param ScHips_                      Input vector of Schmidt numbers for HiPS simulation.
 * \param cantSol                      Input Cantera solution object.
 * \param performReaction_             Input flag for performing chemical reactions
 * \param seed                         Input seed for the random number generator(negative to randomize it).
 */
hips::hips(int nLevels_, 
           double domainLength_, 
           double tau0_, 
           double C_param_, 
           int forceTurb_,
           int nVar_,
           vector<double> &ScHips_,
           shared_ptr<Cantera::Solution> cantSol,
           bool performReaction_,
           int seed) : 

    nLevels(nLevels_), 
    domainLength(domainLength_), 
    tau0(tau0_),
    C_param(C_param_), 
    forceTurb(forceTurb_),       
    ScHips(ScHips_),   
    nVar(nVar_),                       
    LrandSet(true),              
    rand(seed),
    performReaction(performReaction_) { 
    
    gas = cantSol->thermo(); 
    nsp = gas->nSpecies();

    varData.resize(nVar);
    varRho.resize(nparcels);
    varName.resize(nVar);

    bRxr = make_unique<batchReactor_cvode>(cantSol);                  // initialize integrator pointer(cvode)
    //bRxr = make_unique<batchReactor_cantera>(cantSol);                // initialize integrator pointer(cantera)
              
    
    //-------------------------- Set number of parcels, level lengthscales, timescales, and rates, and i_plus && i_batchelor
     
    iEta = nLevels - 3;        // Kolmogorov level; if nLevels = 7, then 0, 1, 2, 3, (4), 5, 6; iEta=4 is the lowest swap level: swap grandchildren of iEta=4 at level 6.
           
    int maxSc = 1.0;

    for (int i=0; i<ScHips.size(); i++)
        maxSc = ScHips[i]>maxSc ? ScHips[i] : maxSc;
    
    if (maxSc > 1.0)
        nLevels += ceil(log(maxSc)/log(4));       // Changing number of levels!
    
    Nm1 = nLevels - 1;
    Nm2 = nLevels - 2;
    Nm3 = nLevels - 3;
    
    // -------------------------- 
    
    nparcels = static_cast<int>(pow(2, Nm1));
    parcelTimes.resize(nparcels,0);
    i_batchelor.resize(nVar,0);
    
    vector<double> levelLengths(nLevels);      // including all levels, but last 2 don't count:
    vector<double> levelTaus(nLevels);         // smallest scale is 2 levels up from bottom
    levelRates   = vector<double>(nLevels);

    for (int i=0; i<nLevels; i++) {
        levelLengths[i] = domainLength * pow(Afac,i);
        levelTaus[i] = tau0 * pow(levelLengths[i]/domainLength, 2.0/3.0) / C_param;
        levelRates[i] = 1.0/levelTaus[i] * pow(2.0,i);
    }
    
    LScHips = ScHips.size() > 0 ? true : false;
    if (LScHips) {                             // correct levels for high Sc (levels > Kolmogorov)
        for (int i=iEta+1; i<nLevels; i++) {
            levelTaus[i] = tau0 *
            pow(levelLengths[iEta]/domainLength, 2.0/3.0) /
            C_param;
            levelRates[i] = 1.0/levelTaus[i] * pow(2.0,i);
        }
    }

    //-------------------------------------------------

    eddyRate_total = 0.0;
    for (int i=0; i<=Nm3; i++)
        eddyRate_total += levelRates[i];
    
    eddyRate_inertial = 0.0;
    for (int i=0; i<=iEta; i++)
        eddyRate_inertial += levelRates[i];
    
    //-------------------
    
    i_plus.resize(nVar);
    for (int k=0; k<nVar; k++) {
        if (ScHips[k] < 1.0)
            i_batchelor[k] = iEta + 1.5*log(ScHips[k])/log(4);
        else if (ScHips[k] > 1.0)
            i_batchelor[k] = iEta +     log(ScHips[k])/log(4);
        else
            i_batchelor[k] = iEta;

        i_plus[k] = ceil(i_batchelor[k]);
    }
    
    //------------------- Set the parcel addresses (index array)
    
    pLoc.resize(nparcels);
    for (int i=0; i<nparcels; i++)
        pLoc[i] = i;
} 
////////////////////////////////////////i///////////////////////////////////////
/**
 * \brief Function to pass the variables, their weights, and their names to the parcels of the tree.  
 * \param v                            Input vector consisting of variables passed to the HiPS tree.
 * \param w                            Input vector containing weights for each flow particle.
 * \param varN                         Input vector containing names corresponding to each variable.
 * \param i                             Index indicating where the data should be assigned in the vectors varData and varName.
 */
void hips::set_varData(std::vector<double> &v, std::vector<double> &w, const std::string &varN, int i) {  

    varData[i] = new vector<double>(projection(v, w));
    varName[i] = varN;
}
///////////////////////////////////////i///////////////////////////////////////
/**
 * \brief Function to pass the variable, their weights, names, and densities to the parcels of the tree. 
 * \param v                              Vector of variables that are passed to the HiPS tree.
 * \param w                              Vector of weights; each flow particle has a weight.
 * \param varN                           Vector of names of the variable.
 * \param rho                            Vector of density; each flow particle has a specific density.
 * \param i                              index indicating where the data should be assigned in the vectors varData and varName.
 * \note This function is overloaded. This version considers particle density.
 */
void hips::set_varData(std::vector<double> &v, std::vector<double> &w, const std::string &varN, const std::vector<double> &rho, int i) {
    
    std::pair<std::vector<double>, std::vector<double>> results = projection(v, w, rho);
    std::vector<double> vh = results.first;
    std::vector<double> rho_h = results.second;

    varData[i] = new std::vector<double>(vh);            
    varRho =  std::vector<double>(rho_h);                 
}

/////////////////////////////////////////i/////////////////////////////////////////////
/**
 * \brief Function for projecting vectors onto a grid.
 * \param v                              Vector to be projected.
 * \param w                              Weight vector.
 * \return Projected vector.
 */
std::vector<double> hips::projection(std::vector<double> &vcfd, std::vector<double> &weight) {
    
    std::vector<double> xc = setGridCfd(weight);
    std::vector<double> xh = setGridHips(nparcels);

    int nc = xc.size() - 1;                 
    int nh = xh.size() - 1; 

    std::vector<double> vh(nh, 0.0);
    int jprev = 0;

    for(int i = 0; i < nh; i++) {
        for(int j = jprev + 1; j <= nc; ++j) {
            if(xc[j] <= xh[i + 1]) {
                double d1 = xc[j] - xc[j - 1];
                double d2 = xc[j] - xh[i];
                double d = std::min(d1, d2);

                vh[i] += vcfd[j - 1] * d;
            } 

        else {
                double d1 = xh[i + 1] - xc[j - 1];
                double d2 = xh[i + 1] - xh[i];
                double d = std::min(d1, d2);

                vh[i] += vcfd[j - 1] * d;
                jprev = j - 1;
                break;
            }
        }
        
        vh[i] /= (xh[i + 1] - xh[i]);           
    }
    return vh;
}

/////////////////////////////////////////i/////////////////////////////////////////////
/**
 * \brief Project vectors onto a grid.
 *
 * This function projects vectors onto a grid. It takes a vector of variables, a weight vector, and a vector of density as input parameters,
 * and returns the projected vector.
 *
 * \param vcfd                          Vector of variables passed to the HiPS tree.
 * \param weight                        Weight vector; each flow particle has a weight.
 * \param density                       Vector of density.
 * \return A pair of vectors representing the projected vector and the density on the grid.
  * \note This function is overloaded. This version considers particle density.
 */
std::pair<std::vector<double>, std::vector<double>> hips::projection(std::vector<double> &vcfd, std::vector<double> &weight, const std::vector<double> &density) {
    
    std::vector<double> xc = setGridCfd(weight);
    std::vector<double> xh = setGridHips(nparcels);

    // Initialize variables
    int nc = xc.size() - 1;
    int nh = xh.size() - 1;

    std::vector<double> vh(nh, 0.0);
    std::vector<double> rho_h(nh, 0.0);
    int jprev = 0;

    // Main loop
    for (int i = 0; i < nh; i++) {
        double total_dx = 0.0;

        for (int j = jprev + 1; j <= nc; ++j) {
            if (xc[j] <= xh[i + 1]) {
                double d = std::min(xc[j] - xc[j - 1], xc[j] - xh[i]);
                total_dx += d;
                rho_h[i] += density[j - 1] * d;
                vh[i] += density[j - 1] * vcfd[j - 1] * d;
            } 

        else {
                double d = std::min(xh[i + 1] - xc[j - 1], xh[i + 1] - xh[i]);
                total_dx += d;
                rho_h[i] += density[j - 1] * d;
                vh[i] += density[j - 1] * vcfd[j - 1] * d;
                jprev = j - 1;
                break;
            }
        }

        // Normalize results
        rho_h[i] /= (xh[i + 1] - xh[i]);
        vh[i] /= rho_h[i] * (xh[i + 1] - xh[i]);
    }

    return {vh, rho_h};
}


/////////////////////////////////////////////////////////////////////////////////
/**
 * \brief Set up a grid for CFD simulation.
 *
 * This function generates a grid for CFD simulation based on a given weight vector.
 *
 * \param w                     Weight vector defining the spacing between grid points.
 * \return Vector representing the grid positions.
 */
std::vector<double> hips::setGridCfd(std::vector<double> &w) {
    
    std::vector<double> pos;                               // Initializing a vector to hold the grid positions
    double posL = 0.0;                                     // Initializing the starting position

    int i = 0;

    while (i <= w.size()) {                                // Generate the grid positions based on the weights
        pos.push_back(posL);                               // Add the current position to the grid
        posL += w[i];                                      // Move to the next position by adding the corresponding weight
        i++;                                              
                                                   
    }

    return pos;                                          // Return the generated grid positions
}

///////////////////////////////////////////////////////////////////////////////
/**
 * \brief Set up a grid for HIPS simulation.
 *
 * This function generates a grid for HIPS simulation with a specified number of grid points.
 *
 * \param N                       Number of grid points.
 * \return Vector representing the grid.
 */
std::vector<double> hips::setGridHips(int N){

    std::vector<double> xh(N + 1);                               // Initialize a vector to hold the grid points
    double step = 1.0 / N;                                       // Calculate the step size

    for(int i = 0; i <= N; i++){                                // Populate the grid with evenly spaced points
        xh[i] = i * step;
    }
        
    return xh;                                                  // Return the generated grid
}

///////////////////////////////////////////////////////////////////////////////
/** The HiPS solver
 * \param tRun                               Input simulation run time
 * \param shouldWriteData                    Set to false by default. If true, data will be written.
 *
 * Sample tee (time of next eddy event)
 * Select and swap subtrees (at current time, not at tee, so that we know who's involved)
 * If the eddy event is at the parcel/micromixing level:
 *     React involved parcels from their current time to tee
 *     Mix the involved parcels (micromixing)
 */

void hips::calculateSolution(const double tRun, bool shouldWriteData) {

    unsigned long long nEddies = 0;                                        // number of eddy events
    int    fileCounter = 0;                                                // number of data files written
    int    iLevel;                                                         // tree level of EE with top at iLevel=0
    int    iTree;                                                          // one of two subtrees involved in swap at iLevel
    double dtEE;                                                           // time increment to next eddy event 
    time = 0.0;                        // initialize simulation time

    sample_hips_eddy(dtEE, iLevel);    // get first EE at time 0+dtEE
    nEddies++;

    while (time+dtEE<=tRun) {
        time += dtEE;
        selectAndSwapTwoSubtrees(iLevel, iTree);
        advanceHips(iLevel, iTree);    // reaction and micromixing (if needed) to t=time

        sample_hips_eddy(dtEE, iLevel);

        nEddies++;
        //cout<<"-------------------------------------"<<nEddies<<endl;
        if(shouldWriteData && nEddies %200000 == 0) writeData(++fileCounter, time);
    }
    time = tRun;
    iLevel = 0; iTree  = 0;
    if(performReaction)
    reactParcels_LevelTree(iLevel, iTree);      // react all parcels up to end time
}


///////////////////////////////////////////////////////////////////////////////
/** Samples stochastic eddy events on the hips tree: time and level.
*   \param dtEE                            Output time increment to next eddy event (EE)
*   \param iLevel                          Output tree level of EE
*/

void hips::sample_hips_eddy(double &dtEE, int &iLevel) {

    static double c1 = 1.0 - pow(2.0, 5.0/3.0*(iEta+1));
    static double c2 = pow(2.0, Nm2) - pow(2.0, iEta+1);
    static double c3 = pow(2.0, iEta+1);

    //--------------- time to next eddy

    double r = rand.getRand();
    dtEE = -log(r)/eddyRate_total;

    //----------------- get eddy level

    r = rand.getRand();

    if ( r <= eddyRate_inertial/eddyRate_total) {     // inertial region
        r = rand.getRand();
        iLevel = ceil(3.0/5.0*log2(1.0-r*c1) - 1.0);
        if (iLevel < 0)    iLevel = 0;
        if (iLevel > iEta) iLevel = iEta;
    }

    else {                                            // "Batchelor" region
        r = rand.getRand();
        iLevel = ceil(log2(r*c2 + c3) - 1.0);
        if(iLevel < iEta+1) iLevel = iEta+1;
        if(iLevel > Nm3) iLevel = Nm3;
    }

    return;
}

///////////////////////////////////////////////////////////////////////////////
/** Function performs eddy events: parcel swaps.
    @param iLevel                          Input  level of the tree for the base of the swap.
    @param iTree                           Output which subtree on the level is selected

    Randomly select a node on iLevel.
    Go down two levels and select nodes 0q and 1r, where q, r are randomly 0 or 1
    Find the starting index of the Q-tree and R-tree to swap and the number of parcels.
    Then swap the cells.

    For a 6 level tree: 0, 1, 2, 3, 4, 5:
    If iLevel = 1, then suppose i=1, 0q = 00 and 1r = 11:
    Then we are swaping 0100** with 0111** or (01|00|**) with (01|11|**)
       or i0qs with i1rs, where i = 01; 0q = 00; 1r = 11; and s = **

    We use bitwise shifts for easy powers of 2.
    The swap is done by adding or subtracting a value (shift),
        which should be equivalent to flipping the swapping the two 0q bits and 1r bits.
                                                                                                              Level
                                                                                                            ---------
                                                    *                                                           0
                                                 /     \
                                              /           \
                                           /                 \
                                        /                       \
                                     /                             \
                                  /                                   \
                               /                                         \
                            /                                               \
                           *                                                (*)  01|0000                        1
                          / \                                               / \
                        /     \                                           /     \
                      /         \                                       /         \
                    /             \                                   /             \
                  /                 \                               /                 \
                /                     \                           /                     \
               *                       *                         *                       *                      2
              / \                     / \                       / \                     / \
            /     \                 /     \                   /     \                 /     \
          /         \             /         \               /         \             /         \
         *           *           *           *            [*] 00|**    *           *          [*] 11|**         3
        / \         / \         / \         / \           / \         / \         / \         / \
       /   \       /   \       /   \       /   \         /   \       /   \       /   \       /   \
      *     *     *     *     *     *     *     *       *     *     *     *     *     *     *     *             4
     / \   / \   / \   / \   / \   / \   / \   / \     / \   / \   / \   / \   / \   / \   / \   / \
    00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15   16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31           5
                                                      ^^^^^^^^^^^                         ^^^^^^^^^^^

*/

void hips::selectAndSwapTwoSubtrees(const int iLevel, int &iTree) {

    iTree = rand.getRandInt((1 << iLevel)-1);
    int zero_q = rand.getRandInt(1);                                    // 0q where q is 0 or 1
    int one_r  = 2 + rand.getRandInt(1);                                // 1r where r is 0 or 1

    int Qstart = (zero_q << (Nm3-iLevel)) + (iTree << (Nm1-iLevel));     // starting index of Q parcels
    int Rstart = (one_r  << (Nm3-iLevel)) + (iTree << (Nm1-iLevel));     // starting index of R parcels
    int nPswap = 1 << (Nm3-iLevel);                                      // number of parcels that will be swapped

    int Qend = Qstart + nPswap;                                          // inclusive indices are Qstart to Qend-1
    int Rend = Rstart + nPswap;                                          // inclusive indices are Rstart to Rend-1
    vector<int> aa(pLoc.begin()+Qstart, pLoc.begin()+Qend);
    copy(pLoc.begin()+Rstart, pLoc.begin()+Rend, pLoc.begin()+Qstart);   // python: pLoc[Qstart:Qend]=pLoc[Rstart:Rend]
    copy(aa.begin(), aa.end(), pLoc.begin()+Rstart);                     // python: pLoc[Rstart:Rend]=aa
}


///////////////////////////////////////////////////////////////////////////////
/** React and mix parcels that are involved in a micromixing process.
  * This is determined by the level and the tree within that level.
  * Might not do anything if the eddy doesn't cause micromixing.
  * \param iLevel                    Input level that the eddy event occurred.
  * \param iTree                     Input root note of the eddy event at iLevel.
  */

void hips::advanceHips(const int iLevel, const int iTree) {
     
    if (forceTurb==2 && iLevel==0)              // forcing for statistically stationary
        forceProfile();

    bool rxnDone = false;                       // react all variables once
    for (int k=0; k<nVar; k++) {                // upon finding first variable needing micromixing
        if ( (iLevel >= i_plus[k]) || 
             (iLevel==i_plus[k]-1 && rand.getRand() <= i_plus[k]-i_batchelor[k]) ) {
                if(!rxnDone && performReaction) {
                   reactParcels_LevelTree(iLevel, iTree);
                   rxnDone = true;
             }
            mixAcrossLevelTree(k, iLevel, iTree);
        }
    }

}
///////////////////////////////////////////////////////////////////////////////
/** React parcels that are involved in a micromixing process.
  * Parcels might react for different amounts of time depending on when they last
  * reacted, which is stored in the parcelTimes array.
  * @param iLevel                    Input level that the eddy event occurred.
  * @param iTree                     Input root note of the eddy event at iLevel.
  */

void hips::reactParcels_LevelTree(const int iLevel, const int iTree) {
    
    int nP  = 1 << (Nm1 - iLevel);
    //int istart = iTree << nP;
    int istart = iTree * nP;
    int iend = istart + nP;
    int ime;
    double dt;
    double h;
    vector<double> y(nsp);

    for (int i=istart; i<iend; i++) {
        ime = pLoc[i];
        //cout<<"pLoc "<<pLoc[i]<<"  gas "<<gas->density()<<"\n\n"<<endl;
        dt = time-parcelTimes[ime];
        h = varData[0][0][ime]; 
        for (int k=0; k<nsp; k++)
            y[k] = varData[k+1][0][ime];

        bRxr->react(h, y, dt);
        varRho[ime] =  bRxr->getDensity();             //Mb
         //cout<<"h    "<<h<<endl;
        varData[0][0][ime] = h;
        for (int k=0; k<nsp; k++)
            varData[k+1][0][ime] = y[k];
    }

}

///////////////////////////////////////////////////////////////////////////////
/** Mix parcels uniformly (using average) at given level and tree
  * Mixing at levels above the lowest enables low Sc variables.
  * @param kVar   \input  variable index to mix (normally a transported var as determined by caller)
  * @param iLevel \input  grandchildren of this iLevel will be mixed
  * @param iTree  \input  at the given iLevel, mix only this subtree
  *
  * Example: For a 5 level tree, we have levels 0, 1, 2, 3, 4 (top to bottom).
  *   
  *          Let iLevel = 2, then nPmix = 2 and we are mixing pairs of parcels at the base of the tree.
  *          iTree can be 0, 1, 2, or 3.
  *          if iTree = 0 we will mix parcels (0,1) and (2,3) as the left subtree and right subtree:
  *                (istart=1, iend=2) and (istart=2, iend=4)
  *          if iTree = 1 we will mix parcels (4,5) and (6,7) with (istart=4,iend=6), (istart=6, iend=8)
  *          if iTree = 2 we will mix parcels (8,9) and (10,11) with (istart=8,iend=10), (istart=10, iend=12)
  *          if iTree = 3 we will mix parcels (12,13) and (14,15) with (istart=12,iend=14), (istart=14, iend=16)
  *         
  *          Let iLevel=1 then we will be mixing groups of four parcels.
  *          iTree can be 0 or 1
  *          if iTree = 0 we will mix parcels (0,1,2,3) and (4,5,6,7) with (istart=0,iend=4), (istart=4, iend=8)
  *          if iTree = 1 we will mix parcels (8,9,10,11) and (12,13,14,15) with (istart=8,iend=12), (istart=12, iend=16)
  *
  * recall: 3 << 4 means 3*2^4 (or 3 = 000011 and 3<<4 = 110000 = 48), that is, we shift the bits left 4 places.
  *          
  * NOTE: BE CAREFUL WITH MIXING SOME SCALARS, LIKE MASS FRACTIONS; CURRENT CODE ASSUMES ALL PARCELS HAVE SAME DENSITY (mixing Yi directly)
  */

void hips::mixAcrossLevelTree(int kVar, const int iLevel, const int iTree) {
    
    int istart;
    int iend;

    int nPmix = 1 << (nLevels - iLevel - 2);                   // number of parcels mixed together

    int ime;

    //---------- Mix left branch of iTree

    istart = iTree << (Nm1-iLevel);  
    iend = istart + nPmix;

    double s = 0;                                               // initialize sum to 0
    for (int i=istart; i<iend; i++) {
        ime = pLoc[i];
        s += varData[kVar][0][ime];
    }
    for (int i=istart; i<iend; i++) {
        ime = pLoc[i];
        varData[kVar][0][ime] = s / nPmix; 
    }

    //--------- Mix right branch of iTree

    istart = iend;
    iend = istart + nPmix;

    s = 0;                   // initialize sum to 0
    for (int i=istart; i<iend; i++) {
        ime = pLoc[i];
        s += varData[kVar][0][ime];
    }
    for (int i=istart; i<iend; i++) {
        ime = pLoc[i];
        varData[kVar][0][ime] = s / nPmix; 

        //cout<<"var Data in mixxx "<<varData[0][0][ime]<<endl;
        
    }   
}

///////////////////////////////////////////////////////////////////////////////
/**
 * @brief Force the HiPS profile to achieve statistical stationarity.
 *
 * This function is meant to demonstrate forcing for simple scalars, such as a mixture fraction variable that varies between 0 and 1.
 * The code within this function is configured to force the left half of parcels to average 0 and the right half of parcels to average 1.
 * 
 * @note This function modifies the data stored in the HiPS profile.
 */
void hips::forceProfile() {
    // Loop through each variable in the HiPS profile
    for (int k = 0; k < varData.size(); k++) {
        double s;                                                    // Temporary variable for summation

        //---------- Force the left half of parcels to average 0 ----------

        for (int i = 0; i < nparcels >> 1; i++)
            s += varData[k][0][pLoc[i]];                             // Calculate the sum of values in the left half of parcels
        
        s /= (nparcels >> 1); // Calculate the average of values in the left half of parcels
        
        for (int i = 0; i < nparcels >> 1; i++)
            varData[k][0][pLoc[i]] += (-s - 0.0);                    // Adjust values in the left half of parcels to achieve an average of 0

        //---------- Force the right half of parcels to average 1 ----------
        s = 0.0;

        for (int i = nparcels >> 1; i < nparcels; i++)
            s += varData[k][0][pLoc[i]];                             // Calculate the sum of values in the right half of parcels
        
        s /= (nparcels >> 1);                                        // Calculate the average of values in the right half of parcels
        
        for (int i = nparcels >> 1; i < nparcels; i++)
            varData[k][0][pLoc[i]] += (-s + 1.0);                    // Adjust values in the right half of parcels to achieve an average of 1
    }
}

///////////////////////////////////////////////////////////////////////////////
/**
 * @brief Write data to a file with a specific index and output time.
 *
 * This function writes data to a file with a filename following a sequential incrementing index pattern 
 * (e.g., Data_00001.dat, Data_00002.dat, etc.). The time of the data file is also written into the file.
 *
 * @param ifile The sequential index used in the file name.
 * @param outputTime The time associated with the data, written into the file.
 */
void hips::writeData(const int ifile, const double outputTime) {
    // Create the "data" directory if it doesn't exist
    system("mkdir -p ../data");

    stringstream ss1;
    string s1;
    ss1 << setfill('0') << setw(5) << ifile;
    ss1 >> s1;                                                                  // Convert the index into a string with leading zeros for filename

    string fname = "../data/Data_" + s1 + ".dat";                               // Construct the filename

    ofstream ofile(fname.c_str());                                              // Open the file for writing

    ofile << "# time = " << outputTime << "\n";                                  // Write the time and number of grid points into the file
    ofile << "# Grid Points = " << nparcels << "\n";

    for (const auto& name : varName)
        ofile << setw(14) << name;                                              // Write variable names into the file
    ofile << scientific << setprecision(10);

    for (int i = 0; i < nparcels; i++) {
        ofile << "\n";                                                          // Start a new line for each parcel

        for (int k = 0; k < nVar; k++)
            ofile << setw(19) << varData[k][0][pLoc[i]];                       // Write data for each variable
    }

    ofile.close();                                                              // Close the file
}

//void hips::writeDataiconst int ifile, const double outputTime) {
//    
//    stringstream ss1;
//    string s1;
//    ss1 << setfill('0') << setw(5) << ifile;
//    ss1 >> s1;                                                                  // Convert the index into a string with leading zeros for filename
//
//    string fname = "../data/Data_" + s1 + ".dat";                               // Construct the filename
//
//    ofstream ofile(fname.c_str());                                              // Open the file for writing
//
//    ofile << "# time = " << outputTime << "\n";                                  // Write the time and number of grid points into the file
//    ofile << "# Grid Points = " << nparcels << "\n";
//
//    for (const auto& name : varName)
//        ofile << setw(14) << name;                                              // Write variable names into the file
//    ofile << scientific << setprecision(10);
//
//    for (int i = 0; i < nparcels; i++) {
//        ofile << "\n";                                                          // Start a new line for each parcel
//
//        for (int k = 0; k < nVar; k++)
//            ofile << setw(19) << varData[k][0][pLoc[i]];                       // Write data for each variable
//    }
//
//    ofile.close();                                                              // Close the file
//}
///////////////////////////////////////////////////////////////////////////////
/**
 * \brief Function for projecting vectors onto a grid.
 * \param vb Vector to be projected back.
 */
std::vector<double> hips::projection_back(std::vector<double> &vb){
    
    int jprev = 0;

    std::vector<double> xh = setGridHips(nparcels);
    std::vector<double> xc = setGridCfd(vb);

    int nh = xh.size() - 1;
    int nc = xc.size() - 1;
    std::vector<double> vc(nc, 0);

    for (int i = 0; i < nc; ++i) {
        for (int j = jprev + 1; j < xh.size(); ++j) {
            double d1 = xh[j] - xh[j - 1];
            double d2 = xh[j] - xc[i];
            double d = (d1 < d2) ? d1 : d2;

            if (xh[j] <= xc[i + 1]) {
                vc[i] += vb[j - 1] * d;
            } else {
                double d1 = xc[i + 1] - xh[j - 1];
                double d2 = xc[i + 1] - xc[i];
                d = (d1 < d2) ? d1 : d2;
                vc[i] += vb[j - 1] * d;
                jprev = j - 1;
                break;
            }
        }
        vc[i] /= (xc[i + 1] - xc[i]);
    }

    return vc;
}    

/////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Retrieve modified data from the Hierarchical Progressive Survey (HiPS) library.
 *
 * This function retrieves modified data from the HiPS library and stores it in the provided vector.
 * Each element of the returned vector contains a vector representing a single data projection.
 *
 * @return A vector of vectors containing modified data retrieved from the HiPS library.
 */
std::vector<std::vector<double>> hips::get_varData(){

    std::vector<std::vector<double>> varDataProjections;                             // Vector to store modified data projections

    for (int i = 0; i < varData.size(); i++) {                                        // Loop through each element of varData and project the data back
        varDataProjections.push_back(projection_back(varData[i][0]));                 // Project the data back and store it in vh
    }

    return varDataProjections;                                                         // Return the vector containing modified data projections
}

////////////////////////////////////////////////////////
 


