#include "hips.h"

#ifdef REACTIONS_ENABLED
    #include "batchReactor_cvode.h"
    #include "batchReactor_cantera.h"
#endif

#include "randomGenerator.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
using namespace std;

////////////////////////////////////////////////////////////////////////////////

int hips::nL = 0;                   // Initialize adjusted number of levels based on the Reynolds number.
double hips::Prob = 0.0;            // Initialize probability value for probability-based solution.
double hips::lStar = 0.0;           // Initialize length of the level associated with the Reynolds number.
double hips::Anew = 0.0;            // Initialize adjusted level lengthscale reduction factor.          

//////////////////////////////////////////////////////////////////////////////

/// @brief Constructor to initialize parameters for creating the HiPS tree.
///
/// This constructor initializes parameters including, eddy rate, forcing turbulence, number of variables, solution obj, etc. 
/// It is particularly useful for initializing the tree multiple times, as it allows for the parameters to be set once and reused. This constructor is typically followed by a call to `set_tree(...)`.
///
/// \param C_param_         Parameter controlling eddy rate.
/// \param forceTurb_       Flag indicating whether to force turbulence.
/// \param nVar_            Number of variables.
/// \param cantSol          Cantera solution object.
/// \param performReaction_ Flag indicating whether to perform chemical reactions.
/// \param seed             Seed for the random number generator (negative value to randomize it).
///
/// \note `bRxr` is a pointer to the integrator object. By default, `batchReactor_cvode` is enabled. To switch to `batchReactor_cantera`, the user needs to uncomment the corresponding section.

//////////////////////////////////////////////////////////////////////////////

hips::hips(double C_param_, 
           int forceTurb_,
           int nVar_,
           bool performReaction_,
         #ifdef REACTIONS_ENABLED
                 shared_ptr<Cantera::Solution> cantSol,
          #endif
         int seed) : 

    C_param(C_param_), 
    forceTurb(forceTurb_),       
    nVar(nVar_),                       
    LrandSet(true),              
    rand(seed),
    performReaction(performReaction_) {
      
    #ifdef REACTIONS_ENABLED
        gas = cantSol->thermo(); 
        nsp = gas->nSpecies();

        // By default, use batchReactor_cvode
        bRxr = make_unique<batchReactor_cvode>(cantSol);

        // Uncomment the following line to switch to batchReactor_cantera
        // bRxr = make_unique<batchReactor_cantera>(cantSol);
    #endif

    // Resize vectors to the number of variables
    varData.resize(nVar);
    varName.resize(nVar);        
} 

/////////////////////////////////////////////////////////////////////////////////////////////////////////

/// @brief Constructor to initialize parameters for creating the HiPS tree.
///
/// This constructor initializes all necessary parameters in a single step, including length scale, time scale, eddy rate, and more. 
/// It is beneficial for scenarios where users need to initialize the tree once, such as for simple mixing simulations.
///
/// \param nLevels_         Number of tree levels.
/// \param domainLength_    Length scale of the domain.
/// \param tau0_            Time scale of the domain.
/// \param C_param_         Parameter to control eddy rate.
/// \param forceTurb_       Flag indicating whether to force turbulence.
/// \param nVar_            Number of variables.
/// \param ScHips_          Vector of Schmidt numbers for HiPS simulation.
/// \param cantSol          Cantera solution object.
/// \param performReaction_ Flag indicating whether to perform chemical reactions.
/// \param seed             Seed for the random number generator (negative value to randomize it).
///
/// \note `bRxr` is a pointer to the integrator object. By default, `batchReactor_cvode` is enabled. To switch to `batchReactor_cantera`, the user needs to uncomment the corresponding section.

/////////////////////////////////////////////////////////////////////////////////////////////////////

hips::hips(int nLevels_, 
           double domainLength_, 
           double tau0_, 
           double C_param_, 
           int forceTurb_,
           int nVar_,
           vector<double> &ScHips_,
           bool performReaction_,
         #ifdef REACTIONS_ENABLED
                 shared_ptr<Cantera::Solution> cantSol,
          #endif
         int seed): 

   nLevels(nLevels_), 
   domainLength(domainLength_), 
   tau0(tau0_),
   C_param(C_param_), 
   forceTurb(forceTurb_),       
   ScHips(ScHips_),   
   nVar(nVar_),                       
   LrandSet(true),              
   rand(seed),
   performReaction(performReaction_){

   #ifdef REACTIONS_ENABLED
        gas = cantSol->thermo(); 
        nsp = gas->nSpecies();

        // By default, use batchReactor_cvode
        bRxr = make_unique<batchReactor_cvode>(cantSol);

        // Uncomment the following line to switch to batchReactor_cantera
        // bRxr = make_unique<batchReactor_cantera>(cantSol);
    #endif

    // Resize vectors to the number of variables
    varData.resize(nVar);
    varName.resize(nVar); 

    set_tree(nLevels, domainLength, tau0, ScHips);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

/// \brief Function to create a tree structure based on the specified parameters.
/// \param nLevels_         Number of levels in the tree.
/// \param domainLength_    Length scale of the domain.
/// \param tau0_            Time scale of the domain.
/// \param ScHips_          Vector of Schmidt numbers for HiPS simulation.
/// \note This function sets up the tree based on the specified number of levels. It is useful when the user knows the number of levels explicitly.

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void hips::set_tree(int nLevels_, double domainLength_, double tau0_, vector<double> &ScHips_){    
 
    nLevels= nLevels_; 
    domainLength = domainLength_; 
    tau0 = tau0_; 
    ScHips = ScHips_; 

 
    if (nLevels == -1)  
        nLevels = nL; 

    iEta = nLevels - 3;                            // Kolmogorov level; if nLevels = 7, then 0, 1, 2, 3, (4), 5, 6; iEta=4 is the lowest swap level: swap grandchildren of iEta=4 at level 6.
           
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
    
    vector<double> levelLengths(nLevels);            // Including all levels, but last 2 don't count:
    vector<double> levelTaus(nLevels);               // Smallest scale is 2 levels up from bottom
    levelRates   = vector<double>(nLevels);

    for (int i=0; i<nLevels; i++) {
        levelLengths[i] = domainLength * pow(Afac,i);
       //levelLengths[i] = domainLength * pow(Anew,i);

        levelTaus[i] = tau0 * pow(levelLengths[i]/domainLength, 2.0/3.0) / C_param;
        levelRates[i] = 1.0/levelTaus[i] * pow(2.0,i);

    }

    LScHips = ScHips.size() > 0 ? true : false;
    if (LScHips) {                                     // Ccorrect levels for high Sc (levels > Kolmogorov)
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
            i_batchelor[k] = iEta + log(ScHips[k])/log(4);
        else
            i_batchelor[k] = iEta;

        i_plus[k] = ceil(i_batchelor[k]);
    }
    
    //------------------- Set the parcel addresses (index array)
   varRho.resize(nparcels);
   Temp.resize(nparcels);
 
    pLoc.resize(nparcels);
    for (int i=0; i<nparcels; i++)
        pLoc[i] = i;
} 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// \brief Function to create a tree structure based on the specified parameters.
/// \param Re_ The Reynolds number.
/// \param approach The method used for setting the number of levels based on the Reynolds number:
///                 - approach "1": This approach is introduced as Rounding to closest level to \f$ i_s^* \f$ for micromixing in the paper.
///                 - approach "2": This approach is introduced as Probability-based solution the paper.
///                 - approach "3": This approach is introduced as Micromixing at level \f$ i \f$ with \f$ \tau_s^* \f$ in the paper.
///                 - The last approach: This approach is considered as Dynamic adjustment of \f$ A \f$ value in the paper.
/// \param domainLength_ Length scale of the domain.
/// \param tau0_ Time scale of the domain.
/// \param ScHips_ Vector of Schmidt numbers for HiPS simulation.
/// \note This function sets up the tree based on the Reynolds number. Users pass the method (approach) and the number of levels.


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hips::set_tree(double Re_, double domainLength_, double tau0_, std::vector<double> &ScHips_, std::string approach = "1") {
    Re = Re_;
    domainLength = domainLength_;
    tau0 = tau0_;
    ScHips = ScHips_;

    double originalLevel = (3.0 / 4) * log(1 / Re) / log(Afac);  // Calculate the original level

    if (approach == "1") {
        int lowerLevel = round(originalLevel);                   // Round the original level to the nearest integer
        nL = lowerLevel + 3;                                     // Set the number of levels for the binary tree structure
    } 
    else if (approach == "2") {
        int lowerLevel = ceil(originalLevel);                    // Ceil the original level to the nearest integer
        int upperLevel = lowerLevel - 1;
        Prob = abs((log(originalLevel) - log(lowerLevel)) / (log(upperLevel) - log(lowerLevel)));  // Calculate the probability 
        nL = lowerLevel + 3;                                     // Set the number of levels for the binary tree structure
    }  
    else if (approach == "3") {
        int lowerLevel = ceil(originalLevel);                    // Ceil the original level to the nearest integer
        lStar = std::pow(Re, -3.0 / 4);                          // Calculate lStar based on Re
        nL = lowerLevel + 3;                                     // Set the number of levels for the binary tree structure
    } 
    else if (approach == "4") {
        int closestLevel = round(originalLevel);                 // Round the original level to the nearest integer
        Anew = exp(-log(Re) / ((4.0 / 3.0) * closestLevel));     // Calculate the new value of parameter A
        nL = closestLevel + 3;                                   // Set the number of levels for the binary tree structure
    }
    else {
        // Handle invalid approach case if needed
        throw std::invalid_argument("Invalid approach specified");
    }

    nLevels = nL;
    iEta = nLevels - 3;                                          // Kolmogorov level 

    int maxSc = 1.0;
    for (const auto &sc : ScHips) {
        maxSc = std::max(maxSc, static_cast<int>(sc));
    }
    if (maxSc > 1.0) {
        nLevels += ceil(log(maxSc) / log(4));                    // Changing number of levels!
    }

    Nm1 = nLevels - 1;
    Nm2 = nLevels - 2;
    Nm3 = nLevels - 3;

    nparcels = static_cast<int>(pow(2, Nm1));
    parcelTimes.resize(nparcels, 0);
    i_batchelor.resize(nVar, 0);

    std::vector<double> levelLengths(nLevels);                   // Including all levels, but last 2 don't count
    std::vector<double> levelTaus(nLevels);                      // Smallest scale is 2 levels up from bottom
    levelRates.resize(nLevels);

    for (int i = 0; i < nLevels; ++i) {
        levelLengths[i] = domainLength * pow(Afac, i);
        // levelLengths[i] = domainLength * pow(Anew, i);

        levelTaus[i] = tau0 * pow(levelLengths[i] / domainLength, 2.0 / 3.0) / C_param;
        levelRates[i] = 1.0 / levelTaus[i] * pow(2.0, i);
    }

    LScHips = !ScHips.empty();
    if (LScHips) {                                              // Correct levels for high Sc (levels > Kolmogorov)
        for (int i = iEta + 1; i < nLevels; ++i) {
            levelTaus[i] = tau0 * pow(levelLengths[iEta] / domainLength, 2.0 / 3.0) / C_param;
            levelRates[i] = 1.0 / levelTaus[i] * pow(2.0, i);
        }
    }
    
    eddyRate_total = 0.0;
    for (int i = 0; i <= Nm3; ++i) {
        eddyRate_total += levelRates[i];
    }

    eddyRate_inertial = 0.0;
    for (int i = 0; i <= iEta; ++i) {
        eddyRate_inertial += levelRates[i];
    }

    i_plus.resize(nVar);
    for (int k = 0; k < nVar; ++k) {
        if (ScHips[k] < 1.0) {
            i_batchelor[k] = iEta + 1.5 * log(ScHips[k]) / log(4);
        } else if (ScHips[k] > 1.0) {
            i_batchelor[k] = iEta + log(ScHips[k]) / log(4);
        } else {
            i_batchelor[k] = iEta;
        }
        i_plus[k] = ceil(i_batchelor[k]);
    }

    varRho.resize(nparcels);
    Temp.resize(nparcels);
    pLoc.resize(nparcels);
    for (int i = 0; i < nparcels; ++i) {
        pLoc[i] = i;
    }
}

///////////////////////////////////////////////////////////////////////////////

/// @brief Function to pass the variables, their weights, and their names to the parcels of the tree.  
/// \param v         Vector consisting of variables passed to the HiPS tree.
/// \param w         Vector containing weights for each flow particle.
/// \param varN      Vector containing names corresponding to each variable.

////////////////////////////////////////////////////////////////////////////////////

void hips::set_varData(std::vector<double> &v, std::vector<double> &w, const std::string &varN) {  
    
    varData[currentIndex] = new vector<double>(projection(v, w));
    varName[currentIndex] = varN;

    currentIndex++; 
}

///////////////////////////////////////i///////////////////////////////////////

/// @brief Function to pass the variable, their weights, names, and densities to the parcels of the tree. 
/// \param v         Vector of variables that are passed to the HiPS tree.
/// \param w         Vector of weights; each flow particle has a weight.
/// \param varN      Vector of names of the variable.
/// \param rho       Vector of density; each flow particle has a specific density.
/// \note This function is overloaded. This version considers particle density.

//////////////////////////////////////////////////////////////////////////////////// 
void hips::set_varData(std::vector<double> &v, std::vector<double> &w, const std::string &varN, const std::vector<double> &rho) {
    
    std::pair<std::vector<double>, std::vector<double>> results = projection(v, w, rho);
    std::vector<double> vh = results.first;
    std::vector<double> rho_h = results.second;

    varData[currentIndex] = new std::vector<double>(vh);            
    varRho =  std::vector<double>(rho_h);
    varName[currentIndex] = varN;
 
    currentIndex++; 

}

/////////////////////////////////////////i/////////////////////////////////////////////

/// @brief This function projects the value in flow particles onto HiPS parcels. If works in cases in which density is constant.  
/// It follows \f[
/// \f[
/// \sum_{i=0}^{\text{Number of FP}} (\phi_{\text{FP}} \, \mathrm{d}x_{\text{FP}})_{i} = \sum_{j=0}^{\text{Number of HP}} (\phi_{\text{HP}} \, \mathrm{d}x_{\text{HP}})_{j}
/// \f]
/// \param vcfd      Vector of variables passed to the HiPS tree.
/// \param weight    Weight vector; each flow particle has a weight.

//////////////////////////////////////////////////////////////////////////////////////////////

std::vector<double> hips::projection(std::vector<double> &vcfd, std::vector<double> &weight) {
    
    xc = setGridCfd(weight);                              // Populate the physical domain for flow particles
    xh = setGridHips(nparcels);                           // Populate the physical domain for hips parcels

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

/// @brief This function projects the value in flow particles onto HiPS parcels.
///
/// \param vcfd      Vector of variables passed to the HiPS tree.
/// \param weight    Weight vector; each flow particle has a weight.
/// \param density   Vector of density.
/// \return A pair of vectors representing the projected vector and the density on the grid.
///
/// \note This function is overloaded. This version considers particle density. 

///////////////////////////////////////////////////////////////////////////////////////////////////

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

/// @brief This function generates a physical domain for flow particles based on the given weights.
/// /// It assumes that the length of the domain is 1, meaning each particle occupies a portion of the domain proportional to its weight.
/// The sum of the portions occupied by all particles should equal 1.
/// \param w                     Weight vector defining the spacing between grid points.
/// \return Vector representing the grid positions for flow particles.
////////////////////////////////////////////////////////////////////////////////////

std::vector<double> hips::setGridCfd(std::vector<double> &w) {
    
    std::vector<double> pos;                               // Initializing a vector to hold the grid positions
    double posL = 0.0;                                     // Initializing the starting position

    int i = 0;

    while (i <= w.size()) {                                // Generate the grid positions based on the weights
        pos.push_back(posL);                               // Add the current position to the grid
        posL += w[i];                                      // Move to the next position by adding the corresponding weight
        i++;                                              
                                                   
    }

    return pos;                                           // Return the generated grid positions
}

///////////////////////////////////////////////////////////////////////////////

/// @brief This function generates a physical domain for HiPS parcels.
/// The size of the HiPS domain matches the size of the physical domain specified in the ``setGridCfd()`` function.
/// All parcels occupy an equal portion of the physical domain, meaning the domain is divided evenly among the parcels.
/// \param N                       Number of grid points.
/// \return A vector representing the HiPS grids.

////////////////////////////////////////////////////////////////////////////////////////////

std::vector<double> hips::setGridHips(int N){

    std::vector<double> xh(N + 1);                               // Initialize a vector to hold the grid points
    double step = 1.0 / N;                                       // Calculate the step size

    for(int i = 0; i <= N; i++){                                // Populate the grid with evenly spaced points
        xh[i] = i * step;
    }
        
    return xh;                                                  // Return the generated grid
}

//////////////////////////////////////////////////////////////////////////////////////////////////

/// The HiPS solver
/// \param tRun                               Input simulation run time
/// \param shouldWriteData                    Set to false by default. If true, data will be written. 
/// Sample tee (time of next eddy event)
/// Select and swap subtrees (at current time, not at tee, so that we know who's involved)
/// If the eddy event is at the parcel/micromixing level:
///     React involved parcels from their current time to tee
///     Mix the involved parcels (micromixing)
/// \note Data is written after a specified number of eddy events. By default, the data is written after "10000" eddy events. Users have the flexibility to adjust this number in the code.

//////////////////////////////////////////////////////////////////////////////////////////////////

void hips::calculateSolution(const double tRun, bool shouldWriteData) {
    
    unsigned long long nEddies = 0;                   // Number of eddy events
    int    fileCounter = 0;                           // Number of data files written
    int    iLevel;                                    // Tree level of EE with top at iLevel=0
    int    iTree;                                     // One of two subtrees involved in swap at iLevel
    dtEE;                                                                 
    time = 0.0;                                       // Initialize simulation time

    sample_hips_eddy(dtEE, iLevel);                   // Get first EE at time 0+dtEE
    nEddies++;

    while (time+dtEE<=tRun) {
        time += dtEE;
        selectAndSwapTwoSubtrees(iLevel, iTree);
        advanceHips(iLevel, iTree);    // reaction and micromixing (if needed) to t=time

        sample_hips_eddy(dtEE, iLevel);

        nEddies++;
        //cout<<"-------------------------------------"<<nEddies<<endl;
        if(shouldWriteData && nEddies %10000 == 0) writeData(++fileCounter, time);
    }
    time = tRun;
    iLevel = 0; iTree  = 0;
    if(performReaction)
    reactParcels_LevelTree(iLevel, iTree);      // react all parcels up to end time

    if(shouldWriteData)
        writeInputParameters();
}

///////////////////////////////////////////////////////////////////////////////

/// Samples stochastic eddy events on the hips tree: time and level.
/// \param dtEE                            Time increment to next eddy event (EE)
/// \param iLevel                          Tree level of EE

////////////////////////////////////////////////////////////////////////////////

void hips::sample_hips_eddy(double &dtEE, int &iLevel) {

    static double c1 = 1.0 - pow(2.0, 5.0/3.0*(iEta+1));
    static double c2 = pow(2.0, Nm2) - pow(2.0, iEta+1);
    static double c3 = pow(2.0, iEta+1);

    //--------------- time to next eddy

    double r = rand.getRand();
    dtEE = -log(r)/eddyRate_total;

    //----------------- get eddy level

    r = rand.getRand();

    if ( r <= eddyRate_inertial/eddyRate_total) {     // Inertial region
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

/// Function performs eddy events: parcel swaps.
/// \param iLevel                          Input  level of the tree for the base of the swap.
/// \param iTree                           Output which subtree on the level is selected
///
///  Randomly select a node on iLevel.
///  Go down two levels and select nodes 0q and 1r, where q, r are randomly 0 or 1
///  Find the starting index of the Q-tree and R-tree to swap and the number of parcels.
///  Then swap the cells.
///
///  For a 6 level tree: 0, 1, 2, 3, 4, 5:
///  If iLevel = 1, then suppose i=1, 0q = 00 and 1r = 11:
///  Then we are swaping 0100** with 0111** or (01|00|**) with (01|11|**)
///     or i0qs with i1rs, where i = 01; 0q = 00; 1r = 11; and s = **
///
///   We use bitwise shifts for easy powers of 2.
///   The swap is done by adding or subtracting a value (shift),
///      which should be equivalent to flipping the swapping the two 0q bits and 1r bits.
///                                                                                                              Level
///                                                                                                            ---------
///                                                    *                                                           0
///                                                 /     \
///                                              /           \
///                                           /                 \
///                                        /                       \
///                                     /                             \
///                                  /                                   \
///                               /                                         \
///                            /                                               \
///                           *                                                (*)  01|0000                        1
///                          / \                                               / \
///                        /     \                                           /     \
///                      /         \                                       /         \
///                    /             \                                   /             \
///                  /                 \                               /                 \
///                /                     \                           /                     \
///               *                       *                         *                       *                      2
///              / \                     / \                       / \                     / \
///            /     \                 /     \                   /     \                 /     \
///          /         \             /         \               /         \             /         \
///         *           *           *           *            [*] 00|**    *           *          [*] 11|**         3
///        / \         / \         / \         / \           / \         / \         / \         / \
///       /   \       /   \       /   \       /   \         /   \       /   \       /   \       /   \
///      *     *     *     *     *     *     *     *       *     *     *     *     *     *     *     *             4
///     / \   / \   / \   / \   / \   / \   / \   / \     / \   / \   / \   / \   / \   / \   / \   / \
///    00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15   16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31           5
///                                                      ^^^^^^^^^^^                         ^^^^^^^^^^^
///

/////////////////////////////////////////////////////////////////////////////////////////////////////

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

/// React and mix parcels that are involved in a micromixing process.
/// This is determined by the level and the tree within that level.
/// Might not do anything if the eddy doesn't cause micromixing.
/// \param iLevel                    Input level that the eddy event occurred.
/// \param iTree                     Input root note of the eddy event at iLevel.

/////////////////////////////////////////////////////////////////////////////////

void hips::advanceHips(const int iLevel, const int iTree) {
     
    if (forceTurb==2 && iLevel==0)              // Forcing for statistically stationary
        forceProfile();

    bool rxnDone = false;                       // React all variables once
    for (int k=0; k<nVar; k++) {                // Upon finding first variable needing micromixing
        if ( (iLevel >= i_plus[k]) || 
              //  (iLevel==i_plus[k]-1 && rand.getRand() <=Prob) ) {
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

/// React parcels that are involved in a micromixing process.
/// Parcels might react for different amounts of time depending on when they last
/// reacted, which is stored in the parcelTimes array.
/// \param iLevel                    Input level that the eddy event occurred.
/// \param iTree                     Input root note of the eddy event at iLevel.

////////////////////////////////////////////////////////////////////////////////

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
        #ifdef REACTIONS_ENABLED
             bRxr->react(h, y, dtEE);
            varRho[ime] =  bRxr->getDensity();             //Mb
            Temp[ime] = bRxr->temperature;
        #endif
        varData[0][0][ime] = h;
        for (int k=0; k<nsp; k++)
            varData[k+1][0][ime] = y[k];
    }

}

///////////////////////////////////////////////////////////////////////////////

/// Mix parcels uniformly (using average) at given level and tree
/// Mixing at levels above the lowest enables low Sc variables.
/// \param kVar     Variable index to mix (normally a transported var as determined by caller)
/// \param iLevel   Grandchildren of this iLevel will be mixed
/// \param iTree    At the given iLevel, mix only this subtree
///
/// Example: For a 5 level tree, we have levels 0, 1, 2, 3, 4 (top to bottom).
///   
///          Let iLevel = 2, then nPmix = 2 and we are mixing pairs of parcels at the base of the tree.
///          iTree can be 0, 1, 2, or 3.
///          if iTree = 0 we will mix parcels (0,1) and (2,3) as the left subtree and right subtree:
///                (istart=1, iend=2) and (istart=2, iend=4)
///          if iTree = 1 we will mix parcels (4,5) and (6,7) with (istart=4,iend=6), (istart=6, iend=8)
///          if iTree = 2 we will mix parcels (8,9) and (10,11) with (istart=8,iend=10), (istart=10, iend=12)
///          if iTree = 3 we will mix parcels (12,13) and (14,15) with (istart=12,iend=14), (istart=14, iend=16)
///         
///          Let iLevel=1 then we will be mixing groups of four parcels.
///          iTree can be 0 or 1
///          if iTree = 0 we will mix parcels (0,1,2,3) and (4,5,6,7) with (istart=0,iend=4), (istart=4, iend=8)
///          if iTree = 1 we will mix parcels (8,9,10,11) and (12,13,14,15) with (istart=8,iend=12), (istart=12, iend=16)
///
/// recall: 3 << 4 means 3*2^4 (or 3 = 000011 and 3<<4 = 110000 = 48), that is, we shift the bits left 4 places.
///          
/// \note BE CAREFUL WITH MIXING SOME SCALARS, LIKE MASS FRACTIONS; CURRENT CODE ASSUMES ALL PARCELS HAVE SAME DENSITY (mixing Yi directly)

////////////////////////////////////////////////////////////////////////////////////////////

void hips::mixAcrossLevelTree(int kVar, const int iLevel, const int iTree) {
    
    int istart;
    int iend;

    int nPmix = 1 << (nLevels - iLevel - 2);                   // Number of parcels mixed together

    int ime;

    //---------- Mix left branch of iTree

    istart = iTree << (Nm1-iLevel);  
    iend = istart + nPmix;

    double s = 0;                                               // Initialize sum to 0
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
        
    }   
}

///////////////////////////////////////////////////////////////////////////////

/// @brief Force the HiPS profile to achieve statistical stationarity.
///
/// This function is meant to demonstrate forcing for simple scalars, such as a mixture fraction variable that varies between 0 and 1.
/// The code within this function is configured to force the left half of parcels to average 0 and the right half of parcels to average 1.
/// 
/// \note This function modifies the data stored in the HiPS profile.

///////////////////////////////////////////////////////////////////////////
 
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

///////////////////////////////////////////////////////////////////////////////////

/// @brief Write data to a file with a specific index and output time.
///
/// This function writes data to a file with a filename following a sequential incrementing index pattern 
/// (e.g., Data_00001.dat, Data_00002.dat, etc.). The time of the data file is also written into the file.
///
/// \param ifile The sequential index used in the file name.
/// \param outputTime The time associated with the data, written into the file.
///////////////////////////////////////////////////////////////////////////////////

void hips::writeData(const int ifile, const double outputTime) {
    // Create the "data" directory if it doesn't exist
    if (system("mkdir -p ../data") != 0) {
        cerr << "Error: Unable to create directory ../data" << endl;
        return;
    }

    stringstream ss;
    ss << setw(5) << setfill('0') << ifile;
    string fileName = "../data/Data_" + ss.str() + ".dat";

    ofstream outputFile(fileName);
    if (!outputFile) {
        cerr << "Error: Unable to open file " << fileName << " for writing" << endl;
        return;
    }

    outputFile << "# time = " << outputTime << "\n";
    outputFile << "# Grid Points = " << nparcels << "\n";
    outputFile << "#temperature" << setw(14);
    for (const auto& name : varName)
        outputFile << "#" << setw(14) << name;
    outputFile << scientific << setprecision(10);

    for (int i = 0; i < nparcels; ++i) {
        outputFile << "\n" << setw(19) << Temp[pLoc[i]];            // Write temp data for each parcel
        for (int k = 0; k < nVar; ++k)
            outputFile << setw(19) << varData[k][0][pLoc[i]];       // Write data for each variable
    }

    outputFile.close();
}

///////////////////////////////////////////////////////////////////////////////

///  @brief Function for projecting vectors onto a grid.
///  \param vb Vector to be projected back.

//////////////////////////////////////////////////////////////////////////////

std::vector<double> hips::projection_back(std::vector<double> &vh) {
    int nh = xh.size() - 1;
    int nc = xc.size() - 1;

    std::vector<double> vc(nc, 0.0);
    int jprev = 0;

    for (double val : vc)
        std::cout << val << " ";
    std::cout << std::endl;

    for (int i = 0; i < nc; ++i) {
        for (int j = jprev + 1; j <= nh; ++j) {
            //std::cout << "first.    i. " << i << "  j. " << j << std::endl;
            if (xh[j] <= xc[i + 1]) {
                double Δ1 = xh[j] - xh[j - 1];
                double Δ2 = xh[j] - xc[i];
                double Δ = std::min(Δ1, Δ2);

                vc[i] += vh[j - 1] * Δ;
            } else {
                double Δ1 = xc[i + 1] - xh[j - 1];
                double Δ2 = xc[i + 1] - xc[i];
                double Δ = std::min(Δ1, Δ2);

                vc[i] += vh[j - 1] * Δ;

                jprev = j - 1;
                break;
            }
        }
        vc[i] /= (xc[i + 1] - xc[i]);
    }
    return vc;
}

/////////////////////////////////////////////////////////////////////////////////////////

/// \brief Retrieve modified data from the Hierarchical Progressive Survey (HiPS) library.
///
/// This function retrieves modified data from the HiPS library and stores it in the provided vector.
/// Each element of the returned vector contains a vector representing a single data projection.
///
/// \return A vector of vectors containing modified data retrieved from the HiPS library.

//////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> hips::get_varData(){

    std::vector<std::vector<double>> varDataProjections;                             // Vector to store modified data projections

    for (int i = 0; i < varData.size(); i++)                                        // Loop through each element of varData and project the data back
        varDataProjections.push_back(projection_back(varData[i][0]));               // Project the data back and store it in vh

    return varDataProjections;                                                      // Return the vector containing modified data projections
}

/////////////////////////////////////////////////////////////////////////////////////

/// \brief Writes the input parameters to a YAML file for post-process simulations.
/// 
/// This function creates a YAML file named "InputParameters.yaml" in the "../data/" directory.
/// The file includes the following simulation parameters:
/// - nLevels: The number of levels in the simulation.
/// - nparcels: The number of parcels used in the simulation.
/// - nVar: The number of variables in the simulation.
/// - varName: A list of variable names.

/////////////////////////////////////////////////////////////////////////////////////

void hips::writeInputParameters() {
    // Construct the YAML file name
    string yamlFileName = "../data/InputParameters.yaml";

    YAML::Node yamlData;                                // Create a YAML node and write the input parameters to it

    yamlData["params"]["nLevels"] = nLevels;            // Adding simulation parameters with comments
    yamlData["params"]["nparcels"] = nparcels;
    yamlData["params"]["nVar"] = nVar;

    YAML::Node varNamesNode;                             // Adding variable names with comments

    for (const auto& name : varName) {
        varNamesNode.push_back(name);
    }
    yamlData["variable_names"] = varNamesNode;


    std::ofstream yamlFile(yamlFileName);                  // Save the YAML node to a file
    if (!yamlFile) {
        cerr << "Error: Unable to open file " << yamlFileName << " for writing" << endl;
        return;
    }

    yamlFile << yamlData;
    yamlFile.close();
}
 


 


