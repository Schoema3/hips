#include "hips.h"

#ifdef REACTIONS_ENABLED
#include "batchReactor_cvode.h"
#include "batchReactor_cantera.h"
#endif

#include "randomGenerator.h"
#include <yaml-cpp/yaml.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Constructor for initializing parameters to create the HiPS tree.
///
/// This constructor sets up the necessary parameters for creating the HiPS tree,
/// including the eddy rate, turbulence forcing, number of variables, and the 
/// Cantera solution object (if reactions are enabled). It is useful for initializing 
/// the tree multiple times, as the parameters are set once and reused. Typically, 
/// this constructor is followed by a call to `set_tree(...)`.
///
/// \param C_param_         Controls the eddy rate.
/// \param forceTurb_       Flag indicating whether to force turbulence (non-zero for true).
/// \param nVar_            Number of variables in the simulation.
/// \param cantSol          Shared pointer to the Cantera solution object (used only if reactions are enabled).
/// \param performReaction_ Flag indicating whether chemical reactions should be performed.
/// \param seed             Seed for the random number generator (if negative, a random seed is generated).
///
/// \note By default, the integrator object `bRxr` is initialized as `batchReactor_cvode`.
///       To switch to `batchReactor_cantera`, uncomment the corresponding line in the code.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

    bRxr = make_unique<batchReactor_cvode>(cantSol);                             // By default, use batchReactor_cvode

    // Uncomment the following line to switch to batchReactor_cantera
    // bRxr = make_unique<batchReactor_cantera>(cantSol);
#endif

    // Resize vectors to the number of variables
    varData.resize(nVar);
    varName.resize(nVar);        
} 

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Constructor for initializing parameters to create the HiPS tree.
///
/// This constructor sets up all necessary parameters in a single step, including the length scale, 
/// time scale, eddy rate, and other key variables. It is particularly useful for scenarios where the 
/// HiPS tree needs to be initialized only once, such as in simple mixing simulations.
///
/// \param nLevels_          Number of levels in the HiPS tree.
/// \param domainLength_     Length scale of the domain.
/// \param tau0_             Time scale of the domain.
/// \param C_param_          Parameter controlling the eddy rate.
/// \param forceTurb_        Flag indicating whether to force turbulence (non-zero for true).
/// \param nVar_             Number of variables in the simulation.
/// \param ScHips_           Vector containing Schmidt numbers for the HiPS simulation.
/// \param cantSol           Shared pointer to the Cantera solution object (used only if reactions are enabled).
/// \param performReaction_  Flag indicating whether chemical reactions should be performed.
/// \param seed              Seed for the random number generator (if negative, a random seed is generated).
///
/// \note By default, the integrator object `bRxr` is initialized as `batchReactor_cvode`.
///       To switch to `batchReactor_cantera`, uncomment the corresponding line in the code.
/////////////////////////////////////////////////////////////////////////////////////////////////////////

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

    bRxr = make_unique<batchReactor_cvode>(cantSol);                                // By default, use batchReactor_cvode

    // Uncomment the following line to switch to batchReactor_cantera
    // bRxr = make_unique<batchReactor_cantera>(cantSol);
#endif

    // Resize vectors to the number of variables
    varData.resize(nVar);
    varName.resize(nVar); 

    set_tree(nLevels, domainLength, tau0, ScHips);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Creates a tree structure based on the specified parameters.
///
/// This function sets up the tree structure according to the provided number of levels, 
/// domain length, time scale, and Schmidt numbers. It adjusts the number of levels if 
/// necessary, based on the maximum Schmidt number, and initializes various internal 
/// parameters used in the HiPS simulation.
///
/// \param nLevels_         Number of levels in the tree. If set to -1, the function uses the value of `nL`.
/// \param domainLength_    Length scale of the domain.
/// \param tau0_            Time scale of the domain.
/// \param ScHips_          Vector of Schmidt numbers for the HiPS simulation.
///
/// \note The function computes and adjusts the number of levels based on the maximum 
///       Schmidt number if it exceeds 1.0. It also calculates various level-specific 
///       properties like length scales, time scales, and rates used in the HiPS model.
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Creates a tree structure based on the specified Reynolds number and approach.
///
/// This function sets up the HiPS tree structure according to the provided Reynolds number (Re), domain length, 
/// time scale, and Schmidt numbers, using a specified approach to determine the number of levels. 
/// The function offers multiple methods to set the number of levels, based on different strategies outlined in the related paper.
///
/// \param Re_              The Reynolds number, used to determine the original level of the tree.
/// \param approach_        The method for setting the number of levels based on the Reynolds number:
///                         - "1": Rounding to the closest level to \f$i_s^*\f$ for micromixing.
///                         - "2": Probability-based solution.
///                         - "3": Micromixing at level \f$i\f$ with \f$\tau_s^*\f$.
///                         - "4": Dynamic adjustment of \f$A\f$ value.
/// \param domainLength_    Length scale of the domain.
/// \param tau0_            Time scale of the domain.
/// \param ScHips_          Vector of Schmidt numbers for the HiPS simulation.
///
/// \note This function adjusts the number of levels based on the Reynolds number and the selected approach, 
///       and initializes various parameters used for the HiPS simulation.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hips::set_tree(double Re_, double domainLength_, double tau0_, std::vector<double> &ScHips_, std::string approach_) {
    Re = Re_;
    domainLength = domainLength_;
    tau0 = tau0_;
    ScHips = ScHips_;
    approach = approach_;

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

    //----------------------------------------------------------------------

    nLevels = nL;
    iEta = nLevels - 3;                                          // Kolmogorov level 

    int maxSc = 1;
    for (const auto &sc : ScHips)
        maxSc = std::max(maxSc, static_cast<int>(sc));
    if (maxSc > 1.0)
        nLevels += ceil(log(maxSc) / log(4));                    // Changing number of levels!

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

        if (approach == "4")                            //Mb---------------
            levelLengths[i] = domainLength * pow(Anew, i);

        levelTaus[i] = tau0 * pow(levelLengths[i] / domainLength, 2.0 / 3.0) / C_param;
        levelRates[i] = 1.0 / levelTaus[i] * pow(2.0, i);

        if (approach == "3") {                         //Mb--------------
            levelTaus[Nm3] = tau0 * pow(lStar / domainLength, 2.0 / 3.0) / C_param;
            levelRates[Nm3] = 1.0 / levelTaus[Nm3] * pow(2.0, Nm3);
        }
    }

    LScHips = !ScHips.empty();
    if (LScHips)                                                // Correct levels for high Sc (levels > Kolmogorov)
        for (int i = iEta + 1; i < nLevels; ++i) {
            levelTaus[i] = tau0 * pow(levelLengths[iEta] / domainLength, 2.0 / 3.0) / C_param;
            levelRates[i] = 1.0 / levelTaus[i] * pow(2.0, i);
        }
   
    //-----------------------------------------------

    eddyRate_total = 0.0;
    for (int i = 0; i <= Nm3; ++i)
        eddyRate_total += levelRates[i];

    eddyRate_inertial = 0.0;
    for (int i = 0; i <= iEta; ++i)
        eddyRate_inertial += levelRates[i];

    //-------------------------------------------------

    i_plus.resize(nVar);
    for (int k = 0; k < nVar; ++k) {
        if (ScHips[k] < 1.0)
            i_batchelor[k] = iEta + 1.5 * log(ScHips[k]) / log(4);
        else if (ScHips[k] > 1.0)
            i_batchelor[k] = iEta + log(ScHips[k]) / log(4);
        else
            i_batchelor[k] = iEta;
        i_plus[k] = ceil(i_batchelor[k]);
    }
 
    //-----------------------------------------------------

    varRho.resize(nparcels);
    Temp.resize(nparcels);
    pLoc.resize(nparcels);
    for (int i = 0; i < nparcels; ++i)
        pLoc[i] = i;
}

////////////////////////////////////////////////////////////////////////////////////
/// \brief Assigns variables, their corresponding weights, and names to the parcels in the HiPS tree.
///
/// This function assigns the given variables, their associated weights, and names to the specific parcels 
/// within the HiPS tree structure. The function performs a weighted projection of the variables onto the parcels.
///
/// \param v         Vector of variables to be assigned to the HiPS tree.
/// \param w         Vector of weights corresponding to each flow particle.
/// \param varN      Name of the variable being assigned, stored as a string.
////////////////////////////////////////////////////////////////////////////////////

void hips::set_varData(std::vector<double> &v, std::vector<double> &w, const std::string &varN) {  
    
    varData[currentIndex] = new vector<double>(projection(v, w));
    varName[currentIndex] = varN;

    currentIndex++; 
}

//////////////////////////////////////////////////////////////////////////////////// 
/// \brief Assigns variables, weights, names, and densities to the parcels in the HiPS tree.
///
/// This overloaded function assigns the specified variables, along with their associated weights,
/// names, and densities, to the parcels within the HiPS tree structure. It considers particle density
/// during the projection process.
///
/// \param v         Vector of variables to be assigned to the HiPS tree.
/// \param w         Vector of weights for each flow particle.
/// \param varN      Name of the variable being assigned, stored as a string.
/// \param rho       Vector of densities for each flow particle.
///
/// \note This function is overloaded to account for particle density in the HiPS tree.
////////////////////////////////////////////////////////////////////////////////////

void hips::set_varData(std::vector<double> &v, std::vector<double> &w, const std::string &varN, const std::vector<double> &rho) {
    
    std::pair<std::vector<double>, std::vector<double>> results = projection(v, w, rho);
    std::vector<double> vh = results.first;
    std::vector<double> rho_h = results.second;

    varData[currentIndex] = new std::vector<double>(vh);            
    varRho = std::vector<double>(rho_h);
    varName[currentIndex] = varN;
 
    currentIndex++; 
}
//////////////////////////////////////////////////////////////////////////////////////////////
/// \brief This function projects the values in flow particles onto HiPS parcels. It works in cases where density is constant.
/// 
/// This function follows the equation:
/// \f[
/// \sum_{i=0}^{\text{Number of Flow Particles (FP)}} (\phi_{\text{FP}} \, \mathrm{d}x_{\text{FP}})_{i} = \sum_{j=0}^{\text{Number of HiPS Parcels (HP)}} (\phi_{\text{HP}} \, \mathrm{d}x_{\text{HP}})_{j}
/// \f]
/// 
/// \param vcfd          Vector of variables passed to the HiPS tree.
/// \param weight        Weight vector; each flow particle has a weight.
/// \return              Vector of values projected onto HiPS parcels.
///////////////////////////////////////////////////////////////////////////////

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

///////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Projects the values from flow particles onto HiPS parcels, considering particle density.
///
/// This function calculates the projection of values from flow particles onto HiPS parcels, 
/// taking into account the density of the particles. The function returns a pair of vectors: 
/// one representing the projected values and one for the densities on the HiPS grid.
///
/// \param vcfd         Vector of variables from flow particles to be projected onto HiPS parcels.
/// \param weight       Weight vector associated with each flow particle.
/// \param density      Vector of densities for each flow particle.
///
/// \return A pair of vectors:
///         - First vector: Projected values on the HiPS parcels.
///         - Second vector: Densities on the HiPS parcels.
///
/// \note This function is overloaded. This version includes particle density in the projection process.
///////////////////////////////////////////////////////////////////////////////////////////////////

std::pair<std::vector<double>, std::vector<double>> hips::projection(std::vector<double> &vcfd, 
                                                                     std::vector<double> &weight, 
                                                                     const std::vector<double> &density) {
    
    std::vector<double> xc = setGridCfd(weight);
    std::vector<double> xh = setGridHips(nparcels);

    // Initialize variables
    int nc = xc.size() - 1;
    int nh = xh.size() - 1;

    std::vector<double> vh(nh, 0.0);
    std::vector<double> rho_h(nh, 0.0);
    int jprev = 0;

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

        // ------------------------ Normalize results

        rho_h[i] /= (xh[i + 1] - xh[i]);
        vh[i] /= rho_h[i] * (xh[i + 1] - xh[i]);
    }
    return {vh, rho_h};
}

/////////////////////////////////////////////////////////////////////////////////
/// \brief Generates a physical domain for flow particles based on the given weights.
///
/// This function assumes that the total length of the domain is 1, and each particle occupies a portion 
/// of the domain proportional to its weight. The sum of all particle portions should equal 1.
///
/// \param w         Weight vector defining the spacing between grid points.
/// \return          Vector representing the grid positions for flow particles.
/////////////////////////////////////////////////////////////////////////////////

std::vector<double> hips::setGridCfd(std::vector<double> &w) {
    
    std::vector<double> pos;      // Initializing a vector to hold the grid positions
    double posL = 0.0;            // Initializing the starting position

    int i = 0;

    while (i <= w.size()) {       // Generate the grid positions based on the weights
        pos.push_back(posL);      // Add the current position to the grid
        posL += w[i];             // Move to the next position by adding the corresponding weight
        i++;                                              
    }
    return pos;                   // Return the generated grid positions
}

///////////////////////////////////////////////////////////////////////////////
/// \brief Generates a physical domain for HiPS parcels.
///
/// This function creates a grid for HiPS parcels, with each parcel occupying an equal portion 
/// of the physical domain. The size of the HiPS domain matches the size of the domain specified 
/// in the `setGridCfd()` function.
///
/// \param N         Number of grid points.
/// \return          Vector representing the HiPS grid positions.
///////////////////////////////////////////////////////////////////////////////

std::vector<double> hips::setGridHips(int N){

    std::vector<double> xh(N + 1);                               // Initialize a vector to hold the grid points
    double step = 1.0 / N;                                       // Calculate the step size

    for(int i = 0; i <= N; i++)                                 // Populate the grid with evenly spaced points
        xh[i] = i * step;
        
    return xh;                                                  // Return the generated grid
}

/////////////////////////////////////////////////////////////////////////////////
/// \brief The HiPS solver for simulating the interaction of parcels within a turbulent flow.
///
/// \param tRun                Simulation run time.
/// \param shouldWriteData     Flag indicating whether to write data during the simulation. Default is false.
///
/// This function simulates the HiPS process by iterating over eddy events until the specified simulation time is reached.
/// It samples the time of the next eddy event (tee), selects and swaps subtrees, and handles reactions and micromixing 
/// at the parcel level if necessary. Data is written after a specified number of eddy events (default is 10,000), but this 
/// number can be adjusted by the user.
///
/// \note Data is written after a specified number of eddy events. By default, data is written after 10,000 eddy events. 
/// Users can modify this threshold within the code.
///////////////////////////////////////////////////////////////////////////////////

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
        advanceHips(iLevel, iTree);                   // reaction and micromixing (if needed) to t=time

        sample_hips_eddy(dtEE, iLevel);

        nEddies++;
        if (shouldWriteData && nEddies %10000 == 0) 
            writeData(++fileCounter, time);
    }
    time = tRun;
    iLevel = 0; iTree  = 0;
    if (performReaction)
    reactParcels_LevelTree(iLevel, iTree);            // react all parcels up to end time

    if (shouldWriteData)
        writeInputParameters();
}

///////////////////////////////////////////////////////////////////////////////
/// \brief Samples stochastic eddy events on the HiPS tree, determining both the time increment to the next eddy event (dtEE) and the tree level of the eddy event (iLevel).
///
/// \param dtEE         Time increment to the next eddy event (EE).
/// \param iLevel       Tree level at which the eddy event occurs.
///
/// This function uses stochastic sampling to determine when and at what level in the HiPS tree an eddy event will occur.
/// The time to the next eddy event is sampled based on the total eddy rate, and the level is determined by whether the event
/// occurs in the inertial or Batchelor region of the turbulence.
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
        if (iLevel < iEta+1) iLevel = iEta+1;
        if (iLevel > Nm3) iLevel = Nm3;
    }
    return;
}

///////////////////////////////////////////////////////////////////////////////
/// Function performs eddy events: parcel swaps.
///
/// \param iLevel         Input  level of the tree for the base of the swap.
/// \param iTree          Output which subtree on the level is selected
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
/// This is determined by the level nd the tree within that level.
/// Might not do anything if the eddy doesn't cause micromixing.
///
/// \param iLevel         Input level that the eddy event occurred.
/// \param iTree          Input root note of the eddy event at iLevel.
/////////////////////////////////////////////////////////////////////////////////

void hips::advanceHips(const int iLevel, const int iTree) {
    if (forceTurb == 2 && iLevel == 0) {
        forceProfile();                                                  // Forcing for statistically stationary
    }

    bool rxnDone = false;                                               // React all variables once
    for (int k = 0; k < nVar; k++) {                                    // Upon finding first variable needing micromixing
        // Combined condition check with approach condition
        if ((iLevel >= i_plus[k]) || 
            (iLevel == i_plus[k] - 1 && rand.getRand() <= i_plus[k] - i_batchelor[k]) || 
            (approach == "2" && iLevel == i_plus[k] - 1 && rand.getRand() <= Prob)) {
                if (!rxnDone && performReaction) {
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
///
/// \param iLevel         Input level that the eddy event occurred.
/// \param iTree          Input root note of the eddy event at iLevel.
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
        if(performReaction) {
                bRxr->react(h, y, dtEE);
                varRho[ime] =  bRxr->getDensity();            
            Temp[ime] = bRxr->temperature;
        }
#endif
        varData[0][0][ime] = h;
        for (int k=0; k<nsp; k++)
            varData[k+1][0][ime] = y[k];
    }
}

///////////////////////////////////////////////////////////////////////////////
/// Mix parcels uniformly (using average) at given level and tree
/// Mixing at levels above the lowest enables low Sc variables.
///
/// \param kVar         Variable index to mix (normally a transported var as determined by caller)
/// \param iLevel       Grandchildren of this iLevel will be mixed
/// \param iTree        At the given iLevel, mix only this subtree
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
/// \brief Force the HiPS profile to achieve statistical stationarity.
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
/// \brief Write data to a file with a specific index and output time.
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

    if (performReaction)
        outputFile << "#temperature" << setw(14);

    for (const auto& name : varName)
        outputFile << "#" << setw(14) << name;
    
    outputFile << scientific << setprecision(10);

    for (int i = 0; i < nparcels; ++i) {
        if (performReaction)
            outputFile << "\n" << setw(19) << Temp[pLoc[i]];  // Write temp data for each parcel if reaction occurs
        else
            outputFile << "\n";
        for (int k = 0; k < nVar; ++k)
            outputFile << setw(19) << varData[k][0][pLoc[i]]; // Write data for each variable
    }

    outputFile.close();
}

///////////////////////////////////////////////////////////////////////////////
/// \brief This function projects the HiPS parcel values back onto the flow particles.
/// 
/// It effectively reverses the projection process to ensure the values in HiPS parcels 
/// are accurately redistributed to the flow particles.
/// 
/// It follows:
/// \f[
/// \sum_{j=0}^{\text{Number of HP}} (\phi_{\text{HP}} \, \mathrm{d}x_{\text{HP}})_{j} = \sum_{i=0}^{\text{Number of FP}} (\phi_{\text{FP}} \, \mathrm{d}x_{\text{FP}})_{i}
/// \f]
/// 
/// \param vh           Vector of values from HiPS parcels to be projected back.
/// \return             Vector of values projected back onto the flow particles.
/// 
/// \note This function is the reverse of the projection function.
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
            if (xh[j] <= xc[i + 1]) {
                double d1 = xh[j] - xh[j - 1];
                double d2 = xh[j] - xc[i];
                double d = std::min(d1, d2);

                vc[i] += vh[j - 1] * d;
            } else {
                double d1 = xc[i + 1] - xh[j - 1];
                double d2 = xc[i + 1] - xc[i];
                double d = std::min(d1, d2);

                vc[i] += vh[j - 1] * d;

                jprev = j - 1;
                break;
            }
        }
        vc[i] /= (xc[i + 1] - xc[i]);
    }
    return vc;
}

/////////////////////////////////////////////////////////////////////////////////////////

/// \brief This function projects the HiPS parcel values and densities back onto the flow particles.
/// 
/// It effectively reverses the projection process to ensure both the values and densities in HiPS parcels 
/// are accurately redistributed to the flow particles.
/// 
/// It follows:
/// \f[
/// \sum_{j=0}^{\text{Number of HP}} (\phi_{\text{HP}} \, \rho_{\text{HP}} \, \mathrm{d}x_{\text{HP}})_{j} = \sum_{i=0}^{\text{Number of FP}} (\phi_{\text{FP}} \, \rho_{\text{FP}} \, \mathrm{d}x_{\text{FP}})_{i}
/// \f]
/// 
/// \param vh           Vector of values from HiPS parcels to be projected back.
/// \param rho_h        Vector of density values from HiPS parcels.
/// \return             A pair of vectors: 
///                     - The first vector contains the values projected back onto the flow particles.
///                     - The second vector contains the corresponding densities for the flow particles.
/// 
/// \note This function is the reverse of the projection function with density.
/////////////////////////////////////////////////////////////////////////////////////////////////////////

std::pair<std::vector<double>, std::vector<double>> hips::projection_back_with_density(std::vector<double> &vh, 
                                                                                       std::vector<double> &rho_h) {
    int nh = xh.size() - 1;
    int nc = xc.size() - 1;

    std::vector<double> vc(nc, 0.0);
    std::vector<double> rho_c(nc, 0.0);
    int jprev = 0;

    for (int i = 0; i < nc; ++i) {
        double total_dx = 0.0;

        for (int j = jprev + 1; j <= nh; ++j) {
            if (xh[j] <= xc[i + 1]) {
                double d = std::min(xh[j] - xh[j - 1], xh[j] - xc[i]);
                total_dx += d;
                rho_c[i] += rho_h[j - 1] * d;
                vc[i] += vh[j - 1] * rho_h[j - 1] * d;
            } else {
                double d = std::min(xc[i + 1] - xh[j - 1], xc[i + 1] - xc[i]);
                total_dx += d;
                rho_c[i] += rho_h[j - 1] * d;
                vc[i] += vh[j - 1] * rho_h[j - 1] * d;

                jprev = j - 1;
                break;
            }
        }

        // Normalize the results
        vc[i] /= rho_c[i];
    }
    return {vc, rho_c};
}

/////////////////////////////////////////////////////////////////////////////////////////
/// \brief This function returns final data from the simulation.
///
/// \note This function is used for integrating HiPS as a subgrid model in CFD simulations.
///
/// \return A vector of vectors containing the final results.
//////////////////////////////////////////////////////////////////////////////////////////

std::vector<std::vector<double>> hips::get_varData(){

    std::vector<std::vector<double>> varDataProjections;               // Vector to store modified data projections

    for (int i = 0; i < varData.size(); i++)                           // Loop through each element of varData and project the data back
        varDataProjections.push_back(projection_back(varData[i][0]));  // Project the data back and store it in vh

    return varDataProjections;                                         // Return the vector containing modified data projections
}
////////////////////////////////////////////////////////////////////////////////////
/// \brief This function returns final data from the simulation, including densities.
///
/// This function processes the HiPS data, including both values and densities, and 
/// projects them back onto the flow particles.
///
/// \note This function is used for integrating HiPS as a subgrid model in CFD simulations, 
///       taking into account both the values and densities of parcels.
///
/// \return A vector of pairs:
///         - Each pair contains two vectors: 
///             - The first vector contains the final results for the values.
///             - The second vector contains the corresponding density results.
///////////////////////////////////////////////////////////////////////////////////

std::vector<std::pair<std::vector<double>, std::vector<double>>> hips::get_varData_with_density() {

    std::vector<std::pair<std::vector<double>, std::vector<double>>> varDataProjections;  // Vector to store data projections with densities

    for (int i = 0; i < varData.size(); i++) {
        // Extract the value and density data
        std::vector<double> vh = varData[i][0];  // Assuming varData[i][0] holds the values
        std::vector<double> rho_h = varRho;      // Assuming varData[i][1] holds the densities

        // Project the data back with density and store the result
        varDataProjections.push_back(projection_back_with_density(vh, rho_h));
    }

    return varDataProjections;  // Return the vector containing modified data projections with densities
}

/////////////////////////////////////////////////////////////////////////////////////
/// \brief This function writes the input parameters to a YAML file for post-processes simulations.
/// 
/// This function creates a YAML file named "InputParameters.yaml" in the "../data/" directory.
///
/// \param nLevels         The number of levels in the simulation.
/// \param nparcels        The number of parcels used in the simulation.
/// \param nVar            The number of variables in the simulation.
/// \param varName         A list of variable names.
///
/// \note Users can add other simulations parameters to this function. 
/////////////////////////////////////////////////////////////////////////////////////

void hips::writeInputParameters() {
    // Construct the YAML file name
    string yamlFileName = "../data/InputParameters.yaml";

    YAML::Node yamlData;               // Create a YAML node and write the input parameters to it

    // Adding simulation parameters with comments
    yamlData["params"]["nLevels"] = nLevels;           
    yamlData["params"]["nparcels"] = nparcels;
    yamlData["params"]["nVar"] = nVar;

    YAML::Node varNamesNode;

    // Adding variable names with comments
    for (const auto& name : varName)
        varNamesNode.push_back(name);
    yamlData["variable_names"] = varNamesNode;


    // Save the YAML node to a file
    std::ofstream yamlFile(yamlFileName);                  
    if (!yamlFile) {
        cerr << "Error: Unable to open file " << yamlFileName << " for writing" << endl;
        return;
    }

    yamlFile << yamlData;
    yamlFile.close();
}
