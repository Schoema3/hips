#include "HiPS.h"

#ifdef REACTIONS_ENABLED
#include "BatchReactor_cvode.h"
#include "BatchReactor_cantera.h"
#endif

#include "RandomGenerator.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

HiPS::HiPS(int nLevels,
           double domainLength,
           double tau0,
           double C_param,
           bool forceTurb,
           int nVar,
           vector<double> &ScHips,
           bool performReaction,
           shared_ptr<void> vcantSol,
           int seed,
           int realization):
    _nLevels(nLevels),
    _domainLength(domainLength),
    _tau0(tau0),
    _C_param(C_param),
    _forceTurb(forceTurb),
    _ScHips(ScHips),
    _nVar(nVar),
    rand(seed),
    _performReaction(performReaction),
    _realization(realization){

    #ifdef REACTIONS_ENABLED
    if(_performReaction) {
        shared_ptr<Cantera::Solution> cantSol = static_pointer_cast<Cantera::Solution>(vcantSol);
        _gas = cantSol->thermo();
        _nsp =    _gas->nSpecies();

        _bRxr = make_shared<BatchReactor_cvode>(cantSol);                                // By default, use BatchReactor_cvode

        // Uncomment the following line to switch to BatchReactor_cantera
        // _bRxr = make_unique<BatchReactor_cantera>(cantSol);

        if(forceTurb)
            throw std::runtime_error("Error: forceTurb should be false if preformReaction is true");
    }
    #endif

    // Resize vectors to the number of variables
    _varData.resize(nVar);
    _varName.resize(nVar);

    set_tree(nLevels, domainLength, tau0);
}

HiPS::HiPS(double C_param,
           bool forceTurb,
           int nVar,
           vector<double> &ScHips,
           bool performReaction,
           shared_ptr<void> vcantSol,
           int seed,
           int realization):
    _C_param(C_param),
    _forceTurb(forceTurb),
    _nVar(nVar),
    _ScHips(ScHips),
    rand(seed),
    _performReaction(performReaction),
    _realization(realization){

    #ifdef REACTIONS_ENABLED
    // Initialize Cantera thermo phase and species count.
    if(_performReaction) {
        shared_ptr<Cantera::Solution> cantSol = static_pointer_cast<Cantera::Solution>(vcantSol);
        _gas = cantSol->thermo();
        _nsp =    _gas->nSpecies();

        // Set up the default batch reactor (cvode).
       // _bRxr = make_shared<BatchReactor_cvode>(cantSol);

        // Uncomment the following line to switch to BatchReactor_cantera.
         _bRxr = make_shared<BatchReactor_cantera>(cantSol);
    }
    #endif

    // Resize vectors to accommodate the number of variables.
    _varData.resize(nVar);
    _varName.resize(nVar);
}

void HiPS::set_tree(int nBaseLevels, double domainLength, double tau0){
 
    _nLevels = nBaseLevels;
    _domainLength = domainLength;
    _tau0 = tau0;

    if (_nLevels == -1)
        _nLevels = _nL;

    _iEta = _nLevels - 3;                                // Kolmogorov level; if _nLevels = 7, then 0, 1, 2, 3, (4), 5, 6; _iEta=4 is the lowest swap level: swap grandchildren of _iEta=4 at level 6.
           
    int maxSc = 1.0;

    for (int i=0; i<_ScHips.size(); i++)
        maxSc = _ScHips[i]>maxSc ? _ScHips[i] : maxSc;
    
    if (maxSc > 1.0)
        _nLevels += ceil(log(maxSc)/log(4));            // Changing number of levels!
    
    _Nm1 = _nLevels - 1;
    _Nm2 = _nLevels - 2;
    _Nm3 = _nLevels - 3;
    
    // -------------------------- 
    
    _nparcels = static_cast<int>(pow(2, _Nm1));
    _parcelTimes.resize(_nparcels,0);
    _i_batchelor.resize(_nVar,0);
    
    vector<double> levelLengths(_nLevels);              // Including all levels, but last 2 don't count:
    vector<double> levelTaus(_nLevels);                 // Smallest scale is 2 levels up from bottom
    _levelRates   = vector<double>(_nLevels);

    for (int i=0; i<_nLevels; i++) {
        levelLengths[i] = domainLength * pow(_Afac,i);
       //levelLengths[i] = domainLength * pow(_Anew,i);

        levelTaus[i] = tau0 * pow(levelLengths[i]/domainLength, 2.0/3.0) / _C_param;
        _levelRates[i] = 1.0/levelTaus[i] * pow(2.0,i);
    }

    _LScHips = _ScHips.size() > 0 ? true : false;
    if (_LScHips) {                                     // Ccorrect levels for high Sc (levels > Kolmogorov)
        for (int i=_iEta+1; i<_nLevels; i++) {
            levelTaus[i] = tau0 *
            pow(levelLengths[_iEta]/domainLength, 2.0/3.0) / _C_param;
            _levelRates[i] = 1.0/levelTaus[i] * pow(2.0,i);
        }
    }

    //-------------------------------------------------

    _eddyRate_total = 0.0;
    for (int i=0; i<=_Nm3; i++)
        _eddyRate_total += _levelRates[i];
    
    _eddyRate_inertial = 0.0;
    for (int i=0; i<=_iEta; i++)
        _eddyRate_inertial += _levelRates[i];
    
    //-------------------
    
    _i_plus.resize(_nVar);
    
    for (int k=0; k<_nVar; k++) {
        if (_ScHips[k] < 1.0)
            _i_batchelor[k] = _iEta + 1.5*log(_ScHips[k])/log(4);
        else if (_ScHips[k] > 1.0)
            _i_batchelor[k] = _iEta + log(_ScHips[k])/log(4);
        else
            _i_batchelor[k] = _iEta;

        _i_plus[k] = ceil(_i_batchelor[k]);
    }
    
    //------------------- Set the parcel addresses (index array)

    _varRho.resize(_nparcels);
    _wPar.assign(_nparcels, 1.0 / _nparcels);
 
    _pLoc.resize(_nparcels);
    for (int i=0; i<_nparcels; i++)
        _pLoc[i] = i;

    _currentIndex = 0;
} 

void HiPS::set_tree(double Re, double domainLength, double tau0, std::string ReApproach){
    _Re = Re;
    _domainLength = domainLength;
    _tau0 = tau0;
    _ReApproach = ReApproach;

    double baseLevelEstimate = (3.0 / 4) * log(1 / _Re) / log(_Afac);                               // Calculate the base tree level estimate (non-integer)
    int baseLevel;

    if (ReApproach == "rounding") {
        baseLevel = round(baseLevelEstimate);                                                     // Round the base level to the nearest integer
    } 
    else if (ReApproach == "probability") {
        baseLevel = ceil(baseLevelEstimate);                                                      // Ceil the base level to the nearest integer
        int previousLevel = baseLevel - 1;
        _probability = baseLevelEstimate - previousLevel;                                                 // Calculate the probability
    } 
    else if (ReApproach == "micromixing") {
        baseLevel = ceil(baseLevelEstimate);                                                      // Ceil the base level to the nearest integer
        _lStar = std::pow(_Re, -3.0 / 4);                                                           // Calculate _lStar based on _Re
    } 
    else if (ReApproach == "dynamic_A") {
        baseLevel = round(baseLevelEstimate);                                                     // Round the base level to the nearest integer
        _Anew = exp(-log(_Re) / ((4.0 / 3.0) * baseLevel));                                         // Calculate the new value of parameter A
    } 
    else {
        throw std::invalid_argument("Invalid ReApproach specified");                                // Handle invalid approach case if needed
    }

    _nL = baseLevel + 3;                                                                           // Set the number of levels for the binary tree structure

    //----------------------------------------------------------------------
    _nLevels = _nL;
    _iEta = _nLevels - 3;  // Kolmogorov level

    int maxSc = 1;
    for (const auto &sc : _ScHips)
        maxSc = std::max(maxSc, static_cast<int>(sc));
    if (maxSc > 1.0)
        _nLevels += ceil(log(maxSc) / log(4));

    _Nm1 = _nLevels - 1;
    _Nm2 = _nLevels - 2;
    _Nm3 = _nLevels - 3;

    _nparcels = static_cast<int>(pow(2, _Nm1));
    _parcelTimes.resize(_nparcels, 0);
    _i_batchelor.resize(_nVar, 0);

    std::vector<double> levelLengths(_nLevels);                                // Including all levels, but last 2 don't count
    std::vector<double> levelTaus(_nLevels);                                   // Smallest scale is 2 levels up from bottom
    _levelRates.resize(_nLevels);

    for (int i = 0; i < _nLevels; ++i) {
        levelLengths[i] = domainLength * pow((ReApproach == "dynamic_A" ? _Anew : _Afac), i);

        levelTaus[i] = tau0 * pow(levelLengths[i] / domainLength, 2.0 / 3.0) / _C_param;
        _levelRates[i] = 1.0 / levelTaus[i] * pow(2.0, i);
    }

    if (ReApproach == "micromixing") {                                          // Adjust rates for micromixing model
        levelTaus[_Nm3] = tau0 * pow(_lStar / domainLength, 2.0 / 3.0) / _C_param;
        _levelRates[_Nm3] = 1.0 / levelTaus[_Nm3] * pow(2.0, _Nm3);
    }

    if (ReApproach == "probability") {                                          // Adjust final mixing rate based on probability
        _levelRates[_Nm3] = _levelRates[_nL - 3] * _probability;
    }

    _LScHips = !_ScHips.empty();                                               // Correct levels for high Sc (levels > Kolmogorov)
    if (_LScHips) {
        for (int i = _iEta + 1; i < _nLevels; ++i) {
            levelTaus[i] = tau0 * pow(levelLengths[_iEta] / domainLength, 2.0 / 3.0) / _C_param;
            _levelRates[i] = 1.0 / levelTaus[i] * pow(2.0, i);
        }
    }

    //-----------------------------------------------------

    _eddyRate_total = 0.0;
    for (int i = 0; i <= _Nm3; ++i)
        _eddyRate_total += _levelRates[i];

    _eddyRate_inertial = 0.0;
    for (int i = 0; i <= _iEta; ++i)
        _eddyRate_inertial += _levelRates[i];

    //-----------------------------------------------------

    _i_plus.resize(_nVar);
    for (int k = 0; k < _nVar; ++k) {
        if (_ScHips[k] < 1.0)
            _i_batchelor[k] = _iEta + 1.5 * log(_ScHips[k]) / log(4);
        else if (_ScHips[k] > 1.0)
            _i_batchelor[k] = _iEta + log(_ScHips[k]) / log(4);
        else
            _i_batchelor[k] = _iEta;
        _i_plus[k] = ceil(_i_batchelor[k]);
    }

    //-----------------------------------------------------

    _varRho.resize(_nparcels);
    _wPar.assign(_nparcels, 1.0 / _nparcels);
    _pLoc.resize(_nparcels);
    for (int i = 0; i < _nparcels; ++i)
        _pLoc[i] = i;

    _currentIndex = 0;
}

void HiPS::set_varData(std::vector<double> &v, std::vector<double> &w, const std::string &varN){
    
    _varData[_currentIndex] = std::make_shared<std::vector<double>>(projection(v, w));
    _varName[_currentIndex] = varN;

    _currentIndex++;
}

void HiPS::set_varData(std::vector<double> &v, std::vector<double> &w, const std::string &varN, const std::vector<double> &rho){

    std::pair<std::vector<double>, std::vector<double>> results = projection(v, w, rho);

    _varData[_currentIndex] = std::make_shared<std::vector<double>>(results.first);
    _varRho = results.second;
    _varName[_currentIndex] = varN;
 
    _currentIndex++;
}

std::vector<double> HiPS::projection(std::vector<double> &vcfd, std::vector<double> &weight){
    
    _xc = setGridCfd(weight);                               // Populate the physical domain for flow particles
    _xh = setGridHips(_nparcels);                            // Populate the physical domain for HiPS parcels

    int nc = _xc.size() - 1;
    int nh = _xh.size() - 1;

    std::vector<double> vh(nh, 0.0);
    int jprev = 0;

    for(int i = 0; i < nh; i++) {
        for(int j = jprev + 1; j <= nc; ++j) {
            if(_xc[j] <= _xh[i + 1]) {
                double d1 = _xc[j] - _xc[j - 1];
                double d2 = _xc[j] - _xh[i];
                // calculation of shortest distance
                // handling if distance is zero but substarction gives comutational error
                double d  = std::min(d1, d2) < 1E-15 ? 0 : std::min(d1, d2);

                vh[i] += vcfd[j - 1] * d;
            } 
            else {
                double d1 = _xh[i + 1] - _xc[j - 1];
                double d2 = _xh[i + 1] - _xh[i];
                // calculation of shortest distance
                // handling if distance is zero but substarction gives comutational error
                double d  = std::min(d1, d2) < 1E-15 ? 0 : std::min(d1, d2);

                vh[i] += vcfd[j - 1] * d;
                jprev = j - 1;
                break;
            }
        }
        vh[i] /= (_xh[i + 1] - _xh[i]);
    }
    return vh;
}

std::pair<std::vector<double>, std::vector<double>>
HiPS::projection(      std::vector<double> &vcfd,
                       std::vector<double> &weight,
                 const std::vector<double> &density){
    // Build CFD and HiPS grids
    _xc = setGridCfd(weight);
    _xh = setGridHips(_nparcels);

    int nc = _xc.size() - 1; // CFD cell count
    int nh = _xh.size() - 1; // HiPS parcel count

    std::vector<double> vh(nh, 0.0);     // parcel-averaged phi
    std::vector<double> rho_h(nh, 0.0);  // parcel-averaged density

    int jprev = 0; // remember where we left off in CFD cells

    for (int i = 0; i < nh; i++) 
    {
        double M = 0.0;     // total mass in this parcel
        double Mphi = 0.0;  // total (mass * phi) in this parcel

        double parcel_length = _xh[i + 1] - _xh[i];

        for (int j = jprev + 1; j <= nc; ++j) 
        {
            // Find geometric overlap between CFD cell and HiPS parcel
            double overlap_start = std::max(_xh[i], _xc[j - 1]);
            double overlap_end   = std::min(_xh[i + 1], _xc[j]);

            double overlap_len = overlap_end - overlap_start;
            if (overlap_len <= 0.0) continue;

            // Effective length scaled by _wPar[i]
            double effective_len = overlap_len * (_wPar[i] / parcel_length);

            // Accumulate mass and mass*phi
            double rho = density[j - 1];
            double phi = vcfd[j - 1];
            M    += rho * effective_len;
            Mphi += rho * phi * effective_len;

            // If CFD cell ends after parcel end  move to next parcel
            if (_xc[j] >= _xh[i + 1]) {
                jprev = j - 1;
                break;
            }
        }

        // Convert total mass  average density and phi
        if (_wPar[i] > 0.0) rho_h[i] = M    / _wPar[i];
        if (       M > 0.0) vh[i]    = Mphi / M;
    }

    return {vh, rho_h};
}

std::vector<double> HiPS::setGridCfd(std::vector<double> &w){

    double sumw = 0.0;
    for(int i=0; i<w.size(); i++)
        sumw += w[i];
   
    std::vector<double> pos;                               // Initializing a vector to hold the grid positions
    double posL = 0.0;                                     // Initializing the starting position

    int i = 0;

    while (i <= w.size()) {                               // Generate the grid positions based on the weights
        pos.push_back(posL);                              // Add the current position to the grid
        posL += w[i]/sumw;                                // Move to the next position by adding the corresponding weight
        i++;                                              
    }
    return pos;                                           // Return the generated grid positions
}

std::vector<double> HiPS::setGridHips(int N){

    std::vector<double> xh(N + 1);                               // Initialize a vector to hold the grid points
    double step = 1.0 / N;                                       // Calculate the step size

    for(int i = 0; i <= N; i++)                                 // Populate the grid with evenly spaced points
        xh[i] = i * step;
        
    return xh;                                                  // Return the generated grid
}

void HiPS::calculateSolution(const double tRun, bool shouldWriteData) {
    
    unsigned long long nEddies = 0;                               // Number of eddy events
    int fileCounter = 0;                                          // Number of data files written
    int iLevel;                                                   // Tree level of EE with top at iLevel=0
    int iTree;                                                    // One of two subtrees involved in swap at iLevel                                           
    _time = 0.0;                                                  // Initialize simulation time
    int lastEddyOutput = 0;                                       // Track last eddy-based output event

    // Apply default values if user hasn't set them
    if (!_useEddyBasedWriting && !_useTimeBasedWriting) {
        _outputIntervalEddy = _DEFAULT_EDDY_INTERVAL;
        _useEddyBasedWriting = true;  // Default to eddy-based writing
    }

    sample_hips_eddy(_dtEE, iLevel);                              // Get first EE at _time 0+_dtEE
    nEddies++;
    _eddyCounter = 0;                                             // Reset eddy counter at start
    _lastOutputTime = 0.0;                                        // Reset last output time

    while (_time + _dtEE <= tRun) {
        _time += _dtEE;
        selectAndSwapTwoSubtrees(iLevel, iTree);
        advanceHips(iLevel, iTree);                              // Reaction and micromixing (if needed) to t=_time

        sample_hips_eddy(_dtEE, iLevel);
        nEddies++;
        _eddyCounter++;

        //  Only check the selected mode (set in example code or by default)
        bool writeByEddy = (_useEddyBasedWriting && _eddyCounter >= lastEddyOutput + _outputIntervalEddy);
        bool writeByTime = (_useTimeBasedWriting && _time - _lastOutputTime >= _outputIntervalTime);

        if (shouldWriteData) {
            if (writeByEddy) {  //  Only write if using eddy-based writing
                writeData(_realization, ++fileCounter, _time);
                lastEddyOutput = _eddyCounter;  //  Update last output event
            }
            else if (writeByTime) {  // Only write if using time-based writing
                 writeData(_realization, ++fileCounter, _time);
                _lastOutputTime = _time;
            }
        }
    }

    // Ensure the final time step completes
    _time = tRun;
    iLevel = 0; 
    iTree = 0;

    if (_performReaction)
        reactParcels_LevelTree(iLevel, iTree);                   // React all parcels up to end time
    saveAllParameters();
}

void HiPS::sample_hips_eddy(double &_dtEE, int &iLevel) {

    static double c1 = 1.0 - pow(2.0, 5.0/3.0*(_iEta+1));
    static double c2 = pow(2.0, _Nm2) - pow(2.0, _iEta+1);
    static double c3 = pow(2.0, _iEta+1);

    //--------------- time to next eddy

    double r = rand.getRand();
    _dtEE = -log(r)/_eddyRate_total;

    //----------------- get eddy level

    r = rand.getRand();

    if ( r <= _eddyRate_inertial/_eddyRate_total) {     // Inertial region
        r = rand.getRand();
        iLevel = ceil(3.0/5.0*log2(1.0-r*c1) - 1.0);
        if (iLevel < 0)    iLevel = 0;
        if (iLevel > _iEta) iLevel = _iEta;
    }

    else {                                            // "Batchelor" region
        r = rand.getRand();
        iLevel = ceil(log2(r*c2 + c3) - 1.0);
        if (iLevel < _iEta+1) iLevel = _iEta+1;
        if (iLevel > _Nm3) iLevel = _Nm3;
    }
    return;
}

void HiPS::selectAndSwapTwoSubtrees(const int iLevel, int &iTree){

    iTree = rand.getRandInt((1 << iLevel)-1);
    int zero_q = rand.getRandInt(1);                                    // 0q where q is 0 or 1
    int one_r  = 2 + rand.getRandInt(1);                                // 1r where r is 0 or 1

    int Qstart = (zero_q << (_Nm3-iLevel)) + (iTree << (_Nm1-iLevel));     // starting index of Q parcels
    int Rstart = (one_r  << (_Nm3-iLevel)) + (iTree << (_Nm1-iLevel));     // starting index of R parcels
    int nPswap = 1 << (_Nm3-iLevel);                                      // number of parcels that will be swapped

    int Qend = Qstart + nPswap;                                          // inclusive indices are Qstart to Qend-1
    int Rend = Rstart + nPswap;                                          // inclusive indices are Rstart to Rend-1
    vector<int> aa(_pLoc.begin()+Qstart, _pLoc.begin()+Qend);
    copy(_pLoc.begin()+Rstart, _pLoc.begin()+Rend, _pLoc.begin()+Qstart); // python: _pLoc[Qstart:Qend]=_pLoc[Rstart:Rend]
    copy(aa.begin(), aa.end(), _pLoc.begin()+Rstart);                     // python: _pLoc[Rstart:Rend]=aa
}

void HiPS::advanceHips(const int iLevel, const int iTree){

    if (_forceTurb && iLevel == 0) {
        forceProfile();                                                  // Forcing for statistically stationary
    }

    bool rxnDone = false;                                               // React all variables once
    for (int k = 0; k < _nVar; k++) {                                    // Upon finding first variable needing micromixing
        // Combined condition check with approach condition
        if ((iLevel >= _i_plus[k]) ||
            (iLevel == _i_plus[k] - 1 && rand.getRand() <= _i_plus[k] - _i_batchelor[k])) {
                if (!rxnDone && _performReaction) {
                    reactParcels_LevelTree(iLevel, iTree);
                    rxnDone = true;
                }
                mixAcrossLevelTree(k, iLevel, iTree);
        }
    }
}

int HiPS::getVariableIndex(const std::string &varName) const{

    auto it = std::find(_varName.begin(), _varName.end(), varName);
    if (it == _varName.end()) {
        throw std::runtime_error("Error: Variable name '" + varName + "' not found.");
    }
    return std::distance(_varName.begin(), it);
}

void HiPS::reactParcels_LevelTree(const int iLevel, const int iTree){

  #ifdef REACTIONS_ENABLED
    // ---- cache indices once
    const int enthalpyIdx = getVariableIndex("enthalpy");
    std::vector<int> yIdx(_nsp);
    for (int k = 0; k < _nsp; ++k) {
        yIdx[k] = getVariableIndex(_gas->speciesName(k)); // map species name -> var index
    }

    const int nP     = 1 << (_Nm1 - iLevel);
    const int istart = iTree * nP;
    const int iend   = istart + nP;

    std::vector<double> y(_nsp);

    for (int i = istart; i < iend; ++i) {
        const int ime = _pLoc[i];
        const double dt = _time - _parcelTimes[ime];

        // Pull current state
        double h = (*_varData[enthalpyIdx])[ime];
        for (int k = 0; k < _nsp; ++k) {
            y[k] = (*_varData[yIdx[k]])[ime];
        }

        // ---- store old density BEFORE chemistry
        const double rho_old = _varRho[ime];

        if (_performReaction) {
            // Advance chemistry; _bRxr updates its state internally
            _bRxr->react(h, y, dt);

            // Get new density from the reactor/EOS
            const double rho_new = _bRxr->getDensity();

            // Keep per-parcel mass m = rho * _wPar * V_tot constant
            if (rho_new > 0.0) {
                _wPar[ime]  *= (rho_old / rho_new);
                _varRho[ime] =  rho_new;
            }

        }

        // Write back enthalpy and species (post-reaction)
        (*_varData[enthalpyIdx])[ime] = h;
        for (int k = 0; k < _nsp; ++k) {
            (*_varData[yIdx[k]])[ime] = y[k];
        }

        _parcelTimes[ime] = _time;
    }

 #endif
}


void HiPS::mixAcrossLevelTree(int kVar, const int iLevel, const int iTree){
    
    int istart;
    int iend;

    int nPmix = 1 << (_nLevels - iLevel - 2);   // Number of parcels mixed together
    int ime;

    //---------- Mix left branch of iTree ----------
    istart = iTree << (_Nm1 - iLevel);
    iend   = istart + nPmix;

    double s    = 0.0;  // sum (value or value*mass)
    double msum = 0.0;  // mass sum (only used if weighted)

    for (int i = istart; i < iend; i++) {
        ime = _pLoc[i];
        if (_performReaction) {
            double m = _varRho[ime] * _wPar[ime];
            s    += (*_varData[kVar])[ime] * m;
            msum += m;
        } else {
            s += (*_varData[kVar])[ime];
        }
    }

    double avg = _performReaction
               ? ((msum > 0.0) ? (s / msum) : 0.0)   
               : (s / nPmix);

    for (int i = istart; i < iend; i++) {
        ime = _pLoc[i];
        (*_varData[kVar])[ime] = avg;
    }

    //---------- Mix right branch of iTree ----------
    istart = iend;
    iend   = istart + nPmix;

    s    = 0.0;
    msum = 0.0;

    for (int i = istart; i < iend; i++) {
        ime = _pLoc[i];
        if (_performReaction) {
            double m = _varRho[ime] * _wPar[ime];
            s    += (*_varData[kVar])[ime] * m;
            msum += m;
        } else {
            s += (*_varData[kVar])[ime];
        }
    }

    avg = _performReaction
        ? ((msum > 0.0) ? (s / msum) : 0.0)
        : (s / nPmix);

    for (int i = istart; i < iend; i++) {
        ime = _pLoc[i];
        (*_varData[kVar])[ime] = avg;
    }
}

void HiPS::forceProfile(){
    // Loop through each variable in the HiPS profile
    for (int k = 0; k < _varData.size(); k++) {
        double s=0;                                                  // Temporary variable for summation

        //---------- Force the left half of parcels to average 0 ----------

        for (int i = 0; i < _nparcels >> 1; i++)
            s += (*_varData[k])[_pLoc[i]];                             // Calculate the sum of values in the left half of parcels
        
        s /= (_nparcels >> 1); // Calculate the average of values in the left half of parcels
        
        for (int i = 0; i < _nparcels >> 1; i++)
            (*_varData[k])[_pLoc[i]] += (-s - 0.0);                    // Adjust values in the left half of parcels to achieve an average of 0

        //---------- Force the right half of parcels to average 1 ----------
        s = 0.0;

        for (int i = _nparcels >> 1; i < _nparcels; i++)
            s += (*_varData[k])[_pLoc[i]];                             // Calculate the sum of values in the right half of parcels
        
        s /= (_nparcels >> 1);                                        // Calculate the average of values in the right half of parcels
        
        for (int i = _nparcels >> 1; i < _nparcels; i++)
            (*_varData[k])[_pLoc[i]] += (-s + 1.0);                    // Adjust values in the right half of parcels to achieve an average of 1
    }
}

void HiPS::writeData(int real, const int ifile, const double outputTime){

    stringstream ss1, ss2;
    string s1, s2;

    // Prepare directory path
    ss1 << "../data/rlz_" << setfill('0') << setw(5) << real;
    ss1 >> s1;

    // Create directories if they don't exist
    #ifdef _WIN32
        system(("mkdir " + s1).c_str());
    #else
        system(("mkdir -p " + s1).c_str());
    #endif

    ss2 << "Data_" << setfill('0') << setw(5) << ifile << ".dat";
    ss2 >> s2;

    string fname = s1 + "/" + s2;
    ofstream ofile(fname.c_str());
    cout << endl << "writing data for time " << outputTime << " to file: " << fname.c_str();

    // Check if file opened successfully
    if (!ofile) {
        cerr << "Error: Unable to open file " << fname << " for writing!" << endl;
        return;
    }

    // Write metadata (header information)
    ofile << "# time = " << outputTime << "\n";
    ofile << "# Grid Points = " << _nparcels << "\n";
        
    // Write column names (include temperature if reactions are enabled)
    if(_performReaction)
        ofile << setw(19) << "# Temp";  // Include temperature column if reactions are enabled

    for (const auto& varN : _varName) {
        ofile << setw(19) << "# " << varN;
    }
    ofile << endl;  // End of the header line

    // Set scientific notation and precision
    ofile << scientific;
    ofile << setprecision(10);

    // Write data
    for (int i = 0; i < _nparcels; i++) {
        // Write temperature first if reactions are enabled
        #ifdef REACTIONS_ENABLED
        if(_performReaction) {
            vector<double> yy(_nsp);
            for(int k=0; k<_nsp; k++)
                yy[k] = (*_varData[k+1])[_pLoc[i]];
            _gas->setMassFractions(yy.data());
            _gas->setState_HP((*_varData[0])[_pLoc[i]], _gas->pressure());
            ofile << setw(19) << _gas->temperature();
        }
        #endif

        // Write variables
        for (int k = 0; k < _nVar; k++) {
            ofile << setw(19) << (*_varData[k])[_pLoc[i]];
        }
        ofile << endl;
    }

    ofile.close();
   // cout << "Data successfully written to: " << fname << endl;
}
   
std::vector<double> HiPS::projection_back(std::vector<double> &vh) {

    int nh = _xh.size() - 1;
    int nc = _xc.size() - 1;

    std::vector<double> vc(nc, 0.0);
    int jprev = 0;

    for (int i = 0; i < nc; ++i) {
        for (int j = jprev + 1; j <= nh; ++j) {
            if (_xh[j] <= _xc[i + 1]) {
                double d1 = _xh[j] - _xh[j - 1];
                double d2 = _xh[j] - _xc[i];
                // calculation of shortest distance
                // handling if distance is zero but substarction gives comutational error
                double d  = std::min(d1, d2) < 1E-15 ? 0 : std::min(d1, d2);

                vc[i] += vh[j - 1] * d;
            } else {
                double d1 = _xc[i + 1] - _xh[j - 1];
                double d2 = _xc[i + 1] - _xc[i];
                // calculation of shortest distance
                // handling if distance is zero but substarction gives comutational error
                double d  = std::min(d1, d2) < 1E-15 ? 0 : std::min(d1, d2);

                vc[i] += vh[j - 1] * d;

                jprev = j - 1;
                break;
            }
        }
        vc[i] /= (_xc[i + 1] - _xc[i]);
    }
    return vc;
}



std::vector<double> HiPS::projection_back_with_density(std::vector<double> &vh,
                                                       std::vector<double> &rho_h,
                                                       std::vector<double> &rho_c){
    const int nh = static_cast<int>(_xh.size()) - 1;  // # HiPS parcels
    const int nc = static_cast<int>(_xc.size()) - 1;  // # CFD cells

    std::vector<double> phi_c(nc, 0.0);   // CFD-side variable (output)
    std::vector<double> M_cfd(nc, 0.0);   // mass in each CFD cell
    std::vector<double> Mphi_cfd(nc, 0.0);// mass*phi in each CFD cell

    int jprev = 0;
    for (int i = 0; i < nc; ++i) {
        // cell geometry (normalized length, since _xc is 0..1)
        const double cell_len = _xc[i + 1] - _xc[i];

        for (int j = jprev + 1; j <= nh; ++j) {
            // geometric overlap between parcel j-1 and CFD cell i
            const double overlap_start = std::max(_xh[j - 1], _xc[i]);
            const double overlap_end   = std::min(_xh[j],     _xc[i + 1]);
            const double overlap_len   = overlap_end - overlap_start;
            if (overlap_len <= 0.0) continue;

            // scale overlap by current parcel volume fraction _wPar
            const double parcel_len    = _xh[j] - _xh[j - 1]; // = 1.0/nh
            const double effective_len = overlap_len * (_wPar[_pLoc[j - 1]] / parcel_len);

            // accumulate mass and mass*phi into CFD cell i
            const double rhoP = rho_h[j - 1];
            const double phiP = vh[j - 1];
            M_cfd[i]    += rhoP * effective_len;
            Mphi_cfd[i] += rhoP * phiP * effective_len;

            // advance parcel index if we reached end of this CFD cell
            if (_xc[i + 1] <= _xh[j]) { jprev = j - 1; break; }
        }

        // recover CFD density and variable
        rho_c[i] = (cell_len > 0.0)        ? (M_cfd[i] / cell_len) : 0.0;
        phi_c[i] = (M_cfd[i] > 1e-300)     ? (Mphi_cfd[i] / M_cfd[i]) : 0.0;
    }
    return phi_c;
}

std::vector<std::vector<double>> HiPS::get_varData(){
    std::vector<std::vector<double>> varDataProjections;

    for (int i = 0; i < _varData.size(); i++) {

        // Reorder using _pLoc
        std::vector<double> vh(_nparcels);
        for (int j = 0; j < _nparcels; j++) {
            vh[j] = (*_varData[i])[_pLoc[j]];
        }

        std::vector<double> vc = projection_back(vh);
        varDataProjections.push_back(vc);
    }

    return varDataProjections;
}

std::pair< std::vector<std::vector<double>>, std::vector<double>>
HiPS::get_varData_with_density(){
    std::vector<std::vector<double>> varDataProjections;
    
    // Reorder _varRho based on _pLoc
    std::vector<double> rho_h(_nparcels);
    for (int i = 0; i < _nparcels; i++) {
        rho_h[i] = _varRho[_pLoc[i]];
    }

    std::vector<double> rho_c = projection_back(rho_h);

    for (int i = 0; i < _varData.size(); i++) {

        // Reorder variable data using _pLoc
        std::vector<double> vh(_nparcels);
        for (int j = 0; j < _nparcels; j++) {
            vh[j] = (*_varData[i])[_pLoc[j]];
        }

        auto vc = projection_back_with_density(vh, rho_h, rho_c);
        varDataProjections.push_back(vc);
    }

    return {varDataProjections, rho_c};
}

void HiPS::setOutputIntervalEddy(int interval){

    _outputIntervalEddy = interval;
    _useEddyBasedWriting = true;  ///< Enables eddy-based writing
    _useTimeBasedWriting = false; ///< Disables time-based writing
}

void HiPS::setOutputIntervalTime(double interval){

    _outputIntervalTime = interval;
    _useTimeBasedWriting = true;  ///< Enables time-based writing
    _useEddyBasedWriting = false; ///< Disables eddy-based writing
}


void HiPS::saveAllParameters(){

    std::string filepath = "../post/parameters.dat";  ///< Output file path for simulation parameters
    std::ofstream file(filepath);

    if (!file) {
        std::cerr << "Error: Could not open " << filepath << " for writing!\n";
        return;
    }

    // Write user-defined input parameters
    file << "nLevels " << _nLevels << "\n";            ///< Number of hierarchical levels in the HiPS model
    file << "domainLength " << _domainLength << "\n";  ///< Length of the computational domain
    file << "tau0 " << _tau0 << "\n";                  ///< Reference eddy turnover time
    file << "C_param " << _C_param << "\n";            ///< Model constant controlling turbulence behavior
    file << "forceTurb " << _forceTurb << "\n";        ///< Flag for forced turbulence (1 = enabled, 0 = disabled)
    file << "nVar " << _nVar << "\n";                  ///< Number of variables tracked in the simulation
    file << "performReaction " << _performReaction << "\n";  ///< Flag indicating whether chemical reactions are simulated
    file << "realization " << _realization << "\n";    ///< Current simulation realization (for multiple runs)

    // Write variable names
    if (!_varName.empty()) {
        file << "varName ";
        for (const std::string &name : _varName) {
            file << name << " ";  ///< Separate variable names by spaces
        }
        file << "\n";
    } else {
        file << "varName (undefined)\n";
    }

    // Write _i_batchelor vector if it exists and has the correct size
    if (!_i_batchelor.empty() && _i_batchelor.size() == _nVar) {
        file << "i_batchelor ";
        for (const auto &val : _i_batchelor) {
            file << val << " ";  ///< Print each element separated by space
        }
        file << "\n";
    } else {
        file << "i_batchelor (undefined or size mismatch)\n";
    }

    file.close();
    //cout << endl << "All parameters saved in: " << filepath << std::endl;
}