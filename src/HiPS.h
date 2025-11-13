#pragma once

#ifdef  REACTIONS_ENABLED
#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "BatchReactor.h"
#endif

#include <iostream>
#include <vector>
#include <string>
#include <vector>
#include <memory>
#include "RandomGenerator.h"


/// \brief Implementation of the HiPS (Hierarchical Parcel Swapping) model.
///
/// This file contains the core implementation of the HiPS model, used for
/// simulating turbulent mixing. It supports optional functionality for chemical
/// reactions, which can be enabled by defining the `REACTIONS_ENABLED` macro.
///
/// Dependencies:
/// - YAML-CPP: Used for reading and parsing configuration files.
/// - Batch reactor classes (`BatchReactor_cvode.h`, `BatchReactor_cantera.h`):
///   Included only when `REACTIONS_ENABLED` is defined, enabling reaction
///   modeling.
/// - Standard C++ libraries: Utilized for input/output, mathematical computations,
///   and data handling.
///
/// \note Define `REACTIONS_ENABLED` at compile time to enable chemical reaction
/// functionality.
class HiPS {

public:
  
    int realization;                                               ///< Number of realizations
    std::vector<std::shared_ptr<std::vector<double>>> varData;     ///< Vector of pointers to vector
    std::vector<double> varRho;                                    ///< Density
    std::vector<int> pLoc;                                         ///< Parcel index array for fast implementation of swaps
    std::vector<double> wPar;                                      ///< Parcel volume fractions


#ifdef REACTIONS_ENABLED
    std::shared_ptr<Cantera::ThermoPhase> gas;                     ///< Shared pointer to a Cantera thermochemistry object
    std::shared_ptr<BatchReactor> bRxr;                            ///< Unique pointer to the integrator object
#endif

    double domainLength;                                           ///< Length of domain (m)
    double tau0;                                                   ///< Integral timescale
    double C_param;                                                ///< Eddy frequency parameter
    
private:

    int nparcels;                                                  ///< Number of parcels
    int currentIndex = 0;                                          ///< Member variable to keep track of current index of variables
    int nLevels;                                                   ///< Number of tree levels
    int nVar;                                                      ///< Number of parcel variables (e.g., h, ysp)
    int nsp;                                                       ///< Number of species
    int Nm1;                                                       ///< nLevels - 1 
    int Nm2;                                                       ///< nLevels - 2                     
    int Nm3;                                                       ///< nLevels - 3 
    int iEta;                                                      ///< Kolmogorov level (needed for variable Sc scalars)
    int nL;                                                        ///< Adjusted number of levels based on the Reynolds number
    bool forceTurb;                                                ///< Forcing function for statistically stationary: -1 = none, 1 = source term, 2 = dir

    bool LScHips;                                                  ///< HiPS schmidt number
    bool performReaction;                                          ///< Flag indicating whether chemical reactions are performed in the simulation
        
    double time;                                                   ///< Current simulation time
    double eddyRate_total;                                         ///< Total rate of all eddies 0 through nLevels-3
    double eddyRate_inertial;                                      ///< Total rate of all eddies 0 through iEta (= eddyRate_total if Sc=1)
    double Afac = 0.5;                                             ///< Level lengthscale reduction factor (0.5)
    double Re;                                                     ///< Reynolds number
    double dtEE;                                                   ///< Time increment to next eddy event
    double Prob;                                                   ///< Probability value for probability-based solution
    double lStar;                                                  ///< Length of the level associated with the Reynolds number
    double Anew;                                                   ///< Adjusted level length scale reduction factor for dynamic adjustment of reduction factor

    RandomGenerator rand;
    
    std::vector<int> i_plus;                                       ///< ceil(i_batchelor)
    std::vector<double> ScHips;                                    ///< Vector containing Schmidt numbers related to each variable
    std::vector<std::string> varName;                              ///< Vector containing the names of parcel variables
    std::vector<double> parcelTimes;                               ///< Current times corresponding to the parcel states
    std::vector<double> levelRates;                                ///< List of eddy event rates at each level
    std::vector<double> i_batchelor;                               ///< Batchelor level for variable Sc scalars; NOTE: double, as in, between levels
    std::vector<double> xc;                                        ///< Vector containing physical domain of flow particles
    std::vector<double> xh;                                        ///< Vector containing physical domain of HiPS parcels
    
    std::string  ReApproach;
    int outputIntervalEddy = 10;                                   ///< Default: write data every 10 eddy events
    double outputIntervalTime = 0.1;                               ///< Default: write data every 0.1s
    int eddyCounter = 0;                                           ///< Counter for eddy events
    double lastOutputTime = 0.0;                                   ///< Last time data was written
    bool useEddyBasedWriting = false;                              ///< Tracks if eddy writing is set
    bool useTimeBasedWriting = false;                              ///< Tracks if time writing is set
    const int DEFAULT_EDDY_INTERVAL = 1000;                        ///< Default: Write every 1000 eddies
    const double DEFAULT_TIME_INTERVAL = 0.1;                      ///< Default: Write every 0.1s

public:

    ////////////////////////////////////////////////////////////////////////////
    /// \brief Constructor for initializing a HiPS object without building the
    /// tree immediately.
    ///
    /// This constructor is designed for simulations where the HiPS tree must be
    /// configured dynamically, such as grid-based simulations with
    /// cell-specific turbulence properties. The tree can be created or updated
    /// later using one of the `set_tree` functions.
    /// \param C_param_         Eddy coefficient controlling mixing rate.
    /// \param forceTurb_       Flag to enforce turbulence activation.
    /// \param nVar_            Number of transported variables.
    /// \param performReaction_ Enables chemical reactions if set to true.
    /// \param cantSol          Cantera solution object (required when
    ///                         REACTIONS_ENABLED).
    /// \param seed             Random seed (negative for random initialization).
    /// \param realization_     Realization index for ensemble or parallel runs.
    ///
    /// \note This constructor does not call any `set_tree()` function. The user
    ///       is responsible for building the tree explicitly by calling either
    ///       version of `set_tree(...)`.
    ///
    /// \note See the full constructor for chemical integrator behavior under
    ///       `REACTIONS_ENABLED`.
    ///
    /// \see HiPS::set_tree() for deferred tree construction.
    ////////////////////////////////////////////////////////////////////////////
    HiPS(double C_param_,
         bool forceTurb_,
         int nVar_,
         std::vector<double> &ScHips_,
         bool performReaction,
         std::shared_ptr<void> vcantSol = nullptr,
         int seed = 10,
         int realization_ = 1);

    ////////////////////////////////////////////////////////////////////////////
    /// \brief Constructor for initializing the full HiPS tree at the time of
    /// object creation.
    ///
    /// This constructor is intended for simulations where the domain structure
    /// and turbulence properties are fixed throughout the run (e.g., standalone
    /// mixing or reacting cases). The HiPS tree is initialized immediately
    /// using the provided parameters, and the structure remains constant
    /// throughout the simulation.
    ///
    /// \param nLevels          Number of levels in the HiPS binary tree (can be
    ///                         adjusted for high Sc).
    /// \param domainLength_    Domain length for defining spatial scales.
    /// \param tau0_            Characteristic time scale for the smallest eddy.
    /// \param C_param_         Eddy coefficient controlling mixing rate.
    /// \param forceTurb_       Flag to enforce turbulence activation.
    /// \param nVar_            Number of transported variables.
    /// \param ScHips_          Vector of Schmidt numbers (one per variable).
    /// \param performReaction_ Enables chemical reactions if set to true.
    /// \param cantSol          Cantera solution object (required when
    ///                         REACTIONS_ENABLED).
    /// \param seed             Random seed (negative for random initialization).
    /// \param realization_     Realization index for ensemble or parallel runs.
    ///
    /// \note This constructor calls `set_tree(nLevels, domainLength, tau0,
    ///       ScHips)` internally to fully build the tree at initialization
    ///       time. This setup is optimal for cases where tree reconfiguration
    ///       is not needed during runtime.
    ///
    /// \note If `REACTIONS_ENABLED` is defined, the default chemical integrator
    ///       is `BatchReactor_cvode`. To use `BatchReactor_cantera` instead,
    ///       uncomment the corresponding line in the constructor code.
    ///
    /// \see HiPS::set_tree() for the internal tree setup logic.
    ////////////////////////////////////////////////////////////////////////////
    HiPS(int nLevels,
         double domainLength_,
         double tau0_,
         double C_param_,
         bool forceTurb_,
         int nVar_,
         std::vector<double> &ScHips_,
         bool performReaction,
         std::shared_ptr<void> vcantSol = nullptr,
         int seed = 10,
         int realization_ = 1);


    ////////////////////////////////////////////////////////////////////////////
    /// \brief Sets up the HiPS tree using explicitly specified tree parameters.
    ///
    /// This method builds the binary HiPS tree using a user-defined number of
    /// levels and physical parameters. It also automatically adjusts the number
    /// of levels to account for high Schmidt numbers, ensuring that micromixing
    /// is properly resolved.
    ///
    /// \param nLevels          Number of levels in the HiPS tree.
    /// \param domainLength_    Domain size for determining eddy length scales.
    /// \param tau0_            Time scale of the smallest eddy (Kolmogorov scale).
    /// \param ScHips_          Vector of Schmidt numbers (one per variable).
    ///
    /// \note This function is used by the full constructor and can also be
    ///       called manually after using the dynamic constructor.
    ///
    /// \warning The number of levels may be increased automatically for large
    ///          Schmidt numbers to ensure accurate scalar mixing across scales.
    ////////////////////////////////////////////////////////////////////////////
    void set_tree(int nLevels, double domainLength_, double tau0_);


    ////////////////////////////////////////////////////////////////////////////
    /// \brief Dynamically builds the HiPS tree based on Reynolds number and a
    /// selected strategy.
    ///
    /// This method configures the HiPS tree using a continuous Reynolds number
    /// and one of several initialization strategies. It is intended for
    /// simulations where local turbulence conditions change over time or space,
    /// requiring the tree to be updated dynamically.
    ///
    /// \param Re_              Reynolds number used to determine base tree level.
    /// \param domainLength_    Domain length for spatial scaling.
    /// \param tau0_            Base time scale for the largest eddy.
    /// \param ScHips_          Vector of Schmidt numbers (one per variable).
    /// \param ReApproach_      Strategy to convert continuous Re to tree level:
    ///                         - "rounding"     → Round to nearest discrete level
    ///                         - "probability"  → Use probabilistic interpolation between levels
    ///                         - "micromixing"  → Use fixed level with adjusted mixing rate
    ///                         - "dynamic_A"    → Adjust geometric scale factor A to fit Re
    ///
    /// \note This method is ideal for Lagrangian simulations using grid cells
    ///       with different Re values. It supports runtime reconfiguration of
    ///       the tree without reinitializing the HiPS object.
    ///
    /// \warning Ensure consistent `ReApproach_` handling across the simulation
    ///          to avoid inconsistencies.
    ///
    /// \see HiPS::set_tree(int, ...) for direct-level setup.
    ////////////////////////////////////////////////////////////////////////////
    void set_tree(double Re_, double domainLength_, double tau0_, std::string ReApproach_ = "rounding");
  
    ////////////////////////////////////////////////////////////////////////////
    /// \brief Assigns variables, their corresponding weights, and names to the
    /// parcels in the HiPS tree.
    ///
    /// This function projects the provided variables onto the parcels within
    /// the HiPS tree structure, using their corresponding weights. It ensures
    /// that each parcel in the tree is associated with the correct variable and
    /// weight, along with a descriptive name for better identification.
    ///
    /// \param v        Vector of variables to be assigned to the parcels in the
    ///                 HiPS tree.
    /// \param w        Vector of weights corresponding to each variable or
    ///                 parcel.
    /// \param varN     String representing the name of the variable being
    ///                 assigned.
    ///
    /// \note The size of `v` and `w` must match to ensure a one-to-one
    ///       correspondence between variables and their weights. The function
    ///       does not perform size validation internally.
    ///
    /// \warning Ensure that the weights in `w` are normalized or appropriately
    ///          scaled, as they directly influence the projection and
    ///          subsequent simulations.
    ////////////////////////////////////////////////////////////////////////////
    void set_varData(      std::vector<double> &v,
                           std::vector<double> &w,
                     const std::string &varN);

    ////////////////////////////////////////////////////////////////////////////
    /// \brief Assigns variables, weights, names, and densities to the parcels
    /// in the HiPS tree.
    ///
    /// This overloaded function assigns the specified variables, along with
    /// their associated weights, names, and densities, to the parcels within
    /// the HiPS tree structure. The function incorporates particle density
    /// during the projection process, ensuring a more accurate representation
    /// of parcel properties in the simulation.
    ///
    /// \param v        Vector of variables to be assigned to the HiPS tree.
    /// \param w        Vector of weights corresponding to each variable or
    ///                 parcel.
    /// \param varN     String representing the name of the variable being
    ///                 assigned.
    /// \param rho      Vector of densities corresponding to each flow particle.
    ///
    /// \note This function is specifically overloaded to account for particle
    ///       density, which enhances the accuracy of the parcel projection.
    ///       Ensure that the size of `v`, `w`, and `rho` are consistent.
    ///
    /// \warning The values in `rho` should be physically meaningful and
    ///          consistent with the simulation's requirements. Improper density
    ///          values may lead to inaccuracies or instabilities in the
    ///          simulation.
    ///          ///////////////////////////////////////////////////////////////
    void set_varData(      std::vector<double> &v,
                           std::vector<double> &w,
                     const std::string &varN,
                     const std::vector<double> &rho);

    ////////////////////////////////////////////////////////////////////////////
    /// \brief Retrieves the final data from the simulation.
    ///
    /// This function returns the final state of the simulation, packaged as a
    /// vector of vectors. It is particularly useful when integrating HiPS as a
    /// subgrid model in CFD simulations, enabling seamless transfer of data for
    /// further analysis or post-processing.
    ///
    /// Projects data from HiPS varData (size nparcels) to number of CFD
    /// particles corresponding to initial set_varData call
    ///
    /// \return A vector of vectors containing the final results, where:
    ///         - Each inner vector represents a specific variable or property.
    ///         - The outer vector contains all variables across parcels.
    ///
    /// \note
    /// - This function is specifically designed for use when HiPS is employed
    ///   as a subgrid model in CFD simulations.
    /// - The structure of the returned data ensures compatibility with CFD
    ///   solvers that require parcel-level data.
    ///
    /// \warning
    /// - Ensure the simulation has reached completion before calling this
    ///   function to avoid incomplete or inconsistent data.
    /// - The returned data structure should be interpreted according to the
    ///   simulation setup and variable ordering.
    ////////////////////////////////////////////////////////////////////////////
    std::vector<std::vector<double>> get_varData();

    std::pair<std::vector<std::vector<double>>, std::vector<double>> get_varData_with_density();
   
    void setOutputIntervalTime(double interval);
    void setOutputIntervalEddy(int interval);

    void calculateSolution(const double tRun, bool shouldWriteData =false);                         // Running simulations 

    void writeData(int real, const int ifile, const double outputTime);                                       // Writing the results for a user-defined number of eddies in the data folder.
    int get_nparcels() const { return nparcels; }
    
    const std::vector<int>& get_pLoc() const { return pLoc; }
    const std::vector<std::shared_ptr<std::vector<double>>>& get_HipsVarData_ptr() const { return varData; } // internal HiPS varData, sized to nparcels

private:

    ////////////////////////////////////////////////////////////////////////////
    /// \brief Projects values from flow particles onto HiPS parcels assuming
    /// constant density.
    ///
    /// This function maps the values of flow particles onto HiPS parcels under
    /// the assumption of constant density. The projection is performed using
    /// the following equation:
    /// \f[
    /// \sum_{i=0}^{\text{Number of Flow Particles (FP)}} (\phi_{\text{FP}} \, \mathrm{d}x_{\text{FP}})_{i} =
    /// \sum_{j=0}^{\text{Number of HiPS Parcels (HP)}} (\phi_{\text{HP}} \, \mathrm{d}x_{\text{HP}})_{j}
    /// \f]
    /// This ensures conservation of properties such as mass or concentration
    /// during the projection.
    ///
    /// \param vcfd     Vector of variables from flow particles to be mapped to
    ///                 HiPS parcels.
    ///
    /// \param weight   Vector of weights, with one weight assigned to each
    ///                 flow particle.
    ///
    /// \return         Vector of projected values for HiPS parcels.
    ///
    /// \note The function assumes constant density throughout the domain. For
    ///       cases with varying density, use an appropriate overloaded function
    ///       or method.
    ///
    /// \warning Ensure that the `vcfd` and `weight` vectors have matching
    ///          sizes, as any discrepancy may result in undefined behavior or
    ///          incorrect projections.
    ////////////////////////////////////////////////////////////////////////////
    std::vector<double> projection(std::vector<double> &vcfd, std::vector<double> &weight);
    
    std::pair<std::vector<double>, std::vector<double>>  projection(std::vector<double> &vcfd, std::vector<double> &weight,                     
                                   const std::vector<double> &density);                             // Perform vector projection flow particles onto HiPS parcels operation with density
    
    std::vector<double> setGridHips(int N);                                                         // Set HiPS grid with a specified number of grid points equal to number of parcels
    std::vector<double> setGridCfd(std::vector<double> &w);                                         // Set CFD grid using provided weight vector
    std::vector<double> projection_back(std::vector<double> &vb);                       

    ////////////////////////////////////////////////////////////////////////////
    /// \brief Projects HiPS parcel values and densities back onto the flow
    /// particles (CFD cells).
    ///
    /// This function reverses the projection process, redistributing both the
    /// values and densities stored in the HiPS parcels back to the flow
    /// particles. It ensures conservation of both the property values and
    /// densities, maintaining consistency between the HiPS parcels and flow
    /// particles.
    ///
    /// ### Conservation Principle:
    /// The projection follows the equation:
    /// \f[
    /// \sum_{j=0}^{\text{Number of HP}} (\phi_{\text{HP}} \, \rho_{\text{HP}} \, \mathrm{d}x_{\text{HP}})_{j} =
    /// \sum_{i=0}^{\text{Number of FP}} (\phi_{\text{FP}} \, \rho_{\text{FP}} \, \mathrm{d}x_{\text{FP}})_{i}
    /// \f]
    /// where:
    /// - \f$\phi_{\text{HP}}\f$: Values in HiPS parcels.
    /// - \f$\rho_{\text{HP}}\f$: Densities in HiPS parcels.
    /// - \f$\mathrm{d}x_{\text{HP}}\f$: Differential volume elements for HiPS parcels.
    /// - \f$\phi_{\text{FP}}\f$: Values in flow particles.
    /// - \f$\rho_{\text{FP}}\f$: Densities in flow particles.
    /// - \f$\mathrm{d}x_{\text{FP}}\f$: Differential volume elements for flow particles.
    ///
    /// \note In the new implementation, geometric overlaps are scaled by
    ///       \f$w_{\text{Par}}/\mathrm{d}x_{\text{HP}}\f$ so that mass is
    ///       preserved when parcel volumes change during chemistry.
    ///
    /// \param vh           Vector of values from HiPS parcels to be projected back.
    /// \param rho_h        Vector of density values from HiPS parcels.
    /// \param rho_c        (output) Vector to receive the densities redistributed to the flow particles.
    /// \return             A vector containing the values projected back onto the flow particles.
    ///
    /// \note
    /// - This function is the reverse of the projection function that includes
    ///   density.
    /// - The input vectors \p vh and \p rho_h should be consistent with the
    ///   HiPS parcel structure and sizes.
    ///
    /// \warning
    /// - Ensure that the HiPS parcels are populated with valid values and
    ///   densities before invoking this function.
    /// - Mismatches in data sizes between HiPS parcels and flow particles may
    ///   lead to inaccurate results.
    ////////////////////////////////////////////////////////////////////////////
    std::vector<double> projection_back_with_density(std::vector<double> &vh,
                                                     std::vector<double> &rho_h,
                                                     std::vector<double> &rho_c);
    



    void sample_hips_eddy(double &dt, int &iLevel);                                                 // Sample HiPS eddy with specified time step and level
    void selectAndSwapTwoSubtrees(const int iLevel, int &iTree);                                    // Select and swap two subtrees in the level tree
    void advanceHips(const int iLevel, const int iTree);                                            // Advancing simulations to do mixing and reaction
   
    int getVariableIndex(const std::string &varName) const;

    void reactParcels_LevelTree(const int iLevel, const int iTree);                                 // Reacting parcels involved in micro-mixing
    void mixAcrossLevelTree(int kVar, const int iMixLevel, const int iTree);                        // Mixing paecels involved in micr0-mixing.
   
    void forceProfile();
    
    void saveAllParameters();                                                                       // Function to save ALL parameters
};