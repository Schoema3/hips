


\mainpage

<!-- #################################################################### -->

# Overview

[HiPS] \cite Kerstein_2013 \cite Kerstein_2014 \cite Kerstein_2021 is an open-source C++ library that functions both as a standalone mixing model and as a sub-grid model for Transported Probability Density Functions.

# Dependencies and Installation

HiPS is designed for Linux-based systems, including macOS and the Linux subsystem for Windows.

## Required Software:
* CMake 3.15+
* C++11

## Optional Software:
* Cantera and Sundials (for reaction support)
* Doxygen (for building documentation)
* Graphviz (for Doxygen)

## Build and Installation Instructions:
1. Create and navigate to a top-level `build` directory.
2. Configure CMake: `cmake ..`
3. Build HiPSlib: `make`
4. Install HiPSlib: `make install`
5. Generate documentation: `make docs`

## CMake Configuration Variables

The default CMake configuration is suitable for users who do not require examples, tests, or documentation immediately. CMake configuration options can be modified by editing the top-level `CMakeLists.txt`, the `CMakeCache.txt` (generated in the `build` directory after running CMake once), or by specifying them on the command line during step 2:

```
cmake -DHIPSLIB_BUILD_EXAMPLES=ON ..
```

The following project-specific CMake variables can be set by the user:

| CMake variable            | Default                       | Description                                                                 |
| ------------------------- | ----------------------------- | --------------------------------------------------------------------------- |
| `CMAKE_INSTALL_PREFIX`     | Top-level project directory    | Installation location                                                       |
| `REACTIONS_ENABLED`        | `OFF`                         | Enable support for chemical reactions                                       |
| `HIPSLIB_BUILD_EXAMPLES`   | `ON`                          | Build HiPS examples                                                         |
| `HIPSLIB_BUILD_DOCS`       | `OFF`                         | Build HiPS documentation with Doxygen                                       |

The `REACTIONS_ENABLED` flag determines if HiPS supports chemical reactions. If set to `ON`, additional libraries like Cantera or Sundials are required. For simple mixing without reactions, set this flag to `OFF`.


# Using HiPS

HiPS consists of three main object classes: `hips`, `batchReactor_cvode`, and `batchReactor_cantera`. These classes act as the computational engine for advancing the reaction system and solving complex differential equations.

## Example Workflow

To use HiPS in C++, include the `hips.h` header file. Interaction with HiPS is primarily done through two constructors:

1. [hips(int nLevels, double domainLength, double tau0, int nVar, int forceTurb, std::vector<double>& ScHips, bool performReaction, std::shared_ptr<Cantera::Solution> cantSol, int seed)](@ref hips(int, double, double, int, int, std::vector<double>&, bool, std::shared_ptr<Cantera::Solution>, int)) for standalone simulations with full parameter initialization.
2. [hips(int nVar, int forceTurb, std::vector<double>& ScHips, bool performReaction, std::shared_ptr<Cantera::Solution> cantSol, int seed)](@ref hips(int, int, std::vector<double>&, bool, std::shared_ptr<Cantera::Solution>, int)) for use as a sub-grid model in CFD simulations, requiring frequent tree structure updates.

To re-initialize the HiPS tree, use the [set_tree(int nLevels, double domainLength, double tau0, std::vector<double>& ScHips)](@ref set_tree(int, double, double, std::vector<double>&)) function. After initialization, pass variable values to the tree with [`set_varData`](@ref set_varData), run simulations with [`calculateSolution`](@ref calculateSolution), and retrieve results with [`get_varData`](@ref get_varData).

## Implementation Guidelines

To implement a HiPS simulation:

1. Compile the code, ensuring all necessary libraries and dependencies are correctly linked.
2. Run the compiled binary to start the simulation with the specified parameters.
3. Monitor the simulation and analyze the output.

Researchers can adjust parameters based on specific simulation needs.

## Examples

Documented example files are available on the [Examples](pages/Examples.md) page.



<!-- Example files are documented on the [Examples](pages/examples.md).-->




