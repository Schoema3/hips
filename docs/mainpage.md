
\mainpage

<!-- #################################################################### -->

# Overview

[HiPS] is an open-source C++ library that serves as both a standalone mixing model and a sub-grid model for Transported Probability Distribution Functions.

# Dependencies and installation

The code is intended to be built and used on Linux-like systems, including MacOS and the Linux subsystem for Windows.

## Required software:
* CMake 3.15+
* C++11
 
## Optional software:
* Cantera and Sundials (for reaction)
* Doxygen (for building documentation)
* Graphviz (for Doxygen)

## Build and installation instructions
1. Create and navigate into a top-level `build` directory
2. Configure CMake: `cmake ..`
3. Build HiPSlib: `make`
4. Install HiPSlib: `make install`
5. Make documentation: `make docs`

## CMake configuration variables

The default CMake configuration should be adequate for users that do not immediately require the examples, tests, or documentation. CMake configuration options can be set by editing the top-level `CMakeLists.txt` file, editing the `CMakeCache.txt` file (generated in the `build` directory after running CMake at least once), or specifying them on the command line during step 2 as follows:
```
cmake -DHIPSLIB_BUILD_EXAMPLES=ON ..
```

The following project-specific CMake configuration variables can be specified by the user; their default values are also indicated.
| CMake variable | Default | Description |
| ----------- | ----- | ------ |
| `CMAKE_INSTALL_PREFIX`   | top-level project directory | Installation location |
| `REACTIONS_ENABLED` | `OFF` | Whether to enable support for chemical reactions within the HiPS library |
| `HIPSLIB_BUILD_EXAMPLES` | `OFF` | Builds HiPS examples |
| `HIPSLIB_BUILD_DOCS`     | `OFF` | Builds HiPS documentation via Doxygen |

The `REACTIONS_ENABLED` flag determines whether the HiPS library supports chemical reactions. When set to   `ON`, the library includes functionality for reactions, requiring additional libraries like Contra or Sundials. When set to `OFF`, only simple mixing is supported, and these additional libraries are not needed. Users can adjust this flag based on whether they require reaction support and have the necessary libraries installed.

![HiPS workflow diagram](Diagram-paper.png)
# Using HiPS

The HiPS library consists of three main object classes that users can interact with: `hips`, `batchReactor_cvode`, and `batchReactor_cantera`. These classes serve as the computational engine responsible for advancing the reaction system through time, solving complex sets of differential equations.

## Example Workflow

To integrate HiPS into C++ code, include the `hips.h` header file. Users interact with HiPS through two main interfaces. The first interface, [hips(int nLevels, double domainLength, double tau0, int nVar, int forceTurb, std::vector<double>& ScHips, bool performReaction, std::shared_ptr<Cantera::Solution> cantSol, int seed)](@ref hips(int, double, double, int, int, std::vector<double>&, bool, std::shared_ptr<Cantera::Solution>, int)), is suited for standalone simulations with comprehensive parameter initialization. The second interface, [hips(int nVar, int forceTurb, std::vector<double>& ScHips, bool performReaction, std::shared_ptr<Cantera::Solution> cantSol, int seed)](@ref hips(int, int, std::vector<double>&, bool, std::shared_ptr<Cantera::Solution>, int)), is optimized for use as a sub-grid model in CFD simulations, where frequent updates to the tree structure are required. This constructor is complemented by the [set_tree(int nLevels, double domainLength, double tau0, std::vector<double>& ScHips)](@ref set_tree(int, double, double, std::vector<double>&)) function, allowing repeated tree initialization, especially based on Reynolds numbers. After initializing the HiPS tree, use the [`set_varData`](@ref set_varData) function to pass variable values (e.g., concentrations, enthalpy) to the tree, the [`calculateSolution`](@ref calculateSolution) function to run the simulations, and the [`get_varData`](@ref get_varData) function to retrieve results. For detailed information on how HiPS simulations work, refer to Section 5.2.

## Implementation Guidelines

To implement a HiPS simulation based on this example in an academic setting, adhere to the following guidelines:

- Comple the code, ensuring all necessary libraries and dependencies are correctly linked.
- Execute the compiled binary to initiate the HiPS simulation with the specified parameters.
- Monitor the simulation progress and analyze the output for relevant information.

Researchers may customize parameters and mechanisms based on specific simulation requirements within the confines of this methodological framework.

## Examples

<!-- Example files are documented on the [Examples](pages/examples.md).-->

HiPS is written in C++ and includes three examples that illustrate its use:

- The example `ex_1.cc` demonstrates a simulation of scalar mixing dynamics, showing how a
basic mixing model behaves. Initially, fluid parcels are split into two groups with values 0 and 1, then mix
to an average value. The base number of levels is 9, with Schmidt numbers of 0.0625 and 16 illustrating the
impact of both high and low Schmidt numbers on mixing. For details on parameter selection, refer to the
[documentation](\ref parameters).


-`ex_2.cc` and `ex_3.cc` illustrate HiPS mixing coupled with chemical
reactions, focusing on combustion processes. The premixed example `ex_2.cc` simulates
turbulent flames using a stoichiometric ethylene/air mixture to initialize a HiPS domain with six levels. In
this setup, 25% of the parcels are pre-combusted, while the remaining 75% are fresh reactants. The Schmidt
number is set to unity for all species, and the domain length scale is 0.01. Experiments are conducted with
two values of \(\tau_0\) (5E-4 and 5E-6) to control mixing rates, impacting ignition and extinction
phenomena. Smaller \(\tau_0\) values correspond to faster mixing rates, enhancing fuel-oxidizer mixing and
promoting ignition, while larger \(\tau_0\) values indicate slower mixing rates, potentially leading to
extinction as combustion reactions struggle to sustain.


The HiPS library supports two integrators for simulating chemical reactions within the framework. The first
integrator, embedded in the `batchReactor_cvode` class, uses the CVODE solver to seamlessly incorporate
chemical kinetics. The second integrator within the `batchReactor_cantera` class offers an alternative
approach, utilizing Cantera's built-in capabilities. Users can choose between these integrators based on
their preferences and requirements, providing flexibility and robustness to efficiently simulate chemical
reactions within the HiPS library.
