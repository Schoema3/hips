
\mainpage

<!-- #################################################################### -->

# Overview

[HiPS](git clone https://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:x-oauth-basic@github.com/BYUignite/hips.git) is an open-source C++ library that is used as a mixing model for Transported Probability Distribution Functions.

# Dependencies and installation

The code is intended to be built and used on Linux-like systems, including MacOS and the Linux subsystem for Windows.

Required software:
* CMake 3.15+
* C++11

Optional software:
* Doxygen (for building documentation)
* graphviz (for Doxygen)

## Build and installation instructions
1. Create and navigate into a top-level `build` directory
2. Configure CMake: `cmake ..`
3. Build HiPSlib: `make`
4. Install HiPSlib: `make install`

## CMake configuration variables
The default CMake configuration should be adequate for users that do not immediately require the examples, tests, or documentation. CMake configuration options can be set by editing the top-level `CMakeLists.txt` file, editing the `CMakeCache.txt` file (generated in the `build` directory after running CMake at least once), or specifying them on the command line during step 2 as follows:
```
cmake -DHIPSLIB_BUILD_EXAMPLES=ON ..
```

The following project-specific CMake configuration variables can be specified by the user; their default values are also indicated.
| CMake variable | Default | Description |
| ----------- | ----- | ------ |
| `CMAKE_INSTALL_PREFIX`   | top-level project directory | Installation location |
| `HIPSLIB_BUILD_EXAMPLES` | `OFF` | Builds HiPS examples |
| `HIPSLIB_BUILD_DOCS`     | `OFF` | Builds HiPS documentation via Doxygen |

# Using HiPS
The HiPS library consists of three main object classes that users can interact with: `hips`, `batchReactor_cvode` and `batchReactor_cantera`. `bachReactor_cvode` and `batchReactor_cantera` are two integrators which serve as the computational engine responsible for advancing the reaction system through time, solving complex sets of differential equations.  

# Example workflow

## 1. HiPS Tree and Constructor Parameters

Define essential parameters for the HiPS simulation, including the number of levels, domain length, time step (tau0), a parameter (C_param), simulation time (tRun), and other pertinent settings.

## 2. Gas Solution Setup (required for reaction)

Initialize the gas solution using Cantera, specifying a mechanism file (e.g., "gri30.yaml"). Extract thermodynamic information and the number of species to facilitate subsequent computations.

## 3. HiPS Tree Creation

Instantiate the HiPS object with specified parameters, encompassing the number of levels, domain length, time step, and additional simulation details.

## 4. Initialize Variables and Parameters

Set up variables and parameters required for the simulation, such as the number of parcels, weights, and variable names.

## 5. Initialize Mixing Fractions

Utilize a designated function to initialize mixing fractions for a predefined number of parcels. The function sets mixing fractions to 0.0 for the first half and 1.0 for the second half.

## 6. Set Initial State in Parcels

Establish the initial state vectors in each parcel using the initialized mixing fractions, weights, and variable names.

## 7. Advance HiPS for Mixing and Reaction

Utilize the HiPS object to advance the simulation for a specified run time (tRun). This involves solving the equations governing particle-laden flows and updating state variables accordingly.

## Implementation Guidelines

To implement a HiPS simulation based on this example in an academic setting, adhere to the following guidelines:

1. Compile the code, ensuring all necessary libraries and dependencies are correctly linked.
2. Execute the compiled binary to initiate the HiPS simulation with the specified parameters.
3. Monitor the simulation progress and analyze the output for relevant information.

Researchers may customize parameters and mechanisms based on specific simulation requirements within the confines of this methodological framework.


 # Examples
HiPS is written in C++ and indludes three examples that illustrate its use. The first example ```ex_1.cc``` shows that how hips is used for simple mixing. ```ex_2.cc``` illustrtares the use of hips in a premixed flame. The third example ```ex_3.cc``` shows the usage of hips in a non-premixed flame. 

