


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

This section outlines how to run simulations using the HiPS library, assuming prior familiarity with the model as described in the Overview and associated publication.

HiPS can be used in two primary modes:
- As a standalone model for scalar mixing and reacting flows,
- As a sub-grid mixing model embedded in a CFD or PDF solver with tree updates at each timestep.

## Basic Workflow

To use HiPS in a simulation:

1. **Initialize** a `hips` object with the required physical parameters.
2. If needed, **reset the tree** using `set_tree()` for each time step (typical in CFD coupling).
3. **Assign scalar values** using `set_varData()`.
4. **Advance the simulation** using `calculateSolution(tRun, shouldWriteData)`.
5. **Extract updated scalar values** with `get_varData()`.

## Constructors

Two constructor options are available:

### Static Tree Mode

```cpp
hips(
  int nLevels,
  double domainLength,
  double tau0,
  int nVar,
  int forceTurb,
  std::vector<double>& ScHips,
  bool performReaction,
  std::shared_ptr<Cantera::Solution> cantSol,
  int seed
);
```

Use this for standalone simulations when the tree structure is fixed.

### CFD-Compatible Mode

```cpp
hips(
  int nVar,
  int forceTurb,
  std::vector<double>& ScHips,
  bool performReaction,
  std::shared_ptr<Cantera::Solution> cantSol,
  int seed
);
```

Pair this with a per-time-step call to:

```cpp
set_tree(int nLevels, double domainLength, double tau0, std::vector<double>& ScHips);
```

## Main Functions

| Function | Purpose |
|----------|---------|
| `set_varData(...)` | Provide scalar values (e.g., species, temperature) to each parcel. |
| `calculateSolution(...)` | Perform mixing and optional reaction from current time to `tRun`. |
| `get_varData()` | Retrieve updated parcel values. |

## Output and Data Access

HiPS can write output internally during the simulation, controlled by:

- `setOutputIntervalEddy(int interval)`: every N eddy events.
- `setOutputIntervalTime(double interval)`: every Δt seconds.

By default, output is written every 1000 eddy events. Output is only written if `shouldWriteData = true` is passed to `calculateSolution()`.

`get_varData()` returns the scalar values of all parcels after mixing and reaction. This is typically called at the end of each time step to extract results for post-processing or for coupling with an external solver.

## Examples

Documented example files are available on the [Examples](pages/Examples.md) page.



<!-- Example files are documented on the [Examples](pages/examples.md).-->




