

# HiPS Example Usage

HiPS is written in C++ and includes three examples that illustrate its use:

## Example 1: Scalar Mixing Dynamics (ex_1.cc)
The example `ex_1.cc` demonstrates a simulation of scalar mixing dynamics, showing how a basic mixing model behaves. Initially, fluid parcels are split into two groups with values 0 and 1, then mix to an average value. The base number of levels is 9, with Schmidt numbers of 0.0625 and 16 illustrating the impact of both high and low Schmidt numbers on mixing. For details on parameter selection, refer to the [documentation](\ref parameters).


## Example 2: Premixed Combustion (ex_2.cc)
The example `ex_2.cc` illustrates HiPS mixing coupled with chemical reactions, focusing on combustion processes. The premixed example simulates turbulent flames using a stoichiometric ethylene/air mixture to initialize a HiPS domain with six levels. In this setup, 25% of the parcels are pre-combusted, while the remaining 75% are fresh reactants. The Schmidt number is set to unity for all species, and the domain length scale is 0.01.


## Example 3: Non-Premixed Combustion (ex_3.cc)
The example `ex_3.cc` also demonstrates HiPS mixing coupled with chemical reactions, but focuses on non-premixed combustion processes. In this setup, the initial condition consists of separate parcels for fuel and oxidizer. The mixing dynamics and ignition/extinction behavior are similar to the premixed case but differ in initialization. The flexibility of HiPS allows users to switch between premixed and non-premixed scenarios by altering the initial configuration.


### Chemical Reaction Integrators
The HiPS library supports two integrators for simulating chemical reactions within the framework. The first integrator, embedded in the `batchReactor_cvode` class, uses the CVODE solver to seamlessly incorporate chemical kinetics. The second integrator within the `batchReactor_cantera` class offers an alternative approach, utilizing Cantera's built-in capabilities. Users can choose between these integrators based on their preferences and requirements, providing flexibility and robustness to efficiently simulate chemical reactions within the HiPS library.

For detailed compilation and execution instructions, please refer to the comments within each example file.

