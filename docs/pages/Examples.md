

# HiPS Example Usage

HiPS is written in C++ and includes three examples that illustrate its use:

## Example 1: Scalar Mixing Dynamics (ex_1.cc)
The example `ex_1.cc` demonstrates a simulation of scalar mixing dynamics, showing how a basic mixing model behaves. Initially, fluid parcels are split into two groups with values 0 and 1, then mix to an average value. The base number of levels is 9, with Schmidt numbers of 0.0625 and 16 illustrating the impact of both high and low Schmidt numbers on mixing. For details on parameter selection, refer to the [documentation](\ref parameters).


## Example 2: Premixed Combustion (ex_2.cc)
The example `ex_2.cc` illustrates HiPS mixing coupled with chemical reactions, focusing on combustion processes. The premixed example simulates turbulent flames using a stoichiometric ethylene/air mixture to initialize a HiPS domain with six levels. In this setup, 25% of the parcels are pre-combusted, while the remaining 75% are fresh reactants. The Schmidt number is set to unity for all species, and the domain length scale is 0.01.


## Example 3: Subgrid mixing (ex_3.cc)
The example `ex_3.cc` demonstrates use of HiPS as a subgrid mixing model for CFD, providing an outline of the code and calls to the user interface.  Simple mixing similar to Example 1 is done.The example assumes a simple configuration with a few CFD grid cells with each cell containing a different number of flow particles. The user code loops over the grid cells, and in each pass through the loop: 
- resets the HiPS tree, 
- assigns variables from the CFD particles to the HiPS parcels, 
- performs the HiPS mixing, 
- transfers the HiPS parcel data back to the CFD particle data.


### Chemical Reaction Integrators
The HiPS library supports two integrators for simulating chemical reactions within the framework. The first integrator, embedded in the `batchReactor_cvode` class, uses the CVODE solver to seamlessly incorporate chemical kinetics. The second integrator within the `batchReactor_cantera` class offers an alternative approach, utilizing Cantera's built-in capabilities. Users can choose between these integrators based on their preferences and requirements, providing flexibility and robustness to efficiently simulate chemical reactions within the HiPS library.

For detailed compilation and execution instructions, please refer to the comments within each example file.

