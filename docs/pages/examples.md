\page examples Examples

HiPS is written in C++ and includes three examples that illustrate its use:

- The example [`ex_1.cc`](\ref ex_1.cc) demonstrates a simulation of scalar mixing dynamics, showing how a basic mixing model behaves. Initially, fluid parcels are split into two groups with values 0 and 1, then mix to an average value. The base number of levels is 9, with Schmidt numbers of 0.0625 and 16 illustrating the impact of both high and low Schmidt numbers on mixing. For details on parameter selection, refer to the [documentation](\ref parameters).


-[`ex_2.cc`](\ref ex_2.cc)  and [`ex_3.cc`](\ref ex_3.cc) Illustrate HiPS mixing coupled with chemical reactions, focusing on combustion processes. The premixed example [`ex_2.cc`](\ref ex_2.cc) simulates turbulent flames using a stoichiometric ethylene/air mixture to initialize a HiPS domain with six levels. In this setup, 25% of the parcels are pre-combusted, while the remaining 75% are fresh reactants. The Schmidt number is set to unity for all species, and the domain length scale is 0.01. Experiments are conducted with two values of \(\tau_0\) (5E-4 and 5E-6) to control mixing rates, impacting ignition and extinction phenomena. Smaller \(\tau_0\) values correspond to faster mixing rates, enhancing fuel-oxidizer mixing and promoting ignition, while larger \(\tau_0\) values indicate slower mixing rates, potentially leading to extinction as combustion reactions struggle to sustain.
 

The HiPS library supports two integrators for simulating chemical reactions within the framework. The first integrator, embedded in the `batchReactor_cvode` class, uses the CVODE solver to seamlessly incorporate chemical kinetics. The second integrator within the `batchReactor_cantera` class offers an alternative approach, utilizing Cantera's built-in capabilities. Users can choose between these integrators based on their preferences and requirements, providing flexibility and robustness to efficiently simulate chemical reactions within the HiPS library.

