\page examples Examples

HiPS is written in C++ and includes three examples that illustrate its use:

- **ex_1.cc:** Demonstrates simple mixing using HiPS, where fluid parcels are initialized with values of 0 and 1, and then mixed to reach an average value.

- **ex_2.cc** and **ex_3.cc:** Illustrate HiPS mixing coupled with chemical reactions, focusing on combustion processes. In these examples, turbulent flames are simulated with unique initialization procedures, showcasing control over ignition and extinction phenomena through parameters like mixing rate ($\tau_0$). Visual insights into the impact of mixing rate on ignition and extinction phenomena are provided.

The HiPS library supports two integrators for simulating chemical reactions within the framework. The first integrator, embedded in the `batchReactor_cvode` class, uses the CVODE solver to seamlessly incorporate chemical kinetics. The second integrator within the `batchReactor_cantera` class offers an alternative approach, utilizing Cantera's built-in capabilities. Users can choose between these integrators based on their preferences and requirements, providing flexibility and robustness to efficiently simulate chemical reactions within the HiPS library.

