# CR3BP MATLAB Library

A MATLAB library I developed for use in university projects and Master thesis at Politecnico di Milano.

## Features

This is a non-exhaustive list of features in the library
- CR3BP orbit propagation
- CR3BP STM propagation
- CR3BP characteristic value computations
- Lagrange point computation
- Jacobian constant estimation
- A set of single shooting Newton differential correctors:
   - variuos symmetric peroidic orbit correctors for L1/L2/L3 LPO and DRO families 
   - a few asymmetric periodic orbit correctors for L4/L5 LPO families and 3D orbits
- Two-level multi-shooting Newton differntial correctors:
   - Method developed by K.C. Howell for NRHOs
   - Alternate method developed by D. Grebow
- Pseudo-arc length continuation methods
- Analytical Monodromy matrix eigenvalue computation and stability
- Analytical first order approximation for LPO
- Analytical thrid order approximation for Halo orbits
- Conversions between non-dimensional rotational to dimensional inertial frames
- Conversions between non-dimensional rotational to Ephemeris frames
- Other frame, vector and keplarian parameter conversions

**Note** I tried to specify any and all references I used in each module. But if you find any missing or if any mistakes in the implementation exists, do contact me on my e-mail: kevincharls.work@mail.com
