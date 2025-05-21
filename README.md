# Adaptive Mesh Self-Consistent Field (ADSCF) Theory Code Package 
ADSCF is a Fortran90-based code package implementing Self-Consistent Field (SCF) theory on adaptive spatial meshes for the efficient and accurate modeling of inhomogeneous polymer systems.

This package is designed for simulations of block copolymer materials, including:
+ Block copolymers in melt or solution
+ End-grafted polymer brushes
+ Thin films and confined systems
across 1D, 2D, and 3D geometries.

### Key Features:
+ Pseudo-spectral method (based on Fast Fourier Transforms) for bulk systems with periodic boundary conditions
+ Finite difference methods (Crank–Nicolson scheme) for real-space simulations with non-periodic boundary conditions
+ Adaptive mesh and contour discretization for thin films and polymer brushes to reduce numerical errors and enhance resolution where needed

This code provides a flexible and efficient framework for simulating complex morphologies in polymeric systems.

For implementation details, numerical schemes, and example results, please refer to our associated publication:

Qiao, L.; Giannakou, M.; Schmid, F. An efficient and accurate SCFT algorithm for Block copolymer
films and brushs using adaptive discretizations. Polymers 2024, 16(9), 1228; https://doi.org/10.3390/polym16091228

Qiao, L.; Vega, D. A.; Schmid, F. Stability and Elasticity of Ultrathin Sphere-Patterned Block Copolymer Films. Macromolecules, 57 (9), 4629-4634 (2024)

If you use ADSCF or a modified version based on ADSCF to publish scientific papers, please kindly cite our paper and acknowledge the usage of our code.  

# License 

Copyright (C) 2024  
Le Qiao (<le.qiao@uni-mainz.de>)  
Marios Giannakou (<mgiannak@uni-mainz.de>)  
Friederike Schmid (<friederike.schmid@uni-mainz.de>)  
Johannes Gutenberg-Universität Mainz, Institute of Physics, Schmid group  

ADSCF is an open-source code package provided free of charge: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation. The GNU General Public License is available at <https://www.gnu.org/licenses/>.



# Files and directories
### The structure of the code package:
Using the 3D code of polymer in bulk as an example:
```bash
├── Bulk_3D_FFT
│   ├── General_Routines
│   ├── Geometry_Modules # Numerical methods used to solve the diffusion equation, in this example FFT is used. 
│   │   └── Bulk
│   │       ├── fftw3.f
│   │       └── mod_fourier.f90
│   ├── Interaction_Modules # Defining the interactions in melt or in implicit/explicit solvent.
│   │   ├── mod_FloryHuggins.f90
│   │   ├── mod_FloryHuggins_Solution.f90
│   │   └── mod_onecomponent_virialexpansion.f90
│   ├── Main_Diblock_Melt # Directory containing the main.f90 script, input parmeters and Makefiles
│   │   ├── input_general
│   │   ├── input_interactions
│   │   ├── main_diblock_melt.f90
│   │   ├── Makefile
│   │   └── mod_global.f90
│   ├── Main_Blend # blend of A and B homopolymers, SCF and DDFT simulations (chain dynamics or external potential dynamics)
│   │   ├── input_general
│   │   ├── input_interactions
│   │   ├── main_blend.f90
│   │   ├── Makefile
│   │   └── mod_global.f90
│   ├── Main_Multiblock # blend of multiblock copolymers, SCF and DDFT simulations (using mobility matrix function as input)
│   │   ├── input_general
│   │   ├── input_interactions
│   │   ├── input_mobilities
│   │   ├── main_multiblock.f90
│   │   ├── Makefile
│   │   └── mod_global.f90
│   └── Polymer_Modules # Defining the type of polymer
│       ├── mod_diblock.f90
│       ├── mod_homopolymer.f90
│       ├── mod_solvent.f90
│       └── mod_triblock.f90
│       └── mod_multiblock.f90
```
### Examples using adaptive schemes:
+ Film_CK_AD: Block copolymer confined in a thin film solved by Crank–Nicolson method using adaptive spatial discretization. The examples are given in 1D and 3D. 

+ Brush_CK_AD: homopolymer brush solved by Crank–Nicolson method using adaptive spatial and contour discretization.   

# Making and running job
Using the 3D code of polymer in bulk as an example: 
```bash
cd Main_Diblock_Melt
make main_diblock_melt  
./main_diblock_melt.exe  
```
