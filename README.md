# ADSCF
ADSCF is a Self-Consistent Field code package designed for calculating the equilibrium properties of block copolymer materials, encompassing block copolymers in melt/solution and grafted polymer brushes in 1-3D. Developed using Fortran90, the code utilizes a pseudo-spectral method based on Fast Fourier transformation (FFT) for bulk calculations with periodic boundary conditions. Additionally, real-space methods, specifically the Finite Difference method (Crank–Nicolson), are employed for polymer systems with sharp interfaces, such as those found in thin films and end-grafted polymers. Notably, we implement an adaptive scheme for spatial and contour discretization in thin film and brush calculations to minimize numerical errors. For further details, please refer to our publication.

# License 

Copyright (C) 2024  
Le Qiao (<le.qiao@uni-mainz.de>)  
Marios Giannakou (<mgiannak@uni-mainz.de>)  
Friederike Schmid (<friederike.schmid@uni-mainz.de>)  
Johannes Gutenberg-Universität Mainz, Institue of Physics, Schmid group  

ADSCF is an open-source code package provided free of charge: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation. The GNU General Public License is available at <https://www.gnu.org/licenses/>.

## Cite us 
If you use ADSCF or a modified version based on ADSCF to publish scientific papers, please kindly cite our paper and acknowledge the usage of our code.  

Qiao, L.; Giannakou, M.;Schmid, F. An efficient and accurate SCFT algorithm for Block copolymer
films and brushs using adaptive discretizations. Polymers, 2023xxxxxxxx (submitted)



# Files and directories
The structure of the script:\
Using the 3D code of polymer in bulk as an example:
```bash
├── Bulk_3D_FFT
│   ├── General_Routines
│   ├── Geometry_Modules # Numerical method used to solve the diffusion equation, in this emaple FFT is used. 
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
│   ├── Main_Blend # blend of A and B homopolymers
│   │   ├── input_general
│   │   ├── input_interactions
│   │   ├── main_blend.f90
│   │   ├── Makefile
│   │   └── mod_global.f90
│   └── Polymer_Modules # Defining the type of polymer
│       ├── mod_diblock.f90
│       ├── mod_homopolymer.f90
│       ├── mod_solvent.f90
│       └── mod_triblock.f90
```
### Examples using adaptive schemes:
Film_CK_AD: Block copolymer confined in a thin film solved by Crank–Nicolson method using adaptive spatial discretization. The examples are given in 1D and 3D. 

Brush_CK_AD: homopolymer brush solved by Crank–Nicolson method using adaptive spatial and contour discretization.   

# Making and running job
Using the 3D code of polymer in bulk as an example: 
```bash
cd Main_Diblock_Melt
make main_diblock_melt  
./main_diblock_melt.exe  
```
