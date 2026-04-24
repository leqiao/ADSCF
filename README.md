# ADSCF — Adaptive Mesh Self-Consistent Field Theory

ADSCF is a Fortran 90 code package implementing Self-Consistent Field (SCF) theory on adaptive spatial meshes for efficient and accurate modeling of inhomogeneous polymer systems.

## Overview

This package supports simulations of:
- Block copolymers in melt or solution
- End-grafted polymer brushes
- Thin films and confined polymer systems

across **1D, 2D, and 3D** geometries.

## Key Features

- **Pseudo-spectral method** (FFT-based) for bulk systems with periodic boundary conditions
- **Finite difference method** (Crank–Nicolson scheme) for real-space simulations with non-periodic boundary conditions
- **Adaptive mesh and contour discretization** for thin films and polymer brushes to reduce numerical errors and enhance spatial resolution where needed
- Supports **SCF and DDFT** (Dynamic Density Functional Theory) simulations

## Dependencies

Before building, ensure the following are installed:

| Dependency | Version | Notes |
|---|---|---|
| Fortran compiler | gfortran ≥ 9.0 or ifort ≥ 19 | `gfortran` recommended |
| [FFTW3](http://www.fftw.org/) | ≥ 3.3 | Required for FFT-based bulk codes |
| GNU Make | ≥ 4.0 | For building with provided Makefiles |
| MPI (optional) | OpenMPI or MPICH | Required only if parallel builds are needed |

Install FFTW3 on Linux:
```bash
sudo apt-get install libfftw3-dev   # Debian/Ubuntu
sudo dnf install fftw-devel         # Fedora/RHEL
```

On macOS (Homebrew):
```bash
brew install fftw
```

## Installation & Usage

Clone the repository:
```bash
git clone https://github.com/leqiao/ADSCF.git
cd ADSCF
```

Build and run using the 3D bulk diblock copolymer code as an example:
```bash
cd Bulk_3D_FFT/Main_Diblock_Melt
make main_diblock_melt
./main_diblock_melt.exe
```

Repeat the same `make <main_name>` + `./<main_name>.exe` pattern for any other `Main_*` directory.


## Repository Structure
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
├── Film_CK_AD/                        # Block copolymer thin films (Crank–Nicolson, adaptive mesh)
│   ├── 1D/                            # 1D film geometry
│   └── 3D/                            # 3D film geometry
└── Brush_CK_AD/                       # Polymer brush (Crank–Nicolson, adaptive spatial + contour)
```


## Examples
### Adaptive-Mesh Examples

| Directory | Description |
|---|---|
| `Film_CK_AD/1D` | Diblock copolymer film in 1D, Crank–Nicolson + adaptive spatial mesh |
| `Film_CK_AD/3D` | Same, extended to 3D geometry |
| `Brush_CK_AD` | Homopolymer brush with adaptive spatial and contour discretization |
# Making and running job
Using the 3D code of polymer in bulk as an example: 
```bash
cd Main_Diblock_Melt
make main_diblock_melt  
./main_diblock_melt.exe  
```
## Citation
If you use ADSCF or a derivative in published work, please cite:

> Qiao, L.; Giannakou, M.; Schmid, F. An efficient and accurate SCFT algorithm for block copolymer films and brushes using adaptive discretizations. *Polymers* **2024**, *16*(9), 1228. https://doi.org/10.3390/polym16091228

> Qiao, L.; Vega, D. A.; Schmid, F. Stability and elasticity of ultrathin sphere-patterned block copolymer films. *Macromolecules* **2024**, *57*(9), 4629–4634. https://doi.org/10.1021/acs.macromol.4c00345

## License

Copyright © 2024  
Le Qiao (<le.qiao@uni-mainz.de>)  
Marios Giannakou (<mgiannak@uni-mainz.de>)  
Friederike Schmid (<friederike.schmid@uni-mainz.de>)  
Johannes Gutenberg-Universität Mainz, Institute of Physics, Schmid group

ADSCF is free software: you can redistribute it and/or modify it under the terms of the **GNU General Public License** as published by the Free Software Foundation (version 3 or later). See [`LICENSE`](LICENSE) for the full text, or visit <https://www.gnu.org/licenses/>.

