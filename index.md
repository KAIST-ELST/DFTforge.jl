---
layout: default
---

# MFT (Magentic force theory)
Calculate spin exchange coupling parameters J from ([OpenMX](http://www.openmx-square.org/))DFT/Wannier hamiltonians via linear-response theory.
We support the following features.
- J in momentum space.
- Orbital resolved J
- Local axis redefinition for orbital resolved J
- Full [OpenMX](http://www.openmx-square.org/) DFT Hamiltonian and Wannier Hamiltonian ([OpenMX](http://www.openmx-square.org/)/Wannier90)

![Logo](Logo.svg)

# Quick-Start

### Install Julia ([https://julialang.org/](https://julialang.org/))

Julia was designed from the beginning for high performance. Julia programs compile to efficient native code for multiple platforms via LLVM.
[See Julia benchmark compared to C, Fortran, Pythons, Matlab & more...](https://julialang.org/benchmarks/).

 * For Linux system: 

[Using Julia auto-installer for Linux](https://github.com/abelsiqueira/jill)
 ```bash
 JULIA_INSTALL=~/opt/bin bash -ci "$(curl –fsSL https://raw.githubusercontent.com/abelsiqueira/jill/master/jill.sh)"
 ```
 
 * For OSX or Windows see: [Julia Download](https://julialang.org/downloads/)

### Install DFTforge
```bash
git clone https://github.com/ElectronicStructureTheory-KAIST/DFTforge.jl/
julia install.jl
```

### Run example

#### G-type AFM (anti ferromagnetic) NiO example
```bash
./example_NiO_OpenMX.sh
```
---

# DFTforge the DFT postprocessing environment
Simplify obtaining Hamiltonian, Eigenvalue, Eigenvector from DFT results and Pre-caching them for parallelized calculations.
Pre-caching results are easy to use when calculating in k and q space.
 * It consists of two parts DFTforge & DFTrefinery.


## DFTforge
The wrapper for calculating Hamiltonian, Eigenvalue, Eigenvector from DFT results.

 * read DFT results from OpenMX, Wannier90(?), EcalJ(?)
 * Calculate Phi, Enk, Hks,

### Etc. functionalities
 * K-point representation in INT (for unique K,Q).
 * Gennerate K-points.
 * Generalised argument parser (support TOML sytle input).



## DFTrefinery
use DFTforge for calcuating various properties
Wrapper module for easy access of DFT-forge, especially easy use of HDF5 cached eigenstate information.

 * Store Eigenvalues & Eigenvectors in HDF5 format (K,Q)
 * Read Stored Eigenvalues & Eigenvectors
 * Simple interface for K space, K,Q space calucations
 ** Generalised parallelized K,Q space calculation function wrapper.


