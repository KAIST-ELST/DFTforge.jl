---
layout: default
---
![Julia1.0](https://img.shields.io/badge/Julia-1.0-blue.svg?longCache=true)  ![Julia1.1](https://img.shields.io/badge/Julia-1.1-blue.svg?longCache=true) ![Julia1.2](https://img.shields.io/badge/Julia-1.2-blue.svg?longCache=true) 



# What is the DFTforge?
The [DFTforge](https://github.com/KAIST-ELST/DFTforge.jl/) is the backend project for MFT code https://github.com/KAIST-ELST/Jx.jl

## DFT postprocessing environment
Simplify obtaining Hamiltonian, Eigenvalue, Eigenvector from DFT results and Pre-caching them for parallelized calculations.
Pre-caching results are easy to use when calculating in k and q space.
 * It consists of two parts DFTforge & DFTrefinery.

### DFTforge
The wrapper for calculating Hamiltonian, Eigenvalue, Eigenvector from DFT results.

 * read DFT results from OpenMX, Wannier90(?), EcalJ(?)
 * Calculate Phi, Enk, Hks,

#### Etc. functionalities
 * K-point representation in INT (for unique K,Q).
 * Gennerate K-points.
 * Generalised argument parser (support TOML sytle input).



### DFTrefinery
use DFTforge for calcuating various properties
Wrapper module for easy access of DFT-forge, especially easy use of HDF5 cached eigenstate information.

 * Store Eigenvalues & Eigenvectors in HDF5 format (K,Q)
 * Read Stored Eigenvalues & Eigenvectors
 * Simple interface for K space, K,Q space calucations
 ** Generalised parallelized K,Q space calculation function wrapper.

# How to install?
1. Install Julia (1.0 or above)

Visit [https://julialang.org/](https://julialang.org/) for details.

2. Excute the [Julia](https://julialang.org/) and type:
```julia
import Pkg
Pkg.add("DFTforge")
```
All done!

