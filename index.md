---
layout: default
---
![Julia1.0](https://img.shields.io/badge/Julia-1.0-blue.svg?longCache=true)  ![Julia1.1](https://img.shields.io/badge/Julia-1.1-blue.svg?longCache=true) ![Julia1.2](https://img.shields.io/badge/Julia-1.2-blue.svg?longCache=true) 



# What is the DFTforge?
The [DFTforge](https://github.com/KAIST-ELST/DFTforge.jl/) is the backend project for MFT code https://github.com/KAIST-ELST/Jx.jl
All codes are based on [Julia](https://julialang.org/).

## DFT postprocessing environment
DFTforge simplifies obtaining Hamiltonian, Eigenvalue, Eigenvector from DFT results.
It also supports caching the Eigenvalue, Eigenvector results for well-parallelized calculations.

### DFTforge
The wrapper for calculating Hamiltonian, Eigenvalue, Eigenvector from DFT results.

 * read DFT results from [OpenMX](http://www.openmx-square.org/), Wannier90(http://www.wannier.org/), EcalJ
 * Calculate Phi, Enk, Hks,
 * Generalised arguments parser (support TOML sytle input).


### DFTrefinery
use DFTforge for calcuating various properties
Wrapper module for easy access of DFTforge, especially easy use of HDF5 cached eigenstate information.

 * Store Eigenvalues & Eigenvectors in HDF5 format (K,Q)
 * Read Stored Eigenvalues & Eigenvectors
 * Simple interface for K space, K,Q space calucations
 ** Generalised parallelized K,Q space calculation function wrapper.

## Why Julia?

Julia was designed from the beginning for high performance. Julia programs compile to efficient native code for multiple platforms via LLVM.
[See Julia benchmark compared to C, Fortran, Pythons, Matlab & more...](https://julialang.org/benchmarks/).


# How to install?

## Install Julia (1.0 or above)
Visit [https://julialang.org/](https://julialang.org/) for details.

## Excute the [Julia](https://julialang.org/) and type:
```julia
import Pkg
Pkg.add("DFTforge")
```
All done!

