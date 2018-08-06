
# QuickStart Guide

## Install
### Install Julia

For Linux system:
```
JULIA_INSTALL=~/opt/bin bash -ci "$(curl –fsSL https://raw.githubusercontent.com/abelsiqueira/jill/master/jill.sh)"
```
For OSX or Windows see: https://julialang.org/downloads/

### Install DFTforge in Julia

```sh
git clone https://github.com/ElectronicStructureTheory-KAIST/DFTforge.jl/
julia install.jl
```

## Run example

 * G-type AFM NiO example
```
./example_NiO_OpenMX.sh
```
---

# DFT postprocessing environment
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
