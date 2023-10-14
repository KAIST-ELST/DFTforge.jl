- Hongkee Yoon Hongkeeyoon@kaist.ac.kr
- 2023.01 works with Julia 1.0 & 1.9
- https://kaist-elst.github.io/DFTforge.jl/

![Julia1.0](https://img.shields.io/badge/Julia-1.0-blue.svg?longCache=true)  ![Julia1.6](https://img.shields.io/badge/Julia-1.6-blue.svg?longCache=true) ~ ![Julia1.9](https://img.shields.io/badge/Julia-1.9-blue.svg?longCache=true)

The backend package for the Jx MFT code [Jx](https://kaist-elst.github.io/Jx.jl/).

# DFTforge: the DFT postprocessing environment
Simplify obtaining Hamiltonian, Eigenvalue, Eigenvector from DFT results, and Pre-caching them for parallelized calculations.
Pre-cached results are easy to use when calculating in k and q space.


## DFTforge
The wrapper package for calculating Hamiltonian, Eigenvalue, Eigenvector from DFT results.

 * read DFT results from OpenMX, Wannier90
 * Calculate Phi, Enk, Hks,

### Etc. functionalities
 * K-point representation in INT (for unique K, Q).
 * Generate K-points.
 * Generalized argument parser (support TOML style input).



[Homepage & How to use](https://kaist-elst.github.io/DFTforge.jl/)

