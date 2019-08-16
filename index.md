---
layout: default
---
![Julia1.0](https://img.shields.io/badge/Julia-1.0-blue.svg?longCache=true)  ![Julia1.0](https://img.shields.io/badge/Julia-0.7-blue.svg?longCache=true) 

##### Please stay tuned. To be officially  released soon via GitHub will be ready soon! (updated: 2018-09-06 )
We are waiting few more Julia packages to be stabilized.

# MFT (Magentic force theory)

![Jx](docs/Logo_text.svg){: height="100px" }



**JX** Calculate spin exchange coupling parameters *J* from [OpenMX](http://www.openmx-square.org/) DFT/Wannier, [Wannier90](http://www.wannier.org/) hamiltonians via linear-response theory.

We support the following features:
- J in momentum space.
- Orbital resolved *J*
- Local axis redefinition for orbital resolved *J*
- Full [OpenMX](http://www.openmx-square.org/) DFT Hamiltonian and Wannier Hamiltonian ([OpenMX](http://www.openmx-square.org/)/[Wannier90](http://www.wannier.org/))

![Logo](docs/Logo.svg)
---
# Quick-Start

### Install Julia ([https://julialang.org/](https://julialang.org/))

Julia was designed from the beginning for high performance. Julia programs compile to efficient native code for multiple platforms via LLVM.
[See Julia benchmark compared to C, Fortran, Pythons, Matlab & more...](https://julialang.org/benchmarks/).

**Currently we support latest Julia 1.0/0.7.**

 * For Linux system:

[Using Julia auto-installer for Linux](https://github.com/abelsiqueira/jill)

 ```bash
echo 'export PATH=~/opt/bin:$PATH' >>~/.profile
echo 'export PATH=~/opt/bin:$PATH' >>~/.bashrc
JULIA_INSTALL=~/opt/bin bash -ci "$(curl –fsSL https://raw.githubusercontent.com/abelsiqueira/jill/master/jill.sh)"
 ```

 * For OSX or Windows see: [Julia Download](https://julialang.org/downloads/)

### Install DFTforge
```bash
git clone github.com/KAIST-ELST/DFTforge.jl
julia install.jl
```

### Run example

#### G-type AFM (anti ferromagnetic) NiO example
```bash
./example_NiO_OpenMX.sh
```



---
# References

[PhysRevB.97.125132](https://doi.org/10.1103/PhysRevB.97.125132)
# Developer
[Hongkee Yoon](https://github.com/bluehope)
# Funding source

![NRF](http://www.nrf.re.kr/resources/img/icon/icon_eng_logo.png) 

---


