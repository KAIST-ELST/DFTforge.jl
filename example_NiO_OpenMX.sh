#!/bin/bash

cd "examples/NiO_G-AFM.OpenMx"

# Unzip dft result
# nio.scf : OpenMX full Hamiltonian info
# nio.HWR : Wannier Hamiltonian from OpenMX
tar xvf nio_dft_result.tar.xz

cd "../../"
# obtain J(q)
julia -p 4 src/Jx_col_spin_exchange.jl  -T examples/NiO_G-AFM.OpenMx/nio_J_openmx.toml
# J(q) -> J(R)
julia  src/Jx_postprocess.jl --cellvectors  3_3_3 --baseatom1 1 --atom2 1,2 --orbital_name all_all  examples/NiO_G-AFM.OpenMx/jx.col.spin_0.0
