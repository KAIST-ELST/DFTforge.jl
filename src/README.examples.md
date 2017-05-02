## Exchange_JQ.jl

In "examples/NiO/4cell/Ni6.0_O5.0-s2p2d2f1-4.25_0.0-4.180-k10"
### OpenMX scfout

julia ~/Dropbox/shared/DFT-forge/source/Exchange_JQ.jl -k 5_5_5 -q 5_5_5 -D OpenMX --atom12 1_1,1_2   nio.scfout

### Wannier example
julia Exchange_JQ.jl -k 5_5_5 -q 5_5_5 -D Wannier90 -W openmxWF  nio.HWR  --atom12 1_1,1_2  --om1 1,2 --om2 1,2 --ommode 1 --omname eg_eg -T nio_J_wannier.toml

--