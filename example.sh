cd examples/NiO/4cell/Ni6.0_O5.0-s2p2d2f1-4.25_0.0-4.180-k10

# obtain J(q)
julia -p 4 ~/.julia/v0.6/DFTforge/src/Spin_Exchange.jl  -T nio_J_openmx.toml

# J(q) -> J(R)
julia ~/.julia/v0.6/DFTforge/src/Jx_postprocess.jl --cellvectors  3_3_3 --atom12 1_1,1_2  --orbital_name all_all --baseatom 1 jq.col.spin_0.0

