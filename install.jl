import Pkg
julia_installed_pkg = Pkg.installed()
#if !haskey(julia_installed_pkg, "TOML" )
#  Pkg.clone("https://github.com/wildart/TOML.jl.git")
#end

#if !haskey(julia_installed_pkg, "DFTforge" )
#  Pkg.clone("https://github.com/KAIST-ELST/DFTforge.jl")
#else
#  Pkg.update("DFTforge")
#end

Pkg.add("Plots")
Pkg.add("JLD2")
Pkg.add("Glob")
Pkg.add("HDF5")
Pkg.add("ColorTypes")
Pkg.add("ProgressMeter")

Pkg.build("HDF5")
Pkg.develop(Pkg.PackageSpec(path=pwd()))
import DFTforge
