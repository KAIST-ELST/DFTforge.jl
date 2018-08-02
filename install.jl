julia_installed_pkg = Pkg.installed()
if !haskey(julia_installed_pkg, "TOML" )
  Pkg.clone("https://github.com/wildart/TOML.jl.git")
end

if !haskey(julia_installed_pkg, "DFTforge" )
  Pkg.clone("https://github.com/bluehope/DFTforge.jl")
end
