###############################################################################
# Hongkee Yoon Hongkeeyoon@kaist.ac.kr
# 2019.05
# https://kaist-elst.github.io/DFTforge.jl/
###############################################################################


__precompile__(true)
#


module DFTforge

export DFTtype, SPINtype
export Hamiltonian_info_type
export k_point_Tuple,k_point_int_Tuple
export cal_colinear_eigenstate,cal_colinear_Hamiltonian
export cal_colinear_eigenstate_as_nc
export cal_nonco_linear_Eigenstate,cal_noncolinear_Hamiltonian
#DFTforge_VERSION = VersionNumber("0.6.4-dev+20180827");
include("UpdateCheck.jl")


#import ..DFTcommon
#using DFTforge.DFTcommon
include("../ext/TOML/src/TOML.jl")

using Distributed
using LinearAlgebra
## to export
using ArgParse
using ProgressMeter
using Statistics
using CSV
using FileIO
using DataFrames
using Glob
#using Plots
export ArgParse,ProgressMeter,Distributed,Statistics,CSV,FileIO,DataFrames,Glob#,Plots
##

#@enum DFTtype OpenMX = 1 Wannier90 = 2
#@enum SPINtype para_type = 1 colinear_type = 2 non_colinear_type = 4
include("DFTcommon.jl")
using .DFTcommon


module OpenMXdata
include("backend/OpenMX_PostCommon.jl")
end
module EcalJdata
include("backend/EcalJ_PostCommon.jl")
end
module Wannierdata
include("backend/Wannier_PostCommon.jl")
end

module Plugins
include("plugins/OpenMX_scfout_update.jl")
end





#type OpenMX_data
#    scf_name::AbstractString
#    scf_r::OpenMXdata.Openmxscf
#    dfttype::DFTforge.DFTtype
#end
export read_dftresult
function read_dftresult(scf_fname::AbstractString, result_file_dict::Dict{AbstractString,AbstractString},
  dfttype::DFTcommon.DFTtype,spin_type::DFTcommon.SPINtype,
  basisTransform_rule::basisTransform_rule_type=basisTransform_rule_type())
    if (DFTcommon.OpenMX == dfttype)
      scf_r = OpenMXdata.read_scf(scf_fname);
      hamiltonian_info = Hamiltonian_info_type(scf_r,dfttype,spin_type,basisTransform_rule)
      return hamiltonian_info;
    elseif (DFTcommon.EcalJ == dfttype)
      scf_r = EcalJdata.read_EcalJ_scf(result_file_dict, spin_type);
      hamiltonian_info = Hamiltonian_info_type(scf_r,dfttype,spin_type,basisTransform_rule)
      return hamiltonian_info;
    end

end

function read_dftresult(wannier_fname::AbstractString,result_file_dict::Dict{AbstractString,AbstractString},
  dfttype::DFTcommon.DFTtype,
  Wannier90_type::DFTcommon.Wannier90type,spin_type::DFTcommon.SPINtype,
    atoms_orbitals_list::Vector{Array{Int64}},
  atomnum::Int,atompos::Array{Float64,2},basisTransform_rule::basisTransform_rule_type=basisTransform_rule_type())
  if (DFTcommon.Wannier90 == dfttype)

    wannier_r = Wannierdata.read_wannier(wannier_fname, result_file_dict, Wannier90_type, spin_type,
      atoms_orbitals_list,atomnum,atompos)
    hamiltonian_info = Hamiltonian_info_type(wannier_r,dfttype,spin_type,basisTransform_rule)
    return hamiltonian_info;
  end
end


function cal_colinear_eigenstate(k_point::k_point_Tuple,
    hamiltonian_info::Hamiltonian_info_type,spin_list=1)
    Eigenstate::Array{Kpoint_eigenstate} = Array{Kpoint_eigenstate}(undef,1);
    dfttype = hamiltonian_info.dfttype;
    if (DFTcommon.OpenMX == dfttype)
      #Eigenstate::DFTforge.Kpoint_eigenstate =
      Eigenstate =
       OpenMXdata.cal_colinear_eigenstate(k_point,hamiltonian_info,spin_list);
    elseif (DFTcommon.EcalJ == dfttype)
      Eigenstate =
        EcalJdata.cal_colinear_eigenstate(k_point,hamiltonian_info,spin_list);
    elseif (DFTcommon.Wannier90 == dfttype)
      Eigenstate = Wannierdata.cal_eigenstate(k_point,hamiltonian_info,spin_list)
    end

    return Eigenstate;
end

function cal_colinear_eigenstate_as_nc(k_point::k_point_Tuple,
    hamiltonian_info::Hamiltonian_info_type)
    #Eigenstate::Array{Kpoint_eigenstate} = Array{Kpoint_eigenstate}();
    dfttype = hamiltonian_info.dfttype;
    if (DFTcommon.OpenMX == dfttype)
      #Eigenstate::DFTforge.Kpoint_eigenstate =
      Eigenstate =
       OpenMXdata.cal_colinear_eigenstate_as_nc(k_point,hamiltonian_info);
       return Eigenstate;
    elseif (DFTcommon.Wannier90 == dfttype)
      Eigenstate = Wannierdata.cal_eigenstate_as_nc(k_point,hamiltonian_info)
      return Eigenstate;
    end

    return Eigenstate;
end


function cal_nonco_linear_Eigenstate(k_point::k_point_Tuple,
    dfttype::DFTcommon.DFTtype,hamiltonian_info::Hamiltonian_info_type)
    if (DFTcommon.OpenMX == dfttype)
      EigenState = OpenMXdata.cal_noncolinear_eigenstate(k_point, hamiltonian_info)
      return EigenState
    elseif (DFTcommon.Wannier90 == dfttype)
      EigenState = Wannierdata.cal_eigenstate(k_point,hamiltonian_info,[1])
      return EigenState
    end
end

function cal_colinear_Hamiltonian(k_point::k_point_Tuple,
  dfttype::DFTcommon.DFTtype,hamiltonian_info::Hamiltonian_info_type,spin=1)

  H::Hamiltonian_type = zeros(Complex_my,2,2);
  if (DFTcommon.OpenMX==dfttype)
    H =  OpenMXdata.colinear_Hamiltonian(k_point,hamiltonian_info,spin)
  elseif (DFTcommon.Wannier90 == dfttype)
    H = Wannierdata.cal_Hamiltonian(k_point,hamiltonian_info,spin)
  end
  return H
end

function cal_noncolinear_Hamiltonian(k_point::k_point_Tuple,
  dfttype::DFTcommon.DFTtype,
  Hmode::DFTcommon.nc_Hamiltonian_selection,hamiltonian_info::Hamiltonian_info_type)
  if (DFTcommon.OpenMX==dfttype)
    H = OpenMXdata.noncolinear_Hamiltonian(k_point,hamiltonian_info,Hmode);
    return H
  elseif (DFTcommon.Wannier90 == dfttype)
    H = Wannierdata.cal_Hamiltonian(k_point, hamiltonian_info, 1)
    return H
  end
end




################################################################################
module DFTrefinery
include("DFTrefinery.jl")
end
end
