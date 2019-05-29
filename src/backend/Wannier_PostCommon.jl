###############################################################################
# Hongkee Yoon Hongkeeyoon@kaist.ac.kr
# 2019.05
# https://kaist-elst.github.io/DFTforge.jl/
###############################################################################


using ..DFTcommon
#using MAT
#using HDF5

export cal_eigenstate,cal_Hamiltonian
#export test_SmallHks

#const OLP_eigen_cut = 1.0e-10;
include("Wannier_read.jl");

function Overlap_Band!(wannier_r::Wannierdatatype,spin::Int,
  Hout::Array{Complex_my,2},
  TotalOrbitalNum::Int64,k1,k2,k3)
  #HWR_mat_list::Array{Array{Complex_my,2}}
  FNAN = length(wannier_r.Hks_R[spin]);
  @assert(size(wannier_r.R_vector_mat[spin])[1] == FNAN)
  k_point::Array{Float_my,1} = [k1,k2,k3];
  for LB_AN = 1:FNAN
    kRn::Float_my = sum(wannier_r.R_vector_mat[spin][LB_AN,:].*k_point);
    Hout[:,:] += wannier_r.Hks_R[spin][LB_AN].* (cos(2.0*pi*kRn)+sin(2.0*pi*kRn)*im);
  end
end
#=
function cal_Hamiltonian(wannier_r::Wannierdatatype,spin::Int)
  return cal_Hamiltonian((0.0,0.0,0.0),wannier_r,spin)
end
=#
function cal_Hamiltonian(k_point::k_point_Tuple,   # wannier_r::Wannierdatatype,
  hamiltonian_info::Hamiltonian_info_type, spin::Int)
  wannier_r = hamiltonian_info.scf_r
  spin_type = hamiltonian_info.spin_type::DFTcommon.SPINtype

  TotalOrbitalNum = sum(wannier_r.Total_NumOrbs[:])::Int
  TotalOrbitalNum2 = TotalOrbitalNum;
  if DFTcommon.non_colinear_type == spin_type
    TotalOrbitalNum2 = TotalOrbitalNum*2;
  end

  orbitalStartIdx_list = zeros(Int,wannier_r.atomnum)
  orbitalStartIdx = 0
  for i = 1:wannier_r.atomnum
      orbitalStartIdx_list[i] = orbitalStartIdx;
      orbitalStartIdx += wannier_r.Total_NumOrbs[i]
  end
  Hout = zeros(Complex_my, TotalOrbitalNum2, TotalOrbitalNum2);

  Overlap_Band!(wannier_r, spin, Hout, TotalOrbitalNum2, k_point[1], k_point[2], k_point[3])
  if hamiltonian_info.basisTransform_rule.orbital_rot_on
    if (DFTcommon.non_colinear_type == spin_type)
      throw(assertionError("Non-collinear spin basis rotation not supported yet "));
    end
    Hout = Heff(Hout,orbitalStartIdx_list, hamiltonian_info.basisTransform_rule, 0.0);
    #println( sum(abs(H2-H)) )
  end

  return Hout;
end

function cal_eigenstate(k_point::k_point_Tuple,hamiltonian_info::Hamiltonian_info_type,spin_list::Array{Int})
  wannier_r = hamiltonian_info.scf_r
  spin_type = hamiltonian_info.spin_type;

  TotalOrbitalNum = sum(wannier_r.Total_NumOrbs[:])::Int;
  TotalOrbitalNum2 = TotalOrbitalNum;
  if DFTcommon.non_colinear_type == spin_type
    TotalOrbitalNum2 = TotalOrbitalNum*2;
  end

  kpoint_common_list = Array{Kpoint_eigenstate}(undef,0);

  for spin in spin_list
    Hout = cal_Hamiltonian(k_point, hamiltonian_info, spin)

    eigvals = zeros(Float_my, TotalOrbitalNum2);
    eigstate = copy(Hout)
    DFTcommon.eigfact_hermitian(eigstate,eigvals)
    kpoint_common = Kpoint_eigenstate(eigstate,eigvals,k_point,Hout);
    #kpoint_common = Kpoint_eigenstate(Hout,eigvals,k_point);
    push!(kpoint_common_list,kpoint_common)
  end
  if DFTcommon.non_colinear_type == spin_type
    return kpoint_common_list[1];
  end
  return kpoint_common_list;
end
