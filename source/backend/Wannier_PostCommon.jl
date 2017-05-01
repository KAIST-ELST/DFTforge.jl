using MAT
using HDF5
using DFTcommon
export cal_eigenstate,cal_Hamiltonian
#export test_SmallHks

#const OLP_eigen_cut = 1.0e-10;
include("Wannier_read.jl");

function Overlap_Band!(wannier_r::Wannierdatatype,spin::Int,
  Hout::Array{Complex_my,2},
  TotalOrbitalNum::Int64,k1,k2,k3)
  #HWR_mat_list::Array{Array{Complex_my,2}}
  FNAN = length(wannier_r.Hks_R[spin]);
  assert(size(wannier_r.R_vector_mat[spin])[1] == FNAN)
  k_point::Array{Float_my,1} = [k1,k2,k3];
  for LB_AN = 1:FNAN
    kRn::Float_my = sum(wannier_r.R_vector_mat[spin][LB_AN,:].*k_point);
    Hout[:,:] += wannier_r.Hks_R[spin][LB_AN].* (cos(2.0*pi*kRn)+sin(2.0*pi*kRn)*im);
  end
end
function cal_Hamiltonian(wannier_r::Wannierdatatype,spin::Int)
  return cal_Hamiltonian((0.0,0.0,0.0),wannier_r,spin)
end
function cal_Hamiltonian(k_point::k_point_Tuple,wannier_r::Wannierdatatype,spin::Int)
  TotalOrbitalNum::Int = sum(wannier_r.Total_NumOrbs[:])
  MPF = zeros(Int,wannier_r.atomnum)

  orbitalStartIdx = 0
  for i = 1:wannier_r.atomnum
      MPF[i] = orbitalStartIdx;
      orbitalStartIdx += wannier_r.Total_NumOrbs[i]
  end
  Hout = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum);
  Overlap_Band!(wannier_r,spin,Hout,TotalOrbitalNum,k_point[1],k_point[2],k_point[3])
  return Hout;
end
function cal_eigenstate(k_point::k_point_Tuple,wannier_r::Wannierdatatype,spin_list::Array{Int})
  TotalOrbitalNum::Int = sum(wannier_r.Total_NumOrbs[:])
  kpoint_common_list = Array{Kpoint_eigenstate}(0);
  for spin in spin_list
    Hout = cal_Hamiltonian(k_point,wannier_r,spin)
    eigvals = zeros(Float_my, TotalOrbitalNum);
    eigstate = copy(Hout)
    DFTcommon.eigfact_hermitian(eigstate,eigvals)
    kpoint_common = Kpoint_eigenstate(eigstate,eigvals,k_point);
    #kpoint_common = Kpoint_eigenstate(Hout,eigvals,k_point);
    push!(kpoint_common_list,kpoint_common)
  end
  return kpoint_common_list;
end
