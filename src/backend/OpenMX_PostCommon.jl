#module OpenMXdata

using MAT
using HDF5
using LinearAlgebra
export update_co_linear_Energy,cal_noncolinear_eigenstate
export test_SmallHks
export colinear_Hamiltonian,noncolinear_Hamiltonian
export cal_colinear_raw_Hamiltonian
const OLP_eigen_cut = 1.0e-10;


using ..DFTcommon

include("OpenMX_read_scf.jl")
#include("read_scf.jl")

#@enum Spin_mode_enum nospin=1 colinear_spin=2 noncolinear_spin=4

#used for all atoms
#global scf_r = Openmxscf
#@enum nc_Hamiltonian_selection nc_allH=0 nc_realH_only=1 nc_imagH_only=2

function Overlap_Band!(OLP::Overlap_type,
    S::Array{Complex{Float_my},2},MP::Array{Int,1},k1,k2,k3,scf_r::Openmxscf)
    #@acc begin
    k_point::Array{Float_my,1} = [k1,k2,k3];
    TotalOrbitalNum::Int = sum(scf_r.Total_NumOrbs[:]);
    #orbitalStartIdx::Int = 0; #각 atom별로 orbital index시작하는 지점

    #@assert(TotalOrbitalNum==orbitalStartIdx);


    #println(TotalOrbitalNum)
    #S_size::Int = orbitalNum + 0;
    S1 = zeros(Float64,TotalOrbitalNum,TotalOrbitalNum);
    S2 = zeros(Float64,TotalOrbitalNum,TotalOrbitalNum);
    @assert(size(S) == size(S1));

    for GA_AN=1:scf_r.atomnum

      atom1_orbitalNum = scf_r.Total_NumOrbs[GA_AN]
      atom1_orbitalStart = MP[GA_AN];
      for LB_AN = 1:scf_r.FNAN[GA_AN]+1 #atom_i is not atom1,2 index

            GB_AN::UInt = scf_r.natn[GA_AN][LB_AN]
            Rn::UInt = 1+scf_r.ncn[GA_AN][LB_AN]
            atom2_orbitalNum::UInt = scf_r.Total_NumOrbs[GB_AN]
            atom2_orbitalStart::UInt = MP[GB_AN];
            kRn::Float_my = sum(scf_r.atv_ijk[Rn,:][:].*k_point);
            si::Float_my = sin(2.0*pi*kRn);
            co::Float_my = cos(2.0*pi*kRn);
            for i = 1:atom1_orbitalNum
                for j =1:atom2_orbitalNum
                    s = OLP[GA_AN][LB_AN][i][j];

                    S1[atom1_orbitalStart+i,atom2_orbitalStart+j] +=  s*co;
                    S2[atom1_orbitalStart+i,atom2_orbitalStart+j] +=  s*si;

                    #S1[atom1_orbitalStart+i,atom2_orbitalStart+j] =  s*co;
                    #S2[atom1_orbitalStart+i,atom2_orbitalStart+j] =  s*si;

                end
            end
        end
    end

    for i = 1:TotalOrbitalNum
      for j = 1:TotalOrbitalNum
          S[i,j] = S1[i,j]+S2[i,j]*im;
      end
    end
    #end
end
function cal_colinear_eigenstate(k_point::k_point_Tuple,hamiltonian_info::Hamiltonian_info_type,spin_list=1)
  scf_r = hamiltonian_info.scf_r;
  ## Overlapmaxitrx S

  TotalOrbitalNum = sum(scf_r.Total_NumOrbs[:]);
  S = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum)
  orbitalStartIdx_list = zeros(Int,scf_r.atomnum)
  orbitalNums = copy(scf_r.Total_NumOrbs);
  #MPF = Array(Int,scf_r.atomnum)
  orbitalStartIdx = 0
  for i = 1:scf_r.atomnum
      orbitalStartIdx_list[i] = orbitalStartIdx;
      orbitalStartIdx += orbitalNums[i]
  end

  Overlap_Band!(scf_r.OLP,S,orbitalStartIdx_list,k_point[1],k_point[2],k_point[3],scf_r);
  # S rotation  & Merge
  if hamiltonian_info.basisTransform_rule.orbital_rot_on
    S = Heff(S,orbitalStartIdx_list,hamiltonian_info.basisTransform_rule,0.0);
  end
  #Sq = sqrtm(S)
  #S2 = inv(Sq);
  S2 = sqrtm_inv(S);

  #H_temp = Array(Array{Complex_my,2},scf_r.SpinP_switch+1);
  kpoint_eigenstate_list =  Array{Kpoint_eigenstate}(undef,0);
  for spin=spin_list #scf_r.SpinP_switch+1
      Coes = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum);
      H = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum);
      Overlap_Band!(scf_r.Hks[spin],H,orbitalStartIdx_list,k_point[1],k_point[2],k_point[3],scf_r);
      # H rotation & Merge
      if hamiltonian_info.basisTransform_rule.orbital_rot_on
        H = Heff(H,orbitalStartIdx_list,hamiltonian_info.basisTransform_rule,0.0);
        #println( sum(abs(H2-H)) )
      end
      Hk_tilta = S2'*(H*S2);
      Eigen_vect = copy(Hk_tilta);
      Eigen_value = zeros(TotalOrbitalNum);

      eigfact_hermitian(Eigen_vect,Eigen_value);
      if (!check_eigmat(Hk_tilta,Eigen_vect,Eigen_value))
          println("C2 : ",k_point,"\tspin: ",spin)
      end
      #
      #Coes = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum);

      #Coes = S3*Ev_H2;
      #kpoint_common::Kpoint_eigenstate = Kpoint_eigenstate(Coes,ko,k_point,H);
      kpoint_common::Kpoint_eigenstate = Kpoint_eigenstate(Eigen_vect, Eigen_value,k_point,Hk_tilta);
      #kpoint_common::Kpoint_eigenstate = Kpoint_eigenstate(Coes,ko,k_point,H);
      push!(kpoint_eigenstate_list,kpoint_common);
  end

  return kpoint_eigenstate_list;
end
function colinear_Hamiltonian_as_nc!(Hout::Array{Complex{Float_my},2},
    H::H_type,#iH::H_type,
    MP,k1::Float64,k2::Float64,k3::Float64,
    Hmode::nc_Hamiltonian_selection,scf_r::Openmxscf)
    # Essentially same as Overlap_Band!
	#println(size(H))
    @assert(2 == size(H)[1]);
    #@assert(3 == size(iH)[1]);

    k_point::Array{Float_my,1} = [k1,k2,k3];
    TotalOrbitalNum::Int = sum(scf_r.Total_NumOrbs[:])
    orbitalStartIdx::Int = 0; #각 atom별로 orbital index시작하는 지점
    MP = Array{Int}(scf_r.atomnum)
    for i = 1:scf_r.atomnum
        MP[i] = orbitalStartIdx;
        orbitalStartIdx += scf_r.Total_NumOrbs[i]
    end
    @assert(TotalOrbitalNum==orbitalStartIdx);

    # non collinear Hamiltonian: the matrix size is 2*TotalOrbitalNum,2*TotalOrbitalNum
    #Hout = zeros(Complex_my,2*TotalOrbitalNum,2*TotalOrbitalNum);
    @assert((2*TotalOrbitalNum,2*TotalOrbitalNum) == size(Hout));

    for GA_AN=1:scf_r.atomnum
        atom1_orbitalNum = scf_r.Total_NumOrbs[GA_AN];
        atom1_orbitalStart = MP[GA_AN];
        for LB_AN = 1:scf_r.FNAN[GA_AN]+1 #atom_i is not atom1,2 index
            GB_AN::UInt = scf_r.natn[GA_AN][LB_AN]
            Rn::UInt = 1+scf_r.ncn[GA_AN][LB_AN]
            atom2_orbitalNum::UInt = scf_r.Total_NumOrbs[GB_AN]
            atom2_orbitalStart::UInt = MP[GB_AN];
            kRn::Float_my = sum(scf_r.atv_ijk[Rn,:][:].*k_point);
            si::Float_my = sin(2.0*pi*kRn);
            co::Float_my = cos(2.0*pi*kRn);
            si_co = co + im*si;
            #
            for i = 1:atom1_orbitalNum
                for j = 1:atom2_orbitalNum
                    RH1 =  H[1][GA_AN][LB_AN][i][j];
                    RH2 =  H[2][GA_AN][LB_AN][i][j];
                    #RH3 =  H[3][GA_AN][LB_AN][i][j];
                    #RH4 =  H[4][GA_AN][LB_AN][i][j];
#=
                    iH1 =  iH[1][GA_AN][LB_AN][i][j];
                    iH2 =  iH[2][GA_AN][LB_AN][i][j];
                    iH3 =  iH[3][GA_AN][LB_AN][i][j];
=#
                    if DFTcommon.nc_realH_only == Hmode
                        iH1 *= 0 ;iH2 *= 0 ;iH3 *= 0 ;
                    elseif DFTcommon.nc_imagH_only == Hmode
                        RH1 *= 0 ;RH2 *= 0 ;RH3 *= 0 ;RH4 *= 0 ;
                    end

                    Hout[atom1_orbitalStart+i,atom2_orbitalStart+j] +=
                    (RH1) * si_co;
                    Hout[atom1_orbitalStart+i+TotalOrbitalNum,atom2_orbitalStart+j+TotalOrbitalNum] +=
                    (RH2) * si_co;
                    #Hout[atom1_orbitalStart+i,atom2_orbitalStart+j+TotalOrbitalNum] +=
                    #(RH3+im*(RH4+iH3)) * si_co;

                end
            end
        end
    end
    # set off-diagnoal part
    #Hout[TotalOrbitalNum+(1:TotalOrbitalNum),1:TotalOrbitalNum] =
    #   conj( Hout'[1:TotalOrbitalNum,TotalOrbitalNum+(1:TotalOrbitalNum)]);
    for i = (1:TotalOrbitalNum)
        for j = TotalOrbitalNum+(1:TotalOrbitalNum)
            Hout[j,i] = conj( Hout[i,j]);
            #Hout[i,j] = conj( Hout[i,j]); # test code
        end
    end
end

function cal_colinear_eigenstate_as_nc(k_point::k_point_Tuple,hamiltonian_info::Hamiltonian_info_type)
    scf_r = hamiltonian_info.scf_r;
    #function nc_update_Energy(k_point_int::k_point_int_Tuple)
    # non collinear Enk and Eigen function \Psi


    #k_point = k_point_int2float(k_point_int);
    ## Common variables
    TotalOrbitalNum = sum(scf_r.Total_NumOrbs[:])
    TotalOrbitalNum3 = TotalOrbitalNum*2;

    ## Overlap matrix S
    S = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum)
    orbitalStartIdx_list = zeros(Int,scf_r.atomnum)
    orbitalStartIdx = 0;
    for i = 1:scf_r.atomnum
        orbitalStartIdx_list[i] = orbitalStartIdx;
        orbitalStartIdx += scf_r.Total_NumOrbs[i]
    end
    Overlap_Band!(scf_r.OLP,S,orbitalStartIdx_list,k_point[1],k_point[2],k_point[3],scf_r);
    #S = S';
    if hamiltonian_info.basisTransform_rule.orbital_rot_on
      S = Heff(S,orbitalStartIdx_list,hamiltonian_info.basisTransform_rule,0.0);
    end

    S2 = sqrtm_inv(S);
    S22 = zeros(Complex_my, 2*TotalOrbitalNum, 2*TotalOrbitalNum)
    up_spin = 1:TotalOrbitalNum;
    dn_spin = TotalOrbitalNum + up_spin;
    S22[up_spin,up_spin] = S2;
    S22[dn_spin,dn_spin] = S2;
    ## non-collinear Eigen funtion (S^(1/2)\Psi) and Eigen values (Enk)
    ##
    #H0 = zeros(Complex_my,2*TotalOrbitalNum,2*TotalOrbitalNum)
    H1 = zeros(Complex_my,2*TotalOrbitalNum,2*TotalOrbitalNum)

    # H_orig is updated
    colinear_Hamiltonian_as_nc!(H1,scf_r.Hks,orbitalStartIdx_list,k_point[1],
    k_point[2],k_point[3],DFTcommon.nc_allH,scf_r);
    # C = M1 Ut H U M1
    if hamiltonian_info.basisTransform_rule.orbital_rot_on
        orbitals_up =   1:TotalOrbitalNum
        orbitals_down =   orbitals_up + TotalOrbitalNum
        Htmp = zeros(Complex_my,2*TotalOrbitalNum,2*TotalOrbitalNum)
        H = H1[orbitals_up,orbitals_up];
        Htmp[orbitals_up,orbitals_up] = Heff(H,orbitalStartIdx_list,hamiltonian_info.basisTransform_rule,0.0);

        H = H1[orbitals_up,orbitals_down];
        Htmp[orbitals_up,orbitals_down] = Heff(H,orbitalStartIdx_list,hamiltonian_info.basisTransform_rule,0.0);

        H = H1[orbitals_down,orbitals_up];
        Htmp[orbitals_down,orbitals_up] = Heff(H,orbitalStartIdx_list,hamiltonian_info.basisTransform_rule,0.0);

        H = H1[orbitals_down,orbitals_down];
        Htmp[orbitals_down,orbitals_down] = Heff(H,orbitalStartIdx_list,hamiltonian_info.basisTransform_rule,0.0);
        H1 = Htmp;
        H = 0;
        #println( sum(abs(H2-H)) )
    end

    # find Eigenvalues
    Hk_tilta = S22' * H1 * S22;
    Eigen_vect = copy(Hk_tilta)
    Eigen_value = zeros(Float_my,2*TotalOrbitalNum);


    eigfact_hermitian(Eigen_vect,Eigen_value);
    if (!check_eigmat(Hk_tilta,Eigen_vect,Eigen_value))
     println("H3 :",k_point)
    end
    kpoint_common = Kpoint_eigenstate(Eigen_vect,Eigen_value,k_point,Hk_tilta);
    return kpoint_common
end


## Reference implimentation wihtout orbital rot & merge
function cal_colinear_eigenstate_deprecated(k_point::k_point_Tuple,scf_r::Openmxscf,spin_list=1)

  ## Overlapmaxitrx S
  TotalOrbitalNum = sum(scf_r.Total_NumOrbs[:]);
  S = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum)
  MPF = zeros(Int,scf_r.atomnum)
  #MPF = Array(Int,scf_r.atomnum)
  orbitalStartIdx = 0
  for i = 1:scf_r.atomnum
      MPF[i] = orbitalStartIdx;
      orbitalStartIdx += scf_r.Total_NumOrbs[i]
  end

  Overlap_Band!(scf_r.OLP,S,MPF,k_point[1],k_point[2],k_point[3],scf_r);

  S2 = copy(S);
  S_eigvals = zeros(Float_my, TotalOrbitalNum);
  eigfact_hermitian(S2,S_eigvals);
  if (!check_eigmat(S,S2,S_eigvals))
    println("S :",k_point)
  end

  M1 = zeros(size(S_eigvals));
  M1[S_eigvals.>OLP_eigen_cut] = 1.0 ./sqrt.(S_eigvals[S_eigvals.>OLP_eigen_cut]);

  ## Coes

  #H_temp = Array(Array{Complex_my,2},scf_r.SpinP_switch+1);
  kpoint_eigenstate_list =  Array{Kpoint_eigenstate}(undef,0);
  for spin=spin_list #scf_r.SpinP_switch+1
      Coes = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum);
      H = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum);
      Overlap_Band!(scf_r.Hks[spin],H,MPF,k_point[1],k_point[2],k_point[3],scf_r);
      # C = M1 S2^t H S2 M1
      # C = M1 U^t H U M1
      #H_temp[spin] = H;
      C::Array{Complex_my,2} = H*S2;
      for jj = 1:TotalOrbitalNum
          C[:,jj]*= M1[jj];
      end
      C2::Array{Complex_my,2} = (S2'*C); #check if conj(S2) is needed
      for ii = 1:TotalOrbitalNum
          C2[ii,:]*= M1[ii];
      end

      C3 = copy(C2);
      ko = zeros(TotalOrbitalNum);
      eigfact_hermitian(C3,ko);
      if (!check_eigmat(C2,C3,ko))
          println("C2 : ",k_point,"\tspin: ",spin)
      end
      #
      Coes = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum);
      S3 = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum);
      for ii=1:TotalOrbitalNum
          S3[:,ii] = S2[:,ii] * M1[ii];
      end
      Coes = S3*C3;
      #kpoint_common::Kpoint_eigenstate = Kpoint_eigenstate(Coes,ko,k_point,H);
      kpoint_common::Kpoint_eigenstate = Kpoint_eigenstate(Coes,ko,k_point,H);
      push!(kpoint_eigenstate_list,kpoint_common);
  end
  return kpoint_eigenstate_list;
end
function colinear_Hamiltonian(k_point::k_point_Tuple,spin::Int,scf_r::Openmxscf)
  #return colinear_Hamiltonian((0.0,0.0,0.0),spin,scf_r);

  TotalOrbitalNum::Int = sum(scf_r.Total_NumOrbs[:])
  MPF = zeros(Int,scf_r.atomnum)
  #MPF = Array(Int,scf_r.atomnum)
  orbitalStartIdx = 0
  for i = 1:scf_r.atomnum
      MPF[i] = orbitalStartIdx;
      orbitalStartIdx += scf_r.Total_NumOrbs[i]
  end
  Hout = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum);
  Overlap_Band!(scf_r.Hks[spin],Hout,MPF,k_point[1],k_point[2],k_point[3],scf_r)
  return Hout;
end

function noncolinear_Hamiltonian!(Hout::Array{Complex{Float_my},2},
    H::H_type,iH::H_type,MP,k1::Float64,k2::Float64,k3::Float64,
    Hmode::DFTcommon.nc_Hamiltonian_selection,scf_r::Openmxscf)
    # Essentially same as Overlap_Band!
	#println(size(H))
    @assert(4 == size(H)[1]);
    @assert(3 == size(iH)[1]);

    k_point::Array{Float_my,1} = [k1,k2,k3];
    TotalOrbitalNum::Int = sum(scf_r.Total_NumOrbs[:])
    orbitalStartIdx::Int = 0; #각 atom별로 orbital index시작하는 지점
    MP = Array{Int}(scf_r.atomnum)
    for i = 1:scf_r.atomnum
        MP[i] = orbitalStartIdx;
        orbitalStartIdx += scf_r.Total_NumOrbs[i]
    end
    @assert(TotalOrbitalNum==orbitalStartIdx);

    # non collinear Hamiltonian: the matrix size is 2*TotalOrbitalNum,2*TotalOrbitalNum
    #Hout = zeros(Complex_my,2*TotalOrbitalNum,2*TotalOrbitalNum);
    @assert((2*TotalOrbitalNum,2*TotalOrbitalNum) == size(Hout));

    for GA_AN=1:scf_r.atomnum
        atom1_orbitalNum = scf_r.Total_NumOrbs[GA_AN];
        atom1_orbitalStart = MP[GA_AN];
        for LB_AN = 1:scf_r.FNAN[GA_AN]+1 #atom_i is not atom1,2 index
            GB_AN::UInt = scf_r.natn[GA_AN][LB_AN]
            Rn::UInt = 1+scf_r.ncn[GA_AN][LB_AN]
            atom2_orbitalNum::UInt = scf_r.Total_NumOrbs[GB_AN]
            atom2_orbitalStart::UInt = MP[GB_AN];
            kRn::Float_my = sum(scf_r.atv_ijk[Rn,:][:].*k_point);
            si::Float_my = sin(2.0*pi*kRn);
            co::Float_my = cos(2.0*pi*kRn);
            si_co = co + im*si;
            #
            for i = 1:atom1_orbitalNum
                for j = 1:atom2_orbitalNum
                    RH1 =  H[1][GA_AN][LB_AN][i][j];
                    RH2 =  H[2][GA_AN][LB_AN][i][j];
                    RH3 =  H[3][GA_AN][LB_AN][i][j];
                    RH4 =  H[4][GA_AN][LB_AN][i][j];

                    iH1 =  iH[1][GA_AN][LB_AN][i][j];
                    iH2 =  iH[2][GA_AN][LB_AN][i][j];
                    iH3 =  iH[3][GA_AN][LB_AN][i][j];
                    if DFTcommon.nc_realH_only == Hmode
                        iH1 *= 0 ;iH2 *= 0 ;iH3 *= 0 ;
                    elseif DFTcommon.nc_imagH_only == Hmode
                        RH1 *= 0 ;RH2 *= 0 ;RH3 *= 0 ;RH4 *= 0 ;
                    end

                    Hout[atom1_orbitalStart+i,atom2_orbitalStart+j] +=
                    (RH1+im*iH1) * si_co;
                    Hout[atom1_orbitalStart+i+TotalOrbitalNum,atom2_orbitalStart+j+TotalOrbitalNum] +=
                    (RH2+im*iH2) * si_co;
                    Hout[atom1_orbitalStart+i,atom2_orbitalStart+j+TotalOrbitalNum] +=
                    (RH3+im*(RH4+iH3)) * si_co;

                end
            end
        end
    end
    # set off-diagnoal part
    #Hout[TotalOrbitalNum+(1:TotalOrbitalNum),1:TotalOrbitalNum] =
    #   conj( Hout'[1:TotalOrbitalNum,TotalOrbitalNum+(1:TotalOrbitalNum)]);
    for i = (1:TotalOrbitalNum)
        for j = TotalOrbitalNum+(1:TotalOrbitalNum)
            Hout[j,i] = conj( Hout[i,j]);
            #Hout[i,j] = conj( Hout[i,j]); # test code
        end
    end
end

function noncolinear_Hamiltonian(k_point, hamiltonian_info::Hamiltonian_info_type,
  Hmode::DFTcommon.nc_Hamiltonian_selection)

  scf_r = hamiltonian_info.scf_r;

  TotalOrbitalNum = sum(scf_r.Total_NumOrbs[:])
  orbitalStartIdx_list = zeros(Int,scf_r.atomnum)
  orbitalStartIdx = 0;
  for i = 1:scf_r.atomnum
      orbitalStartIdx_list[i] = orbitalStartIdx;
      orbitalStartIdx += scf_r.Total_NumOrbs[i]
  end
  H0 = zeros(Complex_my,2*TotalOrbitalNum,2*TotalOrbitalNum)
  noncolinear_Hamiltonian!(H0,scf_r.Hks,scf_r.iHks,orbitalStartIdx_list,k_point[1],
    k_point[2],k_point[3], Hmode, scf_r);
  return H0;
end

function cal_noncolinear_eigenstate(k_point,hamiltonian_info::Hamiltonian_info_type)
  scf_r = hamiltonian_info.scf_r;
  #function nc_update_Energy(k_point_int::k_point_int_Tuple)
  # non collinear Enk and Eigen function \Psi

  #k_point = k_point_int2float(k_point_int);
  ## Common variables
  TotalOrbitalNum = sum(scf_r.Total_NumOrbs[:])

  ## Overlap matrix S
  S = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum)
  MPF = zeros(Int,scf_r.atomnum)
  orbitalStartIdx = 0;
  for i = 1:scf_r.atomnum
      MPF[i] = orbitalStartIdx;
      orbitalStartIdx += scf_r.Total_NumOrbs[i]
  end
  orbitalStartIdx_list = MPF; #TODO: use  orbitalStartIdx_list instead MPF
  Overlap_Band!(scf_r.OLP,S,MPF,k_point[1],k_point[2],k_point[3],scf_r);
  #S = S';
  # S rotation  & Merge
  if hamiltonian_info.basisTransform_rule.orbital_rot_on
    S = Heff(S,orbitalStartIdx_list,hamiltonian_info.basisTransform_rule,0.0);
  end

  ## S2 = (S^1/2)^-1
  #Sq = sqrtm(S)
  #S2 = inv(Sq);
  S2 = sqrtm_inv(S);
  S22 = zeros(Complex_my, 2*TotalOrbitalNum, 2*TotalOrbitalNum)
  up_spin = 1:TotalOrbitalNum;
  dn_spin = TotalOrbitalNum + up_spin;
  S22[up_spin,up_spin] = S2;
  S22[dn_spin,dn_spin] = S2;
  ## non-collinear Eigen funtion (S^(1/2)\Psi) and Eigen values (Enk)
  ##
  H0 = zeros(Complex_my,2*TotalOrbitalNum,2*TotalOrbitalNum)
  # H_orig is updated
  noncolinear_Hamiltonian!(H0,scf_r.Hks,scf_r.iHks,MPF,k_point[1],
  k_point[2],k_point[3],DFTcommon.nc_allH,scf_r);
  # C = M1 Ut H U M1

  if hamiltonian_info.basisTransform_rule.orbital_rot_on
    orbitals_up =   1:TotalOrbitalNum
    orbitals_down =   orbitals_up + TotalOrbitalNum
    Htmp = zeros(Complex_my,2*TotalOrbitalNum,2*TotalOrbitalNum)
    H = H0[orbitals_up,orbitals_up];
    Htmp[orbitals_up,orbitals_up] = Heff(H,orbitalStartIdx_list,hamiltonian_info.basisTransform_rule,0.0);

    H = H0[orbitals_up,orbitals_down];
    Htmp[orbitals_up,orbitals_down] = Heff(H,orbitalStartIdx_list,hamiltonian_info.basisTransform_rule,0.0);

    H = H0[orbitals_down,orbitals_up];
    Htmp[orbitals_down,orbitals_up] = Heff(H,orbitalStartIdx_list,hamiltonian_info.basisTransform_rule,0.0);

    H = H0[orbitals_down,orbitals_down];
    Htmp[orbitals_down,orbitals_down] = Heff(H,orbitalStartIdx_list,hamiltonian_info.basisTransform_rule,0.0);
    H0 = Htmp;
    #println( sum(abs(H2-H)) )
  end

  # find Eigenvalues & Eigenvector
  Hk_tilta = S22' * H0 * S22;
  Eigen_vect = copy(Hk_tilta)
  Eigen_value = zeros(Float_my,2*TotalOrbitalNum);
  eigfact_hermitian(Eigen_vect,Eigen_value);
  if (!check_eigmat(Hk_tilta,Eigen_vect,Eigen_value))
      println("H3 :",k_point)
  end
  kpoint_common = Kpoint_eigenstate(Eigen_vect,Eigen_value,k_point,Hk_tilta);
  return kpoint_common
end

function cal_noncolinear_eigenstate_depericated(k_point,hamiltonian_info::Hamiltonian_info_type)
  scf_r = hamiltonian_info.scf_r;
  #function nc_update_Energy(k_point_int::k_point_int_Tuple)
  # non collinear Enk and Eigen function \Psi

  #k_point = k_point_int2float(k_point_int);
  ## Common variables
  TotalOrbitalNum = sum(scf_r.Total_NumOrbs[:])

  ## Overlap matrix S
  S = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum)
  MPF = zeros(Int,scf_r.atomnum)
  orbitalStartIdx = 0;
  for i = 1:scf_r.atomnum
      MPF[i] = orbitalStartIdx;
      orbitalStartIdx += scf_r.Total_NumOrbs[i]
  end
  orbitalStartIdx_list = MPF; #TODO: use  orbitalStartIdx_list instead MPF
  Overlap_Band!(scf_r.OLP,S,MPF,k_point[1],k_point[2],k_point[3],scf_r);
  #S = S';
  # S rotation  & Merge
  if hamiltonian_info.basisTransform_rule.orbital_rot_on
    S = Heff(S,orbitalStartIdx_list,hamiltonian_info.basisTransform_rule,0.0);
  end

  S2 = copy(S);
  S_eigvals = zeros(Float_my, TotalOrbitalNum);
  # S^1/2
  eigfact_hermitian(S2,S_eigvals);
  if (!check_eigmat(S,S2,S_eigvals))
      println("S :",k_point)
  end
  #  S2 * 1.0/sqrt(S_eigvals[l])
  M1 = zeros(size(S_eigvals))
  M1[S_eigvals.>OLP_eigen_cut] = 1.0 ./sqrt.(S_eigvals[S_eigvals.>OLP_eigen_cut]);

  for j1 = 1:TotalOrbitalNum
      S2[:,j1] *= M1[j1];
  end
  S2 = copy(transpose(S2));

  ## non-collinear Eigen funtion (S^(1/2)\Psi) and Eigen values (Enk)
  ##
  H0 = zeros(Complex_my,2*TotalOrbitalNum,2*TotalOrbitalNum)
  H1 = zeros(Complex_my,2*TotalOrbitalNum,2*TotalOrbitalNum)
  H2 = zeros(Complex_my,2*TotalOrbitalNum,2*TotalOrbitalNum)
  # H_orig is updated
  noncolinear_Hamiltonian!(H0,scf_r.Hks,scf_r.iHks,MPF,k_point[1],
  k_point[2],k_point[3],DFTcommon.nc_allH,scf_r);
  # C = M1 Ut H U M1

  if hamiltonian_info.basisTransform_rule.orbital_rot_on
    orbitals_up =   1:TotalOrbitalNum
    orbitals_down =   orbitals_up + TotalOrbitalNum
    Htmp = zeros(Complex_my,2*TotalOrbitalNum,2*TotalOrbitalNum)
    H = H0[orbitals_up,orbitals_up];
    Htmp[orbitals_up,orbitals_up] = Heff(H,orbitalStartIdx_list,hamiltonian_info.basisTransform_rule,0.0);

    H = H0[orbitals_up,orbitals_down];
    Htmp[orbitals_up,orbitals_down] = Heff(H,orbitalStartIdx_list,hamiltonian_info.basisTransform_rule,0.0);

    H = H0[orbitals_down,orbitals_up];
    Htmp[orbitals_down,orbitals_up] = Heff(H,orbitalStartIdx_list,hamiltonian_info.basisTransform_rule,0.0);

    H = H0[orbitals_down,orbitals_down];
    Htmp[orbitals_down,orbitals_down] = Heff(H,orbitalStartIdx_list,hamiltonian_info.basisTransform_rule,0.0);
    H0 = Htmp;
    #println( sum(abs(H2-H)) )
  end

  H1 = copy(H0);
  ko_all = zeros(Float_my,2*TotalOrbitalNum);

  C = zeros(Complex_my,2*TotalOrbitalNum,2*TotalOrbitalNum)
  # H*U'
  for i1 = 1:2*TotalOrbitalNum
      for j1 = 1:TotalOrbitalNum
          sum_0 = sum(H1[i1,1:TotalOrbitalNum].*S2[j1,:])
          sum_1 = sum(H1[i1,TotalOrbitalNum+(1:TotalOrbitalNum)].*S2[j1,:])
          C[2*j1-1,i1] = sum_0;
          C[2*j1,i1] = sum_1;
      end
  end
  # U*(H*U')
  S22 = conj(S2);
  for j1 = 1:TotalOrbitalNum
      jj1 = 2*j1-1;
      jj2 = 2*j1;
      for i1 = 1:TotalOrbitalNum

          sum_00 = sum(S22[i1,:].*C[jj1,1:TotalOrbitalNum])
          sum_01 = sum(S22[i1,:].*C[jj1,TotalOrbitalNum+(1:TotalOrbitalNum)])

          sum_10 = sum(S22[i1,:].*C[jj2,1:TotalOrbitalNum])
          sum_11 = sum(S22[i1,:].*C[jj2,TotalOrbitalNum+(1:TotalOrbitalNum)])


          H2[jj1,2*i1-1] = sum_00;
          H2[jj1,2*i1] = sum_01;
          H2[jj2,2*i1-1] = sum_10;
          H2[jj2,2*i1] = sum_11;
      end
  end
  # find Eigenvalues

  H3 = copy(transpose(H2));
  #H3 = copy(H2)

  eigfact_hermitian(H3,ko_all);
  if (!check_eigmat(H2,H3,ko_all))
      println("H3 :",k_point)
  end


  H3 = copy(transpose(H3));
  S3 = copy(transpose(S2));
  NC_Es = zeros(Complex_my,2*TotalOrbitalNum,2*TotalOrbitalNum)
  NC_Es2 = zeros(Complex_my,2*TotalOrbitalNum,2*TotalOrbitalNum)
  for j1=1:2*TotalOrbitalNum
      for i1=1:TotalOrbitalNum
          sum_0 = sum(S3[i1,:].*H3[j1,2*(1:TotalOrbitalNum)-1]);
          sum_1 = sum(S3[i1,:].*H3[j1,2*(1:TotalOrbitalNum)]);
          NC_Es[j1,i1] = sum_0;
          NC_Es[j1,i1+TotalOrbitalNum] = sum_1;

          NC_Es2[j1,2*i1-1] = sum_0;
          NC_Es2[j1,i1] = sum_1;
      end
  end
  NC_Es = copy(transpose(NC_Es));
  #    Psi = Psi';
  #kpoint_nc_common = Kpoint_nc_commondata_Type(NC_Es,ko_all,k_point_int);
  #return kpoint_nc_common;
  #kpoint_nc_common = Kpoint_nc_commondata_Type(H3,ko_all,k_point_int);
  #kpoint_nc_common = Kpoint_nc_commondata_debug_Type(NC_Es,ko_all,k_point_int,S2,H0,H1,H2,H3);
  kpoint_common::Kpoint_eigenstate = Kpoint_eigenstate(NC_Es,ko_all,k_point,H0);
  return kpoint_common
end
