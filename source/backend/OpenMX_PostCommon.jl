#module OpenMXdata

using MAT
using HDF5
export update_co_linear_Energy,cal_noncolinear_eigenstate
const OLP_eigen_cut = 1.0e-10;
using DFTcommon
include("OpenMX_read_scf.jl")
#include("read_scf.jl")

#@enum Spin_mode_enum nospin=1 colinear_spin=2 noncolinear_spin=4

#used for all atoms
#global scf_r = Openmxscf

function Overlap_Band!(OLP::Overlap_type,
    S::Array{Complex{Float_my},2},MP::Array{Int32,1},k1,k2,k3,scf_r)
    #@acc begin
    k_point::Array{Float_my,1} = [k1,k2,k3];
    TotalOrbitalNum::Int32 = sum(scf_r.Total_NumOrbs[:]);
    #orbitalStartIdx::Int32 = 0; #각 atom별로 orbital index시작하는 지점

    #assert(TotalOrbitalNum==orbitalStartIdx);


    #println(TotalOrbitalNum)
    #S_size::Int32 = orbitalNum + 0;
    S1 = zeros(Float64,TotalOrbitalNum,TotalOrbitalNum);
    S2 = zeros(Float64,TotalOrbitalNum,TotalOrbitalNum);
    assert(size(S) == size(S1));

    for GA_AN=1:scf_r.atomnum

      atom1_orbitalNum = scf_r.Total_NumOrbs[GA_AN]
      atom1_orbitalStart = MP[GA_AN];
      for LB_AN = 1:scf_r.FNAN[GA_AN]+1 #atom_i is not atom1,2 index

            GB_AN::UInt32 = scf_r.natn[GA_AN][LB_AN]
            Rn::UInt32 = 1+scf_r.ncn[GA_AN][LB_AN]
            atom2_orbitalNum::UInt32 = scf_r.Total_NumOrbs[GB_AN]
            atom2_orbitalStart::UInt32 = MP[GB_AN];
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

function cal_colinear_eigenstate(k_point::k_point_Tuple,scf_r::Openmxscf,spin_list=1)

  ## Overlapmaxitrx S
  TotalOrbitalNum = sum(scf_r.Total_NumOrbs[:]);
  S = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum)
  MPF = zeros(Int32,scf_r.atomnum)
  #MPF = Array(Int32,scf_r.atomnum)
  orbitalStartIdx = 0
  for i = 1:scf_r.atomnum
      MPF[i] = orbitalStartIdx;
      orbitalStartIdx += scf_r.Total_NumOrbs[i]
  end

  Overlap_Band!(scf_r.OLP,S,MPF,k_point[1],k_point[2],k_point[3],scf_r);

  S2 = copy(S);
  S_eigvals = zeros(Float_my, TotalOrbitalNum);
  eigfact_hermitian(S2,S_eigvals);
  if(!check_eigmat(S,S2,S_eigvals))
          println("S :",k_point)
  end

  M1 = zeros(size(S_eigvals));
  M1[S_eigvals.>OLP_eigen_cut] = 1.0 ./sqrt.(S_eigvals[S_eigvals.>OLP_eigen_cut]);

  ## Coes

  #H_temp = Array(Array{Complex_my,2},scf_r.SpinP_switch+1);
  kpoint_eigenstate_list =  Array{Kpoint_eigenstate,1}();
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
      if(!check_eigmat(C2,C3,ko))
          println("C2 : ",k_point,"\tspin: ",spin)
      end
      #
      Coes = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum);
      S3 = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum);
      for ii=1:TotalOrbitalNum
          S3[:,ii] = S2[:,ii] * M1[ii];
      end
      Coes = S3*C3;
      kpoint_common::Kpoint_eigenstate = Kpoint_eigenstate(Coes,ko,k_point);
      push!(kpoint_eigenstate_list,kpoint_common);
  end
  return kpoint_eigenstate_list;
end



@enum nc_Hamiltonian_selection nc_allH=0 nc_realH_only=1 nc_imagH_only=2
function nc_Hamiltonian_basis_Transform!(Hout::Array{Complex{Float_my},2},
    H::Array{Overlap_type,1},iH::Array{Overlap_type,1},
    MP,k1,k2,k3,H_mode::nc_Hamiltonian_selection)
    # Essentially same as Overlap_Band!
    assert(4 == size(H)[1]);
    assert(3 == size(iH)[1]);

    k_point::Array{Float_my,1} = [k1,k2,k3];
    TotalOrbitalNum::Int32 = sum(scf_r.Total_NumOrbs[:])
    orbitalStartIdx::Int32 = 0; #각 atom별로 orbital index시작하는 지점
    MP = Array(Int32,scf_r.atomnum)
    for i = 1:scf_r.atomnum
        MP[i] = orbitalStartIdx;
        orbitalStartIdx += scf_r.Total_NumOrbs[i]
    end
    assert(TotalOrbitalNum==orbitalStartIdx);

    # non collinear Hamiltonian: the matrix size is 2*TotalOrbitalNum,2*TotalOrbitalNum
    #Hout = zeros(Complex_my,2*TotalOrbitalNum,2*TotalOrbitalNum);
    assert((2*TotalOrbitalNum,2*TotalOrbitalNum) == size(Hout));

    for GA_AN=1:scf_r.atomnum
        atom1_orbitalNum = scf_r.Total_NumOrbs[GA_AN];
        atom1_orbitalStart = MP[GA_AN];
        for LB_AN = 1:scf_r.FNAN[GA_AN]+1 #atom_i is not atom1,2 index
            GB_AN::UInt32 = scf_r.natn[GA_AN][LB_AN]
            Rn::UInt32 = 1+scf_r.ncn[GA_AN][LB_AN]
            atom2_orbitalNum::UInt32 = scf_r.Total_NumOrbs[GB_AN]
            atom2_orbitalStart::UInt32 = MP[GB_AN];
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
                    if nc_realH_only == H_mode
                        iH1 *= 0 ;iH2 *= 0 ;iH3 *= 0 ;
                    elseif nc_imagH_only == H_mode
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



function cal_nonco_linear_EigenState(k_point_int::k_point_int_Tuple)
    #function nc_update_Energy(k_point_int::k_point_int_Tuple)
    # non collinear Enk and Eigen function \Psi

    k_point = k_point_int2float(k_point_int);
    ## Common variables
    TotalOrbitalNum = sum(scf_r.Total_NumOrbs[:])

    ## Overlap matrix S
    S = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum)
    MPF = zeros(Int32,scf_r.atomnum)
    Overlap_Band!(scf_r.OLP,S,MPF,k_point[1],k_point[2],k_point[3]);
    #S = S';

    S2 = copy(S);
    S_eigvals = zeros(Float_my, TotalOrbitalNum);
    # S^1/2
    eigfact_hermitian(S2,S_eigvals);
    if(!check_eigmat(S,S2,S_eigvals))
        println("S :",k_point)
    end
    #  S2 * 1.0/sqrt(S_eigvals[l])
    M1 = zeros(size(S_eigvals))
    M1[S_eigvals.>OLP_eigen_cut] = 1.0 ./sqrt(S_eigvals[S_eigvals.>OLP_eigen_cut]);

    for j1 = 1:TotalOrbitalNum
        S2[:,j1] *= M1[j1];
    end
    S2 = S2.';

    ## non-collinear Eigen funtion (S^(1/2)\Psi) and Eigen values (Enk)
    ##
    H0 = zeros(Complex_my,2*TotalOrbitalNum,2*TotalOrbitalNum)
    H1 = zeros(Complex_my,2*TotalOrbitalNum,2*TotalOrbitalNum)
    H2 = zeros(Complex_my,2*TotalOrbitalNum,2*TotalOrbitalNum)
    # H_orig is updated
    nc_Hamiltonian_basis_Transform!(H0,scf_r.Hks,scf_r.iHks,MPF,k_point[1],
    k_point[2],k_point[3],nc_allH);
    # C = M1 Ut H U M1

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

    H2 = H2.' ;
    H3 = copy(H2)

    eigfact_hermitian(H3,ko_all);
    if(!check_eigmat(H2,H3,ko_all))
        println("H3 :",k_point)
    end


    H3 = H3.'
    S3 = copy(S2.');
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
    NC_Es = NC_Es.';
    #    Psi = Psi';
    kpoint_nc_common = Kpoint_nc_commondata_Type(NC_Es,ko_all,k_point_int);
    #kpoint_nc_common = Kpoint_nc_commondata_Type(H3,ko_all,k_point_int);
    #kpoint_nc_common = Kpoint_nc_commondata_debug_Type(NC_Es,ko_all,k_point_int,S2,H0,H1,H2,H3);
    return kpoint_nc_common;
end

function cal_noncolinear_eigenstate(k_point::k_point_Tuple)
  return update_nc_Energy_testing(k_point);
end

function update_nc_Energy_testing(k_point::k_point_Tuple)
  # testing version different from Opnemx Band_MO.c implimentation
  # Currently Eigen values are corrent but, LCAO coef are not matched
  # non collinear Enk and Eigen function \Psi

  #k_point = k_point_int2float(k_point_int);

  ## Common variables
  TotalOrbitalNum = sum(scf_r.Total_NumOrbs[:])

  ## Overlap matrix S
  S = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum)
  MPF = zeros(Int32,scf_r.atomnum)
  Overlap_Band!(scf_r.OLP,S,MPF,k_point[1],k_point[2],k_point[3]);
  #S = S';

  S2 = copy(S);
  S_eigvals = zeros(Float_my, TotalOrbitalNum);
  # S^1/2
  eigfact_hermitian(S2,S_eigvals);
  if(!check_eigmat(S,S2,S_eigvals))
      println("S :",k_point)
  end
  #  S2 * 1.0/sqrt(S_eigvals[l])
  M1 = zeros(size(S_eigvals))
  M1[S_eigvals.>OLP_eigen_cut] = 1.0 ./sqrt(S_eigvals[S_eigvals.>OLP_eigen_cut]);

  for j1 = 1:TotalOrbitalNum
      S2[:,j1] *= M1[j1];
  end

  S_2by2 = zeros(Complex_my,2*TotalOrbitalNum,2*TotalOrbitalNum);

  S_2by2[1:TotalOrbitalNum,1:TotalOrbitalNum] = S2;
  S_2by2[TotalOrbitalNum+(1:TotalOrbitalNum),TotalOrbitalNum+(1:TotalOrbitalNum)] = S2;

  ## non-collinear Eigen funtion (S^(1/2)\Psi) and Eigen values (Enk)
  ##
  H0 = zeros(Complex_my,2*TotalOrbitalNum,2*TotalOrbitalNum)
  # H_orig is updated
  nc_Hamiltonian_basis_Transform!(H0,scf_r.Hks,scf_r.iHks,MPF,k_point[1],
  k_point[2],k_point[3],nc_allH);
  # C = M1 Ut H U M1

  Enks = zeros(Float_my,2*TotalOrbitalNum);
  H2 = S_2by2' * H0 * S_2by2;

  H3 = copy(H2)
  eigfact_hermitian(H3,Enks);
  if(!check_eigmat(H2,H3,Enks))
      println("H3 :",k_point)
  end
  NC_Es = S_2by2' * H3;

  #kpoint_nc_common = Kpoint_nc_commondata_Type(NC_Es,ko_all,k_point_int);
  #kpoint_nc_common = Kpoint_nc_commondata_Type(H3,Enks,k_point_int);
  #kpoint_nc_common = Kpoint_nc_commondata_debug_Type(NC_Es,Enks,k_point_int,S_2by2,H0,H0,H2,H3);
  #return kpoint_nc_common;
	kpoint_common::Kpoint_eigenstate = Kpoint_eigenstate(NC_Es,Enks,k_point);
    return kpoint_common;
end