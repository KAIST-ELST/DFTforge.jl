###############################################################################
# Hongkee Yoon Hongkeeyoon@kaist.ac.kr
# 2019.05
# https://kaist-elst.github.io/DFTforge.jl/
###############################################################################

__precompile__(true)
module DFTcommon
using Distributed

using LinearAlgebra

using ProgressMeter
using Distributed
using Statistics
export ArgParse,ProgressMeter,Distributed,Statistics


#include("../ext/TOML/src/TOML.jl")
#import ..TOML # Pkg.clone("https://github.com/wildart/TOML.jl.git")
#import DFTforge.TOML
using ..TOML

export Kpoint_eigenstate,Kpoint_eigenstate_only
export Complex_my,Float_my,k_point_Tuple,k_point_int_Tuple
export Hamiltonian_type,Hamiltonian_info_type
export k_point_int2float,k_point_float2int,kPoint2BrillouinZone_Tuple,
       kPoint2BrillouinZone_int_Tuple,q_k_unique_points
export kPoint_gen_EquallySpaced,kPoint_gen_GammaCenter
export pwork
export k_point_precision

export orbital_selection_input_Type,orbital_selection_enum
export parse_input,Arg_Inputs,parse_TOML,read_toml
export kPath_band,Arg_Inputs_Band


export eigfact_hermitian,check_eigmat, sqrtm_inv, cal_eigenVectVal
export SPINtype, DFTtype, Wannier90type



################################################################################
# DMFT + DFT
################################################################################
export excute_cmd
export Arg_DMFT_DFT,Arg_MPI

bar_string = "================================================================"
export bar_string;

################################################################################



# #=
@enum SPINtype para_type = 1 colinear_type = 2 non_colinear_type = 4
@enum DFTtype OpenMX = 1 Wannier90 = 2 EcalJ = 3 NULLDFT = -1
@enum Wannier90type OpenMXWF = 1 Wannier90WF = 2 EcalJWF = 3  NULLWANNIER = -1
@enum ChargeAnalysis Lowdin = 1 NAO = 2

export Hartree2cm,cm2meV,cm2meV,Hartree2meV,Hatree2eV,Ry2eV,kBeV
const Hartree2cm = 2.194746*100000.0;
const cm2meV = 0.1240;
const Hartree2meV = Hartree2cm*cm2meV;
const Hatree2eV = 27.21138602;
const Ry2eV = Hatree2eV/2

#const kB=0.000003166813628;  # Boltzmann constant (Hatree/K)
const kBeV = 8.6173303*(10.0^-5);  # Boltzmann constant (eV/K)
export ang2bohr
const ang2bohr = 1.889725989;
const k_point_precision = 10.0^6;
export nc_Hamiltonian_selection
@enum nc_Hamiltonian_selection nc_allH=0 nc_realH_only=1 nc_imagH_only=2

Float_my =  Float64
Complex_my = ComplexF64  #ComplexF64

k_point_Tuple = Tuple{Float64,Float64,Float64};
k_point_int_Tuple = Tuple{Int64,Int64,Int64};
k_point_rational_Tuple = Tuple{Rational{Int64},Rational{Int64},Rational{Int64}};
@enum orbital_selection_enum nomask=0 unmask=1 mask=2

export mpirun_type, excute_cmd, excute_mpi_cmd
@enum mpirun_type mvapich = 1 openmpi = 2


Hamiltonian_type = Array{Complex_my,2};
# =#
export nc_Hamiltonian_selection


include("basisTransform.jl")
export basisTransform_result_type
#using .basisTransform

include("inputHandler.jl")
#using .inputHandler


struct Kpoint_eigenstate_only
    Eigenstate::Array{Complex_my,2};
    Eigenvalues::Array{Float_my,1};
    k_point_int::k_point_Tuple;
    #Hamiltonian::Hamiltonian_type;
end
struct Kpoint_eigenstate
    Eigenstate::Array{Complex_my,2}; #[orbital index, energy index]
    Eigenvalues::Array{Float_my,1};
    k_point_int::k_point_Tuple;
    Hamiltonian::Hamiltonian_type;
end
struct Hamiltonian_info_type
  scf_r;
  dfttype::DFTcommon.DFTtype
  spin_type::SPINtype
  #atomnum_eff::Int
  #orbitalNums_eff::Array{Int}
  #orbital_index_map::Dict{Int,Dict{Int,Int}};
  basisTransform_result::basisTransform_result_type

  basisTransform_rule::basisTransform_rule_type
  function Hamiltonian_info_type(scf_r,dfttype::DFTtype,spin_type::SPINtype)
    atomnum::Int = copy(scf_r.atomnum);
    orbitalNums::Array{Int} = copy(scf_r.Total_NumOrbs);
    basisTransform_rule = basisTransform_rule_type()
    basisTransform_result = basisTransform_init(atomnum,orbitalNums,basisTransform_rule)
    new(scf_r,dfttype,spin_type,
    basisTransform_result,basisTransform_rule);
  end

  function Hamiltonian_info_type(scf_r,dfttype::DFTtype,spin_type::SPINtype,
    basisTransform_rule::basisTransform_rule_type)
    atomnum::Int = copy(scf_r.atomnum);
    orbitalNums::Array{Int} = copy(scf_r.Total_NumOrbs);
    basisTransform_result = basisTransform_init(atomnum,orbitalNums,basisTransform_rule)
    new(scf_r,dfttype,spin_type,
    basisTransform_result,basisTransform_rule);
  end
end





function k_point_int2float(k_point_int::k_point_int_Tuple)
    k_point = k_point_Tuple;
    k_point = (k_point_int[1]/k_point_precision,k_point_int[2]/k_point_precision,
    k_point_int[3]/k_point_precision);
    return k_point;
end
function k_point_float2int(k_point::k_point_Tuple)
    #k_point_int = k_point_int_Tuple;
    k_point_int = (round(Int64,k_point[1]*k_point_precision),round(Int64,k_point[2]*k_point_precision),
    round(Int64,k_point[3]*k_point_precision) );
    return k_point_int;
end
function kPoint2BrillouinZone(k_point::Array{Float64,1})
    # from -0.5 0.5
    # 0.5 is treated as -0.5
    return (rem.(k_point .+ 2.5,1) .- 0.5);
end

function kPoint2BrillouinZone_Tuple(k_point::k_point_Tuple)
    # from -0.5 0.5
    # 0.5 is treated as -0.5
    k_point_array =    kPoint2BrillouinZone([k_point[1],k_point[2],k_point[3]]);
    return (k_point_array[1],k_point_array[2],k_point_array[3])
end
function kPoint2BrillouinZone_int_Tuple(k_point_int::k_point_int_Tuple)
    # from -0.5 0.5
    # 0.5 is treated as -0.5
    k_point_int_array = [k_point_int[1],k_point_int[2],k_point_int[3]]
    k_point_int_array =  (rem.(k_point_int_array+(2.5*k_point_precision),k_point_precision)-0.5*k_point_precision);
    k_point_int_array = round(Int64,k_point_int_array);
    return (k_point_int_array[1],k_point_int_array[2],k_point_int_array[3])
end

function kPoint_gen_GammaCenter(k_point_num)
  k_point_list = Array{k_point_Tuple}(undef,0);

  for kx in (0:(k_point_num[1]-1))/(k_point_num[1]) #- 1/2
      for ky in (0:(k_point_num[2]-1))/(k_point_num[2]) #- 1/2
          for kz in (0:(k_point_num[3]-1))/(k_point_num[3]) #- 1/2
              push!(k_point_list,
                  kPoint2BrillouinZone_Tuple((kx,ky,kz)) );
          end
      end
  end
  k_point_list = sort(unique(k_point_list));
  return k_point_list;
end
function kPoint_gen_EquallySpaced(k_point_num)
  k_point_list = Array{k_point_Tuple}(undef,0);
  kPoint_esp = 10.0^-8;
  for kx in (0.0:1/k_point_num[1]:(1.0-kPoint_esp))
    for ky in (0.0:1/k_point_num[2]:(1.0-kPoint_esp))
      for kz in (0.0:1/k_point_num[3]:(1.0-kPoint_esp))
        push!(k_point_list,
            kPoint2BrillouinZone_Tuple((kx,ky,kz)) );
      end
    end
  end
  k_point_list = sort(unique(k_point_list));
  return k_point_list;
end

function q_k_unique_points(q_point_list,k_point_list)
  kq_point_dict = Dict{k_point_int_Tuple,k_point_Tuple}();
  p = Progress( length(q_point_list),"Computing unique K points from K,K+Q... ");
  p.barglyphs=BarGlyphs("[=> ]")
  p.output = stdout
  for q_point in q_point_list
      for k_point in k_point_list
        kq_point = (q_point[1] + k_point[1],q_point[2] + k_point[2],q_point[3] + k_point[3]) ;
        kq_point = kPoint2BrillouinZone_Tuple(kq_point);
        kq_point_int = k_point_float2int(kq_point);
        kq_point_dict[kq_point_int] = kq_point;
      end
    next!(p);
  end
  GC.gc();
  kq_point_int_list = collect(keys(kq_point_dict));
  kq_point_list = collect(values(kq_point_dict));
  return (kq_point_list,kq_point_int_list);
end


atom12_Tuple = Tuple{Int64,Int64};
struct orbital_selection_input_Type
    orbital_selection1::Array{Int64,1}
    orbital_selection2::Array{Int64,1}

    orbital_selection3::Array{Int64,1}
    orbital_selection4::Array{Int64,1}
    atom12::atom12_Tuple
    orbital_selection_on::Bool

    function orbital_selection_input_Type(orbital_selection1::Array{Int64,1}, orbital_selection2::Array{Int64,1},
      atom12::atom12_Tuple,orbital_selection_on::Bool)

      new(orbital_selection1,orbital_selection2,Array{Int64}(undef,0),Array{Int64}(undef,0),
        atom12,orbital_selection_on)
    end
    function orbital_selection_input_Type(orbital_selection1::Array{Int64,1}, orbital_selection2::Array{Int64,1},
      orbital_selection3::Array{Int64,1}, orbital_selection4::Array{Int64,1},
      atom12::atom12_Tuple,orbital_selection_on::Bool)

      new(orbital_selection1,orbital_selection2,orbital_selection3,orbital_selection4,
        atom12,orbital_selection_on)
    end
end

function orbital_selection_inv(orbital_selection1,atom1_orbitalNum,)
  orbital_selection1_inv = Array{Int64,1}();
  if (nomask == orbital_selection_option)
      if ( 0 < length(orbital_selection1))
          println("INFO: Orbital masking options only works with --ommode 1 or 2")
      end
      orbital_selection1 = Array{Int64,1}();
      orbital_selection_name = "";
  else
      orbital_selection_on = true;
      if (unmask ==  orbital_selection_option) # inverse selection

          orbital_selection1_tmp = collect(1:atom1_orbitalNum);

          for orbit1 in orbital_selection1
              deleteat!(orbital_selection1_tmp, find(orbital_selection1_tmp.==orbit1))
          end

          orbital_selection1 = orbital_selection1_tmp;

      end
      ## orbital_selection1_inv are only used for file saving
      for orbit1 = 1:atom1_orbitalNum
          if (0==length(find( orbital_selection1.== orbit1)))
              push!(orbital_selection1_inv,orbit1);
          end
      end

  end
  return orbital_selection1_inv;
end


function pwork(f,args...)
    worker_list = workers()
    n = length(worker_list)
    results = Array{Any}(undef,n)
    @sync begin
        for (worer_idx,worker_id) in enumerate(worker_list)
            @async begin
            #@everywhere init_all(atom1,atom2)
                results[worer_idx] = remotecall_fetch(f,worker_id,args...)
            end
        end
    end
    return results
end
################################################################################
# Math related function
################################################################################
function eigfact_hermitian(EigVect::Array{Complex_my,2},
    EigVal::Array{Float_my,1})
    for ii = 1:size(EigVect)[1]
        EigVect[ii,ii] = real(EigVect[ii,ii]);
    end
    #EigVect_h = Hermitian(EigVect,:L);
    EigVect_h = Hermitian(EigVect);

    temp = eigen(EigVect_h);
    p = sortperm(temp.values);
    #EigVal[:] = real(temp.values[p])[:];

    EigVal[:] = (temp.values[p])[:];
    EigVect[:] = temp.vectors[:,p][:];
    #ccall((:EigenBand_lapack3,"./lib_jx"),Void,(Ptr{ComplexF64},Ptr{Float64}, Int32,)
    #,EigVect,EigVect,size(EigVect)[1])
end
function check_eigmat(A1,eig_mat,eig_val) #A1 orign  eig_mat # eig_val
    n = size(A1)[1];
    check_eig = zeros(Complex_my,n,)
    for ii = 1:n
        check_eig[ii] =  sum((A1*eig_mat[:,ii])[:] - (eig_mat[:,ii]*eig_val[ii])[:] )
    end
        #A1 = (A1+conj(A1'))/2;
    diff_A = ((A1+conj(A1'))/2 - real(A1)) ;

    #println("")
    if (maximum(abs.(check_eig)) > 10.0^-3)
      println(" sum ",sum(check_eig)," max ",maximum(abs.(check_eig))," Diff ",sum(diff_A),
      " Diff realmax ",maximum(real(diff_A)),
      " Diff imgmax ",maximum(imag(diff_A)));
        return false
    end
    #@assert( maximum(abs(check_eig)) < 10.0^-3);
    return true;
end

function sqrtm_inv(S)
    # S3 = inv(sqrtm(S))
    # Filter small Eigen values(M1) for safe inv
    #const
    OLP_eigen_cut = 10.0^-10.0
    n = size(S)[1]

    S2 = copy(S);
    S_eigvals = zeros(n);
    eigfact_hermitian(S2,S_eigvals);


    M1 = zeros(size(S_eigvals))
    M1[S_eigvals.>OLP_eigen_cut] = 1.0 ./sqrt.(S_eigvals[S_eigvals.>OLP_eigen_cut]);

    M1= 1.0 ./sqrt.(S_eigvals);

    M2 = zeros(n,n)
    for j1 = 1:n
        M2[j1,j1] = M1[j1];
    end

    S3 = S2 * M2 * S2'
    return S3
end

@inline function cal_eigenVectVal(Hk::Array{Complex{Float64},2},S::Array{Complex{Float64},2}; OLP_eigen_cut = 1.0E-7 )
    n = size(S)[1]

    U = copy(S);
    S_eigvals = zeros(n);
    eigfact_hermitian(U,S_eigvals);

    #ni_start =  findfirst(OLP_eigen_cut .< S_eigvals ); # find idx to cut ( if ni_start == 1, nothing should happen)
    M1 = zeros(size(S_eigvals))
    #M1[S_eigvals.>OLP_eigen_cut] = 1.0 ./sqrt.(S_eigvals[S_eigvals.>OLP_eigen_cut]);
    S_eigvals[ S_eigvals .< OLP_eigen_cut ] .= OLP_eigen_cut/2; # negtiave eigen value filter
    selected_S_mask = OLP_eigen_cut .< S_eigvals

    M1= 1.0 ./sqrt.(S_eigvals);

    #S2 = zeros(Complex_my, n,n2)

    S2 = zeros(Complex_my, n, n ) #sum(selected_S_mask))
    for j1 = 1:n # ni_start:n
        if OLP_eigen_cut < S_eigvals[j1]
            S2[:,j1] = M1[j1] * U[:,j1]
        else
            # println(S_eigvals[j1])
        end
    end
    S2 = S2[:,selected_S_mask]

    Hk_tilta_zz = S2' * (Hk * S2)
    Eigen_vect = copy(Hk_tilta_zz)
    Eigen_value = zeros( size(S2)[2]);

    eigfact_hermitian(Eigen_vect,Eigen_value);
    #Psi = U * Eigen_vect; #[orbital index , Enk]
    Psi = U[:,selected_S_mask] * Eigen_vect;
    Hk_tilta_orbitalbasis = U[:,selected_S_mask] * Hk_tilta_zz * adjoint(U[:,selected_S_mask])
    return Psi, Eigen_value, Hk_tilta_orbitalbasis #, U[:,selected_S_mask]
    #return Eigen_vect,Eigen_value;
end

################################################################################
# Excute cmd
################################################################################

function excute_cmd(cmd::String,options,work_dir)
  # TODO: support multiple argument
  run_error_check = false;
  println(cmd)
  try
    cd(work_dir);
    run(`$cmd $options`)
  catch ee
    run_error_check = true;
    println(ee)
  end
  if (run_error_check)
    println("=================================================================")
    println(" Error occuured while:")
    println( cmd)
    println("=================================================================")
    @assert(falses)
  end
end

function excute_mpi_cmd(MPI_paramter::Arg_MPI, cmd::String,
  options,work_dir, mpi_num_core::Int64=-1)
  run_error_check = false;
  println(cmd)
  mpirun = MPI_paramter.mpirun;

  mpi_options = ["-np", MPI_paramter.mpirun_num]
  thread_options = ["-nt", MPI_paramter.thread_num]

  if (isfile(MPI_paramter.mpirun_hostfile))
	  println("mpi file exist")
    if (mvapich == MPI_paramter.mpirun_type)
	    if(mpi_num_core>0)
      mpi_options = ["-f", MPI_paramter.mpirun_hostfile, "-np",mpi_num_core]
	    else
      mpi_options = ["-f", MPI_paramter.mpirun_hostfile]
	    end
    elseif (openmpi == MPI_paramter.mpirun_type)
      mpi_options = ["--hostfile", MPI_paramter.mpirun_hostfile]
    end
  end
  try
    cd(work_dir);
    run(`$mpirun $mpi_options $cmd $options $thread_options`)
  catch ee
    run_error_check = true;
    println(ee)
  end
  if (run_error_check)
    prinltn("=================================================================")
    println(" Error occuured while:")
    println( cmd)
    prinltn("=================================================================")
    @assert(falses)
  end
end

end
