module DFTcommon
using ArgParse
using ProgressMeter
import TOML # Pkg.clone("https://github.com/wildart/TOML.jl.git")

export Kpoint_eigenstate,Kpoint_eigenstate_only
export Complex_my,Float_my,k_point_Tuple,k_point_int_Tuple
export Hamiltonian_type,Hamiltonian_info_type
export k_point_int2float,k_point_float2int,kPoint2BrillouinZone_Tuple,
       kPoint2BrillouinZone_int_Tuple,q_k_unique_points
export kPoint_gen_EquallySpaced,kPoint_gen_GammaCenter
export pwork
export k_point_precision

export orbital_mask_input_Type,orbital_mask_enum
export parse_input,Arg_Inputs,parse_TOML,read_toml
export kPath_band,Arg_Inputs_Band


export eigfact_hermitian,check_eigmat
export SPINtype, DFTtype, Wannier90type



bar_string = "================================================================"
export bar_string;

@enum SPINtype para_type = 1 colinear_type = 2 non_colinear_type = 4
@enum DFTtype OpenMX = 1 Wannier90 = 2 NULLDFT = -1
@enum Wannier90type OpenMXWF = 1 VASPWF = 2 EcalJWF = 3  NULLWANNIER = -1

export Hartree2cm,cm2meV,cm2meV,Hartree2meV,Hatree2eV,kBeV
const Hartree2cm = 2.194746*100000.0;
const cm2meV = 0.1240;
const Hartree2meV = Hartree2cm*cm2meV;
const Hatree2eV = 27.21138602;

#const kB=0.000003166813628;  # Boltzmann constant (Hatree/K)
const kBeV = 8.6173303*(10.0^-5);  # Boltzmann constant (eV/K)
export ang2bohr
const ang2bohr = 1.889725989;
const k_point_precision = 10.0^6;
export nc_Hamiltonian_selection
@enum nc_Hamiltonian_selection nc_allH=0 nc_realH_only=1 nc_imagH_only=2

Float_my =  Float64
Complex_my = Complex128

k_point_Tuple = Tuple{Float64,Float64,Float64};
k_point_int_Tuple = Tuple{Int64,Int64,Int64};
k_point_rational_Tuple = Tuple{Rational{Int64},Rational{Int64},Rational{Int64}};
@enum orbital_mask_enum nomask=0 unmask=1 mask=2

include("basisTransform.jl")
include("inputHandler.jl")

Hamiltonian_type = Array{Complex_my,2};
type Kpoint_eigenstate_only
    Eigenstate::Array{Complex_my,2};
    Eigenvalues::Array{Float_my,1};
    k_point_int::k_point_Tuple;
    #Hamiltonian::Hamiltonian_type;
end
type Kpoint_eigenstate
    Eigenstate::Array{Complex_my,2};
    Eigenvalues::Array{Float_my,1};
    k_point_int::k_point_Tuple;
    Hamiltonian::Hamiltonian_type;
end
type Hamiltonian_info_type
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
    return (rem(k_point+2.5,1)-0.5);
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
    k_point_int_array =  (rem(k_point_int_array+(2.5*k_point_precision),k_point_precision)-0.5*k_point_precision);
    k_point_int_array = round(Int64,k_point_int_array);
    return (k_point_int_array[1],k_point_int_array[2],k_point_int_array[3])
end

function kPoint_gen_GammaCenter(k_point_num)
  k_point_list = Array(k_point_Tuple,0);

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
  k_point_list = Array(k_point_Tuple,0);
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
  p.output = STDOUT
  for q_point in q_point_list
      for k_point in k_point_list
        kq_point = (q_point[1] + k_point[1],q_point[2] + k_point[2],q_point[3] + k_point[3]) ;
        kq_point = kPoint2BrillouinZone_Tuple(kq_point);
        kq_point_int = k_point_float2int(kq_point);
        kq_point_dict[kq_point_int] = kq_point;
      end
    next!(p);
  end
  gc();
  kq_point_int_list = collect(keys(kq_point_dict));
  kq_point_list = collect(values(kq_point_dict));
  return (kq_point_list,kq_point_int_list);
end


atom12_Tuple = Tuple{Int64,Int64};
type orbital_mask_input_Type
    orbital_mask1::Array{Int64,1}
    orbital_mask2::Array{Int64,1}

    orbital_mask3::Array{Int64,1}
    orbital_mask4::Array{Int64,1}
    atom12::atom12_Tuple
    orbital_mask_on::Bool

    function orbital_mask_input_Type(orbital_mask1::Array{Int64,1}, orbital_mask2::Array{Int64,1},
      atom12::atom12_Tuple,orbital_mask_on::Bool)

      new(orbital_mask1,orbital_mask2,Array{Int64}(0),Array{Int64}(0),
        atom12,orbital_mask_on)
    end
    function orbital_mask_input_Type(orbital_mask1::Array{Int64,1}, orbital_mask2::Array{Int64,1},
      orbital_mask3::Array{Int64,1}, orbital_mask4::Array{Int64,1},
      atom12::atom12_Tuple,orbital_mask_on::Bool)

      new(orbital_mask1,orbital_mask2,orbital_mask3,orbital_mask4,
        atom12,orbital_mask_on)
    end
end

function orbital_mask_inv(orbital_mask1,atom1_orbitalNum,)
  orbital_mask1_inv = Array{Int64,1}();
  if (nomask == orbital_mask_option)
      if ( 0 < length(orbital_mask1))
          println("INFO: Orbital masking options only works with --ommode 1 or 2")
      end
      orbital_mask1 = Array{Int64,1}();
      orbital_mask_name = "";
  else
      orbital_mask_on = true;
      if (unmask ==  orbital_mask_option) # inverse selection

          orbital_mask1_tmp = collect(1:atom1_orbitalNum);

          for orbit1 in orbital_mask1
              deleteat!(orbital_mask1_tmp, find(orbital_mask1_tmp.==orbit1))
          end

          orbital_mask1 = orbital_mask1_tmp;

      end
      ## orbital_mask1_inv are only used for file saving
      for orbit1 = 1:atom1_orbitalNum
          if (0==length(find( orbital_mask1.== orbit1)))
              push!(orbital_mask1_inv,orbit1);
          end
      end

  end
  return orbital_mask1_inv;
end


function pwork(f,args)
    worker_list = workers()
    n = length(worker_list)
    results = Array{Any}(n)
    @sync begin
        for (worer_idx,worker_id) in enumerate(worker_list)
            @async begin
            #@everywhere init_all(atom1,atom2)
                results[worer_idx] = remotecall_fetch(f,worker_id,args)
            end
        end
    end
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

    temp = eigfact(EigVect_h);
    p = sortperm(temp.values);
    #EigVal[:] = real(temp.values[p])[:];

    EigVal[:] = (temp.values[p])[:];
    EigVect[:] = temp.vectors[:,p][:];
    #ccall((:EigenBand_lapack3,"./lib_jx"),Void,(Ptr{Complex128},Ptr{Float64}, Int32,)
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
    if (maximum(abs(check_eig)) > 10.0^-3)
      println(" sum ",sum(check_eig)," max ",maximum(abs(check_eig))," Diff ",sum(diff_A),
      " Diff realmax ",maximum(real(diff_A)),
      " Diff imgmax ",maximum(imag(diff_A)));
        return false
    end
    #assert( maximum(abs(check_eig)) < 10.0^-3);
    return true;
end

end
