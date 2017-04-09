module DFTcommon
using ArgParse
using ProgressMeter
import TOML # Pkg.clone("https://github.com/wildart/TOML.jl.git")

export Kpoint_eigenstate,Complex_my,Float_my,k_point_Tuple,k_point_int_Tuple
export Hamiltonian_type
export k_point_int2float,k_point_float2int,kPoint2BrillouinZone_Tuple,
       kPoint2BrillouinZone_int_Tuple,kPoint_gen_GammaCenter,q_k_unique_points
export pwork
export k_point_precision

export parse_input,Arg_Inputs


export eigfact_hermitian,check_eigmat

export Hartree2cm,cm2meV,cm2meV,Hartree2meV,Hatree2eV,kB
const Hartree2cm = 2.194746*100000.0;
const cm2meV = 0.1240;
const Hartree2meV = Hartree2cm*cm2meV;
const Hatree2eV = 27.2114;
const kB=0.000003166813628;  # Boltzman constant (Hatree/K)

const k_point_precision = 10.0^6;


Float_my =  Float64
Complex_my = Complex128

k_point_Tuple = Tuple{Float64,Float64,Float64};
k_point_int_Tuple = Tuple{Int64,Int64,Int64};
k_point_rational_Tuple = Tuple{Rational{Int64},Rational{Int64},Rational{Int64}};

type Kpoint_eigenstate
    Eigenstate::Array{Complex_my,2};
    Eigenvalues::Array{Float_my,1};
    k_point_int::k_point_Tuple;
end
Hamiltonian_type = Array{Complex_my,2};


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

  for kx in (1:k_point_num[1])/(1+k_point_num[1]) - 1/2
      for ky in (1:k_point_num[2])/(1+k_point_num[2]) - 1/2
          for kz in (1:k_point_num[3])/(1+k_point_num[3]) - 1/2
              push!(k_point_list,
                  kPoint2BrillouinZone_Tuple((kx,ky,kz)) );
          end
      end
  end
  return k_point_list;
end

function q_k_unique_points(q_point_list,k_point_list)
  kq_point_dict = Dict{k_point_int_Tuple,k_point_Tuple}();
  p = Progress( length(q_point_list),"Computing unique K points from K,K+Q... ");
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

@enum orbital_mask_enum nomask=0 unmask=1 mask=2
type orbital_mask_input_Type
    orbital_mask1::Array{Int64,1}
    orbital_mask2::Array{Int64,1}
    #atom12::atom12_Tuple
    orbital_mask_on::Bool
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
    if(maximum(abs.(check_eig)) > 10.0^-3)
      println(" sum ",sum(check_eig)," max ",maximum(abs(check_eig))," Diff ",sum(diff_A),
      " Diff realmax ",maximum(real(diff_A)),
      " Diff imgmax ",maximum(imag(diff_A)));
        return false
    end
    #assert( maximum(abs(check_eig)) < 10.0^-3);
    return true;
end

################################################################################
# Input parse
################################################################################
type Arg_Inputs
  scf_name::AbstractString
  atom1::Int
  atom2::Int
  k_point_num::Array{Int,1}
  q_point_num::Array{Int,1}
  orbital_mask1
  orbital_mask2
  orbital_mask_option::orbital_mask_enum
  orbital_mask_name
  hdftmpdir
  ChemP_delta_ev::Float64
  Arg_Inputs() = new("",-1,-1,[2,2,2],[2,2,2],
    Array{Int64,1}(),Array{Int64,1}(),nomask,"All","",0.0)
end

function read_toml(toml_fname)
  toml_inputs =  TOML.parse(readstring(toml_fname));

end


function parse_input(args)

    input::Arg_Inputs =  Arg_Inputs()

    s = ArgParseSettings("Example 2 for argparse.jl: " *  # description
                      "flags, options help, " *
                      "required arguments.")

    @add_arg_table s begin
        "--TOMLinput","-T"
        help = "input.toml file ex:) nio.toml "
        "--DFTtype","-D"
        help = "openmx, openmxWF, vaspWF "
        "--atom12"
        help = "target atom1&2 ex:) 1_2"     # used by the help screen
        "--kpoint", "-k"
        #action = :store_true   # this makes it a flag
        help = "k_point ex:) 5_5_5"
        "--qpoint", "-q"
        #action = :store_true   # this makes it a flag
        help = "q_point ex:) 5_5_5"
        "--om1"
        help = "orbital_mask1 ex) 1,2,3"
        "--om2"
        help = "orbital_mask1 ex) 1,2,3"
        "--ommode"
        help = "0 (default: no orbital control) 1(selected orbitals are unmasked) 2(selected orbitals are masked)"
        arg_type = Int
        "--omname"
        help = "obital masking name ex) d_d "
        "--hdftmpdir"
        help = "specify hdf shared mem tmpdir dir (should be visible to all processes). Default: root_dir/jq/*.hdf5"
        "--chempdelta"
        help = "Chemical potential shift in eV(default 0.0 eV)"
        arg_type = Float64
        "scfname"
        help = "scf file name ex:) nio.scfout"
        required = true        # makes the argument mandatory
    end

    parsed_args = parse_args(args, s)
    println("Parsed args:",parsed_args)
    for (key,val) in parsed_args
        #println("  $key  =>  $(repr(val))")
        #println("  $key  =>  $val")

        if(key == "atom12" && (Void != typeof(val)))
            atom12 =   map(x->parse(Int64,x),split(val,"_"))
            input.atom1 = atom12[1];
            input.atom2 = atom12[2];
        end
        if(key == "scfname")
            if(isfile(val) && ".scfout" == splitext(val)[2])
                input.scf_name  = val;
                #println("scf file:$scf_name")
            end
        end
        if(key == "kpoint" && typeof(val) <: AbstractString)
            k_point_num_tmp =   map(x->parse(Int64,x),split(val,"_"))
            if(3 == length(k_point_num_tmp))
                input.k_point_num = k_point_num_tmp;
            else
                println("k point should be like 10_10_10")
            end
        end
        if(key == "qpoint" && typeof(val) <: AbstractString)
            q_point_num_tmp =   map(x->parse(Int64,x),split(val,"_"))
            if(3 == length(q_point_num_tmp))
                input.q_point_num = q_point_num_tmp;
            else
                println("q point should be like 10_10_10")
            end
        end
        if(key =="om1" && typeof(val) <: AbstractString)
            println(val)
            input.orbital_mask1 = parse_int_list(val)
        end
        if(key =="om2" && typeof(val) <: AbstractString)
            input.orbital_mask2 = parse_int_list(val)
        end
        if("ommode"==key)
            if(1 == val)
               input.orbital_mask_option = unmask
            elseif(2 == val)
               input.orbital_mask_option = mask
            end
        end
        if("omname"==key)
            println(val)
            input.orbital_mask_name = val;
        end
        if("hdftmpdir" ==key)
            if (nothing!=val)
                if isdir(val)
                    input.hdftmpdir = val;
                else
                    println(" $val is not directory")
                end
            end
        end
        if("chempdelta" == key)
            if(nothing != val)
                input.ChemP_delta_ev = val;
                ChemP_delta = val/Hatree2eV;
            end
        end
    end
    return input;
end


end
