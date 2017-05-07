module DFTcommon
using ArgParse
using ProgressMeter
import TOML # Pkg.clone("https://github.com/wildart/TOML.jl.git")

export Kpoint_eigenstate,Complex_my,Float_my,k_point_Tuple,k_point_int_Tuple
export Hamiltonian_type
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

const k_point_precision = 10.0^6;
export nc_Hamiltonian_selection
@enum nc_Hamiltonian_selection nc_allH=0 nc_realH_only=1 nc_imagH_only=2

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
atom12_Tuple = Tuple{Int64,Int64};
type orbital_mask_input_Type
    orbital_mask1::Array{Int64,1}
    orbital_mask2::Array{Int64,1}
    atom12::atom12_Tuple
    orbital_mask_on::Bool
end

function orbital_mask_inv(orbital_mask1,atom1_orbitalNum,)
  orbital_mask1_inv = Array{Int64,1}();
  if(nomask == orbital_mask_option)
      if( 0 < length(orbital_mask1))
          println("INFO: Orbital masking options only works with --ommode 1 or 2")
      end
      orbital_mask1 = Array{Int64,1}();
      orbital_mask_name = "";
  else
      orbital_mask_on = true;
      if(unmask ==  orbital_mask_option) # inverse selection

          orbital_mask1_tmp = collect(1:atom1_orbitalNum);

          for orbit1 in orbital_mask1
              deleteat!(orbital_mask1_tmp, find(orbital_mask1_tmp.==orbit1))
          end

          orbital_mask1 = orbital_mask1_tmp;

      end
      ## orbital_mask1_inv are only used for file saving
      for orbit1 = 1:atom1_orbitalNum
          if(0==length(find( orbital_mask1.== orbit1)))
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

################################################################################
# Input parse
################################################################################
type Wannier_OptionalInfo
  atomnum::Int
  atompos::Array{Float64,2}
  atoms_orbitals_list::Array{Array{Int},1}
  Wannier_OptionalInfo() = new(0,zeros(0,3),Array{Array{Int},1}(0));
end

type Arg_Inputs
  result_file::AbstractString
  atom1::Int
  atom2::Int
  atom12_list::Vector{Tuple{Int64,Int64}}
  k_point_num::Array{Int,1}
  q_point_num::Array{Int,1}
  orbital_mask1_list::Array{Array{Int}}
  orbital_mask1_names::Array{AbstractString,1}
  orbital_mask2_list::Array{Array{Int}}
  orbital_mask2_names::Array{AbstractString,1}
  orbital_mask_option::orbital_mask_enum
  #orbital_mask_name


  Hmode::nc_Hamiltonian_selection
  hdftmpdir
  ChemP_delta_ev::Float64
  TOMLinput::AbstractString
  k_point_list::Vector{k_point_Tuple};
  DFT_type::DFTtype
  Wannier90_type::Wannier90type
  Wannier_Optional_Info::Wannier_OptionalInfo
  #Optinal::Dict{AbstractString,Any};
  spin_type::SPINtype;
  Optional::Dict{AbstractString,Any};
  Arg_Inputs() = new("",-1,-1,[(-1,-1)],[2,2,2],[2,2,2],
    convert(Array{Array{Int}}, [[]]),["all"],
    convert(Array{Array{Int}}, [[]]),["all"],
    nomask,nc_allH,"",0.0,"",
    Vector{k_point_Tuple}(0),NULLDFT,NULLWANNIER,Wannier_OptionalInfo(),
    colinear_type,
    Dict{AbstractString,Any}())
end

type kPath_band
  K_point_list::Array{k_point_Tuple}
  K_start_point_name::AbstractString
  K_end_point_name::AbstractString
end
type Arg_Inputs_Band
  bandplot::Bool
  K_point_groups::Array{kPath_band}
  Arg_Inputs_Band() = new(true,Array{kPath_band}(0))
end
function read_toml(toml_fname)
  toml_inputs =  TOML.parse(readstring(toml_fname));

end
function parse_int_list(num_str)
    a = split(num_str,",")
    intarray = [];
    for parsedV in a
        if (contains(parsedV,":"))
            tmp = split(parsedV,":");
            int_list = collect(parse(Int,tmp[1]):parse(Int,tmp[2]));
            for intv in int_list
              push!(intarray,intv);
            end
        else
            intv = parse(Int,parsedV)
            push!(intarray,intv);
        end
    end
    intarray = unique(intarray);
    return intarray;
end
function parse_Kpath(kPoint_toml,kPoint_step_num)
  K_point_groups = Array{kPath_band}(0);
  #kPoint_toml = toml_inputs["bandplot"]["kPath_list"]
  for (i,v) in enumerate(kPoint_toml)
    kPoint_step_num_interal = kPoint_step_num;
    k_point_start = (convert(Array{Float64,1}, v[1]))
    k_point_end = (convert(Array{Float64,1}, v[2]))
    assert(3==length(k_point_start) && 3==length(k_point_end))
    K_start_point_name = ""
    K_end_point_name = ""
    if (length(v)>=3 && (typeof(v[3]) <: AbstractString) )
      K_start_point_name = v[3][1];
      K_end_point_name = v[3][2];
      if (length(v)>=4 && (Int == typeof(v[4])))
        if (v[3]>0)
          kPoint_step_num_interal = v[4]
        end
      end
    end
    kPoint_steps = (k_point_end-k_point_start)/kPoint_step_num_interal
    K_point_list = Array{k_point_Tuple}(0);
    for steps_i in 0:kPoint_step_num_interal-1
      steps_f = steps_i*1.0;
      kPoint = (k_point_start[1] + steps_f*kPoint_steps[1]
      ,k_point_start[2]  + steps_f*kPoint_steps[2]
      ,k_point_start[3]  + steps_f*kPoint_steps[3] );
      push!(K_point_list,kPoint)
    end
    push!(K_point_groups,
      kPath_band(K_point_list,K_start_point_name,K_end_point_name))
  end
  return K_point_groups;
end
function parse_TOML(toml_file,input::Arg_Inputs)
  if (isfile(toml_file))
    toml_inputs = TOML.parse(readstring(input.TOMLinput))
    #println(toml_inputs)
    if (haskey(toml_inputs,"result_file"))
      input.result_file = toml_inputs["result_file"]
    end
    if (haskey(toml_inputs,"DFTtype"))
      DFT_type::AbstractString = toml_inputs["DFTtype"]
      if ( lowercase("OpenMX") == lowercase(DFT_type) )
        input.DFT_type = OpenMX
      elseif ( lowercase("Wannier90") == lowercase(DFT_type) )
        input.DFT_type = Wannier90
      end
    end
    if (haskey(toml_inputs,"spintype"))
      val = toml_inputs["spintype"];
      if ("para" == lowercase(val))
        input.spin_type = para_type
      elseif ("co_spin" == lowercase(val))
        input.spin_type = colinear_type
      elseif ("nc_spin" == lowercase(val))
        input.spin_type = non_colinear_type
      end
    end
    if (haskey(toml_inputs,"k_point_num"))
      k_point_num = convert(Array{Int,1},toml_inputs["k_point_num"]);
      assert(3 == length(k_point_num))
      input.k_point_num = k_point_num
    end
    if (haskey(toml_inputs,"q_point_num"))
      q_point_num = convert(Array{Int,1},toml_inputs["q_point_num"]);
      assert(3 == length(q_point_num))
      input.q_point_num = q_point_num
    end
    if (haskey(toml_inputs,"Wannier90type"))
      result_type::AbstractString = toml_inputs["Wannier90type"]
      if ( lowercase("openmx") == lowercase(result_type) )
        input.Wannier90_type = OpenMXWF;
      elseif ( lowercase("vasp") == lowercase(result_type))
        input.Wannier90_type = VASPWF;
      elseif ( lowercase("ecalj") == lowercase(result_type))
        input.Wannier90_type = EcalJWF;
      elseif ( lowercase("wien2k") == lowercase(result_type))

      end
    end
    if (haskey(toml_inputs,"atom12"))
      input.atom12_list = Vector{Tuple{Int64,Int64}}(0);
      for (k,v) in enumerate(toml_inputs["atom12"])
        push!(input.atom12_list,(v[1],v[2]));
      end
    end
    if (haskey(toml_inputs,"wannier_optional"))
      wannier_options = Wannier_OptionalInfo()

      atomnum = toml_inputs["wannier_optional"]["atomnum"]
      atompos = toml_inputs["wannier_optional"]["atompos"]
      atompos2 = convert(Array{Array{Float64},1},atompos)
      atompos = zeros(atomnum,3);
      for i in 1:atomnum
        atompos[i,:] = atompos2[i]
      end
      atoms_orbitals_list2 = toml_inputs["wannier_optional"]["atoms_orbitals_list"]
      atoms_orbitals_list = convert(Array{Array{Int},1},atoms_orbitals_list2)

      wannier_options.atomnum = atomnum;
      wannier_options.atompos = atompos;
      wannier_options.atoms_orbitals_list = atoms_orbitals_list;
      println(wannier_options)
      input.Wannier_Optional_Info = wannier_options
    end
    if (haskey(toml_inputs,"orbitals"))
      #println(toml_inputs["orbitals"])
      if (haskey(toml_inputs["orbitals"],"orbitalselection"))
        if ("on" == lowercase(toml_inputs["orbitals"]["orbitalselection"]))
          input.orbital_mask_option = unmask
          #input.orbital_mask_name = "masked"
          if (haskey(toml_inputs["orbitals"],"orbital_mask_option"))
            if ("unmask" == lowercase(toml_inputs["orbitals"]["orbital_mask_option"]))
              input.orbital_mask_option = unmask
            elseif ("mask" == lowercase(toml_inputs["orbitals"]["orbital_mask_option"]))
              input.orbital_mask_option = mask
            end
          end
          #if (haskey(toml_inputs["orbitals"],"orbital_mask_name"))
          #  input.orbital_mask_name = toml_inputs["orbitals"]["orbital_mask_name"]
          #end
          input.orbital_mask1_list = convert(Array{Array{Int}}, toml_inputs["orbitals"]["orbital_mask1_list"])
          input.orbital_mask1_names = split(toml_inputs["orbitals"]["orbital_mask1_names"],r"\[|\]|,",keep=false)
          input.orbital_mask2_list = convert(Array{Array{Int}},toml_inputs["orbitals"]["orbital_mask2_list"])
          input.orbital_mask2_names = split(toml_inputs["orbitals"]["orbital_mask2_names"],r"\[|\]|,",keep=false)
        end
      end
    end

    if (haskey(toml_inputs,"bandplot"))
      bandplot = true;
      kPoint_step_num = 20; # default steps
      if (haskey(toml_inputs["bandplot"],"bandplot"))
        bandplot = toml_inputs["bandplot"]["bandplot"];
      end
      if bandplot
        if (haskey(toml_inputs["bandplot"],"kPoint_step"))
          if (0 < toml_inputs["bandplot"]["kPoint_step"] )
            kPoint_step_num =  toml_inputs["bandplot"]["kPoint_step"]
          end
        end
        if (haskey(toml_inputs["bandplot"],"kPath_list"))
          K_point_groups::Array{kPath_band} =
            parse_Kpath(toml_inputs["bandplot"]["kPath_list"],kPoint_step_num);
          input.Optional["bandplot"] = K_point_groups
        end
      end
    end
  end
  input_checker(input);
  return input;
end
function parse_input(args,input::Arg_Inputs)

    #input::Arg_Inputs =  Arg_Inputs()

    s = ArgParseSettings("Example 2 for argparse.jl: " *  # description
                      "flags, options help, " *
                      "required arguments.")

    @add_arg_table s begin
        "--TOMLinput","-T"
        help = "input.toml file ex:) nio.toml "
        "--DFTtype","-D"
        help = " [OpenMX, Wannier90]"
        "--Wannier90type","-W"
        help =" Specify Wannier function type. Use with -D Wannier90 option. ex:) openmxWF [openmxWF, vaspWF, ecaljWF, wien2kWF] (openmxWF only supported for now)."
        "--atom12"
        help = "target atom1&2 ex:) 1_1,1_2,2_5"     # used by the help screen
        "--kpoint", "-k"
        #action = :store_true   # this makes it a flag
        help = "k_point ex:) 5_5_5"
        "--qpoint", "-q"
        #action = :store_true   # this makes it a flag
        help = "q_point ex:) 5_5_5"
        "--om1"
        help = "orbital_mask1_list ex) [[],[9,10,11,12,13],[9]] <= [all,d-only,z2]"
        "--om1names"
        help = "obital mask1 names ex) [all,d,z2]]"

        "--om2"
        help = "orbital_mask2_list ex) [[],[9,10],[11]] <= [all,eg,xy]"
        "--om2names"
        help = "obital mask2 names ex) [all,eg,xy]"

        "--soc"
        help = "spin orbitcoupling ex:) nc_allH=0 nc_realH_only=1 nc_imagH_only=2 "
        arg_type = Int
        "--ommode"
        help = "0 (default: no orbital control) 1(selected orbitals are unmasked) 2(selected orbitals are masked)"
        arg_type = Int
        #"--omname"
        #help = "obital masking name ex) d_d "
        "--hdftmpdir"
        help = "Specify hdf shared mem tmpdir dir (should be visible to all processes). Default: root_dir/jq/*.hdf5"
        "--chempdelta"
        help = "Chemical potential shift in eV(default 0.0 eV)"
        arg_type = Float64
        "--spintype","-s"
        help = "Spin type [para,cospin,ncspin] "
        "result_file"
        help = "Result file name ex:) nio.scfout (OpenMX scf), nio.HWR (OpenMX wannier)"
        required = false        # makes the argument mandatory
    end

    parsed_args = parse_args(args, s)
    #println("Parsed args:",parsed_args)
    for (key,val) in parsed_args
        #println("  $key  =>  $(repr(val))")
        #println("  $key  =>  $val")

        if (key == "atom12" && (Void != typeof(val)))
            atom_str_list = split(val,",")
            atom12_list = Vector{Tuple{Int64,Int64}}();
            for (atom_i, atom12_str) in enumerate(atom_str_list)
              if contains(atom12_str,"_")
                atom12 =   map(x->parse(Int64,x),split(atom12_str,"_"))
                push!(atom12_list,(atom12[1],atom12[2]));
              end
            end
            input.atom1 = atom12_list[1][1]
            input.atom2 = atom12_list[1][2]
            input.atom12_list = atom12_list
            #input.atom1 = atom12[1];
            #input.atom2 = atom12[2];
        end
        if (key == "result_file")
            if (typeof(val) <:AbstractString)
              #if (isfile(val) && ".scfout" == splitext(val)[2])
                  input.result_file  = val;
                  #println("scf file:$scf_name")
              #end
            end
        end
        if (key == "DFTtype" && typeof(val) <: AbstractString)
          DFT_type::AbstractString = val
          if ( lowercase("OpenMX") == lowercase(DFT_type) )
            input.DFT_type = OpenMX
          elseif ( lowercase("Wannier90") == lowercase(DFT_type) )
            input.DFT_type = Wannier90
          end
        end
        if (key == "Wannier90type" && typeof(val) <: AbstractString)
          result_type::AbstractString = val
          if ( lowercase("openmx") == lowercase(result_type) )
            input.Wannier90_type = OpenMXWF;
          elseif ( lowercase("vasp") == lowercase(result_type))
            input.Wannier90_type = VASPWF;
          elseif ( lowercase("ecalj") == lowercase(result_type))
            input.Wannier90_type = EcalJWF;
          elseif ( lowercase("wien2k") == lowercase(result_type))

          end
        end
        if (key == "TOMLinput" && typeof(val) <: AbstractString)
          # Check if file name endswith ".toml"
          input.TOMLinput = val;
        end
        if (key == "kpoint" && typeof(val) <: AbstractString)
            k_point_num_tmp =   map(x->parse(Int64,x),split(val,"_"))
            if(3 == length(k_point_num_tmp))
                input.k_point_num = k_point_num_tmp;
            else
                println("k point should be like 10_10_10")
            end
        end
        if (key == "qpoint" && typeof(val) <: AbstractString)
            q_point_num_tmp =   map(x->parse(Int64,x),split(val,"_"))
            if(3 == length(q_point_num_tmp))
                input.q_point_num = q_point_num_tmp;
            else
                println("q point should be like 10_10_10")
            end
        end
        if (key=="soc" && typeof(val) <: Int)
          # nc_allH=0 nc_realH_only=1 nc_imagH_only=2
          input.Hmode = val;
        end
        if (key =="om1" && typeof(val) <: AbstractString)
            #input.orbital_mask1 = parse_int_list(val)
            orbital_mask1_list = string("orbital_mask1_list = ",val);
            v = TOML.parse(orbital_mask1_list);
            input.orbital_mask1_list = convert(Array{Array{Int}},v["orbital_mask1_list"])
        end
        if (key =="om1names" && typeof(val) <: AbstractString)
            #println(val)
            input.orbital_mask1_names = split(val,r"\[|\]|,",keep=false)
        end
        if (key =="om2" && typeof(val) <: AbstractString)
            #input.orbital_mask2 = parse_int_list(val)
            orbital_mask2_list = string("orbital_mask2_list = ",val);
            v = TOML.parse(orbital_mask2_list);
            input.orbital_mask2_list =  convert(Array{Array{Int}},v["orbital_mask2_list"])
        end
        if (key =="om2names" && typeof(val) <: AbstractString)
            #println(val)
            input.orbital_mask2_names = split(val,r"\[|\]|,",keep=false)
        end
        if ("ommode"==key)
            if(1 == val)
               input.orbital_mask_option = unmask
            elseif(2 == val)
               input.orbital_mask_option = mask
            end
        end
        if ("omname"==key)
            #println(val)
            input.orbital_mask_name = val;
        end
        if ("hdftmpdir" ==key)
            if (nothing!=val)
                if isdir(val)
                    input.hdftmpdir = val;
                else
                    println(" $val is not directory")
                end
            end
        end
        if ("chempdelta" == key)
            if (nothing != val)
                input.ChemP_delta_ev = val;
            end
        end
        if ("spintype" == key && typeof(val) <: AbstractString)
          if ("para" == lowercase(val) )
            input.spin_type = para_type
          elseif ("co_spin" == lowercase(val))
            input.spin_type = colinear_type
          elseif ("nc_spin" == lowercase(val))
            input.spin_type = non_colinear_type
          else
            println(" Spintype must be one of [para, co_spin, nc_spin].")
            exit(1);
          end
        end

    end

    return input;
end
function input_checker(input::Arg_Inputs)
  exit_programe = false;
  if (""==input.result_file)
    # no result file
    println(" NO RESULT FILE SPECIFIED. TRY -h OPTION FOR HELP.")
    exit(1);
  end
  # Check Wannier90 properties
  if (Wannier90 == input.DFT_type)
    if (NULLWANNIER == input.Wannier90_type)
      println(" Set Wannier90type with -W option. TRY -h OPTION FOR HELP.")
      exit_programe = true;
    end
    if (OpenMXWF == input.Wannier90_type)
      if (input.Wannier_Optional_Info.atomnum <= 0)
        println(" For OpenMX wannier function, TOML file is need. Set atom position infomation in TOML. TRY -h OPTION FOR HELP.")
        exit_programe = true;
      end
    end
  end
  # Check orbital_selection properties
  if (nomask != input.orbital_mask_option)
    if (0 == length(input.orbital_mask1_list)|| 0 == length(input.orbital_mask2_list) )
      println(" orbital_mask1_list or orbital_mask2 is not set. ");
      exit_programe = true;
    end
    if (length(input.orbital_mask1_list) != length(input.orbital_mask1_names))
      println(" Length of orbital_mask1_list  orbital_mask1_names is not same. ");
      exit_programe = true;
    end
    if (length(input.orbital_mask2_list) != length(input.orbital_mask2_names))
      println(" Length of orbital_mask1_list  orbital_mask2_names is not same. ");
      exit_programe = true;
    end
  end

  if exit_programe
    println(DFTcommon.bar_string) # print ====...====
    println(" Exiting programe. Please set informations" )
    println(DFTcommon.bar_string) # print ====...====
    exit(1)
  end
end


end
