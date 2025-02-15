###############################################################################
# Hongkee Yoon Hongkeeyoon@kaist.ac.kr
# 2019.05
# https://kaist-elst.github.io/DFTforge.jl/
###############################################################################



################################################################################
# Input parser
################################################################################

using ArgParse

export nc_Hamiltonian_selection
export Arg_MPI
#@enum orbital_mask_enum nomask=0 unmask=1 mask=2
#@enum nc_Hamiltonian_selection nc_allH=0 nc_realH_only=1 nc_imagH_only=2

mutable struct Wannier_OptionalInfo
  atomnum::Int
  atompos::Array{Float64,2}
  atoms_orbitals_list::Array{Array{Int},1}
  Wannier_OptionalInfo() = new(0,zeros(0,3),Array{Array{Int},1}(undef,0));
end

mutable struct Arg_Inputs
  result_file::AbstractString
  result_file_dict::Dict{AbstractString,AbstractString}
  atom1::Int
  atom2::Int
  atom12_list::Vector{Tuple{Int64,Int64}}
  k_point_num::Array{Int,1}
  q_point_num::Array{Int,1}
  orbital_selection1_list::Array{Array{Int}}
  orbital_selection1_names::Array{AbstractString,1}
  orbital_selection2_list::Array{Array{Int}}
  orbital_selection2_names::Array{AbstractString,1}

  orbital_selection3_list::Array{Array{Int}}
  orbital_selection3_names::Array{AbstractString,1}
  orbital_selection4_list::Array{Array{Int}}
  orbital_selection4_names::Array{AbstractString,1}

  orbital_selection_option::orbital_selection_enum
  #orbital_selection_name


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
  Arg_Inputs() = new("",Dict{AbstractString,AbstractString}(),
  -1,-1,[(-1,-1)],[2,2,2],[2,2,2],
    convert(Array{Array{Int}}, [[]]),["all"],
    convert(Array{Array{Int}}, [[]]),["all"],
    convert(Array{Array{Int}}, [[]]),["all"], #orbital_selection3_list,orbital_selection3_names
    convert(Array{Array{Int}}, [[]]),["all"], #orbital_selection4_list,orbital_selection4_names
    nomask,nc_allH,"",0.0,"",
    Vector{k_point_Tuple}(undef,0),NULLDFT,NULLWANNIER,Wannier_OptionalInfo(),
    colinear_type,
    Dict{AbstractString,Any}())
end

struct kPath_band
  K_point_list::Array{k_point_Tuple}
  K_start_point_name::AbstractString
  K_end_point_name::AbstractString
end
struct Arg_Inputs_Band
  bandplot::Bool
  K_point_groups::Array{kPath_band}
  Arg_Inputs_Band() = new(true,Array{kPath_band}(undef,0))
end

struct Arg_DMFT_DFT
  start_iter::Int
  max_iter::Int
  dmft_executable::String

  openmx_input::String
  openmx_DM_executable::String
  openmx_thread_num::Int
  openmx_POST_executable::String

  charge_mixing::Float64
end

######### added by TJ
const eVToKelvin = 11604.5250061657
function SetOrDefault( parser, item_name, item_default_val )
    try
        return parser[item_name]
    catch
        return item_default_val
    end
end

struct Arg_DMFT_
  Calculation_mode::String
  mpi_prefix::String
  DMFT_loop_N::Int64
  KgridNum::Array{Int64,1}
  RgridNum::Array{Int64,1}
  iWgridCut::Float64
  # NImFreq::Int64
  Solver_green_Cut::Float64
  # beta::Float64
  Temperature::Float64
  Corr_atom_Ind::Array{Int64,1}
  Corr_orbital_Ind::Array{Array{Int64,1},1}
  Corr_atom_equiv::Array{Int64,1}
  Spin_type::Int64

  Mix_selfE::Float64
  init_bias::Float64
  smth_step::Int64
  cal_susc::Bool
  compute_EF::Bool

  DMFT_solver::String
  imp_dc_type::String
  imp_U::Array{Float64,1}
  imp_J::Array{Float64,1}
  imp_dc::Array{Float64,1}
  imp_int_type::Array{String,1}
  imp_int_parameterisation::Array{String,1}
  imp_block::Array{Array{Int64,2},1}
  imp_lev_shift::Array{Array{Float64,2},1}
  imp_Measure_time::Float64
  imp_Thermal_time::Float64

  basis_transform::String
  consider_degeneracy::Bool
  green_basis::String
  green_legendre_cutoff::Int64

  EDMFTF_MonteCar_step::Int64
  EDMFTF_warmup_step::Int64
  EDMFTF_GlobalFlip::Int64
  EDMFTF_tsample::Int64
  EDMFTF_nom::Int64
  EDMFTF_PChangeOrder::Float64

  SelfE_file::String
  chem_int::Float64
end


struct Arg_DFT_Jx
  Calculation_mode::String
  qmode::Int64
  Neighbors_cut::Int64
  KgridNum::Array{Int64,1}
  RgridNum::Array{Int64,1}
  iWgridCut::Float64
  Temperature::Float64
  Corr_atom_Ind::Array{Int64,1}
  Corr_orbital_Ind::Array{Array{Int64,1},1}
  Corr_atom_equiv::Array{Int64,1}
  beta::Float64
  NImFreq::Int64
end

struct Arg_GWEDMFT
  Dimension::Int64
  UnitCellVector::Array{Float64,2}
  OrbitalN::Int64
  AtomPosition::Array{Array{Float64,1},1}
  OrbitalInd::Array{Array{Int64,1},1}
  Corr_atom_Ind::Array{Int64,1}
  Corr_orbital_Ind::Array{Array{Int64,1},1}
  Corr_atom_equiv::Array{Int64,1}


  spintype::Int64

  OneBody::String
  Coulomb::String

  TaugridNum::Int64
  iWgridCut::Float64
  KgridNum::Array{Int64,1}
  RgridNum::Array{Int64,1}

  MaxLoopNum::Int64
  Temperature::Float64
  EleNum::Float64

  ThresholdChemP::Float64
  ThresholdGreenf::Float64

  GreenfOutAllk::Bool
  LocalSelfEOut::Bool
  LogOutOn::Bool

  GWMixingP::Float64
  GWMixingS::Float64

  DMFTMixingP::Float64
  DMFTMixingS::Float64

  Calculation_mode::String

  MushiftInDMFT_1stloop::Float64

  init_bias::Float64
  delta::Float64

  mpi_prefix::String
  compute_EF::Bool
  DMFT_solver::String
  Solver_green_Cut::Float64
  
  imp_half_filled::Array{Bool,1}
  imp_U::Array{Float64,1}
  imp_J::Array{Float64,1}
  imp_int_type::Array{String,1}
  imp_int_parameterisation::Array{String,1}
  imp_block::Array{Array{Int64,2},1}
  imp_lev_shift::Array{Array{Float64,2},1}
  imp_Measure_time::Float64
  imp_Thermal_time::Float64

  basis_transform::String
  consider_degeneracy::Bool
  green_basis::String
  green_legendre_cutoff::Int64

end

######### added by TJ





struct Arg_MPI
  mpirun_type::mpirun_type
  mpirun::String
  mpirun_num::Int
  mpirun_hostfile::String
  thread_num::Int
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
  K_point_groups = Array{kPath_band}(undef,0);
  #kPoint_toml = toml_inputs["bandplot"]["kPath_list"]
  for (i,v) in enumerate(kPoint_toml)
    kPoint_step_num_interal = kPoint_step_num;
    k_point_start = (convert(Array{Float64,1}, v[1]))
    k_point_end = (convert(Array{Float64,1}, v[2]))
    @assert(3==length(k_point_start) && 3==length(k_point_end))
    K_start_point_name = ""
    K_end_point_name = ""
    if (length(v)>=3 && (typeof(v[3][1]) <: AbstractString) )
      println(v)
      K_start_point_name = v[3][1];
      K_end_point_name = v[3][2];
      if (length(v)>=4 && (Int == typeof(v[4])))
        if (v[3]>0)
          kPoint_step_num_interal = v[4]
        end
      end
    end
    kPoint_steps = (k_point_end-k_point_start)/kPoint_step_num_interal
    K_point_list = Array{k_point_Tuple}(undef,0);
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
function detect_file(result_file,toml_realpath)
  found = false;
  if isfile(joinpath(pwd(),result_file))
    found = true;
    result_file = joinpath(pwd(),result_file)
  elseif isfile(result_file)
    found = true;
    result_file = result_file
  else
    toml_dir =  dirname(toml_realpath)
    if isfile(joinpath(toml_dir,result_file))
      result_file = joinpath(toml_dir,result_file)
      found = true;
    else

    end
  end
  return (found,result_file);
end
function parse_TOML(toml_file,input::Arg_Inputs)

  if (isfile(toml_file))
    toml_inputs = TOML.parse(read(input.TOMLinput,String))
    toml_realpath = normpath(input.TOMLinput); #realpath(input.TOMLinput); 
    toml_dir =  dirname(toml_realpath)
    println(" TOML file: ",toml_realpath)
    if (haskey(toml_inputs,"DFTtype"))
      DFT_type::AbstractString = toml_inputs["DFTtype"]
      if ( lowercase("OpenMX") == lowercase(DFT_type) )
        input.DFT_type = OpenMX
      elseif ( lowercase("EcalJ") == lowercase(DFT_type) )
        input.DFT_type = EcalJ
      elseif ( lowercase("Wannier90") == lowercase(DFT_type) || lowercase("Wannier") == lowercase(DFT_type) )
        input.DFT_type = Wannier90
      end
    end
    if (haskey(toml_inputs,"Wanniertype"))
      result_type::AbstractString = toml_inputs["Wanniertype"]
      if ( lowercase("openmx") == lowercase(result_type) )
        input.Wannier90_type = OpenMXWF;
      elseif ( lowercase("Wannier90") == lowercase(result_type))
        input.Wannier90_type = Wannier90WF;
      elseif ( lowercase("ecalj") == lowercase(result_type))
        input.Wannier90_type = EcalJWF;
      #elseif ( lowercase("wien2k") == lowercase(result_type))
      end
    end

    if (haskey(toml_inputs,"HamiltonianType"))
      if haskey(toml_inputs,"DFTtype")
        println("If `HamiltonianType` keyword is used, remove `DFTtype` and `Wanniertype` keyword.")
        exit(0);
      end
      Hamiltonian_type::AbstractString = toml_inputs["HamiltonianType"]
      if ( lowercase("OpenMX") == lowercase(Hamiltonian_type) )
        input.DFT_type = OpenMX
      elseif ( lowercase("OpenMXWannier") == lowercase(Hamiltonian_type))
        input.DFT_type = Wannier90
        input.Wannier90_type = OpenMXWF;
      elseif ( lowercase("Wannier90") == lowercase(Hamiltonian_type))
        input.DFT_type = Wannier90
        input.Wannier90_type = Wannier90WF;
      elseif ( lowercase("ecalj") == lowercase(Hamiltonian_type))
        input.DFT_type = EcalJ
      elseif ( lowercase("ecaljWannier") == lowercase(Hamiltonian_type))
        input.DFT_type = Wannier90
        input.Wannier90_type = EcalJWF;
      #elseif ( lowercase("wien2k") == lowercase(result_type))
      elseif ( lowercase("PlainwaveLobster") == lowercase(Hamiltonian_type))
        input.DFT_type = PlainwaveLobster
        #elseif ( lowercase("wien2k") == lowercase(result_type))
      elseif (lowercase("Wien2kDMFT") == lowercase(Hamiltonian_type))
        input.DFT_type = Wien2kDMFT
      else
        println("Hamiltonian_type ", Hamiltonian_type, " is not valid. Select from OpenMX, OpenMXWannier, Wannier90, ecalj, ecaljWannier, PlainwaveLobster, Wien2kDMFT")
      end
    end

    # Input files for OpenMX outputs

    if (OpenMX == input.DFT_type ||
      (Wannier90 == input.DFT_type &&  OpenMXWF ==  input.Wannier90_type) )
      ## OpenMX
      if (haskey(toml_inputs,"result_file"))
        result_file =  toml_inputs["result_file"]
        println(result_file) # TODO: remove it
        #if !isfile(result_file)
        toml_dir =  dirname(toml_realpath)
        if isfile(joinpath(pwd(),result_file))
          result_file = joinpath(pwd(),result_file)
        elseif isfile(joinpath(toml_dir,result_file))
          result_file = joinpath(toml_dir,result_file)
        end
        #end
        println(result_file)
        if isfile(result_file)
          input.result_file = result_file
        else
          print("File not found?")
        end
      end
    elseif ( EcalJ == input.DFT_type )
      if (haskey(toml_inputs,"result_file"))

          result_file_list_input =  toml_inputs["result_file"]
          if (length(result_file_list_input) > 1)
            result_file_1 = result_file_list_input[1]
            input.result_file = split(split(result_file_1,".")[1],"_")[1];
          end
          result_file_list = Array{AbstractString}(undef,0);
          for (k,v) in enumerate(result_file_list_input)
            (exits_check_dat,result_file_path) = detect_file(v , toml_realpath)
            println(result_file_path)
            if (exits_check_dat)
              result_file_dir  = dirname(result_file_path)
              result_file_path = joinpath(result_file_dir,v);
              push!(result_file_list,result_file_path);
            else
              println("File not found  ",v)
              @assert(false)
            end
          end
          result_file_dict = Dict(
          "result_file_up" => result_file_list[1],
          "result_file_down" => result_file_list[2])
          input.result_file_dict = result_file_dict;
          println(result_file_dict)
      end
    elseif ((Wannier90 == input.DFT_type &&  Wannier90WF ==  input.Wannier90_type))
      ## Wannier90
      #result_file =  toml_inputs["result_file"];
      if (haskey(toml_inputs,"result_file"))

        result_file_list_input =  toml_inputs["result_file"]
        if (length(result_file_list_input) > 1)
          result_file_1 = result_file_list_input[1]
          input.result_file = split(split(result_file_1,".")[1],"_")[1];
        end
        result_file_list = Array{AbstractString}(undef,0);
        for (k,v) in enumerate(result_file_list_input)
          (exits_check_dat,result_file_path) = detect_file(v * "_hr.dat", toml_realpath)
          println(result_file_path)
          (exits_check_win,result_file_path) = detect_file(v * ".win", toml_realpath)
          println(result_file_path)
          if (exits_check_dat && exits_check_win)
            result_file_dir  = dirname(result_file_path)
            result_file_path = joinpath(result_file_dir,v);
            push!(result_file_list,result_file_path);
          else
            println("File not found  ",v)
            @assert(false)
          end
        end
        
        if (haskey(toml_inputs,"spintype"))
          if (lowercase(toml_inputs["spintype"]) != "nc_spin")
            result_file_dict = Dict(
              "result_file_up" => result_file_list[1],
              "result_file_down" => result_file_list[2])
            input.result_file_dict = result_file_dict;
            println(result_file_dict)
            #result_file_dict = Dict
          else
            #input.result_file = result_file_list_input[1]
            result_file_dict = Dict(
              "result_file" => result_file_list_input[1])
              
            input.result_file_dict = result_file_dict;
            println(result_file_dict)
            #println(input.result_file)

          end

        end

      end
    elseif ((Wannier90 == input.DFT_type &&  EcalJWF ==  input.Wannier90_type))
      ## Wannier90
      #result_file =  toml_inputs["result_file"];
      if (haskey(toml_inputs,"result_file"))

        result_file_list_input =  toml_inputs["result_file"]
        if (length(result_file_list_input) > 1)
          result_file_1 = result_file_list_input[1]
          input.result_file = split(split(result_file_1,".")[1],"_")[1];
        end
        result_file_list = Array{AbstractString}(undef,0);
        for (k,v) in enumerate(result_file_list_input)
          (exits_check,result_file_path) = detect_file(v, toml_realpath)
          if (exits_check)
            push!(result_file_list,result_file_path);
          else
            println("File not found  ",v)
            @assert(false)
          end

        end
        result_file_dict = Dict(
        "result_file_up" => result_file_list[1],
        "result_file_down" => result_file_list[2])
        input.result_file_dict = result_file_dict;
        println(result_file_dict)
        #result_file_dict = Dict
      end
    elseif ((PlainwaveLobster == input.DFT_type))
      # PlainwaveLobster
      println("PlainwaveLobster ")
      if (haskey(toml_inputs,"result_file"))
        result_file_list_input =  toml_inputs["result_file"]
        if (length(result_file_list_input) > 1)
          #result_file_1 = result_file_list_input[1]
          #input.result_file = split(split(result_file_1,".")[1],"_")[1];
        end
        result_file_list = Array{AbstractString}(undef,0);
        for (k,v) in enumerate(result_file_list_input)
          (exits_check,result_file_path) = detect_file(v, toml_realpath)
          if (exits_check)
            push!(result_file_list,result_file_path);
          else
            println("File not found  ",v)
            @assert(false)
          end
        end
        result_file_dict = Dict(
        "RealSpaceHamiltonians" => result_file_list[1],
        "RealSpaceOverlaps" => result_file_list[2])
        input.result_file_dict = result_file_dict;
        println(result_file_dict)
      end

    elseif ((Wien2kDMFT == input.DFT_type))
      # Wien2K + DMFT_rutgers & json output patched
      if (haskey(toml_inputs,"result_file"))
        result_file =  toml_inputs["result_file"]
        println(result_file) # TODO: remove it
        #if !isfile(result_file)
        toml_dir =  dirname(toml_realpath)
        if isfile(joinpath(pwd(),result_file))
          result_file = joinpath(pwd(),result_file)
        elseif isfile(joinpath(toml_dir,result_file))
          result_file = joinpath(toml_dir,result_file)
        end
        #end
        println(result_file)
        if isfile(result_file)
          input.result_file = result_file
        else
          print("File not found?")
        end
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
      @assert(3 == length(k_point_num))
      input.k_point_num = k_point_num
    end
    if (haskey(toml_inputs,"q_point_num"))
      q_point_num = convert(Array{Int,1},toml_inputs["q_point_num"]);
      @assert(3 == length(q_point_num))
      input.q_point_num = q_point_num
    end

    if (haskey(toml_inputs,"atom12"))
      input.atom12_list = Vector{Tuple{Int64,Int64}}(undef,0);
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

    if (haskey(toml_inputs,"lobster_optional"))
      wannier_options = Wannier_OptionalInfo()

      atomnum = toml_inputs["lobster_optional"]["atomnum"]
      atompos = toml_inputs["lobster_optional"]["atompos"]
      atompos2 = convert(Array{Array{Float64},1},atompos)
      atompos = zeros(atomnum,3);
      for i in 1:atomnum
        atompos[i,:] = atompos2[i]
      end
      #atoms_orbitals_list2 = toml_inputs["lobster_optional"]["atoms_orbitals_list"]
      #atoms_orbitals_list = convert(Array{Array{Int},1},atoms_orbitals_list2)

      wannier_options.atomnum = atomnum;
      wannier_options.atompos = atompos;
      #wannier_options.atoms_orbitals_list = atoms_orbitals_list;
      println(wannier_options)
      input.Wannier_Optional_Info = wannier_options

      tv = zeros(3,3)
      a = toml_inputs["lobster_optional"]["tv"]
      println(a[1],typeof(a[1]))
      tv[1,:] = a[1]; tv[2,:] = a[2]; tv[3,:] = a[3]
      input.Optional["tv"] = tv;
      input.Optional["ChemP"] = toml_inputs["lobster_optional"]["ChemP"]
    end

    if (haskey(toml_inputs,"bandselection"))
      input.Optional["band_selection"] =  toml_inputs["bandselection"]["band_selection"]::Bool
      input.Optional["band_selection_boundary"] =  toml_inputs["bandselection"]["band_selection_boundary"]
    end
    if (haskey(toml_inputs,"orbitals"))
      #println(toml_inputs["orbitals"])
      if (haskey(toml_inputs["orbitals"],"orbitalselection"))
        if (toml_inputs["orbitals"]["orbitalselection"])
          input.orbital_selection_option = unmask
          #input.orbital_selection_name = "masked"


          if (haskey(toml_inputs["orbitals"],"orbital_mask_option"))

            if ("unmask" == lowercase(toml_inputs["orbitals"]["orbital_mask_option"]))
              input.orbital_selection_option = unmask
            elseif ("mask" == lowercase(toml_inputs["orbitals"]["orbital_mask_option"]))
              input.orbital_selection_option = mask
            end


          end

          ## The option keyword 'orbital_mask_option' -> renamed to 'orbital_selection_option' 201905
          if (haskey(toml_inputs["orbitals"],"orbital_mask1_list"))
            input.orbital_selection1_list = convert(Array{Array{Int}}, toml_inputs["orbitals"]["orbital_mask1_list"])
            input.orbital_selection1_names = split(toml_inputs["orbitals"]["orbital_mask1_names"],r"\[|\]|,",keepempty=false)
            input.orbital_selection2_list = convert(Array{Array{Int}},toml_inputs["orbitals"]["orbital_mask2_list"])
            input.orbital_selection2_names = split(toml_inputs["orbitals"]["orbital_mask2_names"],r"\[|\]|,",keepempty=false)
          end

          if (haskey(toml_inputs["orbitals"],"orbital_selection1_list"))
            input.orbital_selection1_list = convert(Array{Array{Int}}, toml_inputs["orbitals"]["orbital_selection1_list"])
            input.orbital_selection1_names = split(toml_inputs["orbitals"]["orbital_selection1_names"],r"\[|\]|,",keepempty=false)
            input.orbital_selection2_list = convert(Array{Array{Int}},toml_inputs["orbitals"]["orbital_selection2_list"])
            input.orbital_selection2_names = split(toml_inputs["orbitals"]["orbital_selection2_names"],r"\[|\]|,",keepempty=false)
          end
          ##

          if (haskey(toml_inputs["orbitals"],"orbital_selection_option"))
            if ("unmask" == lowercase(toml_inputs["orbitals"]["orbital_selection_option"]))
              input.orbital_selection_option = unmask
            elseif ("mask" == lowercase(toml_inputs["orbitals"]["orbital_selection_option"]))
              input.orbital_selection_option = mask
            end

          end
          #if (haskey(toml_inputs["orbitals"],"orbital_mask_name"))
          #  input.orbital_mask_name = toml_inputs["orbitals"]["orbital_mask_name"]
          #end

          if (haskey(toml_inputs["orbitals"],"orbital_selection3_list"))
            input.orbital_selection3_list = convert(Array{Array{Int}}, toml_inputs["orbitals"]["orbital_selection3_list"])
            input.orbital_selection3_names = split(toml_inputs["orbitals"]["orbital_selection3_names"],r"\[|\]|,",keepempty=false)
          end
          if (haskey(toml_inputs["orbitals"],"orbital_selection4_list"))
            input.orbital_selection4_list = convert(Array{Array{Int}}, toml_inputs["orbitals"]["orbital_selection4_list"])
            input.orbital_selection4_names = split(toml_inputs["orbitals"]["orbital_selection4_names"],r"\[|\]|,",keepempty=false)
          end

        end
      end
    end

    if (haskey(toml_inputs,"energywindow"))
        if (haskey(toml_inputs["energywindow"],"energywindow"))
            energywindow_on = toml_inputs["energywindow"]["energywindow"]
            if energywindow_on
                #input.Optional["energywindow"]  = Dict{String,Array{Array{Float64,1},1}}()
                input.Optional["energywindow"]  = Dict{String,Any}()
                if (!haskey(toml_inputs["energywindow"],"energywindow_name"))
                    println(" In [energywindow] the energywindow_name is missing "); @assert(false);
                end
                input.Optional["energywindow"]["energywindow_name"] = toml_inputs["energywindow"]["energywindow_name"]


                input.Optional["energywindow"]["energywindow_all_list"]  = Array{Array{Float64,1},1}()
                input.Optional["energywindow"]["energywindow_1_list"]  = Array{Array{Float64,1},1}()
                input.Optional["energywindow"]["energywindow_2_list"]  = Array{Array{Float64,1},1}()
                input.Optional["energywindow"]["energywindow_3_list"]  = Array{Array{Float64,1},1}()
                input.Optional["energywindow"]["energywindow_4_list"]  = Array{Array{Float64,1},1}()

                if (haskey(toml_inputs["energywindow"],"energywindow_all_list"))
                    energywindow_all_list = convert(Array{Array{Float64,1},1},  toml_inputs["energywindow"]["energywindow_all_list"])
                    for eRange in energywindow_all_list
                        if (2 != length(eRange))
                            println(" energywindow required lowerBound & upperBound  ", eRange);
                            @assert(false);
                        end
                    end
                    input.Optional["energywindow"]["energywindow_all_list"] = energywindow_all_list;
                end

                if (haskey(toml_inputs["energywindow"],"energywindow_1_list"))
                    energywindow_1_list = convert(Array{Array{Float64,1},1},  toml_inputs["energywindow"]["energywindow_1_list"])
                    for eRange in energywindow_1_list
                        if (2 != length(eRange))
                            println(" energywindow required lowerBound & upperBound  ", eRange);
                            @assert(false);
                        end
                    end
                    input.Optional["energywindow"]["energywindow_1_list"] = energywindow_1_list;
                end

                if (haskey(toml_inputs["energywindow"],"energywindow_2_list"))
                    energywindow_2_list = convert(Array{Array{Float64,1},1},  toml_inputs["energywindow"]["energywindow_2_list"])
                    for eRange in energywindow_2_list
                        if (2 != length(eRange))
                            println(" energywindow required lowerBound & upperBound  ", eRange);
                            @assert(false);
                        end
                    end
                    input.Optional["energywindow"]["energywindow_2_list"] = energywindow_2_list;
                end

                if (haskey(toml_inputs["energywindow"],"energywindow_3_list"))
                    energywindow_3_list = convert(Array{Array{Float64,1},1},  toml_inputs["energywindow"]["energywindow_3_list"])
                    for eRange in energywindow_3_list
                        if (2 != length(eRange))
                            println(" energywindow required lowerBound & upperBound  ", eRange);
                            @assert(false);
                        end
                    end
                    input.Optional["energywindow"]["energywindow_3_list"] = energywindow_3_list;
                end

                if (haskey(toml_inputs["energywindow"],"energywindow_4_list"))
                    energywindow_4_list = convert(Array{Array{Float64,1},1},  toml_inputs["energywindow"]["energywindow_4_list"])
                    for eRange in energywindow_4_list
                        if (2 != length(eRange))
                            println(" energywindow required lowerBound & upperBound  ", eRange);
                            @assert(false);
                        end
                    end
                    input.Optional["energywindow"]["energywindow_4_list"] = energywindow_4_list;
                end

            end
        end

    end

    if (haskey(toml_inputs,"orbital_reassign") || haskey(toml_inputs,"jeff_projection_on"))

      orbital_rot_on = false
      orbital_rot_rules = Dict{Int,orbital_rot_type}()
      if (haskey(toml_inputs["orbital_reassign"],"orbital_rot_on"))
        orbital_rot_on = toml_inputs["orbital_reassign"]["orbital_rot_on"]
      end

      jeff_projection_on = false
      jeff_orbitals = Nothing
      if (haskey(toml_inputs["orbital_reassign"],"jeff_projection_on"))
        jeff_projection_on = toml_inputs["orbital_reassign"]["jeff_projection_on"];
        if jeff_projection_on
          orbital_rot_on = true
        end
      end

      if orbital_rot_on || jeff_projection_on
        if (haskey(toml_inputs["orbital_reassign"],"d_orbital_rot"))
          #orbital_rot_d_dict = Dict{Int,orbital_rot_d_type}()
          for (k,v) in enumerate(toml_inputs["orbital_reassign"]["d_orbital_rot"])
            # orbital_rot_d_type
            atom1 = convert(Int, v[1][1])
            Z_vect = convert(Vector{Float64},v[2]);
            X_vect = convert(Vector{Float64},v[3]);
            d_orbital_list = convert(Array{Array{Int}},v[4+0:4+4]);
            orbital_rot_rules[atom1] =
              orbital_rot_type(atom1,Z_vect,X_vect,d_orbital_list,5);
          end
            #TODO: migrate into Arg_Inputs
          #input.Optional["orbital_rot_d_dict"] = orbital_rot_d_dict;
        end

        orbital_downfold_on = false;
        if (haskey(toml_inputs["orbital_reassign"],"orbital_downfold_on"))
          orbital_downfold_on = toml_inputs["orbital_reassign"]["orbital_downfold_on"]
        end
        orbital_downfold_rules = Dict{Int,orbital_downfold_type}()
        keep_unmerged_atoms = true;
        keep_unmerged_orbitals = true;
        if orbital_downfold_on
          if (haskey(toml_inputs["orbital_reassign"],"keep_unmerged_atoms"))
            keep_unmerged_atoms = toml_inputs["orbital_reassign"]["keep_unmerged_atoms"]
          end
          if (haskey(toml_inputs["orbital_reassign"],"keep_unmerged_orbitals"))
            keep_unmerged_orbitals = toml_inputs["orbital_reassign"]["keep_unmerged_orbitals"]
          end
          if (haskey(toml_inputs["orbital_reassign"],"orbital_merge_downfolding"))


            for (k,v) in enumerate(toml_inputs["orbital_reassign"]["orbital_merge_downfolding"])
              # orbital_rot_d_type
              atom1 = convert(Int, v[1][1])
              orbital_list = convert(Array{Array{Int}},v[2:end]);
              @assert(issorted( map(v->v[1],orbital_list) ));
              orbital_downfold_rules[atom1] =
                orbital_downfold_type(atom1,orbital_list);
            end
              #TODO: migrate into Arg_Inputs

          end

        
        end

        position_merge_on = false
        postion_merge_rule_list = Nothing
        if (haskey(toml_inputs["orbital_reassign"],"position_merge"))
          position_merge_on = toml_inputs["orbital_reassign"]["position_merge"]
          if position_merge_on
            postion_merge_rule_list = toml_inputs["orbital_reassign"]["position_merge_rule"]
            Array{Union{Array{Float64,1}, Array{Int64,1}},1} == typeof(postion_merge_rule_list)

            # [new atom index],[prev atom indices],[new x,new y,new z]
            for rule_i in 1:length(postion_merge_rule_list)
              postion_merge_rule = postion_merge_rule_list[rule_i]
              @assert(typeof(postion_merge_rule[1]) <: Int)
              @assert(typeof(postion_merge_rule[2]) <: Array{Int,1})
              @assert(typeof(postion_merge_rule[3]) <: Array{Float64,1})
              @assert(3 == length(postion_merge_rule[3]))
              typeof(postion_merge_rule[1])
            end
            sort!(postion_merge_rule_list, by = x -> x[1])

          end
        end

        orbital_custom_transform_on = false;
        custom_transfrom_rules = Nothing
        if (haskey(toml_inputs["orbital_reassign"],"custom_transform"))
          orbital_custom_transform_on = toml_inputs["orbital_reassign"]["custom_transform"] 

          if orbital_custom_transform_on
            custom_transfrom_rules = Array{custom_transfrom_type,1}(undef,0)
            custom_rule = toml_inputs["orbital_reassign"]["custom_transfrom_rule"] # TODO: MY 

            for rule_i in 1:length(custom_rule)
              print(custom_rule[rule_i][1]," ",custom_rule[rule_i][2])
              custom_rule_filename =  custom_rule[rule_i][1][1]
              custom_rule_varible =  custom_rule[rule_i][1][2]
          
              custom_rule_remap = custom_rule[rule_i][2]
              custom_U = Nothing
              print(custom_rule_filename)
              if isfile(abspath(custom_rule_filename))
                  custom_U = include(abspath(custom_rule_filename))
              elseif isfile(joinpath(toml_dir, custom_rule_filename))
                  custom_U = include(joinpath(toml_dir, custom_rule_filename))
              else
                  println(" custom_rule_filename ",custom_rule_filename ," dose not exist.")
                  @assert(false)
              end
              
          
              if Nothing != custom_U
                  U_matrix = custom_U[custom_rule_varible]
                  orbital_remap_list = Array{remap_type,1}(undef,0)
                  display(custom_U[custom_rule_varible])
                  # print(custom_U["testf"](4,2))
                  display(custom_rule_remap)
                  for (remap_i,reamp) in enumerate(custom_rule_remap)
                    # (atom1, orbital)
                    push!(orbital_remap_list, remap_type(reamp[1], reamp[2])) # , reamp[3]))
                  end
                  push!(custom_transfrom_rules, 
                    custom_transfrom_type(copy(U_matrix),orbital_remap_list))
              end
            end
        
          end
          # input.Optional["custom_orbital_tranformfule"] = custom_rule
        end

        basisTransform = basisTransform_rule_type(orbital_rot_on,orbital_rot_rules,orbital_downfold_on,
          orbital_downfold_rules,keep_unmerged_atoms,keep_unmerged_orbitals, 
          postion_merge_rule_list, custom_transfrom_rules);
          
        #=
        jeff_orbitals = Nothing
        if jeff_projection_on
          jeff_orbitals = toml_inputs["orbital_reassign"]["jeff_orbitals"];
          # Check   
          jeff_projection_num = length(jeff_orbitals)
          for Jeff_i in 1:jeff_projection_num
            jeff_orbital_cluster = jeff_orbitals[Jeff_i]
              @assert(4 == length(jeff_orbital_cluster))

              for i in 1:4
                  @assert(4 == length(jeff_orbital_cluster[i]))  # [[atom],[dxy],[dxz],[dyz]]
                  @assert(1 == length(jeff_orbital_cluster[i][1]))  # [[atom],[dxy],[dxz],[dyz]]
              end
          end
          #input.Optional["jeff_orbitals"] = jeff_orbitals;
        end

        #basisTransform = basisTransform_rule_type(orbital_rot_on,orbital_rot_rules,orbital_merge_on,
        #  orbital_merge_rules,keep_unmerged_atoms,keep_unmerged_orbitals);
        basisTransform = basisTransform_rule_type(orbital_rot_on, orbital_rot_rules, orbital_downfold_on,
        orbital_downfold_rules, keep_unmerged_atoms, keep_unmerged_orbitals, jeff_orbitals); # TODO: Jeff
        =#
        input.Optional["basisTransform"] = basisTransform;
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
    if (haskey(toml_inputs,"soc"))
      val = toml_inputs["soc"]
      if (0==val)
        input.Hmode = nc_allH;
      elseif 1==val
        input.Hmode = nc_realH_only;
      elseif 2==val
        input.Hmode = nc_imagH_only;
      else
        println(" soc must be on of number ex:)[nc_allH=0 nc_realH_only=1 nc_imagH_only=2].")
        exit(1);
      end
    end




    ############# added by TJ #############
    if ( haskey(toml_inputs,"DMFT_Jx") || haskey(toml_inputs,"KaiEDJ") )
       
       header_dmft_jx = "DMFT_Jx"
       if( haskey(toml_inputs,"KaiEDJ") )
          header_dmft_jx = "KaiEDJ"
       end
       toml_input_header_dmft_jx = toml_inputs[header_dmft_jx]

       Calculation_mode = "DMFT"
       
       if (haskey(toml_input_header_dmft_jx,"Calculation_mode"))
          Calculation_mode = toml_input_header_dmft_jx["Calculation_mode"]
       end

       if (Calculation_mode == "DMFT")
          KgridNum = [5,5,5]
          RgridNum = [5,5,5]
          # iWgridCut = 30
          Solver_green_Cut = 30
          imp_Measure_time = 20;
          imp_Thermal_time = 1;
          mpi_prefix="mpirun -np 4"
          Mix_selfE = 0.5;
          DMFT_loop_N = 100;
          init_bias =0.0;
          basis_transform = "false"
          consider_degeneracy = false
          smth_step = 10;
          cal_susc = false
          green_basis ="matsubara"
          green_legendre_cutoff = 50
          imp_dc_type ="nominal"
          DMFT_solver ="ComCTQMC"
          SelfE_file = ""
          chem_int = 0.0


          EDMFTF_MonteCar_step = 1000000
          EDMFTF_warmup_step = 100000
          EDMFTF_GlobalFlip = 100000
          EDMFTF_tsample = 100
          EDMFTF_nom = 100
          EDMFTF_PChangeOrder = 0.9

          compute_EF = true

          if (haskey(toml_input_header_dmft_jx,"compute_EF"))
            compute_EF = toml_input_header_dmft_jx["compute_EF"]
          end

          if (haskey(toml_input_header_dmft_jx,"EDMFTF_MonteCar_step"))
            EDMFTF_MonteCar_step = toml_input_header_dmft_jx["EDMFTF_MonteCar_step"]
          end

          if (haskey(toml_input_header_dmft_jx,"EDMFTF_warmup_step"))
            EDMFTF_warmup_step = toml_input_header_dmft_jx["EDMFTF_warmup_step"]
          end

          if (haskey(toml_input_header_dmft_jx,"EDMFTF_GlobalFlip"))
            EDMFTF_GlobalFlip = toml_input_header_dmft_jx["EDMFTF_GlobalFlip"]
          end

          if (haskey(toml_input_header_dmft_jx,"EDMFTF_tsample"))
            EDMFTF_tsample = toml_input_header_dmft_jx["EDMFTF_tsample"]
          end

          if (haskey(toml_input_header_dmft_jx,"EDMFTF_nom"))
            EDMFTF_nom = toml_input_header_dmft_jx["EDMFTF_nom"]
          end

          if (haskey(toml_input_header_dmft_jx,"EDMFTF_PChangeOrder"))
            EDMFTF_PChangeOrder = toml_input_header_dmft_jx["EDMFTF_PChangeOrder"]
          end





          if (haskey(toml_input_header_dmft_jx,"DMFT_solver"))
            DMFT_solver = toml_input_header_dmft_jx["DMFT_solver"]
          end

          if (haskey(toml_input_header_dmft_jx,"DMFT_loop_N"))
            DMFT_loop_N = toml_input_header_dmft_jx["DMFT_loop_N"]
          end

          if (haskey(toml_input_header_dmft_jx,"basis_transform"))
            basis_transform = toml_input_header_dmft_jx["basis_transform"]
          end

          if (haskey(toml_input_header_dmft_jx,"consider_degeneracy"))
            consider_degeneracy = toml_input_header_dmft_jx["consider_degeneracy"]
          end


          if (haskey(toml_input_header_dmft_jx,"KgridNum"))
            KgridNum = toml_input_header_dmft_jx["KgridNum"]
          end
          if (haskey(toml_input_header_dmft_jx,"RgridNum"))
            RgridNum = toml_input_header_dmft_jx["RgridNum"]
          end

          beta          = SetOrDefault( toml_input_header_dmft_jx, "beta", NaN )
          if (isnan(beta))
            Temperature   = SetOrDefault( toml_input_header_dmft_jx, "Temperature", 90.66035161070313 )
          else
            Temperature   = 1. / beta * eVToKelvin
          end

          NImFreq       = SetOrDefault( toml_input_header_dmft_jx, "NImFreq", NaN )
          if (isnan(NImFreq))
            iWgridCut = SetOrDefault( toml_input_header_dmft_jx, "iWgridCut", NaN )
          else
            iWgridCut = (2*NImFreq)*pi/beta
          end
          # if (haskey(toml_input_header_dmft_jx,"iWgridCut"))
          #   iWgridCut = toml_input_header_dmft_jx["iWgridCut"]
          # end

          if (haskey(toml_input_header_dmft_jx,"Solver_green_Cut"))
            Solver_green_Cut = toml_input_header_dmft_jx["Solver_green_Cut"]
          end

          if (haskey(toml_input_header_dmft_jx,"mpi_prefix"))
            mpi_prefix = toml_input_header_dmft_jx["mpi_prefix"]
          end
          if (haskey(toml_input_header_dmft_jx,"smth_step"))
            smth_step = toml_input_header_dmft_jx["smth_step"]
          end
          if (haskey(toml_input_header_dmft_jx,"cal_susc"))
            cal_susc = toml_input_header_dmft_jx["cal_susc"]
          end

          if (haskey(toml_input_header_dmft_jx,"green_basis"))
            green_basis = toml_input_header_dmft_jx["green_basis"]
          end
          if (haskey(toml_input_header_dmft_jx,"green_legendre_cutoff"))
            green_legendre_cutoff = toml_input_header_dmft_jx["green_legendre_cutoff"]
          end




          if (haskey(toml_input_header_dmft_jx,"imp_Measure_time"))
            imp_Measure_time = toml_input_header_dmft_jx["imp_Measure_time"]
          end
          if (haskey(toml_input_header_dmft_jx,"imp_Thermal_time"))
            imp_Thermal_time = toml_input_header_dmft_jx["imp_Thermal_time"]
          end

          if (haskey(toml_input_header_dmft_jx,"Mix_selfE"))
            Mix_selfE = toml_input_header_dmft_jx["Mix_selfE"]
          end

          if (haskey(toml_input_header_dmft_jx,"init_bias"))
            init_bias = toml_input_header_dmft_jx["init_bias"]
          end


          Corr_atom_Ind = toml_input_header_dmft_jx["Corr_atom_Ind"]
          Corr_orbital_Ind = toml_input_header_dmft_jx["Corr_orbital_Ind"]
          Corr_atom_equiv = toml_input_header_dmft_jx["Corr_atom_equiv"]

          Spin_type = toml_input_header_dmft_jx["DMFT_Spin_type"]

          if (haskey(toml_input_header_dmft_jx,"imp_dc_type"))
             imp_dc_type = toml_input_header_dmft_jx["imp_dc_type"]
          end
          imp_U=[]
          imp_J=[]
          imp_dc=[]
          imp_int_type =[]
          imp_int_parameterisation =[]
          imp_block=[]
          imp_lev_shift=[]
          for i in unique(abs.(Corr_atom_equiv))
             str=string("imp",i,"_U")
             push!(imp_U, toml_input_header_dmft_jx[str])

             str=string("imp",i,"_J")
             push!(imp_J,toml_input_header_dmft_jx[str])

             str=string("imp",i,"_dc")
             push!(imp_dc,toml_input_header_dmft_jx[str])

             str=string("imp",i,"_int_type")
             push!(imp_int_type,toml_input_header_dmft_jx[str])

             str=string("imp",i,"_int_parameterisation")
             push!(imp_int_parameterisation,toml_input_header_dmft_jx[str])

             str=string("imp",i,"_block")
             Ind_2d_block=zeros(Int64, size( toml_input_header_dmft_jx[str] )[1], size( toml_input_header_dmft_jx[str][1] )[1] )
             for j = 1:size( toml_input_header_dmft_jx[str] )[1]
               Ind_2d_block[j,:] = toml_input_header_dmft_jx[str][j]
             end
             push!(imp_block,Ind_2d_block)


             str2=string("imp",i,"_lev_shift")
             lev_shift_2d_block=zeros(Float64, size( toml_input_header_dmft_jx[str] )[1], size( toml_input_header_dmft_jx[str][1] )[1] )
             if (haskey(toml_input_header_dmft_jx,str2))
                 for j = 1:size( toml_input_header_dmft_jx[str2] )[1]
                    lev_shift_2d_block[j,:] = toml_input_header_dmft_jx[str2][j]
                 end
             end
             push!(imp_lev_shift,lev_shift_2d_block)

          end
          imp_U=convert(Array{Float64,1},imp_U)
          imp_J=convert(Array{Float64,1},imp_J)
          imp_dc=convert(Array{Float64,1},imp_dc)
          imp_int_type=convert(Array{String,1},imp_int_type)
          imp_int_parameterisation=convert(Array{String,1},imp_int_parameterisation)
          imp_block = convert(Array{Array{Int64,2},1}, imp_block)
          imp_lev_shift = convert(Array{Array{Float64,2},1}, imp_lev_shift)
          input.Optional[header_dmft_jx] = Arg_DMFT_(Calculation_mode,mpi_prefix,DMFT_loop_N, KgridNum, RgridNum, iWgridCut, Solver_green_Cut, Temperature, Corr_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, Spin_type, Mix_selfE,init_bias,smth_step,cal_susc,compute_EF,DMFT_solver,imp_dc_type, imp_U,imp_J,imp_dc,imp_int_type,imp_int_parameterisation,imp_block, imp_lev_shift, imp_Measure_time, imp_Thermal_time, basis_transform,consider_degeneracy, green_basis, green_legendre_cutoff, EDMFTF_MonteCar_step, EDMFTF_warmup_step, EDMFTF_GlobalFlip, EDMFTF_tsample, EDMFTF_nom, EDMFTF_PChangeOrder,SelfE_file,chem_int)


       elseif(Calculation_mode == "Jx-DMFT" || Calculation_mode == "Jx-DMFT-private" || Calculation_mode=="DMFT+MFT")
          KgridNum = [5,5,5]
          RgridNum = [5,5,5]
          # iWgridCut = 30
          Solver_green_Cut = 30
          imp_Measure_time = 20;
          imp_Thermal_time = 1;
          mpi_prefix="mpirun -np 4"
          Mix_selfE = 0.5;
          DMFT_loop_N = 100;
          init_bias =0.0;
          consider_degeneracy = false
          basis_transform = "false"
          smth_step = 10;
          cal_susc = false
          green_basis ="matsubara"
          green_legendre_cutoff = 50
          imp_dc_type ="nominal"
          DMFT_solver = "ComCTQMC"
          SelfE_file = ""
          chem_int = 0.0

          EDMFTF_MonteCar_step = 1000000
          EDMFTF_warmup_step = 100000
          EDMFTF_GlobalFlip = 100000
          EDMFTF_tsample = 100
          EDMFTF_nom = 100
          EDMFTF_PChangeOrder = 0.9

          compute_EF = true

          if (haskey(toml_input_header_dmft_jx,"compute_EF"))
            compute_EF = toml_input_header_dmft_jx["compute_EF"]
          end


          if (haskey(toml_input_header_dmft_jx,"EDMFTF_MonteCar_step"))
            EDMFTF_MonteCar_step = toml_input_header_dmft_jx["EDMFTF_MonteCar_step"]
          end

          if (haskey(toml_input_header_dmft_jx,"EDMFTF_warmup_step"))
            EDMFTF_warmup_step = toml_input_header_dmft_jx["EDMFTF_warmup_step"]
          end

          if (haskey(toml_input_header_dmft_jx,"EDMFTF_GlobalFlip"))
            EDMFTF_GlobalFlip = toml_input_header_dmft_jx["EDMFTF_GlobalFlip"]
          end

          if (haskey(toml_input_header_dmft_jx,"EDMFTF_tsample"))
            EDMFTF_tsample = toml_input_header_dmft_jx["EDMFTF_tsample"]
          end

          if (haskey(toml_input_header_dmft_jx,"EDMFTF_nom"))
            EDMFTF_nom = toml_input_header_dmft_jx["EDMFTF_nom"]
          end

          if (haskey(toml_input_header_dmft_jx,"EDMFTF_PChangeOrder"))
            EDMFTF_PChangeOrder = toml_input_header_dmft_jx["EDMFTF_PChangeOrder"]
          end


          if (haskey(toml_input_header_dmft_jx,"DMFT_solver"))
            DMFT_solver = toml_input_header_dmft_jx["DMFT_solver"]
          end

          if (haskey(toml_input_header_dmft_jx,"SelfE_file"))
            SelfE_file = toml_input_header_dmft_jx["SelfE_file"]
          end

          if (haskey(toml_input_header_dmft_jx,"chem_int"))
            chem_int = toml_input_header_dmft_jx["chem_int"]
          end


          if (haskey(toml_input_header_dmft_jx,"DMFT_loop_N"))
            DMFT_loop_N = toml_input_header_dmft_jx["DMFT_loop_N"]
          end

          if (haskey(toml_input_header_dmft_jx,"basis_transform"))
            basis_transform = toml_input_header_dmft_jx["basis_transform"]
          end

          if (haskey(toml_input_header_dmft_jx,"consider_degeneracy"))
            consider_degeneracy = toml_input_header_dmft_jx["consider_degeneracy"]
          end



          if (haskey(toml_input_header_dmft_jx,"KgridNum"))
            KgridNum = toml_input_header_dmft_jx["KgridNum"]
          end
          if (haskey(toml_input_header_dmft_jx,"RgridNum"))
            RgridNum = toml_input_header_dmft_jx["RgridNum"]
          end

          beta          = SetOrDefault( toml_input_header_dmft_jx, "beta", NaN )
          if (isnan(beta))
            Temperature   = SetOrDefault( toml_input_header_dmft_jx, "Temperature", 90.66035161070313 )
          else
            Temperature   = 1. / beta * eVToKelvin
          end

          NImFreq       = SetOrDefault( toml_input_header_dmft_jx, "NImFreq", NaN )
          if (isnan(NImFreq))
            iWgridCut = SetOrDefault( toml_input_header_dmft_jx, "iWgridCut", NaN )
          else
            iWgridCut = (2*NImFreq)*pi/beta
          end

          if (haskey(toml_input_header_dmft_jx,"Solver_green_Cut"))
            Solver_green_Cut = toml_input_header_dmft_jx["Solver_green_Cut"]
          end

          if (haskey(toml_input_header_dmft_jx,"mpi_prefix"))
            mpi_prefix = toml_input_header_dmft_jx["mpi_prefix"]
          end
          if (haskey(toml_input_header_dmft_jx,"smth_step"))
            smth_step = toml_input_header_dmft_jx["smth_step"]
          end
          if (haskey(toml_input_header_dmft_jx,"cal_susc"))
            cal_susc = toml_input_header_dmft_jx["cal_susc"]
          end

          if (haskey(toml_input_header_dmft_jx,"green_basis"))
            green_basis = toml_input_header_dmft_jx["green_basis"]
          end

          if (haskey(toml_input_header_dmft_jx,"green_legendre_cutoff"))
            green_legendre_cutoff = toml_input_header_dmft_jx["green_legendre_cutoff"]
          end

          if (haskey(toml_input_header_dmft_jx,"imp_Measure_time"))
            imp_Measure_time = toml_input_header_dmft_jx["imp_Measure_time"]
          end
          if (haskey(toml_input_header_dmft_jx,"imp_Thermal_time"))
            imp_Thermal_time = toml_input_header_dmft_jx["imp_Thermal_time"]
          end

          if (haskey(toml_input_header_dmft_jx,"Mix_selfE"))
            Mix_selfE = toml_input_header_dmft_jx["Mix_selfE"]
          end

          if (haskey(toml_input_header_dmft_jx,"init_bias"))
            init_bias = toml_input_header_dmft_jx["init_bias"]
          end
          Corr_atom_Ind = toml_input_header_dmft_jx["Corr_atom_Ind"]
          Corr_orbital_Ind = toml_input_header_dmft_jx["Corr_orbital_Ind"]
          Corr_atom_equiv = toml_input_header_dmft_jx["Corr_atom_equiv"]

          Spin_type = toml_input_header_dmft_jx["DMFT_Spin_type"]

          if (haskey(toml_input_header_dmft_jx,"imp_dc_type"))
             imp_dc_type = toml_input_header_dmft_jx["imp_dc_type"]
          end

          imp_U=[]
          imp_J=[]
          imp_dc=[]
          imp_int_type =[]
          imp_int_parameterisation =[]
          imp_block=[]
          imp_lev_shift=[]
          for i in unique(abs.(Corr_atom_equiv))
             str=string("imp",i,"_U")
             push!(imp_U, toml_input_header_dmft_jx[str])

             str=string("imp",i,"_J")
             push!(imp_J,toml_input_header_dmft_jx[str])

             str=string("imp",i,"_dc")
             push!(imp_dc,toml_input_header_dmft_jx[str])

             str=string("imp",i,"_int_type")
             push!(imp_int_type,toml_input_header_dmft_jx[str])

             str=string("imp",i,"_int_parameterisation")
             push!(imp_int_parameterisation,toml_input_header_dmft_jx[str])

             str=string("imp",i,"_block")
             Ind_2d_block=zeros(Int64, size( toml_input_header_dmft_jx[str] )[1], size( toml_input_header_dmft_jx[str][1] )[1] )
             for j = 1:size( toml_input_header_dmft_jx[str] )[1]
               Ind_2d_block[j,:] = toml_input_header_dmft_jx[str][j]
             end
             push!(imp_block,Ind_2d_block)
             str2=string("imp",i,"_lev_shift")
             lev_shift_2d_block=zeros(Float64, size( toml_input_header_dmft_jx[str] )[1], size( toml_input_header_dmft_jx[str][1] )[1] )
             if (haskey(toml_input_header_dmft_jx,str2))
                 for j = 1:size( toml_input_header_dmft_jx[str2] )[1]
                    lev_shift_2d_block[j,:] = toml_input_header_dmft_jx[str2][j]
                 end
             end
             push!(imp_lev_shift,lev_shift_2d_block)


          end
          imp_U=convert(Array{Float64,1},imp_U)
          imp_J=convert(Array{Float64,1},imp_J)
          imp_dc=convert(Array{Float64,1},imp_dc)
          imp_int_type=convert(Array{String,1},imp_int_type)
          imp_int_parameterisation=convert(Array{String,1},imp_int_parameterisation)
          imp_block = convert(Array{Array{Int64,2},1}, imp_block)
          imp_lev_shift = convert(Array{Array{Float64,2},1}, imp_lev_shift)
          input.Optional[header_dmft_jx] = Arg_DMFT_(Calculation_mode,mpi_prefix,DMFT_loop_N, KgridNum, RgridNum, iWgridCut, Solver_green_Cut, Temperature, Corr_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, Spin_type, Mix_selfE,init_bias,smth_step,cal_susc,compute_EF,DMFT_solver,imp_dc_type, imp_U,imp_J,imp_dc,imp_int_type,imp_int_parameterisation,imp_block, imp_lev_shift, imp_Measure_time, imp_Thermal_time, basis_transform,consider_degeneracy, green_basis, green_legendre_cutoff, EDMFTF_MonteCar_step, EDMFTF_warmup_step, EDMFTF_GlobalFlip, EDMFTF_tsample, EDMFTF_nom, EDMFTF_PChangeOrder,SelfE_file,chem_int)



       elseif (Calculation_mode == "Jx0" || Calculation_mode == "MFT" || Calculation_mode == "Spinwave"|| Calculation_mode == "Magnon")

          KgridNum = [5,5,5]
          RgridNum = [5,5,5]
          # iWgridCut = 30
          # Temperature = 500
          qmode = 1
          Neighbors_cut = 300

          if (haskey(toml_input_header_dmft_jx,"qmode"))
            qmode = toml_input_header_dmft_jx["qmode"]
          end

          if (haskey(toml_input_header_dmft_jx,"Neighbors_cut"))
            Neighbors_cut = toml_input_header_dmft_jx["Neighbors_cut"]
          end

          beta          = SetOrDefault( toml_input_header_dmft_jx, "beta", 0 )
          if (iszero(beta))
            Temperature   = SetOrDefault( toml_input_header_dmft_jx, "Temperature", 90.66035161070313 )
          else
            Temperature   = 1. / beta * eVToKelvin
          end

          NImFreq       = SetOrDefault( toml_input_header_dmft_jx, "NImFreq", 0 )
          if (iszero(NImFreq))
            iWgridCut = SetOrDefault( toml_input_header_dmft_jx, "iWgridCut", 30 )
          else
            iWgridCut = (2*NImFreq)*pi/beta
          end

          if (haskey(toml_input_header_dmft_jx,"KgridNum"))
            KgridNum = toml_input_header_dmft_jx["KgridNum"]
          end
          if (haskey(toml_input_header_dmft_jx,"RgridNum"))
            RgridNum = toml_input_header_dmft_jx["RgridNum"]
          end

          if (haskey(toml_input_header_dmft_jx,"iWgridCut"))
            iWgridCut = toml_input_header_dmft_jx["iWgridCut"]
          end


          Corr_atom_Ind = toml_input_header_dmft_jx["Corr_atom_Ind"]
          Corr_orbital_Ind = toml_input_header_dmft_jx["Corr_orbital_Ind"]
          Corr_atom_equiv = toml_input_header_dmft_jx["Corr_atom_equiv"]


          input.Optional[header_dmft_jx] = Arg_DFT_Jx(Calculation_mode, qmode, Neighbors_cut,KgridNum, RgridNum,iWgridCut,Temperature, Corr_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, beta, NImFreq)
       end

    end


    if (haskey(toml_inputs,"GW_EDMFT"))
       TaugridNum = 2000
       iWgridCut = 20.0
       KgridNum = [6,6,6]
       RgridNum = [6,6,6]

       OneBody = "OneBody.in"
       Coulomb = "Coulomb.in"

       MaxLoopNum = 20

       ThresholdChemP = 0.000000001
       ThresholdGreenf = 0.0001

       GreenfOutAllk = false
       LocalSelfEOut = true
       LogOutOn = false

       GWMixingP = 0.5
       GWMixingS = 0.5

       DMFTMixingP = 0.5
       DMFTMixingS = 0.5

       MushiftInDMFT_1stloop = 0.0

       init_bias = 0.0
       delta = 0.0
       Solver_green_Cut = iWgridCut

       mpi_prefix="mpirun -np 4"
       compute_EF = true
       DMFT_solver ="ComCTQMC"
       imp_Measure_time = 20;
       imp_Thermal_time = 1;

       basis_transform = "false"
       consider_degeneracy = false
       green_basis ="matsubara"
       green_legendre_cutoff = 50  



       Dimension = toml_inputs["GW_EDMFT"]["Dimension"]
       dum_cell_vec = zeros(3,3)
       for i=1:3
          for ii =1:3
             dum_cell_vec[i,ii] = toml_inputs["GW_EDMFT"]["UnitCellVector"][i][ii]
          end
       end
       UnitCellVector = dum_cell_vec

       OrbitalN = toml_inputs["GW_EDMFT"]["OrbitalN"]

       AtomPosition = toml_inputs["GW_EDMFT"]["AtomPosition"]
       OrbitalInd = toml_inputs["GW_EDMFT"]["OrbitalInd"]
       Corr_atom_Ind = toml_inputs["GW_EDMFT"]["Corr_atom_Ind"]
       Corr_orbital_Ind = toml_inputs["GW_EDMFT"]["Corr_orbital_Ind"]
       Corr_atom_equiv = toml_inputs["GW_EDMFT"]["Corr_atom_equiv"]
       spintype = toml_inputs["GW_EDMFT"]["spintype"]

       if (haskey(toml_inputs["GW_EDMFT"],"OneBody"))
          OneBody = toml_inputs["GW_EDMFT"]["OneBody"]
       end
       if (haskey(toml_inputs["GW_EDMFT"],"Coulomb"))
          Coulomb = toml_inputs["GW_EDMFT"]["Coulomb"]
       end

       if (haskey(toml_inputs["GW_EDMFT"],"TaugridNum"))
          TaugridNum = toml_inputs["GW_EDMFT"]["TaugridNum"]
       end
       if (haskey(toml_inputs["GW_EDMFT"],"iWgridCut"))
          iWgridCut = toml_inputs["GW_EDMFT"]["iWgridCut"]
       end
       if (haskey(toml_inputs["GW_EDMFT"],"KgridNum"))
          KgridNum = toml_inputs["GW_EDMFT"]["KgridNum"]
       end
       if (haskey(toml_inputs["GW_EDMFT"],"RgridNum"))
          RgridNum = toml_inputs["GW_EDMFT"]["RgridNum"]
       end

       if (haskey(toml_inputs["GW_EDMFT"],"MaxLoopNum"))
          MaxLoopNum = toml_inputs["GW_EDMFT"]["MaxLoopNum"]
       end

       Temperature = toml_inputs["GW_EDMFT"]["Temperature"]
       EleNum = toml_inputs["GW_EDMFT"]["EleNum"]

       if (haskey(toml_inputs["GW_EDMFT"],"ThresholdChemP"))
          ThresholdChemP = toml_inputs["GW_EDMFT"]["ThresholdChemP"]
       end
       if (haskey(toml_inputs["GW_EDMFT"],"ThresholdGreenf"))
          ThresholdGreenf = toml_inputs["GW_EDMFT"]["ThresholdGreenf"]
       end


       if (haskey(toml_inputs["GW_EDMFT"],"GreenfOutAllk"))
          GreenfOutAllk = toml_inputs["GW_EDMFT"]["GreenfOutAllk"]
       end
       if (haskey(toml_inputs["GW_EDMFT"],"LocalSelfEOut"))
          LocalSelfEOut = toml_inputs["GW_EDMFT"]["LocalSelfEOut"]
       end
       if (haskey(toml_inputs["GW_EDMFT"],"LogOutOn"))
          LogOutOn = toml_inputs["GW_EDMFT"]["LogOutOn"]
       end


       if (haskey(toml_inputs["GW_EDMFT"],"GWMixingP"))
          GWMixingP = toml_inputs["GW_EDMFT"]["GWMixingP"]
       end
       if (haskey(toml_inputs["GW_EDMFT"],"GWMixingS"))
          GWMixingS = toml_inputs["GW_EDMFT"]["GWMixingS"]
       end
       if (haskey(toml_inputs["GW_EDMFT"],"DMFTMixingP"))
          DMFTMixingP = toml_inputs["GW_EDMFT"]["DMFTMixingP"]
       end
       if (haskey(toml_inputs["GW_EDMFT"],"DMFTMixingS"))
          DMFTMixingS = toml_inputs["GW_EDMFT"]["DMFTMixingS"]
       end

       Calculation_mode = toml_inputs["GW_EDMFT"]["Calculation_mode"]

       if (haskey(toml_inputs["GW_EDMFT"],"MushiftInDMFT_1stloop"))
          MushiftInDMFT_1stloop = toml_inputs["GW_EDMFT"]["MushiftInDMFT_1stloop"]
       end




       if (haskey(toml_inputs["GW_EDMFT"],"init_bias"))
          init_bias = toml_inputs["GW_EDMFT"]["init_bias"]
       end
       if (haskey(toml_inputs["GW_EDMFT"],"delta"))
          delta = toml_inputs["GW_EDMFT"]["delta"]
       end
       if (haskey(toml_inputs["GW_EDMFT"],"Solver_green_Cut"))
          Solver_green_Cut = toml_inputs["GW_EDMFT"]["Solver_green_Cut"]
       end


       if (haskey(toml_inputs["GW_EDMFT"],"mpi_prefix"))
          mpi_prefix = toml_inputs["GW_EDMFT"]["mpi_prefix"]
       end

       if (haskey(toml_inputs["GW_EDMFT"],"compute_EF"))
          compute_EF = toml_inputs["GW_EDMFT"]["compute_EF"]
       end

       if (haskey(toml_inputs["GW_EDMFT"],"DMFT_solver"))
          DMFT_solver = toml_inputs["GW_EDMFT"]["DMFT_solver"]
       end

       if (haskey(toml_inputs["GW_EDMFT"],"imp_Measure_time"))
          imp_Measure_time = toml_inputs["GW_EDMFT"]["imp_Measure_time"]
       end
       if (haskey(toml_inputs["GW_EDMFT"],"imp_Thermal_time"))
          imp_Thermal_time = toml_inputs["GW_EDMFT"]["imp_Thermal_time"]
       end

       if (haskey(toml_inputs["GW_EDMFT"],"basis_transform"))
          basis_transform = toml_inputs["GW_EDMFT"]["basis_transform"]
       end

       if (haskey(toml_inputs["GW_EDMFT"],"consider_degeneracy"))
          consider_degeneracy = toml_inputs["GW_EDMFT"]["consider_degeneracy"]
       end

       if (haskey(toml_inputs["GW_EDMFT"],"green_basis"))
          green_basis = toml_inputs["GW_EDMFT"]["green_basis"]
       end

       if (haskey(toml_inputs["GW_EDMFT"],"green_legendre_cutoff"))
          green_legendre_cutoff = toml_inputs["GW_EDMFT"]["green_legendre_cutoff"]
       end


       imp_half_filled = []
       imp_U=[]
       imp_J=[]
       imp_int_type =[]
       imp_int_parameterisation =[]
       imp_block=[]
       imp_lev_shift=[]
       for i in unique(abs.(Corr_atom_equiv))
          push!(imp_half_filled, false)
          
          str=string("imp",i,"_U")
          push!(imp_U, toml_inputs["GW_EDMFT"][str])

          str=string("imp",i,"_J")
          push!(imp_J,toml_inputs["GW_EDMFT"][str])


          str=string("imp",i,"_int_type")
          push!(imp_int_type,toml_inputs["GW_EDMFT"][str])

          str=string("imp",i,"_int_parameterisation")
          push!(imp_int_parameterisation,toml_inputs["GW_EDMFT"][str])

          str=string("imp",i,"_block")
          Ind_2d_block=zeros(Int64, size( toml_inputs["GW_EDMFT"][str] )[1], size( toml_inputs["GW_EDMFT"][str][1] )[1] )
          for j = 1:size( toml_inputs["GW_EDMFT"][str] )[1]
             Ind_2d_block[j,:] = toml_inputs["GW_EDMFT"][str][j]
          end
          push!(imp_block,Ind_2d_block)
          str2=string("imp",i,"_lev_shift")
          lev_shift_2d_block=zeros(Float64, size( toml_inputs["GW_EDMFT"][str] )[1], size( toml_inputs["GW_EDMFT"][str][1] )[1] )
          if (haskey(toml_inputs["GW_EDMFT"],str2))
              for j = 1:size( toml_inputs["GW_EDMFT"][str2] )[1]
                 lev_shift_2d_block[j,:] = toml_inputs["GW_EDMFT"][str2][j]
              end
          end
          push!(imp_lev_shift,lev_shift_2d_block)
       end

       if (haskey(toml_inputs["GW_EDMFT"],"imp_half_filled"))
          imp_half_filled = toml_inputs["GW_EDMFT"]["imp_half_filled"]
       end

       imp_half_filled = convert(Array{Bool,1},imp_half_filled)
       imp_U=convert(Array{Float64,1},imp_U)
       imp_J=convert(Array{Float64,1},imp_J)
       imp_int_type=convert(Array{String,1},imp_int_type)
       imp_int_parameterisation=convert(Array{String,1},imp_int_parameterisation)
       imp_block = convert(Array{Array{Int64,2},1}, imp_block)
       imp_lev_shift = convert(Array{Array{Float64,2},1}, imp_lev_shift)


       input.Optional["GW_EDMFT"] = Arg_GWEDMFT(Dimension, UnitCellVector, OrbitalN, AtomPosition, OrbitalInd, Corr_atom_Ind, Corr_orbital_Ind, Corr_atom_equiv, spintype, OneBody, Coulomb, TaugridNum, iWgridCut, KgridNum, RgridNum, MaxLoopNum, Temperature, EleNum, ThresholdChemP, ThresholdGreenf, GreenfOutAllk, LocalSelfEOut, LogOutOn, GWMixingP, GWMixingS, DMFTMixingP, DMFTMixingS, Calculation_mode, MushiftInDMFT_1stloop, init_bias, delta, mpi_prefix, compute_EF, DMFT_solver,  Solver_green_Cut, imp_half_filled, imp_U, imp_J, imp_int_type,  imp_int_parameterisation,  imp_block,  imp_lev_shift,  imp_Measure_time,  imp_Thermal_time, basis_transform, consider_degeneracy,  green_basis,  green_legendre_cutoff)
    end


    ############# added by TJ #############














    #Input parser for the DMFT + DFT
    if (haskey(toml_inputs,"MPI"))


      mpirun_type_str = lowercase(toml_inputs["MPI"]["mpirun_type"])
      mpirun_type = mvapich;
      if ("mvapich" == mpirun_type_str)
        mpirun_type = mvapich;
      elseif ("openmpi" == mpirun_type_str)
        mpirun_type = openmpi;
      end

      mpirun = "mpirun"
      if (haskey(toml_inputs["MPI"], "mpirun"))
        mpirun =  toml_inputs["MPI"]["mpirun"]
      end

      mpirun_num = 4
      input.Optional["MPI"] = Arg_MPI(mpirun_type, mpirun, mpirun_num, "", 1)
      if haskey(input.Optional,"MPI_HOSTFILE")
        MPI_hostfile_input = input.Optional["MPI_HOSTFILE"];
        MPI_hostfile =  abspath(MPI_hostfile_input)
        # input.Optional["MPI"].mpirun_hostfile = input.Optional["MPI_HOSTFILE"];
        if (isfile(MPI_hostfile))
           input.Optional["MPI"].mpirun_hostfile = MPI_hostfile
        end
      end
      if haskey(toml_inputs["MPI"], "thread_num")
        input.Optional["MPI"].thread_num = toml_inputs["MPI"]["thread_num"];
      end

    end

    
    if (haskey(toml_inputs,"DMFT"))

        start_iter = 0;
        if (haskey(toml_inputs["DMFT"],"start_iter"))
          start_iter = toml_inputs["DMFT"]["start_iter"];
          @assert(start_iter>=0);
        end
        max_iter = convert(Int64,toml_inputs["DMFT"]["max_iter"]);
        dmft_executable = toml_inputs["DMFT"]["dmft_executable"]::String
        charge_mixing = 0.3;
        if (haskey(toml_inputs["DMFT"],"charge_mixing"))
          charge_mixing   = toml_inputs["DMFT"]["charge_mixing"]::Float64
        end
        # DMFT.DFT
        openmx_input = toml_inputs["DMFT"]["DFT"]["openmx_input"]
        openmx_input_full = joinpath(toml_dir,"dft",openmx_input)

        if (isfile(openmx_input_full))

        else
          println("openmx_input dft/*.dat dose not exists", openmx_input_full)
          @assert(false)
        end

        openmx_DM_executable = toml_inputs["DMFT"]["DFT"]["openmx_DM_executable"]
        openmx_POST_executable = toml_inputs["DMFT"]["DFT"]["openmx_POST_executable"]
        openmx_thread_num = toml_inputs["DMFT"]["DFT"]["openmx_DM_executable_thread_num"]


        input.Optional["DMFT_DFT"] = Arg_DMFT_DFT(start_iter, max_iter, dmft_executable,
         openmx_input_full, openmx_DM_executable, openmx_thread_num, openmx_POST_executable, charge_mixing);
    end

  end
  input_checker(input);
  return input;
end


function parse_input(args,input::Arg_Inputs)

    #input::Arg_Inputs =  Arg_Inputs()

    s = ArgParseSettings("Input help page: " *  # description
                      "flags, options help, " *
                      "required arguments.")

    @add_arg_table! s begin
        "--TOMLinput","-T"
        help = "input.toml file ex:) nio.toml "
        "--MPI_HOSTFILE"
        help = "MPI_HOSTFILE: file containing the host names "
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
        help = "orbital_selection1_list ex) [[],[9,10,11,12,13],[9]] <= [all,d-only,z2]"
        "--om1names"
        help = "obital mask1 names ex) [all,d,z2]"

        "--om2"
        help = "orbital_selection2_list ex) [[],[9,10],[11]] <= [all,eg,xy]"
        "--om2names"
        help = "obital mask2 names ex) [all,eg,xy]"

        "--om3"
        help = "orbital_selection1_list ex) [[],[9,10,11,12,13] <= [all,d-only]"
        "--om3names"
        help = "obital mask1 names ex) [all,d]"

        "--om4"
        help = "orbital_selection2_list ex) [[],[9,10,11,12,13] <= [all,d-only]"
        "--om4names"
        help = "obital mask2 names ex) [all,d]"

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
        help = "Spin type [para, co_spin, nc_spin] "
        "result_file"
        help = "Result file name ex:) nio.scfout (OpenMX scf), nio.HWR (OpenMX wannier)"
        required = false        # makes the argument mandatory
    end

    parsed_args = parse_args(args, s)
    #println("Parsed args:",parsed_args)
    for (key,val) in parsed_args
        #println("  $key  =>  $(repr(val))")
        #println("  $key  =>  $val")

        if (key == "atom12" && (Nothing != typeof(val)))
            atom_str_list = split(val,",")
            atom12_list = Vector{Tuple{Int64,Int64}}();
            for (atom_i, atom12_str) in enumerate(atom_str_list)
              if occursin("_",atom12_str)
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
              result_file = val;
              if !isfile(result_file)
                if isfile(joinpath(pwd(),result_file))
                  result_file = joinpath(pwd(),result_file)
                end
              end
              if isfile(result_file)
                input.result_file  = val;
              end
            end
        end
        if (key == "DFTtype" && typeof(val) <: AbstractString)
          DFT_type::AbstractString = val
          if ( lowercase("OpenMX") == lowercase(DFT_type) )
            input.DFT_type = OpenMX
        elseif ( lowercase("Wannier90") == lowercase(DFT_type) || lowercase("Wannier") == lowercase(DFT_type) )
            input.DFT_type = Wannier90
          end
        end
        if (key == "Wannier90type" && typeof(val) <: AbstractString)
          result_type::AbstractString = val
          if ( lowercase("openmx") == lowercase(result_type) )
            input.Wannier90_type = OpenMXWF;
          elseif ( lowercase("Wannier90") == lowercase(result_type))
            input.Wannier90_type = Wannier90WF;
          elseif ( lowercase("ecalj") == lowercase(result_type))
            input.Wannier90_type = EcalJWF;
          elseif ( lowercase("wien2k") == lowercase(result_type))

          end
        end
        if (key == "TOMLinput" && typeof(val) <: AbstractString)
          # Check if file name endswith ".toml"
          input.TOMLinput = val;
        end

        if (key == "MPI_HOSTFILE" && typeof(val) <: AbstractString)
          if (isfile(val))
            input.Optional["MPI_HOSTFILE"] = val;
          else
            println("MPI_HOSTFILE ",val, " dose not exits!")
          end
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
          if (0==val)
            input.Hmode = nc_allH;
          elseif 1==val
            input.Hmode = nc_realH_only;
          elseif 2==val
            input.Hmode = nc_imagH_only;
          else
            println(" soc must be on of number ex:)[nc_allH=0 nc_realH_only=1 nc_imagH_only=2].")
            exit(1);
          end


        end
        if (key =="om1" && typeof(val) <: AbstractString)
            #input.orbital_mask1 = parse_int_list(val)
            orbital_selection1_list = string("orbital_selection1_list = ",val);
            v = TOML.parse(orbital_selection1_list);
            input.orbital_selection1_list = convert(Array{Array{Int}},v["orbital_selection1_list"])
        end
        if (key =="om1names" && typeof(val) <: AbstractString)
            #println(val)
            input.orbital_selection1_names = split(val,r"\[|\]|,",keepempty=false)
        end
        if (key =="om2" && typeof(val) <: AbstractString)
            #input.orbital_selection2 = parse_int_list(val)
            orbital_selection2_list = string("orbital_selection2_list = ",val);
            v = TOML.parse(orbital_selection2_list);
            input.orbital_selection2_list =  convert(Array{Array{Int}},v["orbital_selection2_list"])
        end
        if (key =="om2names" && typeof(val) <: AbstractString)
            #println(val)
            input.orbital_selection2_names = split(val,r"\[|\]|,",keepempty=false)
        end

        if (key =="om3" && typeof(val) <: AbstractString)
            #input.orbital_selection1 = parse_int_list(val)
            orbital_selection3_list = string("orbital_selection3_list = ",val);
            v = TOML.parse(orbital_selection3_list);
            input.orbital_selection3_list = convert(Array{Array{Int}},v["orbital_selection3_list"])
        end
        if (key =="om3names" && typeof(val) <: AbstractString)
            #println(val)
            input.orbital_selection3_names = split(val,r"\[|\]|,",keepempty=false)
        end
        if (key =="om4" && typeof(val) <: AbstractString)
            #input.orbital_selection2 = parse_int_list(val)
            orbital_selection4_list = string("orbital_selection4_list = ",val);
            v = TOML.parse(orbital_selection4_list);
            input.orbital_selection4_list =  convert(Array{Array{Int}},v["orbital_selection4_list"])
        end
        if (key =="om4names" && typeof(val) <: AbstractString)
            #println(val)
            input.orbital_selection4_names = split(val,r"\[|\]|,",keepempty=false)
        end

        if ("ommode"==key)
            if(1 == val)
               input.orbital_selection_option = unmask
            elseif(2 == val)
               input.orbital_selection_option = mask
            end
        end
        if ("omname"==key)
            #println(val)
            input.orbital_selection_name = val;
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
  if (""==input.result_file && 0 == length(input.result_file_dict) )
    # no result file
    println(" NO RESULT FILE SPECIFIED. TRY -h OPTION FOR HELP.")
    exit_programe = true;
    #exit(1);
  end
  if !isfile(input.result_file) && (0 == length(input.result_file_dict)) && (0 == length(input.result_file_dict))
    println(" result file is not found ", input.result_file );
  end
  for (k,v) in input.result_file_dict
    if (!isfile(v))
      println(" result file is not found ",k,"\t",v );
    end
  end
  if (DFTcommon.OpenMX == input.DFT_type)
  elseif (DFTcommon.EcalJ == input.DFT_type)
  elseif (DFTcommon.Wannier90 == input.DFT_type)
  elseif (DFTcommon.PlainwaveLobster == input.DFT_type)
  elseif (DFTcommon.Wien2kDMFT == input.DFT_type )
  else
    println(" Set DFTresult type with -D option. TRY -h OPTION FOR HELP.")
    exit_programe = true;
  end

  if (0==length(input.atom12_list))
    println(" WARNNING: Set atom12 list with -atom12 option. TRY -h OPTION FOR HELP.")
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
  if (nomask != input.orbital_selection_option)
    if (0 == length(input.orbital_selection1_list)|| 0 == length(input.orbital_selection2_list) )
      println(" orbital_selection1_list or orbital_selection2 is not set. ");
      exit_programe = true;
    end
    if (length(input.orbital_selection1_list) != length(input.orbital_selection1_names))
      println(" Length of orbital_selection1_list  orbital_selection1_names is not same. ");
      exit_programe = true;
    end
    if (length(input.orbital_selection2_list) != length(input.orbital_selection2_names))
      println(" Length of orbital_selection1_list  orbital_mask2_names is not same. ");
      exit_programe = true;
    end
  end

  if exit_programe
    sleep(2);
    println(DFTcommon.bar_string) # print ====...====
    println(" Exiting programe. Please set informations" )
    println(DFTcommon.bar_string) # print ====...====
    exit(1)
  end
end
#end
