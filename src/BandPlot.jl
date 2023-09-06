import Plots
using Distributed
using ProgressMeter
import DFTforge
using DFTforge.DFTrefinery
using DFTforge.DFTcommon

X_VERSION = VersionNumber("0.4.0-dev+20170515");

if 1 == myid()
  println(" X_VERSION: ",X_VERSION)
end

@everywhere import DFTforge
@everywhere using DFTforge.DFTrefinery
@everywhere using DFTforge.DFTcommon
##############################################################################
## 1. Read INPUT
## 1.1 Set Default values
## 1.2 Read input from argument & TOML file
## 1.3 Set values from intput (arg_input)
## 1.4 Set caluations type and ouput folder
##############################################################################

hdftmpdir = ""
## 1.1 Set Default values
#orbital_selection1 = Array{Int64,1}();
#orbital_selection2 = Array{Int64,1}();
orbital_selection1_list = Array{Array{Int}}(undef,0);
orbital_selection1_names = Array{AbstractString}(undef,0);
orbital_selection2_list = Array{Array{Int}}(undef,0);
orbital_selection2_names = Array{AbstractString}(undef,0);

orbital_selection3_list = Array{Array{Int}}(undef,0);
orbital_selection3_names = Array{AbstractString}(undef,0);
orbital_selection4_list = Array{Array{Int}}(undef,0);
orbital_selection4_names = Array{AbstractString}(undef,0);

orbital_selection_option = DFTcommon.nomask;
orbital_selection_on = false;

band_selection_on = false;
band_selection_upper = 0;
band_selection_lower = 0;

k_point_num = [3,3,3]
q_point_num = [3,3,3]
ChemP_delta_ev = 0.0
DFT_type = DFTcommon.OpenMX

## 1.2 Read input from argument & TOML file
arg_input = DFTcommon.Arg_Inputs();
arg_input = parse_input(ARGS,arg_input)
#arg_input.TOMLinput = "vasp_lobster/ver_3.2.0/test.toml" # Debug

arg_input = parse_TOML(arg_input.TOMLinput,arg_input)
# let argument override
arg_input = parse_input(ARGS,arg_input)

## 1.3 Set values from intput (arg_input)
DFT_type = arg_input.DFT_type
Wannier90_type = arg_input.Wannier90_type
spin_type = arg_input.spin_type

result_file = arg_input.result_file
result_file_dict = arg_input.result_file_dict;

ChemP_delta_ev = arg_input.ChemP_delta_ev
 # k,q point num
k_point_num = arg_input.k_point_num
q_point_num = arg_input.q_point_num
 # atom 12
#atom1 = arg_input.atom1;
#atom2 = arg_input.atom2;
atom12_list = arg_input.atom12_list;
hdftmpdir = arg_input.hdftmpdir;

# orbital mask
orbital_selection_option = arg_input.orbital_selection_option;

orbital_selection1_list = arg_input.orbital_selection1_list;
orbital_selection1_names = arg_input.orbital_selection1_names;
orbital_selection2_list = arg_input.orbital_selection2_list;
orbital_selection2_names = arg_input.orbital_selection2_names;

orbital_selection3_list = arg_input.orbital_selection3_list;
orbital_selection3_names = arg_input.orbital_selection3_names;
orbital_selection4_list = arg_input.orbital_selection4_list;
orbital_selection4_names = arg_input.orbital_selection4_names;

println(orbital_selection1_list," ",orbital_selection1_names)
println(orbital_selection2_list," ",orbital_selection2_names)
println(orbital_selection3_list," ",orbital_selection3_names)
println(orbital_selection4_list," ",orbital_selection4_names)
@assert(length(orbital_selection1_list) == length(orbital_selection1_names));
@assert(length(orbital_selection2_list) == length(orbital_selection2_names));
@assert(length(orbital_selection3_list) == length(orbital_selection3_names));
@assert(length(orbital_selection4_list) == length(orbital_selection4_names));
# Band selection
if haskey(arg_input.Optional,"band_selection")
  band_selection_on =  arg_input.Optional["band_selection"]
  band_selection_lower =  arg_input.Optional["band_selection_boundary"][1]
  band_selection_upper =  arg_input.Optional["band_selection_boundary"][2]
end

if ((DFTcommon.unmask == orbital_selection_option) || (DFTcommon.mask == orbital_selection_option) )
 #orbital_selection_name = arg_input.orbital_selection_name
 orbital_selection_on = true
end

# orbital orbital_reassign
basisTransform_rule = basisTransform_rule_type()
if haskey(arg_input.Optional,"basisTransform")
  basisTransform_rule = arg_input.Optional["basisTransform"]
end

println(DFTcommon.bar_string) # print ====...====
println("atom12_list: ",atom12_list)
println("q_point_num ",q_point_num, "\tk_point_num ",k_point_num)
println(string("DFT_type ",DFT_type))
println(string("orbital_selection_option ",orbital_selection_option))
println("mask1list ",orbital_selection1_list,"\tmask2list ",orbital_selection2_list)
println("basisTransform", basisTransform_rule)

## 1.4 Set caluations type and ouput folder
cal_type = "bandplot" # xq, ...

if (DFTcommon.Wannier90 == DFT_type)
  cal_type = string(cal_type,".wannier")
end

root_dir = dirname(result_file)
result_file_only = splitext(basename(result_file))
cal_name = result_file_only[1];
jq_output_dir =  joinpath(root_dir,string(cal_type,"_" ,ChemP_delta_ev))
if (!isdir(jq_output_dir))
  mkdir(jq_output_dir)
end
if ("" == hdftmpdir || !isdir(hdftmpdir) )
  hdftmpdir = jq_output_dir
end
hdf_cache_name = joinpath(hdftmpdir,string(cal_name,".hdf5"))
println(hdf_cache_name)
println(DFTcommon.bar_string) # print ====...====


##############################################################################
## 2. Calculate & Store k,q points information
## 2.1 Set Input info
## 2.2 Generate k,q points
## 2.3 Calculate Eigenstate & Store Eigenstate into file in HDF5 format
## 2.4 Send Eigenstate info to child processes
##############################################################################

## 2.1 Set Input info
#scf_r = [];
#basisTransform_rule = basisTransform_rule_type()


hamiltonian_info = [];
if (DFTcommon.OpenMX == DFT_type || DFTcommon.EcalJ == DFT_type)

  hamiltonian_info = set_current_dftdataset(result_file, result_file_dict, DFT_type, spin_type,basisTransform_rule)
elseif (DFTcommon.Wannier90 == DFT_type)
  atomnum = arg_input.Wannier_Optional_Info.atomnum
  atompos = arg_input.Wannier_Optional_Info.atompos
  atoms_orbitals_list = arg_input.Wannier_Optional_Info.atoms_orbitals_list

  #hamiltonian_info = DFTforge.read_dftresult(result_file,DFT_type,Wannier90_type,atoms_orbitals_list,atomnum,atompos,basisTransform_rule)
  hamiltonian_info = set_current_dftdataset(result_file, result_file_dict, DFT_type,Wannier90_type,spin_type,atoms_orbitals_list,atomnum,atompos,basisTransform_rule)
  #hamiltonian_info = set_current_dftdataset(scf_r, DFT_type, spin_type,basisTransform_rule)
elseif (DFTcommon.PlainwaveLobster == DFT_type)
  print("PlainwaveLobster")
  atomnum = arg_input.Wannier_Optional_Info.atomnum
  atompos = arg_input.Wannier_Optional_Info.atompos
  tv  = arg_input.Optional["tv"]
  optionalInfo = Dict{AbstractString,Any}()
  optionalInfo["tv"] = tv
  optionalInfo["atomnum"] = atomnum
  optionalInfo["atompos"] = atompos
  optionalInfo["ChemP"] = arg_input.Optional["ChemP"]
  println("optionalInfo[ChemP] ",optionalInfo["ChemP"])
  println(tv)
  hamiltonian_info = set_current_dftdataset(result_file, result_file_dict, DFT_type, spin_type,basisTransform_rule,optionalInfo)
end


DFTforge.pwork(set_current_dftdataset,(hamiltonian_info, 1));

#H = DFTforge.cal_colinear_Hamiltonian(k_point_test, DFT_type, hamiltonian_info, 1)
H_type = Array{Array{Array{ComplexF64,2}}}
function Overlap_Band!(H::H_type, R_vector_mat, spin::Int,
    Hout::Array{Complex_my,2},
    TotalOrbitalNum::Int64,k1,k2,k3)
    #HWR_mat_list::Array{Array{Complex_my,2}}
    FNAN = length(H[spin]);
    @assert(size(R_vector_mat[spin])[1] == FNAN)
    k_point::Array{Float_my,1} = [k1,k2,k3];
    #println(FNAN," spin ", 1)
    for LB_AN = 1:FNAN
        kRn::Float_my = -sum(R_vector_mat[spin][LB_AN,:].*k_point);
        Hout[:,:] += H[spin][LB_AN].* (cos(2.0*pi*kRn)+sin(2.0*pi*kRn)*im);
    end
end

TotalOrbitalNum = sum(hamiltonian_info.scf_r.Total_NumOrbs)
lobster_r = hamiltonian_info.scf_r

S = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum);
#k_point = k_point_test


############################
# 2.2.1 Get kpoint path
############################

K_point_groups = arg_input.Optional["bandplot"]
K_point_groups[1].K_point_list

K_point_list_all = Array{k_point_Tuple}(undef,0);

for (i,v) in enumerate(K_point_groups)
  map(x -> push!(K_point_list_all,x) ,v.K_point_list)
end


K_point_list_all_mat =  zeros(length(K_point_list_all),3)
K_point_name_list = Array{Tuple{Int64,String} }(undef,0)

Cnt = 1
for (i,v) in enumerate(K_point_groups)
    
    for (i2,K_point) in enumerate(v.K_point_list)
        if(1 == i2)
            push!(K_point_name_list,(Cnt,v.K_start_point_name))
        end
        K_point_cat = hamiltonian_info.scf_r.rv * [K_point[1],K_point[2],K_point[3]];
        K_point_cat = [K_point[1],K_point[2],K_point[3]]' * hamiltonian_info.scf_r.rv;
        K_point_list_all_mat[Cnt,:] = K_point_cat;
        
        if(i2 == length(v.K_point_list) && i == length(K_point_groups) )
            push!(K_point_name_list,(Cnt,v.K_end_point_name))
        end
        global Cnt += 1
    end
end

# Convert K_point distances
K_point_dist_list = zeros(length(K_point_list_all))

for i in 2:length(K_point_list_all)
    dist = sum( (K_point_list_all_mat[i-1,:] -     K_point_list_all_mat[i,:]).^2)
    K_point_dist_list[i] = sqrt(dist);
end
#
K_point_dist_list = cumsum(K_point_dist_list);

#K_point_name_list


K_point_tick_pos = map(x -> K_point_dist_list[x[1]],K_point_name_list)
K_point_tick_name = map(x -> x[2],K_point_name_list)
print(K_point_tick_pos," ",K_point_tick_name)


## 2.3 Calculate Eigenstate & Store Eigenstate into file in HDF5 format
eigenstate_cache = cachecal_all_Qpoint_eigenstats(K_point_list_all,hdf_cache_name);
## 2.4 Send Eigenstate info to child processes
DFTforge.pwork(cacheset,eigenstate_cache)
#tic();
DFTforge.pwork(cacheread_lampup,K_point_list_all)
#toc();


if (DFTcommon.non_colinear_type != eigenstate_cache.spin_type)


    cache_index = 1;
    eigenvalues_spin_up = zeros(length(K_point_list_all),(eigenstate_cache.TotalOrbitalNum)) .- 100;
    eigenvalues_spin_down = zeros(length(K_point_list_all),(eigenstate_cache.TotalOrbitalNum)) .- 100;
    #tic()
    for (i,k_point) in enumerate(K_point_list_all)
        #eigenstate = cal_colinear_eigenstate(v,[1,2]);
        eigenstate_k_up = cacheread_eigenstate(k_point,1,cache_index)
        eigenstate_k_dn = cacheread_eigenstate(k_point,2,cache_index)
        
        #println(length(eigenstate_k_up.Eigenvalues))
        #println(size(eigenstate_k_down.Eigenvalues))
        energy_idx_num_up = length(eigenstate_k_up.Eigenvalues)
        energy_idx_num_dn = length(eigenstate_k_dn.Eigenvalues)
        eigenvalues_spin_up[i,1:energy_idx_num_up] = eigenstate_k_up.Eigenvalues;
        eigenvalues_spin_down[i,1:energy_idx_num_dn] = eigenstate_k_dn.Eigenvalues;
    end
    #toc()

    Plots.plot!(K_point_dist_list, eigenvalues_spin_up .- get_ChempP() ,xaxis=("K-path"),yaxis=("eV",(-5,5))   ,    xticks=(K_point_tick_pos, K_point_tick_name),color=Plots.RGBA(1,0,0,1),lab="",ms=2)
    Plots.plot!(K_point_dist_list, eigenvalues_spin_down .- get_ChempP()  ,xaxis=("K-path"),yaxis=("eV",(-2.5,2.5))   ,    xticks=(K_point_tick_pos, K_point_tick_name),color=Plots.RGBA(0,0,1,1),lab="",ms=2)
end

Plots.savefig(joinpath(jq_output_dir,string("all_bands",".pdf") ))