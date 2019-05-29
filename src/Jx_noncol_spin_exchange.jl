################################################################################
using ProgressMeter
import DFTforge
using DFTforge.DFTrefinery
using DFTforge.DFTcommon
#import MAT
using Statistics
using Distributed;

X_VERSION = VersionNumber("0.1.0-dev.test.20170701");
println(" X_VERSION: ",X_VERSION)
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
#orbital_mask1 = Array{Int64,1}();
#orbital_mask2 = Array{Int64,1}();
orbital_mask1_list = Array{Array{Int}}(undef,0);
orbital_mask1_names = Array{AbstractString}(undef,0);
orbital_mask2_list = Array{Array{Int}}(undef,0);
orbital_mask2_names = Array{AbstractString}(undef,0);

orbital_mask3_list = Array{Array{Int}}(undef,0);
orbital_mask3_names = Array{AbstractString}(undef,0);
orbital_mask4_list = Array{Array{Int}}(undef,0);
orbital_mask4_names = Array{AbstractString}(undef,0);

orbital_mask_option = DFTcommon.nomask;
orbital_mask_on = false;

k_point_num = [3,3,3]
q_point_num = [3,3,3]
ChemP_delta_ev = 0.0
DFT_type = DFTcommon.OpenMX

## 1.2 Read input from argument & TOML file
#ARGS=[""]
arg_input = DFTcommon.Arg_Inputs();
arg_input = parse_input(ARGS,arg_input)
#arg_input.TOMLinput = "nio_J_wannier.toml" # Debug
#arg_input.TOMLinput = "nio_J_openmx.toml" # Debug
#arg_input.TOMLinput = "../examples/NiO/16cell_NC/s2p2d1_TMO.toml"
#arg_input.TOMLinput = "/home/users1/bluehope/work_local/NC_testset/NiO/U5.5/AF2_16cell/Ni6.0_O5.0-s2p2d2f1-4.180-k10/90_0/s2p2d1_TMO.toml"
#arg_input.TOMLinput = "../examples/NiO/16cell_NC_FM/s2p2d1_TMO.toml"
#arg_input.TOMLinput = "/home/users1/bluehope/work_local/NC_testset/Fe/s2p2d2f1/90_90_90_90_noSoC/s2p2d2f1_TM.toml"

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
atom1 = arg_input.atom1;
atom2 = arg_input.atom2;
atom12_list = arg_input.atom12_list;
hdftmpdir = arg_input.hdftmpdir;

Hmode = arg_input.Hmode;
# orbital mask
orbital_mask_option = arg_input.orbital_mask_option;

orbital_mask1_list = arg_input.orbital_mask1_list;
orbital_mask1_names = arg_input.orbital_mask1_names;
orbital_mask2_list = arg_input.orbital_mask2_list;
orbital_mask2_names = arg_input.orbital_mask2_names;

orbital_mask3_list = arg_input.orbital_mask3_list;
orbital_mask3_names = arg_input.orbital_mask3_names;
orbital_mask4_list = arg_input.orbital_mask4_list;
orbital_mask4_names = arg_input.orbital_mask4_names;

println(orbital_mask1_list," ",orbital_mask1_names)
println(orbital_mask2_list," ",orbital_mask2_names)
println(orbital_mask3_list," ",orbital_mask3_names)
println(orbital_mask4_list," ",orbital_mask4_names)
@assert(length(orbital_mask1_list) == length(orbital_mask1_names));
@assert(length(orbital_mask2_list) == length(orbital_mask2_names));
@assert(length(orbital_mask3_list) == length(orbital_mask3_names));
@assert(length(orbital_mask4_list) == length(orbital_mask4_names));

if ((DFTcommon.unmask == orbital_mask_option) || (DFTcommon.mask == orbital_mask_option) )
 #orbital_mask_name = arg_input.orbital_mask_name
 orbital_mask_on = true
end


# orbital orbital_reassign
basisTransform_rule = basisTransform_rule_type()
if haskey(arg_input.Optional,"basisTransform")
 basisTransform_rule = arg_input.Optional["basisTransform"]
end
println(DFTcommon.bar_string) # print ====...====
println(atom12_list)
println("q_point_num ",q_point_num, "\tk_point_num ",k_point_num)
println(string("DFT_type ",DFT_type))
println(string("orbital_mask_option ",orbital_mask_option))
println("mask1list ",orbital_mask1_list,"\tmask2list ",orbital_mask2_list)
println("basisTransform", basisTransform_rule)

## 1.4 Set caluations type and ouput folder
cal_type = "jq.opposite" # xq, ...
cal_type = "jq.same" # xq, ...
if (DFTcommon.nc_allH == Hmode)
  cal_type = string(cal_type,".nc_allH")
elseif (DFTcommon.nc_realH_only == Hmode)
  cal_type = string(cal_type,".nc_realH_only")
elseif (DFTcommon.nc_imagH_only == Hmode)
  cal_type = string(cal_type,".nc_imagH_only")
end

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
hamiltonian_info = [];
if (DFTcommon.OpenMX == DFT_type || DFTcommon.EcalJ == DFT_type)
  #scf_r = set_current_dftdataset(result_file, result_file_dict, DFTcommon.OpenMX, DFTcommon.non_colinear_type)
  hamiltonian_info = set_current_dftdataset(result_file,result_file_dict, DFT_type, spin_type,basisTransform_rule)
elseif (DFTcommon.Wannier90 == DFT_type)
  atomnum = arg_input.Wannier_Optional_Info.atomnum
  atompos = arg_input.Wannier_Optional_Info.atompos
  atoms_orbitals_list = arg_input.Wannier_Optional_Info.atoms_orbitals_list

  #scf_r = DFTforge.read_dftresult(result_file, result_file_dict, DFT_type,Wannier90_type,atoms_orbitals_list,atomnum,atompos)
  #scf_r = set_current_dftdataset(scf_r, result_file_dict, DFT_type, DFTcommon.non_colinear_type)
  hamiltonian_info = set_current_dftdataset(result_file,result_file_dict,
    DFT_type,Wannier90_type,spin_type,atoms_orbitals_list,atomnum,atompos,basisTransform_rule)
end

DFTforge.pwork(set_current_dftdataset,(hamiltonian_info, 1));

## 2.2 Generate k,q points
k_point_list = kPoint_gen_GammaCenter(k_point_num);
q_point_list = kPoint_gen_GammaCenter(q_point_num);

#k_point_list = kPoint_gen_EquallySpaced(k_point_num);
#q_point_list = kPoint_gen_EquallySpaced(q_point_num);

(kq_point_list,kq_point_int_list) = q_k_unique_points(q_point_list,k_point_list)
println(string(" kq_point_list ",length(kq_point_list)))
println(string(" q point ",length(q_point_list) ," k point ",length(k_point_list)))

## 2.3 Calculate Eigenstate & Store Eigenstate into file in HDF5 format
eigenstate_cache = cachecal_all_Qpoint_eigenstats(kq_point_list,hdf_cache_name);
gc();

## 2.4 Send Eigenstate info to child processes
DFTforge.pwork(cacheset,eigenstate_cache)
#tic();
@time DFTforge.pwork(cacheread_lampup,kq_point_list)
#toc();
################################################################################

##############################################################################
## 3. Setup extra infos (orbital, chemp shift)
##############################################################################
@everywhere function init_orbital_mask(orbital_mask_input::orbital_mask_input_Type)
    global orbital_mask1,orbital_mask2,orbital_mask_on
    global orbital_mask3,orbital_mask4
    orbital_mask1 = Array{Int64,1}();
    orbital_mask2 = Array{Int64,1}();

    orbital_mask3 = Array{Int64,1}();
    orbital_mask4 = Array{Int64,1}();
    if (orbital_mask_input.orbital_mask_on)
        orbital_mask1 = orbital_mask_input.orbital_mask1;
        orbital_mask2 = orbital_mask_input.orbital_mask2;

        orbital_mask3 = orbital_mask_input.orbital_mask3;
        orbital_mask4 = orbital_mask_input.orbital_mask4;
        orbital_mask_on = true;
    else
        orbital_mask_on = false;
    end
    #println(orbital_mask1)
end
@everywhere function init_variables(input)
  global ChemP_delta_ev,Hmode,theta_phi_list;
  Input_ChemP_delta_ev = input[1];
  Input_Hmode          = input[2];
  theta_phi_list       = input[3];

  ChemP_delta_ev = Input_ChemP_delta_ev;
  Hmode = Input_Hmode;
end


@everywhere #import MAT
#DFTforge.pwork(Init_SmallHks,(atom1,atom2))
#=
@everywhere function init_Hks(input)
  global Hks_0;
    result_index = input[1]
    Hmode = input[2]
    Hks_0  = cal_Hamiltonian(1,result_index,Hmode)
  end
DFTforge.pwork(init_Hks,(1,Hmode))
=#

##############################################################################
## Physical properties calculation define section
## 4. Magnetic exchange function define
## 4.1 Do K,Q sum
## 4.2 reduce K,Q to Q space
##############################################################################
num_return = 15;
## Test codes

## 4. Magnetic exchange function define
@everywhere function Magnetic_Exchange_J_noncolinear(input::Job_input_kq_atom_list_Type)
  global orbital_mask1,orbital_mask2,orbital_mask_on
  global orbital_mask3,orbital_mask4

  global Hks_0;
  global Hmode;
  global ChemP_delta_ev;
  #global SmallHks;
  #global Hks_0;
  ############################################################################
  ## Accessing Data Start
  ############################################################################
  # Common input info
  num_return = 15;
  k_point::DFTforge.k_point_Tuple =  input.k_point
  kq_point::DFTforge.k_point_Tuple =  input.kq_point
  spin_type::DFTforge.SPINtype = input.spin_type;

  #atom1::Int = input.atom12[1][1];
  #atom2::Int = input.atom12[1][2];
  atom12_list::Vector{Tuple{Int64,Int64}} = input.atom12_list;
  result_mat = zeros(Complex_my,num_return,length(atom12_list))

  #atom1::Int = input.atom1;
  #atom2::Int = input.atom2;

  result_index =  input.result_index;
  cache_index = input.cache_index;

  # Get Chemp, E_temp
  TotalOrbitalNum = cacheread(cache_index).TotalOrbitalNum;
  TotalOrbitalNum2 = TotalOrbitalNum;
  if (DFTcommon.non_colinear_type == spin_type)
    TotalOrbitalNum2 = 2*TotalOrbitalNum;
  end

  ChemP = get_dftdataset(result_index).scf_r.ChemP + ChemP_delta_ev;
  E_temp = get_dftdataset(result_index).scf_r.E_Temp;

  # Get EigenState
  #eigenstate_k_up::Kpoint_eigenstate  = cacheread_eigenstate(k_point,1,cache_index)
  eigenstate_kq = cacheread_eigenstate(kq_point,1,cache_index)
  eigenstate_k  = cacheread_eigenstate(k_point,1,cache_index)
  #eigenstate_kq_down::Kpoint_eigenstate = cacheread_eigenstate(kq_point,2,cache_index)

  # Get Hamiltonian
  #Hks_updown = cacheread_Hamiltonian(1,cache_index)
  #Hks_updown = cal_Hamiltonian(1,1)
  #Hks_updown = copy(Hks_0);

  #Hks_updown = cacheread_Hamiltonian(1,cache_index)
  #Hks_updown = copy(Hks_0);
  #@assert( sum(real(Hks_updown - Hks_down2)) < 10.0^-4 )
  Hks_updown_k  = cacheread_Hamiltonian(k_point, Hmode,cache_index)::Hamiltonian_type;
  Hks_updown_kq  = cacheread_Hamiltonian(kq_point, Hmode,cache_index)::Hamiltonian_type;
  @assert(size(Hks_updown_k)[1] == TotalOrbitalNum2);
  #Hks_k_down::Hamiltonian_type = cacheread_Hamiltonian(k_point,2,cache_index)
  #Hks_kq_down::Hamiltonian_type = cacheread_Hamiltonian(kq_point,2,cache_index)

  (orbitalStartIdx_list,orbitalNums) = cacheread_atomsOrbital_lists(cache_index)

  #En_k_up::Array{Float_my,1} = eigenstate_k_up.Eigenvalues;
  Em_kq = eigenstate_kq.Eigenvalues;
  En_k  = eigenstate_k.Eigenvalues;
  #Em_kq_down::Array{Float_my,1} = eigenstate_kq_down.Eigenvalues;

  #Es_n_k_up::Array{Complex_my,2} = eigenstate_k_up.Eigenstate;
  Es_m_kq  = eigenstate_kq.Eigenstate;
  Es_n_k  = eigenstate_k.Eigenstate;

  #for mask2 in orbital_mask2
  #    Es_m_kq_atom2[mask2,:]=0.0
  #end
  #Plots.heatmap(real(Es_n_k))
  #Es_m_kq_down::Array{Complex_my,2} = eigenstate_kq_down.Eigenstate;
  Fftn_k  = 1.0./(exp.( ((En_k)  - ChemP)/(kBeV*E_temp)) + 1.0 );
  Fftm_kq = 1.0./(exp.( ((Em_kq) - ChemP)/(kBeV*E_temp)) + 1.0 );
  # Index convention: dFnk[nk,mkq]
  dFnk_Fmkq  =
    Fftn_k*ones(1,TotalOrbitalNum2)  - ones(TotalOrbitalNum2,1)*Fftm_kq' ;
  # Index convention: Enk_Emkq[nk,mkq]
  Enk_Emkq =
    En_k[:]*ones(1,TotalOrbitalNum2) - ones(TotalOrbitalNum2,1)*Em_kq[:]' ;
  Enk_Emkq += im*0.0001;

  # Index convention: part1[nk,mkq]
  part1 = dFnk_Fmkq./(-Enk_Emkq);

  for (atom12_i,atom12) in enumerate(atom12_list)
    atom1 = atom12[1]
    atom2 = atom12[2]

    Es_m_kq_atom1 = copy(Es_m_kq);
    Es_m_kq_atom2 = copy(Es_m_kq);
    Es_n_k_atom1 = copy(Es_n_k);
    Es_n_k_atom2 = copy(Es_n_k);
    if (length(orbital_mask1)>0)
      orbital_mask1_tmp = collect(1:orbitalNums[atom1]);
      for orbit1 in orbital_mask1
          deleteat!(orbital_mask1_tmp, find(orbital_mask1_tmp.==orbit1))
      end
      Es_n_k_atom1[orbitalStartIdx_list[atom1]+orbital_mask1_tmp,:]=0.0; # up
      Es_n_k_atom1[TotalOrbitalNum + orbitalStartIdx_list[atom1]+orbital_mask1_tmp,:]=0.0; # down
    end
    if (length(orbital_mask2)>0)
      orbital_mask2_tmp = collect(1:orbitalNums[atom2]);
      for orbit2 in orbital_mask2
          deleteat!(orbital_mask2_tmp, find(orbital_mask2_tmp.==orbit2))
      end
      Es_m_kq_atom2[orbitalStartIdx_list[atom2]+orbital_mask2_tmp,:]=0.0; # up
      Es_m_kq_atom2[TotalOrbitalNum + orbitalStartIdx_list[atom2]+orbital_mask2_tmp,:]=0.0; # down
    end



    atom1_orbits_up   = orbitalStartIdx_list[atom1] + (1:orbitalNums[atom1])
    atom1_orbits_down = TotalOrbitalNum + atom1_orbits_up
    atom2_orbits_up   = orbitalStartIdx_list[atom2] + (1:orbitalNums[atom2])
    atom2_orbits_down = TotalOrbitalNum + atom2_orbits_up

    atom1_orbits_updown = [atom1_orbits_up;atom1_orbits_down];
    atom2_orbits_updown = [atom2_orbits_up;atom2_orbits_down];
    atom1_orbits_downup = [atom1_orbits_down;atom1_orbits_up];
    atom2_orbits_downup = [atom2_orbits_down;atom2_orbits_up];


    atom1_orbits_up_rel = (1:orbitalNums[atom1]);
    atom1_orbits_down_rel = orbitalNums[atom1]+(1:orbitalNums[atom1]);
    atom2_orbits_up_rel = (1:orbitalNums[atom2]);
    atom2_orbits_down_rel = orbitalNums[atom2]+(1:orbitalNums[atom2]);
    #Plots.heatmap(real(Hks_updown))

    ############################################################################
    ## Accessing Data End
    ############################################################################
    ## Do auctual calucations
    Hks_updown_k_atom1  = copy(Hks_updown_k[atom1_orbits_updown,atom1_orbits_updown]);
    Hks_updown_k_atom1_flip = copy(Hks_updown_k_atom1);

    # spin mirror flip
    # rho_flip[up,up]     =  rho_orig[down,down]
    # rho_flip[down,down] =  rho_orig[up,up]
    # rho_flip[up,down]   = -rho_orig[up,down]
    # rho_flip[down,up]   = -rho_orig[down,up]

    # spin 1 flip
    Hks_updown_k_atom1_flip[atom1_orbits_up_rel, atom1_orbits_up_rel] = Hks_updown_k_atom1[atom1_orbits_down_rel, atom1_orbits_down_rel]
    Hks_updown_k_atom1_flip[atom1_orbits_down_rel, atom1_orbits_down_rel] = Hks_updown_k_atom1[atom1_orbits_up_rel, atom1_orbits_up_rel]

    Hks_updown_k_atom1_flip[atom1_orbits_down_rel, atom1_orbits_up_rel] = -Hks_updown_k_atom1[atom1_orbits_down_rel, atom1_orbits_up_rel]
    Hks_updown_k_atom1_flip[atom1_orbits_up_rel, atom1_orbits_down_rel] = -Hks_updown_k_atom1[atom1_orbits_up_rel, atom1_orbits_down_rel]

    # spin 1 H_rot - H_orig
    Hks_updown_k_atom1_diff = Hks_updown_k_atom1_flip - Hks_updown_k_atom1;

    # spin 2 flip
    Hks_updown_kq_atom2 = copy(Hks_updown_kq[atom2_orbits_updown,atom2_orbits_updown]);
    Hks_updown_kq_atom2_flip = copy(Hks_updown_kq_atom2);

    Hks_updown_kq_atom2_flip[atom1_orbits_up_rel, atom1_orbits_up_rel] = Hks_updown_kq_atom2[atom1_orbits_down_rel, atom1_orbits_down_rel]
    Hks_updown_kq_atom2_flip[atom1_orbits_down_rel, atom1_orbits_down_rel] = Hks_updown_kq_atom2[atom1_orbits_up_rel, atom1_orbits_up_rel]

    Hks_updown_kq_atom2_flip[atom1_orbits_down_rel,atom1_orbits_up_rel] = -Hks_updown_kq_atom2[atom1_orbits_down_rel,atom1_orbits_up_rel]
    Hks_updown_kq_atom2_flip[atom1_orbits_up_rel,atom1_orbits_down_rel] = -Hks_updown_kq_atom2[atom1_orbits_up_rel,atom1_orbits_down_rel]

    # spin 2 H_rot - H_orig
    Hks_updown_kq_atom2_diff = Hks_updown_kq_atom2_flip - Hks_updown_kq_atom2;


    Vz_1 = 0.5*(Hks_updown_k_atom1_diff[atom1_orbits_up_rel,atom1_orbits_up_rel] - Hks_updown_k_atom1_diff[atom1_orbits_down_rel,atom1_orbits_down_rel])/2.0
    #Voff_1 = Hks_updown[atom1_orbits_up,atom1_orbits_down]
    #Voff_1 = conj(Hks_updown_k[atom1_orbits_up,atom1_orbits_down]) + Hks_updown_k[atom1_orbits_down,atom1_orbits_up]
    Vx_1 = 0.5*(Hks_updown_k_atom1_diff[atom1_orbits_up_rel, atom1_orbits_down_rel] + Hks_updown_k_atom1_diff[atom1_orbits_down_rel, atom1_orbits_up_rel])/2.0;
    Vy_1 = 0.5*(-Hks_updown_k_atom1_diff[atom1_orbits_up_rel, atom1_orbits_down_rel] + Hks_updown_k_atom1_diff[atom1_orbits_down_rel, atom1_orbits_up_rel])/(2.0*im);

    #Vx_1 = 0.5*real(Hks_updown_k_atom1_diff[atom1_orbits_up_rel, atom1_orbits_down_rel]);
    #Vy_1 = 0.5*imag(Hks_updown_k_atom1_diff[atom1_orbits_up_rel, atom1_orbits_down_rel])*(1*im)
    #Plots.heatmap(real(Vz_1))
    #Plots.heatmap(real(Vy_1))

    Vz_2 = 0.5*(Hks_updown_kq_atom2_diff[atom2_orbits_up_rel, atom2_orbits_up_rel] - Hks_updown_kq_atom2_diff[atom2_orbits_down_rel, atom2_orbits_down_rel])/2.0;
    #Voff_2 = conj(Hks_updown_kq[atom2_orbits_up,atom2_orbits_down]) + Hks_updown_kq[atom2_orbits_down,atom2_orbits_up]
    # Voff_2 = (Hks_updown[atom2_orbits_up,atom2_orbits_down])
           #+ Hks_updown[atom2_orbits_down,atom2_orbits_up]

    Vx_2 = 0.5*(Hks_updown_kq_atom2_diff[atom2_orbits_up_rel, atom2_orbits_down_rel] + Hks_updown_kq_atom2_diff[atom2_orbits_down_rel, atom2_orbits_up_rel])/2.0;
    Vy_2 = 0.5*(-Hks_updown_kq_atom2_diff[atom2_orbits_up_rel, atom2_orbits_down_rel] + Hks_updown_kq_atom2_diff[atom2_orbits_down_rel, atom2_orbits_up_rel] )/(2.0*im);

    #Vx_2 = 0.5*real(Hks_updown_kq_atom2_diff[atom2_orbits_up_rel, atom2_orbits_down_rel]);
    #Vy_2 = -0.5*imag(Hks_updown_kq_atom2_diff[atom2_orbits_up_rel, atom2_orbits_down_rel])*(1*im)

    #Plots.heatmap(real(part1))

    G1V1_z = (Es_n_k_atom1[atom1_orbits_down,:]'  *Vz_1* Es_m_kq_atom1[atom1_orbits_up,:]);
    G2V2_z = (Es_m_kq_atom2[atom2_orbits_up,:]'   *Vz_2.' * Es_n_k_atom2[atom2_orbits_down,:]);
    #G2V2_z = (Es_m_kq[atom2_orbits_up,:]' * Vz_2 * Es_n_k[atom2_orbits_down,:]);
    G1V1_x = (Es_n_k_atom1[atom1_orbits_down,:]'  *Vx_1* Es_m_kq_atom1[atom1_orbits_up,:]);
    G2V2_x = (Es_m_kq_atom2[atom2_orbits_up,:]'   *Vx_2.' * Es_n_k_atom2[atom2_orbits_down,:]);

    G1V1_y = (Es_n_k_atom1[atom1_orbits_down,:]'  *Vy_1* Es_m_kq_atom1[atom1_orbits_up,:]);
    G2V2_y = (Es_m_kq_atom2[atom2_orbits_up,:]'   *Vy_2.' * Es_n_k_atom2[atom2_orbits_down,:]);


    G1V1_z2 = (Es_n_k_atom1[atom1_orbits_up,:]'   * Vz_1   * Es_m_kq_atom1[atom1_orbits_down,:]);
    G2V2_z2 = (Es_m_kq_atom2[atom2_orbits_up,:]' * -Vz_2.'  * Es_n_k_atom2[atom2_orbits_down,:]);

    G1V1_x2 = (Es_n_k_atom1[atom1_orbits_up,:]'   * Vx_1   * Es_m_kq_atom1[atom1_orbits_down,:]);
    G2V2_x2 = (Es_m_kq_atom2[atom2_orbits_up,:]' * -Vx_2.'  * Es_n_k_atom2[atom2_orbits_down,:]);

    G1V1_y2 = (Es_n_k_atom1[atom1_orbits_up,:]'   * Vy_1   * Es_m_kq_atom1[atom1_orbits_down,:]);
    G2V2_y2 = (Es_m_kq_atom2[atom2_orbits_up,:]' * Vy_2.'  * Es_n_k_atom2[atom2_orbits_down,:]);

#=
    G1V1_all = (Es_n_k_atom1[atom1_orbits_updown,:]'  *Hks_updown_k_atom1_diff* Es_m_kq_atom1[atom1_orbits_updown,:]);
    G2V2_all = (Es_m_kq_atom2[atom2_orbits_updown,:]'   *Hks_updown_kq_atom2_diff.' * Es_n_k_atom2[atom2_orbits_updown,:]);
    =#
    #G2V2_z = (Es_m_kq[atom2_orbits_up,:]' * Vz_2 * Es_n_k[atom2_orbits_down,:]);

    # Index convention: J_ij_(xyz,xyz)[nk,mkq]
    J_ij_xx =  0.5* sum(part1.* G1V1_x .* transpose(G2V2_x) );
    J_ij_xy =  0.5* sum(part1.* G1V1_x .* transpose(G2V2_y) );
    J_ij_xz =  0.5* sum(part1.* G1V1_x .* transpose(G2V2_z) );

    J_ij_yx =  0.5* sum(part1.* G1V1_y .* transpose(G2V2_x) );
    J_ij_yy =  0.5* sum(part1.* G1V1_y .* transpose(G2V2_y) );
    J_ij_yz =  0.5* sum(part1.* G1V1_y .* transpose(G2V2_z) );

    J_ij_zx =  0.5* sum(part1.* G1V1_z .* transpose(G2V2_x) );
    J_ij_zy =  0.5* sum(part1.* G1V1_z .* transpose(G2V2_y) );
    J_ij_zz =  0.5* sum(part1.* G1V1_z .* transpose(G2V2_z) );
    #=
    J_ij_all = sum(part1.* G1V1_all .* transpose(G2V2_all) );
    =#
    J_ij_xx2 =  0.5* sum(part1.* G1V1_x2 .* transpose(G2V2_x2) );
    J_ij_yy2 =  0.5* sum(part1.* G1V1_y2 .* transpose(G2V2_y2) );
    J_ij_zz2 =  0.5* sum(part1.* G1V1_z2 .* transpose(G2V2_z2) );


    X_ij = sum(part1.*
    (Es_n_k[atom1_orbits_down,:]' * Es_m_kq[atom1_orbits_up,:]) .*
    transpose(Es_m_kq[atom2_orbits_up,:]' * Es_n_k[atom2_orbits_down,:]) );


    X_ij_charge_updown = sum(part1.*
    (Es_n_k[atom1_orbits_up,:]' * Es_m_kq[atom1_orbits_up,:]) .*
    transpose(Es_m_kq[atom2_orbits_down,:]' * Es_n_k[atom2_orbits_down,:]) );

    X_ij_charge_downup = sum(part1.*
    (Es_n_k[atom1_orbits_down,:]' * Es_m_kq[atom1_orbits_down,:]) .*
    transpose(Es_m_kq[atom2_orbits_up,:]' * Es_n_k[atom2_orbits_up,:]) );



    result_subset = [
    J_ij_xx,J_ij_xy,J_ij_xz,
    J_ij_yx,J_ij_yy,J_ij_yz,
    J_ij_zx,J_ij_zy,J_ij_zz, J_ij_xx2, J_ij_yy2, J_ij_zz2, X_ij, X_ij_charge_updown, X_ij_charge_downup];
    result_mat[:,atom12_i] = result_subset
    #=
    println(atom12)
    println(  result_mat[:,atom12_i])
    MAT.matwrite(string("/home/users1/bluehope/Dropbox/shared/DFT-forge/source/debug/",atom12_i,".mat")
      ,Dict("Vz_1" =>Vz_1
      ,"Vy_1" => Vy_1
      ,"Vx_1" => Vx_1
      ,"Voff_1" => Voff_1
      ,"Hks_updown" => Hks_updown
      ,"atom1_orbits_up" => collect(atom1_orbits_up)
      ,"atom1_orbits_down" => collect(atom1_orbits_down)
      ,"result_subset" => result_subset
      #,"k_point" => k_point
      #,"kq_point" => kq_point
      )
      );
      =#
  end

  return result_mat #sum(J_ij[:]);
end


num_return = 14;
## Test codes

## 4. Magnetic exchange function define
@everywhere function Magnetic_Exchange_J_noncolinear_fixed_v1(input::Job_input_kq_atom_list_Type)
  global orbital_mask1,orbital_mask2,orbital_mask_on
  global orbital_mask3,orbital_mask4

  global Hks_0;
  global Hmode;
  global ChemP_delta_ev;
  #global SmallHks;
  #global Hks_0;
  ############################################################################
  ## Accessing Data Start
  ############################################################################
  # Common input info
  num_return = 14;
  k_point::DFTforge.k_point_Tuple =  input.k_point
  kq_point::DFTforge.k_point_Tuple =  input.kq_point
  spin_type::DFTforge.SPINtype = input.spin_type;

  #atom1::Int = input.atom12[1][1];
  #atom2::Int = input.atom12[1][2];
  atom12_list::Vector{Tuple{Int64,Int64}} = input.atom12_list;
  result_mat = zeros(Complex_my,num_return,length(atom12_list))

  #atom1::Int = input.atom1;
  #atom2::Int = input.atom2;

  result_index =  input.result_index;
  cache_index = input.cache_index;

  # Get Chemp, E_temp
  TotalOrbitalNum = cacheread(cache_index).TotalOrbitalNum;
  TotalOrbitalNum2 = TotalOrbitalNum;
  if (DFTcommon.non_colinear_type == spin_type)
    TotalOrbitalNum2 = 2*TotalOrbitalNum;
  end

  ChemP = get_dftdataset(result_index).scf_r.ChemP + ChemP_delta_ev;
  E_temp = get_dftdataset(result_index).scf_r.E_Temp;

  # Get EigenState
  #eigenstate_k_up::Kpoint_eigenstate  = cacheread_eigenstate(k_point,1,cache_index)
  eigenstate_kq = cacheread_eigenstate(kq_point,1,cache_index)
  eigenstate_k  = cacheread_eigenstate(k_point,1,cache_index)
  #eigenstate_kq_down::Kpoint_eigenstate = cacheread_eigenstate(kq_point,2,cache_index)

  # Get Hamiltonian
  #Hks_updown = cacheread_Hamiltonian(1,cache_index)
  #Hks_updown = cal_Hamiltonian(1,1)
  #Hks_updown = copy(Hks_0);

  #Hks_updown = cacheread_Hamiltonian(1,cache_index)
  #Hks_updown = copy(Hks_0);
  #@assert( sum(real(Hks_updown - Hks_down2)) < 10.0^-4 )
  Hks_updown_k  = cacheread_Hamiltonian(k_point, Hmode,cache_index)::Hamiltonian_type;
  Hks_updown_kq  = cacheread_Hamiltonian(kq_point, Hmode,cache_index)::Hamiltonian_type;
  @assert(size(Hks_updown_k)[1] == TotalOrbitalNum2);
  #Hks_k_down::Hamiltonian_type = cacheread_Hamiltonian(k_point,2,cache_index)
  #Hks_kq_down::Hamiltonian_type = cacheread_Hamiltonian(kq_point,2,cache_index)

  (orbitalStartIdx_list,orbitalNums) = cacheread_atomsOrbital_lists(cache_index)

  #En_k_up::Array{Float_my,1} = eigenstate_k_up.Eigenvalues;
  Em_kq = eigenstate_kq.Eigenvalues;
  En_k  = eigenstate_k.Eigenvalues;
  #Em_kq_down::Array{Float_my,1} = eigenstate_kq_down.Eigenvalues;

  #Es_n_k_up::Array{Complex_my,2} = eigenstate_k_up.Eigenstate;
  Es_m_kq  = eigenstate_kq.Eigenstate;
  Es_n_k  = eigenstate_k.Eigenstate;

  #for mask2 in orbital_mask2
  #    Es_m_kq_atom2[mask2,:]=0.0
  #end
  #Plots.heatmap(real(Es_n_k))
  #Es_m_kq_down::Array{Complex_my,2} = eigenstate_kq_down.Eigenstate;
  Fftn_k  = 1.0./(exp.( ((En_k)  - ChemP)/(kBeV*E_temp)) + 1.0 );
  Fftm_kq = 1.0./(exp.( ((Em_kq) - ChemP)/(kBeV*E_temp)) + 1.0 );
  # Index convention: dFnk[nk,mkq]
  dFnk_Fmkq  =
    Fftn_k*ones(1,TotalOrbitalNum2)  - ones(TotalOrbitalNum2,1)*Fftm_kq' ;
  # Index convention: Enk_Emkq[nk,mkq]
  Enk_Emkq =
    En_k[:]*ones(1,TotalOrbitalNum2) - ones(TotalOrbitalNum2,1)*Em_kq[:]' ;
  Enk_Emkq += im*0.0001;

  # Index convention: part1[nk,mkq]
  part1_nk_mkq = dFnk_Fmkq./(-Enk_Emkq);

  for (atom12_i,atom12) in enumerate(atom12_list)
    atom1 = atom12[1]
    atom2 = atom12[2]

    Es_m_kq_atom1 = copy(Es_m_kq);
    Es_m_kq_atom2 = copy(Es_m_kq);
    Es_n_k_atom1 = copy(Es_n_k);
    Es_n_k_atom2 = copy(Es_n_k);
    if (length(orbital_mask1)>0)
      orbital_mask1_tmp = collect(1:orbitalNums[atom1]);
      for orbit1 in orbital_mask1
          deleteat!(orbital_mask1_tmp, find(orbital_mask1_tmp.==orbit1))
      end
      Es_n_k_atom1[orbitalStartIdx_list[atom1]+orbital_mask1_tmp,:]=0.0; # up
      Es_n_k_atom1[TotalOrbitalNum + orbitalStartIdx_list[atom1]+orbital_mask1_tmp,:]=0.0; # down
    end
    if (length(orbital_mask2)>0)
      orbital_mask2_tmp = collect(1:orbitalNums[atom2]);
      for orbit2 in orbital_mask2
          deleteat!(orbital_mask2_tmp, find(orbital_mask2_tmp.==orbit2))
      end
      Es_m_kq_atom2[orbitalStartIdx_list[atom2]+orbital_mask2_tmp,:]=0.0; # up
      Es_m_kq_atom2[TotalOrbitalNum + orbitalStartIdx_list[atom2]+orbital_mask2_tmp,:]=0.0; # down
    end



    atom1_orbits_up   = orbitalStartIdx_list[atom1] + (1:orbitalNums[atom1])
    atom1_orbits_down = TotalOrbitalNum + atom1_orbits_up
    atom2_orbits_up   = orbitalStartIdx_list[atom2] + (1:orbitalNums[atom2])
    atom2_orbits_down = TotalOrbitalNum + atom2_orbits_up

    atom1_orbits_updown = [atom1_orbits_up;atom1_orbits_down];
    atom2_orbits_updown = [atom2_orbits_up;atom2_orbits_down];
    atom1_orbits_downup = [atom1_orbits_down;atom1_orbits_up];
    atom2_orbits_downup = [atom2_orbits_down;atom2_orbits_up];


    atom1_orbits_up_rel = (1:orbitalNums[atom1]);
    atom1_orbits_down_rel = orbitalNums[atom1]+(1:orbitalNums[atom1]);
    atom2_orbits_up_rel = (1:orbitalNums[atom2]);
    atom2_orbits_down_rel = orbitalNums[atom2]+(1:orbitalNums[atom2]);
    #Plots.heatmap(real(Hks_updown))

    ############################################################################
    ## Accessing Data End
    ############################################################################
    ## Do auctual calucations
    #Hks_updown_k_atom1  = copy(Hks_updown_k[atom1_orbits_updown,atom1_orbits_updown]);
    #Hks_updown_k_atom1_flip = copy(Hks_updown_k_atom1);

    # spin mirror flip
    # rho_flip[up,up]     =  rho_orig[down,down]
    # rho_flip[down,down] =  rho_orig[up,up]
    # rho_flip[up,down]   = -rho_orig[up,down]
    # rho_flip[down,up]   = -rho_orig[down,up]

    Es_n_k_atom1_structual  = Array{Array{Complex_my,2},2}(2,1);
    Es_m_kq_atom1_structual  = Array{Array{Complex_my,2},2}(2,1);

    Es_n_k_atom2_structual = Array{Array{Complex_my,2},2}(2,1);
    Es_m_kq_atom2_structual = Array{Array{Complex_my,2},2}(2,1);

    Es_n_k_atom1_structual[1,1] = Es_n_k_atom1[atom1_orbits_up,:];
    Es_n_k_atom1_structual[2,1] = Es_n_k_atom1[atom1_orbits_down,:];

    Es_m_kq_atom1_structual[1,1] = Es_m_kq_atom1[atom1_orbits_up,:];
    Es_m_kq_atom1_structual[2,1] = Es_m_kq_atom1[atom1_orbits_down,:];

    Es_n_k_atom2_structual[1,1] = Es_n_k_atom2[atom2_orbits_up,:];
    Es_n_k_atom2_structual[2,1] = Es_n_k_atom2[atom2_orbits_down,:];

    Es_m_kq_atom2_structual[1,1] = Es_m_kq_atom2[atom2_orbits_up,:];
    Es_m_kq_atom2_structual[2,1] = Es_m_kq_atom2[atom2_orbits_down,:];

    H_k_atom1 = Array{Array{Complex_my,2},2}(2,2);
    H_kq_atom2 = Array{Array{Complex_my,2},2}(2,2);
    H_k_atom1[1,1] =  Hks_updown_k[atom1_orbits_up,atom1_orbits_up];
    H_k_atom1[2,2] =  Hks_updown_k[atom1_orbits_down,atom1_orbits_down];
    H_k_atom1[1,2] =  Hks_updown_k[atom1_orbits_up,atom1_orbits_down];
    H_k_atom1[2,1] =  Hks_updown_k[atom1_orbits_down,atom1_orbits_up];

    H_kq_atom2[1,1] =  Hks_updown_kq[atom2_orbits_up,atom2_orbits_up];
    H_kq_atom2[2,2] =  Hks_updown_kq[atom2_orbits_down,atom2_orbits_down];
    H_kq_atom2[1,2] =  Hks_updown_kq[atom2_orbits_up,atom2_orbits_down];
    H_kq_atom2[2,1] =  Hks_updown_kq[atom2_orbits_down,atom2_orbits_up];

    V_k_atom1_x = Array{Array{Complex_my,2},2}(2,2);
    V_k_atom1_x[1,1] = H_k_atom1[1,2] - H_k_atom1[2,1];
    V_k_atom1_x[2,2] = -(H_k_atom1[1,2] - H_k_atom1[2,1]);
    V_k_atom1_x[1,2] = H_k_atom1[1,1] - H_k_atom1[2,2];
    V_k_atom1_x[2,1] = -(H_k_atom1[1,1] - H_k_atom1[2,2]);

    V_kq_atom2_x = Array{Array{Complex_my,2},2}(2,2);
    V_kq_atom2_x[1,1] = H_kq_atom2[1,2] - H_kq_atom2[2,1];
    V_kq_atom2_x[2,2] = -(H_kq_atom2[1,2] - H_kq_atom2[2,1]);
    V_kq_atom2_x[1,2] = H_kq_atom2[1,1] - H_kq_atom2[2,2];
    V_kq_atom2_x[2,1] = -(H_kq_atom2[1,1] - H_kq_atom2[2,2]);

    V_k_atom1_y = Array{Array{Complex_my,2},2}(2,2);
    V_k_atom1_y[1,1] = im * (H_k_atom1[1,2] + H_k_atom1[2,1]);
    V_k_atom1_y[2,2] = -im * (H_k_atom1[1,2] + H_k_atom1[2,1]);
    V_k_atom1_y[1,2] = -im*(H_k_atom1[1,1] - H_k_atom1[2,2]);
    V_k_atom1_y[2,1] = -im*(H_k_atom1[1,1] - H_k_atom1[2,2]);

    V_kq_atom2_y = Array{Array{Complex_my,2},2}(2,2);
    V_kq_atom2_y[1,1] = im * (H_kq_atom2[1,2] + H_kq_atom2[2,1]);
    V_kq_atom2_y[2,2] = -im * (H_kq_atom2[1,2] + H_kq_atom2[2,1]);
    V_kq_atom2_y[1,2] = -im* (H_kq_atom2[1,1] - H_kq_atom2[2,2]);
    V_kq_atom2_y[2,1] = -im* (H_kq_atom2[1,1] - H_kq_atom2[2,2]);

    V_k_atom1_z = Array{Array{Complex_my,2},2}(2,2);
    V_k_atom1_z[1,1] = zeros(Complex_my,orbitalNums[atom1],orbitalNums[atom1])
    V_k_atom1_z[2,2] = zeros(Complex_my,orbitalNums[atom1],orbitalNums[atom1])
    V_k_atom1_z[1,2] = -2* H_k_atom1[1,2];
    V_k_atom1_z[2,1] = 2* H_k_atom1[2,1];

    V_kq_atom2_z = Array{Array{Complex_my,2},2}(2,2);
    V_kq_atom2_z[1,1] = zeros(Complex_my,orbitalNums[atom2],orbitalNums[atom2])
    V_kq_atom2_z[2,2] = zeros(Complex_my,orbitalNums[atom2],orbitalNums[atom2])
    V_kq_atom2_z[1,2] = -2* H_kq_atom2[1,2];
    V_kq_atom2_z[2,1] = 2* H_kq_atom2[2,1];


    Es_mkq_Vk_Es_nk_atom1_x  =  (Es_n_k_atom1_structual' * V_k_atom1_x * Es_m_kq_atom1_structual )[1,1];
    Es_nk_Vkq_Es_mkq_atom2_x =  (Es_m_kq_atom2_structual' * V_kq_atom2_x * Es_n_k_atom2_structual)[1,1];

    Es_mkq_Vk_Es_nk_atom1_y  =  (Es_n_k_atom1_structual' * V_k_atom1_y * Es_m_kq_atom1_structual )[1,1];
    Es_nk_Vkq_Es_mkq_atom2_y =  (Es_m_kq_atom2_structual' * V_kq_atom2_y * Es_n_k_atom2_structual)[1,1];

    Es_mkq_Vk_Es_nk_atom1_z  =  (Es_n_k_atom1_structual' * V_k_atom1_z * Es_m_kq_atom1_structual )[1,1];
    Es_nk_Vkq_Es_mkq_atom2_z =  (Es_m_kq_atom2_structual' * V_kq_atom2_z * Es_n_k_atom2_structual)[1,1];

# zz Test
    # col
    Es_mkq_Vk_Es_nk_atom1_z_test_col = Es_n_k_atom1_structual[1,1]'  * H_k_atom1[1,2] * Es_m_kq_atom1_structual[2,1];
    Es_nk_Vkq_Es_mkq_atom2_z_test_col = Es_m_kq_atom2_structual[2,1]'  * V_kq_atom2_z[1,2] * Es_n_k_atom2_structual[1,1];
    # non-col
    Es_mkq_Vk_Es_nk_atom1_z_test_noncol = Es_n_k_atom1_structual[1,1]'  * H_k_atom1[1,2] * Es_m_kq_atom1_structual[2,1];
    Es_nk_Vkq_Es_mkq_atom2_z_test_noncol = Es_m_kq_atom2_structual[1,1]'  * V_kq_atom2_z[1,2] * Es_n_k_atom2_structual[2,1];

#=
    G1V1_all = (Es_n_k_atom1[atom1_orbits_updown,:]'  *Hks_updown_k_atom1_diff* Es_m_kq_atom1[atom1_orbits_updown,:]);
    G2V2_all = (Es_m_kq_atom2[atom2_orbits_updown,:]'   *Hks_updown_kq_atom2_diff.' * Es_n_k_atom2[atom2_orbits_updown,:]);
    =#
    #G2V2_z = (Es_m_kq[atom2_orbits_up,:]' * Vz_2 * Es_n_k[atom2_orbits_down,:]);

    # Index convention: J_ij_(xyz,xyz)[nk,mkq]
    J_ij_xx =  -1/8.0 * 0.5* sum(part1_nk_mkq.* Es_mkq_Vk_Es_nk_atom1_x.' .* Es_nk_Vkq_Es_mkq_atom2_x );
    J_ij_xy =  -1/8.0 * 0.5* sum(part1_nk_mkq.* Es_mkq_Vk_Es_nk_atom1_x.' .* Es_nk_Vkq_Es_mkq_atom2_y );
    J_ij_xz =  -1/8.0 * 0.5* sum(part1_nk_mkq.* Es_mkq_Vk_Es_nk_atom1_x.' .* Es_nk_Vkq_Es_mkq_atom2_z );

    J_ij_yx =  -1/8.0 * 0.5* sum(part1_nk_mkq.* Es_mkq_Vk_Es_nk_atom1_y.' .* Es_nk_Vkq_Es_mkq_atom2_x );
    J_ij_yy =  -1/8.0 * 0.5* sum(part1_nk_mkq.* Es_mkq_Vk_Es_nk_atom1_y.' .* Es_nk_Vkq_Es_mkq_atom2_y );
    J_ij_yz =  -1/8.0 * 0.5* sum(part1_nk_mkq.* Es_mkq_Vk_Es_nk_atom1_y.' .* Es_nk_Vkq_Es_mkq_atom2_z );

    J_ij_zx =  -1/8.0 * 0.5* sum(part1_nk_mkq.* Es_mkq_Vk_Es_nk_atom1_z.' .* Es_nk_Vkq_Es_mkq_atom2_x );
    J_ij_zy =  -1/8.0 * 0.5* sum(part1_nk_mkq.* Es_mkq_Vk_Es_nk_atom1_z.' .* Es_nk_Vkq_Es_mkq_atom2_y );
    J_ij_zz =  -1/8.0 * 0.5* sum(part1_nk_mkq.* Es_mkq_Vk_Es_nk_atom1_z.' .* Es_nk_Vkq_Es_mkq_atom2_z );


    J_ij_zz_col =  1/2.0 * 0.5* sum(part1_nk_mkq.* Es_mkq_Vk_Es_nk_atom1_z_test_col.' .* Es_nk_Vkq_Es_mkq_atom2_z_test_col );
    J_ij_zz_noncol =  1/2.0 * 0.5* sum(part1_nk_mkq.* Es_mkq_Vk_Es_nk_atom1_z_test_noncol.' .* Es_nk_Vkq_Es_mkq_atom2_z_test_noncol );



    X_ij = sum(part1_nk_mkq.*
    (copy(Es_n_k_atom1_structual') * Es_m_kq_atom1_structual)[1,1] .*
    transpose(copy(Es_m_kq_atom2_structual')  * Es_n_k_atom2_structual)[1,1] );


    X_ij_charge_updown = sum(part1_nk_mkq.*
    (Es_n_k[atom1_orbits_up,:]' * Es_m_kq[atom1_orbits_up,:]) .*
    transpose(Es_m_kq[atom2_orbits_down,:]' * Es_n_k[atom2_orbits_down,:]) );

    X_ij_charge_downup = sum(part1_nk_mkq.*
    (Es_n_k[atom1_orbits_down,:]' * Es_m_kq[atom1_orbits_down,:]) .*
    transpose(Es_m_kq[atom2_orbits_up,:]' * Es_n_k[atom2_orbits_up,:]) );



    result_subset = [
    J_ij_xx,J_ij_xy,J_ij_xz,
    J_ij_yx,J_ij_yy,J_ij_yz,
    J_ij_zx,J_ij_zy,J_ij_zz,
    J_ij_zz_col, J_ij_zz_noncol,
    X_ij, X_ij_charge_updown, X_ij_charge_downup];
    result_mat[:,atom12_i] = result_subset
    #=
    println(atom12)
    println(  result_mat[:,atom12_i])
    MAT.matwrite(string("/home/users1/bluehope/Dropbox/shared/DFT-forge/source/debug/",atom12_i,".mat")
      ,Dict("Vz_1" =>Vz_1
      ,"Vy_1" => Vy_1
      ,"Vx_1" => Vx_1
      ,"Voff_1" => Voff_1
      ,"Hks_updown" => Hks_updown
      ,"atom1_orbits_up" => collect(atom1_orbits_up)
      ,"atom1_orbits_down" => collect(atom1_orbits_down)
      ,"result_subset" => result_subset
      #,"k_point" => k_point
      #,"kq_point" => kq_point
      )
      );
      =#
  end

  return result_mat #sum(J_ij[:]);
end



#
theta_list = 0.0:0.1:1
phi_list = 0.0:0.2:2
#global theta_phi_num

theta_phi_num = 0;
theta_phi_list = Array{Array{Float64},1}()
for theta in theta_list
  for phi in phi_list
    theta_phi_num += 1;
    push!(theta_phi_list,[theta,phi])
  end
end

#

  ## 4.1 Do K,Q sum
  # for orbital_mask1_list,orbital_mask2_list combinations



#println(DFTcommon.bar_string) # print ====...====
#println(orbital_mask_name," mask1 ",orbital_mask1,"\tmask2 ",orbital_mask2)
#=
import Plots
Plots.plotly()
Plots.plot(rand(30,3))
Plots.savefig("test.html")
=#
# #=
# setup extra info
#=
orbital_mask_input = orbital_mask_input_Type(orbital_mask1_list[1],orbital_mask1_list[1],(-1,-1),false)
DFTforge.pwork(init_orbital_mask,orbital_mask_input)
DFTforge.pwork(init_variables,(ChemP_delta_ev, Hmode, theta_phi_list))
k_point = (0.0,0.0,0.0)
kq_point = (0.0,0.0,0.0)
spin_type = DFTcommon.non_colinear_type
atom12 = atom12_list
input = Job_input_kq_atom_list_Type(k_point,kq_point,spin_type,atom12)
import Plots
Plots.plotly()
@everywhere import Plots
@everywhere Plots.plotly()

=#

#=


global num_return = num_return;
=#
#num_return = 10
@everywhere function Magnetic_Exchange_J_noncolinear_v2(input::Job_input_kq_atom_list_Type)
  global Hks_0;
  global Hmode;
  #global SmallHks;
  #global Hks_0;
  ############################################################################
  ## Accessing Data Start
  ############################################################################
  # Common input info
  #num_return = 10;
  global num_return;
  global theta_phi_list;
  num_return = 10
  k_point =  input.k_point::DFTforge.k_point_Tuple
  kq_point =  input.kq_point::DFTforge.k_point_Tuple
  spin_type = input.spin_type::DFTforge.SPINtype;

  #atom1::Int = input.atom12[1][1];
  #atom2::Int = input.atom12[1][2];
  atom12_list = input.atom12_list::Vector{Tuple{Int64,Int64}};
  result_mat = zeros(Complex_my,num_return,length(theta_phi_list),length(atom12_list));

  #atom1::Int = input.atom1;
  #atom2::Int = input.atom2;

  result_index =  input.result_index::Int;
  cache_index = input.cache_index::Int;

  # Get Chemp, E_temp
  TotalOrbitalNum = cacheread(cache_index).TotalOrbitalNum::Int;
  TotalOrbitalNum2 = TotalOrbitalNum;
  if (DFTcommon.non_colinear_type == spin_type)
    TotalOrbitalNum2 = 2*TotalOrbitalNum;
  end

  ChemP = get_dftdataset(result_index).scf_r.ChemP;
  E_temp = get_dftdataset(result_index).scf_r.E_Temp;

  # Get EigenState
  #eigenstate_k_up::Kpoint_eigenstate  = cacheread_eigenstate(k_point,1,cache_index)
  eigenstate_kq = cacheread_eigenstate(kq_point,1,cache_index);
  eigenstate_k  = cacheread_eigenstate(k_point,1,cache_index);
  #eigenstate_kq_down::Kpoint_eigenstate = cacheread_eigenstate(kq_point,2,cache_index)

  # Get Hamiltonian
  #Hks_updown = cacheread_Hamiltonian(1,cache_index)
  #Hks_updown = cal_Hamiltonian(1,1)
  #Hks_updown = copy(Hks_0);

  #Hks_updown = cacheread_Hamiltonian(1,cache_index)
  #Hks_updown = copy(Hks_0);
  #@assert( sum(real(Hks_updown - Hks_down2)) < 10.0^-4 )
  Hks_updown_k  = cacheread_Hamiltonian(k_point, Hmode,cache_index)::Hamiltonian_type;
  Hks_updown_kq  = cacheread_Hamiltonian(kq_point, Hmode,cache_index)::Hamiltonian_type;
  @assert(size(Hks_updown_k)[1] == TotalOrbitalNum2);
  #Hks_k_down::Hamiltonian_type = cacheread_Hamiltonian(k_point,2,cache_index)
  #Hks_kq_down::Hamiltonian_type = cacheread_Hamiltonian(kq_point,2,cache_index)

  (orbitalStartIdx_list,orbitalNums) = cacheread_atomsOrbital_lists(cache_index);

  #En_k_up::Array{Float_my,1} = eigenstate_k_up.Eigenvalues;
  Em_kq = eigenstate_kq.Eigenvalues;
  En_k  = eigenstate_k.Eigenvalues;
  #Em_kq_down::Array{Float_my,1} = eigenstate_kq_down.Eigenvalues;

  #Es_n_k_up::Array{Complex_my,2} = eigenstate_k_up.Eigenstate;
  Es_m_kq  = eigenstate_kq.Eigenstate;
  Es_n_k  = eigenstate_k.Eigenstate;

  #for mask2 in orbital_mask2
  #    Es_m_kq_atom2[mask2,:]=0.0
  #end
  #Plots.heatmap(real(Es_n_k))
  #Es_m_kq_down::Array{Complex_my,2} = eigenstate_kq_down.Eigenstate;
  Fftn_k  = 1.0./(exp.( ((En_k)  - ChemP)/(kBeV*E_temp)) + 1.0 );
  Fftm_kq = 1.0./(exp.( ((Em_kq) - ChemP)/(kBeV*E_temp)) + 1.0 );
  # Index convention: dFnk[nk,mkq]
  dFnk_Fmkq  =
    Fftn_k*ones(1,TotalOrbitalNum2)  - ones(TotalOrbitalNum2,1)*Fftm_kq' ;
  # Index convention: Enk_Emkq[nk,mkq]
  Enk_Emkq =
    En_k[:]*ones(1,TotalOrbitalNum2) - ones(TotalOrbitalNum2,1)*Em_kq[:]' ;
  Enk_Emkq += im*0.0001;

  # Index convention: part1[nk,mkq]
  part1 = dFnk_Fmkq./(-Enk_Emkq);

  for (atom12_i,atom12) in enumerate(atom12_list)
  #  atom12_i =1;
    atom12 = atom12_list[atom12_i]

    atom1 = atom12[1]
    atom2 = atom12[2]

    Es_m_kq_atom1 = copy(Es_m_kq);
    Es_m_kq_atom2 = copy(Es_m_kq);
    Es_n_k_atom1 = copy(Es_n_k);
    Es_n_k_atom2 = copy(Es_n_k);
    if (length(orbital_mask1)>0)
      orbital_mask1_tmp = collect(1:orbitalNums[atom1]);
      for orbit1 in orbital_mask1
          deleteat!(orbital_mask1_tmp, find(orbital_mask1_tmp.==orbit1))
      end
      Es_n_k_atom1[orbitalStartIdx_list[atom1]+orbital_mask1_tmp,:]=0.0; # up
      Es_n_k_atom1[TotalOrbitalNum + orbitalStartIdx_list[atom1]+orbital_mask1_tmp,:]=0.0; # down
    end
    if (length(orbital_mask2)>0)
      orbital_mask2_tmp = collect(1:orbitalNums[atom2]);
      for orbit2 in orbital_mask2
          deleteat!(orbital_mask2_tmp, find(orbital_mask2_tmp.==orbit2))
      end
      Es_m_kq_atom2[orbitalStartIdx_list[atom2]+orbital_mask2_tmp,:]=0.0; # up
      Es_m_kq_atom2[TotalOrbitalNum + orbitalStartIdx_list[atom2]+orbital_mask2_tmp,:]=0.0; # down
    end



    atom1_orbits_up   = orbitalStartIdx_list[atom1] + (1:orbitalNums[atom1])
    atom1_orbits_down = TotalOrbitalNum + atom1_orbits_up
    atom2_orbits_up   = orbitalStartIdx_list[atom2] + (1:orbitalNums[atom2])
    atom2_orbits_down = TotalOrbitalNum + atom2_orbits_up

    atom1_orbits_updown = [atom1_orbits_up;atom1_orbits_down];
    atom2_orbits_updown = [atom2_orbits_up;atom2_orbits_down];
    atom1_orbits_up_rel = (1:orbitalNums[atom1]);
    atom1_orbits_down_rel = orbitalNums[atom1]+(1:orbitalNums[atom1]);
    atom2_orbits_up_rel = (1:orbitalNums[atom2]);
    atom2_orbits_down_rel = orbitalNums[atom2]+(1:orbitalNums[atom2]);
    #Plots.heatmap(real(Hks_updown))

    ############################################################################
    ## Accessing Data End
    ############################################################################
    ## Do auctual calucations

    #=
    Vz_1 = 0.5*(Hks_updown_k[atom1_orbits_up,atom1_orbits_up] - Hks_updown_k[atom1_orbits_down,atom1_orbits_down]);
    #Plots.heatmap(real(Vz_1));Plots.savefig("fig/Vz_1_r.html");
    #Voff_1 = Hks_updown[atom1_orbits_up,atom1_orbits_down]
    Voff_1 = conj(Hks_updown_k[atom1_orbits_up,atom1_orbits_down]) + Hks_updown_k[atom1_orbits_down,atom1_orbits_up] ;
    Vx_1 = 0.5*(Voff_1 + conj(Voff_1))/2.0;
    Vy_1 = 0.5*(Voff_1 - conj(Voff_1))/(2.0*im);
    #Plots.heatmap(real(Vz_1))
    #Plots.heatmap(real(Vy_1))


    Vz_2 = 0.5*(Hks_updown_kq[atom2_orbits_up,atom2_orbits_up] - Hks_updown_kq[atom2_orbits_down,atom2_orbits_down]) ;
    #Plots.heatmap(real(Vz_2));Plots.savefig("fig/Vz_2_r.html");
    Voff_2 = conj(Hks_updown_kq[atom2_orbits_up,atom2_orbits_down]) + Hks_updown_kq[atom2_orbits_down,atom2_orbits_up];
    # Voff_2 = (Hks_updown[atom2_orbits_up,atom2_orbits_down])
           #+ Hks_updown[atom2_orbits_down,atom2_orbits_up]

    Vx_2 = 0.5*(Voff_2 + conj(Voff_2))/2.0;
    Vy_2 = 0.5*(Voff_2 - conj(Voff_2))/(2.0*im);
    =#

    #theta_list = 0.0:0.05:1.0
    #phi_list = 0.0:0.05:2.0
    #global theta_list# = 0.0:0.1:1
    #global phi_list# = 0.0:0.1:2

    results_tmp = zeros(ComplexF64,length(theta_phi_list),num_return);
    # #=
    for (key_index,v) in enumerate(theta_phi_list)
      theta = v[1];
      phi = v[2]
    # =#
      #key_index = 1
      #theta = 0.0;
      #phi = 0.1;
      #println(k," ",theta," ",phi)
      if true
        Hks_updown_k_atom1 = copy(Hks_updown_k[atom1_orbits_updown,atom1_orbits_updown]);
        Hks_updown_k_atom2 = copy(Hks_updown_k[atom2_orbits_updown,atom2_orbits_updown]);
        Rot_theta_phi_1 = RotPauliThetaPhi(pi*theta, pi*phi)
        #Rot_ZXZ_2 = RotPauliZXZ(-pi*.2,-pi*0.3,-pi*.15)

        Rot_theta_phi_2 = inv(Rot_theta_phi_1)
        Rot_theta_phi_2 = copy(Rot_theta_phi_1)
        #Rot_theta_phi_2 = RotPauliThetaPhi(-pi*0.25,-pi*0.0);
        @assert(1.0 - norm(Rot_theta_phi_1) < 0.001);
        #=
        Rot_ZYZ = RotPauliZ(pi*1.0)*RotPauliY2(pi*0.5)*RotPauliZ(pi*0.2)
        norm(Rot_ZYZ)
        Rot_ZXZ = RotPauliZ(pi*1.0)*RotPauliX(pi*0.5)*RotPauliZ(pi*1.0)
        norm(Rot_ZXZ)
        =#
        Hks_updown_k_atom1 = ncHamiltonian_rot!(Hks_updown_k_atom1, atom1_orbits_up_rel,Rot_theta_phi_1);
        Hks_updown_k_atom1_diff = Hks_updown_k_atom1 - Hks_updown_k[atom1_orbits_updown,atom1_orbits_updown];
        Hks_updown_k_atom1_diff *= 0.5;
        sum(abs(Hks_updown_k_atom1_diff[:]))
        Hks_updown_k_atom2 = ncHamiltonian_rot!(Hks_updown_k_atom2, atom2_orbits_up_rel,Rot_theta_phi_2);
        Hks_updown_k_atom2_diff = Hks_updown_k_atom2 - Hks_updown_k[atom2_orbits_updown,atom2_orbits_updown];
        Hks_updown_k_atom2_diff *= 0.5;
        sum(abs(Hks_updown_k_atom2_diff[:]))


        #=
        Plots.heatmap(real(Hks_updown_k[atom1_orbits_updown,atom1_orbits_updown]));Plots.savefig("fig/Hks_updown_k_r.html")
        Plots.heatmap(real(Hks_updown_k[atom1_orbits_up,atom1_orbits_down]));Plots.savefig("fig/Hks_updown_k_off_r.html")
        Plots.heatmap(real(Hks_updown_k_atom1));Plots.savefig("fig/Hks_updown_k_atom1_r.html")
        Plots.heatmap(real(Hks_updown_k_atom1_diff)/2);Plots.savefig("fig/Hks_updown_k_atom1_diff_r.html")
        Plots.heatmap(real(Hks_updown_k_atom2_diff)/2);Plots.savefig("fig/Hks_updown_k_atom2_diff_r.html")

        =#
        Vz_1 = 0.5*(Hks_updown_k_atom1_diff[atom1_orbits_up_rel,atom1_orbits_up_rel] - Hks_updown_k_atom1_diff[atom1_orbits_down_rel,atom1_orbits_down_rel])/2.0;
        #Voff_1 = conj(Hks_updown_k_atom1_diff[atom1_orbits_up_rel,atom1_orbits_down_rel]) + Hks_updown_k_atom1_diff[atom1_orbits_down_rel,atom1_orbits_up_rel];

        #Voff_1 = conj(Hks_updown_k_atom1_diff[atom1_orbits_up_rel,atom1_orbits_down_rel]);

        Vx_1 = 0.5*((Hks_updown_k_atom1_diff[atom1_orbits_up_rel,atom1_orbits_down_rel]) + Hks_updown_k_atom1_diff[atom1_orbits_down_rel,atom1_orbits_up_rel])/2.0;
        Vy_1 = 0.5*(-(Hks_updown_k_atom1_diff[atom1_orbits_up_rel,atom1_orbits_down_rel]) + Hks_updown_k_atom1_diff[atom1_orbits_down_rel,atom1_orbits_up_rel])/(2.0*im);
        #=
        Plots.heatmap(real(Vx_1));Plots.savefig("fig/Vx_1.html")
        Plots.heatmap(real(Vy_1));Plots.savefig("fig/Vy_1.html")
        Plots.heatmap(real(Vz_1));Plots.savefig("fig/Vz_1.html")
        =#
        Vz_2 = 0.5*(Hks_updown_k_atom2_diff[atom2_orbits_up_rel,atom2_orbits_up_rel] - Hks_updown_k_atom2_diff[atom2_orbits_down_rel,atom2_orbits_down_rel])/2.0;
        #Voff_2 = conj(Hks_updown_k_atom2_diff[atom2_orbits_up_rel,atom2_orbits_down_rel]) + Hks_updown_k_atom2_diff[atom2_orbits_down_rel,atom2_orbits_up_rel];

        #Voff_2 = conj(Hks_updown_k_atom2_diff[atom2_orbits_up_rel,atom2_orbits_down_rel])

        Vx_2 = 0.5*((Hks_updown_k_atom2_diff[atom2_orbits_up_rel,atom2_orbits_down_rel]) + Hks_updown_k_atom2_diff[atom2_orbits_down_rel,atom2_orbits_up_rel])/2.0;
        Vy_2 = 0.5*(-(Hks_updown_k_atom2_diff[atom2_orbits_up_rel,atom2_orbits_down_rel]) + Hks_updown_k_atom2_diff[atom2_orbits_down_rel,atom2_orbits_up_rel])/(2.0*im);
        #=
        Plots.heatmap(real(Vx_2));Plots.savefig("fig/Vx_2.html")
        Plots.heatmap(real(Vy_2));Plots.savefig("fig/Vy_2.html")
        Plots.heatmap(real(Vz_2));Plots.savefig("fig/Vz_2.html")
        =#
        #Plots.heatmap(real(part1))

        G1V1_z = (Es_n_k_atom1[atom1_orbits_down,:]'  *Vz_1* Es_m_kq_atom1[atom1_orbits_up,:]);
        #Plots.heatmap(real(G1V1_z));Plots.savefig("fig/G1V1_z_r.html")
        G2V2_z = (Es_m_kq_atom2[atom2_orbits_up,:]'   *Vz_2* Es_n_k_atom2[atom2_orbits_down,:]);
        #G2V2_z = (Es_m_kq[atom2_orbits_up,:]' * Vz_2 * Es_n_k[atom2_orbits_down,:]);
        G1V1_x = (Es_n_k_atom1[atom1_orbits_down,:]'  *Vx_1* Es_m_kq_atom1[atom1_orbits_up,:]);
        G2V2_x = (Es_m_kq_atom2[atom2_orbits_up,:]'   *Vx_2* Es_n_k_atom2[atom2_orbits_down,:]);

        G1V1_y = (Es_n_k_atom1[atom1_orbits_down,:]'  *Vy_1* Es_m_kq_atom1[atom1_orbits_up,:]);
        G2V2_y = (Es_m_kq_atom2[atom2_orbits_up,:]'   *Vy_2* Es_n_k_atom2[atom2_orbits_down,:]);

        # Index convention: J_ij_(xyz,xyz)[nk,mkq]
        J_ij_xx =  0.5* sum(part1.* G1V1_x .* transpose(G2V2_x) );
        J_ij_xy =  0.5* sum(part1.* G1V1_x .* transpose(G2V2_y) );
        J_ij_xz =  0.5* sum(part1.* G1V1_x .* transpose(G2V2_z) );

        J_ij_yx =  0.5* sum(part1.* G1V1_y .* transpose(G2V2_x) );
        J_ij_yy =  0.5* sum(part1.* G1V1_y .* transpose(G2V2_y) );
        J_ij_yz =  0.5* sum(part1.* G1V1_y .* transpose(G2V2_z) );

        J_ij_zx =  0.5* sum(part1.* G1V1_z .* transpose(G2V2_x) );
        J_ij_zy =  0.5* sum(part1.* G1V1_z .* transpose(G2V2_y) );
        J_ij_zz =  0.5* sum(part1.* G1V1_z .* transpose(G2V2_z) )

        X_ij = sum(part1.* (Es_n_k[atom1_orbits_down,:]' * Es_m_kq[atom1_orbits_up,:]) .* transpose(Es_m_kq[atom2_orbits_up,:]' * Es_n_k[atom2_orbits_down,:]) );
        result_subset = [J_ij_xx,J_ij_xy,J_ij_xz,  J_ij_yx,J_ij_yy,J_ij_yz,  J_ij_zx,J_ij_zy,J_ij_zz, X_ij];
        result_mat[:,key_index,atom12_i] = result_subset
        #

        # #=
        # =#
      end
    end
      #result_subset*1000.0
      #println(sum(result_subset[1:9]*1000.0))
  end


  return result_mat #sum(J_ij[:]);
end




## 4.1 Do K,Q sum
# for orbital_mask1_list,orbital_mask2_list combinations
for (orbital1_i,orbital_mask1) in enumerate(orbital_mask1_list)
  for (orbital2_i,orbital_mask2) in enumerate(orbital_mask2_list)

    for (orbital3_i,orbital_mask3) in enumerate(orbital_mask3_list)
      for (orbital4_i,orbital_mask4) in enumerate(orbital_mask4_list)
        orbital_mask_input = orbital_mask_input_Type(orbital_mask1,orbital_mask2,
          orbital_mask3,orbital_mask4,(-1,-1),false)
        if (orbital_mask_on)
            orbital_mask_input = orbital_mask_input_Type(orbital_mask1,orbital_mask2,
            orbital_mask3,orbital_mask4,(-1,-1),true)
        end
        orbital_mask_name = orbital_mask1_names[orbital1_i]*"_"*orbital_mask2_names[orbital2_i]*"__"*
          orbital_mask3_names[orbital3_i]*"__"*orbital_mask4_names[orbital4_i];
        println(DFTcommon.bar_string) # print ====...====
        println(orbital_mask_name," mask1 ",orbital_mask1,"\tmask2 ",orbital_mask2,
        "\tmask3 ",orbital_mask3,"\tmask4 ",orbital_mask4)

        # setup extra info
        DFTforge.pwork(init_orbital_mask,orbital_mask_input)
        DFTforge.pwork(init_variables,(ChemP_delta_ev,Hmode,theta_phi_list))

        #(X_Q_nc,X_Q_mean_nc) = Qspace_Ksum_atomlist_parallel_nc(Magnetic_Exchange_J_noncolinear_v2,
        #q_point_list,k_point_list,atom12_list,theta_phi_list,num_return)
        (X_Q_nc,X_Q_mean_nc) = Qspace_Ksum_atomlist_parallel(Magnetic_Exchange_J_noncolinear_fixed_v1,
        q_point_list,k_point_list,atom12_list,num_return)
        #println(typeof(X_Q_mean_nc))
        println(DFTcommon.bar_string) # print ====...====
        println("[1:Jxx, 2:Jxy, 3:Jyy]")
        println("[4:Jyx, 5:Jyy, 6:Jyz]")
        println("[7:Jzx, 8:Jzy, 9:Jzz] 10:X")
        ## 4.2 reduce K,Q to Q space
        # Average X_Q results

        #= v2
        Xij_Q_mean_matlab = Array(Array{Complex_my,1},num_return,length(theta_phi_list),length(atom12_list));
        for (atom12_i,atom12) in enumerate(atom12_list)
          atom1 = atom12[1];
          atom2 = atom12[2];
          for xyz_i = 1:num_return
            for theta_phi_i = 1:length(theta_phi_list)
              Xij_Q_mean_matlab[xyz_i,theta_phi_i,atom12_i] = zeros(Complex_my,length(q_point_list));
              for (q_i,q_point) in enumerate(q_point_list)
                q_point_int = k_point_float2int(q_point);
                Xij_Q_mean_matlab[xyz_i,theta_phi_i,atom12_i][q_i] =
                  X_Q_mean_nc[xyz_i,theta_phi_i,atom12_i][q_point_int];
              end
            end
            println(string(" Gamma point J [",atom1,",",atom2,"] ",xyz_i," :",
              1000.0*mean(Xij_Q_mean_matlab[xyz_i,:,atom12_i][:])," meV"))
          end
        end
        =#
		Xij_Q_mean_matlab = Array{Array{Complex_my,1}}(num_return,length(atom12_list));
        for (atom12_i,atom12) in enumerate(atom12_list)
          atom1 = atom12[1];
          atom2 = atom12[2];
          for xyz_i = 1:num_return

            Xij_Q_mean_matlab[xyz_i,atom12_i] = zeros(Complex_my,length(q_point_list));
            for (q_i,q_point) in enumerate(q_point_list)
              q_point_int = k_point_float2int(q_point);
              Xij_Q_mean_matlab[xyz_i,atom12_i][q_i] =
                X_Q_mean_nc[xyz_i,atom12_i][q_point_int];
            end

            println(string(" Gamma point J [",atom1,",",atom2,"] ",xyz_i," :",
              1000.0*mean(Xij_Q_mean_matlab[xyz_i,atom12_i][:])," meV"))
          end
        end
        ###############################################################################
        ## 5. Save results and clear hdf5 file
        ## 5.1 Prepaire infomations for outout
        ## 5.2 Write to MAT
        ###############################################################################
        optionalOutputDict = Dict{AbstractString,Any}()
        optionalOutputDict["num_return"] = num_return;
        optionalOutputDict["VERSION_Spin_Exchange"] = string(X_VERSION);

        optionalOutputDict["theta_phi_list"] = theta_phi_list;
        #=
        export2mat_K_Q_nc(Xij_Q_mean_matlab,theta_phi_list,hamiltonian_info,q_point_list,k_point_list,atom12_list,
        orbital_mask_on,orbital_mask1,orbital_mask2,orbital_mask3,orbital_mask4,ChemP_delta_ev,
        optionalOutputDict,
        jq_output_dir,cal_name,
        orbital_mask_name,cal_type);
        =#
        export2mat_K_Q(Xij_Q_mean_matlab,hamiltonian_info,q_point_list,k_point_list,atom12_list,
        orbital_mask_on,orbital_mask1,orbital_mask2,orbital_mask3,orbital_mask4,ChemP_delta_ev,
        optionalOutputDict,
        jq_output_dir,cal_name,
        orbital_mask_name,cal_type);
      end
    end
  end
end
## 5.3 Cleanup HDF5 cache file
println(DFTcommon.bar_string) # print ====...====
println("hdf_cache_name:",hdf_cache_name)
if (isfile(hdf_cache_name))
  rm(hdf_cache_name)
  file = open(hdf_cache_name,"w");
  close(file);
end
