################################################################################

using ProgressMeter
import DFTforge
using DFTforge.DFTrefinery
using DFTcommon
import MAT


X_VERSION = VersionNumber("0.1.0-dev+20170401");
println("X_VERSION: ",X_VERSION)

@everywhere import DFTforge
@everywhere using DFTforge.DFTrefinery
@everywhere using DFTcommon

#scf_name = "../examples/NiO/4cell/Ni6.0_O5.0-s2p2d2f1-4.25_0.0-4.180-k10/nio.scfout"
#scf_name = "/home/users1/bluehope/work_local/GaVS/LDA_FS_rhom/GaV4S8.scfout"
#scf_name = "/home/users1/bluehope/work_local/GaVS/LDA_FS_cubic/GaV4S8.scfout"
#scf_name = "/home/users1/bluehope/work/Jx/X0/FeTe/global/FeTe_PM_root2xroot2_AF_s2p2d2f1_expstr_U0.0.LDA/feTe.scfout"
#scf_name = "/home/users1/bluehope/work/Jx/X0/SRO/1x1/SRO_ortho.scfout"

##############################################################################
hdftmpdir = ""

orbital_mask1 = Array{Int64,1}();
orbital_mask2 = Array{Int64,1}();
orbital_mask_name = "";
orbital_mask_option = DFTcommon.nomask;
orbital_mask_on = false;
atom1 = -1;
atom2 = -1;
# #=
if true
arg_input = parse_input(ARGS)
scf_name = arg_input.scf_name
ChemP_delta_ev = arg_input.ChemP_delta_ev
k_point_num = arg_input.k_point_num
q_point_num = arg_input.q_point_num
atom1 = arg_input.atom1;
atom2 = arg_input.atom2;
hdftmpdir = arg_input.hdftmpdir;
orbital_mask1 = arg_input.orbital_mask1;
orbital_mask2 = arg_input.orbital_mask2;
end
# =#

###
# TOML override (should be merged to DFTcommon.jl)
###
if false
import TOML
#toml_inputs = TOML.parse(readstring("nio_J.toml"))
toml_inputs = TOML.parse(readstring("Fe_J.toml"))
atom1 = toml_inputs["atom12"][1][1];
atom2 = toml_inputs["atom12"][1][2];
scf_name = toml_inputs["scf_fname"];
k_point_num = [3,3,3]
q_point_num = [3,3,3]
#k_point_num = [9,9,9]
#q_point_num = [9,9,9]
ChemP_delta_ev = 0.0
end

###############################################################################

println((atom1,atom2))

root_dir = dirname(scf_name)
scf_name_only = splitext(basename(scf_name))
cal_name = scf_name_only[1];

jq_output_dir =  joinpath(root_dir,string("jq_test_", ChemP_delta_ev))
if (!isdir(jq_output_dir))
  mkdir(jq_output_dir)
end
if ("" == hdftmpdir || !isdir(hdftmpdir) )
  hdftmpdir = jq_output_dir
end
hdf_cache_name = joinpath(hdftmpdir,string(cal_name,".hdf5"))
println(hdf_cache_name)

##############################################################################
scf_test = DFTforge.OpenMXdata.read_scf(scf_name);

scf_r = set_current_dftdataset(scf_name, DFTforge.OpenMX, DFTforge.colinear_type)
DFTforge.pwork(set_current_dftdataset,(scf_name, DFTforge.OpenMX, DFTforge.colinear_type,1));
#scf_r = set_current_dftdataset(scf_name, DFTforge.OpenMX, DFTforge.non_colinear_type)
#DFTforge.pwork(set_current_dftdataset,(scf_name, DFTforge.OpenMX, DFTforge.non_colinear_type,1));
##


#k_point_num = [2,2,2]
#q_point_num = [2,2,2]

k_point_list = kPoint_gen_GammaCenter(k_point_num);
q_point_list = kPoint_gen_GammaCenter(q_point_num);

#k_point_list = kPoint_gen_EquallySpaced(k_point_num);
#q_point_list = kPoint_gen_EquallySpaced(q_point_num);

#q_point_list = unique(q_point_list);
(kq_point_list,kq_point_int_list) = q_k_unique_points(q_point_list,k_point_list)
println(string(" kq_point_list ",length(kq_point_list)))
println(string(" q point ",length(q_point_list) ," k point ",length(k_point_list)))

eigenstate_cache = cachecal_all_Qpoint_eigenstats(kq_point_list,hdf_cache_name);
gc();
DFTforge.pwork(cacheset,eigenstate_cache)
tic();
DFTforge.pwork(cacheread_lampup,kq_point_list)
toc();
#DFTforge.pwork(cacheread_eigenstate,1,(0.0,0.0,0.0))
##############################################################################
##############################################################################
@everywhere function init_orbital_mask(orbital_mask_input::orbital_mask_input_Type)
    global orbital_mask1,orbital_mask2,orbital_mask_on
    orbital_mask1 = Array{Int64,1}();
    orbital_mask2 = Array{Int64,1}();
    if (orbital_mask_input.orbital_mask_on)
        orbital_mask1 = orbital_mask_input.orbital_mask1;
        orbital_mask2 = orbital_mask_input.orbital_mask2;
        orbital_mask_on = true;
    else
        orbital_mask_on = false;
    end
    #println(orbital_mask1)
end
@everywhere function init_variables(Input_ChemP_delta_ev)
  global ChemP_delta_ev
  ChemP_delta_ev = Input_ChemP_delta_ev;
end
orbital_mask_input = orbital_mask_input_Type(orbital_mask1,orbital_mask2,(-1,-1),false)
if ((DFTcommon.unmask == orbital_mask_option) || (DFTcommon.mask == orbital_mask_option) )
    orbital_mask_input = orbital_mask_input_Type(orbital_mask1,orbital_mask2,(-1,-1),true)
end
DFTforge.pwork(init_orbital_mask,orbital_mask_input)
DFTforge.pwork(init_variables,ChemP_delta_ev)
##############################################################################

@everywhere function Init_SmallHks(atom12)
    atom1 = atom12[1];
    atom2 = atom12[2];
    global SmallHks;
    result_index = 1;
    SmallHks =  DFTforge.OpenMXdata.test_SmallHks(atom1,atom2,get_dftdataset(result_index).scf_r);
end


DFTforge.pwork(Init_SmallHks,(atom1,atom2))

@everywhere function Magnetic_Exchange_J_colinear(input::Job_input_kq_atom_Type)
    global orbital_mask1,orbital_mask2,orbital_mask_on
    global ChemP_delta_ev
    #global SmallHks;
    ############################################################################
    ## Accessing Data Start
    ############################################################################
    # Common input info
    k_point::DFTforge.k_point_Tuple =  input.k_point
    kq_point::DFTforge.k_point_Tuple =  input.kq_point
    spin_type::DFTforge.SPINtype = input.spin_type;

    atom1::Int = input.atom1;
    atom2::Int = input.atom2;

    result_index::Int = input.result_index;
    cache_index::Int = input.cache_index;

    # Get Chemp, E_temp
    TotalOrbitalNum = cacheread(cache_index).TotalOrbitalNum;
    TotalOrbitalNum2 = TotalOrbitalNum;
    if (DFTforge.non_colinear_type == spin_type)
      TotalOrbitalNum2 = 2*TotalOrbitalNum;
    end

    ChemP = get_dftdataset(result_index).scf_r.ChemP;
    E_temp = get_dftdataset(result_index).scf_r.E_Temp;

    # Get EigenState
    #eigenstate_k_up::Kpoint_eigenstate  = cacheread_eigenstate(k_point,1,cache_index)
    eigenstate_kq_up::Kpoint_eigenstate = cacheread_eigenstate(kq_point,1,cache_index)
    eigenstate_k_down::Kpoint_eigenstate  = cacheread_eigenstate(k_point,2,cache_index)
    #eigenstate_kq_down::Kpoint_eigenstate = cacheread_eigenstate(kq_point,2,cache_index)

    # Get Hamiltonian
    Hks_up::Hamiltonian_type = cacheread_Hamiltonian(1,cache_index)
    Hks_down::Hamiltonian_type = cacheread_Hamiltonian(2,cache_index)
    (orbitalStartIdx,orbitalNums) = cacheread_atomsOrbital_lists(cache_index)

    atom1_orbitals = orbitalStartIdx[atom1]+(1:orbitalNums[atom1]);
    atom2_orbitals = orbitalStartIdx[atom2]+(1:orbitalNums[atom2]);
    ############################################################################
    ## Accessing Data End
    ############################################################################
    ## Do auctual calucations
    #En_k_up::Array{Float_my,1} = eigenstate_k_up.Eigenvalues;
    Em_kq_up::Array{Float_my,1} = eigenstate_kq_up.Eigenvalues;
    En_k_down::Array{Float_my,1} = eigenstate_k_down.Eigenvalues;
    #Em_kq_down::Array{Float_my,1} = eigenstate_kq_down.Eigenvalues;

    #Es_n_k_up::Array{Complex_my,2} = eigenstate_k_up.Eigenstate;
    Es_m_kq_up::Array{Complex_my,2} = eigenstate_kq_up.Eigenstate;
    Es_n_k_down::Array{Complex_my,2} = eigenstate_k_down.Eigenstate;
    #Es_m_kq_down::Array{Complex_my,2} = eigenstate_kq_down.Eigenstate;

    Es_n_k_down_atom1 = copy(Es_n_k_down);
    Es_m_kq_up_atom1  = copy(Es_m_kq_up);
    Es_m_kq_up_atom2  = copy(Es_m_kq_up);
    Es_n_k_down_atom2 = copy(Es_n_k_down);
    for mask1 in orbital_mask1
        Es_n_k_down_atom1[mask1,:]=0.0;
    end
    for mask2 in orbital_mask2
        Es_m_kq_up_atom2[mask2,:]=0.0
    end

    Fftn_k::Array{Float_my,1}  = 1.0./(exp( ((En_k_down)  - ChemP)/(kB*E_temp)) + 1.0 );
    Fftm_kq::Array{Float_my,1} = 1.0./(exp( ((Em_kq_up) - ChemP)/(kB*E_temp)) + 1.0 );

    dFnk_Fmkq::Array{Float_my,2} =
      Fftn_k*ones(1,TotalOrbitalNum2)  - ones(TotalOrbitalNum2,1)*Fftm_kq[:]' ;

    Enk_Emkq::Array{Complex_my,2} =
      En_k_down[:]*ones(1,TotalOrbitalNum2) - ones(TotalOrbitalNum2,1)*Em_kq_up[:]' ;
    #Enk_Emkq += im*0.001;

    V1 = 0.5 *(Hks_up[atom1_orbitals,atom1_orbitals] - Hks_down[atom1_orbitals,atom1_orbitals])
    V2 = 0.5 *(Hks_up[atom2_orbitals,atom2_orbitals] - Hks_down[atom2_orbitals,atom2_orbitals])

    atom1_orbitals_rel = 1:orbitalNums[atom1];
    atom2_orbitals_rel = 1:orbitalNums[atom2];
    atom2_orbitals_rel2 = orbitalNums[atom1]+atom2_orbitals_rel;

    V1_2 = 0.5 *( SmallHks[1][atom1_orbitals_rel,atom1_orbitals_rel] -
      SmallHks[2][atom1_orbitals_rel,atom1_orbitals_rel] );
    V2_2 = 0.5 *( SmallHks[1][atom2_orbitals_rel2,atom2_orbitals_rel2] -
        SmallHks[2][atom2_orbitals_rel2,atom2_orbitals_rel2] );

    VV1::Array{Complex_my,2} = (Es_n_k_down_atom1[atom1_orbitals,:]' * V1 * Es_m_kq_up_atom1[atom1_orbitals,:]);
    VV2::Array{Complex_my,2} = (Es_m_kq_up_atom2[atom2_orbitals,:]' * V2 * Es_n_k_down_atom2[atom2_orbitals,:]);
    Vi_Vj = transpose(VV1).*VV2;


    J_ij::Array{Complex_my,2} =  0.5./(-Enk_Emkq).*dFnk_Fmkq .* Vi_Vj ;
    #return sum(J_ij[:])*Hartree2cm;
    return sum(J_ij[!isnan(J_ij)] )*Hartree2cm;
  end

##############################
#DEBUG
##############################
##############################
#DEBUG END
##############################


###############################################################################
# Average X_Q results
###############################################################################
#X_Q::Dict{k_point_int_Tuple,Array{Complex_my,1}} = X_succeptability();
X_Q = Qspace_Ksum_atom_parallel(Magnetic_Exchange_J_colinear,q_point_list,k_point_list,atom1,atom2)
println(typeof(X_Q))
X_Q_matlab = zeros(Complex128,length(q_point_list),1);
for (q_i,q_point) in enumerate(q_point_list)
  q_point_int = k_point_float2int(q_point);
  X_Q_matlab[q_i] = mean(X_Q[q_point_int])
end
println(string(" X_Q sum ",sum(X_Q_matlab)*cm2meV));
###############################################################################
# Prepaire infomations for outout
###############################################################################
q_point_int_list = Array{k_point_int_Tuple,1}();
k_point_int_list = Array{k_point_int_Tuple,1}();
for (q_i,q_point) in enumerate(q_point_list)
  push!(q_point_int_list,k_point_float2int(q_point))
end
for (k_i,k_point) in enumerate(k_point_int_list)
  push!(k_point_int_list,k_point_float2int(k_point))
end

q_point_int_list_matlab = reinterpret(Int64,q_point_int_list,(3,length(q_point_int_list)))';
k_point_int_list_matlab = reinterpret(Int64,k_point_int_list,(3,length(k_point_int_list)))';

###############################################################################
# Write to MAT
###############################################################################
tv = get_dftdataset().scf_r.tv;
rv = get_dftdataset().scf_r.rv;
Gxy = get_dftdataset().scf_r.Gxyz;
atom_num = get_dftdataset().scf_r.atomnum;

#jq_output_file = "test.mat"
f_name = string(cal_name,"_meshk_",atom1,"_",atom2,"_[all]","_ChemPdelta_",ChemP_delta_ev);
if (orbital_mask_on)
    println(" ", orbital_mask1_inv," : ",orbital_mask2_inv)
    f_name = string(cal_name,"_meshk_",atom1,"_",atom2);
    mask_name = string("_atom1m_[",join(orbital_mask1_inv,","),
    "]_atom2m_[",join(orbital_mask2_inv,","),"]");
    f_name = string(f_name,mask_name,"_[",orbital_mask_name,"]","_ChemPdelta_",ChemP_delta_ev);
end
result_fname = string("jq_",f_name,".mat");


println(jq_output_dir)
jq_output_file = joinpath(jq_output_dir,result_fname)

MAT.matwrite(jq_output_file,Dict("Jij_Q_matlab" =>X_Q_matlab
  #,"Jij_Q_K" => Jij_Q_K_matlab
  ,"q_point_list" => q_point_int_list_matlab
  ,"k_point_list" => k_point_int_list_matlab
  ,"k_point_precision" => k_point_precision
  ,"tv" => tv
  ,"rv" => rv
  ,"Gxyz" => scf_r.Gxyz
  ,"atomnum" => scf_r.atomnum
  ,"atom1" => atom1
  ,"atom2" => atom2
  ,"scf_name" => scf_name
  ,"orbital_mask1" => orbital_mask1
  ,"orbital_mask2" => orbital_mask2
  #,"Jij_history" => cal_history_dat["Jij_history"]
  ,"orbital_mask_on" => orbital_mask_on
  #,"orbital_mask1_inv" => orbital_mask1_inv
  #,"orbital_mask2_inv" => orbital_mask2_inv
  ,"ChemP_delta" => 0.0
  ,"X_VERSION" => string(X_VERSION)
  ));

###############################################################################
# Clear HDF5 cache
###############################################################################
println("hdf_cache_name:",hdf_cache_name)
if (isfile(hdf_cache_name))
  rm(hdf_cache_name)
  file = open(hdf_cache_name,"w");
  close(file);
end
