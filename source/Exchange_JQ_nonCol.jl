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

scf_name = "../examples/NiO/4cell/Ni6.0_O5.0-s2p2d2f1-4.25_0.0-4.180-k10/nio.scfout"
#scf_name = "/home/users1/bluehope/work_local/GaVS/LDA_FS_rhom/GaV4S8.scfout"
#scf_name = "/home/users1/bluehope/work_local/GaVS/LDA_FS_cubic/GaV4S8.scfout"
#scf_name = "/home/users1/bluehope/work/Jx/X0/FeTe/global/FeTe_PM_root2xroot2_AF_s2p2d2f1_expstr_U0.0.LDA/feTe.scfout"
#scf_name = "/home/users1/bluehope/work/Jx/X0/SRO/1x1/SRO_ortho.scfout"

##############################################################################
orbital_mask1 = Array{Int64,1}();
orbital_mask2 = Array{Int64,1}();
orbital_mask_name = "";
orbital_mask_option = DFTcommon.nomask;
orbital_mask_on = false;
atom1 = -1;
atom2 = -1;
atom12_list = Vector{Tuple{Int64,Int64}}();

if true
  arg_input = parse_input(ARGS)
  scf_name = arg_input.scf_name
  ChemP_delta_ev = arg_input.ChemP_delta_ev
  k_point_num = arg_input.k_point_num
  q_point_num = arg_input.q_point_num
  #atom1 = arg_input.atom1;
  #atom2 = arg_input.atom2;
  atom12_list = arg_input.atom12_list;
  hdftmpdir = arg_input.hdftmpdir;
end
###
# TOML override (should be merged to DFTcommon.jl)
###
#=
import TOML
toml_inputs = TOML.parse(readstring("nio_J.toml"))
atom1 = toml_inputs["atom12"][1][1];
atom2 = toml_inputs["atom12"][1][2];
scf_name = toml_inputs["scf_fname"];
#k_point_num = [3,3,3]
#q_point_num = [3,3,3]
k_point_num = [5,5,5]
q_point_num = [5,5,5]
=#
ChemP_delta_ev = 0.0

###############################################################################

println(atom12_list)

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

#scf_r = set_current_dftdataset(scf_name, DFTforge.OpenMX, DFTforge.colinear_type)
#DFTforge.pwork(set_current_dftdataset,(scf_name, DFTforge.OpenMX, DFTforge.colinear_type,1));
scf_r = set_current_dftdataset(scf_name, DFTforge.OpenMX, DFTforge.non_colinear_type)
DFTforge.pwork(set_current_dftdataset,(scf_name, DFTforge.OpenMX, DFTforge.non_colinear_type,1));
##


#k_point_num = [2,2,2]
#q_point_num = [2,2,2]

k_point_list = kPoint_gen_GammaCenter(k_point_num);
q_point_list = kPoint_gen_GammaCenter(q_point_num);
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
orbital_mask_input = orbital_mask_input_Type(orbital_mask1,orbital_mask2,(atom1,atom2),false)
if ((DFTcommon.unmask == orbital_mask_option) || (DFTcommon.mask == orbital_mask_option) )
    orbital_mask_input = orbital_mask_input_Type(orbital_mask1,orbital_mask2,(atom1,atom2),true)
end
pwork(init_orbital_mask,orbital_mask_input);

#@everywhere function Init_SmallHks(atom12)
#    atom1 = atom12[1];
#    atom2 = atom12[2];
#    global SmallHks;
#    result_index = 1;
#    SmallHks =  DFTforge.OpenMXdata.test_SmallHks(atom1,atom2,get_dftdataset(result_index).scf_r);
#end


#DFTforge.pwork(Init_SmallHks,(atom1,atom2))

@everywhere function Magnetic_Exchange_J_noncolinear(input::Job_input_kq_atom_list_Type)
    #global SmallHks;
    ############################################################################
    ## Accessing Data Start
    ############################################################################
    # Common input info
    k_point::DFTforge.k_point_Tuple =  input.k_point
    kq_point::DFTforge.k_point_Tuple =  input.kq_point
    spin_type::DFTforge.SPINtype = input.spin_type;

    #atom1::Int = input.atom12[1][1];
    #atom2::Int = input.atom12[1][2];
    atom12_list::Vector{Tuple{Int64,Int64}} = input.atom12_list;
    result_mat = zeros(Complex_my,10,length(atom12_list))

    #atom1::Int = input.atom1;
    #atom2::Int = input.atom2;

    result_index =  input.result_index;
    cache_index = input.cache_index;

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
    eigenstate_kq = cacheread_eigenstate(kq_point,1,cache_index)
    eigenstate_k  = cacheread_eigenstate(k_point,1,cache_index)
    #eigenstate_kq_down::Kpoint_eigenstate = cacheread_eigenstate(kq_point,2,cache_index)

    # Get Hamiltonian
    Hks_updown = cacheread_Hamiltonian(1,cache_index)
    #Hks_down = cacheread_Hamiltonian(2,cache_index)
    (orbitalStartIdx,orbitalNums) = cacheread_atomsOrbital_lists(cache_index)

    #En_k_up::Array{Float_my,1} = eigenstate_k_up.Eigenvalues;
    Em_kq = eigenstate_kq.Eigenvalues;
    En_k  = eigenstate_k.Eigenvalues;
    #Em_kq_down::Array{Float_my,1} = eigenstate_kq_down.Eigenvalues;

    #Es_n_k_up::Array{Complex_my,2} = eigenstate_k_up.Eigenstate;
    Es_m_kq  = eigenstate_kq.Eigenstate;
    Es_n_k  = eigenstate_k.Eigenstate;
    #Plots.heatmap(real(Es_n_k))
    #Es_m_kq_down::Array{Complex_my,2} = eigenstate_kq_down.Eigenstate;
    Fftn_k  = 1.0./(exp( ((En_k)  - ChemP)/(kB*E_temp)) + 1.0 );
    Fftm_kq = 1.0./(exp( ((Em_kq) - ChemP)/(kB*E_temp)) + 1.0 );
    dFnk_Fmkq  =
      Fftn_k*ones(1,TotalOrbitalNum2)  - ones(TotalOrbitalNum2,1)*Fftm_kq' ;

    Enk_Emkq =
      En_k[:]*ones(1,TotalOrbitalNum2) - ones(TotalOrbitalNum2,1)*Em_kq[:]' ;
    Enk_Emkq += im*0.0001;

    part1 = dFnk_Fmkq./(-Enk_Emkq);

    for (atom12_i,atom12) in enumerate(atom12_list)
      atom1 = atom12[1]
      atom2 = atom12[2]

      atom1_orbits_up   = orbitalStartIdx[atom1] + (1:orbitalNums[atom1])
      atom1_orbits_down = TotalOrbitalNum + atom1_orbits_up
      atom2_orbits_up   = orbitalStartIdx[atom2] + (1:orbitalNums[atom2])
      atom2_orbits_down = TotalOrbitalNum + atom2_orbits_up
      #Plots.heatmap(real(Hks_updown))

      ############################################################################
      ## Accessing Data End
      ############################################################################
      ## Do auctual calucations
      Vz_1 = 0.5*(Hks_updown[atom1_orbits_up,atom1_orbits_up] -
        Hks_updown[atom1_orbits_down,atom1_orbits_down])
      Voff_1 = Hks_updown[atom1_orbits_up,atom1_orbits_down]
      Vx_1 = 1*(Voff_1 + conj(Voff_1))/2.0;
      Vy_1 = 1*(Voff_1 - conj(Voff_1))/(2.0*im);
      #Plots.heatmap(real(Vz_1))
      #Plots.heatmap(real(Vy_1))

      Vz_2 = 0.5*(Hks_updown[atom2_orbits_up,atom2_orbits_up] -
        Hks_updown[atom2_orbits_down,atom2_orbits_down])
      Voff_2 = Hks_updown[atom2_orbits_up,atom2_orbits_down]
      Vx_2 = 1*(Voff_2 + conj(Voff_2))/2.0;
      Vy_2 = 1*(Voff_2 - conj(Voff_2))/(2.0*im);

      #Plots.heatmap(real(part1))

      G1V1_z = (Es_n_k[atom1_orbits_down,:]'  *Vz_1* Es_m_kq[atom1_orbits_up,:]);
      G2V2_z = (Es_m_kq[atom2_orbits_up,:]'   *Vz_2* Es_n_k[atom2_orbits_down,:]);
      #G2V2_z = (Es_m_kq[atom2_orbits_up,:]' * Vz_2 * Es_n_k[atom2_orbits_down,:]);
      G1V1_x = (Es_n_k[atom1_orbits_down,:]'  *Vx_1* Es_m_kq[atom1_orbits_up,:]);
      G2V2_x = (Es_m_kq[atom2_orbits_up,:]'   *Vx_2* Es_n_k[atom2_orbits_down,:]);

      G1V1_y = (Es_n_k[atom1_orbits_down,:]'  *Vy_1* Es_m_kq[atom1_orbits_up,:]);
      G2V2_y = (Es_m_kq[atom2_orbits_up,:]'   *Vy_2* Es_n_k[atom2_orbits_down,:]);
      #J_ij_zz::Array{Complex_my,2} =  0.5* part1.* transpose(G1V1_z) .* G2V2_z * Hartree2cm;

      J_ij_xx =  0.5* sum(part1.* transpose(G1V1_x) .* G2V2_x )* Hartree2cm;
      J_ij_xy =  0.5* sum(part1.* transpose(G1V1_x) .* G2V2_y )* Hartree2cm;
      J_ij_xz =  0.5* sum(part1.* transpose(G1V1_x) .* G2V2_z )* Hartree2cm;

      J_ij_yx =  0.5* sum(part1.* transpose(G1V1_y) .* G2V2_x )* Hartree2cm;
      J_ij_yy =  0.5* sum(part1.* transpose(G1V1_y) .* G2V2_y )* Hartree2cm;
      J_ij_yz =  0.5* sum(part1.* transpose(G1V1_y) .* G2V2_z )* Hartree2cm;

      J_ij_zx =  0.5* sum(part1.* transpose(G1V1_z) .* G2V2_x )* Hartree2cm;
      J_ij_zy =  0.5* sum(part1.* transpose(G1V1_z) .* G2V2_y )* Hartree2cm;
      J_ij_zz =  0.5* sum(part1.* transpose(G1V1_z) .* G2V2_z )* Hartree2cm;

      X_ij = sum(part1.*
      transpose(Es_n_k[atom1_orbits_down,:]' * Es_m_kq[atom1_orbits_up,:]) .*
      (Es_m_kq[atom2_orbits_up,:]' * Es_n_k[atom2_orbits_down,:]));

      result_mat[:,atom12_i] = [
      J_ij_xx,J_ij_xy,J_ij_xz,
      J_ij_yx,J_ij_yy,J_ij_yz,
      J_ij_zx,J_ij_zy,J_ij_zz, X_ij];
    end
    return result_mat #sum(J_ij[:]);
  end



##############################
#DEBUG
##############################
#cacheset(eigenstate_cache)

##############################
#DEBUG END
##############################


###############################################################################
# Average X_Q results
###############################################################################
#X_Q::Dict{k_point_int_Tuple,Array{Complex_my,1}} = X_succeptability();
#atom12_list =[(atom1,atom2)];
(X_Q_nc,X_Q_mean_nc) = Qspace_Ksum_atomlist_parallel_nc(Magnetic_Exchange_J_noncolinear,
q_point_list,k_point_list,atom12_list)
println(typeof(X_Q_mean_nc))
#X_Q_matlab = zeros(Complex128,length(q_point_list),10,length(atom12_list));
#for (q_i,q_point) in enumerate(q_point_list)
  #q_point_int = k_point_float2int(q_point);
  #X_Q_matlab[q_i] = mean(X_Q[q_point_int])
#end
Xij_Q_mean_matlab = Array(Array{Complex_my,1},10,length(atom12_list));
for (atom12_i,atom12) in enumerate(atom12_list)
  atom1 = atom12[1];
  atom2 = atom12[2];
  for xyz_i = 1:10
    Xij_Q_mean_matlab[xyz_i,atom12_i] = zeros(Complex_my,length(q_point_list));
    for (q_i,q_point) in enumerate(q_point_list)
      q_point_int = k_point_float2int(q_point);
      Xij_Q_mean_matlab[xyz_i,atom12_i][q_i] =
        X_Q_mean_nc[xyz_i,atom12_i][q_point_int];
    end
  end
end


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
for (atom12_i,atom12) in enumerate(atom12_list)
  atom1 = atom12[1];
  atom2 = atom12[2];

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

  MAT.matwrite(jq_output_file,Dict("Jij_Q_matlab" =>Xij_Q_mean_matlab[:,atom12_i]
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
end
###############################################################################
# Clear HDF5 cache
###############################################################################
println("hdf_cache_name:",hdf_cache_name)
if (isfile(hdf_cache_name))
  rm(hdf_cache_name)
  file = open(hdf_cache_name,"w");
  close(file);
end
