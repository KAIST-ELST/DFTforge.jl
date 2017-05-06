using ProgressMeter
import DFTforge
using DFTforge.DFTrefinery
using DFTcommon
import MAT
X_VERSION = VersionNumber("0.2.0-dev+20170503");
println("X_VERSION: ",X_VERSION)

@everywhere import DFTforge
@everywhere using DFTforge.DFTrefinery
@everywhere using DFTcommon
##############################################################################
## 1. Read INPUT
## 1.1 Set Default values
## 1.2 Read input from argument & TOML file
## 1.3 Set values from intput (arg_input)
## 1.4 Set caluations type and ouput folder
##############################################################################
hdftmpdir = ""
## 1.1 Set Default values
orbital_mask1 = Array{Int64,1}();
orbital_mask2 = Array{Int64,1}();
orbital_mask_name = "";
orbital_mask_option = DFTcommon.nomask;
orbital_mask_on = false;

k_point_num = [3,3,3]
q_point_num = [3,3,3]
ChemP_delta_ev = 0.0
DFT_type = DFTcommon.OpenMX

## 1.2 Read input from argument & TOML file
arg_input = DFTcommon.Arg_Inputs();
arg_input = parse_input(ARGS,arg_input)
#arg_input.TOMLinput = "nio_J_wannier.toml" # Debug
#arg_input.TOMLinput = "nio_J_openmx.toml" # Debug
arg_input = parse_TOML(arg_input.TOMLinput,arg_input)
# let argument override
arg_input = parse_input(ARGS,arg_input)

## 1.3 Set values from intput (arg_input)
DFT_type = arg_input.DFT_type
Wannier90_type = arg_input.Wannier90_type

result_file = arg_input.result_file
ChemP_delta_ev = arg_input.ChemP_delta_ev
 # k,q point num
k_point_num = arg_input.k_point_num
q_point_num = arg_input.q_point_num
 # atom 12
atom1 = arg_input.atom1;
atom2 = arg_input.atom2;
atom12_list = arg_input.atom12_list;
hdftmpdir = arg_input.hdftmpdir;

 # orbital mask
orbital_mask_option = arg_input.orbital_mask_option;
if ((DFTcommon.unmask == orbital_mask_option) || (DFTcommon.mask == orbital_mask_option) )
  orbital_mask_name = arg_input.orbital_mask_name
  orbital_mask1 = arg_input.orbital_mask1;
  orbital_mask2 = arg_input.orbital_mask2;
  orbital_mask_on = true
end

println("===================================================")
println(atom12_list)
println(string("DFT_type ",DFT_type))
println(string("orbital_mask_option ",orbital_mask_option,"\t",orbital_mask_name))

## 1.4 Set caluations type and ouput folder
cal_type = "jq" # xq, ...

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
println("===================================================")
##############################################################################
## 2. Calculate & Store k,q points information
## 2.1 Set Input info
## 2.2 Generate k,q points
## 2.3 Calculate Eigenstate & Store Eigenstate into file in HDF5 format
## 2.4 Send Eigenstate info to child processes
##############################################################################

## 2.1 Set Input info
scf_r = [];
if (DFTcommon.OpenMX == DFT_type)
  scf_r = set_current_dftdataset(result_file, DFTcommon.OpenMX, DFTcommon.colinear_type)
elseif (DFTcommon.Wannier90 == DFT_type)
  atomnum = arg_input.Wannier_Optional_Info.atomnum
  atompos = arg_input.Wannier_Optional_Info.atompos
  atoms_orbitals_list = arg_input.Wannier_Optional_Info.atoms_orbitals_list

  scf_r = DFTforge.read_dftresult(result_file,DFT_type,Wannier90_type,atoms_orbitals_list,atomnum,atompos)
  scf_r = set_current_dftdataset(scf_r, DFT_type, DFTcommon.colinear_type)
end

DFTforge.pwork(set_current_dftdataset,(scf_r, DFT_type, DFTcommon.colinear_type,1));

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
tic();
DFTforge.pwork(cacheread_lampup,kq_point_list)
toc();
################################################################################

##############################################################################
## 3. Setup extra infos (orbital, chemp shift)
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
if (orbital_mask_on)
    orbital_mask_input = orbital_mask_input_Type(orbital_mask1,orbital_mask2,(-1,-1),true)
end
DFTforge.pwork(init_orbital_mask,orbital_mask_input)
DFTforge.pwork(init_variables,ChemP_delta_ev)


##############################################################################
## Physical properties calculation define section
## 4. Magnetic exchange function define
## 4.1 Do K,Q sum
## 4.2 reduce K,Q to Q space
##############################################################################

## 4. Magnetic exchange function define
@everywhere function Magnetic_Exchange_J_colinear(input::Job_input_kq_atom_list_Type)
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

  atom12_list::Vector{Tuple{Int64,Int64}} = input.atom12_list;
  result_mat = zeros(Complex_my,1,length(atom12_list))
  #atom1::Int = input.atom1;
  #atom2::Int = input.atom2;

  result_index::Int = input.result_index;
  cache_index::Int = input.cache_index;

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
  eigenstate_kq_up::Kpoint_eigenstate = cacheread_eigenstate(kq_point,1,cache_index)
  eigenstate_k_down::Kpoint_eigenstate  = cacheread_eigenstate(k_point,2,cache_index)
  #eigenstate_kq_down::Kpoint_eigenstate = cacheread_eigenstate(kq_point,2,cache_index)

  # Get Hamiltonian
  Hks_up::Hamiltonian_type = cacheread_Hamiltonian(1,cache_index)
  Hks_down::Hamiltonian_type = cacheread_Hamiltonian(2,cache_index)
  (orbitalStartIdx,orbitalNums) = cacheread_atomsOrbital_lists(cache_index)

  #En_k_up::Array{Float_my,1} = eigenstate_k_up.Eigenvalues;
  Em_kq_up::Array{Float_my,1} = eigenstate_kq_up.Eigenvalues;
  En_k_down::Array{Float_my,1} = eigenstate_k_down.Eigenvalues;
  #Em_kq_down::Array{Float_my,1} = eigenstate_kq_down.Eigenvalues;

  #Es_n_k_up::Array{Complex_my,2} = eigenstate_k_up.Eigenstate;
  Es_m_kq_up::Array{Complex_my,2} = eigenstate_kq_up.Eigenstate;
  Es_n_k_down::Array{Complex_my,2} = eigenstate_k_down.Eigenstate;
  #Es_m_kq_down::Array{Complex_my,2} = eigenstate_kq_down.Eigenstate;

  for (atom12_i,atom12) in enumerate(atom12_list)
    atom1 = atom12[1]
    atom2 = atom12[2]

    atom1_orbitals = orbitalStartIdx[atom1]+(1:orbitalNums[atom1]);
    atom2_orbitals = orbitalStartIdx[atom2]+(1:orbitalNums[atom2]);
    ############################################################################
    ## Accessing Data End
    ############################################################################

    # mask oribtals
    Es_n_k_down_atom1 = copy(Es_n_k_down);
    Es_m_kq_up_atom1  = copy(Es_m_kq_up);
    Es_m_kq_up_atom2  = copy(Es_m_kq_up);
    Es_n_k_down_atom2 = copy(Es_n_k_down);
    if (length(orbital_mask1)>0)
      orbital_mask1_tmp = collect(1:orbitalNums[atom1]);
      for orbit1 in orbital_mask1
          deleteat!(orbital_mask1_tmp, find(orbital_mask1_tmp.==orbit1))
      end
      Es_n_k_down_atom1[orbitalStartIdx[atom1]+orbital_mask1_tmp,:]=0.0;
    end
    if (length(orbital_mask2)>0)
      orbital_mask2_tmp = collect(1:orbitalNums[atom2]);
      for orbit2 in orbital_mask2
          deleteat!(orbital_mask2_tmp, find(orbital_mask2_tmp.==orbit2))
      end
      Es_m_kq_up_atom2[orbitalStartIdx[atom2]+orbital_mask2_tmp,:]=0.0;
    end

    ## Do auctual calucations
    Fftn_k::Array{Float_my,1}  = 1.0./(exp( ((En_k_down)  - ChemP)/(kBeV*E_temp)) + 1.0 );
    Fftm_kq::Array{Float_my,1} = 1.0./(exp( ((Em_kq_up) - ChemP)/(kBeV*E_temp)) + 1.0 );

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

    #V1_2 = 0.5 *( SmallHks[1][atom1_orbitals_rel,atom1_orbitals_rel] -
    #  SmallHks[2][atom1_orbitals_rel,atom1_orbitals_rel] );
    #V2_2 = 0.5 *( SmallHks[1][atom2_orbitals_rel2,atom2_orbitals_rel2] -
    #   SmallHks[2][atom2_orbitals_rel2,atom2_orbitals_rel2] );

    VV1::Array{Complex_my,2} = (Es_n_k_down_atom1[atom1_orbitals,:]' * V1 * Es_m_kq_up_atom1[atom1_orbitals,:]);
    VV2::Array{Complex_my,2} = (Es_m_kq_up_atom2[atom2_orbitals,:]' * V2 * Es_n_k_down_atom2[atom2_orbitals,:]);

    #VV1::Array{Complex_my,2} = (Es_n_k_down_atom1[atom1_orbitals,:]' * Es_m_kq_up_atom1[atom1_orbitals,:]);
    #VV2::Array{Complex_my,2} = (Es_m_kq_up_atom2[atom2_orbitals,:]'  * Es_n_k_down_atom2[atom2_orbitals,:]);

    Vi_Vj = transpose(VV1).*VV2;

    J_ij::Array{Complex_my,2} =  0.5./(-Enk_Emkq).*dFnk_Fmkq .* Vi_Vj ;
    #return sum(J_ij[:])*Hartree2cm;
    #return sum(J_ij[!isnan(J_ij)] )*Hartree2cm;
    result_mat[1,atom12_i] = sum(J_ij[!isnan(J_ij)] );
  end
  return result_mat
end

num_return = 1;

## 4.1 Do K,Q sum
(X_Q_nc,X_Q_mean_nc) = Qspace_Ksum_atomlist_parallel(Magnetic_Exchange_J_colinear,
q_point_list,k_point_list,atom12_list,num_return)
#println(typeof(X_Q_mean_nc))
println("===================================================")
## 4.2 reduce K,Q to Q space
# Average X_Q results
Xij_Q_mean_matlab = Array(Array{Complex_my,1},num_return,length(atom12_list));
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
    println(string(" Gamma point J [",atom1,",",atom2,"]: ",
      1000.0*mean(Xij_Q_mean_matlab[xyz_i,atom12_i][:])," meV"))
  end
end


###############################################################################
## 5. Save results and clear hdf5 file
## 5.1 Prepaire infomations for outout
## 5.2 Write to MAT
## 5.3 Cleanup HDF5 cache file
###############################################################################

## 5.1 Prepaire infomations for outout
q_point_int_list = Array{k_point_int_Tuple,1}();
k_point_int_list = Array{k_point_int_Tuple,1}();
for (q_i,q_point) in enumerate(q_point_list)
  push!(q_point_int_list,k_point_float2int(q_point))
end
for (k_i,k_point) in enumerate(k_point_list)
  push!(k_point_int_list,k_point_float2int(k_point))
end

q_point_int_list_matlab = reinterpret(Int64,q_point_int_list,(3,length(q_point_int_list)))';
k_point_int_list_matlab = reinterpret(Int64,k_point_int_list,(3,length(k_point_int_list)))';





tv = get_dftdataset().scf_r.tv;
rv = get_dftdataset().scf_r.rv;
Gxy = get_dftdataset().scf_r.Gxyz;
atom_num = get_dftdataset().scf_r.atomnum;
println(jq_output_dir)
#jq_output_file = "test.mat"
for (atom12_i,atom12) in enumerate(atom12_list)
  atom1 = atom12[1];
  atom2 = atom12[2];

  f_name = string(cal_name,"_meshk_",atom1,"_",atom2,"_[all]","_ChemPdelta_",ChemP_delta_ev);
  if (orbital_mask_on)
      #println(" ", orbital_mask1_inv," : ",orbital_mask2_inv)
      f_name = string(cal_name,"_meshk_",atom1,"_",atom2);
      #mask_name = string("_atom1m_[",join(orbital_mask1_inv,","),
      #"]_atom2m_[",join(orbital_mask2_inv,","),"]");
      mask_name = string("_atom1m_[",",", "]_atom2m_[",",","]");
      f_name = string(f_name,mask_name,"_[",orbital_mask_name,"]","_ChemPdelta_",ChemP_delta_ev);
  end
  result_fname = string("jq_",f_name,".mat");

  jq_output_file = joinpath(jq_output_dir,result_fname)

## 5.2 Write to MAT
  println(jq_output_file)
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
    ,"cal_name" => cal_name
    ,"orbital_mask1" => orbital_mask1
    ,"orbital_mask2" => orbital_mask2
    #,"Jij_history" => cal_history_dat["Jij_history"]
    ,"orbital_mask_on" => orbital_mask_on
    #,"orbital_mask1_inv" => orbital_mask1_inv
    #,"orbital_mask2_inv" => orbital_mask2_inv
    ,"ChemP_delta" => ChemP_delta_ev
    ,"X_VERSION" => string(X_VERSION)
    ));
end

## 5.3 Cleanup HDF5 cache file
println("hdf_cache_name:",hdf_cache_name)
if (isfile(hdf_cache_name))
  rm(hdf_cache_name)
  file = open(hdf_cache_name,"w");
  close(file);
end
