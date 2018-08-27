
##import MAT
import FileIO
function export2mat_K_Q(Xij_Q_mean_matlab,hamiltonian_info,
  q_point_list::Array{k_point_Tuple},
  k_point_list::Array{k_point_Tuple},atom12_list::Vector{Tuple{Int64,Int64}},
  orbital_mask_on::Bool,orbital_mask1,orbital_mask2,ChemP_delta_ev::Float64,
  optionalOutputDict::Dict{AbstractString,Any},
  jq_output_dir::AbstractString,cal_name::AbstractString,
  orbital_mask_name::AbstractString,
  cal_type::AbstractString)
  export2mat_K_Q(Xij_Q_mean_matlab,hamiltonian_info,
      q_point_list,
      k_point_list,atom12_list,
      orbital_mask_on,
      orbital_mask1,orbital_mask2,Array{Int}(undef,0),Array{Int}(undef,0),
      ChemP_delta_ev,
      optionalOutputDict,
      jq_output_dir,cal_name,
      orbital_mask_name,
      cal_type)

end
function export2mat_K_Q(Xij_Q_mean_matlab,hamiltonian_info,
  q_point_list::Array{k_point_Tuple},
  k_point_list::Array{k_point_Tuple},atom12_list::Vector{Tuple{Int64,Int64}},
  orbital_mask_on::Bool,
  orbital_mask1,orbital_mask2,orbital_mask3,orbital_mask4,
  ChemP_delta_ev::Float64,
  optionalOutputDict::Dict{AbstractString,Any},
  jq_output_dir::AbstractString,cal_name::AbstractString,
  orbital_mask_name::AbstractString,
  cal_type::AbstractString)


  export2mat_K_Q(Xij_Q_mean_matlab,hamiltonian_info,
    q_point_list,
    k_point_list,atom12_list,
    orbital_mask_on,
    orbital_mask1,orbital_mask2,orbital_mask3,orbital_mask4,
    [],
    ChemP_delta_ev,
    optionalOutputDict,
    jq_output_dir,cal_name,
    orbital_mask_name,
    cal_type)
end
function export2mat_K_Q(Xij_Q_mean_matlab,hamiltonian_info,
  q_point_list::Array{k_point_Tuple},
  k_point_list::Array{k_point_Tuple},atom12_list::Vector{Tuple{Int64,Int64}},
  orbital_mask_on::Bool,
  orbital_mask1,orbital_mask2,orbital_mask3,orbital_mask4,
  energywindow_all1234_list,
  ChemP_delta_ev::Float64,
  optionalOutputDict::Dict{AbstractString,Any},
  jq_output_dir::AbstractString,cal_name::AbstractString,
  orbital_mask_name::AbstractString,
  cal_type::AbstractString)

  scf_r = hamiltonian_info.scf_r;
## .1 Prepaire infomations for outout
  q_point_int_list = Array{k_point_int_Tuple,1}();
  k_point_int_list = Array{k_point_int_Tuple,1}();
  for (q_i,q_point) in enumerate(q_point_list)
    push!(q_point_int_list,k_point_float2int(q_point))
  end
  for (k_i,k_point) in enumerate(k_point_list)
    push!(k_point_int_list,k_point_float2int(k_point))
  end

  q_point_int_list_matlab = collect(reshape(reinterpret(Int64,q_point_int_list),(3,length(q_point_int_list)))');
  k_point_int_list_matlab = collect(reshape(reinterpret(Int64,k_point_int_list),(3,length(k_point_int_list)))');

  tv = scf_r.tv;
  rv = scf_r.rv;
  Gxy = scf_r.Gxyz;
  atom_num = scf_r.atomnum;
  println(jq_output_dir)
  #jq_output_file = "test.mat"
  for (atom12_i,atom12) in enumerate(atom12_list)
    atom1 = atom12[1];
    atom2 = atom12[2];

    f_name = string(cal_name,"_atomij_",atom1,"_",atom2,"_[all_all]","_ChemPdelta_",ChemP_delta_ev);
    if (orbital_mask_on)
        #println(" ", orbital_mask1_inv," : ",orbital_mask2_inv)
        f_name = string(cal_name,"_atomij_",atom1,"_",atom2);
        #mask_name = string("_atom1m_[",join(orbital_mask1_inv,","),
        #"]_atom2m_[",join(orbital_mask2_inv,","),"]");
        mask_name = string("_atom1m_[",",", "]_atom2m_[",",","]");
        f_name = string(f_name,mask_name,"_[",orbital_mask_name,"]","_ChemPdelta_",ChemP_delta_ev);
    end
    #result_fname = string(cal_type,"_",f_name,".mat");
    result_fname = string(cal_type,"_",f_name,".jld2");


    jq_output_file = joinpath(jq_output_dir,result_fname)

  ## .2 Write to MAT
    println(basename(jq_output_file))
    outputDict = Dict("Jij_Q_matlab" =>Xij_Q_mean_matlab[:,atom12_i]
      #,"Jij_Q_K" => Jij_Q_K_matlab
      ,"q_point_list" => q_point_int_list_matlab
      ,"k_point_list" => k_point_int_list_matlab
      ,"k_point_precision" => DFTcommon.k_point_precision
      ,"tv" => tv
      ,"rv" => rv
      ,"Gxyz" => scf_r.Gxyz
      ,"atomnum" => scf_r.atomnum
      ,"atom1" => atom1
      ,"atom2" => atom2
      ,"cal_name" => cal_name
      ,"orbital_mask1" => orbital_mask1
      ,"orbital_mask2" => orbital_mask2
      ,"orbital_mask3" => orbital_mask3
      ,"orbital_mask4" => orbital_mask4
      ,"energywindow_all1234_list" => energywindow_all1234_list
      #,"Jij_history" => cal_history_dat["Jij_history"]
      ,"orbital_mask_on" => orbital_mask_on
      #,"orbital_mask1_inv" => orbital_mask1_inv
      #,"orbital_mask2_inv" => orbital_mask2_inv
      ,"ChemP_delta" => ChemP_delta_ev
      #,"X_VERSION" => string(X_VERSION)
      ,"DFTforge_VERSION" => string(DFTforge.get_DFTforge_VERSION())
      )
  ## Add optional infos
    for (key,value) in optionalOutputDict
      outputDict[key] = value;
    end
    #println(keys(outputDict))
    #println(map(x-> typeof(x), values(outputDict)))
  ## Write
    #MAT.matwrite(jq_output_file,outputDict);
    FileIO.save(jq_output_file,outputDict);
  end
end


function export2mat_K_Q_nc(Xij_Q_mean_matlab,theta_phi_list,hamiltonian_info,
  q_point_list::Array{k_point_Tuple},
  k_point_list::Array{k_point_Tuple},atom12_list::Vector{Tuple{Int64,Int64}},
  orbital_mask_on::Bool,
  orbital_mask1,orbital_mask2,orbital_mask3,orbital_mask4,
  ChemP_delta_ev::Float64,
  optionalOutputDict::Dict{AbstractString,Any},
  jq_output_dir::AbstractString,cal_name::AbstractString,
  orbital_mask_name::AbstractString,
  cal_type::AbstractString)

  scf_r = hamiltonian_info.scf_r;
## .1 Prepaire infomations for outout
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

  tv = scf_r.tv;
  rv = scf_r.rv;
  Gxy = scf_r.Gxyz;
  atom_num = scf_r.atomnum;
  println(jq_output_dir)
  #jq_output_file = "test.mat"
  for (atom12_i,atom12) in enumerate(atom12_list)
    atom1 = atom12[1];
    atom2 = atom12[2];

    f_name = string(cal_name,"_meshk_",atom1,"_",atom2,"_[all_all]","_ChemPdelta_",ChemP_delta_ev);
    if (orbital_mask_on)
        #println(" ", orbital_mask1_inv," : ",orbital_mask2_inv)
        f_name = string(cal_name,"_meshk_",atom1,"_",atom2);
        #mask_name = string("_atom1m_[",join(orbital_mask1_inv,","),
        #"]_atom2m_[",join(orbital_mask2_inv,","),"]");
        mask_name = string("_atom1m_[",",", "]_atom2m_[",",","]");
        f_name = string(f_name,mask_name,"_[",orbital_mask_name,"]","_ChemPdelta_",ChemP_delta_ev);
    end
    #result_fname = string(cal_type,"_",f_name,".mat");
    result_fname = string(cal_type,"_",f_name,".jld2");

    jq_output_file = joinpath(jq_output_dir,result_fname)

  ## .2 Write to MAT
    println(basename(jq_output_file))
    outputDict = Dict("Jij_Q_matlab" =>Xij_Q_mean_matlab[:,:,atom12_i]
      #,"Jij_Q_K" => Jij_Q_K_matlab
      ,"q_point_list" => q_point_int_list_matlab
      ,"k_point_list" => k_point_int_list_matlab
      ,"k_point_precision" => DFTcommon.k_point_precision
      ,"tv" => tv
      ,"rv" => rv
      ,"Gxyz" => scf_r.Gxyz
      ,"atomnum" => scf_r.atomnum
      ,"atom1" => atom1
      ,"atom2" => atom2
      ,"cal_name" => cal_name
      ,"orbital_mask1" => orbital_mask1
      ,"orbital_mask2" => orbital_mask2
      ,"orbital_mask3" => orbital_mask3
      ,"orbital_mask4" => orbital_mask4
      #,"Jij_history" => cal_history_dat["Jij_history"]
      ,"orbital_mask_on" => orbital_mask_on
      #,"orbital_mask1_inv" => orbital_mask1_inv
      #,"orbital_mask2_inv" => orbital_mask2_inv
      ,"ChemP_delta" => ChemP_delta_ev
      #,"X_VERSION" => string(X_VERSION)
      ,"DFTforge_VERSION" => string(DFTforge.get_DFTforge_VERSION())
      )
  ## Add optional infos
    for (key,value) in optionalOutputDict
      outputDict[key] = value;
    end
    #println(keys(outputDict))
  ## Write
    #MAT.matwrite(jq_output_file,outputDict);
    FileIO.matwrite(jq_output_file,outputDict);
    end
end
