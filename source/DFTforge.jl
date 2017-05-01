__precompile__(false)

module DFTforge
using DFTcommon
export DFTtype, SPINtype
export k_point_Tuple,k_point_int_Tuple
export cal_colinear_eigenstate,cal_colinear_Hamiltonian
export cal_nonco_linear_Eigenstate,cal_noncolinear_Hamiltonian


#@enum DFTtype OpenMX = 1 Wannier90 = 2
@enum SPINtype para_type = 1 colinear_type = 2 non_colinear_type = 4

module OpenMXdata
include("backend/OpenMX_PostCommon.jl")
end
module Wannierdata
include("backend/Wannier_PostCommon.jl")
end



#type OpenMX_data
#    scf_name::AbstractString
#    scf_r::OpenMXdata.Openmxscf
#    dfttype::DFTforge.DFTtype
#end

function read_dftresult(scf_r::AbstractString, dfttype::DFTcommon.DFTtype)
    if (DFTcommon.OpenMX == dfttype)
      scf_r = OpenMXdata.read_scf(scf_r);
      return scf_r;
    end

end


function read_dftresult(wannier_fname::AbstractString,dfttype::DFTcommon.DFTtype,
  typeof_wannier::AbstractString,atoms_orbitals_list::Vector{Array{Int64}},
  atomnum::Int,atompos::Array{Float64,2})
  if (DFTcommon.Wannier90 == dfttype)
    wannier_r = Wannierdata.read_wannier(wannier_fname,typeof_wannier,
      atoms_orbitals_list,atomnum,atompos)
    return wannier_r;
  end
end




function cal_colinear_eigenstate(k_point::k_point_Tuple,
    dfttype::DFTcommon.DFTtype,scf_r,spin_list=1)
    Eigenstate::Array{Kpoint_eigenstate,1} = Array{Kpoint_eigenstate,1}();
    if (DFTcommon.OpenMX == dfttype)
      #Eigenstate::DFTforge.Kpoint_eigenstate =
      Eigenstate =
       OpenMXdata.cal_colinear_eigenstate(k_point,scf_r,spin_list);
    elseif (DFTcommon.Wannier90 == dfttype)
      Eigenstate = Wannierdata.cal_eigenstate(k_point,scf_r,spin_list)
    end

    return Eigenstate;
end

function cal_nonco_linear_Eigenstate(k_point::k_point_Tuple,
    dfttype::DFTcommon.DFTtype,scf_r)
    if (OpenMX == dfttype)
      EigenState = OpenMXdata.cal_noncolinear_eigenstate(k_point, scf_r)
      return EigenState
    elseif (Wannier90 == dfttype)
      Eigenstate = Wannierdata.cal_eigenstate(k_point,scf_r,[1])
      return EigenState
    end
end

function cal_colinear_Hamiltonian(dfttype::DFTcommon.DFTtype,scf_r,spin=1)

  H::Hamiltonian_type = zeros(Complex_my,2,2);
  if (DFTcommon.OpenMX==dfttype)
    H =  OpenMXdata.colinear_Hamiltonian(spin,scf_r)
  elseif (DFTcommon.Wannier90 == dfttype)
    H = Wannierdata.cal_Hamiltonian(scf_r,spin)
  end
  return H
end

function cal_noncolinear_Hamiltonian(dfttype::DFTcommon.DFTtype,
  Hmode::DFTcommon.nc_Hamiltonian_selection,scf_r)
  if (DFTcommon.OpenMX==dfttype)
    H = OpenMXdata.noncolinear_Hamiltonian(scf_r,Hmode);
    return H
  elseif (DFTcommon.Wannier90 == dfttype)
    H = Wannierdata.cal_Hamiltonian(k_point,scf_r,spin)
    return H
  end
end




################################################################################
module DFTrefinery
using DFTforge
using DFTcommon
using HDF5
using ProgressMeter

export set_current_dftdataset,cal_colinear_eigenstate,cal_Hamiltonian,get_dftdataset
export Job_input_Type,Job_input_kq_Type,Job_input_kq_atom_Type,Job_input_kq_atom_list_Type

export cachecal_all_Qpoint_eigenstats,cacheset,cacheread_eigenstate,cacheread,
  cacheread_lampup,cacheread_atomsOrbital_lists,cacheread_Hamiltonian
export get_ChempP
export Qspace_Ksum_parallel,Qspace_Ksum_atom_parallel,
  Kspace_parallel,Qspace_Ksum_atomlist_parallel_nc

  type DFTdataset
    dfttype::DFTcommon.DFTtype
    scf_r
    spin_type::SPINtype
    orbitalStartIdx::Array{Int,1}
    orbitalNums::Array{Int,1}
  end

  type Eigenstate_hdf5
    hdf_cache_name::AbstractString
    fid_hdf::HDF5.HDF5File
    q_points::Array{k_point_Tuple}
    q_points_int::Array{k_point_int_Tuple};
    q_points_intdic::Dict{k_point_int_Tuple,Int};
    spin_type::SPINtype
    TotalOrbitalNum::Int
    dftresult::DFTdataset
    Eigenvect_real::Array{Float64,4}
    Eigenvect_imag::Array{Float64,4}
    Eigenvalues::Array{Float64,3}
    Hamiltonian_real::Array{Float64,3}
    Hamiltonian_imag::Array{Float64,3}


    Eigenstate_hdf5(hdf_cache_name,fid_hdf,q_points,q_points_int,
      q_points_intdic,spin_type,TotalOrbitalNum,dftresult) =
    new(hdf_cache_name,fid_hdf,q_points,q_points_int,
      q_points_intdic,spin_type,TotalOrbitalNum,dftresult);
  end

  type Job_input_Type
    k_point::k_point_Tuple
    spin_type::SPINtype
    result_index::Int
    Job_input_Type(k_point,spin_type) = new(k_point,spin_type,1)
    Job_input_Type(k_point,spin_type,result_index) = new(k_point,spin_type,result_index)
  end

  type Job_input_kq_Type
    k_point::k_point_Tuple
    kq_point::k_point_Tuple
    spin_type::SPINtype
    result_index::Int
    cache_index::Int
    Job_input_kq_Type(k_point,kq_point,spin_type) =
      new(k_point,kq_point,spin_type,1,1)
    Job_input_kq_Type(k_point,kq_point,spin_type,result_index,cache_index) =
      new(k_point,kq_point,spin_type,result_index,cache_index)
  end
  type Job_input_kq_atom_Type
    k_point::k_point_Tuple
    kq_point::k_point_Tuple
    spin_type::SPINtype
    atom1::Int
    atom2::Int
    result_index::Int
    cache_index::Int
    Job_input_kq_atom_Type(k_point,kq_point,spin_type,atom1,atom2) =
      new(k_point,kq_point,spin_type,atom1,atom2,1,1)
    Job_input_kq_atom_Type(k_point,kq_point,spin_type,atom1,atom2,result_index,cache_index) =
      new(k_point,kq_point,spin_type,atom1,atom2,result_index,cache_index)
  end

  type Job_input_kq_atom_list_Type
    k_point::k_point_Tuple
    kq_point::k_point_Tuple
    spin_type::SPINtype
    atom12_list::Vector{Tuple{Int64,Int64}}

    result_index::Int
    cache_index::Int
    Job_input_kq_atom_list_Type(k_point,kq_point,spin_type,atom12) =
      new(k_point,kq_point,spin_type,atom12,1,1)
    Job_input_kq_atom_list_Type(k_point,kq_point,spin_type,atom12,result_index,cache_index) =
      new(k_point,kq_point,spin_type,atom12,result_index,cache_index)
  end

  global dftresult = Array{DFTdataset}();
  global eigenstate_list =  Array{Eigenstate_hdf5}(); #cached Eigenstates

  function set_current_dftdataset(scf_name::AbstractString,
    dfttype::DFTcommon.DFTtype,spin_type::SPINtype,result_index=1)
    if (DFTcommon.OpenMX == dfttype)
      # Read SCF and Set as current dftdata
      scf_r = DFTforge.OpenMXdata.read_scf(scf_name);
      orbitalStartIdx = zeros(Int,scf_r.atomnum);
      orbitalIdx::Int = 0; #각 atom별로 orbital index시작하는 지점
      orbitalNums = zeros(Int,scf_r.atomnum)
      for i = 1:scf_r.atomnum
          orbitalStartIdx[i] = orbitalIdx;
          orbitalIdx += scf_r.Total_NumOrbs[i]
          orbitalNums[i] = scf_r.Total_NumOrbs[i];
      end

      #assert(0 == scf_r.SpinP_switch - spin_type);
      dftresult[result_index] =
      DFTdataset(dfttype, scf_r,spin_type,
        orbitalStartIdx,orbitalNums);
      return scf_r;
    end
  end
  function set_current_dftdataset(scf_r,
    dfttype::DFTcommon.DFTtype,spin_type::SPINtype,result_index=1)
    if (DFTcommon.OpenMX == dfttype)
      # Read SCF and Set as current dftdata
      #scf_r = DFTforge.OpenMXdata.read_scf(scf_name);
      orbitalStartIdx = zeros(Int,scf_r.atomnum);
      orbitalIdx = 0; #각 atom별로 orbital index시작하는 지점
      orbitalNums = zeros(Int,scf_r.atomnum)
      for i = 1:scf_r.atomnum
          orbitalStartIdx[i] = orbitalIdx;
          orbitalIdx += scf_r.Total_NumOrbs[i]
          orbitalNums[i] = scf_r.Total_NumOrbs[i];
      end

      #assert(0 == scf_r.SpinP_switch - spin_type);
      dftresult[result_index] =
      DFTdataset(dfttype, scf_r,spin_type,
        orbitalStartIdx,orbitalNums);
      return scf_r;
    elseif (DFTcommon.Wannier90 == dfttype)
      orbitalStartIdx = zeros(Int,scf_r.atomnum);
      orbitalIdx = 0; #각 atom별로 orbital index시작하는 지점
      orbitalNums = zeros(Int,scf_r.atomnum)
      for i = 1:scf_r.atomnum
          orbitalStartIdx[i] = orbitalIdx;
          orbitalIdx += scf_r.Total_NumOrbs[i]
          orbitalNums[i] = scf_r.Total_NumOrbs[i];
      end
      #assert(0 == scf_r.SpinP_switch - spin_type);
      dftresult[result_index] =
      DFTdataset(dfttype, scf_r,spin_type,
        orbitalStartIdx,orbitalNums);
      return scf_r;
    end
  end
  function set_current_dftdataset(input)

    set_current_dftdataset(input[1],input[2],input[3],input[4])
  end
  function get_dftdataset(result_index=1)
    global dftresult
    return dftresult[result_index]
  end
  function cal_colinear_eigenstate(Kpoint::k_point_Tuple, spin_list=1,result_index=1)
    global dftresult
    return DFTforge.cal_colinear_eigenstate(Kpoint,
    dftresult[result_index].dfttype,dftresult[result_index].scf_r,spin_list )
  end
  function cal_nonco_linear_Eigenstate(Kpoint::k_point_Tuple,result_index=1)
    global dftresult
    return DFTforge.cal_nonco_linear_Eigenstate(Kpoint,
    dftresult[result_index].dfttype,dftresult[result_index].scf_r);
  end

  # for pmap
  function cal_eigenstate(input::Job_input_Type,result_index=1)
    # specfify spin type is required

    if (DFTforge.para_type == input.spin_type)
      kpoint_eigenstate_list = Array{Kpoint_eigenstate}();
      push!(kpoint_eigenstate_list,
      cal_colinear_eigenstate(input.k_point,1,input.result_index) );
      return kpoint_eigenstate_list
    elseif (DFTforge.colinear_type ==  input.spin_type)
      return cal_colinear_eigenstate(input.k_point,[1,2],input.result_index)
    elseif (DFTforge.non_colinear_type == input.spin_type)
      return cal_nonco_linear_Eigenstate(input.k_point,input.result_index)
    end
  end
  function cal_Hamiltonian(spin=1,result_index=1,Hmode::DFTcommon.nc_Hamiltonian_selection=DFTcommon.nc_allH)
    global dftresult;
    spin_type = dftresult[result_index].spin_type;
    dfttype::DFTtype = dftresult[result_index].dfttype;
    if (DFTforge.para_type == spin_type || DFTforge.colinear_type == spin_type)
      return cal_colinear_Hamiltonian(dfttype,dftresult[result_index].scf_r,spin);
    elseif (DFTforge.non_colinear_type == spin_type)
      return cal_noncolinear_Hamiltonian(dfttype,Hmode,dftresult[result_index].scf_r)
    end

  end


  function cachecal_all_Qpoint_eigenstats(q_point_list::Array{k_point_Tuple},
    hdf_cache_name,
    result_index=1,cache_index=1)
    global dftresult;

    Total_q_point_num = length(q_point_list)
    TotalOrbitalNum = get_TotalOrbitalNum(result_index);
    spin_type = dftresult[result_index].spin_type;
    spin_dim  = Int(spin_type)
    println(string("spin_dim ",spin_dim))


    TotalOrbitalNum2 = TotalOrbitalNum;
    if (DFTforge.non_colinear_type == spin_type)
      TotalOrbitalNum2 = 2*TotalOrbitalNum;
    end
    println((TotalOrbitalNum,TotalOrbitalNum2))

    #print(spin_dim )
    fid_hdf = h5open(hdf_cache_name,"w");
    hdf5_eigenstate_real = d_create(fid_hdf,"Eigenvect_real",datatype(Float64),
    dataspace(TotalOrbitalNum2,TotalOrbitalNum2, spin_dim, Total_q_point_num));

    hdf5_eigenstate_imag = d_create(fid_hdf,"Eigenvect_imag",datatype(Float64),
    dataspace(TotalOrbitalNum2,TotalOrbitalNum2, spin_dim, Total_q_point_num));

    hdf5_eigenvalues = d_create(fid_hdf,"Eigenvalues",datatype(Float64),
    dataspace(TotalOrbitalNum2, spin_dim, Total_q_point_num));

    hdf5_hamiltonian_real = d_create(fid_hdf,"Hamiltonian_real",datatype(Float64),
    dataspace(TotalOrbitalNum2,TotalOrbitalNum2, spin_dim));

    hdf5_hamiltonian_imag = d_create(fid_hdf,"Hamiltonian_imag",datatype(Float64),
    dataspace(TotalOrbitalNum2,TotalOrbitalNum2, spin_dim));
    # Write hamiltonian

    job_list = Array(Job_input_Type,0)
    q_points_int = Array{k_point_int_Tuple}(Total_q_point_num);
    q_points_intdic = Dict{k_point_int_Tuple,Int}();
    for (index,q) in enumerate(q_point_list)
      k_point = (q[1],q[2],q[3]);
      push!(job_list,Job_input_Type(k_point,spin_type,result_index));

      q_points_int[index] = k_point_float2int(kPoint2BrillouinZone_Tuple(k_point));
      q_points_intdic[q_points_int[index]] = index;
    end

    batch_size = 2*nprocs();
    cnt = 1;

    p = Progress(floor(Int, length(q_point_list)/batch_size),
    "Computing Eigenstates(q)...");

    while cnt <= Total_q_point_num
      # pmap
      start_idx = cnt;
      end_idx = minimum([cnt+batch_size-1,Total_q_point_num]);
      temp = pmap(cal_eigenstate,job_list[start_idx:end_idx]);
      ii = 1;
	    if (DFTforge.colinear_type == spin_type || DFTforge.para_type == spin_type)
        for jj = start_idx:end_idx

          hdf5_eigenstate_real[:,:,1,jj] = real(temp[ii][1].Eigenstate);
          hdf5_eigenstate_imag[:,:,1,jj] = imag(temp[ii][1].Eigenstate);
          hdf5_eigenvalues[:,1,jj] = temp[ii][1].Eigenvalues;
          if (DFTforge.colinear_type == spin_type )
            hdf5_eigenstate_real[:,:,2,jj] = real(temp[ii][2].Eigenstate);
            hdf5_eigenstate_imag[:,:,2,jj] = imag(temp[ii][2].Eigenstate);
            hdf5_eigenvalues[:,2,jj] = temp[ii][2].Eigenvalues;
          end
          ii += 1;
	      end
      elseif (DFTforge.non_colinear_type == spin_type)
        for jj = start_idx:end_idx
          hdf5_eigenstate_real[:,:,1,jj] = real(temp[ii].Eigenstate);
          hdf5_eigenstate_imag[:,:,1,jj] = imag(temp[ii].Eigenstate);
          hdf5_eigenvalues[:,1,jj] = temp[ii].Eigenvalues;
          ii += 1;
        end
      end
      # write to hdf5
      cnt = end_idx + 1;
      next!(p)
    end

    H::Hamiltonian_type = cal_Hamiltonian(1,result_index);

    println(size(H))
    hdf5_hamiltonian_real[:,:,1] = real(H);
    hdf5_hamiltonian_imag[:,:,1] = imag(H);
    if (DFTforge.colinear_type == spin_type )
      H = cal_Hamiltonian(2,result_index);
      hdf5_hamiltonian_real[:,:,2] = real(H);
      hdf5_hamiltonian_imag[:,:,2] = imag(H);
    end
    flush(fid_hdf);
    close(fid_hdf);
    #fid_hdf = h5open(hdf_cache_name,"r");
    #q_points_int = Array{k_point_int_Tuple}(Total_q_point_num);


    eigenstate_cache =  Eigenstate_hdf5(hdf_cache_name,fid_hdf,q_point_list,
      q_points_int,q_points_intdic,
      dftresult[result_index].spin_type,TotalOrbitalNum,
      dftresult[result_index]);

    eigenstate_list[cache_index] = eigenstate_cache;
    hdf5_eigenstate_real = [];
    hdf5_eigenstate_imag = [];
    hdf5_eigenvalues = [];
    gc();
    return eigenstate_cache;
  end

  function cacheset(eigenstate_cache::Eigenstate_hdf5,cache_index=1)
    global eigenstate_list
    fid_hdf = h5open(eigenstate_cache.hdf_cache_name,"r");

    eigenstate_list[cache_index] = eigenstate_cache;

    eigenstate_list[cache_index].Eigenvect_real = readmmap(fid_hdf["Eigenvect_real"]);
    eigenstate_list[cache_index].Eigenvect_imag = readmmap(fid_hdf["Eigenvect_imag"]);
    eigenstate_list[cache_index].Eigenvalues    = readmmap(fid_hdf["Eigenvalues"]);
    eigenstate_list[cache_index].Hamiltonian_real    = readmmap(fid_hdf["Hamiltonian_real"]);
    eigenstate_list[cache_index].Hamiltonian_imag    = readmmap(fid_hdf["Hamiltonian_imag"]);

    #eigenstate_cache.fid_hdf = fid_hdf;

    gc();
  end
  function cacheread(cache_index=1)
    global eigenstate_list;
    return eigenstate_list[cache_index]
  end
  function cacheread_eigenstate(k_point::k_point_Tuple,spin,cache_index=1)
    global eigenstate_list;
    global orbital_mask1,orbital_mask2,orbital_mask_on

    k_point_int = k_point_float2int(kPoint2BrillouinZone_Tuple(k_point));

    # non spin :
    # collinear spin :
    # * spin = 1 down spin
    # * spin = 2 up spin
    # non-collienar spin :
    spin_type = eigenstate_list[cache_index].spin_type
    q_index = -1;
    TotalOrbitalNum = eigenstate_list[cache_index].TotalOrbitalNum;
    TotalOrbitalNum2 = TotalOrbitalNum
    if (DFTforge.non_colinear_type == spin_type)
      TotalOrbitalNum2 = 2*TotalOrbitalNum;
    end
    Eigenstate = zeros(Complex_my,TotalOrbitalNum2,TotalOrbitalNum2);
    Eigenvalues = zeros(Float_my,TotalOrbitalNum2);

    if (haskey(eigenstate_list[cache_index].q_points_intdic, k_point_int))
      q_index =  eigenstate_list[cache_index].q_points_intdic[k_point_int];
    else
      error(string("cacheread_eigenstate ",k_point," Not Found Exception"))
    end
    if (q_index>=1)
      if (DFTforge.para_type == spin_type)
          Eigenstate[:,:] = eigenstate_list[cache_index].Eigenvect_real[:,:,1,q_index] +
          im * eigenstate_list[cache_index].Eigenvect_imag[:,:,1,q_index];
          Eigenvalues[:] = eigenstate_list[cache_index].Eigenvalues[:,1,q_index];
      elseif (DFTforge.colinear_type == spin_type)
          Eigenstate[:,:] = eigenstate_list[cache_index].Eigenvect_real[:,:,spin,q_index] +
          im * eigenstate_list[cache_index].Eigenvect_imag[:,:,spin,q_index];
          Eigenvalues[:] = eigenstate_list[cache_index].Eigenvalues[:,spin,q_index];
      elseif (DFTforge.non_colinear_type == spin_type)
        Eigenstate[:,:] = eigenstate_list[cache_index].Eigenvect_real[:,:,1,q_index] +
        im * eigenstate_list[cache_index].Eigenvect_imag[:,:,1,q_index];
        Eigenvalues[:] = eigenstate_list[cache_index].Eigenvalues[:,1,q_index];
      end
    end

    return Kpoint_eigenstate(Eigenstate,Eigenvalues,k_point);
  end
  function cacheread_Hamiltonian(spin=1,cache_index=1)
    global eigenstate_list;
    TotalOrbitalNum = eigenstate_list[cache_index].TotalOrbitalNum;
    TotalOrbitalNum2 = TotalOrbitalNum
    spin_type =  eigenstate_list[cache_index].spin_type;
    if (DFTforge.non_colinear_type == spin_type)
      TotalOrbitalNum2 = 2*TotalOrbitalNum;
    end
    Hamiltonian = zeros(Complex_my,TotalOrbitalNum2,TotalOrbitalNum2);
    Hamiltonian[:,:] = eigenstate_list[cache_index].Hamiltonian_real[:,:,spin] +
      im * eigenstate_list[cache_index].Hamiltonian_imag[:,:,spin]

    return Hamiltonian;
  end
  function cacheread_atomsOrbital_lists(cache_index=1)
    global eigenstate_list;
    return (eigenstate_list[cache_index].dftresult.orbitalStartIdx,
    eigenstate_list[cache_index].dftresult.orbitalNums)
  end
  function cacheread_lampup(q_point_list::Array{k_point_Tuple},cache_index=1)
    for q_point in q_point_list
      cacheread_eigenstate(q_point,cache_index)
    end
  end

  function get_TotalOrbitalNum(result_index=1)
    global dftresult
    TotalOrbitalNum = sum(dftresult[result_index].scf_r.Total_NumOrbs);
    #print(TotalOrbitalNum)
    #print(typeof(TotalOrbitalNum))
    assert(Int == typeof(TotalOrbitalNum));
    return TotalOrbitalNum
  end
  function get_ChempP(result_index=1)
    global dftresult
    ChemP::Float64 = dftresult[result_index].scf_r.ChemP;
    return ChemP;
  end

  function cacheread(cache_name::AbstractString)

  end
  function clear_cache(cache_index=1)

  end

  function clear_dftdataset()
    dftresult = Array{DFTdataset}();
  end

################################################################################
# Genernal K space, K-Q space
################################################################################
  function Qspace_Ksum_parallel(kq_function,q_point_list,k_point_list,
    result_index=1,cache_index=1)
    batch_size = 2*nprocs();
    cnt = 1
    spin_type = get_dftdataset(result_index).spin_type;

    Q_ksum = Dict{k_point_int_Tuple,Array{Complex_my,1}}();

    p = Progress( round(Int, length(q_point_list)/4),
      string("Computing  (Q:",length(q_point_list),", K:",length(k_point_list),")...") );
    for (q_i,q_point) in enumerate(q_point_list)
      q_point_int = k_point_float2int(q_point);

      job_list = Array(Job_input_kq_Type,0)
      for k_point in k_point_list
        kq_point = (q_point[1] + k_point[1],q_point[2] + k_point[2],q_point[3] + k_point[3]) ;
        kq_point = kPoint2BrillouinZone_Tuple(kq_point);
        kq_point_int = k_point_float2int(kq_point);
        push!(job_list,Job_input_kq_Type(k_point,kq_point,spin_type));
      end
      temp = pmap(kq_function,job_list);

      Q_ksum[q_point_int] = vcat(temp...);
      ## End of each q_point
      if (1==rem(q_i,4))
        next!(p)
      end
      if (1==rem(q_i,50))
        @everywhere gc()
      end
    end
    return Q_ksum;
  end
  function Qspace_Ksum_atom_parallel(kq_function,q_point_list,k_point_list,atom1,atom2,
    result_index=1,cache_index=1)
    batch_size = 2*nprocs();
    cnt = 1
    spin_type = get_dftdataset(result_index).spin_type;

    Q_ksum = Dict{k_point_int_Tuple,Array{Complex_my,1}}();

    p = Progress( round(Int, length(q_point_list)/6),
      string("Computing  (Q:",length(q_point_list),", K:",length(k_point_list),")...") );
    for (q_i,q_point) in enumerate(q_point_list)
      q_point_int = k_point_float2int(q_point);

      job_list = Array(Job_input_kq_atom_Type,0)
      for k_point in k_point_list
        kq_point = (q_point[1] + k_point[1],q_point[2] + k_point[2],q_point[3] + k_point[3]) ;
        kq_point = kPoint2BrillouinZone_Tuple(kq_point);
        kq_point_int = k_point_float2int(kq_point);
        push!(job_list,Job_input_kq_atom_Type(k_point,kq_point,spin_type,atom1,atom2));
      end
      temp = pmap(kq_function,job_list);

      Q_ksum[q_point_int] = vcat(temp...);
      ## End of each q_point
      if (1==rem(q_i,6))
        next!(p)
      end
      if (1==rem(q_i,50))
        @everywhere gc()
      end
    end
    return Q_ksum;
  end
  function Qspace_Ksum_atomlist_parallel_nc(kq_function,q_point_list,k_point_list,
    atom12_list::Vector{Tuple{Int64,Int64}},
    result_index=1,cache_index=1)

    batch_size = 2*nprocs();
    cnt = 1
    spin_type = get_dftdataset(result_index).spin_type;

    Xij_Q = Array(Dict{k_point_int_Tuple,Array{Complex_my,1}},10,length(atom12_list));
    Xij_Q_mean = Array(Dict{k_point_int_Tuple,Complex_my},10,length(atom12_list));
    for xyz_ij = 1:10
        for atom12_i = 1:length(atom12_list)
            Xij_Q[xyz_ij,atom12_i] = Dict{k_point_int_Tuple,Array{Complex_my,1}}();
            Xij_Q_mean[xyz_ij,atom12_i] = Dict{k_point_int_Tuple,Complex_my}();
        end
    end
    #Q_ksum = Dict{k_point_int_Tuple,Array{Complex_my,1}}();

    p = Progress( round(Int, length(q_point_list)/6),
      string("Computing  (Q:",length(q_point_list),", K:",length(k_point_list),")...") );
    for (q_i,q_point) in enumerate(q_point_list)
      q_point_int = k_point_float2int(q_point);

      job_list = Array(Job_input_kq_atom_list_Type,0)
      for k_point in k_point_list
        kq_point = (q_point[1] + k_point[1],q_point[2] + k_point[2],q_point[3] + k_point[3]) ;
        kq_point = kPoint2BrillouinZone_Tuple(kq_point);
        kq_point_int = k_point_float2int(kq_point);
        push!(job_list,Job_input_kq_atom_list_Type(k_point,kq_point,spin_type,atom12_list));
      end
      X_temp = pmap(kq_function,job_list);
      for xyz_ij = 1:10
        for atom12_i = 1:length(atom12_list)
          tmp = zeros(Complex_my,length(k_point_list))
          for ii = 1:length(k_point_list)
              tmp[ii] = X_temp[ii][xyz_ij,atom12_i];
          end
          Xij_Q[xyz_ij,atom12_i][q_point_int] = copy(tmp);
          Xij_Q_mean[xyz_ij,atom12_i][q_point_int] = mean(tmp);
        end
      end

      #Q_ksum[q_point_int] = vcat(temp...);
      ## End of each q_point
      if (1==rem(q_i,6))
        next!(p)
      end
      if (1==rem(q_i,50))
        @everywhere gc()
      end
    end
    return (Xij_Q,Xij_Q_mean);
  end

  end


end
