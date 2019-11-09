###############################################################################
# Hongkee Yoon Hongkeeyoon@kaist.ac.kr
# 2019.05
# https://kaist-elst.github.io/DFTforge.jl/
###############################################################################


__precompile__(true)
using DFTforge
using ..DFTcommon
using HDF5
using ProgressMeter
using Distributed

# Julia 1.0 support
using LinearAlgebra
using Mmap
using Statistics
#
include("resultExport.jl")
export export2mat_K_Q,export2mat_K_Q_nc

export set_current_dftdataset,cal_colinear_eigenstate,cal_nonco_linear_Eigenstate,
cal_Hamiltonian,get_dftdataset
export Job_input_Type,Job_input_kq_Type,Job_input_kq_atom_Type,Job_input_kq_atom_list_Type

export cachecal_all_Qpoint_eigenstats,cacheset,cacheread_eigenstate,cacheread,
cacheread_lampup,cacheread_atomsOrbital_lists,cacheread_Hamiltonian
export cachecal_all_Qpoint_eigenstats_as_nc,cacheread_eigenstate_as_nc,
cacheread_lampup_as_nc,cacheread_atomsOrbital_lists,cacheread_Hamiltonian_as_nc
export get_ChempP
export Qspace_Ksum_parallel,Qspace_Ksum_atom_parallel,
Kspace_parallel,Qspace_Ksum_atomlist_parallel,Qspace_Ksum_atomlist_parallel_nc

#=
type DFTdataset
  dfttype::DFTcommon.DFTtype
  scf_r
  spin_type::SPINtype
  orbitalStartIdx::Array{Int,1}
  orbitalNums::Array{Int,1}
end
=#

mutable struct Eigenstate_hdf5
  hdf_cache_name::AbstractString
  fid_hdf::HDF5.HDF5File
  q_points::Array{k_point_Tuple}
  q_points_int::Array{k_point_int_Tuple};
  q_points_intdic::Dict{k_point_int_Tuple,Int};
  spin_type::SPINtype
  TotalOrbitalNum::Int
  dftresult::Hamiltonian_info_type
  Energy_idx_num::Array{Int64,2}
  Eigenvect_real::Array{Float64,4}
  Eigenvect_imag::Array{Float64,4}
  Eigenvalues::Array{Float64,3}
  Hamiltonian_real::Array{Float64,4}
  Hamiltonian_imag::Array{Float64,4}


  Eigenstate_hdf5(hdf_cache_name,fid_hdf,q_points,q_points_int,
    q_points_intdic,spin_type,TotalOrbitalNum,dftresult) =
  new(hdf_cache_name,fid_hdf,q_points,q_points_int,
    q_points_intdic,spin_type,TotalOrbitalNum,dftresult);
end

struct Job_input_Type
  k_point::k_point_Tuple
  spin_type::SPINtype
  Hmode::DFTcommon.nc_Hamiltonian_selection
  result_index::Int
  Job_input_Type(k_point, spin_type) = new(k_point, spin_type, DFTcommon.nc_allH, 1)
  Job_input_Type(k_point, spin_type, result_index) = new(k_point, spin_type, DFTcommon.nc_allH, result_index)

  Job_input_Type(k_point, spin_type, Hmode::DFTcommon.nc_Hamiltonian_selection) = new(k_point, spin_type, Hmode, 1)
  Job_input_Type(k_point, spin_type, Hmode::DFTcommon.nc_Hamiltonian_selection, result_index) = new(k_point, spin_type, Hmode, result_index)
end

struct Job_input_kq_Type
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
struct Job_input_kq_atom_Type
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

struct Job_input_kq_atom_list_Type
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

global dftresult = Dict{Int64, Hamiltonian_info_type}();
global eigenstate_list =  Dict{Int64, Eigenstate_hdf5}(); #cached Eigenstates

function set_current_dftdataset(scf_name::AbstractString,result_file_dict::Dict{AbstractString,AbstractString},
  dfttype::DFTcommon.DFTtype,spin_type::SPINtype,
    basisTransform_rule::basisTransform_rule_type=basisTransform_rule_type(),result_index=1)
  global dftresult;
  if (DFTcommon.OpenMX == dfttype || DFTcommon.EcalJ == dfttype)
    # Read SCF and Set as current dftdata
    #scf_r = DFTforge.OpenMXdata.read_scf(scf_name);
    hamiltonian_info = read_dftresult(scf_name,result_file_dict,dfttype,spin_type,basisTransform_rule)
    #=
    orbitalStartIdx = zeros(Int,scf_r.atomnum);
    orbitalIdx::Int = 0; #각 atom별로 orbital index시작하는 지점
    orbitalNums = zeros(Int,scf_r.atomnum)
    for i = 1:scf_r.atomnum
        orbitalStartIdx[i] = orbitalIdx;
        orbitalIdx += scf_r.Total_NumOrbs[i]
        orbitalNums[i] = scf_r.Total_NumOrbs[i];
    end

    #@assert(0 == scf_r.SpinP_switch - spin_type);
    dftresult[result_index] =
    DFTdataset(dfttype, scf_r,spin_type,
      orbitalStartIdx,orbitalNums);

    =#
    dftresult[result_index] = hamiltonian_info;
    return hamiltonian_info;
  end
end
function set_current_dftdataset(wannier_fname::AbstractString,result_file_dict::Dict{AbstractString,AbstractString},
  dfttype::DFTcommon.DFTtype,
    Wannier90_type::DFTcommon.Wannier90type,spin_type::DFTcommon.SPINtype,
      atoms_orbitals_list::Vector{Array{Int64}},
    atomnum::Int,atompos::Array{Float64,2},basisTransform_rule::basisTransform_rule_type=basisTransform_rule_type(),
    result_index=1)
  global dftresult;
  if (DFTcommon.Wannier90 ==dfttype)
    hamiltonian_info = read_dftresult(wannier_fname, result_file_dict, dfttype,
        Wannier90_type,spin_type,
          atoms_orbitals_list,
        atomnum,atompos,basisTransform_rule)
    dftresult[result_index] = hamiltonian_info;
    return hamiltonian_info;
  end
end
function set_current_dftdataset(scf_r,
  dfttype::DFTcommon.DFTtype,spin_type::DFTcommon.SPINtype,result_index=1)
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

    #@assert(0 == scf_r.SpinP_switch - spin_type);
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
    #@assert(0 == scf_r.SpinP_switch - spin_type);
    dftresult[result_index] =
    DFTdataset(dfttype, scf_r,spin_type,
      orbitalStartIdx,orbitalNums);
    return scf_r;
  end
end
function set_current_dftdataset(input)
  hamiltonian_info::Hamiltonian_info_type = input[1]
  result_index::Int = input[2]

  dftresult[result_index] = hamiltonian_info;
  #set_current_dftdataset(input[1],input[2],input[3],input[4])
end
function get_dftdataset(result_index=1)
  global dftresult
  return dftresult[result_index]
end
function cal_colinear_eigenstate(Kpoint::k_point_Tuple, spin_list=1,result_index=1)
  global dftresult
  return DFTforge.cal_colinear_eigenstate(Kpoint,dftresult[result_index],spin_list)
end
function cal_colinear_eigenstate_as_nc(Kpoint::k_point_Tuple, result_index=1)
  global dftresult
  return DFTforge.cal_colinear_eigenstate_as_nc(Kpoint,dftresult[result_index])
end
function cal_nonco_linear_Eigenstate(Kpoint::k_point_Tuple,result_index=1)
  global dftresult
  return DFTforge.cal_nonco_linear_Eigenstate(Kpoint,
  dftresult[result_index].dfttype,dftresult[result_index]);
end


## for pmap functions
# for pmap
function cal_eigenstate(input::Job_input_Type,result_index=1)
  # specfify spin type is required

  if (DFTcommon.para_type == input.spin_type)
    return cal_colinear_eigenstate(input.k_point,[1],input.result_index)
  elseif (DFTcommon.colinear_type ==  input.spin_type)
    return cal_colinear_eigenstate(input.k_point,[1,2],input.result_index)
  elseif (DFTcommon.non_colinear_type == input.spin_type)
    return cal_nonco_linear_Eigenstate(input.k_point,input.result_index)
  end
end

function cal_eigenstate_as_nc(input::Job_input_Type,result_index=1)
  # specfify spin type is required

  if (DFTcommon.para_type == input.spin_type)
    #=
    kpoint_eigenstate_list = Array{Kpoint_eigenstate}();
    push!(kpoint_eigenstate_list,
    cal_colinear_eigenstate(input.k_point,1,input.result_index) );
    return kpoint_eigenstate_list
    =#
  elseif (DFTcommon.colinear_type ==  input.spin_type)
    return cal_colinear_eigenstate_as_nc(input.k_point,input.result_index)
  elseif (DFTcommon.non_colinear_type == input.spin_type)
    return cal_nonco_linear_Eigenstate(input.k_point,input.result_index)
  end
end

#for pmap
function cal_noncolinear_Hamiltonian(input::Job_input_Type)
  @assert(DFTcommon.non_colinear_type == input.spin_type)
  global dftresult
  result_index = input.result_index;
  return DFTforge.cal_noncolinear_Hamiltonian(input.k_point,
  dftresult[result_index].dfttype, input.Hmode, dftresult[result_index]);
end


function cachecal_all_Qpoint_eigenstats(q_point_list::Array{k_point_Tuple},
  hdf_cache_name,
  result_index=1,cache_index=1)
  global dftresult;

  Total_q_point_num = length(q_point_list)
  TotalOrbitalNum = get_TotalOrbitalNum(result_index);
  spin_type = dftresult[result_index].spin_type;

  spin_dim  = 1 # Int(spin_type)
  TotalOrbitalNum2 = TotalOrbitalNum;
  if (DFTcommon.non_colinear_type == spin_type)
    TotalOrbitalNum2 = 2*TotalOrbitalNum;
    spin_dim = 3;
  elseif (DFTcommon.colinear_type == spin_type)
    spin_dim = 2;
  elseif (DFTcommon.para_type == spin_type)
    spin_dim = 1;
  end
  println(string("spin_dim ",spin_dim))
  println((TotalOrbitalNum,TotalOrbitalNum2))

  #print(spin_dim )
  fid_hdf = h5open(hdf_cache_name,"w");
  hdf5_energy_idx_num = d_create(fid_hdf,"Energy_idx_num",datatype(Int64),
  dataspace(spin_dim, Total_q_point_num));

  hdf5_eigenstate_real = d_create(fid_hdf,"Eigenvect_real",datatype(Float64),
  dataspace(TotalOrbitalNum2,TotalOrbitalNum2, spin_dim, Total_q_point_num));

  hdf5_eigenstate_imag = d_create(fid_hdf,"Eigenvect_imag",datatype(Float64),
  dataspace(TotalOrbitalNum2,TotalOrbitalNum2, spin_dim, Total_q_point_num));

  hdf5_eigenvalues = d_create(fid_hdf,"Eigenvalues",datatype(Float64),
  dataspace(TotalOrbitalNum2, spin_dim, Total_q_point_num));

  hdf5_hamiltonian_real = d_create(fid_hdf,"Hamiltonian_real",datatype(Float64),
  dataspace(TotalOrbitalNum2,TotalOrbitalNum2, spin_dim, Total_q_point_num));

  hdf5_hamiltonian_imag = d_create(fid_hdf,"Hamiltonian_imag",datatype(Float64),
  dataspace(TotalOrbitalNum2,TotalOrbitalNum2, spin_dim,Total_q_point_num ));


  # Write hamiltonian

  job_list = Array{Job_input_Type}(undef,0)
  q_points_int = Array{k_point_int_Tuple}(undef,Total_q_point_num);
  q_points_intdic = Dict{k_point_int_Tuple,Int}();
  for (index,q) in enumerate(q_point_list)
    k_point = (q[1],q[2],q[3]);
    push!(job_list,Job_input_Type(k_point,spin_type,result_index));

    q_points_int[index] = k_point_float2int(kPoint2BrillouinZone_Tuple(k_point));
    q_points_intdic[q_points_int[index]] = index;
  end

  batch_size = 5*nprocs();
  cnt = 1;

  p = Progress(ceil(Int, 1.0+length(q_point_list)/batch_size),
  "Computing Eigenstates(q)...");
  p.barglyphs=BarGlyphs("[=> ]")
  p.output = stdout
  while cnt <= Total_q_point_num
    # pmap
    start_idx = cnt;
    end_idx = minimum([cnt+batch_size-1,Total_q_point_num]);

    temp = pmap(cal_eigenstate,job_list[start_idx:end_idx]);

    ii = 1;
    if (DFTcommon.colinear_type == spin_type || DFTcommon.para_type == spin_type)
      for jj = start_idx:end_idx
        # energy_idx_num could be smaller then TotalOrbitalNum2 when overlap is small
        energy_idx_num = length(temp[ii][1].Eigenvalues)

        hdf5_eigenstate_real[:,1:energy_idx_num,1,jj] = real(temp[ii][1].Eigenstate);
        hdf5_eigenstate_imag[:,1:energy_idx_num,1,jj] = imag(temp[ii][1].Eigenstate);
        hdf5_eigenvalues[1:energy_idx_num,1,jj] = temp[ii][1].Eigenvalues;

        hdf5_hamiltonian_real[:,:,1,jj] = real(temp[ii][1].Hamiltonian);
        hdf5_hamiltonian_imag[:,:,1,jj] = imag(temp[ii][1].Hamiltonian);

        hdf5_energy_idx_num[1,jj] = energy_idx_num;

        if (DFTcommon.colinear_type == spin_type )
          energy_idx_num = length(temp[ii][2].Eigenvalues)

          hdf5_eigenstate_real[:,1:energy_idx_num,2,jj] = real(temp[ii][2].Eigenstate);
          hdf5_eigenstate_imag[:,1:energy_idx_num,2,jj] = imag(temp[ii][2].Eigenstate);
          hdf5_eigenvalues[1:energy_idx_num,2,jj] = temp[ii][2].Eigenvalues;

          hdf5_hamiltonian_real[:,:,2,jj] = real(temp[ii][2].Hamiltonian);
          hdf5_hamiltonian_imag[:,:,2,jj] = imag(temp[ii][2].Hamiltonian);

          hdf5_energy_idx_num[2,jj] = energy_idx_num;
        end

        ii += 1;
      end
    elseif (DFTcommon.non_colinear_type == spin_type)
      for jj = start_idx:end_idx
        energy_idx_num = length(temp[ii][1].Eigenvalues)

        hdf5_eigenstate_real[:,1:energy_idx_num,1,jj] = real(temp[ii].Eigenstate);
        hdf5_eigenstate_imag[:,1:energy_idx_num,1,jj] = imag(temp[ii].Eigenstate);
        hdf5_eigenvalues[1:energy_idx_num,1,jj] = temp[ii].Eigenvalues;

        hdf5_hamiltonian_real[:,:,1,jj] = real(temp[ii].Hamiltonian);
        hdf5_hamiltonian_imag[:,:,1,jj] = imag(temp[ii].Hamiltonian);

        ii += 1;
      end
    end

    #for non_colinear_type spin
    # H of (real only, image only)
    if (DFTcommon.non_colinear_type == spin_type)
      job_list_nc_realH_only = Array{Job_input_Type}(undef,0)
      job_list_nc_imagH_only = Array{Job_input_Type}(undef,0)
      for (index,job_item) in enumerate(job_list[start_idx:end_idx])
        push!(job_list_nc_realH_only,Job_input_Type(job_item.k_point ,spin_type, DFTcommon.nc_realH_only, result_index));
        push!(job_list_nc_imagH_only,Job_input_Type(job_item.k_point ,spin_type, DFTcommon.nc_imagH_only, result_index));
      end
      tmp_realH_only = pmap(cal_noncolinear_Hamiltonian, job_list_nc_realH_only);
      tmp_imagH_only = pmap(cal_noncolinear_Hamiltonian, job_list_nc_imagH_only);

      ii = 1;
      for jj = start_idx:end_idx
        energy_idx_num = length(temp[ii][1].Eigenvalues)

        hdf5_hamiltonian_real[:,1:energy_idx_num,2,jj] = real(tmp_realH_only[ii]);
        hdf5_hamiltonian_imag[:,1:energy_idx_num,2,jj] = imag(tmp_realH_only[ii]);

        hdf5_hamiltonian_real[:,1:energy_idx_num,3,jj] = real(tmp_imagH_only[ii]);
        hdf5_hamiltonian_imag[:,1:energy_idx_num,3,jj] = imag(tmp_imagH_only[ii]);

        ii += 1;
      end
    end

    # write to hdf5
    cnt = end_idx + 1;
    next!(p)
  end
  next!(p)
  #=
  #H::Hamiltonian_type = cal_Hamiltonian(1,result_index);

  #println(size(H))
  hdf5_hamiltonian_real[:,:,1] = real(H);
  hdf5_hamiltonian_imag[:,:,1] = imag(H);
  if (DFTcommon.colinear_type == spin_type )
    H = cal_Hamiltonian(2,result_index);
    hdf5_hamiltonian_real[:,:,2] = real(H);
    hdf5_hamiltonian_imag[:,:,2] = imag(H);
  end
  =#
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
  GC.gc();
  return eigenstate_cache;
end

function cachecal_all_Qpoint_eigenstats_as_nc(q_point_list::Array{k_point_Tuple},
  hdf_cache_name,
  result_index=1,cache_index=1)
  println(" Deprecated function ")
  exit(0)
  global dftresult;

  Total_q_point_num = length(q_point_list)
  TotalOrbitalNum = get_TotalOrbitalNum(result_index);
  spin_type = dftresult[result_index].spin_type;
  spin_dim  = Int(spin_type)
  println(string("spin_dim ",spin_dim))
  if (DFTcommon.para_type == spin_type)
      println(" Para Hamiltonian will not be treated as nc.")
      println(" Try to use 'cachecal_all_Qpoint_eigenstats' instead. ")
      @assert(false)
  end

  TotalOrbitalNum2 = TotalOrbitalNum;
  TotalOrbitalNum3 = TotalOrbitalNum*2;
  if (DFTcommon.non_colinear_type == spin_type)
    TotalOrbitalNum2 = 2*TotalOrbitalNum;
  end
  println((TotalOrbitalNum,TotalOrbitalNum2,TotalOrbitalNum3))

  #print(spin_dim )
  fid_hdf = h5open(hdf_cache_name,"w");
  hdf5_eigenstate_real = d_create(fid_hdf,"Eigenvect_real",datatype(Float64),
  dataspace(TotalOrbitalNum3,TotalOrbitalNum3, spin_dim, Total_q_point_num));

  hdf5_eigenstate_imag = d_create(fid_hdf,"Eigenvect_imag",datatype(Float64),
  dataspace(TotalOrbitalNum3,TotalOrbitalNum3, spin_dim, Total_q_point_num));

  hdf5_eigenvalues = d_create(fid_hdf,"Eigenvalues",datatype(Float64),
  dataspace(TotalOrbitalNum3, spin_dim, Total_q_point_num));

  hdf5_hamiltonian_real = d_create(fid_hdf,"Hamiltonian_real",datatype(Float64),
  dataspace(TotalOrbitalNum3,TotalOrbitalNum3, spin_dim, Total_q_point_num));

  hdf5_hamiltonian_imag = d_create(fid_hdf,"Hamiltonian_imag",datatype(Float64),
  dataspace(TotalOrbitalNum3,TotalOrbitalNum3, spin_dim,Total_q_point_num ));
  # Write hamiltonian

  job_list = Array{Job_input_Type}(undef,0)
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

  p = Progress(ceil(Int, 1.0+length(q_point_list)/batch_size),
  "Computing Eigenstates as nc(q)...");
  p.barglyphs=BarGlyphs("[=> ]")
  p.output = stdout
  while cnt <= Total_q_point_num
    # pmap
    start_idx = cnt;
    end_idx = minimum([cnt+batch_size-1,Total_q_point_num]);

    temp = pmap(cal_eigenstate_as_nc,job_list[start_idx:end_idx]);

    ii = 1;
    if ((DFTcommon.colinear_type == spin_type) ||  (DFTcommon.non_colinear_type == spin_type))
      for jj = start_idx:end_idx
        hdf5_eigenstate_real[:,:,1,jj] = real(temp[ii].Eigenstate);
        hdf5_eigenstate_imag[:,:,1,jj] = imag(temp[ii].Eigenstate);
        hdf5_eigenvalues[:,1,jj] = temp[ii].Eigenvalues;

        hdf5_hamiltonian_real[:,:,1,jj] = real(temp[ii].Hamiltonian);
        hdf5_hamiltonian_imag[:,:,1,jj] = imag(temp[ii].Hamiltonian);

        ii += 1;
      end
    end

    #for non_colinear_type spin
    # H of (real only, image only)
    if (DFTcommon.non_colinear_type == spin_type)
      job_list_nc_realH_only = Array{Job_input_Type}(undef,0)
      job_list_nc_imagH_only = Array{Job_input_Type}(undef,0)
      for (index,job_item) in enumerate(job_list[start_idx:end_idx])
        push!(job_list_nc_realH_only,Job_input_Type(job_item.k_point ,spin_type, DFTcommon.nc_realH_only, result_index));
        push!(job_list_nc_imagH_only,Job_input_Type(job_item.k_point ,spin_type, DFTcommon.nc_imagH_only, result_index));
      end
      tmp_realH_only = pmap(cal_noncolinear_Hamiltonian, job_list_nc_realH_only);
      tmp_imagH_only = pmap(cal_noncolinear_Hamiltonian, job_list_nc_imagH_only);

      ii = 1;
      for jj = start_idx:end_idx

        hdf5_hamiltonian_real[:,:,2,jj] = real(tmp_realH_only[ii]);
        hdf5_hamiltonian_imag[:,:,2,jj] = imag(tmp_realH_only[ii]);

        hdf5_hamiltonian_real[:,:,3,jj] = real(tmp_imagH_only[ii]);
        hdf5_hamiltonian_imag[:,:,3,jj] = imag(tmp_imagH_only[ii]);

        ii += 1;
      end
    end

    # write to hdf5
    cnt = end_idx + 1;
    next!(p)
  end
  next!(p)
  #=
  #H::Hamiltonian_type = cal_Hamiltonian(1,result_index);

  #println(size(H))
  hdf5_hamiltonian_real[:,:,1] = real(H);
  hdf5_hamiltonian_imag[:,:,1] = imag(H);
  if (DFTcommon.colinear_type == spin_type )
    H = cal_Hamiltonian(2,result_index);
    hdf5_hamiltonian_real[:,:,2] = real(H);
    hdf5_hamiltonian_imag[:,:,2] = imag(H);
  end
  =#
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
  GC.gc();
  return eigenstate_cache;
end

function cacheset(eigenstate_cache::Eigenstate_hdf5,cache_index=1)
  global eigenstate_list
  fid_hdf = h5open(eigenstate_cache.hdf_cache_name,"r");

  eigenstate_list[cache_index] = eigenstate_cache;

  eigenstate_list[cache_index].Energy_idx_num = readmmap(fid_hdf["Energy_idx_num"]);
  eigenstate_list[cache_index].Eigenvect_real = readmmap(fid_hdf["Eigenvect_real"]);
  eigenstate_list[cache_index].Eigenvect_imag = readmmap(fid_hdf["Eigenvect_imag"]);
  eigenstate_list[cache_index].Eigenvalues    = readmmap(fid_hdf["Eigenvalues"]);
  eigenstate_list[cache_index].Hamiltonian_real    = readmmap(fid_hdf["Hamiltonian_real"]);
  eigenstate_list[cache_index].Hamiltonian_imag    = readmmap(fid_hdf["Hamiltonian_imag"]);

  #eigenstate_cache.fid_hdf = fid_hdf;

  GC.gc();
end
function cacheread(cache_index=1)
  global eigenstate_list;
  return eigenstate_list[cache_index]
end
function cacheread_kpoint2q_index(k_point::k_point_Tuple,cache_index=1)
  global eigenstate_list;
  k_point_int = k_point_float2int(kPoint2BrillouinZone_Tuple(k_point));
  q_index = -1;
  if (haskey(eigenstate_list[cache_index].q_points_intdic, k_point_int))
    q_index =  eigenstate_list[cache_index].q_points_intdic[k_point_int];
  else
    error(string("cacheread_eigenstate ",k_point," Not Found Exception"))
  end
  return q_index;
end
function cacheread_eigenstate(k_point::k_point_Tuple,spin,cache_index=1)
  global eigenstate_list;
  global orbital_selection1,orbital_selection2,orbital_selection_on

  #k_point_int = k_point_float2int(kPoint2BrillouinZone_Tuple(k_point));

  # non spin :
  # collinear spin :
  # * spin = 1 down spin
  # * spin = 2 up spin
  # non-collienar spin :
  spin_type = eigenstate_list[cache_index].spin_type
  q_index = cacheread_kpoint2q_index(k_point,cache_index);
  TotalOrbitalNum = eigenstate_list[cache_index].TotalOrbitalNum;
  TotalOrbitalNum2 = TotalOrbitalNum
  if (DFTcommon.non_colinear_type == spin_type)
    TotalOrbitalNum2 = 2*TotalOrbitalNum;
  end

  # energy_idx_num could be smaller then TotalOrbitalNum2 when overlap is small
  energy_idx_num = TotalOrbitalNum2 #length(temp[ii][1].Eigenvalues)
  if (1 <= q_index)
    energy_idx_num = eigenstate_list[cache_index].Energy_idx_num[1,q_index];
  end

  Eigenstate = zeros(Complex_my, TotalOrbitalNum2, energy_idx_num);
  Eigenvalues = zeros(Float_my, energy_idx_num);

  if (q_index>=1)
    if (DFTcommon.para_type == spin_type)
        Eigenstate[:,:] = eigenstate_list[cache_index].Eigenvect_real[:,1:energy_idx_num,1,q_index] +
        im * eigenstate_list[cache_index].Eigenvect_imag[:,1:energy_idx_num,1,q_index];
        Eigenvalues[:] = eigenstate_list[cache_index].Eigenvalues[1:energy_idx_num,1,q_index];
    elseif (DFTcommon.colinear_type == spin_type)
        Eigenstate[:,:] = eigenstate_list[cache_index].Eigenvect_real[:,1:energy_idx_num,spin,q_index] +
        im * eigenstate_list[cache_index].Eigenvect_imag[:,1:energy_idx_num,spin,q_index];
        Eigenvalues[:] = eigenstate_list[cache_index].Eigenvalues[1:energy_idx_num,spin,q_index];
    elseif (DFTcommon.non_colinear_type == spin_type)
      Eigenstate[:,:] = eigenstate_list[cache_index].Eigenvect_real[:,1:energy_idx_num,1,q_index] +
      im * eigenstate_list[cache_index].Eigenvect_imag[:,1:energy_idx_num,1,q_index];
      Eigenvalues[:] = eigenstate_list[cache_index].Eigenvalues[1:energy_idx_num,1,q_index];
    end
  end

  return Kpoint_eigenstate_only(Eigenstate,Eigenvalues,k_point);
end

function cacheread_eigenstate_as_nc(k_point::k_point_Tuple,spin,cache_index=1)
  global eigenstate_list;
  global orbital_selection1,orbital_selection2,orbital_selection_on

  #k_point_int = k_point_float2int(kPoint2BrillouinZone_Tuple(k_point));

  # non spin :
  # collinear spin :
  # * spin = 1 down spin
  # * spin = 2 up spin
  # non-collienar spin :
  spin_type = eigenstate_list[cache_index].spin_type
  q_index = cacheread_kpoint2q_index(k_point,cache_index);
  TotalOrbitalNum = eigenstate_list[cache_index].TotalOrbitalNum;
  TotalOrbitalNum2 = TotalOrbitalNum
  TotalOrbitalNum3 = TotalOrbitalNum*2
  if (DFTcommon.non_colinear_type == spin_type)
    TotalOrbitalNum2 = 2*TotalOrbitalNum;
  end
  if (DFTcommon.para_type == spin_type)
      println(" Para Hamiltonian will not be treated as nc.")
      println(" Try to use 'cacheread_eigenstate' instead. ")
      @assert(false)
  end
  Eigenstate = zeros(Complex_my,TotalOrbitalNum3,TotalOrbitalNum3);
  Eigenvalues = zeros(Float_my,TotalOrbitalNum3);

  if (q_index>=1)
    if (DFTcommon.para_type == spin_type)
        Eigenstate[:,:] = eigenstate_list[cache_index].Eigenvect_real[:,:,1,q_index] +
        im * eigenstate_list[cache_index].Eigenvect_imag[:,:,1,q_index];
        Eigenvalues[:] = eigenstate_list[cache_index].Eigenvalues[:,1,q_index];
    elseif (DFTcommon.colinear_type == spin_type)
        Eigenstate[:,:] = eigenstate_list[cache_index].Eigenvect_real[:,:,spin,q_index] +
        im * eigenstate_list[cache_index].Eigenvect_imag[:,:,spin,q_index];
        Eigenvalues[:] = eigenstate_list[cache_index].Eigenvalues[:,spin,q_index];
    elseif (DFTcommon.non_colinear_type == spin_type)
      Eigenstate[:,:] = eigenstate_list[cache_index].Eigenvect_real[:,:,1,q_index] +
      im * eigenstate_list[cache_index].Eigenvect_imag[:,:,1,q_index];
      Eigenvalues[:] = eigenstate_list[cache_index].Eigenvalues[:,1,q_index];
    end
  end

  return Kpoint_eigenstate_only(Eigenstate,Eigenvalues,k_point);
end

function cacheread_Hamiltonian(k_point::k_point_Tuple,Hmode::DFTcommon.nc_Hamiltonian_selection, cache_index=1)
  spin = 1;
  if (DFTcommon.nc_allH == Hmode)
    spin = 1;
  elseif (DFTcommon.nc_realH_only == Hmode)
    spin = 2;
  elseif (DFTcommon.nc_imagH_only == Hmode)
    spin = 3;
  end
  return cacheread_Hamiltonian(k_point,spin,cache_index);
end
function cacheread_Hamiltonian(k_point::k_point_Tuple,spin::Int=1,cache_index=1)
  global eigenstate_list;

  # non spin :
  # collinear spin :
  # * spin = 1 down spin
  # * spin = 2 up spin
  # non-collienar spin :
  spin_type = eigenstate_list[cache_index].spin_type
  q_index = cacheread_kpoint2q_index(k_point,cache_index);
  TotalOrbitalNum = eigenstate_list[cache_index].TotalOrbitalNum;
  TotalOrbitalNum2 = TotalOrbitalNum
  if (DFTcommon.non_colinear_type == spin_type)
    TotalOrbitalNum2 = 2*TotalOrbitalNum;
  end

  Hamiltonian = zeros(Complex_my,TotalOrbitalNum2,TotalOrbitalNum2);

  if (q_index>=1)
    if (DFTcommon.para_type == spin_type)
      Hamiltonian[:,:] = eigenstate_list[cache_index].Hamiltonian_real[:,:,1,q_index] +
      im * eigenstate_list[cache_index].Hamiltonian_imag[:,:,1,q_index];

    elseif (DFTcommon.colinear_type == spin_type)
      Hamiltonian[:,:] = eigenstate_list[cache_index].Hamiltonian_real[:,:,spin,q_index] +
      im * eigenstate_list[cache_index].Hamiltonian_imag[:,:,spin,q_index];

    elseif (DFTcommon.non_colinear_type == spin_type)
      Hamiltonian[:,:] = eigenstate_list[cache_index].Hamiltonian_real[:,:,spin,q_index] +
      im * eigenstate_list[cache_index].Hamiltonian_imag[:,:,spin,q_index];

    end
  end
  return Hamiltonian;
end
function cacheread_Hamiltonian_as_nc(k_point::k_point_Tuple,Hmode::DFTcommon.nc_Hamiltonian_selection, cache_index=1)
    spin = 1;
    spin_type = eigenstate_list[cache_index].spin_type
    if (DFTcommon.non_colinear_type == spin_type)
        if (DFTcommon.nc_allH == Hmode)
          spin = 1;
        elseif (DFTcommon.nc_realH_only == Hmode)
          spin = 2;
        elseif (DFTcommon.nc_imagH_only == Hmode)
          spin = 3;
        end
    end
    return cacheread_Hamiltonian_as_nc(k_point,spin,cache_index);
end

function cacheread_Hamiltonian_as_nc(k_point::k_point_Tuple,spin::Int=1,cache_index=1)
  global eigenstate_list;

  # non spin :
  # collinear spin :
  # * spin = 1 down spin
  # * spin = 2 up spin
  # non-collienar spin :
  spin_type = eigenstate_list[cache_index].spin_type
  q_index = cacheread_kpoint2q_index(k_point,cache_index);
  TotalOrbitalNum = eigenstate_list[cache_index].TotalOrbitalNum;
  TotalOrbitalNum2 = TotalOrbitalNum
  TotalOrbitalNum3 = TotalOrbitalNum*2
  if (DFTcommon.non_colinear_type == spin_type)
    TotalOrbitalNum2 = 2*TotalOrbitalNum;
  end
  if (DFTcommon.para_type == spin_type)
    println(" Para Hamiltonian will not be treated as nc.")
    println(" Try to use 'cacheread_Hamiltonian' instead. ")
    @assert(false)
  end

  Hamiltonian = zeros(Complex_my,TotalOrbitalNum3,TotalOrbitalNum3);

  if (q_index>=1)
    if (DFTcommon.para_type == spin_type)
      Hamiltonian[:,:] = eigenstate_list[cache_index].Hamiltonian_real[:,:,1,q_index] +
      im * eigenstate_list[cache_index].Hamiltonian_imag[:,:,1,q_index];

    elseif (DFTcommon.colinear_type == spin_type)
      Hamiltonian[:,:] = eigenstate_list[cache_index].Hamiltonian_real[:,:,spin,q_index] +
      im * eigenstate_list[cache_index].Hamiltonian_imag[:,:,spin,q_index];

    elseif (DFTcommon.non_colinear_type == spin_type)
      Hamiltonian[:,:] = eigenstate_list[cache_index].Hamiltonian_real[:,:,spin,q_index] +
      im * eigenstate_list[cache_index].Hamiltonian_imag[:,:,spin,q_index];

    end
  end
  return Hamiltonian;
end

function cacheread_atomsOrbital_lists(cache_index=1)
  global eigenstate_list;
  return (eigenstate_list[cache_index].dftresult.basisTransform_result.orbitalStartIdx_list,
  eigenstate_list[cache_index].dftresult.basisTransform_result.orbitalNums)
end
function cacheread_lampup(q_point_list::Array{k_point_Tuple},cache_index=1)
  global eigenstate_list;
  spin_type = eigenstate_list[cache_index].spin_type
  for q_point in q_point_list
    cacheread_eigenstate(q_point,1,cache_index)
  end
  if (DFTcommon.colinear_type == spin_type)
    for q_point in q_point_list
      cacheread_eigenstate(q_point,2,cache_index)
    end
  end
  # Read Hamiltonian
  for q_point in q_point_list
    cacheread_Hamiltonian(q_point,1,cache_index)
  end
  if (DFTcommon.colinear_type == spin_type)
    for q_point in q_point_list
      cacheread_Hamiltonian(q_point,2,cache_index)
    end
  end
  if (DFTcommon.non_colinear_type == spin_type)
    for q_point in q_point_list
      cacheread_Hamiltonian(q_point,1,cache_index)
      cacheread_Hamiltonian(q_point,2,cache_index)
      cacheread_Hamiltonian(q_point,3,cache_index)
    end
  end
end
function cacheread_lampup_as_nc(q_point_list::Array{k_point_Tuple},cache_index=1)
  global eigenstate_list;
  spin_type = eigenstate_list[cache_index].spin_type
  for q_point in q_point_list
    cacheread_eigenstate_as_nc(q_point,1,cache_index)
  end
  #=
  if (DFTcommon.colinear_type == spin_type)
    for q_point in q_point_list
      cacheread_eigenstate(q_point,2,cache_index)
    end
  end
  =#
  # Read Hamiltonian
  for q_point in q_point_list
    cacheread_Hamiltonian_as_nc(q_point,1,cache_index)
  end
  #=
  if (DFTcommon.colinear_type == spin_type)
    for q_point in q_point_list
      cacheread_Hamiltonian(q_point,2,cache_index)
    end
  end
  =#
  if (DFTcommon.non_colinear_type == spin_type)
    for q_point in q_point_list
      cacheread_Hamiltonian_as_nc(q_point,1,cache_index)
      cacheread_Hamiltonian_as_nc(q_point,2,cache_index)
      cacheread_Hamiltonian_as_nc(q_point,3,cache_index)
    end
  end
end

function get_TotalOrbitalNum(result_index=1)
  global dftresult
  TotalOrbitalNum = sum(dftresult[result_index].scf_r.Total_NumOrbs);
  #print(TotalOrbitalNum)
  #print(typeof(TotalOrbitalNum))
  @assert(Int == typeof(TotalOrbitalNum));
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
      @everywhere GC.gc()
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
      @everywhere GC.gc()
    end
  end
  return Q_ksum;
end
function Qspace_Ksum_atomlist_parallel(kq_function,q_point_list,k_point_list,
  atom12_list::Vector{Tuple{Int64,Int64}},num_return=1,
  result_index=1,cache_index=1)
  batch_size = 2*nprocs();
  if (20 > batch_size)
    batch_size = 20;
  end
  cnt = 1
  spin_type = get_dftdataset(result_index).spin_type;

  Xij_Q = Array{Dict{k_point_int_Tuple,Array{Complex_my,1}}}(undef,num_return,length(atom12_list));
  Xij_Q_mean = Array{Dict{k_point_int_Tuple,Complex_my}}(undef,num_return,length(atom12_list));
  for xyz_ij = 1:num_return
      for atom12_i = 1:length(atom12_list)
          Xij_Q[xyz_ij,atom12_i] = Dict{k_point_int_Tuple, Array{Complex_my,1}}();
          Xij_Q_mean[xyz_ij,atom12_i] = Dict{k_point_int_Tuple, Complex_my}();
      end
  end
  #Q_ksum = Dict{k_point_int_Tuple,Array{Complex_my,1}}();
  println(DFTcommon.bar_string) # print ====...====
  p = Progress( ceil(Int, length(q_point_list)/10),
    string("Computing  (Q:",length(q_point_list),", K:",length(k_point_list),")...") );
  p.barglyphs=BarGlyphs("[=> ]")
  p.output = stdout
  for (q_i,q_point) in enumerate(q_point_list)
    q_point_int = k_point_float2int(q_point);

    job_list = Array{Job_input_kq_atom_list_Type}(undef,0)
    for k_point in k_point_list
      kq_point = (q_point[1] + k_point[1],q_point[2] + k_point[2],q_point[3] + k_point[3]) ;
      kq_point = kPoint2BrillouinZone_Tuple(kq_point);
      kq_point_int = k_point_float2int(kq_point);
      push!(job_list,Job_input_kq_atom_list_Type(k_point,kq_point,spin_type,atom12_list));
    end
    X_temp = pmap(kq_function,job_list);
    for xyz_ij = 1:num_return
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
    if (1==rem(q_i,10))
      next!(p)
    end
    if (1==rem(q_i,50))
      @everywhere GC.gc()
    end
  end
  return (Xij_Q,Xij_Q_mean);

end
