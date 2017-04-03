module DFTforge
using DFTcommon
export DFTtype, SPINtype
export k_point_Tuple,k_point_int_Tuple,cal_colinear_eigenstate

@enum DFTtype OpenMX = 1 Wannier90 = 2
@enum SPINtype para_type = 1 colinear_type = 2 non_colinear_type = 4

module OpenMXdata
include("backend/OpenMX_PostCommon.jl")
end



type OpenMX_data
    scf_name::AbstractString
    scf_r::OpenMXdata.Openmxscf
    dfttype::DFTforge.DFTtype
end

function read_dftresult(scf_name::AbstractString, dfttype::DFTforge.DFTtype)
    if(OpenMX == dfttype)
      scf_r = OpenMXdata.read_scf(scf_name);
      return scf_r;
    elseif (Wannier90 == dfttype)

    end

end




function cal_colinear_eigenstate(k_point::k_point_Tuple,
    dfttype::DFTforge.DFTtype,scf_r,spin_list=1)
    Eigenstate::Array{Kpoint_eigenstate,1} = Array{Kpoint_eigenstate,1}();
    if(OpenMX==dfttype)
      #Eigenstate::DFTforge.Kpoint_eigenstate =
      Eigenstate =
       OpenMXdata.cal_colinear_eigenstate(k_point,scf_r,spin_list);
    end

    return Eigenstate;
end

function cal_nonco_linear_Eigenstate()
end

################################################################################
module DFTrefinery
using DFTforge
using DFTcommon
using HDF5
using ProgressMeter


export set_current_dftdataset,cal_colinear_eigenstate,get_dftdataset
export Job_input_Type,Job_input_kq_Type
export cachecal_all_Qpoint_eigenstats,cacheset,cacheread_eigenstate
export get_ChempP

  type DFTdataset
    dfttype::DFTtype
    scf_r
    spin_type::SPINtype
  end

  type Eigenstate_hdf5
    hdf_cache_name::AbstractString
    fid_hdf::HDF5.HDF5File
    q_points::Array{k_point_Tuple}
    q_points_int::Array{k_point_int_Tuple};
    q_points_intdic::Dict{k_point_int_Tuple,Int64};
    spin_type::SPINtype
    TotalOrbitalNum::Int64
    dftresult::DFTdataset
    Eigenvect_real::Array{Float64,4}
    Eigenvect_imag::Array{Float64,4}
    Eigenvalues::Array{Float64,3}

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
    Job_input_kq_Type(k_point,kq_point,spin_type) =
      new(k_point,kq_point,spin_type,1)
    Job_input_kq_Type(k_point,kq_point,spin_type,result_index) =
      new(k_point,kq_point,spin_type,result_index)
  end

  global dftresult = Array{DFTdataset}();
  global eigenstate_list =  Array{Eigenstate_hdf5}();

  function set_current_dftdataset(scf_name::AbstractString,
    dfttype::DFTtype,spin_type::SPINtype,result_index=1)
    if(DFTforge.OpenMX == dfttype)
      # Read SCF and Set as current dftdata
      scf_r = DFTforge.OpenMXdata.read_scf(scf_name);
      #assert(0 == scf_r.SpinP_switch - spin_type);
      dftresult[result_index] =
      DFTdataset(dfttype, scf_r,spin_type);
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

  # for pmap
  function cal_eigenstate(input::Job_input_Type)
    # specfify spin type is required

    if(DFTforge.para_type == input.spin_type)
      kpoint_eigenstate_list = Array{Kpoint_eigenstate}();
      push!(kpoint_eigenstate_list,
      cal_colinear_eigenstate(input.k_point,1,input.result_index) );
      return kpoint_eigenstate_list
    elseif(DFTforge.colinear_type ==  input.spin_type)
      return cal_colinear_eigenstate(input.k_point,[1,2],input.result_index)
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

    TotalOrbitalNum2 = TotalOrbitalNum;

    if(DFTforge.non_colinear_type == spin_type)
      TotalOrbitalNum2 = 2*TotalOrbitalNum;
    end

    #print(spin_dim )
    fid_hdf = h5open(hdf_cache_name,"w");
    hdf5_eigenstate_real = d_create(fid_hdf,"Eigenvect_real",datatype(Float64),
    dataspace(TotalOrbitalNum2,TotalOrbitalNum2, spin_dim, Total_q_point_num));

    hdf5_eigenstate_imag = d_create(fid_hdf,"Eigenvect_imag",datatype(Float64),
    dataspace(TotalOrbitalNum2,TotalOrbitalNum2, spin_dim, Total_q_point_num));

    hdf5_eigenvalues = d_create(fid_hdf,"Eigenvalues",datatype(Float64),
    dataspace(TotalOrbitalNum2, spin_dim, Total_q_point_num));

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
    p = Progress(floor(Int64, length(q_point_list)/batch_size),
    "Computing Eigenstates(q)...");

    while cnt < Total_q_point_num
      # pmap
      start_idx = cnt;
      end_idx = minimum([cnt+batch_size-1,Total_q_point_num]);
      temp = pmap(cal_eigenstate,job_list[start_idx:end_idx]);
      ii = 1;
      for jj = start_idx:end_idx

        hdf5_eigenstate_real[:,:,1,jj] = real(temp[ii][1].Eigenstate);
        hdf5_eigenstate_imag[:,:,1,jj] = imag(temp[ii][1].Eigenstate);
        hdf5_eigenvalues[:,1,jj] = temp[ii][1].Eigenvalues;
        if(DFTforge.colinear_type == spin_type )
          hdf5_eigenstate_real[:,:,2,jj] = real(temp[ii][2].Eigenstate);
          hdf5_eigenstate_imag[:,:,2,jj] = imag(temp[ii][2].Eigenstate);
          hdf5_eigenvalues[:,2,jj] = temp[ii][2].Eigenvalues;
        end
        ii += 1;
      end

      # write to hdf5
      cnt = end_idx + 1;
      next!(p)
    end
    close(fid_hdf);
    fid_hdf = h5open(hdf_cache_name,"r");
    q_points_int = Array{k_point_int_Tuple}(Total_q_point_num);


    eigenstate_cache =  Eigenstate_hdf5(hdf_cache_name,fid_hdf,q_point_list,
      q_points_int,q_points_intdic,
      dftresult[result_index].spin_type,TotalOrbitalNum,
      dftresult[result_index]);

    eigenstate_list[cache_index] = eigenstate_cache;
    return eigenstate_cache;
  end

  function cacheset(eigenstate_cache::Eigenstate_hdf5,cache_index=1)
    global eigenstate_list
    fid_hdf = h5open(eigenstate_cache.hdf_cache_name,"r");


    eigenstate_cache.Eigenvect_real = readmmap(fid_hdf["Eigenvect_real"]);
    eigenstate_cache.Eigenvect_imag = readmmap(fid_hdf["Eigenvect_imag"]);
    eigenstate_cache.Eigenvalues     = readmmap(fid_hdf["Eigenvalues"]);
    eigenstate_cache.fid_hdf = fid_hdf;
    eigenstate_list[cache_index] = eigenstate_cache;

  end
  function cacheread_eigenstate(k_point::k_point_Tuple,spin=1,cache_index=1)
    global eigenstate_list;
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
    if(DFTforge.non_colinear_type == spin_type)
      TotalOrbitalNum2 = 2*TotalOrbitalNum;
    end
    Eigenstate = zeros(Complex_my,TotalOrbitalNum2,TotalOrbitalNum2);
    Eigenvalues = zeros(Float_my,TotalOrbitalNum2);

    if(haskey(eigenstate_list[cache_index].q_points_intdic, k_point_int))
      q_index =  eigenstate_list[cache_index].q_points_intdic[k_point_int];
    end
    if(q_index>=1)
      if(DFTforge.para_type == spin_type)
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

  function get_TotalOrbitalNum(result_index=1)
    global dftresult
    TotalOrbitalNum = sum(dftresult[result_index].scf_r.Total_NumOrbs);
    #print(TotalOrbitalNum)
    #print(typeof(TotalOrbitalNum))
    assert(Int64 == typeof(TotalOrbitalNum));
    return TotalOrbitalNum
  end
  function get_ChempP(result_index=1)
    global dftresult
    ChemP::Float64 = dftresult[result_index].scf_r.ChemP;
    return ChemP;
  end

  function cacheread(cache_name::AbstractString)

  end
  function clear_cache(result_index=1)

  end

  function clear_dftdataset()
    dftresult = Array{DFTdataset}();
  end
  end

end
