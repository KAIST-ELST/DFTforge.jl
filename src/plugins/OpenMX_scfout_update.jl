using ProgressMeter
import DFTforge
#using DFTforge.DFTrefinery
using ..DFTcommon

export update_scfout;

function update_scfout(orginal_scfout_fname::String, updated_scfout_fname::String, dm_fname::String, mixing::Float64=0.5)
  if !(isfile(orginal_scfout_fname))
    println(" Orignal scfout dose not exists ",orginal_scfout_fname)
    @assert(falses)
  end
  if !(isfile(dm_fname))
    println(" update DM dose not exists ",dm_fname)
    @assert(falses)
  end
  scf_r = DFTforge.OpenMXdata.read_scf(orginal_scfout_fname);

  ##############################################################################
  # Orbital index generation
  ##############################################################################
  scf_r.SpinP_switch
  spin_type =  DFTcommon.colinear_type
  if (3 == scf_r.SpinP_switch)
    spin_type = DFTcommon.non_colinear_type;
  end

  orbitalStartIdx_list = zeros(Int,scf_r.atomnum)
  orbitalNums = copy(scf_r.Total_NumOrbs);
  #MPF = Array(Int,scf_r.atomnum)
  orbitalStartIdx = 0
  for i = 1:scf_r.atomnum
      orbitalStartIdx_list[i] = orbitalStartIdx;
      orbitalStartIdx += orbitalNums[i]
  end
  ##############################################################################
  #  Dictionary generation
  #  [lmn (cell index) ij (atom index)] to LB_AN
  ##############################################################################
  lmnij2LBAN =  Dict{Tuple{Int32,Int32,Int32,Int32,Int32},Int}()

  for GA_AN=1:scf_r.atomnum
    atom1_orbitalNum = scf_r.Total_NumOrbs[GA_AN];
    atom1_orbitalStart = orbitalStartIdx_list[GA_AN];
    for LB_AN = 1:scf_r.FNAN[GA_AN]+1 #atom_i is not atom1,2 index
        GB_AN = scf_r.natn[GA_AN][LB_AN]
        Rn = 1+scf_r.ncn[GA_AN][LB_AN]
        atom2_orbitalNum::UInt = scf_r.Total_NumOrbs[GB_AN]
        atom2_orbitalStart::UInt = orbitalStartIdx_list[GB_AN];
        #kRn = sum(scf_r.atv_ijk[Rn,:][:].*k_point)
        #println(scf_r.atv_ijk[Rn,:])
        #lmnij_key = vcat(scf_r.atv_ijk[Rn,1:3],GB_AN,GB_AN)
        lmnij_key = (scf_r.atv_ijk[Rn,1],scf_r.atv_ijk[Rn,2],scf_r.atv_ijk[Rn,3],GA_AN,GB_AN)
        lmnij2LBAN[lmnij_key] = LB_AN;

    end
  end
  ##############################################################################
  # Read DM matrix (updated part)
  ##############################################################################
  csv_result = readdlm(dm_fname);
  line_num = size(csv_result)[1]

  DM_delta_orig = 0.0;
  dup_cnt =  Dict{Tuple{Int32,Int32,Int32,Int32,Int32},Float64}()

  ##############################################################################
  # Updated DM matrix
  # TODO: iDM is not updated yet
  ##############################################################################
  #tic()
  for line_index = 1:line_num
    #line_index = 888746;
    lmn_cell_vector = convert(Array{Int32},csv_result[line_index,1:3])
    atom_i = convert(Int32,csv_result[line_index,4])
    atom_j = convert(Int32,csv_result[line_index,5])
    lmnij_key = (lmn_cell_vector[1],lmn_cell_vector[2],lmn_cell_vector[3],atom_i,atom_j)
    if (haskey(lmnij2LBAN,lmnij_key))
      LB_AN = lmnij2LBAN[lmnij_key]

      DM_ij = csv_result[line_index,8]
      iDM_ij= csv_result[line_index,9]
      atom_i_orbital = convert(Int64,csv_result[line_index,6])
      atom_j_orbital = convert(Int64,csv_result[line_index,7])

      spin_i_updown = atom_i_orbital % 2
      spin_j_updown = atom_j_orbital % 2

      atom_i_orbital_rel = convert(Int64,(atom_i_orbital + spin_i_updown) / 2);
      atom_j_orbital_rel = convert(Int64,(atom_j_orbital + spin_j_updown) / 2);


      #@assert(atom_j_orbital_rel <= orbitalStartIdx_list[atom_j]);
      spin = 1;
      if (DFTcommon.colinear_type == spin_type )
        @assert(spin_i_updown == spin_j_updown)
        spin = spin_i_updown + 1
      elseif (DFTcommon.non_colinear_type == spin_type)
        if (0 == spin_i_updown && 0 == spin_j_updown)
          spin = 1;
        elseif (1 == spin_i_updown && 1 == spin_j_updown)
          spin = 2;
        elseif (1 == spin_i_updown && 0 == spin_j_updown)
          spin = 3; # TODO Check required
        elseif (1 == spin_i_updown && 1 == spin_j_updown)
          spin = 4;
        end

      end
      if (atom_j!=scf_r.natn[atom_i][LB_AN])
        println(line_index," a1")
        println(lmnij_key)
        println(spin," ",atom_i," ",LB_AN," ",atom_i_orbital_rel," ",atom_j_orbital_rel," ",atom_j);
        @assert(false)
      end
      if (atom_i_orbital_rel <= orbitalNums[atom_i]) && (atom_j_orbital_rel <= orbitalNums[atom_j])
      else
        println(line_index, " a2")
        println(spin," ",atom_i," ",LB_AN," ",atom_i_orbital_rel," ",atom_j_orbital_rel," ",atom_j);
        @assert(false)
      end
      try

        delta =  abs(scf_r.DM[spin][atom_i][LB_AN][atom_i_orbital_rel][atom_j_orbital_rel] - DM_ij);
        DM_delta_orig += delta
	#simple mixing
	prev_DM = scf_r.DM[spin][atom_i][LB_AN][atom_i_orbital_rel][atom_j_orbital_rel]
        mixed_DM   = 0.5*DM_ij + 0.5*prev_DM
        scf_r.DM[spin][atom_i][LB_AN][atom_i_orbital_rel][atom_j_orbital_rel] =  mixed_DM
	#
        dup_key = (spin,atom_i,LB_AN,atom_i_orbital_rel,atom_j_orbital_rel);
        if (!haskey(dup_cnt,dup_key))
          dup_cnt[dup_key] = delta ;#delta
        else
          #dup_cnt[dup_key] += delta ;#delta
        end

      catch ee
        println(ee)
        println(line_index," a3")
        println(spin," ",atom_i," ",LB_AN," ",atom_i_orbital_rel," ",atom_j_orbital_rel," ",atom_j);
        @assert(false)
      end
      #scf_r.iDM[spin][atom_i][LB_AN][atom_i_orbital_rel][atom_j_orbital_rel] = iDM_ij; # TODO: Waitfor OpenMX tobe patched

      #println( LB_AN)
    end
  end
  #toc()
  sum(0 .== csv_result[:,7])
  println(" DM delta ",DM_delta_orig)

  DM_delta_orig2 = 0.0
  for (k,v) in dup_cnt
    DM_delta_orig2 += v;
  end
  DM_delta_orig2

  ##############################################################################
  # Write new SCFOUT
  # TODO: iDM is not updated yet
  ##############################################################################

  DFTforge.OpenMXdata.write_scf(scf_r, updated_scfout_fname)
  #println(scf_r.ChemP)
  scf_r_test = DFTforge.OpenMXdata.read_scf(updated_scfout_fname);
  scf_r_test.ChemP
  #scf_r_test = DFTforge.OpenMXdata.read_scf("/home/users1/bluehope/work_local/scf_read_OpenMX_DMFT/LaFeAsO/PM/laFeAsO.scfout");

  ##############################################################################
  # check the difference between scf_r & scf_r_test (wirten data)
  # TODO: iDM is not updated yet
  ##############################################################################

  DM_delta = 0.0;
  for GA_AN=1:scf_r_test.atomnum
      atom1_orbitalNum = scf_r_test.Total_NumOrbs[GA_AN];
      atom1_orbitalStart = orbitalStartIdx_list[GA_AN];
      for LB_AN = 1:scf_r_test.FNAN[GA_AN]+1 #atom_i is not atom1,2 index
          GB_AN::UInt = scf_r_test.natn[GA_AN][LB_AN]
          Rn::UInt = 1+scf_r_test.ncn[GA_AN][LB_AN]
          @assert(GB_AN == scf_r.natn[GA_AN][LB_AN])
          @assert(Rn == 1+scf_r.ncn[GA_AN][LB_AN])
          atom2_orbitalNum::UInt = scf_r.Total_NumOrbs[GB_AN]
          atom2_orbitalStart::UInt = orbitalStartIdx_list[GB_AN];

          #
          for spin_i = 1:scf_r_test.SpinP_switch+1
            for i = 1:atom1_orbitalNum
                for j = 1:atom2_orbitalNum
                    DM_old =  scf_r.DM[spin_i][GA_AN][LB_AN][i][j];
                    DM_new =  scf_r_test.DM[spin_i][GA_AN][LB_AN][i][j];
                    DM_delta += abs(DM_new-DM_old);
                end
            end
          end
      end
  end
  if (DM_delta >0.001)
    println(" Warnning scfout may not be written corrently [",updated_scfout_fname,"]")
  end

  return DM_delta_orig;
end
