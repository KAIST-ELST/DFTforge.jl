using ProgressMeter
import DFTforge
using DFTforge.DFTrefinery
using DFTcommon

export update_scfout,rotation_item;
type  rotation_item
  atom_i::Int
  atom_j::Int
  theta::Float64
  phi::Float64
end

function update_scfout_rotation(orginal_scfout_fname::String, updated_scfout_fname::String,
  rotation_dict::Dict{Int,rotation_item})
  if !(isfile(orginal_scfout_fname))
    println(" Orignal scfout dose not exists ",orginal_scfout_fname)
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
  for GA_AN=1:scf_r.atomnum
    atom1_orbitalNum = scf_r.Total_NumOrbs[GA_AN];
    atom1_orbitalStart = orbitalStartIdx_list[GA_AN];
    for LB_AN = 1:scf_r.FNAN[GA_AN]+1 #atom_i is not atom1,2 index
      GB_AN = scf_r.natn[GA_AN][LB_AN]
      Rn = 1+scf_r.ncn[GA_AN][LB_AN]
      cell_vector  = scf_r.atv_ijk[Rn,:]
      #if (atom_i == GA_AN && atom_j == GB_AN)
      #if (atom_i == GA_AN && atom_j == GB_AN)


      if (haskey(rotation_dict, GA_AN)) # && (0==cell_vector[1] && 0==cell_vector[2] && 0==cell_vector[3]) )
        theta = rotation_dict[GA_AN].theta;
        phi = rotation_dict[GA_AN].phi;
        Rot_theta_phi_mat = RotPauliThetaPhi(pi*theta,pi*phi)


        println(cell_vector)
        atom2_orbitalNum::UInt = scf_r.Total_NumOrbs[GB_AN]
        atom2_orbitalStart::UInt = orbitalStartIdx_list[GB_AN];
        println("GA_AN ",GA_AN," GB_AN " ,GB_AN)
        println(Rot_theta_phi_mat)


        for i = 1:atom1_orbitalNum
          for j = 1:atom2_orbitalNum
            DM_pauli_part = zeros(2,2);


            DM_pauli_part[1,1] =  scf_r.DM[1][GA_AN][LB_AN][i][j];
            DM_pauli_part[2,2] =  scf_r.DM[2][GA_AN][LB_AN][i][j];
            DM_pauli_part[1,2] =  scf_r.DM[3][GA_AN][LB_AN][i][j];
            DM_pauli_part[2,1] =  scf_r.DM[4][GA_AN][LB_AN][i][j];
            n = (DM_pauli_part[1,1] + DM_pauli_part[2,2])/2
            m = DM_pauli_part;
            m[1,1] -= n;
            m[2,2] -= n;

            m2  = (Rot_theta_phi_mat')*m*Rot_theta_phi_mat;
            DM_pauli_part_new = m2;
            DM_pauli_part_new[1,1] += n;
            DM_pauli_part_new[2,2] += n;

            scf_r.DM[1][GA_AN][LB_AN][i][j] = DM_pauli_part_new[1,1];
            scf_r.DM[2][GA_AN][LB_AN][i][j] = DM_pauli_part_new[2,2];
            scf_r.DM[3][GA_AN][LB_AN][i][j] = DM_pauli_part_new[1,2];
            scf_r.DM[4][GA_AN][LB_AN][i][j] = DM_pauli_part_new[2,1];

          end
        end
      end


    end
  end
  ##############################################################################
  # Write new SCFOUT
  # TODO: iDM is not updated yet
  ##############################################################################

  DFTforge.OpenMXdata.write_scf(scf_r, updated_scfout_fname)
  #println(scf_r.ChemP)
  scf_r_test = DFTforge.OpenMXdata.read_scf(orginal_scfout_fname);
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
    #println(" Warnning scfout may not be written corrently [",updated_scfout_fname,"]")
  end
  println("DM_delta ",DM_delta)

  return DM_delta;

  #if (DM_delta >0.001)
  #  println(" Warnning scfout may not be written corrently [",updated_scfout_fname,"]")
  #end

  #return DM_delta_orig;
end
