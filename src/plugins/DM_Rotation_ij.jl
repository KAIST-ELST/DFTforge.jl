#include("OpenMX_scfout_update.jl")
include("OpenMX_scfout_dm_replace.jl")
using ProgressMeter

import DFTforge
using DFTforge.DFTrefinery
using DFTcommon

function excute_cmd(cmd,work_dir)
  run_error_check = false;
  println(cmd)
  try
    cd(work_dir);
    run(cmd)
  catch ee
    run_error_check = true;
    println(ee)
  end
  if (run_error_check)
    println("=================================================================")
    println(" Error occuured while:")
    println( cmd)
    println("=================================================================")
    @assert(falses)
  end
end


orginal_scfout_fname = "../examples/NiO/16cell_NC_FM/nio.orig.scfout"

orginal_scfout_fname = "nio.orig.scfout"
#orginal_scfout_fname = "../examples/NiO/16cell_NC_rot/nio.orig.scfout"
#arg_input.TOMLinput = "../examples/NiO/16cell_NC_FM/s2p2d1_TMO.toml"
#openmx_fname = "../examples/NiO/16cell_NC_rot/NiO.dat"
openmx_fname = "NiO.dat"


work_dir = dirname(orginal_scfout_fname)
work_dir = pwd()
base_scfout_name = basename(orginal_scfout_fname)
# run time options TODO: set with TOML
mpirun_option = "-np 4"
mpirun = "/home/users1/bluehope/opt/local/mvapich2/2.2/bin/mpirun"
#mpirun_parallel = mpirun*mpirun_option
#mpirun_threadonly = mpirun*" -np 1"

# binary info TODO: set with TOML
#dmft_executable = "/home/users1/bluehope/work_local/DMFT.DFT.scf_read_OpenMX_DMFT/bin/dmft"
openmx_DM_executable = "/home/users1/bluehope/work_local/DMFT.DFT.scf_read_OpenMX_DMFT/bin/openmx3.8.DM"
openmx_POST_executable = "/home/users1/bluehope/work_local/DMFT.DFT.scf_read_OpenMX_DMFT/bin/analysis_example"
rel_dft = "dft"


theta_i = 0.0;
phi_i   = 0.0;
theta_j = 0.25;
phi_j   = 0.0;

atom12 = [1,5]

rotation_dict = Dict{Int,rotation_item}()
for atom_i in 1:8
  rotation_dict[atom_i] = rotation_item(atom_i,-1,1.0,0.0)
end
#rotation_dict[5] = rotation_item(5,-1,0.0,0.0)
#rotation_dict[6] = rotation_item(6,-1,0.0,0.0)

theta_phi_name = string(atom12[1],"_",theta_i,"_",phi_i,"_",atom12[2],"_",theta_j,"_",phi_j);

println("theta_phi_name: ",theta_phi_name)
updated_scfout_fname = joinpath(work_dir,"nio.scfout")
println(updated_scfout_fname)
update_scfout_rotation(orginal_scfout_fname, updated_scfout_fname, rotation_dict)
#scf_r = DFTforge.OpenMXdata.read_scf(orginal_scfout_fname);
#DFTforge.OpenMXdata.write_scf(scf_r, updated_scfout_fname);
#scf_r_test = DFTforge.OpenMXdata.read_scf(updated_scfout_fname);
#=
DM_delta = 0.0;
orbitalStartIdx_list = zeros(Int,scf_r.atomnum)
orbitalNums = copy(scf_r.Total_NumOrbs);
#MPF = Array(Int,scf_r.atomnum)
orbitalStartIdx = 0
for i = 1:scf_r.atomnum
    orbitalStartIdx_list[i] = orbitalStartIdx;
    orbitalStartIdx += orbitalNums[i]
end

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
println("ChemP ",scf_r.ChemP," ",scf_r_test.ChemP)
scf_r.Catomnum
scf_r_test.Catomnum
scf_r_test.TCpyCell
scf_r.rv - scf_r_test.rv
=#

dft_work_dir = work_dir;
println(" dft_work_dir: ",dft_work_dir)

excute_cmd(`$mpirun -np 1 $openmx_DM_executable $openmx_fname -nt 8`,dft_work_dir)
