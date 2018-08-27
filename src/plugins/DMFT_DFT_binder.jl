################################################################################
# DMFT + DFT self consistent calculation
# Bind J.H. Sim  DMFT with OpenMX DFT
################################################################################
using ProgressMeter
import DFTforge
using DFTforge.DFTrefinery
using DFTcommon
using DFTforge.Plugins
################################################################################
# Input Variable & Run time Variable
################################################################################

## 1.2 Read input from argument & TOML file
arg_input = DFTcommon.Arg_Inputs();
arg_input = parse_input(ARGS,arg_input)
#arg_input.TOMLinput = "nio_J_openmx.toml" # Debug
arg_input = parse_TOML(arg_input.TOMLinput,arg_input)
# let argument override
arg_input = parse_input(ARGS,arg_input)

# Input paramters from TOML file
DMFT_DFT_paramter = arg_input.Optional["DMFT_DFT"]::Arg_DMFT_DFT;
start_scfout_fname = arg_input.result_file;

systemName = splitext(basename(arg_input.result_file))[1];
openmx_input           = DMFT_DFT_paramter.openmx_input;
openmx_DM_executable   = DMFT_DFT_paramter.openmx_DM_executable;
openmx_POST_executable = DMFT_DFT_paramter.openmx_POST_executable;

dmft_executable        = DMFT_DFT_paramter.dmft_executable;
max_iter               = DMFT_DFT_paramter.max_iter;
start_iter             = DMFT_DFT_paramter.start_iter

charge_mixing  = DMFT_DFT_paramter.charge_mixing
# Set MPI parmaters (the patched Openmx runs with -np 1)
DMFT_MPI_paramter = arg_input.Optional["MPI"]::Arg_MPI;
DFT_MPI_paramter = deepcopy(DMFT_MPI_paramter);
DFT_MPI_paramter.mpirun_num = 1;
DFT_MPI_paramter.thread_num = DMFT_DFT_paramter.openmx_thread_num;

rel_dft = "dft"
################################################################################
# Input Check
################################################################################
toml_realpath = realpath(arg_input.TOMLinput);
toml_dir =  dirname(toml_realpath);

dmft_work_dir = toml_dir
dft_work_dir = joinpath(toml_dir,"dft")

println(" dmft_work_dir:", dmft_work_dir);
println(" dft_work_dir:", dft_work_dir);

input_error_check  = false;
if !(isfile(start_scfout_fname))
  println("scfout not exist :",start_scfout_fname)
  input_error_check = true;
end
if !(isdir(dft_work_dir))
  println(" make ",dft_work_dir)
  println(" put *.dat with scf.restart  onDM && scf.maxIter  0")
  input_error_check = true;
end
if !(isfile(openmx_input))
  println(" Openmx input dosenot exists :",openmx_input)
  input_error_check = true;
end
if input_error_check
  println(" error occurred")
  @assert(false)
end

function check_openmxInput(openmx_input_fname::String, systemName)
  #TODO: Check below criterias
  # 1. scf.restart                onDM
  # 2. scf.maxIter                0 #5000        # default=40
  # 3. HS.fileout                  on
  # 4. System.Name                 $systemName
  openmx_input_lines =readlines(openmx_input_fname)
  num_lines = length(openmx_input_lines);
  for line_idx = 1:num_lines
    line = openmx_input_lines[line_idx]
    line_splited = split(line);
    if 0 < length(line_splited)
      key_str = line_splited[1];
      if ("scf.restart" == key_str)
        val = line_splited[2];
        @assert("onDM" == val);
      elseif ("scf.maxIter" == key_str)
        val = line_splited[2];
        @assert("0" == val);
      elseif ("HS.fileout" == key_str)
        val = line_splited[2];
        @assert("on" == lowercase(val));
      elseif ("System.Name" == key_str)
        val = line_splited[2];
        @assert(systemName == val);
      end
    end
  end
end
check_openmxInput(openmx_input, systemName)

################################################################################
# Generate internal variables
################################################################################
last_scfout_fname = start_scfout_fname;
base_scfout_name = splitext(basename(start_scfout_fname))[1];
Base.sendfile(start_scfout_fname,joinpath(dft_work_dir, string(base_scfout_name,".orig.scfout")) ); # make backup (H.scfout -> H.orig.scfout)

DM_tot_history = zeros(max_iter+1,1)
run_error_check = false;
################################################################################
# DMFT + DFT Start
################################################################################
cnt = 1;
for iter = start_iter:start_iter+max_iter
  #iter = 1
  println(" DMFT_DFT: iter ",iter)
# 1. read DFT & Gen Hks output for DMFT input
## excute analysis_example
  excute_cmd(openmx_POST_executable, start_scfout_fname, dft_work_dir)                   #run analysis_example
  Base.sendfile(joinpath(dft_work_dir,"Hk.HWR"), joinpath(dmft_work_dir,"Hk.HWR") )      #cp Hk.HWR
  Base.sendfile(joinpath(dft_work_dir,"Hk.HWR"), joinpath(dft_work_dir,string("Hk.",iter,".HWR") )) #history
  Base.sendfile(joinpath(dft_work_dir,"OverlapMatrix.HWR"),joinpath(dmft_work_dir,"OverlapMatrix.HWR"))  #cp OverlapMatrix.HWR
  Base.sendfile(joinpath(dft_work_dir,"OverlapMatrix.HWR"), joinpath(dft_work_dir,string("OverlapMatrix.",iter,".HWR") )) #backup

# 2. Run DMFT
  excute_mpi_cmd(DMFT_MPI_paramter, dmft_executable, ["Hk.HWR",string(iter)], dmft_work_dir)

# 3. Read updated DM from DMFT & Writefor OpenmMX (scfout)
  dm_fname = joinpath(dmft_work_dir,"NumMat_DFT.HWR");

  #TODO: mixing maybe added
  current_scfout_fname = joinpath(dft_work_dir, string(base_scfout_name,".",iter,".scfout")); # ex) H.1.scfout   
  Base.sendfile(last_scfout_fname, current_scfout_fname) # backup with number  ,last_scfout_fname = dft_work_dir/systemName.scfout
  dm_delta = update_scfout(current_scfout_fname, start_scfout_fname, dm_fname, charge_mixing);

# 4. Run Openmx for new Wannier
  excute_mpi_cmd(DFT_MPI_paramter,openmx_DM_executable,openmx_input,dft_work_dir,1)
  Base.sendfile(joinpath(dft_work_dir,string(systemName,".out")), joinpath(dft_work_dir,string(systemName,".",iter,".out") )) #backup
# Last. Check the converence -> exit
# Read input (on the fly stop)
# Gen new SCFOUT for DFT to gen new Hks_R
# goto 1.
  DM_tot_history[cnt] = dm_delta
# Check converge
  println(" dm_delta ",dm_delta)
  println(" DM_tot_history ",DM_tot_history)

  cnt += 1
end


################################################################################
# DMFT + DFT Finished
################################################################################
println(" DM_tot_history ",DM_tot_history)
