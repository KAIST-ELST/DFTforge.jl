###############################################################################
# Hongkee Yoon Hongkeeyoon@kaist.ac.kr
# 2019.05
# https://kaist-elst.github.io/DFTforge.jl/
###############################################################################

using ..DFTcommon

#export Wannierdata
export read_wannier

struct Wannierdatatype
  atomnum::Int
  Total_NumOrbs::Array{Int,1}
  SpinP_switch::Int

  tv::Array{Float64,2}
  rv::Array{Float64,2}
  Gxyz::Array{Float64,2}

  Hks_R_raw::Array{Array{Array{Complex_my,2}}}
  Hks_R::Array{Array{Array{Complex_my,2}}}
  R_vector_mat::Array{Array{Int,2}};
  num_degen_list::Array{Int};
  ChemP::Float64

  wannierOrbital2atomGrouped::Array{Int,1}

  E_Temp::Float64;
  spin_type::SPINtype
end

function read_wannier(wannier_fname::AbstractString,result_file_dict::Dict{AbstractString,AbstractString},
  Wannier90_type::Wannier90type, spin_type::SPINtype,
  atoms_orbitals_list::Vector{Array{Int64}},atomnum::Int,atompos::Array{Float64,2})
  @assert(atomnum == length(atoms_orbitals_list))
  if (DFTcommon.OpenMXWF == Wannier90_type)
    if (DFTcommon.para_type == spin_type || DFTcommon.colinear_type == spin_type)
      return read_wannier_OpenMX_ParaCol_linear(wannier_fname, atoms_orbitals_list, atomnum, atompos, spin_type);
    elseif  DFTcommon.non_colinear_type == spin_type
      return read_wannier_OpenMX_nonCol_linear(wannier_fname, atoms_orbitals_list, atomnum, atompos, spin_type);
    else
      throw(assertionError("Unexpected spin type ",spin_type));
    end
  elseif (DFTcommon.EcalJWF  == Wannier90_type)
    return read_wannier_EcalJ(result_file_dict, atoms_orbitals_list, atomnum, atompos, spin_type);
  elseif (DFTcommon.Wannier90WF == Wannier90_type)
    return read_wannier_Wannier90(result_file_dict, atoms_orbitals_list, atomnum, atompos, spin_type);
  end
end

function read_wannier_OpenMX_ParaCol_linear(wannier_fname,atoms_orbitals_list::Vector{Array{Int64}},
  atomnum::Int,atompos::Array{Float64,2}, spin_type::SPINtype)
  @assert(atomnum == length(atoms_orbitals_list))
  #@assert( spin_type ) #TODO: spin_type check
  Total_NumOrbs = Array{Int,1}(undef,0);
  for i in 1:length(atoms_orbitals_list)
    push!(Total_NumOrbs,length(atoms_orbitals_list[i]));
  end
  TotalOrbitalNum = sum(Total_NumOrbs)

  f = open(wannier_fname)
  lines = readlines(f)
  close(f)

  tv = zeros(3,3)
  tv[1,:] = map(x->parse(Float64,x),split(lines[5]))
  tv[2,:] = map(x->parse(Float64,x),split(lines[6]))
  tv[3,:] = map(x->parse(Float64,x),split(lines[7]))
  tv = tv/ang2bohr; # all angstrong

  rv = 2*pi*inv(collect(tv')); # Check required

  ChemP = parse(Float64,split(lines[9])[end])* Hatree2eV
  SpinP_switch  = parse(Int64,split(lines[8])[end])
  # read HWR
  current_line = 10
  start_linenums = Array{Int}(undef,0);
  num_degen_list = Array{Int}(undef,0);
  while (current_line <= length(lines))
    if (startswith(lines[current_line],"R"))
      #println(lines[current_line])
      degen = parse(Int64, split(chomp(lines[current_line]))[end]);
      push!(start_linenums, current_line);
      push!(num_degen_list, degen);
    end
    current_line+= 1
  end
  if (length(start_linenums) > 1)
    linePer_R =  start_linenums[2:end] .- start_linenums[1:end-1] .- 1;
    @assert(0 == sum(false .==(linePer_R[1] .== linePer_R)))
    @assert(Float64(Int64(sqrt(linePer_R[1]))) == sqrt(linePer_R[1]))
    @assert(TotalOrbitalNum == Int64(sqrt(linePer_R[1])));
  end

  @assert(sum(Total_NumOrbs)==TotalOrbitalNum);
  wannierOrbital2atomGrouped = zeros(Int,TotalOrbitalNum)
  cnt = 1;
  for i in 1:length(atoms_orbitals_list)
    for j in 1:length(atoms_orbitals_list[i])
      wannierOrbital2atomGrouped[cnt] = atoms_orbitals_list[i][j];
      cnt+=1;
    end
  end

  R_vector_mat = Array{Array{Int,2}}(undef,SpinP_switch)
  #R_vector_mat = zeros(length(start_linenums),3)
  Hks_R = Array{Array{Array{Complex_my,2}}}(undef,SpinP_switch)
  for spin in 1:SpinP_switch
    Hks_R[spin] = Array{Array{Complex_my,2}}(undef,0);
    R_vector_mat[spin] = zeros(convert(Int,length(start_linenums)/SpinP_switch),3);
  end
  spin = 1;

  cnt = 1;
  for (i,start_linenum) in enumerate(start_linenums)
    if ((length(start_linenums)/2) +1 == i )
      spin = 2;
      cnt = 1;
    end
    R_vector =  map(x->parse(Int64,x),split(lines[start_linenum])[3:5])
    R_vector_mat[spin][cnt,:] = R_vector

    HWR_mat = zeros(Complex_my, TotalOrbitalNum,TotalOrbitalNum);
    rel_line = 0;
    for x in 1:TotalOrbitalNum
      for y in 1:TotalOrbitalNum
        rel_line += 1
            #println(split(chomp(lines[start_linenum+rel_line]) ))
    		tmp = map(x->parse(Float64,x),split(chomp(lines[start_linenum+rel_line]))[3:4]);
    		#tmp = map(x->parse(Float64,x),split(chomp(lines[start_linenums+rel_line]))[3:4])
    		value = tmp[1] + tmp[2]*im;
    		HWR_mat[x,y] = value;
	    end
    end
    HWR_mat *= Hatree2eV/Float64(num_degen_list[i]);
    push!(Hks_R[spin],HWR_mat);
    cnt += 1;
  end
  # re arange wannier group same atoms orbitals
  Hks_R_grouped = Array{Array{Array{Complex_my,2}}}(undef,SpinP_switch)
  for spin in 1:SpinP_switch
    Hks_R_grouped[spin] = Array{Array{Complex_my,2}}(undef,0);
  end
  for spin in 1:SpinP_switch
    for i in 1:length(Hks_R[spin])

      HWR_mat = zeros(Complex_my, TotalOrbitalNum,TotalOrbitalNum);
      HWR_mat[:,:] = Hks_R[spin][i][wannierOrbital2atomGrouped,wannierOrbital2atomGrouped]
      push!(Hks_R_grouped[spin],HWR_mat);
    end
  end
  Gxyz = zeros(atomnum,3);
  for i in 1:atomnum
    Gxyz[i,:] += atompos[i,1]*tv[1,:] ;
    Gxyz[i,:] += atompos[i,2]*tv[2,:] ;
    Gxyz[i,:] += atompos[i,3]*tv[3,:] ;
  end

  Wannier_info = Wannierdatatype(atomnum, Total_NumOrbs, SpinP_switch, tv, rv,
    Gxyz, Hks_R, Hks_R_grouped, R_vector_mat, num_degen_list, ChemP,
    wannierOrbital2atomGrouped, 300.0, spin_type)
  return Wannier_info;
end

function read_wannier_OpenMX_nonCol_linear(wannier_fname,atoms_orbitals_list::Vector{Array{Int64}},
  atomnum::Int,atompos::Array{Float64,2}, spin_type::DFTcommon.SPINtype)
  @assert(atomnum == length(atoms_orbitals_list))
  #@assert( spin_type ) #TODO: spin_type check
  Total_NumOrbs = Array{Int,1}(undef,0);
  for i in 1:length(atoms_orbitals_list)
    push!(Total_NumOrbs,length(atoms_orbitals_list[i]));
  end
  TotalOrbitalNum = sum(Total_NumOrbs)
  TotalOrbitalNum2 = TotalOrbitalNum * 2;
  f = open(wannier_fname)
  lines = readlines(f)
  close(f)

  tv = zeros(3,3)
  tv[1,:] = map(x->parse(Float64,x),split(lines[5]))
  tv[2,:] = map(x->parse(Float64,x),split(lines[6]))
  tv[3,:] = map(x->parse(Float64,x),split(lines[7]))

  tv = tv/ang2bohr; # all angstrong
  rv = 2*pi*inv(tv'); # Check required

  ChemP = parse(Float64,split(lines[9])[end])*DFTcommon.Hatree2eV
  SpinP_switch  = parse(Int64,split(lines[8])[end])
  @assert(1 == SpinP_switch); #
  # read HWR
  current_line = 10
  start_linenums = Array{Int}(undef,0);
  num_degen_list = Array{Int}(undef,0);
  while (current_line <= length(lines))
    if (startswith(lines[current_line],"R"))
      #println(lines[current_line])
      degen = parse(Int64, split(chomp(lines[current_line]))[end]);
      push!(start_linenums, current_line);
      push!(num_degen_list, degen);
    end
    current_line+= 1
  end
  if (length(start_linenums) > 1)
    linePer_R =  start_linenums[2:end] - start_linenums[1:end-1] - 1;
    @assert(0 == sum(false .==(linePer_R[1] .== linePer_R)))
    @assert(Float64(Int64(sqrt(linePer_R[1]))) == sqrt(linePer_R[1]))
    @assert(TotalOrbitalNum2 == Int64(sqrt(linePer_R[1])))
  end

  @assert(2 * TotalOrbitalNum == TotalOrbitalNum2); # non-collinear spin
  wannierOrbital2atomGrouped = zeros(Int,TotalOrbitalNum2)
  cnt = 1;
  for i in 1:length(atoms_orbitals_list)
    for j in 1:length(atoms_orbitals_list[i])
      wannierOrbital2atomGrouped[cnt] = atoms_orbitals_list[i][j]; # non-collinear spin up
      wannierOrbital2atomGrouped[TotalOrbitalNum+cnt] = TotalOrbitalNum + atoms_orbitals_list[i][j]; # non-collinear spin down
      cnt+=1;
    end
  end
  #println(wannierOrbital2atomGrouped)

  R_vector_mat = Array{Array{Int,2}}(SpinP_switch)
  #R_vector_mat = zeros(length(start_linenums),3)
  Hks_R = Array{Array{Array{Complex_my,2}}}(SpinP_switch)
  for spin in 1:SpinP_switch
    Hks_R[spin] = Array{Array{Complex_my,2}}(undef,0);
    R_vector_mat[spin] = zeros(convert(Int,length(start_linenums)/SpinP_switch),3);
  end
  spin = 1;

  cnt = 1;
  for (i,start_linenum) in enumerate(start_linenums)
    if ((length(start_linenums)/2) +1 == i )
      spin = 2;
      cnt = 1;
    end
    R_vector =  map(x->parse(Int64,x),split(lines[start_linenum])[3:5])
    R_vector_mat[spin][cnt,:] = R_vector

    HWR_mat = zeros(Complex_my, TotalOrbitalNum2,TotalOrbitalNum2);
    rel_line = 0;
    for x in 1:TotalOrbitalNum2
      for y in 1:TotalOrbitalNum2
        rel_line += 1

        #println(split(chomp(lines[start_linenum+rel_line]) ))
    		tmp = map(x->parse(Float64,x),split(chomp(lines[start_linenum+rel_line]))[3:4]);
        tmp_ab = map(x->parse(Int64,x),split(chomp(lines[start_linenum+rel_line]))[1:2]);

        @assert(tmp_ab[1] == x && tmp_ab[2] == y )
    		value = tmp[1] + tmp[2]*im;
    		HWR_mat[x,y] = value;
	    end
    end
    HWR_mat *= DFTcommon.Hatree2eV/Float64(num_degen_list[i]);
    push!(Hks_R[spin], HWR_mat);
    cnt += 1;
  end
  # re arange wannier group same atoms orbitals
  Hks_R_grouped = Array{Array{Array{Complex_my,2}}}(SpinP_switch)
  for spin in 1:SpinP_switch
    Hks_R_grouped[spin] = Array{Array{Complex_my,2}}(undef,0);
  end
  for spin in 1:SpinP_switch
    for i in 1:length(Hks_R[spin])

      HWR_mat = zeros(Complex_my, TotalOrbitalNum2, TotalOrbitalNum2);
      HWR_mat[:,:] = Hks_R[spin][i][wannierOrbital2atomGrouped, wannierOrbital2atomGrouped]
      push!(Hks_R_grouped[spin], HWR_mat);
    end
  end
  Gxyz = zeros(atomnum,3);
  for i in 1:atomnum
    Gxyz[i,:] += atompos[i,1]*tv[1,:] ;
    Gxyz[i,:] += atompos[i,2]*tv[2,:] ;
    Gxyz[i,:] += atompos[i,3]*tv[3,:] ;
  end

  Wannier_info = Wannierdatatype(atomnum, Total_NumOrbs, SpinP_switch, tv, rv,
    Gxyz, Hks_R, Hks_R_grouped, R_vector_mat, num_degen_list, ChemP,
    wannierOrbital2atomGrouped, 300.0, spin_type)
  return Wannier_info;
end


function read_wannier_Wannier90(result_file_dict::Dict{AbstractString,AbstractString},
  atoms_orbitals_list::Vector{Array{Int64}},
  atomnum::Int,atompos::Array{Float64,2}, spin_type::DFTcommon.SPINtype)
  if (DFTcommon.colinear_type == spin_type)
      Wannier_info_up = read_wannier_Wannier90_internal(result_file_dict["result_file_up"],
        atoms_orbitals_list,
        atomnum, atompos, spin_type)
      Wannier_info_down = read_wannier_Wannier90_internal(result_file_dict["result_file_down"],
        atoms_orbitals_list,
        atomnum, atompos, spin_type)

      Total_NumOrbs = Wannier_info_up.Total_NumOrbs;
      @assert(Wannier_info_up.Total_NumOrbs == Wannier_info_down.Total_NumOrbs);
      @assert(Wannier_info_up.num_degen_list == Wannier_info_down.num_degen_list)
      tv = Wannier_info_up.tv;
      @assert(Wannier_info_up.tv ==  Wannier_info_down.tv);
      rv = Wannier_info_up.rv;
      @assert(Wannier_info_up.rv == Wannier_info_down.rv);
      Gxyz = Wannier_info_up.Gxyz;
      @assert(Wannier_info_up.Gxyz == Wannier_info_down.Gxyz)
      Hks_R_raw = Array{Array{Array{Complex_my,2}}}(undef,2)
      Hks_R_raw[1] = Wannier_info_up.Hks_R_raw[1];
      Hks_R_raw[2] = Wannier_info_down.Hks_R_raw[1];

      Hks_R = Array{Array{Array{Complex_my,2}}}(undef,2)
      Hks_R[1] = Wannier_info_up.Hks_R[1];
      Hks_R[2] = Wannier_info_down.Hks_R[1];


      #R_vector_mat = Wannier_info_up.R_vector_mat;
      R_vector_mat = Array{Array{Int,2}}(undef,2)
      @assert(Wannier_info_up.R_vector_mat == Wannier_info_down.R_vector_mat);
      R_vector_mat[1]  = Wannier_info_up.R_vector_mat[1]
      R_vector_mat[2]  = Wannier_info_down.R_vector_mat[1]

      ChemP = Wannier_info_up.ChemP;
      #@assert(0.0 == ChemP); # In EcalJ Wannier
      @assert(Wannier_info_up.ChemP == Wannier_info_down.ChemP)

      wannierOrbital2atomGrouped = Wannier_info_up.wannierOrbital2atomGrouped;
      @assert(Wannier_info_up.wannierOrbital2atomGrouped == Wannier_info_down.wannierOrbital2atomGrouped)

      Wannier_info = Wannierdatatype(atomnum, Total_NumOrbs, 2, tv, rv,
          Gxyz, Hks_R_raw, Hks_R, R_vector_mat, Wannier_info_up.num_degen_list, ChemP,
          wannierOrbital2atomGrouped, 300.0, spin_type)
      return Wannier_info;
  elseif (DFTcommon.para_type == spin_type)
      return  read_wannier_Wannier90_internal(result_file_dict["result_file_up"],
          atoms_orbitals_list,
          atomnum, atompos, spin_type)
   end

end

function read_wannier_Wannier90_internal(wannier_fname::AbstractString,
  atoms_orbitals_list::Vector{Array{Int64}},
  atomnum::Int, atompos::Array{Float64,2}, spin_type::DFTcommon.SPINtype)

  @assert(atomnum == length(atoms_orbitals_list))
  @assert(atomnum == size(atompos)[1]);

  #TODO: Read ChemP froms somewhere



  Total_NumOrbs = Array{Int,1}(undef,0);
  for i in 1:length(atoms_orbitals_list)
    push!(Total_NumOrbs,length(atoms_orbitals_list[i]));
  end
  # read wannier90 "win" file

  # Read Cell Vector
  f = open((wannier_fname * ".win"))
  println(wannier_fname * ".win")
  wannier90_win_file = readlines(f)
  close(f)
  wannier90_win_file_lowercase = map(x->strip(lowercase(x)), wannier90_win_file);
  # Read cell vector
  cell_vect_start_line = findfirst("begin unit_cell_cart" .== wannier90_win_file_lowercase ) + 1
  cell_vect_end_line = findfirst("end unit_cell_cart" .== wannier90_win_file_lowercase ) -1
  #println(cell_vect_start_line,cell_vect_end_line)
  @assert(cell_vect_end_line > cell_vect_start_line)
  @assert(3 <= (cell_vect_end_line - cell_vect_start_line + 1)  )

  tv = zeros(3,3)
  tv[1,:] = map(x->parse(Float64,x),split(wannier90_win_file_lowercase[cell_vect_end_line-2]))
  tv[2,:] = map(x->parse(Float64,x),split(wannier90_win_file_lowercase[cell_vect_end_line-1]))
  tv[3,:] = map(x->parse(Float64,x),split(wannier90_win_file_lowercase[cell_vect_end_line-0]))

  rv = 2*pi*inv(collect(tv')); # Check required

  #ChemP = 8.1400
  # Read Chemical potentail
  chemp_line = 0
  try
    chemp_line = findfirst( 0 .< map(x-> sum("fermi_energy".==split(x,['=',' '],keepempty=false)), wannier90_win_file_lowercase));
    if (nothing==chemp_line)
        println(" No fermi_energy found in " * wannier_fname * ".win ex) fermi_energy = 1.234 !eV ")
        @assert(true)
    end
  catch
    println(" No fermi_energy found in " * wannier_fname * ".win ex) fermi_energy = 1.234 !eV ")
  end
  ChemP = parse(Float64, split(wannier90_win_file_lowercase[chemp_line],['=',' '],keepempty=false)[2])
  println(wannier90_win_file_lowercase[chemp_line])

  num_wann_line  = findfirst( 0 .< map(x-> sum("num_wann".==split(x,['=',' '],keepempty=false)), wannier90_win_file_lowercase));
  num_wann = parse(Int64, split(wannier90_win_file_lowercase[num_wann_line],['=',' '],keepempty=false)[2]);

  # Read position vector
  #=
  pos_vect_start_line = findfirst(lowercase("Begin Projections") .== wannier90_win_file_lowercase) + 1
  pos_vect_end_line = findfirst(lowercase("End Projections") .== wannier90_win_file_lowercase) - 1

  postion_list = [];
  for line_num = pos_vect_start_line:pos_vect_end_line
    postion_raw_str = wannier90_win_file[line_num];
    postion_raw_str = strip(postion_raw_str);
    if !startswith(postion_raw_str,"#")
        postion_1 = map(x->parse(Float64,x),
        split(split(postion_raw_str,":")[1],[',','='])[2:4])
        push!(postion_list, postion_1)
    end
  end
  # Check if atom infomations are correct
  check_fail = false;
  for (i,v) in enumerate(atoms_orbitals_list)
      if (sum(abs.(atompos[i,:] - postion_list[i])) > 10.0^-4.0)
          println(" read_wannier_Wannier90 atom ",i," position not matched ",position_1,"  ",atompos[i,:])
          check_fail = true;
      end
  end
  if check_fail
    println(" read_wannier_Wannier90 atom position not matched ")
    println(DFTcommon.bar_string) # print ====...====
    @assert(true);
  end
  =#

  # Read wannier H
  f = open((wannier_fname * "_hr.dat"))
  wannier90_hr_file = readlines(f)
  close(f)

  TotalOrbitalNum = parse(Int64,split(wannier90_hr_file[2])[1])
  @assert(num_wann == TotalOrbitalNum); # Check orbital numbers between .win and _hr.dat file
  num_Rvector = parse(Int64,split(wannier90_hr_file[3])[1])

  # Read num_degen_list
  num_degen_linenum = Int64(ceil(num_Rvector/15))
  num_degen_list = Array{Int}(undef,0);
  for num_degen_lineidx = 1:num_degen_linenum
    tmp = map( x->parse(Int64,x),split(wannier90_hr_file[4+num_degen_lineidx-1 ]));
    append!(num_degen_list,tmp)
  end
  num_degen_list

  # wannierOrbital2atomGrouped
  wannierOrbital2atomGrouped = zeros(Int,TotalOrbitalNum)
  cnt = 1;
  for i in 1:length(atoms_orbitals_list)
    for j in 1:length(atoms_orbitals_list[i])
      wannierOrbital2atomGrouped[cnt] = atoms_orbitals_list[i][j];
      cnt+=1;
    end
  end

  # Read wannier H

  start_linenum = 4 + num_degen_linenum ;
  end_linenum = length(wannier90_hr_file)

  @assert((end_linenum - start_linenum+1) == TotalOrbitalNum*TotalOrbitalNum*num_Rvector);

  SpinP_switch = 1
  spin = 1
  R_vector_mat = Array{Array{Int,2}}(undef,SpinP_switch)
  Hks_R = Array{Array{Array{Complex_my,2}}}(undef,SpinP_switch)
  for spin in 1:SpinP_switch
    Hks_R[spin] = Array{Array{ComplexF64,2}}(undef,0);
    R_vector_mat[spin] = zeros(num_Rvector,3);
  #end
    for R_vect_idx = 1:num_Rvector
      start_cell_vect =  start_linenum + (R_vect_idx-1)*TotalOrbitalNum^2;
      HWR_mat = zeros(Complex_my, TotalOrbitalNum,TotalOrbitalNum);
      rel_line = 0;

      R_vector =  map(x->parse(Int,x),split(chomp(wannier90_hr_file[start_cell_vect+0]))[1:3]);
      R_vector_mat[spin][R_vect_idx,:] = R_vector;
      for y in 1:TotalOrbitalNum
          for x in 1:TotalOrbitalNum
              tmp = map(x->parse(Float64,x),split(chomp(wannier90_hr_file[start_cell_vect+rel_line]))[6:7]);
              value = tmp[1] + tmp[2]*im;
      		HWR_mat[x,y] = value;

              tmp = map(x->parse(Int,x),split(chomp(wannier90_hr_file[start_cell_vect+rel_line]))[1:5]);
              @assert(tmp[1:3] == R_vector) # Check if same R_vector
              #println(tmp)
              @assert(tmp[4] == x && tmp[5] == y)

              rel_line += 1;
          end
      end
      HWR_mat *= 1.0 # Transform to eV unit
      push!(Hks_R[spin],HWR_mat);
    end
  end

  # re arange wannier group same atoms orbitals
  Hks_R_grouped = Array{Array{Array{Complex_my,2}}}(undef,SpinP_switch)
  for spin in 1:SpinP_switch
    Hks_R_grouped[spin] = Array{Array{Complex_my,2}}(undef,0);
  end
  for spin in 1:SpinP_switch
    for i in 1:length(Hks_R[spin])

      HWR_mat = zeros(Complex_my, TotalOrbitalNum,TotalOrbitalNum);
      HWR_mat[:,:] = Hks_R[spin][i][wannierOrbital2atomGrouped,wannierOrbital2atomGrouped]
      push!(Hks_R_grouped[spin],HWR_mat);
    end
  end
  Gxyz = zeros(atomnum,3);
  for i in 1:atomnum
    Gxyz[i,:] += atompos[i,1]*tv[1,:] ;
    Gxyz[i,:] += atompos[i,2]*tv[2,:] ;
    Gxyz[i,:] += atompos[i,3]*tv[3,:] ;
  end



  Wannier_info = Wannierdatatype(atomnum, Total_NumOrbs, SpinP_switch, tv, rv,
    Gxyz, Hks_R, Hks_R_grouped, R_vector_mat, num_degen_list, ChemP,
    wannierOrbital2atomGrouped, 300.0, spin_type)
  return Wannier_info;
end

function read_wannier_EcalJ(result_file_dict::Dict{AbstractString,AbstractString},
  atoms_orbitals_list::Vector{Array{Int64}},
  atomnum::Int,atompos::Array{Float64,2}, spin_type::DFTcommon.SPINtype)

  if (DFTcommon.colinear_type == spin_type)
    @assert(2 == length(result_file_dict) );
    Wannier_info_up = read_wannier_EcalJInternal(result_file_dict["result_file_up"],
      atoms_orbitals_list,
      atomnum, atompos, spin_type)

    Wannier_info_down = read_wannier_EcalJInternal(result_file_dict["result_file_down"],
      atoms_orbitals_list,
      atomnum, atompos, spin_type)

    Total_NumOrbs = Wannier_info_up.Total_NumOrbs;
    @assert(Wannier_info_up.Total_NumOrbs == Wannier_info_down.Total_NumOrbs);
    @assert(Wannier_info_up.num_degen_list == Wannier_info_down.num_degen_list)
    tv = Wannier_info_up.tv;
    @assert(Wannier_info_up.tv ==  Wannier_info_down.tv);
    rv = Wannier_info_up.rv;
    @assert(Wannier_info_up.rv == Wannier_info_down.rv);
    Gxyz = Wannier_info_up.Gxyz;
    @assert(Wannier_info_up.Gxyz == Wannier_info_down.Gxyz)
    Hks_R_raw = Array{Array{Array{Complex_my,2}}}(undef,2)
    Hks_R_raw[1] = Wannier_info_up.Hks_R_raw[1];
    Hks_R_raw[2] = Wannier_info_down.Hks_R_raw[1];

    Hks_R = Array{Array{Array{Complex_my,2}}}(undef,2)
    Hks_R[1] = Wannier_info_up.Hks_R[1];
    Hks_R[2] = Wannier_info_down.Hks_R[1];


    #R_vector_mat = Wannier_info_up.R_vector_mat;
    R_vector_mat = Array{Array{Int,2}}(undef,2)
    @assert(Wannier_info_up.R_vector_mat == Wannier_info_down.R_vector_mat);
    R_vector_mat[1]  = Wannier_info_up.R_vector_mat[1]
    R_vector_mat[2]  = Wannier_info_down.R_vector_mat[1]

    ChemP = Wannier_info_up.ChemP;
    #@assert(0.0 == ChemP); # In EcalJ Wannier
    @assert(Wannier_info_up.ChemP == Wannier_info_down.ChemP)

    wannierOrbital2atomGrouped = Wannier_info_up.wannierOrbital2atomGrouped;
    @assert(Wannier_info_up.wannierOrbital2atomGrouped == Wannier_info_down.wannierOrbital2atomGrouped)

    Wannier_info = Wannierdatatype(atomnum, Total_NumOrbs, 2, tv, rv,
        Gxyz, Hks_R_raw, Hks_R, R_vector_mat, Wannier_info_up.num_degen_list, ChemP,
        wannierOrbital2atomGrouped, 300.0, spin_type)
    return Wannier_info;
  elseif (DFTcommon.para_type == spin_type)
    return  read_wannier_EcalJInternal(result_file_dict["result_file_up"],
        atoms_orbitals_list,
        atomnum, atompos, spin_type)
  end
end

function read_wannier_EcalJInternal(wannier_fname::AbstractString,
  atoms_orbitals_list::Vector{Array{Int64}},
  atomnum::Int, atompos::Array{Float64,2}, spin_type::DFTcommon.SPINtype)

  @assert(atomnum == length(atoms_orbitals_list))
  @assert(atomnum == size(atompos)[1]);


  Total_NumOrbs = Array{Int,1}(undef,0);
  for i in 1:length(atoms_orbitals_list)
    push!(Total_NumOrbs,length(atoms_orbitals_list[i]));
  end

  f = open(wannier_fname)
  lines = readlines(f)
  close(f)

  fractional_scale = 1.0 #1.224860;
  fractional_scale = parse(Float64, lines[2])

  TotalOrbitalNum = parse(Int64,split(lines[6])[1]);
  println("sum(Total_NumOrbs):",sum(Total_NumOrbs)," TotalOrbitalNum:", TotalOrbitalNum)
  @assert(sum(Total_NumOrbs)==TotalOrbitalNum);
  num_Rvector = parse(Int64,split(lines[6])[2]);
  total_lines = parse(Int64,split(lines[6])[3]);
  @assert(total_lines == TotalOrbitalNum^2 * num_Rvector);

  tv = zeros(3,3)
  tv[1,:] = map(x->parse(Float64,x),split(lines[3]))
  tv[2,:] = map(x->parse(Float64,x),split(lines[4]))
  tv[3,:] = map(x->parse(Float64,x),split(lines[5]))

  tv *= fractional_scale;

  rv = 2*pi*inv(tv'); # Check required

  ChemP = parse(Float64,split(lines[7])[end]) #* DFTcommon.Hatree2eV/2 # TODO:Check for Wannier Chem potential parse(Float64,split(lines[9])[end])*DFTcommon.Hatree2eV
  SpinP_switch  = 1 #parse(Int64,split(lines[6])[end])

  # Check if atom infomations are correct
  check_fail = false;
  for (i,v) in enumerate(atoms_orbitals_list)
    position_1 = map(x->parse(Float64,x), split(lines[9+v[1]])[end-2:end]);
    #position_1 = map(x-> rem(x,1.0), position_1 + 2.0)
    if (sum(abs.(atompos[i,:] - position_1)) > 10.0^-4.0)
      println(" read_wannier_EcalJ orbital index ",i," position not matched Wannier file:",position_1," input Toml file:",atompos[i,:])
      check_fail = true;
    end
    for (i2,v2) in enumerate(v)
      position_2 = map(x->parse(Float64,x), split(lines[9+v[i2] ])[end-2:end]);
      #position_2 = map(x-> rem(x,1.0), position_2 + 2.0)
      if (sum(abs.(atompos[i,:] - position_2)) > 10.0^-4.0)
        println(" read_wannier_EcalJ atom ",i2," position not matched Wannier file:",position_2," input Toml file:",atompos[i,:])
        check_fail = true;
      end
    end
  end
  if check_fail
    println(" read_wannier_EcalJ atom position not matched ")
    println(DFTcommon.bar_string) # print ====...====
    @assert(true);
  end
  # Read wannier H
  R_vector_mat = Array{Array{Int,2}}(undef,SpinP_switch)

  Hks_R = Array{Array{Array{Complex_my,2}}}(undef,SpinP_switch)
  wannierOrbital2atomGrouped = zeros(Int,TotalOrbitalNum)

  ## Check R vectors
  # TODO: Spin up down
  start_linenums = Array{Int}(undef,0);
  num_degen_list = Array{Int}(undef,0);
  start_line = 8 + TotalOrbitalNum + 2;
  R_vector_catesian = map(x->parse(Float64,x),split(lines[start_line])[4:6]) * fractional_scale;
  degen = parse(Int64, split(chomp(lines[start_line]))[7]);
  R_frac =  rv * R_vector_catesian[:];
  R_int = round.(Int64,R_frac)
  # store current R_vect
  prev_R_vect = R_vector_catesian;
  push!(start_linenums, start_line);
  push!(num_degen_list, degen )

  @assert(sum(abs.(R_frac - R_int)) < 10.0^-4.0); # Check if cell vector is Int
  for current_line in start_line:length(lines)
    R_vector_catesian = map(x->parse(Float64,x),split(lines[current_line])[4:6]) * fractional_scale;
    degen = parse(Int64, split(chomp(lines[current_line]))[7]);

    if sum(abs.(prev_R_vect-R_vector_catesian)) < 10.0^-4.0
    else
      # new R vector
      push!(start_linenums,current_line);
      push!(num_degen_list, degen )

      R_frac =  (rv) * R_vector_catesian[:]/(2*pi);
      R_int = round.(Int64,R_frac)
      if (sum(abs.(R_frac - R_int)) > 10.0^-4.0)
        println(" Cartesian R vector  ",R_vector_catesian," did not convert to fractional R", R_frac," @line# ",current_line)
        @assert(false); # Check if cell vector is Int
      end
      prev_R_vect = R_vector_catesian;
    end
  end


  cnt = 1;
  for i in 1:length(atoms_orbitals_list)
    for j in 1:length(atoms_orbitals_list[i])
      wannierOrbital2atomGrouped[cnt] = atoms_orbitals_list[i][j];
      cnt+=1;
    end
  end
  #for spin in 1:SpinP_switch
    Hks_R[1] = Array{Array{Complex_my,2}}(undef,0);
    R_vector_mat[1] = zeros(length(start_linenums),3);

  #end
  spin = 1;

  cnt = 1;
  @assert(length(start_linenums)== length(num_degen_list))
  for (i,start_linenum) in enumerate(start_linenums)

    R_vector_catesian =  map(x->parse(Float64,x),split(lines[start_linenum])[4:6]) * fractional_scale
    R_frac =  (rv) * R_vector_catesian[:]/(2*pi);
    R_int = round.(Int64,R_frac)
    if (sum(abs.(R_frac - R_int)) > 10.0^-4.0)
      println(" Cartesian R vector  ",R_vector_catesian," did not convert to fractional R", R_frac," @line# ",start_linenum)
      @assert(false); # Check if cell vector is Int
    end
    R_vector_mat[spin][cnt,:] = R_int;

    HWR_mat = zeros(Complex_my, TotalOrbitalNum,TotalOrbitalNum);
    rel_line = 0;
    for x in 1:TotalOrbitalNum
      for y in 1:TotalOrbitalNum

        #println(split(chomp(lines[start_linenum+rel_line]) ))
        tmp = map(x->parse(Int64,x),split(chomp(lines[start_linenum+rel_line]))[8:9]);
        if ( (tmp[1] != x) || (tmp[2] != y) )
          println(" Orbital index l1,l2 (",x,",",y,") not matched (",tmp[1],",",tmp[2],")")
          @assert(false);
        end

    		tmp = map(x->parse(Float64,x),split(chomp(lines[start_linenum+rel_line]))[10:11]);
        rel_line += 1
    		#tmp = map(x->parse(Float64,x),split(chomp(lines[start_linenums+rel_line]))[3:4])
    		value = tmp[1] + tmp[2]*im;
    		HWR_mat[x,y] = value;
	    end
    end
    #HWR_mat *= DFTcommon.Hatree2eV;
    HWR_mat *=1.0/Float64(num_degen_list[i]);
    push!(Hks_R[spin],HWR_mat);
    cnt += 1;
  end

  # re arange wannier group same atoms orbitals
  Hks_R_grouped = Array{Array{Array{Complex_my,2}}}(undef,SpinP_switch)
  for spin in 1:SpinP_switch
    Hks_R_grouped[spin] = Array{Array{Complex_my,2}}(undef,0);
  end
  for spin in 1:SpinP_switch
    for i in 1:length(Hks_R[spin])

      HWR_mat = zeros(Complex_my, TotalOrbitalNum,TotalOrbitalNum);
      HWR_mat[:,:] = Hks_R[spin][i][wannierOrbital2atomGrouped,wannierOrbital2atomGrouped]
      push!(Hks_R_grouped[spin],HWR_mat);
    end
  end
  Gxyz = zeros(atomnum,3);
  for i in 1:atomnum
    Gxyz[i,:] += atompos[i,1]*tv[1,:] ;
    Gxyz[i,:] += atompos[i,2]*tv[2,:] ;
    Gxyz[i,:] += atompos[i,3]*tv[3,:] ;
  end

  Wannier_info = Wannierdatatype(atomnum, Total_NumOrbs, SpinP_switch, tv, rv,
    Gxyz, Hks_R, Hks_R_grouped, R_vector_mat, num_degen_list, ChemP,
    wannierOrbital2atomGrouped, 300.0, spin_type)
  return Wannier_info;

end
