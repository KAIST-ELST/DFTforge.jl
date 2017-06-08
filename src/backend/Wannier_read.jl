using DFTcommon
export Wannierdata
export read_wannier

immutable Wannierdatatype
  atomnum::Int
  Total_NumOrbs::Array{Int,1}
  SpinP_switch::Int

  tv::Array{Float64,2}
  rv::Array{Float64,2}
  Gxyz::Array{Float64,2}

  Hks_R_raw::Array{Array{Array{Complex_my,2}}}
  Hks_R::Array{Array{Array{Complex_my,2}}}
  R_vector_mat::Array{Array{Int,2}};
  ChemP::Float64

  wannierOrbital2atomGrouped::Array{Int,1}

  E_Temp::Float64;
end

function read_wannier(wannier_fname::AbstractString,result_file_dict::Dict{AbstractString,AbstractString},
  Wannier90_type::DFTcommon.Wannier90type, spin_type::DFTcommon.SPINtype,
  atoms_orbitals_list::Vector{Array{Int64}},atomnum::Int,atompos::Array{Float64,2})
  assert(atomnum == length(atoms_orbitals_list))
  if (DFTcommon.OpenMXWF == Wannier90_type)
    return read_wannier_OpenMX(wannier_fname,atoms_orbitals_list,atomnum,atompos);
  elseif (DFTcommon.EcalJWF  == Wannier90_type)
    return read_wannier_EcalJ(result_file_dict, spin_type, atoms_orbitals_list,atomnum,atompos);
  end
end

function read_wannier_OpenMX(wannier_fname,atoms_orbitals_list::Vector{Array{Int64}},
  atomnum::Int,atompos::Array{Float64,2})
  assert(atomnum == length(atoms_orbitals_list))
  Total_NumOrbs = Array{Int,1}(0);
  for i in 1:length(atoms_orbitals_list)
    push!(Total_NumOrbs,length(atoms_orbitals_list[i]));
  end

  f = open(wannier_fname)
  lines = readlines(f)
  close(f)

  tv = zeros(3,3)
  tv[1,:] = map(x->parse(Float64,x),split(lines[5]))
  tv[2,:] = map(x->parse(Float64,x),split(lines[6]))
  tv[3,:] = map(x->parse(Float64,x),split(lines[7]))


  rv = 2*pi*inv(tv'); # Check required

  ChemP = parse(Float64,split(lines[9])[end])*DFTcommon.Hatree2eV
  SpinP_switch  = parse(Int64,split(lines[8])[end])
  # read HWR
  current_line = 10
  start_linenums = Array{Int}(0);
  while (current_line <= length(lines))
    if (startswith(lines[current_line],"R"))
      #println(lines[current_line])
      push!(start_linenums,current_line)
    end
    current_line+= 1
  end
  linePer_R =  start_linenums[2:end] - start_linenums[1:end-1] - 1;
  assert(0 == sum(false .==(linePer_R[1] .== linePer_R)))
  assert(Float64(Int64(sqrt(linePer_R[1]))) == sqrt(linePer_R[1]))
  TotalOrbitalNum = Int64(sqrt(linePer_R[1]));
  assert(sum(Total_NumOrbs)==TotalOrbitalNum);
  wannierOrbital2atomGrouped = zeros(Int,TotalOrbitalNum)
  cnt = 1;
  for i in 1:length(atoms_orbitals_list)
    for j in 1:length(atoms_orbitals_list[i])
      wannierOrbital2atomGrouped[cnt] = atoms_orbitals_list[i][j];
      cnt+=1;
    end
  end

  R_vector_mat = Array{Array{Int,2}}(SpinP_switch)
  #R_vector_mat = zeros(length(start_linenums),3)
  Hks_R = Array{Array{Array{Complex_my,2}}}(SpinP_switch)
  for spin in 1:SpinP_switch
    Hks_R[spin] = Array{Array{Complex_my,2}}(0);
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
            #println(split(chop(lines[start_linenum+rel_line]) ))
    		tmp = map(x->parse(Float64,x),split(chop(lines[start_linenum+rel_line]))[3:4]);
    		#tmp = map(x->parse(Float64,x),split(chop(lines[start_linenums+rel_line]))[3:4])
    		value = tmp[1] + tmp[2]*im;
    		HWR_mat[x,y] = value;
	    end
    end
    HWR_mat *= DFTcommon.Hatree2eV;
    push!(Hks_R[spin],HWR_mat);
    cnt += 1;
  end
  # re arange wannier group same atoms orbitals
  Hks_R_grouped = Array{Array{Array{Complex_my,2}}}(SpinP_switch)
  for spin in 1:SpinP_switch
    Hks_R_grouped[spin] = Array{Array{Complex_my,2}}(0);
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

  Wannier_info = Wannierdatatype(atomnum,Total_NumOrbs,SpinP_switch,tv,rv,
    Gxyz,Hks_R,Hks_R_grouped,R_vector_mat,ChemP,
    wannierOrbital2atomGrouped,300.0)
  return Wannier_info;
end
function read_wannier_EcalJ(result_file_dict::Dict{AbstractString,AbstractString},
  spin_type::DFTcommon.SPINtype,
  atoms_orbitals_list::Vector{Array{Int64}},
  atomnum::Int,atompos::Array{Float64,2})

  if (DFTcommon.colinear_type == spin_type)
    assert(2 == length(result_file_dict) );
    Wannier_info_up = read_wannier_EcalJInternal(result_file_dict["result_file_up"],
      atoms_orbitals_list,
      atomnum,atompos)

    Wannier_info_down = read_wannier_EcalJInternal(result_file_dict["result_file_down"],
      atoms_orbitals_list,
      atomnum,atompos)

    Total_NumOrbs = Wannier_info_up.Total_NumOrbs;
    assert(Wannier_info_up.Total_NumOrbs == Wannier_info_down.Total_NumOrbs);
    tv = Wannier_info_up.tv;
    assert(Wannier_info_up.tv ==  Wannier_info_down.tv;)
    rv = Wannier_info_up.rv;
    assert(Wannier_info_up.rv == Wannier_info_down.rv);
    Gxyz = Wannier_info_up.Gxyz;
    assert(Wannier_info_up.Gxyz == Wannier_info_down.Gxyz)
    Hks_R_raw = Array{Array{Array{Complex_my,2}}}(2)
    Hks_R_raw[1] = Wannier_info_up.Hks_R_raw[1];
    Hks_R_raw[2] = Wannier_info_down.Hks_R_raw[1];

    Hks_R = Array{Array{Array{Complex_my,2}}}(2)
    Hks_R[1] = Wannier_info_up.Hks_R[1];
    Hks_R[2] = Wannier_info_down.Hks_R[1];


    #R_vector_mat = Wannier_info_up.R_vector_mat;
    R_vector_mat = Array{Array{Int,2}}(2)
    assert(Wannier_info_up.R_vector_mat == Wannier_info_down.R_vector_mat);
    R_vector_mat[1]  = Wannier_info_up.R_vector_mat[1]
    R_vector_mat[2]  = Wannier_info_down.R_vector_mat[1]

    ChemP = Wannier_info_up.ChemP;
    assert(0.0 == ChemP); # In EcalJ Wannier
    assert(Wannier_info_up.ChemP == Wannier_info_down.ChemP)

    wannierOrbital2atomGrouped = Wannier_info_up.wannierOrbital2atomGrouped;
    assert(Wannier_info_up.wannierOrbital2atomGrouped == Wannier_info_down.wannierOrbital2atomGrouped)

    Wannier_info = Wannierdatatype(atomnum,Total_NumOrbs,2,tv,rv,
        Gxyz,Hks_R_raw,Hks_R,R_vector_mat,ChemP,
        wannierOrbital2atomGrouped,300.0)
    return Wannier_info;
  elseif (DFTcommon.para_type == spin_type)
    return  read_wannier_EcalJInternal(result_file_dict["result_file_up"],
        atoms_orbitals_list,
        atomnum,atompos)
  end
end

function read_wannier_EcalJInternal(wannier_fname::AbstractString,
  atoms_orbitals_list::Vector{Array{Int64}},
  atomnum::Int,atompos::Array{Float64,2})

  assert(atomnum == length(atoms_orbitals_list))
  assert(atomnum == size(atompos)[1]);


  Total_NumOrbs = Array{Int,1}(0);
  for i in 1:length(atoms_orbitals_list)
    push!(Total_NumOrbs,length(atoms_orbitals_list[i]));
  end

  f = open(wannier_fname)
  lines = readlines(f)
  close(f)

  TotalOrbitalNum = parse(Int64,split(lines[5])[1]);
  assert(sum(Total_NumOrbs)==TotalOrbitalNum);
  num_Rvector = parse(Int64,split(lines[5])[2]);
  total_lines = parse(Int64,split(lines[5])[3]);
  assert(total_lines == TotalOrbitalNum^2 * num_Rvector);

  tv = zeros(3,3)
  tv[1,:] = map(x->parse(Float64,x),split(lines[2]))
  tv[2,:] = map(x->parse(Float64,x),split(lines[3]))
  tv[3,:] = map(x->parse(Float64,x),split(lines[4]))

  rv = 2*pi*inv(tv'); # Check required

  ChemP = 0.0 # TODO:Check for Wannier Chem potential parse(Float64,split(lines[9])[end])*DFTcommon.Hatree2eV
  SpinP_switch  = 1 #parse(Int64,split(lines[6])[end])

  # Check if atom infomations are correct
  check_fail = false;
  for (i,v) in enumerate(atoms_orbitals_list)
    position_1 = map(x->parse(Float64,x), split(lines[7+v[1]])[end-2:end]);
    #position_1 = map(x-> rem(x,1.0), position_1 + 2.0)
    if (sum(abs(atompos[i,:] - position_1)) > 10.0^-4.0)
      println(" read_wannier_EcalJ atom ",i," position not matched ",position_1,"  ",atompos[i,:])
      check_fail = true;
    end
    for (i2,v2) in enumerate(v)
      position_2 = map(x->parse(Float64,x), split(lines[7+v[i2] ])[end-2:end]);
      #position_2 = map(x-> rem(x,1.0), position_2 + 2.0)
      if (sum(abs(atompos[i,:] - position_2)) > 10.0^-4.0)
        println(" read_wannier_EcalJ atom ",i2," position not matched ",position_2,"  ",atompos[i,:])
        check_fail = true;
      end
    end
  end
  if check_fail
    println(" read_wannier_EcalJ atom position not matched ")
    println(DFTcommon.bar_string) # print ====...====
    assert(true);
  end
  # Read wannier H
  R_vector_mat = Array{Array{Int,2}}(SpinP_switch)
  Hks_R = Array{Array{Array{Complex_my,2}}}(SpinP_switch)
  wannierOrbital2atomGrouped = zeros(Int,TotalOrbitalNum)

  ## Check R vectors
  # TODO: Spin up down
  start_linenums = Array{Int}(0);
  start_line = 7 + TotalOrbitalNum + 1;
  R_vector_catesian = map(x->parse(Float64,x),split(lines[start_line])[4:6]);
  R_frac =  rv * R_vector_catesian[:];
  R_int = round(Int64,R_frac)
  # store current R_vect
  prev_R_vect = R_vector_catesian;
  push!(start_linenums, start_line);

  assert(sum(abs(R_frac - R_int)) < 10.^-4.0); # Check if cell vector is Int
  for current_line in start_line:length(lines)
    R_vector_catesian = map(x->parse(Float64,x),split(lines[current_line])[4:6]);
    if sum(abs(prev_R_vect-R_vector_catesian)) < 10.^-4.0
    else
      # new R vector
      push!(start_linenums,current_line);

      R_frac =  (rv) * R_vector_catesian[:]/(2*pi);
      R_int = round(Int64,R_frac)
      if (sum(abs(R_frac - R_int)) > 10.^-4.0)
        println(" Cartesian R vector  ",R_vector_catesian," did not convert to fractional R", R_frac," @line# ",current_line)
        assert(false); # Check if cell vector is Int
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
    Hks_R[1] = Array{Array{Complex_my,2}}(0);
    R_vector_mat[1] = zeros(length(start_linenums),3);
  #end
  spin = 1;

  cnt = 1;
  for (i,start_linenum) in enumerate(start_linenums)

    R_vector_catesian =  map(x->parse(Float64,x),split(lines[start_linenum])[4:6])
    R_frac =  (rv) * R_vector_catesian[:]/(2*pi);
    R_int = round(Int64,R_frac)
    if (sum(abs(R_frac - R_int)) > 10.^-4.0)
      println(" Cartesian R vector  ",R_vector_catesian," did not convert to fractional R", R_frac," @line# ",start_linenum)
      assert(false); # Check if cell vector is Int
    end
    R_vector_mat[spin][cnt,:] = R_int;

    HWR_mat = zeros(Complex_my, TotalOrbitalNum,TotalOrbitalNum);
    rel_line = 0;
    for x in 1:TotalOrbitalNum
      for y in 1:TotalOrbitalNum

        #println(split(chop(lines[start_linenum+rel_line]) ))
        tmp = map(x->parse(Int64,x),split(chop(lines[start_linenum+rel_line]))[7:8]);
        if ( (tmp[1] != x) || (tmp[2] != y) )
          println(" Orbital index l1,l2 (",x,",",y,") not matched (",tmp[1],",",tmp[2],")")
          assert(false);
        end

    		tmp = map(x->parse(Float64,x),split(chop(lines[start_linenum+rel_line]))[9:10]);
        rel_line += 1
    		#tmp = map(x->parse(Float64,x),split(chop(lines[start_linenums+rel_line]))[3:4])
    		value = tmp[1] + tmp[2]*im;
    		HWR_mat[x,y] = value;
	    end
    end
    #HWR_mat *= DFTcommon.Hatree2eV;
    push!(Hks_R[spin],HWR_mat);
    cnt += 1;
  end

  # re arange wannier group same atoms orbitals
  Hks_R_grouped = Array{Array{Array{Complex_my,2}}}(SpinP_switch)
  for spin in 1:SpinP_switch
    Hks_R_grouped[spin] = Array{Array{Complex_my,2}}(0);
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

  Wannier_info = Wannierdatatype(atomnum,Total_NumOrbs,SpinP_switch,tv,rv,
    Gxyz,Hks_R,Hks_R_grouped,R_vector_mat,ChemP,
    wannierOrbital2atomGrouped,300.0)
  return Wannier_info;

end
