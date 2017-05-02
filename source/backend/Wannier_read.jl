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

function read_wannier(wannier_fname::AbstractString,typeof_wannier::AbstractString,
  atoms_orbitals_list::Vector{Array{Int64}},atomnum::Int,atompos::Array{Float64,2})
  assert(atomnum == length(atoms_orbitals_list))
  return read_wannier_OpenMX(wannier_fname,atoms_orbitals_list,atomnum,atompos);
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

  ChemP = parse(Float64,split(lines[9])[end])
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
  if 2 == SpinP_switch
    #println(length(start_linenums))
    assert(0==rem(length(start_linenums),2))
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
