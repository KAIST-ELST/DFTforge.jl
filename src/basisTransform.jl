###############################################################################
# Hongkee Yoon Hongkeeyoon@kaist.ac.kr
# 2019.05
# https://kaist-elst.github.io/DFTforge.jl/
###############################################################################



#module basisTransform
using Rotations # Pkg.add("Rotations")
using DataStructures
#using DFTcommon
using LinearAlgebra

export rot_matrixByZXaxisbase,rot_D_orbital
export basisTransform_rule_type,orbital_rot_type
export rot_basis!,Heff
export ncHamiltonian_rot!
export RotPauliX,RotPauliY,RotPauliZ,RotPauliZYZ,RotPauliZX,RotPauliZXZ,RotPauliThetaPhi

################################################################################


function RotPauliX(theta)
	rot = zeros(ComplexF64  ,2,2);
	rot[1,1] = cos(theta/2);
	rot[2,2] = cos(theta/2);
	rot[1,2] = -im*(sin(theta/2));
	rot[2,1] = -im*(sin(theta/2));
	return rot;
end

function RotPauliY(theta)
	rot = zeros(ComplexF64  ,2,2);
	rot[1,1] = cos(theta/2);
	rot[2,2] = cos(theta/2);
	rot[1,2] = (sin(theta/2));
	rot[2,1] = -(sin(theta/2));
	return rot;
end

function RotPauliZ(theta)
	rot = zeros(ComplexF64  ,2,2);
	rot[1,1] = exp(-im*(theta/2))
	rot[2,2] = exp(im*(theta/2))
	rot[1,2] = 0;
	rot[2,1] = 0;
	return rot;
end

function RotPauliZYZ(alpha::Float_my,beta::Float_my,gamma::Float_my)
  rotZYZ = RotPauliZ(alpha)*RotPauliY(beta)*RotPauliZ(gamma);
  return rotZYZ;
end
function RotPauliZXZ(alpha::Float_my,beta::Float_my,gamma::Float_my)
  rotZYZ = RotPauliZ(alpha)*RotPauliX(beta)*RotPauliZ(gamma);
  return rotZYZ;
end

function RotPauliZX(theta::Float_my,phi::Float_my)
  rot_theta_phi = RotPauliZ(theta)*RotPauliX(phi);
  return rot_theta_phi;
end

function RotPauliThetaPhi(theta::Float_my,phi::Float_my)
	cos_theta = cos(theta/2.0);
	sin_theta = sin(theta/2.0);
	exp_nphi = exp(-im*0.5*phi);
	exp_phi = exp(im*0.5*phi);
	rot = zeros(ComplexF64  ,2,2);
	rot[1,1] = cos_theta*exp_nphi;
	rot[2,2] = cos_theta*exp_phi;
	rot[1,2] = -sin_theta * exp_nphi;
	rot[2,1] = sin_theta * exp_phi;
	return rot
end

function ncHamiltonian_rot!(H::Array{Complex_my,2},orbital_indexs,
  theta::Float_my, phi::Float_my)
  TotalOrbitalNum2 =  size(H)[1];
  @assert(size(H)[2] == TotalOrbitalNum2);
  @assert(TotalOrbitalNum2 % 2 == 0); # Dimension must be even Number
  TotalOrbitalNum = convert(Int64, TotalOrbitalNum2/2);
  @assert(maximum(orbital_indexs) <= TotalOrbitalNum);
  rot_theta_phi = RotPauliZ(phi)*RotPauliX(theta);
  @assert(2==size(rot_theta_phi)[1]);@assert(2==size(rot_theta_phi)[2]);

  Rot_mat_allbasis = eye(Complex_my,TotalOrbitalNum2,TotalOrbitalNum2);
  #Rot_mat_allbasis[TotalOrbitalNum+1:TotalOrbitalNum2, 1:TotalOrbitalNum] = eye(Complex_my,TotalOrbitalNum,TotalOrbitalNum);
  #Rot_mat_allbasis[1:TotalOrbitalNum, TotalOrbitalNum+1:TotalOrbitalNum2] = eye(Complex_my,TotalOrbitalNum,TotalOrbitalNum);
  for oribtal_index in orbital_indexs
    Rot_mat_allbasis[oribtal_index,oribtal_index] = rot_theta_phi[1,1];
    Rot_mat_allbasis[TotalOrbitalNum+oribtal_index,TotalOrbitalNum+oribtal_index] = rot_theta_phi[2,2];
    #Rot_mat_allbasis[oribtal_index,TotalOrbitalNum+oribtal_index] = rot_theta_phi[1,2];
    #Rot_mat_allbasis[TotalOrbitalNum+oribtal_index,oribtal_index] = rot_theta_phi[2,1];
  end
  H[:,:]=(Rot_mat_allbasis')*H*(Rot_mat_allbasis);
end

function ncHamiltonian_rot!(H::Array{Complex_my,2},orbital_indexs,
  rot_pauli_mat::Array{ComplexF64 ,2})
  TotalOrbitalNum2 =  size(H)[1];
  @assert(size(H)[2] == TotalOrbitalNum2);
  @assert(TotalOrbitalNum2 % 2 == 0); # Dimension must be even Number
  TotalOrbitalNum = convert(Int64, TotalOrbitalNum2/2);
  @assert(maximum(orbital_indexs) <= TotalOrbitalNum);
  #rot_zyz = RotPauliZYZ(alpha, beta, gamma);
  @assert(2==size(rot_pauli_mat)[1]);@assert(2==size(rot_pauli_mat)[2]);

  Rot_mat_allbasis = eye(Complex_my,TotalOrbitalNum2,TotalOrbitalNum2);
	Rot_mat_allbasis[TotalOrbitalNum+1:TotalOrbitalNum2, 1:TotalOrbitalNum] = eye(Complex_my,TotalOrbitalNum,TotalOrbitalNum);
  Rot_mat_allbasis[1:TotalOrbitalNum, TotalOrbitalNum+1:TotalOrbitalNum2] = eye(Complex_my,TotalOrbitalNum,TotalOrbitalNum);
  for oribtal_index in orbital_indexs
    Rot_mat_allbasis[oribtal_index,oribtal_index] = rot_pauli_mat[1,1];
    Rot_mat_allbasis[TotalOrbitalNum+oribtal_index,TotalOrbitalNum+oribtal_index] = rot_pauli_mat[2,2];
    Rot_mat_allbasis[oribtal_index,TotalOrbitalNum+oribtal_index] = rot_pauli_mat[1,2];
    Rot_mat_allbasis[TotalOrbitalNum+oribtal_index,oribtal_index] = rot_pauli_mat[2,1];
  end
  H[:,:] = (Rot_mat_allbasis')*H*Rot_mat_allbasis;
end


struct orbital_rot_type
  atom1::Int
  orbital_rot::Array{Float64,2}
  orbital_list::Array{Array{Int}}
  #atomnum::Int
  #Orbital_rot::Array{Float64,2}
  #orbital_list::Array{Array{Int}} #index of d-orbitals z^2, x2y2, xy, xz, yz
  orbital_num::Int # for d-orbital 5
  duplicated_orbital_num::Int # ex) In OpenMX LCAO with s2p2d2f1 : d-orbital 2
  #Z_vect::Vector{Float64} # for debug
  #X_vect::Vector{Float64} # for debug

  function orbital_rot_type(atomnum,Z_vect,X_vect,d_orbital_list,orbital_num::Int)
    @assert(3==length(Z_vect))
    @assert(3==length(X_vect))
    @assert(orbital_num==length(d_orbital_list))
    duplicated_orbital_num = length(d_orbital_list[1])
    for (k,v) in enumerate(d_orbital_list)
      @assert(length(v) == duplicated_orbital_num);
    end
    R = rot_matrixByZXaxisbase(Z_vect,X_vect);
    orbital_rot_d = rot_D_orbital(R);
    new(atomnum,orbital_rot_d,d_orbital_list,orbital_num,duplicated_orbital_num)
	#,Z_vect,X_vect)
  end
  function orbital_rot_type(atomnum,R::Array{Float64,2},d_orbital_list,orbital_num::Int)
	#@assert(3==length(Z_vect))
	#@assert(3==length(X_vect))
	@assert(orbital_num==length(d_orbital_list))
	duplicated_orbital_num = length(d_orbital_list[1])
	for (k,v) in enumerate(d_orbital_list)
	  @assert(length(v) == duplicated_orbital_num);
	end
	#R = rot_matrixByZXaxisbase(Z_vect,X_vect);
	orbital_rot_d = rot_D_orbital(R);
	new(atomnum,orbital_rot_d,d_orbital_list,orbital_num,duplicated_orbital_num)#,
	#[0.0,0.0,1.0],[1.0,0.0,0.0])
  end
end


function rot_basis!(mat_original,orbitalStartIdx::Array{Int},
  orbital_rot_rules::Dict{Int,orbital_rot_type})
  for (k,v) in orbital_rot_rules
    atom1 = v.atom1;
    orbital_num = v.orbital_num;
    duplicated_orbital_num = v.duplicated_orbital_num;
    @assert(k == atom1);
    orbital_rot = v.orbital_rot;
    @assert(Array{Float64,2} == typeof(orbital_rot));
    for duplicated_orbital_i in 1:duplicated_orbital_num
      rel_orbital_list = map(x->x[duplicated_orbital_i], v.orbital_list);
      atom_orbitals = orbitalStartIdx[atom1] .+ rel_orbital_list
      # gen orbital_rot_full
      #orbital_rot_full = eye(size(mat_original)[1]);
	  orbital_rot_full = Matrix(1.0I, size(mat_original)[1], size(mat_original)[1]);
      orbital_rot_full[atom_orbitals,atom_orbitals] = orbital_rot
      # actual rotation
      # Aroted = U * A * A'
      mat_original[:,:] = orbital_rot_full'*mat_original*(orbital_rot_full);
    end
  end
  #return mat_original;
end

function rot_matrixByZXaxisbase(Z_vect::Vector{Float64},X_vect::Vector{Float64})
  Z_vect = normalize(Z_vect)
  X_vect = normalize(X_vect)
  Y_vect = normalize(cross(Z_vect,X_vect))
  R = zeros(3,3)

  R[1,:] = X_vect
  R[2,:] = Y_vect
  R[3,:] = Z_vect


  return collect(R');
end
#=
function rot_matrixByZYXangle(theta_Z::Float64,theta_Y::Float64::Float64)
  Rz1 = RotZ(theta_Z)
  Ry =  RotY(theta_Y)
  Rz2 = RotZ(phi)
  R = Rz2 * Ry * Rz1
  R = convert(Array{Float64,2},R)

end
=#

function rot_D_orbital(R::Array{Float64,2})
  # From Jae-Hoon Sim's orbital rotation code
  diff = ( abs(opnorm(R,2)-1.0)) # Julia 0.6 diff = ( abs(norm(R)-1.0))
  if !( diff < 10.0^-4.0)
	  println(" new global X, Z are not orthorgonal ", diff)
	  @assert(false)
  end

  h = 2.0
  s3 = sqrt(3.0)
  DR = Matrix(1.0I,5,5);

  DR[1,1]=3/h*(R[3,3])^2.0-1/h
  DR[1,2]=s3/h*(R[3,1]^2.0-R[3,2]^2.0)
  DR[1,3]=s3*(R[3,1]*R[3,2])
  DR[1,4]=s3*(R[3,1]*R[3,3])
  DR[1,5]=s3*(R[3,2]*R[3,3])

  DR[2,1]=s3/h*(R[1,3]^2.0-R[2,3]^2.0)
  DR[2,2]=1/h*(R[1,1]^2.0-R[2,1]^2.0-R[1,2]^2.0+R[2,2]^2.0)
  DR[2,3]=(R[1,1]*R[1,2]-R[2,1]*R[2,2])
  DR[2,4]=(R[1,1]*R[1,3]-R[2,1]*R[2,3])
  DR[2,5]=(R[1,2]*R[1,3]-R[2,2]*R[2,3])

  DR[3,1]=s3*(R[1,3]*R[2,3])
  DR[3,2]=(R[1,1]*R[2,1]-R[1,2]*R[2,2])
  DR[3,3]=R[1,1]*R[2,2]+R[1,2]*R[2,1]
  DR[3,4]=R[1,1]*R[2,3]+R[1,3]*R[2,1]
  DR[3,5]=R[1,2]*R[2,3]+R[1,3]*R[2,2]

  DR[4,1]=s3*(R[1,3]*R[3,3])
  DR[4,2]=(R[1,1]*R[3,1]-R[1,2]*R[3,2])
  DR[4,3]=R[1,1]*R[3,2]+R[3,1]*R[1,2]
  DR[4,4]=R[1,1]*R[3,3]+R[1,3]*R[3,1]
  DR[4,5]=R[1,2]*R[3,3]+R[1,3]*R[3,2]

  DR[5,1]=s3*(R[2,3]*R[3,3])
  DR[5,2]=(R[2,1]*R[3,1]-R[2,2]*R[3,2])
  DR[5,3]=R[2,1]*R[3,2]+R[3,1]*R[2,2]
  DR[5,4]=R[2,1]*R[3,3]+R[3,1]*R[2,3]
  DR[5,5]=R[2,2]*R[3,3]+R[2,3]*R[3,2]
  @assert( abs(opnorm(DR,2)-1.0) < 10.0^-4.0); # julia 0.6 assert( abs(norm(DR)-1.0) < 10.0^-4.0);

  return DR;
end


struct orbital_merge_type
  atom1::Int
  rel_orbital2merge::Array{Array{Int}}
end

struct basisTransform_rule_type
  orbital_rot_on::Bool;
  orbital_rot_rules::Dict{Int,orbital_rot_type}; # orbital_rot_d::orbital_rot_d_type

  orbital_merge_on::Bool;
  orbital_merge_rules::Dict{Int,orbital_merge_type}
  keep_unmerged_atoms::Bool
  keep_unmerged_orbitals::Bool
  function basisTransform_rule_type()
    new(false,Dict{Int,orbital_rot_type}(),false,
      Dict{Int,orbital_merge_type}(),true,true)
  end
  function basisTransform_rule_type(orbital_rot_on,orbital_rot_rules,orbital_merge_on,
    orbital_merge_rules,keep_unmerged_atoms,keep_unmerged_orbitals)
    new(orbital_rot_on,orbital_rot_rules,orbital_merge_on,
      orbital_merge_rules,keep_unmerged_atoms,keep_unmerged_orbitals)
  end
end


export orbital_merge_type,basisTransform_rule_type,basisTransform_result_type

struct basisTransform_result_type
  atomnum::Int
  orbitalNums::Array{Int}
  orbitalStartIdx_list::Array{Int}

  orbital_index_orig2new::SortedDict{Int,Dict{Int,Int}};
  survieved_orbitals_dict::SortedDict{Int,Array{Int}};
  unsurvieved_orbitals_dict::SortedDict{Int,Array{Int}};

  function basisTransform_result_type(atomnum::Int,  orbitalNums::Array{Int},
    orbital_index_orig2new::SortedDict{Int,Dict{Int,Int}},
    survieved_orbitals_dict::SortedDict{Int,Array{Int}},
    unsurvieved_orbitals_dict::SortedDict{Int,Array{Int}})

    orbitalStartIdx_list = Array{Int}(undef,atomnum);
    orbitalStartIdx = 0;
    @assert(length(orbitalNums) == atomnum);
    for (k,v) in enumerate(orbitalNums)
      orbitalStartIdx_list[k] = orbitalStartIdx
      orbitalStartIdx += v;
    end
    @assert(length(orbitalStartIdx_list) == atomnum);
    println(" orbitalStartIdx_list ", orbitalStartIdx_list)
    #TODO: atomnum,orbitalNums,orbitalStartIdx_list could be changed if orbitals are merged
    new(atomnum,orbitalNums,orbitalStartIdx_list,
    orbital_index_orig2new,survieved_orbitals_dict,unsurvieved_orbitals_dict)
  end

end

export basisTransform_init

function basisTransform_init(atomnum::Int,orbitalNums::Array{Int},basisTransform::basisTransform_rule_type)
  println(" basisTransform_init ")
  # Default setting
  @assert(atomnum > 0)
  @assert(length(orbitalNums) == atomnum)
  atomnum_eff = copy(atomnum)
  orbitalNums_eff = copy(orbitalNums)

  orbital_index_orig2new = Dict{Int,Dict{Int,Int}}();
  survieved_orbitals_dict = Dict{Int,Array{Int}}();
  unsurvieved_orbitals_dict = Dict{Int,Array{Int}}();
  for (k,atom1_orig) in enumerate(1:atomnum)
    orbital_index_orig2new[atom1_orig] = Dict( i => 1 for i = 1:orbitalNums[atom1_orig]);
    survieved_orbitals_dict[atom1_orig] = collect(1:orbitalNums[atom1_orig]);
    unsurvieved_orbitals_dict[atom1_orig] = Array{Int}(undef,0);
  end
  #
  if basisTransform.orbital_merge_on

    orbital_index_new = Dict{Int,Int}();
    orbital_index_orig2new = Dict{Int,Dict{Int,Int}}();
    survieved_orbitals_dict = Dict{Int,Array{Int}}();
    unsurvieved_orbitals_dict = Dict{Int,Array{Int}}();

    atom_survived = Array{Int}(undef,0);

    for (k,orbital_merge_rule) in basisTransform.orbital_merge_rules
      atom1_orig = orbital_merge_rule.atom1;
      push!(atom_survived,atom1_orig)

      atom1_index_orig2new = Dict{Int,Int}();#Dict( i => i for i = 1:orbitalNums[atom1_orig]);
      # unspecified orbitals : survieved_orbitals (set 1)
      atom1_orbitals_all_dict =  Dict( i => 1 for i = 1:orbitalNums[atom1_orig]);
      if basisTransform.keep_unmerged_orbitals
      else
        # unspecified orbitals : unsurvieved_orbitals (set 0)
        atom1_orbitals_all_dict =  Dict( i => 0 for i = 1:orbitalNums[atom1_orig]);
      end

      # mark survieved_orbitals and unsurvived orbitals
      # all atom1 orbital are survieved only unsurvived orbitals are marked
      for (k2,rel_orbital2merge) in enumerate(orbital_merge_rule.rel_orbital2merge)
        # representative orbital for unsurvieved_orbitals
        atom1_orbitals_all_dict[rel_orbital2merge[1]] = 2;

        for rel_orbital in rel_orbital2merge[2:end]
          atom1_orbitals_all_dict[rel_orbital] = -1; # unsurvieved_orbitals
        end
      end
      # group  survieved_orbitals & unsurvieved_orbitals
      survieved_orbitals = Array{Int}(undef,0);
      unsurvieved_orbitals = Array{Int}(undef,0);
      atom1_orbital_cnt = 0;
      for atom1_rel_oribtal in 1:orbitalNums[atom1_orig]
        if atom1_orbitals_all_dict[atom1_rel_oribtal] > 0
          atom1_orbital_cnt += 1
          atom1_index_orig2new[atom1_rel_oribtal] = atom1_orbital_cnt;
          push!(survieved_orbitals,atom1_rel_oribtal);
        else
          push!(unsurvieved_orbitals,atom1_rel_oribtal);
        end
      end
      # unsurvieved_orbitals are maped in to representative orbitals
      for (k2,rel_orbital2merge) in enumerate(orbital_merge_rule.rel_orbital2merge)
        for rel_orbital in rel_orbital2merge[2:end]
          atom1_index_orig2new[rel_orbital] = atom1_index_orig2new[rel_orbital2merge[1]]
        end
      end
      @assert(length(keys(atom1_index_orig2new))  == orbitalNums[atom1_orig])
      #atom1_index_orig2new[rel_orbital] = rel_orbital2merge[1];

      orbital_index_orig2new[atom1_orig] = atom1_index_orig2new;
      survieved_orbitals_dict[atom1_orig] = survieved_orbitals
      unsurvieved_orbitals_dict[atom1_orig] = unsurvieved_orbitals
    end

    for atom1_orig in 1:atomnum
      if !haskey(orbital_index_orig2new,atom1_orig)
        survieved_orbitals = collect(1:orbitalNums[atom1_orig]);
        unsurvieved_orbitals = Array{Int}(undef,0);
        atom1_index_orig2new = Dict(i=>i for i in 1:orbitalNums[atom1_orig]);
        if !basisTransform.keep_unmerged_atoms
          survieved_orbitals = Array{Int}(undef,0);
          unsurvieved_orbitals = collect(1:orbitalNums[atom1_orig]);
          atom1_index_orig2new = Dict( i => -1 for i in 1:orbitalNums[atom1_orig])
        end

        orbital_index_orig2new[atom1_orig] = atom1_index_orig2new;
        survieved_orbitals_dict[atom1_orig] = survieved_orbitals
        unsurvieved_orbitals_dict[atom1_orig] = unsurvieved_orbitals
      end
    end

    atomnum_eff = length(orbital_index_new);
    orbitalNums_eff = Array{Int}(undef,0);
    for (k,v) in orbital_index_new
      push!(orbitalNums_eff,v)
    end
  end
  orbital_index_orig2new = SortedDict(orbital_index_orig2new);
  survieved_orbitals_dict = SortedDict(survieved_orbitals_dict);
  unsurvieved_orbitals_dict = SortedDict(unsurvieved_orbitals_dict);

  println(atomnum_eff, orbitalNums_eff)

  return basisTransform_result_type(atomnum_eff,orbitalNums_eff,
    orbital_index_orig2new,survieved_orbitals_dict,unsurvieved_orbitals_dict)
end

function Heff(H,orbitalStartIdx,basisTransform_rule::basisTransform_rule_type,w)
    #FHeff[H00_, H01_, H11_, H10_, w_] := (H00 + H01.Inverse[w - H11].H10)
    (TotalOrbitalNum_1,TotalOrbitalNum_2)=size(H)
    @assert(TotalOrbitalNum_1 == TotalOrbitalNum_2)
    TotalOrbitalNum2 = TotalOrbitalNum_1;


    ## Rot orbitals
    if (basisTransform_rule.orbital_rot_on)
      rot_basis!(H,orbitalStartIdx,basisTransform_rule.orbital_rot_rules)
    end
    Heff = H;
    #=
    ## Merge orbitals

    unselected_orbitals = collect(1:TotalOrbitalNum2);
    survieved_orbitals = Array{Int}(undef,0);
    for (k,v) in enumerate(basisTransform.orbital_merge_rules)
        atomnum = v.atomnum
        rel_orbital2merge = v.rel_orbital2merge;
        orbitals2merge::Array{Int} = orbitalNums[atomnum] + rel_orbital2merge
        push!(survieved_orbitals,orbitals2merge);
        deleteat!(unselected_orbitals,v);
    end
    @assert( 0 < length(survieved_orbitals))

    ##

    H00 = H[survieved_orbitals,survieved_orbitals];
    H11 = H[unselected_orbitals,unselected_orbitals];
    H01 = H[survieved_orbitals,unselected_orbitals];
    H10 = H[unselected_orbitals,survieved_orbitals];

    Heff = H00 + H01 * (inv(H11)-w)* H10
    =#
    return Heff
end

struct energyCut_type
  E_low_window::Float_my
  E_upper_window::Float_my
end

#end
