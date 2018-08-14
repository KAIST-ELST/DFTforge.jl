using ..DFTcommon

export EcalJscf
export read_EcalJ_scf

struct EcalJ_H
    orbital_cell_orbitals::Array{Int64,2};
    orbital_cell_indexs::Array{Int64,2};

    orbital_cell_Degen::Array{Int64,1}
    orbital_cell_H::Array{Complex_my,1};
    orbital_cell_S::Array{Complex_my,1};
end

struct EcalJscf
    atomnum::Int32
    Total_NumOrbs::Array{Int,1}
    SpinP_switch::Int

    tv::Array{Float64,2}
    rv::Array{Float64,2}
    qlat::Array{Float64,2}
    Gxyz::Array{Float64,2}

    ChemP::Float64


    E_Temp::Float64;
    spin_type::SPINtype

    EcalJ_H_list::Array{EcalJ_H,1}

    plat::Array{Float64,2}
end

function cal_colinear_eigenstate(k_point_frac_input,hamiltonian_info,spin_list);
# In EcalJ k_point_frac, k_point_cartesian should be seperated
# Internally k_point_cartesian is used
    scf_r = hamiltonian_info.scf_r;
    k_point_frac = [k_point_frac_input[1],k_point_frac_input[2],k_point_frac_input[3]];


    k_point_cartesian = scf_r.qlat'* k_point_frac;
    #k_point_cartesian = (k_point_frac'* scf_r.qlat)[:];
    TotalOrbitalNum = sum(scf_r.Total_NumOrbs);
    plat = scf_r.plat;

    kpoint_eigenstate_list =  Array{Kpoint_eigenstate}(0);
    for spin=spin_list
        Hk = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum)
        S  = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum)

        num_Hamil_lines = length(scf_r.EcalJ_H_list[spin].orbital_cell_H)
        @fastmath @simd for idx = 1:num_Hamil_lines

            alpha =  scf_r.EcalJ_H_list[spin].orbital_cell_orbitals[idx,1]
            beta  =  scf_r.EcalJ_H_list[spin].orbital_cell_orbitals[idx,2]
            #phase = (1./orbital_cell_Degen[idx])*exp(-im* 2.0*pi * sum( k_point .* (plat * orbital_cell_indexs[idx,3:5])'  ) )

            #phase = (1./scf_r.EcalJ_H_list[spin].orbital_cell_Degen[idx])*exp(-im* 2.0*pi *
            #    sum( k_point_cartesian .* ( plat * scf_r.EcalJ_H_list[spin].orbital_cell_indexs[idx,:]  )  ) )
            phase = (1. /scf_r.EcalJ_H_list[spin].orbital_cell_Degen[idx])*exp(-im* 2.0*pi *
                sum( k_point_cartesian .* ( plat[1,:] * scf_r.EcalJ_H_list[spin].orbital_cell_indexs[idx,1] +
                 plat[2,:] * scf_r.EcalJ_H_list[spin].orbital_cell_indexs[idx,2] +
                 plat[3,:] * scf_r.EcalJ_H_list[spin].orbital_cell_indexs[idx,3] )) )

            Hk[alpha,beta] += scf_r.EcalJ_H_list[spin].orbital_cell_H[idx]*phase;
            S[alpha,beta] += scf_r.EcalJ_H_list[spin].orbital_cell_S[idx]*phase
        end

        Hk *= Ry2eV;
#=
        # Cal S and S^(-1/2)
        # With Cut for small overlap
        S2 = copy(S);
        OLP_eigen_cut = 10.^-14.0
        OLP_eigen_cut = 10.^-20.0
        S_eigvals = zeros(Float_my, TotalOrbitalNum);
        # S^1/2
        eigfact_hermitian(S2,S_eigvals);
        if (!check_eigmat(S,S2,S_eigvals))
            println("S :",k_point)
        end
        #  S2 * 1.0/sqrt(S_eigvals[l])
        M1 = zeros(size(S_eigvals))
        M1[S_eigvals.>OLP_eigen_cut] = 1.0 ./sqrt.(S_eigvals[S_eigvals.>OLP_eigen_cut]);

        for j1 = 1:TotalOrbitalNum
            S2[:,j1] *= M1[j1];
        end
        #S2 = S2.';
=#
        # Wihtout Cut
        ##=
        Sq = sqrtm(S)
        S2 = inv(Sq);
        # =#

        Hk_tilta = S2' *Hk * S2
        Eigen_vect = copy(Hk_tilta)
        #Eigen_vect = copy( H )
        Eigen_value = zeros(TotalOrbitalNum);
        eigfact_hermitian(Eigen_vect,Eigen_value);

        Coes = S2 * Eigen_vect;
        #kpoint_common::Kpoint_eigenstate = Kpoint_eigenstate(Coes, Eigen_value, k_point_frac_input, Hk);
        kpoint_common::Kpoint_eigenstate = Kpoint_eigenstate(Eigen_vect, Eigen_value, k_point_frac_input, Hk_tilta);
        push!(kpoint_eigenstate_list,kpoint_common);

    end
    return kpoint_eigenstate_list;
end

function read_EcalJ_scf(result_file_dict::Dict{AbstractString,AbstractString},
    spin_type::DFTcommon.SPINtype)
    if (DFTcommon.colinear_type == spin_type)
        EcalJ_info_up = read_EcalJ_scf_interal(result_file_dict["result_file_up"])
        EcalJ_info_dn = read_EcalJ_scf_interal(result_file_dict["result_file_down"])

        Total_NumOrbs = EcalJ_info_up.Total_NumOrbs;
        assert(EcalJ_info_up.Total_NumOrbs == EcalJ_info_dn.Total_NumOrbs);
        #assert(EcalJ_info_up.EcalJ_H_list[1].orbital_cell_Degen
        #    == EcalJ_info_dn.EcalJ_H_list[1].orbital_cell_Degen)

        #tv = EcalJ_info_up.tv;
        assert(EcalJ_info_up.tv ==  EcalJ_info_dn.tv;)
        #rv = EcalJ_info_up.rv;
        assert(EcalJ_info_up.rv == EcalJ_info_dn.rv);
        #Gxyz = EcalJ_info_up.Gxyz;
        assert(EcalJ_info_up.Gxyz == EcalJ_info_dn.Gxyz)

        EcalJ_H_list = Array{EcalJ_H,1}(0);
        push!(EcalJ_H_list, deepcopy(EcalJ_info_up.EcalJ_H_list[1]) )
        push!(EcalJ_H_list, deepcopy(EcalJ_info_dn.EcalJ_H_list[1]) )

        EcalJ_info_updn = EcalJscf(EcalJ_info_up.atomnum,
        EcalJ_info_up.Total_NumOrbs, 2, EcalJ_info_up.tv, EcalJ_info_up.rv,
        EcalJ_info_up.qlat, EcalJ_info_up.Gxyz, EcalJ_info_up.ChemP,
            EcalJ_info_up.E_Temp, DFTcommon.colinear_type,
            EcalJ_H_list, EcalJ_info_up.plat )
        return EcalJ_info_updn;
    end


end

function read_EcalJ_scf_interal(scf_name::String)
    f = open(scf_name)
    lines = readlines(f)
    close(f)
    function search_next_delimiter(search_line)
        while(search_line < length(lines)  )
        if( contains(lines[search_line], "===="))
            break;
        end
        search_line+=1
        end
        return search_line
    end
    ChemP = parse(Float64,split(lines[2])[end])

    alat = parse(Float64,split(lines[4])[end])
    plat = zeros(3,3)
    plat[1,:] = map(x->parse(Float64,x),split(lines[6]))
    plat[2,:] = map(x->parse(Float64,x),split(lines[7]))
    plat[3,:] = map(x->parse(Float64,x),split(lines[8]))
    plat

    #tv =   plat * alat *  (1.0/ang2bohr) # borh -> ang
    tv =   plat * alat/ang2bohr
    rv = 2*pi*inv(tv') # Check required
    rv_inv = inv(tv')

    #tmp = map(x->parse(Int64,x),split(lines[10])[[2,4]]  )
    #atomnum = tmp[1]
    #TotalOrbitalNum = tmp[2]

    # Read Total_NumOrbs
    ib_table = Array{Int64}(0)

    #start_line = next_delimiter+1
    start_line = 12
    start_line = 10
    next_delimiter = search_next_delimiter(start_line)
    for current_line = start_line:next_delimiter-1
        tmp = map(x->parse(Int64,x),split(lines[current_line]))
        append!(ib_table,tmp)
    end
    atomnum = length(unique( ib_table))
    #assert( atomnum == length(unique( ib_table))) # Check atom num

    # Generate Total_NumOrbs
    Total_NumOrbs = Array{Int64}(atomnum)
    for atom_i = 1:atomnum
        Total_NumOrbs[atom_i] = sum( ib_table .== atom_i)
    end
    TotalOrbitalNum = sum(Total_NumOrbs);
    assert(TotalOrbitalNum == sum(Total_NumOrbs)) # Check TotalOrbitalNum
    println(Total_NumOrbs," ",TotalOrbitalNum)
    # npair
    start_line = next_delimiter + 1
    next_delimiter = search_next_delimiter(start_line)
    # qlat
    start_line =  next_delimiter + 1
    qlat = zeros(3,3)
    qlat[1,:] = map(x->parse(Float64,x),split(lines[start_line]))
    qlat[2,:] = map(x->parse(Float64,x),split(lines[start_line+1]))
    qlat[3,:] = map(x->parse(Float64,x),split(lines[start_line+2]))
    assert(sum(abs.(qlat - rv_inv * alat/ang2bohr )) < 10.0^-4)
    # Read atom Pos
    start_line = start_line+4
    #ib_table = Array{Int64}(0)


    next_delimiter = search_next_delimiter(start_line)
    atomnum_tmp = next_delimiter - start_line
    assert(atomnum == atomnum_tmp)

    Gxyz_cartesian = zeros(atomnum,3)
    Gxyz_frac = zeros(atomnum,3)
    cnt = 1
    for current_line = start_line:next_delimiter-1
        tmp = map(x->parse(Float64,x),split(lines[current_line]))
        Gxyz_cartesian[cnt,:] = tmp
        cnt += 1
    end
    Gxyz_cartesian

    for atom_i = 1:atomnum
        Gxyz_frac[atom_i,:] = qlat*Gxyz_cartesian[atom_i,:];
        #TODO: Check Gxyz_frac[atom_i,:] = inv(tv)*Gxyz_cartesian[atom_i,:];  or inv(tv')
    end
    # Read T vector numes
    start_line = next_delimiter+1
    next_delimiter = search_next_delimiter(start_line)

    # Scan all Cell vectors
    num_Hamil_lines = length(lines)- start_line + 1
    #orbital_cell_indexs =  zeros(Int64,num_Hamil_lines,5);
    orbital_cell_orbitals =  zeros(Int64,num_Hamil_lines,2);
    orbital_cell_indexs =  zeros(Int64,num_Hamil_lines,3);

    orbital_cell_Degen =  zeros(Int64,num_Hamil_lines);
    orbital_cell_H =  zeros(Complex_my, num_Hamil_lines);
    orbital_cell_S =  zeros(Complex_my, num_Hamil_lines);
    cnt = 1;
    for current_line = start_line:start_line+num_Hamil_lines - 1
        tmp = map(x->parse(Int64,x),split(lines[current_line])[1:6])

        orbital_cell_orbitals[cnt,:] = tmp[1:2];
        orbital_cell_indexs[cnt,:] = tmp[3:5];

        orbital_cell_Degen[cnt,:] = tmp[6]

        tmp2 = map(x->parse(Float64,x),split(lines[current_line])[7:10])

        orbital_cell_H[cnt] = tmp2[1] + im* tmp2[2]
        orbital_cell_S[cnt] = tmp2[3] + im* tmp2[4]
        cnt +=1;
    end
    cnt
    # Filterout small H & S 1.0E-14
    filter_small_H_S = .! ((abs.(orbital_cell_H) .< 1.0E-13) .& (abs.(orbital_cell_S) .< 1.0E-13))
    println(" Filtered ",
        (100.0*(1.0- (sum(filter_small_H_S)) / length(orbital_cell_H)) ),"% of H,S datas from #",
        length(orbital_cell_H),"  (Filtered abs(H) abs(S) are smaller than 1.0E-13)" )
    orbital_cell_orbitals = orbital_cell_orbitals[filter_small_H_S,:];
    orbital_cell_indexs = orbital_cell_indexs[filter_small_H_S,:];

    orbital_cell_Degen = orbital_cell_Degen[filter_small_H_S];
    orbital_cell_H = orbital_cell_H[filter_small_H_S];
    orbital_cell_S = orbital_cell_S[filter_small_H_S];

    # Store to EcalJ_H_list

    EcalJ_H_list = Array{EcalJ_H,1}(0);
    push!(EcalJ_H_list,EcalJ_H( orbital_cell_orbitals,orbital_cell_indexs,orbital_cell_Degen,orbital_cell_H,orbital_cell_S ) )

    E_Temp = 300.0
    EcalJ_info = EcalJscf(atomnum,Total_NumOrbs,1,tv,rv,qlat,Gxyz_cartesian,ChemP,
        E_Temp, DFTcommon.para_type,
        EcalJ_H_list, plat )
    return EcalJ_info;
end

function cal_colinear_eigenstate(k_point::k_point_Tuple,hamiltonian_info::Hamiltonian_info_type)

end
