###############################################################################
# Hongkee Yoon Hongkeeyoon@kaist.ac.kr
# 2019.05
# https://kaist-elst.github.io/DFTforge.jl/
###############################################################################


__precompile__(true)
using ..DFTcommon
using LinearAlgebra
#using DMFTscf
import JSON
export read_lobster, cal_Hamiltonian, cal_colinear_eigenstate

H_type = Array{Array{Array{ComplexF64,2}}}
struct Lobsterdatatype
    atomnum::Int32
    SpinP_switch::Int32

    #atv::Array{Float64,2}
    Total_NumOrbs::Array{Int32,1}


    tv::Array{Float64,2}
    rv::Array{Float64,2}

    Gxyz::Array{Float64,2}

    R_vector_mat::Array{Array{Int,2}}
    Hks_R::H_type

    R_vector_Overlap_mat::Array{Array{Int,2}}
    OverlapS_R::H_type


    ChemP::Float64
    spin_type::SPINtype
    E_Temp::Float64
end


function Overlap_Band!(H::H_type, R_vector_mat, spin::Int,
    Hout::Array{Complex_my,2},
    TotalOrbitalNum::Int64,k1,k2,k3)
    #HWR_mat_list::Array{Array{Complex_my,2}}
    FNAN = length(H[spin]);
    @assert(size(R_vector_mat[spin])[1] == FNAN)
    k_point::Array{Float_my,1} = [k1,k2,k3];
    #println(FNAN," spin ", 1)
    for LB_AN = 1:FNAN
        kRn::Float_my = -sum(R_vector_mat[spin][LB_AN,:].*k_point);
        Hout[:,:] += H[spin][LB_AN].* (cos(2.0*pi*kRn)+sin(2.0*pi*kRn)*im);
    end
end

function cal_Hamiltonian(k_point::k_point_Tuple,   # lobster_r::Wannierdatatype,
    hamiltonian_info::Hamiltonian_info_type, spin::Int)
    lobster_r = hamiltonian_info.scf_r
    spin_type = hamiltonian_info.spin_type::DFTcommon.SPINtype

    TotalOrbitalNum = sum(lobster_r.Total_NumOrbs[:])::Int
    TotalOrbitalNum2 = TotalOrbitalNum;
    if DFTcommon.non_colinear_type == spin_type
        TotalOrbitalNum2 = TotalOrbitalNum*2;
    end

    orbitalStartIdx_list = zeros(Int,lobster_r.atomnum)
    orbitalStartIdx = 0
    for i = 1:lobster_r.atomnum
        orbitalStartIdx_list[i] = orbitalStartIdx;
        orbitalStartIdx += lobster_r.Total_NumOrbs[i]
    end
    Hout = zeros(Complex_my, TotalOrbitalNum2, TotalOrbitalNum2);

    Overlap_Band!(lobster_r.Hks_R, lobster_r.R_vector_mat, spin, Hout, TotalOrbitalNum2, k_point[1], k_point[2], k_point[3])
    if hamiltonian_info.basisTransform_rule.orbital_rot_on
        if (DFTcommon.non_colinear_type == spin_type)
            throw(assertionError("Non-collinear spin basis rotation not supported yet "));
        end
        Hout = Heff(Hout,orbitalStartIdx_list, hamiltonian_info.basisTransform_rule, hamiltonian_info.basisTransform_result, 0.0);
    #println( sum(abs(H2-H)) )
    end
    return Hout;
end
function cal_colinear_eigenstate(k_point::k_point_Tuple,hamiltonian_info::Hamiltonian_info_type,spin_list::Array{Int})
    lobster_r = hamiltonian_info.scf_r
    spin_type = hamiltonian_info.spin_type;

    ## Overlapmaxitrx S

    TotalOrbitalNum = sum(lobster_r.Total_NumOrbs[:]);
    S = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum)
    orbitalStartIdx_list = zeros(Int,lobster_r.atomnum)
    orbitalNums = copy(lobster_r.Total_NumOrbs);
    #MPF = Array(Int,lobster_r.atomnum)
    orbitalStartIdx = 0
    for i = 1:lobster_r.atomnum
        orbitalStartIdx_list[i] = orbitalStartIdx;
        orbitalStartIdx += orbitalNums[i]
    end

    Overlap_Band!(lobster_r.OverlapS_R, lobster_r.R_vector_Overlap_mat, 1, S, TotalOrbitalNum, k_point[1], k_point[2], k_point[3])

    # S rotation  & Merge
    if hamiltonian_info.basisTransform_rule.orbital_rot_on
        S = Heff(S,orbitalStartIdx_list,hamiltonian_info.basisTransform_rule,
                hamiltonian_info.basisTransform_result,0.0);
    end
    #Sq = sqrtm(S)
    #S2 = inv(Sq);
    kpoint_eigenstate_list =  Array{Kpoint_eigenstate}(undef,0);
    for spin=spin_list #lobster_r.SpinP_switch+1

        H = zeros(Complex_my,TotalOrbitalNum,TotalOrbitalNum);
        H[:,:] = cal_Hamiltonian(k_point, hamiltonian_info, spin)

        Eigen_vect::Array{Complex{Float64},2},Eigen_value::Array{Float64,1}, Hk_tilta::Array{Complex{Float64},2} = cal_eigenVectVal(H,S, OLP_eigen_cut=1e-4)

        kpoint_common::Kpoint_eigenstate = Kpoint_eigenstate(Eigen_vect, Eigen_value,k_point,Hk_tilta); # Lowdin
        #kpoint_common::Kpoint_eigenstate = Kpoint_eigenstate(Coes,ko,k_point,H);
        push!(kpoint_eigenstate_list,kpoint_common);
    end

    return kpoint_eigenstate_list;
end

function read_lobster(lobster_dirname::AbstractString,spin_type::SPINtype,
    atomnum::Int64, atompos::Array{Float64,2}, tv::Array{Float64,2}, ChemP::Float64)

    rv = 2*pi*inv((tv')); # Check required

    Gxyz = zeros(atomnum,3);
    for i in 1:atomnum
      Gxyz[i,:] += atompos[i,1]*tv[1,:] ;
      Gxyz[i,:] += atompos[i,2]*tv[2,:] ;
      Gxyz[i,:] += atompos[i,3]*tv[3,:] ;
    end
    ###################################################
    if(isdirpath(lobster_dirname))
        println("lobster_dirname ",lobster_dirname)
    end
    function search_next_delimiter(search_line, lines, target; maxscan = -1)
        scan_cnt = 0
        found = false
        while(search_line < length(lines)  )
            if (occursin(target , lines[search_line]))
                found = true
                break;
            end
            search_line+=1
            scan_cnt += 1
            if (0 < maxscan & maxscan < scan_cnt )
                print("Reached scan_cnt line ", search_line)
                break;
            end
        end
        if found
            return search_line
        end
        return -1 # not found
    end
    ###################################################
    lobsterout_fname = joinpath(lobster_dirname,"lobsterout")
    #readlines(lobsterout_fname)
    # 5 detecting used PAW program... VASP


    hamiltonMatrices_fname = joinpath(lobster_dirname,"RealSpaceHamiltonians.lobster")
    overlapMatrices_fname = joinpath(lobster_dirname,"RealSpaceOverlaps.lobster")

    hamiltonMatrices_fname_text  = readlines(hamiltonMatrices_fname);

    println(" Total Hamiltonian length: ",length(hamiltonMatrices_fname_text))
    for i in 1:9
        println(hamiltonMatrices_fname_text[i])
    end
    if DFTcommon.non_colinear_type == spin_type
        # TODO
    end
    SpinP_switch = 2 # TODO: currently support only for colinear type

    # num_R_vector
    tmp_str_list = split(hamiltonMatrices_fname_text[1])
    num_R_vector = parse(Int64, tmp_str_list[4])
    R_vector_cnt = zeros(Int,SpinP_switch)
    println("num_R_vector ",num_R_vector)
    # TotalOrbitalNum
    tmp_str_list = split(hamiltonMatrices_fname_text[2])
    TotalOrbitalNum   = parse( Int64, tmp_str_list[end-2])
    TotalOrbitalNum_0 = parse( Int64, tmp_str_list[end])
    @assert(TotalOrbitalNum == TotalOrbitalNum_0)
    println("TotalOrbitalNum ", TotalOrbitalNum)

    # Total_NumOrbs
    tmp_str_list = split(hamiltonMatrices_fname_text[9])[2:end]

    @assert(length(tmp_str_list) == TotalOrbitalNum)
    #Atom_names = Array{String,1}()
    Atom_names = [split(i,"_")[1] for i in tmp_str_list]
    Orbital_names = [split(i,"_")[2:end] for i in tmp_str_list]
    Atom_unique_names = unique(Atom_names)

    atomnum = length(Atom_unique_names)
    Total_NumOrbs = zeros(Int, atomnum)
    for (k,v) in enumerate(tmp_str_list)
        #println(k," ", v)
        atom_name = split(v,"_")[1]
        atom_i = findfirst(x -> x == atom_name, Atom_unique_names)

        Total_NumOrbs[atom_i] += 1
    end


    println("Total_NumOrbs ",Total_NumOrbs)
    @assert(sum(Total_NumOrbs) == TotalOrbitalNum)


    ###################################################
    ## HK_R matrix
    ###################################################

    R_vector_mat = Array{Array{Int,2}}(undef,SpinP_switch)
    Hks_R = Array{Array{Array{ComplexF64,2}}}(undef,SpinP_switch)

    for spin in 1:SpinP_switch
        Hks_R[spin]      = Array{Array{ComplexF64,2}}(undef,num_R_vector);
        R_vector_mat[spin] = zeros(Int,num_R_vector,3)
    end

    current_line = 1
    while(current_line < length(hamiltonMatrices_fname_text))
        next_delimiter = search_next_delimiter(current_line, hamiltonMatrices_fname_text, "Real-space Hamiltonian for spin")

        # Detect spin
        tmp_str_list = split(hamiltonMatrices_fname_text[next_delimiter])
        spin = parse(Int, tmp_str_list[end])
        ##println("spin ",spin)

        # Detect R vector
        tmp_str_list = split(hamiltonMatrices_fname_text[next_delimiter+1])
        R_vector = zeros(Int,3)
        R_vector[1] = parse(Int,tmp_str_list[end-2])
        R_vector[2] = parse(Int,tmp_str_list[end-1])
        R_vector[3] = parse(Int,tmp_str_list[end-0])


        ##println("R_vector ", R_vector)

        current_R_vector_i = (R_vector_cnt[spin] +=1)
        R_vector_mat[spin][current_R_vector_i,:] = R_vector

        Hks_R[spin][current_R_vector_i] = zeros(ComplexF64, TotalOrbitalNum, TotalOrbitalNum)
        for i in 1:TotalOrbitalNum
            # Hk_R
            tmp_str_list_real = split(hamiltonMatrices_fname_text[next_delimiter+4+i]);
            tmp_str_list_imag = split(hamiltonMatrices_fname_text[next_delimiter+4+TotalOrbitalNum+3+ i])

            ham_real =  parse.(Float64, tmp_str_list_real[2:end])
            ham_imag =  parse.(Float64, tmp_str_list_imag[2:end])

            Hks_R[spin][current_R_vector_i][i,:] = ham_real + 1.0im*ham_imag

        end
        current_line = next_delimiter+4+TotalOrbitalNum+3+ TotalOrbitalNum
    end
    hamiltonMatrices_fname_text = []
    ###################################################
    ## Read Overlap matrix
    ###################################################
    overlapMatrices_fname_text  = readlines(overlapMatrices_fname);
    println(" ===================== ")
    println(" Total Overlap Matrix length: ",length(overlapMatrices_fname_text))
    for i in 1:9
        println(overlapMatrices_fname_text[i])
    end

    # num_R_vector_Overlap
    tmp_str_list = split(overlapMatrices_fname_text[1])
    num_R_vector_Overlap = parse(Int64, tmp_str_list[4])
    R_vector_Overlap_cnt = zeros(Int,SpinP_switch)
    println("num_R_vector_Overlap ",num_R_vector_Overlap)

    R_vector_Overlap_mat = Array{Array{Int,2}}(undef,SpinP_switch)
    OverlapS_R = Array{Array{Array{ComplexF64,2}}}(undef,SpinP_switch)

    for spin in 1:SpinP_switch
        OverlapS_R[spin]           = Array{Array{ComplexF64,2}}(undef,num_R_vector_Overlap);
        R_vector_Overlap_mat[spin] = zeros(Int,num_R_vector_Overlap,3)
    end

    current_line = 1
    next_delimiter = search_next_delimiter(current_line, overlapMatrices_fname_text,
            "Real-space Overlap for spin", maxscan=10)

    if -1 == next_delimiter
        # Lobster after 4.0: Overlap matrix of spin 1,2 is same duplicated infor is reduced
        next_delimiter = search_next_delimiter(current_line, overlapMatrices_fname_text,
            "Real-space Overlap", maxscan=10)

            while(current_line < length(overlapMatrices_fname_text))
                next_delimiter = search_next_delimiter(current_line, overlapMatrices_fname_text,
                    "Real-space Overlap", maxscan=10)
                
        
                # Detect spin
                #tmp_str_list = split(overlapMatrices_fname_text[next_delimiter])
                #spin = parse(Int, tmp_str_list[end])
                spin = 1
                ##println("spin ",spin)
        
                # Detect R vector
                tmp_str_list = split(overlapMatrices_fname_text[next_delimiter+1])
                R_vector = zeros(Int,3)
                R_vector[1] = parse(Int,tmp_str_list[end-2])
                R_vector[2] = parse(Int,tmp_str_list[end-1])
                R_vector[3] = parse(Int,tmp_str_list[end-0])
        
        
                ##println("R_vector ", R_vector)
        
                current_R_vector_i = (R_vector_Overlap_cnt[spin] +=1)
                R_vector_Overlap_mat[spin][current_R_vector_i,:] = R_vector
        
                OverlapS_R[spin][current_R_vector_i]  = zeros(ComplexF64, TotalOrbitalNum, TotalOrbitalNum)
                for i in 1:TotalOrbitalNum
                    # OverlapMatix
                    tmp_str_list_real = split(overlapMatrices_fname_text[next_delimiter+4+i]);
                    print(" DEBUG ",next_delimiter+4+i, " " ,tmp_str_list_real)
                    #tmp_str_list_imag = split(overlapMatrices_fname_text[next_delimiter+4+TotalOrbitalNum+3+ i])
        
                    overlap_real =  parse.(Float64, tmp_str_list_real[2:end])
                    #overlap_imag =  parse.(Float64, tmp_str_list_imag[2:end])
        
                    OverlapS_R[spin][current_R_vector_i][i,:] = overlap_real # + 1.0im*overlap_imag
                end
                current_line = next_delimiter+4+TotalOrbitalNum+3 # + TotalOrbitalNum
            end
            if 2==SpinP_switch
                OverlapS_R[2] = OverlapS_R[1]
            end
    else
        # Lobster before 4.0 : spin 1,2 is written
        while(current_line < length(overlapMatrices_fname_text))
            next_delimiter = search_next_delimiter(current_line, overlapMatrices_fname_text,
                "Real-space Overlap for spin", maxscan=10)
            
    
            # Detect spin
            tmp_str_list = split(overlapMatrices_fname_text[next_delimiter])
            spin = parse(Int, tmp_str_list[end])
            ##println("spin ",spin)
    
            # Detect R vector
            tmp_str_list = split(overlapMatrices_fname_text[next_delimiter+1])
            R_vector = zeros(Int,3)
            R_vector[1] = parse(Int,tmp_str_list[end-2])
            R_vector[2] = parse(Int,tmp_str_list[end-1])
            R_vector[3] = parse(Int,tmp_str_list[end-0])
    
    
            ##println("R_vector ", R_vector)
    
            current_R_vector_i = (R_vector_Overlap_cnt[spin] +=1)
            R_vector_Overlap_mat[spin][current_R_vector_i,:] = R_vector
    
            OverlapS_R[spin][current_R_vector_i]  = zeros(ComplexF64, TotalOrbitalNum, TotalOrbitalNum)
            for i in 1:TotalOrbitalNum
                # OverlapMatix
                tmp_str_list_real = split(overlapMatrices_fname_text[next_delimiter+4+i]);
                tmp_str_list_imag = split(overlapMatrices_fname_text[next_delimiter+4+TotalOrbitalNum+3+ i])
    
                overlap_real =  parse.(Float64, tmp_str_list_real[2:end])
                overlap_imag =  parse.(Float64, tmp_str_list_imag[2:end])
    
                OverlapS_R[spin][current_R_vector_i][i,:] = overlap_real + 1.0im*overlap_imag
            end
            current_line = next_delimiter+4+TotalOrbitalNum+3+ TotalOrbitalNum
        end
    
    end
    

    overlapMatrices_fname_text = []
    ###################################################
    #ChemP = 0.0
    Wannier_info = Lobsterdatatype(atomnum,  SpinP_switch,Total_NumOrbs, tv, rv,
    Gxyz, R_vector_mat, Hks_R, R_vector_Overlap_mat, OverlapS_R, ChemP,
    spin_type, 300.0)
    return Wannier_info;
end
