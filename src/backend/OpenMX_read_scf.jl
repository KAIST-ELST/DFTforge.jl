
#using ArrayViews

#scf_name = "/home/users1/bluehope/work/Jx/NiO/4cell/AF2/Ni6.0_O5.0-s2p2d2f1-5.5_0.0-4.180-k10/nio.scfout";
const Hatree2eV = 27.2114;
macro init_OLPmat(mat)
:(
    for ct_AN=1:atomnum
        TNO1 = Total_NumOrbs[ct_AN];
        $mat[ct_AN] =  Array(Array{Array{Float64,},},FNAN[ct_AN]+1);

        for h_AN=1:FNAN[ct_AN]+1
            $mat[ct_AN][h_AN] =  Array(Array{Float64,},TNO1);
            Gh_AN = natn[ct_AN][h_AN];
            TNO2 = Total_NumOrbs[Gh_AN];
            for i=1:TNO1
                $mat[ct_AN][h_AN][i] = Array(Float64,TNO2);
            end
        end
    end)
end
macro read_OLPmat(mat)
    :(
        for ct_AN=1:atomnum
            TNO1 = Total_NumOrbs[ct_AN];

            #read(f,Float64)
            for h_AN=1:FNAN[ct_AN]+1
                Gh_AN = natn[ct_AN][h_AN];
                TNO2 = Total_NumOrbs[Gh_AN];
                for i=1:TNO1
                    $mat[ct_AN][h_AN][i] = read(f,Float64,TNO2)
                end
            end
        end

    )
    #:(println($mat[1][1][1]);)
end

macro init_Hamil(mat,spin,atomnum)
    :(
    for ct_AN=1:$atomnum
        TNO1 = Total_NumOrbs[ct_AN];
        $mat[$spin][ct_AN] =
        Array(Array{Array{Float64,},},FNAN[ct_AN]+1);

        for h_AN=1:FNAN[ct_AN]+1
            $mat[$spin][ct_AN][h_AN] =
                    Array(Array{Float64,},TNO1);

            Gh_AN = natn[ct_AN][h_AN];
            TNO2 = Total_NumOrbs[Gh_AN];

            for i=1:TNO1
                $mat[$spin][ct_AN][h_AN][i] =
                Array(Float64,TNO2);
            end
        end
    end
    )
end

macro read_Hamil(mat,spin,atomnum)
    :(
        for ct_AN=1:$atomnum
        TNO1 = Total_NumOrbs[ct_AN];

        #read(f,Float64)
        for h_AN=1:FNAN[ct_AN]+1
            Gh_AN = natn[ct_AN][h_AN];
            TNO2 = Total_NumOrbs[Gh_AN];
            for i=1:TNO1
                $mat[$spin][ct_AN][h_AN][i] = read(f,Float64,TNO2)*Hatree2eV
                #for j =1:TNO2
                #    Hks[spin][ct_AN][h_AN][i][j] = read(f,Float64)
                #end
            end
        end
    end
    )
end


H_type = Array{Array{Array{Array{Array{Float64}}}}};
Overlap_type =  Array{Array{Array{Array{Float64}}}};



immutable Openmxscf
    atomnum::Int32
    SpinP_switch::Int32
    Catomnum::Int32
    Latomnum::Int32
    Ratomnum::Int32
    TCpyCell::Int32

    atv::Array{Float64,2}
    atv_ijk::Array{Int32,2}
    Total_NumOrbs::Array{Int32,1}
    FNAN::Array{Int32,1}

    natn::Array{Array{Int32},1}
    ncn::Array{Array{Int32},1}

    tv::Array{Float64,2}
    rv::Array{Float64,2}

    Gxyz::Array{Float64,2}

    Hks::H_type
    iHks::H_type

    OLP::Overlap_type
    OLPpox::Overlap_type
    OLPpoy::Overlap_type
    OLPpoz::Overlap_type
    DM::H_type

    Solver::Int32

    ChemP::Float64
    E_Temp::Float64
    #dipole_moment_core = zeros(3,)
    #dipole_moment_core[1] = d_vec[3];
    #dipole_moment_core[2] = d_vec[4];
    #dipole_moment_core[3] = d_vec[5];
    #dipole_moment_background = zeros(3,)
    #dipole_moment_background[1] = d_vec[6];
    #dipole_moment_background[2] = d_vec[7];
    #dipole_moment_background[3] = d_vec[8];
    #Valence_Electrons = d_vec[9];
    #Total_SpinS = d_vec[10];
end

function read_scf(scf_name::AbstractString)
    # atomnum

    f = open(scf_name,"r")
    atomnum =  read(f,Int32)
    SpinP_switch =  read(f,Int32)
    Catomnum =  read(f,Int32)
    Latomnum =  read(f,Int32)
    Ratomnum =  read(f,Int32)
    TCpyCell =  read(f,Int32)

    # atv
    atv = zeros(Float64,TCpyCell+1,3)
    for Rn=1:TCpyCell+1
        read(f,Float64)
        for ii=1:3
            atv[Rn,ii] = read(f,Float64) / DFTcommon.ang2bohr;
        end
    end


    # atv_ijk
    atv_ijk = zeros(Int32,TCpyCell+1,3)
    for Rn=1:TCpyCell+1
        read(f,Int32)
        for ii=1:3
            atv_ijk[Rn,ii] = read(f,Int32)
        end
    end
    #Total_NumOrbs
    Total_NumOrbs = zeros(Int32,atomnum,)
    for ii=1:atomnum
        Total_NumOrbs[ii] = read(f,Int32)
    end
    #read(f,Int32)
    #println(Total_NumOrbs)
    #FNAN
    FNAN = zeros(Int32,atomnum+1,)
    for ii=1:atomnum
        FNAN[ii] = read(f,Int32)
    end
    #println(FNAN)


    #natn[atomnum+1][FNAN[ct_AN]+1];
    #ncn[atomnum+1][FNAN[ct_AN]+1];
    natn =  Array{Array{Int32,}}(atomnum)
    ncn =  Array{Array{Int32,}}(atomnum)
    for ii=1:atomnum
        natn[ii] = zeros(Int32,FNAN[ii]+1)
        natn[ii] = read(f,Int32,FNAN[ii]+1)
        #natn[ii] = zeros(Int32,FNAN[ii])
        #read(f,Int32);
        #natn[ii] = read(f,Int32,FNAN[ii])
    end
    #println(natn)
    for ii=1:atomnum
        ncn[ii] = zeros(Int32,FNAN[ii]+1)
        ncn[ii] = read(f,Int32,FNAN[ii]+1)
        #ncn[ii] = zeros(Int32,FNAN[ii])
        #read(f,Int32);
        #ncn[ii] = read(f,Int32,FNAN[ii])

    end

    # natn[][]:
    # grobal index of neighboring atoms of an atom ct_AN
    # tv[4][4]: -> tv[3][3]
    # unit cell vectors in Bohr
    tv = zeros(Float64,3,3)
    for ii=1:3
        read(f,Float64)
        for jj =1:3
            tv[ii,jj] = read(f,Float64)  / DFTcommon.ang2bohr
        end
    end

    # unit cell vectors in Bohr

    rv = zeros(Float64,3,3)
    for ii=1:3
        read(f,Float64)
        for jj =1:3
            rv[ii,jj] = read(f,Float64)  * DFTcommon.ang2bohr
        end
    end

    if(sum((2*pi*inv(tv') .- rv)[:]) > 10.0^-10)
        println("SCF reading error")
    end

    #  Gxyz[][1-3]:
    #  atomic coordinates in Bohr
    Gxyz = zeros(Float64,atomnum,3)
    for ii=1:atomnum
        read(f,Float64)
        for jj =1:3
            Gxyz[ii,jj] = read(f,Float64)/DFTcommon.ang2bohr
        end
    end


    #     allocation of arrays:
    #
    #     Kohn-Sham Hamiltonian
    #
    #      dooble Hks[SpinP_switch+1]
    #                [atomnum+1]
    #                [FNAN[ct_AN]+1]
    #                [Total_NumOrbs[ct_AN]]
    #                [Total_NumOrbs[h_AN]];
    #
    #     Overlap matrix
    #
    #      dooble OLP[atomnum+1]
    #                [FNAN[ct_AN]+1]
    #                [Total_NumOrbs[ct_AN]]
    #                [Total_NumOrbs[h_AN]];
    #
    #     Overlap matrix with position operator x, y, z
    #
    #      dooble OLPpox,y,z
    #                  [atomnum+1]
    #                  [FNAN[ct_AN]+1]
    #                  [Total_NumOrbs[ct_AN]]
    #                  [Total_NumOrbs[h_AN]];
    #
    #     Density matrix
    #
    #      dooble DM[SpinP_switch+1]
    #               [atomnum+1]
    #               [FNAN[ct_AN]+1]
    #               [Total_NumOrbs[ct_AN]]
    #                [Total_NumOrbs[h_AN]];
    Hks = Array{Array{Array{Array{Array{Float64}}}}}(SpinP_switch+1)
    iHks = Array{Array{Array{Array{Array{Float64}}}}}(3)

    OLP = Array{Array{Array{Array{Float64}}}}(atomnum)
    OLPpox = Array{Array{Array{Array{Float64}}}}(atomnum);
    OLPpoy = Array{Array{Array{Array{Float64}}}}(atomnum);
    OLPpoz = Array{Array{Array{Array{Float64}}}}(atomnum);
    DM = Array{Array{Array{Array{Array{Float64}}}}}(SpinP_switch+1)

    for spin = 1:SpinP_switch+1
        Hks[spin] = Array{Array{Array{Array{Float64}}}}(atomnum)
        @init_Hamil(Hks,spin,atomnum)
    end
    for spin = 1:3
        iHks[spin] = Array{Array{Array{Array{Float64}}}}(atomnum)
        @init_Hamil(iHks,spin,atomnum)
    end

    # Hamiltonian matrix
    for spin = 1:SpinP_switch+1
        @read_Hamil(Hks,spin,atomnum)
    end
    if 3 == SpinP_switch # non-collinear
        for spin = 1:3#SpinP_switch+1
            @read_Hamil(iHks,spin,atomnum)
        end
    end

    # OLP OLPx OLPy OLPz
    @init_OLPmat(OLP)
    @init_OLPmat(OLPpox)
    @init_OLPmat(OLPpoy)
    @init_OLPmat(OLPpoz)

    @read_OLPmat(OLP)
    @read_OLPmat(OLPpox)
    @read_OLPmat(OLPpoy)
    @read_OLPmat(OLPpoz)

    #for spin = 1:3
    #    iHks[spin] = Array(Array{Array{Array{Float64,},},},atomnum)
    #    @init_Hamil(iHks,spin)
    #end

    for spin = 1:SpinP_switch+1
        DM[spin] = Array{Array{Array{Array{Float64}}}}(atomnum)
        @init_Hamil(DM,spin,atomnum)
    end
    for spin = 1:SpinP_switch+1
        @read_Hamil(DM,spin,atomnum)
    end
    #############
    # Solver
    #############
    Solver = read(f,Int32)
    d_vec = read(f,Float64,10)
    ChemP  = d_vec[1] * Hatree2eV;
    E_Temp = d_vec[2];
    dipole_moment_core = zeros(3,)
    dipole_moment_core[1] = d_vec[3];
    dipole_moment_core[2] = d_vec[4];
    dipole_moment_core[3] = d_vec[5];
    dipole_moment_background = zeros(3,)
    dipole_moment_background[1] = d_vec[6];
    dipole_moment_background[2] = d_vec[7];
    dipole_moment_background[3] = d_vec[8];
    Valence_Electrons = d_vec[9];
    Total_SpinS = d_vec[10];


    #
    ########################
    # input file
    ########################
    num_lines = read(f,Int32)

    #reading input file is not implimented
    close(f)

    result = Openmxscf(    atomnum,
    SpinP_switch,
    Catomnum,
    Latomnum,
    Ratomnum,
    TCpyCell,

    atv,
    atv_ijk,
    Total_NumOrbs,
    FNAN,

    natn,
    ncn,

    tv,
    rv,

    Gxyz,

    Hks,
    iHks,

    OLP,
    OLPpox,
    OLPpoy,
    OLPpoz,
    DM,

    Solver,

    ChemP,
    E_Temp)
    #println("scf_name: ",scf_name)
    return result
end
