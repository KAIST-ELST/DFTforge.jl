###############################################################################
# Hongkee Yoon Hongkeeyoon@kaist.ac.kr
# 2019.05
# https://kaist-elst.github.io/DFTforge.jl/
###############################################################################

#using ArrayViews

const Hatree2eV = 27.2114;
function init_OLPmat!(mat, atomnum, Total_NumOrbs, FNAN, natn)
    for ct_AN=1:atomnum
        TNO1 = Total_NumOrbs[ct_AN];
        mat[ct_AN] =  Array{Array{Array{Float64}}}(undef,FNAN[ct_AN]+1);

        for h_AN=1:FNAN[ct_AN]+1
            mat[ct_AN][h_AN] =  Array{Array{Float64}}(undef,TNO1);
            Gh_AN = natn[ct_AN][h_AN];
            TNO2 = Total_NumOrbs[Gh_AN];
            for i=1:TNO1
                mat[ct_AN][h_AN][i] = Array{Float64}(undef,TNO2);
            end
        end
    end
end

function read_OLPmat!(f, mat, atomnum, Total_NumOrbs, FNAN, natn)
  for ct_AN=1:atomnum
      TNO1 = Total_NumOrbs[ct_AN];
      for h_AN=1:FNAN[ct_AN]+1
          Gh_AN = natn[ct_AN][h_AN];
          TNO2 = Total_NumOrbs[Gh_AN];
          for i=1:TNO1
              #mat[ct_AN][h_AN][i] = read(f,Float64,TNO2)
              read!(f,mat[ct_AN][h_AN][i])
          end
      end
  end
end

function  write_OLPmat(f, mat, scalefactor, scf_openmx)

    natn = scf_openmx.natn;
    Total_NumOrbs = scf_openmx.Total_NumOrbs
    FNAN = scf_openmx.FNAN
    atomnum = scf_openmx.atomnum

    for ct_AN=1:atomnum
        TNO1 = Total_NumOrbs[ct_AN];
        for h_AN=1:FNAN[ct_AN]+1
            Gh_AN = natn[ct_AN][h_AN];
            TNO2 = Total_NumOrbs[Gh_AN];
            for i=1:TNO1
                #$mat[ct_AN][h_AN][i] = read(f,Float64,TNO2)
                for j = 1:TNO2
                  write(f,convert(Float64,mat[ct_AN][h_AN][i][j]))
                end
            end
        end
    end


    #:(println($mat[1][1][1]);)
end

function init_Hamil!(mat, spin, atomnum, Total_NumOrbs, FNAN, natn)
    for ct_AN=1:atomnum
        TNO1 = Total_NumOrbs[ct_AN];
        mat[spin][ct_AN] =
        Array{Array{Array{Float64}}}(undef,FNAN[ct_AN]+1);

        for h_AN=1:FNAN[ct_AN]+1
            mat[spin][ct_AN][h_AN] =
                    Array{Array{Float64}}(undef,TNO1);

            Gh_AN = natn[ct_AN][h_AN];
            TNO2 = Total_NumOrbs[Gh_AN];

            for i=1:TNO1
                mat[spin][ct_AN][h_AN][i] =
                Array{Float64}(undef,TNO2);
            end
        end
    end

end

function read_Hamil(f, mat, spin, atomnum, Total_NumOrbs, FNAN, natn, scalefactor)
    for ct_AN=1:atomnum
      TNO1 = Total_NumOrbs[ct_AN];
        for h_AN=1:FNAN[ct_AN]+1
            Gh_AN = natn[ct_AN][h_AN];
            TNO2 = Total_NumOrbs[Gh_AN];
            for i=1:TNO1
                #mat[spin][ct_AN][h_AN][i] = read(f,Float64,TNO2)*scalefactor
                read!(f,mat[spin][ct_AN][h_AN][i])
                mat[spin][ct_AN][h_AN][i] *= scalefactor
                #for j =1:TNO2
                #    Hks[spin][ct_AN][h_AN][i][j] = read(f,Float64)
                #end
            end
        end
    end
end

function write_Hamil(f, spin , mat, scalefactor, scf_openmx)

  natn = scf_openmx.natn;
  Total_NumOrbs = scf_openmx.Total_NumOrbs
  FNAN = scf_openmx.FNAN
  atomnum = scf_openmx.atomnum

  for ct_AN=1:atomnum
    TNO1 = Total_NumOrbs[ct_AN];
    for h_AN=1:FNAN[ct_AN]+1
      Gh_AN = natn[ct_AN][h_AN];
      TNO2 = Total_NumOrbs[Gh_AN];
      for i=1:TNO1
          #$mat[$spin][ct_AN][h_AN][i] = read(f,Float64,TNO2)*Hatree2eV
          for j = 1:TNO2
            write(f,convert(Float64,mat[spin][ct_AN][h_AN][i][j]*scalefactor))
          end
          #for j =1:TNO2
          #    Hks[spin][ct_AN][h_AN][i][j] = read(f,Float64)
          #end
      end
    end
  end
end


H_type = Array{Array{Array{Array{Array{Float64}}}}};
Overlap_type =  Array{Array{Array{Array{Float64}}}};



struct Openmxscf
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
    natn = Array{Array{Int32}}(undef,atomnum)
    ncn =  Array{Array{Int32}}(undef,atomnum)
    for ii=1:atomnum
        natn[ii] = zeros(Int32,FNAN[ii]+1)
        #natn[ii] = read(f,Int32,FNAN[ii]+1)
        read!(f,natn[ii])

        #natn[ii] = zeros(Int32,FNAN[ii])
        #read(f,Int32);
        #natn[ii] = read(f,Int32,FNAN[ii])
    end
    #println(natn)
    for ii=1:atomnum
        ncn[ii] = zeros(Int32,FNAN[ii]+1)
        #ncn[ii] = read(f,Int32,FNAN[ii]+1)
        read!(f,ncn[ii])
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
    #println(tv,rv)
    if (sum((2*pi*inv( copy(transpose(tv)) ) .- rv)[:]) > 10.0^-10)
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
    Hks = Array{Array{Array{Array{Array{Float64}}}}}(undef,SpinP_switch+1)
    iHks = Array{Array{Array{Array{Array{Float64}}}}}(undef,3)

    OLP = Array{Array{Array{Array{Float64}}}}(undef,atomnum)
    OLPpox = Array{Array{Array{Array{Float64}}}}(undef,atomnum);
    OLPpoy = Array{Array{Array{Array{Float64}}}}(undef,atomnum);
    OLPpoz = Array{Array{Array{Array{Float64}}}}(undef,atomnum);
    DM = Array{Array{Array{Array{Array{Float64}}}}}(undef,SpinP_switch+1)

    for spin = 1:SpinP_switch+1
        Hks[spin] = Array{Array{Array{Array{Float64}}}}(undef,atomnum)
        init_Hamil!(Hks, spin, atomnum, Total_NumOrbs, FNAN, natn)
    end
    for spin = 1:3
        iHks[spin] = Array{Array{Array{Array{Float64}}}}(undef,atomnum)
        init_Hamil!(iHks, spin, atomnum, Total_NumOrbs, FNAN, natn)
    end

    # Hamiltonian matrix
    for spin = 1:SpinP_switch+1
        read_Hamil(f, Hks, spin,atomnum, Total_NumOrbs, FNAN, natn, Hatree2eV)
    end
    if 3 == SpinP_switch # non-collinear
        for spin = 1:3#SpinP_switch+1
            read_Hamil(f, iHks, spin, atomnum, Total_NumOrbs, FNAN, natn, Hatree2eV)
        end
    end

    # OLP OLPx OLPy OLPz
    init_OLPmat!(OLP, atomnum, Total_NumOrbs, FNAN, natn)
    init_OLPmat!(OLPpox, atomnum, Total_NumOrbs, FNAN, natn)
    init_OLPmat!(OLPpoy, atomnum, Total_NumOrbs, FNAN, natn)
    init_OLPmat!(OLPpoz, atomnum, Total_NumOrbs, FNAN, natn)

    read_OLPmat!(f, OLP, atomnum, Total_NumOrbs, FNAN, natn)
    read_OLPmat!(f, OLPpox, atomnum, Total_NumOrbs, FNAN, natn)
    read_OLPmat!(f, OLPpoy, atomnum, Total_NumOrbs, FNAN, natn)
    read_OLPmat!(f, OLPpoz, atomnum, Total_NumOrbs, FNAN, natn)

    #for spin = 1:3
    #    iHks[spin] = Array(Array{Array{Array{Float64,},},},atomnum)
    #    init_Hamil!(iHks,spin)
    #end

    for spin = 1:SpinP_switch+1
        DM[spin] = Array{Array{Array{Array{Float64}}}}(undef,atomnum)
        init_Hamil!(DM, spin, atomnum, Total_NumOrbs,FNAN, natn)
    end
    for spin = 1:SpinP_switch+1
        read_Hamil(f, DM, spin, atomnum, Total_NumOrbs, FNAN, natn, 1.0)
    end
    #############
    # Solver
    #############
    Solver = read(f,Int32)
    #d_vec = read(f,Float64,10)
    d_vec = Array{Float64}(undef,10);
    read!(f,d_vec)
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

function write_scf(scf_openmx::Openmxscf,scf_name::AbstractString)
  # atomnum

  f = open(scf_name,"w")
  #atomnum =  read(f,Int32)
  write(f,convert(Int32, scf_openmx.atomnum));
  #SpinP_switch =  read(f,Int32)
  write(f,convert(Int32, scf_openmx.SpinP_switch));
  #Catomnum =  read(f,Int32)
  write(f,convert(Int32, scf_openmx.Catomnum));
  #Latomnum =  read(f,Int32)
  write(f,convert(Int32, scf_openmx.Latomnum));
  #Ratomnum =  read(f,Int32)
  write(f,convert(Int32, scf_openmx.Ratomnum))
  #TCpyCell =  read(f,Int32)
  write(f,convert(Int32, scf_openmx.TCpyCell))

  # atv
  #atv = zeros(Float64,TCpyCell+1,3)
  for Rn=1:scf_openmx.TCpyCell+1
      #read(f,Float64)
      write(f,convert(Float64, 0.0));
      for ii=1:3
          #atv[Rn,ii] = read(f,Float64) / DFTcommon.ang2bohr;
          write(f,convert(Float64, scf_openmx.atv[Rn,ii] * DFTcommon.ang2bohr));
      end
  end


  # atv_ijk
  #atv_ijk = zeros(Int32,TCpyCell+1,3)
  for Rn=1:scf_openmx.TCpyCell+1
      #read(f,Int32)
      write(f,convert(Int32,0))
      for ii=1:3
          #atv_ijk[Rn,ii] = read(f,Int32)
          write(f,convert(Int32,scf_openmx.atv_ijk[Rn,ii]))
      end
  end

  #Total_NumOrbs
  #Total_NumOrbs = zeros(Int32,atomnum,)
  for ii=1:scf_openmx.atomnum
      #Total_NumOrbs[ii] = read(f,Int32)
      write(f,convert(Int32,scf_openmx.Total_NumOrbs[ii]))
  end
  #read(f,Int32)
  #println(Total_NumOrbs)

  #FNAN
  #FNAN = zeros(Int32,atomnum+1,)
  for ii=1:scf_openmx.atomnum
      #FNAN[ii] = read(f,Int32)
      write(f,convert(Int32, scf_openmx.FNAN[ii]))
  end
  #println(FNAN)


  #natn[atomnum+1][FNAN[ct_AN]+1];
  #ncn[atomnum+1][FNAN[ct_AN]+1];

  #natn =  Array{Array{Int32,}}(atomnum)
  #ncn =  Array{Array{Int32,}}(atomnum)
  for ii=1:scf_openmx.atomnum
      #natn[ii] = zeros(Int32,FNAN[ii]+1)
      #natn[ii] = read(f,Int32,FNAN[ii]+1)
      for j=1:scf_openmx.FNAN[ii]+1
        write(f,convert(Int32, scf_openmx.natn[ii][j]));
      end
  end
  #println(natn)
  for ii=1:scf_openmx.atomnum
      #ncn[ii] = zeros(Int32,FNAN[ii]+1)
      #ncn[ii] = read(f,Int32,FNAN[ii]+1)
      for j=1:scf_openmx.FNAN[ii]+1
        write(f,convert(Int32, scf_openmx.ncn[ii][j]));
      end
  end

  # natn[][]:
  # grobal index of neighboring atoms of an atom ct_AN
  # tv[4][4]: -> tv[3][3]
  # unit cell vectors in Bohr
  #tv = zeros(Float64,3,3)
  for ii=1:3
      #read(f,Float64)
      write(f, convert(Float64, 0.0)) # null write
      for jj =1:3
          #tv[ii,jj] = read(f,Float64)  / DFTcommon.ang2bohr
          write(f,convert(Float64, scf_openmx.tv[ii,jj]* DFTcommon.ang2bohr));
      end
  end

  # unit cell vectors in Bohr

  #rv = zeros(Float64,3,3)
  for ii=1:3
      #read(f,Float64)
      write(f, convert(Float64, 0.0)) # null write
      for jj =1:3
          #rv[ii,jj] = read(f,Float64)  * DFTcommon.ang2bohr
          write(f, convert(Float64, scf_openmx.rv[ii,jj] / DFTcommon.ang2bohr))
      end
  end

  if (sum((2*pi*inv(scf_openmx.tv') .- scf_openmx.rv)[:]) > 10.0^-10)
      println("SCF writting error")
  end

  #  Gxyz[][1-3]:
  #  atomic coordinates in Bohr
  #Gxyz = zeros(Float64,atomnum,3)
  for ii=1:scf_openmx.atomnum
      #read(f,Float64)
      write(f, convert(Float64, 0.0)) # null write
      for jj =1:3
          #Gxyz[ii,jj] = read(f,Float64)/DFTcommon.ang2bohr
          write(f, convert(Float64, scf_openmx.Gxyz[ii,jj] * DFTcommon.ang2bohr))
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

  #=
  Hks = Array{Array{Array{Array{Array{Float64}}}}}(SpinP_switch+1)
  iHks = Array{Array{Array{Array{Array{Float64}}}}}(3)

  OLP = Array{Array{Array{Array{Float64}}}}(atomnum)
  OLPpox = Array{Array{Array{Array{Float64}}}}(atomnum);
  OLPpoy = Array{Array{Array{Array{Float64}}}}(atomnum);
  OLPpoz = Array{Array{Array{Array{Float64}}}}(atomnum);
  DM = Array{Array{Array{Array{Array{Float64}}}}}(SpinP_switch+1)
 =#

#=
  for spin = 1:SpinP_switch+1
      Hks[spin] = Array{Array{Array{Array{Float64}}}}(atomnum)
      init_Hamil!(Hks,spin,atomnum)
  end
=#
#=
  for spin = 1:3
      iHks[spin] = Array{Array{Array{Array{Float64}}}}(atomnum)
      #init_Hamil!(iHks,spin,atomnum)
  end
=#

  # Hamiltonian matrix
  for spin = 1:scf_openmx.SpinP_switch+1
      #@read_Hamil(Hks,spin,atomnum)
      write_Hamil(f, spin, scf_openmx.Hks, 1.0/Hatree2eV, scf_openmx)
  end
  if 3 == scf_openmx.SpinP_switch # non-collinear
      for spin = 1:3#SpinP_switch+1
          #@read_Hamil(iHks,spin,atomnum)
          write_Hamil(f, spin, scf_openmx.iHks, 1.0/Hatree2eV, scf_openmx)
      end
  end

  # OLP OLPx OLPy OLPz

  #=
  @init_OLPmat(OLP)
  @init_OLPmat(OLPpox)
  @init_OLPmat(OLPpoy)
  @init_OLPmat(OLPpoz)
=#
#=
  @read_OLPmat(OLP)
  @read_OLPmat(OLPpox)
  @read_OLPmat(OLPpoy)
  @read_OLPmat(OLPpoz)
=#
  write_OLPmat(f, scf_openmx.OLP, 1.0, scf_openmx)
  write_OLPmat(f, scf_openmx.OLPpox, 1.0, scf_openmx)
  write_OLPmat(f, scf_openmx.OLPpoy, 1.0, scf_openmx)
  write_OLPmat(f, scf_openmx.OLPpoz, 1.0, scf_openmx)
  #for spin = 1:3
  #    iHks[spin] = Array(Array{Array{Array{Float64,},},},atomnum)
  #    init_Hamil!(iHks,spin)
  #end
#=
  for spin = 1:SpinP_switch+1
      DM[spin] = Array{Array{Array{Array{Float64}}}}(atomnum)
      init_Hamil!(DM,spin,atomnum)
  end
=#
  for spin = 1:scf_openmx.SpinP_switch+1
      #@read_Hamil(DM,spin,atomnum)
      write_Hamil(f, spin, scf_openmx.DM, 1.0, scf_openmx)
  end
  #############
  # Solver
  #############
  #Solver = read(f,Int32)
  write(f, convert(Int32, scf_openmx.Solver));
  #d_vec = read(f,Float64,10)

  #ChemP  = d_vec[1] * Hatree2eV;
  write(f, convert(Float64, scf_openmx.ChemP/Hatree2eV));
  #E_Temp = d_vec[2];
  write(f, convert(Float64, scf_openmx.E_Temp));
  #=
  dipole_moment_core = zeros(3,)
  dipole_moment_core[1] = d_vec[3];
  dipole_moment_core[2] = d_vec[4];
  dipole_moment_core[3] = d_vec[5];
  =#
  write(f, convert(Float64, 0.0));
  write(f, convert(Float64, 0.0));
  write(f, convert(Float64, 0.0));
  #=
  dipole_moment_background = zeros(3,)
  dipole_moment_background[1] = d_vec[6];
  dipole_moment_background[2] = d_vec[7];
  dipole_moment_background[3] = d_vec[8];
  =#

  write(f, convert(Float64, 0.0));
  write(f, convert(Float64, 0.0));
  write(f, convert(Float64, 0.0));
  #Valence_Electrons = d_vec[9];
  write(f, convert(Float64, 0.0));
  #Total_SpinS = d_vec[10];
  write(f, convert(Float64, 0.0));


  #
  ########################
  # input file
  ########################
  #num_lines = read(f,Int32)
  write(f, convert(Int32, 0));
  #reading input file is not implimented
  close(f)

#=
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
  =#
end
