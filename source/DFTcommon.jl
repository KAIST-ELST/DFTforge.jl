module DFTcommon
export Kpoint_eigenstate,Complex_my,Float_my,k_point_Tuple,k_point_int_Tuple
export k_point_int2float,k_point_float2int,kPoint2BrillouinZone_Tuple
export pwork

export eigfact_hermitian,check_eigmat

const Hartree2cm = 2.194746*100000.0;
const cm2meV = 0.1240;
const Hartree2meV = Hartree2cm*cm2meV;
const kB=0.000003166813628;  # Boltzman constant (Hatree/K)
const k_point_precision = 10.0^6;


Float_my =  Float64
Complex_my = Complex128

k_point_Tuple = Tuple{Float64,Float64,Float64};
k_point_int_Tuple = Tuple{Int64,Int64,Int64};

type Kpoint_eigenstate
    Eigenstate::Array{Complex_my,2};
    Eigenvalues::Array{Float_my,1};
    k_point_int::k_point_Tuple;
end


function k_point_int2float(k_point_int::k_point_int_Tuple)
    k_point = k_point_Tuple;
    k_point = (k_point_int[1]/k_point_precision,k_point_int[2]/k_point_precision,
    k_point_int[3]/k_point_precision);
    return k_point;
end
function k_point_float2int(k_point::k_point_Tuple)
    #k_point_int = k_point_int_Tuple;
    k_point_int = (round(Int64,k_point[1]*k_point_precision),round(Int64,k_point[2]*k_point_precision),
    round(Int64,k_point[3]*k_point_precision) );
    return k_point_int;
end
function kPoint2BrillouinZone(k_point::Array{Float64,1})
    # from -0.5 0.5
    # 0.5 is treated as -0.5
    return (rem(k_point+2.5,1)-0.5);
end

function kPoint2BrillouinZone_Tuple(k_point::k_point_Tuple)
    # from -0.5 0.5
    # 0.5 is treated as -0.5
    k_point_array =    kPoint2BrillouinZone([k_point[1],k_point[2],k_point[3]]);
    return (k_point_array[1],k_point_array[2],k_point_array[3])
end
function kPoint2BrillouinZone_int_Tuple(k_point_int::k_point_int_Tuple)
    # from -0.5 0.5
    # 0.5 is treated as -0.5
    k_point_int_array = [k_point_int[1],k_point_int[2],k_point_int[3]]
    k_point_int_array =  (rem(k_point_int_array+(2.5*k_point_precision),k_point_precision)-0.5*k_point_precision);
    k_point_int_array = round(Int64,k_point_int_array);
    return (k_point_int_array[1],k_point_int_array[2],k_point_int_array[3])
end


@enum orbital_mask_enum nomask=0 unmask=1 mask=2
type orbital_mask_input_Type
    orbital_mask1::Array{Int64,1}
    orbital_mask2::Array{Int64,1}
    #atom12::atom12_Tuple
    orbital_mask_on::Bool
end

function pwork(f,args)
    worker_list = workers()
    n = length(worker_list)
    results = Array{Any}(n)
    @sync begin
        for (worer_idx,worker_id) in enumerate(worker_list)
            @async begin
            #@everywhere init_all(atom1,atom2)
                results[worer_idx] = remotecall_fetch(f,worker_id,args)
            end
        end
    end
end
################################################################################
# Math related function
################################################################################
function eigfact_hermitian(EigVect::Array{Complex_my,2},
    EigVal::Array{Float_my,1})
    for ii = 1:size(EigVect)[1]
        EigVect[ii,ii] = real(EigVect[ii,ii]);
    end
    #EigVect_h = Hermitian(EigVect,:L);
    EigVect_h = Hermitian(EigVect);

    temp = eigfact(EigVect_h);
    p = sortperm(temp.values);
    #EigVal[:] = real(temp.values[p])[:];

    EigVal[:] = (temp.values[p])[:];
    EigVect[:] = temp.vectors[:,p][:];
    #ccall((:EigenBand_lapack3,"./lib_jx"),Void,(Ptr{Complex128},Ptr{Float64}, Int32,)
    #,EigVect,EigVect,size(EigVect)[1])
end
function check_eigmat(A1,eig_mat,eig_val) #A1 orign  eig_mat # eig_val
    n = size(A1)[1];
    check_eig = zeros(Complex_my,n,)
    for ii = 1:n
        check_eig[ii] =  sum((A1*eig_mat[:,ii])[:] - (eig_mat[:,ii]*eig_val[ii])[:] )
    end
        #A1 = (A1+conj(A1'))/2;
    diff_A = ((A1+conj(A1'))/2 - real(A1)) ;

    #println("")
    if(maximum(abs.(check_eig)) > 10.0^-3)
      println(" sum ",sum(check_eig)," max ",maximum(abs(check_eig))," Diff ",sum(diff_A),
      " Diff realmax ",maximum(real(diff_A)),
      " Diff imgmax ",maximum(imag(diff_A)));
        return false
    end
    #assert( maximum(abs(check_eig)) < 10.0^-3);
    return true;
end


end
