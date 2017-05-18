using Rotations # Pkg.add("Rotations")
export rot_matrixByZXaxisbase,rot_D_orbital
export orbital_rot_d_type

type orbital_rot_d_type
  atomnum::Int
  Z_vect::Vector{Float64}
  X_vect::Vector{Float64}
  d_orbital_list::Array{Array{Int}} #index of d-orbitals z^2, x2y2, xy, xz, yz

  function orbital_rot_d_type(atomnum,Z_vect,X_vect,d_orbital_list)
    assert(3==length(Z_vect))
    assert(3==length(X_vect))
    assert(5==length(d_orbital_list))
    new(atomnum,Z_vect,X_vect,d_orbital_list)
  end
end

function rot_matrixByZXaxisbase(Z_vect::Vector{Float64},X_vect::Vector{Float64})
  #Z_vect = normalize([0.0 , 0.0, 2.0])
  #X_vect = normalize([1.0 , 1.0, 0.0])
  Y_vect = normalize(cross(Z_vect,X_vect))
  R = zeros(3,3)

  R[1,:] = X_vect
  R[2,:] = Y_vect
  R[3,:] = Z_vect
  return R;
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
  assert( abs(norm(R)-1.0) < 10.0^-6.0);

  h = 2.0
  s3 = sqrt(3.0)
  DR = eye(5,5)

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
  assert( abs(norm(DR)-1.0) < 10.0^-6.0);

  return DR;
end
