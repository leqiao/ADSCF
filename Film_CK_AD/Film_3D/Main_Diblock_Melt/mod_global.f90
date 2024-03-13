!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    GLOBAL VARIABLES (SYSTEM DIMENSIONS ETC.)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module Global
!
  implicit none
  double precision, parameter :: pi = 3.141592653589793d0
  double precision, parameter :: delta_z = 0.2d0
  double precision, parameter :: LambdaA = -60.0d0 ! choose 0 to turn surface interaction off
  double precision, parameter :: LambdaB = -5.0d0
  integer, parameter          :: monomer_types=2
  integer,parameter           :: gridx=20, gridy=35, gridb=20, gridz=100
  integer                     :: int_scheme, mesh_type
  double precision            :: sizex, sizey, sizez, volume, ds, dx, dy
  double precision, dimension(0:gridz)  :: dz, zr !discretization

end module Global
