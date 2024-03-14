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
!
  integer, parameter          :: monomer_types=2
!
  integer,parameter           :: gridx=86, gridy=1, gridz=75
!  integer,parameter           :: gridx=86, gridy=1, gridz=1
  integer,parameter           :: ngrid = gridx*gridy*gridz
  double precision            :: ds 
  double precision            :: sizex, sizey, sizez, volume, dvol
!
  integer, parameter          :: mixing_dim=10 ! should at least be 2

end module Global
