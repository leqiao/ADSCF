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
  integer,parameter           :: gridx=64, gridy=16, gridz=1
  double precision            :: ds 
  double precision            :: sizex, sizey, sizez, volume, dvol

end module Global
