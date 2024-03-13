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
  double precision, parameter :: pi = 3.14159265358979311599796346854d0
  integer, parameter          :: monomer_types=1  
  integer,parameter           :: gridx=1, gridy=1, gridz=1000
  integer                     :: int_scheme, mesh_type ,delta_rep ,stat_job ,segments , ds_scheme, wall_scheme,eps0
  double precision            :: sizex, sizey, sizez, ds, volume, gra_pt , eps ,LambdaA  
  double precision, dimension(0:gridz)  :: dz, zr !discretization
  double precision            :: monomer_volume(0:monomer_types)   ! monomer volumes
  double precision            :: ww, delta_z ! interaction parameters 
end module Global
