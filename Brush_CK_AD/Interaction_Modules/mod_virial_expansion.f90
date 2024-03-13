!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!        ROUTINES FOR ONE COMPONENT VIRIAL EXPANSION
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

module Interactions_Virial_Expansion
!
  use Global
  use Adaptive
!
  implicit none
!
  private
!
  integer, parameter     :: components = monomer_types     ! solvent is component 0
!
  public :: Initialize_Interactions
  public :: Calculate_Mean_Ext_Fields
  public :: External_Energy
  public :: Interaction_Energy
  public :: Check_Volume_Fraction
!
contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           INITIALIZE INTERACTIONS
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Initialize_Interactions
!
    character*30 :: file_name
    integer      :: io_error
    integer      :: component, component2
!
    file_name='input_interactions'
    open(unit=20,file=file_name,status='old',action='read',iostat=io_error) 
!
   
    if(io_error == 0) then
      read(20,*) ww
      read(20,*)
      read(20,*) monomer_volume(0)
      read(20,*) monomer_volume(1)
    else
      print*, 'Error', io_error, ' while trying to open', file_name 
    end if


  end subroutine Initialize_Interactions
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           CALCULATE INTERACTION ENERGY
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  double precision function Interaction_Energy(density)
!
    double precision, dimension(gridx,gridy,0:gridz,components), intent(in) :: density
!
    double precision, dimension(gridx,gridy,0:gridz,components) :: volume_fraction

    integer          :: component, component2 ,i



    call integrator(density(:,:,:,1)**2,Interaction_Energy)

    interaction_energy= Interaction_Energy * ww * 0.5d0

 
  end function Interaction_Energy


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           CHECK THAT SUM OF AVERAGE VOLUME FRACTIONS GIVES ONE.
!		    (FOR USE IN THE CANONICAL ENSEMBLE)
!           THIN FILMS: COUNT ONLY ACCESSIBLE GRID POINTS
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Check_Volume_Fraction(density)
!
    double precision, dimension(gridx,gridy,0:gridz,components), intent(in) :: density
!
    double precision, dimension(gridx,gridy,0:gridz,components) :: volume_fraction
    double precision  :: average_volume_fraction
    integer           :: component, i
!
    do component = 1,components
      volume_fraction(:,:,:,component) = density(:,:,:,component)*monomer_volume(component)
    end do
!
    average_volume_fraction = sum(volume_fraction) / dble(gridx*gridy*gridz)

  end subroutine Check_Volume_Fraction


!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           CALCULATE MEAN FIELDS WITH WALL-INTERACTION
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Calculate_Mean_Ext_Fields(density,field)
!
    double precision, dimension(gridx,gridy,0:gridz,components), intent(in) :: density
    double precision, dimension(gridx,gridy,0:gridz,components), intent(out):: field
!
    double precision, dimension(gridx,gridy,0:gridz,components)               :: volume_fraction
    double precision, dimension(gridx,gridy,0:gridz,components)             :: surface_potential


    integer          :: component, component2

    do component = 1,components
      volume_fraction(:,:,:,component) = density(:,:,:,component)*monomer_volume(component)
    end do

    Call Calculate_Potential(surface_potential)
    
 

    do component = 1,components

      field(:,:,:,component) =  ( ww * density(:,:,:,component)+surface_potential(:,:,:,component) )

    end do

 
  end subroutine Calculate_Mean_Ext_Fields
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           CALCULATE EXTERNAL ENERGY
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  double precision function External_Energy(density)
!
    double precision, dimension(gridx,gridy,0:gridz,components), intent(in) :: density
    double precision, dimension(gridx,gridy,0:gridz,components) :: volume_fraction
    double precision, dimension(gridx,gridy,0:gridz,components) :: surface_potential
    integer               :: component
!
!
    Call Calculate_Potential(surface_potential)
!
    do component = 1,components
      volume_fraction(:,:,:,component) = density(:,:,:,component)*monomer_volume(component)
    end do

    call integrator(sum(surface_potential*volume_fraction,dim=4),external_energy)

  end function External_Energy
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!                 EXTERNAL FIELD
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Calculate_Potential(surface_potential)
!
    double precision, dimension(0:gridz) :: potential_shape
    double precision, dimension(gridx,gridy,0:gridz,components) :: surface_potential
    double precision :: z0
    integer :: i

    potential_shape = 0.0d0

    z0 = sizez - delta_z
    do i = 0,gridz
       if ( (zr(i) .le. delta_z) .and. ((wall_scheme .eq. 0) .or. (wall_scheme .eq. 2) )) then
          potential_shape(i) = 1.d0 + cos(pi*zr(i)/delta_z)
       elseif ((zr(i) .gt. delta_z) .and. (zr(i) .lt. z0) ) then
          potential_shape(i) = 0.d0
       elseif ((zr(i) .ge. z0) .and. (zr(i).le. sizez) .and. ((wall_scheme .eq. 1) .or. (wall_scheme .eq. 2) )) then
          potential_shape(i) = 1.d0 + cos(pi*(sizez-zr(i))/delta_z)
       end if
       surface_potential(:,:,i,1) = LambdaA * potential_shape(i)
    end do


    open(unit=72,file='test_potential.dat', status='replace',action='write')
        do i=0,gridz
          write(72, *) i, '  ', zr(i) , '   ' ,surface_potential(1,1,i,1)
        end do
    close(unit=72)

    return

  end subroutine Calculate_Potential



end module Interactions_Virial_Expansion
