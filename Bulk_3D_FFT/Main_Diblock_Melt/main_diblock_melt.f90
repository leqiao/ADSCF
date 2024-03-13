!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    MAIN PROGRAM (DIBLOCK COPOLYMER MELT)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
program main_diblock_melt
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    DECLARATIONS 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  use Global
  use Molecule_Diblock_Copolymer
  use Interactions_FloryHuggins
  use mtmod
!
  implicit none
!
! GLOBAL PARAMETERS: DECLARED / DEFINED IN MODULE GLOBAL
!
! double precision, parameter :: pi = 3.141592653589793d0
! integer, parameter  :: gridx, gridy, gridz 		! lattice points in each direction
! double precision    :: sizex, sizey, sizez     	! system size
! double precision    :: volume	               		! system volume
! double precision    :: dvol 				! volume element
! double precision    :: ds		                ! discretization along chain
!
  type(diblock_copolymer) :: diblock
!
  double precision    ::  C      ! dimensionless polymer concentration
  double precision    ::  mu     ! chemical potential of polymer

  double precision    ::  energy

  integer             ::  itermax_scf    ! maximum number of iterations
  integer             ::  iterations     ! actual number of iterations
  double precision    ::  accuracy_scf   ! desired accuracy in iterations
  double precision    ::  accuracy       ! actual accuracy
  double precision    ::  lambda         ! mixing parameters
  integer             ::  mixing_type    ! type of mixing
!      -------------  1: simple mixing, 2: lambda mixing, 3: anderson mixing ----------

  integer             ::  io_error,x,y,z,iter,i

  double precision, dimension(gridx,gridy,gridz,monomer_types) :: density
  double precision, dimension(gridx,gridy,gridz,monomer_types) :: newdensity
  double precision, dimension(gridx,gridy,gridz,monomer_types) :: field
  double precision, dimension(gridx,gridy,gridz,monomer_types) :: newfield, dfield

  double precision, dimension(gridx,gridy,gridz,monomer_types) :: saddledensity
  double precision, dimension(gridx,gridy,gridz,monomer_types) :: saddlefield
!  double precision :: interaction_energy

  double precision, dimension(gridx,gridy,gridz,monomer_types) :: density_ref
  double precision, dimension(gridx,gridy,gridz,monomer_types) :: field_ref
  integer          :: dum
  double precision :: dumA

  logical          :: qexist
  character*24     :: file_name
  double precision :: block_lengthA, block_lengthB
  integer          :: polymer_label, monomer_typeA, monomer_typeB
  integer          :: polymer_ensemble   ! canonical:1, grand canonical:2
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	     INITIALIZATIONS
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! ------ READ INPUT PARAMETERS 
!
file_name='input_general'
call read_input_parameters(file_name)
! ------ GENERAL INITIALIZATIONS 
!
call Initialize
!
! ------ INITIALIZE DIBLOCK COPOLYMER
!
polymer_label = 1
monomer_typeA = 1
monomer_typeB = 2
block_lengthA = 0.7d0
block_lengthB = 0.3d0
!    polymer_ensemble  = 1  ! canonical
polymer_ensemble  = 2  ! grand canonical
call Initialize_Diblock(diblock, polymer_label, &
                 (/monomer_typeA, monomer_typeB/),&
                 (/block_lengthA, block_lengthB/),&
                 polymer_ensemble,C,mu)
!
! ------ INITIALIZE INTERACTIONS
!
call Initialize_Interactions
!
! ------ CHECKS
!
! call Calculate_Densities(field,density)
! call Check_Volume_Fraction(density)  ! for fully canonical with Flory Huggins interactions

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   DETERMINATION OF SADDLE POINT
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                LOAD PREFERENCE DESITY IF YOU HAVE ONE
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! open(unit=3,file='initial_dens.dat',status='old',action='read',iostat=io_error)
! do z=1,gridz
!     do y=1,gridy
!         do x=1,gridx
!             read(3,*) dum, dum, dum, density(x,y,z,1), density(x,y,z,2)
!         end do
!     end do
! end do
! close(unit=3)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!            IF NO EXSITING PREFERENCE DESITY, THEN GUESS ONE FROM BELOW FUNCTION
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
do x = 1,gridx
    do z = 1,gridz
        density(x,:,z,1) = block_lengthA*(1.0d0 + 0.5d0*cos(dble(x)/dble(gridx)*2*pi)&
            *cos(dble(z)/dble(gridz)*2*pi))
        density(x,:,z,2) = block_lengthB*(1.0d0 - 0.5d0*cos(dble(x)/dble(gridx)*2*pi)&
            *cos(dble(z)/dble(gridz)*2*pi))
        density(x,:,:,1) = diblock%rho*block_lengthA*(1.0d0 + 0.5d0*cos(dble(x)/dble(gridx)*pi))
        density(x,:,:,2) = diblock%rho*block_lengthB*(1.0d0 - 0.5d0*cos(dble(x)/dble(gridx)*pi))
    end do
end do

!---- TEST VARY SIZE LENGTH
!  call Calculate_Mean_Fields(density,field)
call Calculate_Mean_Ext_Fields(density,field)
!
! -------- ITERATION LOOP
!
write(30,*) 'FIND SADDLE POINT'
write(30,*) 
!
mixing_type = 2
newdensity = density
newfield   = field
 
do iter = 1,itermax_scf

    call Calculate_Densities(field,newdensity)
    !
    call Check_Accuracy(density,newdensity,accuracy)

      call Calculate_Energy                                   ! FOR TESTING
      print*, iter, energy, accuracy, diblock%rho, diblock%mu ! FOR TESTING

    if ( .not.(accuracy>accuracy_scf) ) then     !double negation to deal with NaNs
      exit
    end if

    !    call Calculate_Mean_Fields(newdensity,newfield)
    call Calculate_Mean_Ext_Fields(newdensity,newfield)
    dfield = newfield - field
    Call Mix_Fields(mixing_type,field,dfield,lambda,iter)
    density=newdensity

end do
iterations = iter

saddledensity = newdensity
saddlefield   = newfield

call Check_Volume_Fraction(saddledensity)  ! for fully canonical with Flory Huggins interactions

call Calculate_Energy

if (.not.(accuracy < accuracy_scf)) then  ! DOUBLE NEGATION DEALS WITH NaNs
    print*, 'Caution: SCF saddle point not found within ',iterations,' iteration steps!'&
       ,' mu_scf: ', diblock%mu
    call finish
    stop 1
end if
!
print*, '   Saddle point found after', iterations, 'iterations'
print*, '   Accuracy: ', accuracy, ' Energy: ', energy, ' C: ', diblock%rho,  &
      ' mu_scf: ', diblock%mu
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   OUTPUT OF SADDLE POINT RESULTS
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  open(unit=60,file='saddle_dens.dat', status='replace',action='write',iostat=io_error)
  open(unit=61,file='saddle_field.dat', status='replace',action='write',iostat=io_error)
!
  open(unit=70, file='energy_f07-1p.dat', status='unknown', action='write', iostat=io_error)
!
  open (unit=71, file='density_f07-1p.dat', status='unknown', action='write', iostat=io_error)

    do z=1,gridz
       do y=1,gridy
          do x=1,gridx
             write(60, *) x, '   ' , y, '   ',z, '   ', (saddledensity(x,y,z,i),i=1,monomer_types)
             write(61, *) x, '   ' , y, '   ',z, '   ', (saddlefield(x,y,z,i),i=1,monomer_types)
          end do
          write(60,*)
       end do
    end do
    close(unit=60)
    close(unit=61)
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!                 WRITE ENERGY INTO FILE
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    write(70,*) sizex, '    ', sizey, '   ',sizez, '    ', Energy, '    ', diblock%rho   ! Grand Canonical
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!        WRITE DENSITIES INTO FILE
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!        do z=1,gridz
!           do x=1,gridx
!              write(71,*) sizex, '    ', x, '    ', sizez, '   ', z, '    ', (saddledensity(x,1,z,i),i=1,monomer_types)
!           end do
!           write(71,*)
!        end do
!
!
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  close(unit=70)
  close(unit=71)
  close(unit=72)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   FINISH (CLOSE FILES ETC.)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
call finish
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   END OF MAIN PROGRAM. NOW INNER SUBROUTINES
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   SUBROUTINE TO READ INPUT PARAMETERS  ! THIS DEPENDS ON SYSTEM !!!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine read_input_parameters(file_name)
    
    character*24, intent(in) :: file_name
    integer                  :: io_error

    open(unit=20,file=file_name,status='old',action='read',iostat=io_error)

      if(io_error ==0) then

        read(20,*) sizex
        read(20,*) sizey
        read(20,*) sizez
        read(20,*)
        read(20,*) C
        read(20,*) mu
        read(20,*)
        read(20,*) ds
        read(20,*)
        read(20,*) itermax_scf
        read(20,*) accuracy_scf
        read(20,*)
        read(20,*) lambda

      else
        write(*,*) 'Error', io_error , ' while trying to open ', file_name
      end if

    close(unit=20)
  
  end subroutine read_input_parameters
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   SUBROUTINE TO CALCULATE DENSITIES  ! THIS DEPENDS ON SYSTEM !!!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Calculate_Densities(field,density)
!
   use Propagate
!
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in)  :: field
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(out) :: density

    double precision, dimension(gridx,gridy,gridz,monomer_types) :: expfield

    integer :: monomer_type,i

    call Calculate_Expfield(field,expfield)
    call Propagate_Diblock(diblock,expfield(:,:,:,diblock%monomer_type))

    density = 0.0d0
    do i = 1,2
      monomer_type = diblock%monomer_type(i)
      density(:,:,:,monomer_type) = density(:,:,:,monomer_type) + diblock%density(:,:,:,i)
    end do

    density(:,:,:,diblock%monomer_type) = diblock%density

    return

  end subroutine Calculate_Densities
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!          SUBROUTINE TO CALCULATE ENERGY  ! THIS DEPENDS ON SYSTEM !!!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Calculate_Energy
!
    double precision              :: chain_energy

    print*, diblock%energy
    print*, dvol*sum(saddledensity*saddlefield)
    print*, dvol
    chain_energy = diblock%energy - dvol*sum(saddledensity*saddlefield)
    energy = chain_energy + interaction_energy(density)

    return

  end subroutine Calculate_Energy
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   SUBROUTINE TO MIX FIELDS IN SCF ITERATION
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Mix_Fields(mixing_type,field,dfield,lambda,iter)
!
   use Iterate
!
     implicit none

     integer, intent(in)               :: mixing_type, iter
     double precision, intent(in)      :: lambda
      
     double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: field
     double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: dfield
!
     integer, parameter :: nvar = gridx*gridy*gridz*monomer_types
     double precision, dimension(nvar,0:mixing_dim) :: var,dvar 

     integer :: x,y,z,n,monomer_type

      n=0
      do monomer_type = 1,monomer_types 
        do z = 1,gridz
         do y = 1,gridy
          do x = 1,gridx
            n=n+1
             var(n,0)  = field(x,y,z,monomer_type)
             dvar(n,0) = dfield(x,y,z,monomer_type)
          end do
         end do
        end do
      end do

      mixing: select case (mixing_type)
           case (1)
              call mixing_simple(var,dvar,nvar,lambda,iter)
           case (2)
              call mixing_lambda(var,dvar,nvar,lambda,iter)
           case (3)
              call mixing_anderson(var,dvar,nvar,lambda,iter)
           case default
              print*, 'Invalid mixing option in find_saddle_point!'
              stop
      end select mixing

      n=0
      do monomer_type = 1,monomer_types 
        do z = 1,gridz
         do y = 1,gridy
          do x = 1,gridx
            n=n+1
             field(x,y,z,monomer_type) = var(n,0)  
             dfield(x,y,z,monomer_type) = dvar(n,0) 
          end do
         end do
        end do
      end do

    return

  end subroutine Mix_Fields
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           CALCULATE ACHIEVED ACCURACY
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Check_Accuracy(density,newdensity,accuracy)

    integer :: x,y,z
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in) ::    density
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in) :: newdensity
    double precision, intent(out)       :: accuracy

    accuracy=sum((density-newdensity)*(density-newdensity))

    return
  end subroutine check_accuracy
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	     T H E    E N D 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end program main_diblock_melt
