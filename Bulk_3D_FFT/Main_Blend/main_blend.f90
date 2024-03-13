!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    MAIN PROGRAM (AB HOMOPOLYMER BLEND)
!           EXTERNAL POTENTIAL DYNAMICS SIMULATION
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
program main_blend
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    DECLARATIONS 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  use Global
  use Molecule_Homopolymer
  use Interactions_FloryHuggins
  use mtmod
!
  implicit none
!
! GLOBAL PARAMETERS: DECLARED / DEFINED IN MODULE GLOBAL
!
! double precision, parameter :: pi = 3.141592653589793d0
!
! integer, parameter          :: monomer_types=2
!
! integer, parameter  :: gridx, gridy, gridz ! lattice points in each direction
! double precision    :: sizex, sizey, sizez ! system size
! double precision    :: volume	             ! system volume
! double precision    :: dvol 		     ! volume element
! double precision    :: ds		     ! discretization along chain
!
  integer, parameter  :: polymer_types=monomer_types

  type(homopolymer), dimension(polymer_types) :: polymer
 
  double precision :: Ctotal    ! total dimensionless concentration of polymers

  double precision, dimension(polymer_types) :: mu  ! chemical potential of polymer 
  double precision, dimension(polymer_types) :: C   ! dimensionless concentration of polymer 
  double precision, dimension(polymer_types) :: number_of_chains
  double precision, dimension(polymer_types) :: chain_length

  ! double precision    ::  energy, interaction_energy
  double precision    ::  energy

  integer             ::  number_of_timesteps ! number of timesteps
  double precision    ::  timestep            ! length of timestep
  
  double precision, dimension(gridx,gridy,gridz,monomer_types) :: density
  double precision, dimension(gridx,gridy,gridz,monomer_types) :: field

  logical          :: qexist, qread, qepd
  character*24     :: file_name
  integer          :: polymer_label, monomer_type
  integer          :: polymer_ensemble    ! canonical:1, grand canonical:2
  integer          :: qinit   ! 1:sharp interface, 2: random densities, 3: random fields
  double precision :: dx, dy, dz

  integer          :: io_error, x,y,z,step, rgx,rgy,rgz, i
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
!
! ------ GENERAL INITIALIZATIONS 
!
  call Initialize

  C = number_of_chains/volume
  Ctotal = sum(C)
!
  dx = sizex/dble(gridx)
  dy = sizey/dble(gridy)
  dz = sizez/dble(gridz)
!
! ------ INITIALIZE HOMOPOLYMERS
!
  do polymer_label = 1,2
    monomer_type = polymer_label
    polymer_ensemble  = 1  ! canonical
!   polymer_ensemble  = 2  ! grand canonical
    call Initialize_Homopolymer(polymer(monomer_type), polymer_label, &
                     monomer_type, chain_length(polymer_label), &
                     polymer_ensemble, C(polymer_label), mu(polymer_label))
  end do
!
! ------ INITIALIZE INTERACTIONS
!
  call Initialize_Interactions
!
! ------ PREPARE OUTPUT FILES
!
  inquire(file='check.dat', exist=qexist)
  if (qexist) then
    open(unit=30,file='check.dat', status='old',action='write',position='append',iostat=io_error)
  else
    open(unit=30,file='check.dat', status='new',action='write',iostat=io_error)
  end if
  if(io_error /=0) then                   
        write(*,*) 'Error', io_error , ' while trying to open check.dat'
  end if
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   INITIAL CONFIGURATION
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  rgx = nint(chain_length(1)/dx)
  rgy = nint(chain_length(1)/dy)
  rgz = nint(chain_length(1)/dz)
  if (qread) then
    open (unit=33, file='cnf.in', status='old', action='read', iostat=io_error)
    if(io_error /=0) then                   
      write(*,*) 'Error', io_error , ' while trying to open cnf.in'
    end if
    read (33,*) field
    close (unit=33)
  else
    if (qinit.eq.1) then
       do z=1,gridz
         do y=1,gridy
           do x=1,gridx/2+1
             density(x,y,z,1) = 1.
             density(x,y,z,2) = 0.
           end do
           do x=gridx/2+1,gridx
             density(x,y,z,1) = 0.
             density(x,y,z,2) = 1.
           end do
         end do
       end do
       call Get_Conjugate_Fields(density,field)
     else if (qinit .eq. 2) then
       do i = 1,monomer_types
         do z=1,gridz
           do y=1,gridy
              do x=1,gridx
                density(x,y,z,i) = grnd()     
              end do
           end do
         end do
       end do
      density = Ctotal*(1.0d0 + (2.0d0*density -1.0d0)*0.001d0)
      call Get_Conjugate_Fields(density,field)
     else 
       do i = 1,monomer_types
         do z=1,gridz
           do y=1,gridy
              do x=1,gridx
                field(x,y,z,i) = grnd()     
              end do
           end do
         end do
       end do
    end if
  end if
!
  call Calculate_Densities(field,density)
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   EXTERNAL POTENTIAL DYNAMICS SIMULATION / SCF CALCULATION
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  if (qepd) then

    write(30,*) 'BEGIN EPD SIMULATION'
    write(30,*) 

    do step = 1,number_of_timesteps

       call EPD_Modified_Timestep_Euler(field,density,timestep)
!      call EPD_Timestep_Euler(field,timestep)
!      call EPD_Timestep_RungeKutta(field,timestep)
!
!      if (((step/10)*10).eq. step) then                  ! FOR TESTING
!       call Calculate_Energy                             ! FOR TESTING
!       print*, step, energy, polymer%rho, polymer%mu     ! FOR TESTING
!      end if                                             ! FOR TESTING

    end do

 else

    write(30,*) 'BEGIN SCF CALCULATION'
    write(30,*) 

    call Find_SCF_Saddle_Point(density,field)
    call Calculate_Densities(field,density)

  end if

  call Calculate_Energy

  if (qepd) then
    write(30,*) '   Final values after ', number_of_timesteps, 'steps'
    print*, ' Final values after ', number_of_timesteps, 'steps'
  else
    write(30,*) '   Final values '
    print*, ' Final values '

  end if
  write(30,*) '   Energy: ', energy, ' C: ', polymer%rho, ' mu_scf: ', polymer%mu
  write(30,*) 
!
  print*, '   Energy: ', energy, ' C: ', polymer%rho, ' mu_scf: ', polymer%mu
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   OUTPUT OF FINAL CONFIGURATION
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!  call Calculate_Densities(field,density)

      open(unit=60,file='dens.dat', status='replace',action='write',iostat=io_error)
      open(unit=61,file='field.dat', status='replace',action='write',iostat=io_error)
        do z=1,gridz
           do y=1,gridy
              do x=1,gridx
                 write(60, *) x*dx, '   ' , y*dy, '   ',z*dz, '   ', (density(x,y,z,i),i=1,monomer_types)
                 write(61, *) x*dx, '   ' , y*dy, '   ',z*dz, '   ', (field(x,y,z,i),i=1,monomer_types)
              end do
             write(60, *) 
             write(61, *) 
           end do
        end do
      close(unit=60)
      close(unit=61)
!
      open(unit=33, file='cnf.out', status = 'replace', action='write', iostat=io_error)
        write (33,*) field
      close (unit=33)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   FINISH (CLOSE FILES ETC.)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  close (unit=30)
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
        do polymer_label = 1,polymer_types
          read(20,*) number_of_chains(polymer_label)
          read(20,*) mu (polymer_label)
          read(20,*) chain_length(polymer_label)
          read(20,*)
        end do
        read(20,*) ds
        read(20,*)
        read(20,*) number_of_timesteps
        read(20,*) timestep
        read(20,*)
        read(20,*) qread  ! qread = .true. : read from file cnf.in
        read(20,*) qinit  ! only for qread=.false.  
            ! 1=interface, 2=random densities, 3=random fields
        read(20,*) qepd   ! qepd = .true. : EPD simulation; otherweise SCF calculation

      else
        write(*,*) 'Error', io_error , ' while trying to open ', file_name
      end if

    close(unit=20)
  
    return

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

    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in)  :: field
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(out) :: density

    double precision, dimension(gridx,gridy,gridz,monomer_types) :: expfield

    integer :: monomer_type, polymer_type, i

    call Calculate_Expfield(field,expfield)
    do i = 1,polymer_types
      call Propagate_Homopolymer(polymer(i),expfield(:,:,:,polymer(i)%monomer_type))
    end do
 
    density = 0.0d0
    do i = 1,polymer_types
       monomer_type = polymer(i)%monomer_type
       density(:,:,:,monomer_type) = density(:,:,:,monomer_type) + polymer(i)%density(:,:,:)
    end do
 
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
    double precision :: chain_energy
    integer :: i

    chain_energy = - dvol*sum(density*field)
    do i = 1,polymer_types
       chain_energy = chain_energy + polymer(i)%energy 
    end do
    energy = chain_energy + interaction_energy(density)                

    return

  end subroutine Calculate_Energy
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   SUBROUTINE TO CALCULATE FIELDS CORRESPONDING TO GIVEN DENSITIES  
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Get_Conjugate_Fields(density,field)
!
   use Propagate

    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: field
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in)    :: density

    double precision, dimension(gridx,gridy,gridz,monomer_types) :: dfield
    double precision, dimension(gridx,gridy,gridz,monomer_types) :: density0, newdensity

    double precision :: lambda, accuracy, accuracy_scf, small
    integer :: mixing_type, iter_status, iter, iter_max

    accuracy_scf = 0.01
    mixing_type = 3
    iter_max = 100000
    lambda = 0.5

!    small = 1.d-6
!    where (density > small)
!      density0 = density
!    elsewhere
!      density0 = small
!    end where
!
    iter_status = 1

    do iter = 1,iter_max

      call Calculate_Densities(field,newdensity)
      dfield = newdensity - density
      call Check_Accuracy(dfield,accuracy)
!      if ( ( (iter/100)*100) .eq. iter) then
!        print*, iter, accuracy
!      end if

      if (.not.(accuracy > accuracy_scf) ) then  ! double negation to deal with NaNs
!        PRINT*, 'exit Get_Conjugate_Field after', iter, 'steps '
!        PRINT*, 'Accuracy: ', accuracy
        exit
      end if
      call Mix_Fields(mixing_type,field,dfield,lambda,iter_status)
 
    end do
    
    if (.not.(accuracy < accuracy_scf)) then  ! double negation to deal with NaNs
      print*, 'Caution: Conjugate fields not found within ',iter,' iteration steps!'
      print*, 'Accuracy: ', accuracy
    end if
 
    return

  end subroutine Get_Conjugate_Fields
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   SUBROUTINE TO FIND SCF SADDLE POINT
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Find_SCF_Saddle_Point(density,field)
!
   use Propagate
!   use Interactions_FloryHuggins

    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: field
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: density

    double precision, dimension(gridx,gridy,gridz,monomer_types) :: dfield
    double precision, dimension(gridx,gridy,gridz,monomer_types) :: newdensity, newfield

    double precision :: lambda, accuracy, accuracy_scf
    integer :: mixing_type, iter_status, iter, iter_max

    accuracy_scf = 0.0000001
    mixing_type = 2
    iter_max = 100000
    lambda = 0.5

    iter_status = 1

    do iter = 1,iter_max

      call Calculate_Densities(field,newdensity)
      call Calculate_Mean_Fields(newdensity,newfield)
      dfield = newfield - field
      call Check_Accuracy(dfield,accuracy)
!      if ( ( (iter/100)*100) .eq. iter) then
!        print*, iter, accuracy
!      end if

      if (.not.(accuracy > accuracy_scf) ) then  ! double negation to deal with NaNs
        print*, 'SCF saddle point found after ',iter,' iteration steps!'
        print*, 'Accuracy: ', accuracy
        exit
      end if
      call Mix_Fields(mixing_type,field,dfield,lambda,iter_status)
 
    end do
    
    if (.not.(accuracy < accuracy_scf)) then  ! double negation to deal with NaNs
      print*, 'Caution: SCF saddle point not found within ',iter,' iteration steps!'
      print*, 'Accuracy: ', accuracy
    end if

    density = newdensity
 
    return

  end subroutine Find_SCF_Saddle_Point
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!          SUBROUTINE TO MIX FIELDS IN SCF ITERATION
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Mix_Fields(mixing_type,field,dfield,lambda,iter_status)
!
   use Iterate
!
     implicit none

     integer, intent(in)          :: mixing_type
     double precision, intent(in) :: lambda
     integer, intent(inout)       :: iter_status

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
              call mixing_simple(var,dvar,nvar,lambda,iter_status)
           case (2)
              call mixing_lambda(var,dvar,nvar,lambda,iter_status)
           case (3)
              call mixing_anderson(var,dvar,nvar,lambda,iter_status)
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

      iter_status = iter_status + 1

    return

  end subroutine Mix_Fields
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           CALCULATE ACHIEVED ACCURACY
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Check_Accuracy(dfield,accuracy)

    integer :: x,y,z
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in) :: dfield
    double precision, intent(out)       :: accuracy

    accuracy=sum(abs(dfield))*dvol

    return
  end subroutine check_accuracy
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   EPD TIME STEP 
!          USING FORWARD EULER METHOD
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine EPD_Timestep_Euler (field,timestep)
!
     implicit none

     double precision, intent(in) :: timestep
      
     double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: field

     double precision, dimension(gridx,gridy,gridz,monomer_types) :: dfield
!
     call Calculate_Field_Derivative(field,dfield)

     field = field + timestep*dfield

    return

  end subroutine EPD_Timestep_Euler
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   MODIFIED EPD TIME STEP 
!          USING FORWARD EULER METHOD
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine EPD_Modified_Timestep_Euler (field,density,timestep)
!
   use Fourier_fftw3

     implicit none

     double precision, intent(in) :: timestep
      
     double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: field
     double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: density

     double precision, dimension(gridx,gridy,gridz,monomer_types) :: newfield
     double precision, dimension(gridx,gridy,gridz,monomer_types) :: mulocal
     double precision, dimension(gridx,gridy,gridz,monomer_types) :: olddensity, newdensity
!
     double precision, dimension(gridx,gridy,gridz,monomer_types,3) :: dfield
     double precision, dimension(gridx,gridy,gridz,monomer_types,3) :: ddensity
!
     double precision, dimension(gridx,gridy,gridz)   :: dummy 
     double complex, dimension(gridx/2+1,gridy,gridz) :: dummyfourier, dummyfourier0
!
     double precision :: eps
     integer :: x,y,z,n,monomer_type, spacedim

!
     eps = 0.0001

     call Calculate_Densities(field,olddensity)
     call Calculate_Mean_Fields(olddensity,newfield)
     mulocal = newfield - field

     do monomer_type = 1, monomer_types
       dummy = mulocal(:,:,:,monomer_type) 
       call Real2Fourier(dummy,dummyfourier0)
       do spacedim = 1,3
         dummyfourier = nabla(:,:,:,spacedim)*dummyfourier0
         call Fourier2Real(dummy,dummyfourier)
         dfield(:,:,:,monomer_type,spacedim) = dummy
       end do
     end do

     do spacedim = 1,3
       newfield = field - eps*dfield(:,:,:,:,spacedim)
       call Calculate_Densities(newfield,newdensity)
       ddensity(:,:,:,:,spacedim) = newdensity - olddensity
       do monomer_type = 1, monomer_types
         dummy = ddensity(:,:,:,monomer_type,spacedim) 
         call Real2Fourier(dummy,dummyfourier)
         dummyfourier = nabla(:,:,:,spacedim)*dummyfourier
         call Fourier2Real(dummy,dummyfourier)
         ddensity(:,:,:,monomer_type,spacedim) = dummy
       end do
     end do

!    must distinguish between olddensity and density, 
!       otherwise small density increments have no effect at all!
!
     density = density + sum(ddensity,dim=5) * timestep/eps
     call Get_Conjugate_Fields(density,field)

    return

  end subroutine EPD_Modified_Timestep_Euler
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   EPD TIME STEP 
!          USING FOURTH ORDER RUNGE-KUTTA METHOD
!
!          NB: DOES NOT SEEM SUPERIOR TO EULER METHOD
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine EPD_Timestep_RungeKutta (field,timestep)
!
     implicit none

     double precision, intent(in) :: timestep
      
     double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: field
!
     double precision, dimension(gridx,gridy,gridz,monomer_types) :: dummy_field
     double precision, dimension(gridx,gridy,gridz,monomer_types) :: df1,df2,df3,df4
!
     call Calculate_Field_Derivative(field,df1)
     dummy_field = field + 0.5d0*timestep * df1
     call Calculate_Field_Derivative(dummy_field,df2)
     dummy_field = field + 0.5d0*timestep * df2
     call Calculate_Field_Derivative(dummy_field,df3)
     dummy_field = field + timestep * df3
     call Calculate_Field_Derivative(dummy_field,df4)

     field = field + timestep*(df1/6.0d0+df2/3.0d0+df3/3.0d0+df4/6.0d0)

    return

  end subroutine EPD_Timestep_RungeKutta
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   CALCULATE FIELD DERIVATIVE FOR EPD TIME INTEGRATION
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Calculate_Field_Derivative (field,dfield)
!
   use Fourier_fftw3
!
     implicit none

     double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in)  :: field
     double precision, dimension(gridx,gridy,gridz,monomer_types), intent(out) :: dfield

     double precision, dimension(gridx,gridy,gridz,monomer_types) :: density
     double precision, dimension(gridx,gridy,gridz,monomer_types) :: newfield
!
     double precision, dimension(gridx,gridy,gridz)   :: dummy 
     double complex, dimension(gridx/2+1,gridy,gridz) :: dummyfourier
!
     integer :: x,y,z,n,monomer_type

     call Calculate_Densities(field,density)
     call Calculate_Mean_Fields(density,newfield)
     dfield = newfield - field

     do monomer_type = 1, monomer_types
       dummy = dfield(:,:,:,monomer_type) 
       call Real2Fourier(dummy,dummyfourier)
       dummyfourier = - laplace*dummyfourier
       call Fourier2Real(dummy,dummyfourier)
       dfield(:,:,:,monomer_type) = dummy
     end do

    return

  end subroutine Calculate_Field_Derivative
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	     T H E    E N D 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end program main_blend
