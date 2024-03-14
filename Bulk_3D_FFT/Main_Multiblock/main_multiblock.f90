!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    MAIN PROGRAM (MULTIBLOCK COPOLYMER MELT)
!
!       ALTERNATIVELY
!       
!          SCF CALCULATION OR
!          DYNAMIC DENSITY FUNCTIONAL CALCULATION WITH IMPOSED MOBILITY MATRIX
!
!          References (for DDFT): 
!              S. Mantha, S. Qi, F. Schmid, Macromolecules 53, 3409 (2020)
!                 doi:10.1021/acs.macromol.0c00130
!              F. Schmid, B. Li, Polymers 12, 2205 (2020)
!                 doi:10.3390/polym12102205
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
program main_multiblock
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    DECLARATIONS 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  use Global
  use Molecule_Multiblock_Copolymer
  use Interactions_FloryHuggins
  use mtmod
!
  implicit none
!
! GLOBAL PARAMETERS: DECLARED / DEFINED IN MODULE GLOBAL
!
! double precision, parameter :: pi = 3.141592653589793d0
! integer, parameter  :: gridx, gridy, gridz 		! lattice points in each direction
! integer, parameter  :: ngrid                      ! total number of lattice points
! double precision    :: sizex, sizey, sizez     	! system size
! double precision    :: volume	               		! system volume
! double precision    :: dvol 				    ! volume element
! double precision    :: ds		                ! discretization along chain
! integer, paramter   :: mixing_dim             ! Mixing dimension for Anderson mixing
!
  type(multiblock_copolymer) :: multiblock
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

  double precision    ::  wcutoff        ! cutoff for DDFT calculations (numerical stability)

  double precision, dimension(gridx,gridy,gridz,monomer_types) :: density
  double precision, dimension(gridx,gridy,gridz,monomer_types) :: field
!
! REMOVE LATER
  double precision, dimension(gridx,gridy,gridz,monomer_types) :: newdensity, newfield, dfield
! REMOVE LATER

!  double precision :: interaction_energy

  double precision, dimension(gridx,gridy,gridz,monomer_types) :: potential

  double precision, dimension(gridx/2+1,gridy,gridz,monomer_types,monomer_types) :: mobility_qq


  logical           :: qexist, qread, qddft, qdmin
  character*24      :: file_name
  integer           :: polymer_label
  integer           :: polymer_ensemble   ! canonical:1, grand canonical:2
  integer           :: qinit   ! 1: cosine, 2: random densities, 3: random fields
!  integer,parameter :: nblocks = 10 
   integer,parameter :: nblocks = 6
!  integer,parameter :: nblocks = 7
!  integer,parameter :: nblocks = 2
  double precision, dimension(nblocks) :: block_length
  integer, dimension(nblocks)          :: block_monomer_type
!
  integer          :: number_of_timesteps ! number of time steps
  double precision :: timestep            ! length of time step
  double precision :: start_time          ! start time

  double precision :: dx, dy, dz, dum
  integer          :: io_error, x,y,z, step, i,i0, gcfstat
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
!
! ------ INITIALIZE MULTIBLOCK COPOLYMER
!
    polymer_label = 1
!    block_monomer_type = (/1,2/)
      block_monomer_type = (/1,2,1,2,1,2/)
!     block_monomer_type = (/1,2,1,2,1,2,1/)
!    block_monomer_type = (/1,2,1,2,1,2,1,2,1,2/)
!    block_length = (/0.5d0,0.5d0/)
     block_length = (/0.5d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0/)
!    block_length = (/0.1d0,0.1d0,0.1d0,0.4d0,0.1d0,0.1d0,0.1d0/)
!    block_length = (/0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0,0.1d0/)
    
    polymer_ensemble  = 1  ! canonical
!    polymer_ensemble  = 2  ! grand canonical
    call Initialize_Multiblock(multiblock, polymer_label, nblocks, &
                     block_monomer_type, block_length, &
                     polymer_ensemble,C,mu)
!
! ------ INITIALIZE INTERACTIONS
!
  call Initialize_Interactions
!
! ------ READ AND INITIALIZE MOBILITIES 
!
  if (qddft) then
    file_name='input_mobilities'
    call Read_Mobilities(file_name)
  end if
!
! ------ CHECKS
!
  call Calculate_Densities(field,density)

  call Check_Volume_Fraction(density)  ! for fully canonical with Flory Huggins interactions

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
!      INITIAL CONFIGURATION
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  if (qread) then
     open (unit=33, file='cnf.in', status='old', action='read', iostat=io_error)
     if(io_error /=0) then
       write(*,*) 'Error', io_error , ' while trying to open cnf.in'
     end if
     read (33,*) field
     close (unit=33)
     call Calculate_Densities(field,density)
   else
    if (qinit.eq.1) then
      do x = 1,gridx
         do z = 1,gridz
           density(x,:,z,1) = multiblock%rho &
!            * ( block_length(1)+block_length(3)+block_length(5) ) &
             *  block_length(1) &
!             * ( 1.0d0 + 0.2d0*cos(dble(x)/dble(gridx)*2*pi) )
             * ( 1.0d0 - 0.2d0*cos(dble(x)/dble(gridx)*2*pi)*cos(dble(z)/dble(gridz)*2*pi) )
           density(x,:,z,2) = 1.0d0 - density(x,:,z,1)
         end do
       end do
!       call Get_Conjugate_Fields(density,field,1.d-6,100000,3,gcfstat)
       call Calculate_Mean_Fields(density,field)
       call Calculate_Densities(field,density)
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

        density(:,:,:,1) = multiblock%rho &
!           *(block_length(1)+block_length(3)+block_length(5)) &
           *  block_length(1) &
           *(1.0d0 + (2.0d0*density(:,:,:,1) -1.0d0)*0.01d0)
        density(:,:,:,2) = 1.0d0 - density(:,:,:,1)
!             --- densities not yet normalized
       call Get_Conjugate_Fields(density,field,0.001d0,10000,3,gcfstat)
       call Calculate_Densities(field,density)
       call Get_Conjugate_Fields(density,field,0.001d0,100000,3,gcfstat)
 
!             --- continue with normalized densities
     else if (qinit .eq. 3) then
        do i = 1,monomer_types
          do z=1,gridz
            do y=1,gridy
               do x=1,gridx
                 field(x,y,z,i) = grnd()
               end do
            end do
          end do
        end do
       call Calculate_Densities(field,density)
     else
       print*, 'Initialization not well-defined!'
       stop
     end if
   end if
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   OUTPUT OF INITIAL CONFIGURATION
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      open(unit=60,file='initial_dens.dat', status='replace',action='write',iostat=io_error)
      open(unit=61,file='initial_field.dat', status='replace',action='write',iostat=io_error)
        dx = sizex/float(gridx)
        dy = sizey/float(gridy)
        dz = sizez/float(gridz)
        do x=1,gridx
           do y=1,gridy
              do z=1,gridz
                 write(60, *) x*dx, '   ' , y*dy, '   ',z*dz, '   ', &
                             (density(x,y,z,i),i=1,monomer_types)
                 write(61, *) x*dx, '   ' , y*dy, '   ',z*dz, '   ', &
                             (field(x,y,z,i),i=1,monomer_types)
              end do
           end do
           write(60, *)
           write(61, *)
        end do
      close(unit=60)
      close(unit=61)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!      DDFT SIMULATION / SCF CALCULATION
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  if (qddft) then

    open(unit=10,file='profiles.dat', status='replace',action='write',iostat=io_error)
    open(unit=40,file='energy.dat', status='replace',action='write',iostat=io_error)

    print*, 'BEGIN DDFT SIMULATION'
    write(30,*) 'BEGIN DDFT SIMULATION'
    write(30,*)

    qdmin = .true.

    do step = 1,number_of_timesteps

      call DDFT_Timestep_Euler(density,field,mobility_qq,timestep,qdmin)
!
      if (((step/100)*100).eq. step) then                 
       call Calculate_Energy                                
       write(40,*) start_time + step*timestep, energy, multiblock%rho, multiblock%mu  
       print*, start_time + step*timestep, energy, multiblock%rho, multiblock%mu  
      end if                                                
!      
      if (((step/2000)*2000).eq. step) then                 ! FOR TESTING
       call Calculate_Energy                                ! FOR TESTING
       write(30,*) step, energy, multiblock%rho, multiblock%mu  ! FOR TESTING
       call Check_Volume_Fraction(density)                  ! FOR TESTING
      end if                                                ! FOR TESTING
!
      i = nint(start_time/timestep) + step
      i0 = nint(0.1/timestep)
      if (((i/i0)*i0).eq. i) then                   ! FOR TESTING
        dx = sizex/float(gridx)
        dy = sizey/float(gridy)
        dz = sizez/float(gridz)
        do x=1,gridx
           do y=1,gridy
              do z=1,gridz
                 write(10, *) start_time + step*timestep,'  ',x*dx,'  ',y*dy,'  ',z*dz,'  ', &
                             (density(x,y,z,i),i=1,monomer_types)
              end do
           end do
           write(10, *)
        end do
!
        open(unit=33, file='cnf.out', status = 'replace', action='write', iostat=io_error)
           write (33,*) field
        close (unit=33)

      end if


    end do
    close(unit=10)
    close(unit=40)

 else

    print*, 'BEGIN SCF CALCULATION'
    write(30,*) 'BEGIN SCF CALCULATION'
    write(30,*)

    call Find_SCF_Saddle_Point(density,field)
    call Calculate_Densities(field,density)

  end if

  call Calculate_Energy

  if (qddft) then
    write(30,*) '   Final values after ', number_of_timesteps, 'steps'
    print*, ' Final values after ', number_of_timesteps, 'steps'
  else
    write(30,*) '   Final values '
    print*, ' Final values '

  end if
  write(30,*) '   Energy: ', energy, ' C: ', multiblock%rho, ' mu_scf: ', multiblock%mu
  write(30,*)
!
  print*, '   Energy: ', energy, energy/volume, ' C: ', multiblock%rho, ' mu_scf: ', multiblock%mu
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   OUTPUT OF FINAL CONFIGURATION
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      open(unit=60,file='final_dens.dat', status='replace',action='write',iostat=io_error)
      open(unit=61,file='final_field.dat', status='replace',action='write',iostat=io_error)
        write(60,*) "# Energy : ", energy
        write(60,*) 
        write(61,*) "# Energy : ", energy
        write(61,*) 
        dx = sizex/float(gridx)
        dy = sizey/float(gridy)
        dz = sizez/float(gridz)
        do x=1,gridx
           do y=1,gridy
              do z=1,gridz
                 write(60, *) x*dx, '   ' , y*dy, '   ',z*dz, '   ', &
                             (density(x,y,z,i),i=1,monomer_types)
                 write(61, *) x*dx, '   ' , y*dy, '   ',z*dz, '   ', &
                             (field(x,y,z,i),i=1,monomer_types)
              end do
           end do
           write(60, *)
           write(61, *)
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
  subroutine Read_Input_Parameters(file_name)
    
    implicit none

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
        read(20,*) start_time  ! start time
        read(20,*) number_of_timesteps
        read(20,*) timestep
        read(20,*)
        read(20,*) qread  ! qread = .true. : read from file cnf.in
        read(20,*) qinit  ! only for qread=.false.  
            ! 1=cosine, 2=random densities, 3=random fields
        read(20,*) qddft  ! qepd = .true. : EPD simulation; otherweise SCF calculation
        read(20,*)
        read(20,*) wcutoff

      else
        write(*,*) 'Error', io_error , ' while trying to open ', file_name
      end if
!
    close(unit=20)

  end subroutine Read_Input_Parameters
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   SUBROUTINE TO READ AND INITIALIZE MOBILITIES  ! THIS DEPENDS ON SYSTEM !!!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Read_Mobilities(file_name)
    
    use Fourier_fftw3

    implicit none

    character*24, intent(in) :: file_name

    integer, parameter       :: nqmax=1000
    double precision, dimension(nqmax,monomer_types,monomer_types) :: mobility_table

    double precision :: dum, dq, qq, q, small
    integer          :: io_error, nq, iq, x, y, z, iqcut
!
!   ------------ read mobilities 
!
    small = 1.d-6

    open(unit=20,file=file_name,status='old',action='read',iostat=io_error)
      if(io_error ==0) then
        read(20,*) nq
          if (nq.gt.nqmax) then
            stop 'too many entries in mobility file -> adjust nqmax'
          end if
        read(20,*) dq
        read(20,*)
        read(20,*)
        read(20,*)
        do iq = 1, nq
           read(20,*) dum , mobility_table(iq,1,1)  &
                     , mobility_table(iq,1,2),  mobility_table(iq,2,2)
           if (abs(dum -iq*dq) .ge. small ) then
             print*, 'q values in mobility are not consistent', dum, iq*dq
             stop 
           end if
        end do
        mobility_table(:,2,1) = mobility_table(:,1,2)
      else
        write(*,*) 'Error', io_error , ' while trying to open ', file_name
      end if
    close(unit=20)
!
!   --------------- Calculate mobility_qq = Lambda(q)*q^2
!

    do z = 1,gridz
     do y = 1,gridy
      do x = 1,gridx/2+1
         qq = real(-laplace(x,y,z), kind(0.0d0))
         q  = dsqrt(qq)
         iq = int(q/dq)
         if (iq .eq. 0) then
           mobility_qq(x,y,z,:,:) = mobility_table(1,:,:)*q/dq
         else if (iq .ge. nq) then
           mobility_qq(x,y,z,:,:) = mobility_table(nq,:,:)
         else 
           mobility_qq(x,y,z,:,:)  &
              = (   mobility_table(iq,:,:)  *(q  -  iq*dq) &
                  + mobility_table(iq+1,:,:)*((iq+1)*dq-q) ) / dq
         end if
         mobility_qq(x,y,z,:,:) = min(qq*mobility_qq(x,y,z,:,:),wcutoff)
      end do
     end do
    end do
  
  end subroutine Read_Mobilities
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
    implicit none

    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in)  :: field
    double precision, dimension(gridx,gridy,gridz,monomer_types), intent(out) :: density

    double precision, dimension(gridx,gridy,gridz,monomer_types) :: expfield

    integer :: monomer_type,i

    call Calculate_Expfield(field,expfield)
    call Propagate_Multiblock(multiblock,expfield)

    density = 0.0d0
    do i = 1,nblocks
      monomer_type = multiblock%monomer_type(i)
      density(:,:,:,monomer_type) = density(:,:,:,monomer_type) + multiblock%density(:,:,:,i)
    end do

!    density(:,:,:,diblock%monomer_type) = diblock%density

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

    chain_energy = multiblock%energy - dvol*sum(density*field)
    energy = chain_energy + interaction_energy(density)

    return

  end subroutine Calculate_Energy
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!      SUBROUTINE TO CALCULATE FIELDS CORRESPONDING TO GIVEN DENSITIES  
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Get_Conjugate_Fields(density,field,accuracy_scf,itermax_scf,mixing_type,gcfstat)
!
   use Propagate
 
     implicit none

     double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: field
     double precision, dimension(gridx,gridy,gridz,monomer_types), intent(in)    :: density
     integer, intent(out) :: gcfstat
 
     double precision, dimension(gridx,gridy,gridz,monomer_types) :: dfield
     double precision, dimension(gridx,gridy,gridz,monomer_types) :: density0, newdensity

     double precision, intent(in) :: accuracy_scf
     integer, intent(in) :: itermax_scf, mixing_type
 
     double precision :: lambda, accuracy, small
     integer :: iter_status, iter
 
!    accuracy_scf = 0.001
!    itermax_scf = 100000
!    mixing_type = 3

     lambda = 0.5
 
!    small = 1.d-6
!    where (density > small)
!      density0 = density
!    elsewhere
!      density0 = small
!    end where
!
     iter_status = 1
 
     do iter = 1,itermax_scf
 
       call Calculate_Densities(field,newdensity)
       dfield = newdensity - density
       call Check_Accuracy(dfield,accuracy)
!        if ( ( (iter/100)*100) .eq. iter) then
!           print*, iter, accuracy
!        end if

       if (.not.(accuracy > accuracy_scf*0.5) ) then  ! double negation to deal with NaNs
!        PRINT*, 'exit Get_Conjugate_Field after', iter, 'steps '
!        PRINT*, 'Accuracy: ', accuracy
         exit
       end if
       iter_status = min(iter,mixing_dim+1)
       call Mix_Fields(mixing_type,field,dfield,lambda,iter_status)

     end do

     if (.not.(accuracy < accuracy_scf*0.5)) then  ! double negation to deal with NaNs
       if (.not.(accuracy <  accuracy_scf)) then  ! double negation to deal with NaNs
         print*, 'Warning: Conjugate fields not found within ',iter,' iteration steps!'
         print*, 'Accuracy: ', accuracy
         gcfstat = 2 ! iteration not successful
         if (.not.(accuracy <  0.1d0)) then  ! double negation to deal with NaNs
           stop
         end if
       else
         print*, 'Caution: Finding conjugate fields takes long!'
         print*, 'Accuracy after ', iter, 'iteration steps: ', accuracy
         gcfstat = 1 ! iteration slow
       end if
     else
       gcfstat = 0 ! iteration successful
!       print*, 'Conjugate fields found within ',iter,' iteration steps!'
!       print*, 'Accuracy: ', accuracy
     end if

     return

  end subroutine Get_Conjugate_Fields
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!      SUBROUTINE TO FIND SCF SADDLE POINT
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
   subroutine Find_SCF_Saddle_Point(density,field)
!
    use Propagate
!   use Interactions_FloryHuggins
 
     implicit none

     double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: field
     double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: density
 
     double precision, dimension(gridx,gridy,gridz,monomer_types) :: dfield
     double precision, dimension(gridx,gridy,gridz,monomer_types) :: newdensity, newfield
 
     double precision :: lambda, accuracy, accuracy_scf
     integer :: mixing_type, iter_status, iter, iter_max
 
     accuracy_scf = 0.00001
     mixing_type = 3
     iter_max = 1000000
     lambda = 0.5
 
     iter_status = 1
 
     do iter = 1,iter_max
 
       call Calculate_Densities(field,newdensity)
       call Calculate_Mean_Fields(newdensity,newfield)
       dfield = newfield - field
       call Check_Accuracy(dfield,accuracy)
      if ( ( (iter/100)*100) .eq. iter) then
        print*, iter, accuracy
      end if
 
       if (.not.(accuracy > accuracy_scf) ) then  ! double negation to deal with NaNs
         print*, 'SCF saddle point found after ',iter,' iteration steps!'
         print*, 'Accuracy: ', accuracy
         exit
       end if
       iter_status = min(iter,1000)
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
!	   SUBROUTINE TO MIX FIELDS IN SCF ITERATION
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Mix_Fields(mixing_type,field,dfield,lambda,iter_status)
!
   use Iterate
!
     implicit none

     integer, intent(in)               :: mixing_type, iter_status
     double precision, intent(in)      :: lambda
      
     double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: field
     double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: dfield
!
     integer, parameter :: nvar = gridx*gridy*gridz*monomer_types
     double precision, dimension(nvar,0:mixing_dim), save :: var,dvar 

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

!    accuracy=sum(Abs(dfield))*dvol
    accuracy=maxval(Abs(dfield))

    return
  end subroutine check_accuracy
!
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!      DDFT TIME STEP
!          USING FORWARD EULER METHOD
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
   subroutine DDFT_Timestep_Euler (density,field,mobility_qq,timestep,qdmin)
 !
    use Fourier_fftw3

      implicit none

      double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: density
      double precision, dimension(gridx,gridy,gridz,monomer_types), intent(inout) :: field

      double precision, dimension(gridx/2+1,gridy,gridz,monomer_types,monomer_types), intent(in) :: mobility_qq
      double precision, intent(in) :: timestep
      logical, intent(in)          :: qdmin

      double precision, dimension(gridx,gridy,gridz,monomer_types)   :: newfield
      double precision, dimension(gridx,gridy,gridz,monomer_types)   :: mulocal
      double precision, dimension(gridx,gridy,gridz,monomer_types)   :: ddensity

      double complex, dimension(gridx/2+1,gridy,gridz,monomer_types) :: mu_q
      double complex, dimension(gridx/2+1,gridy,gridz,monomer_types) :: ddensity_q

      double precision, dimension(gridx,gridy,gridz)   :: dummy
      double complex, dimension(gridx/2+1,gridy,gridz) :: dummyfourier
!
      double precision :: dsmall
      integer :: x,y,z,n,monomer_type,type1,type2, gcfstat
      double precision :: offset(monomer_types), dum
!
!     Local potential in Real space
!
      call Get_Conjugate_Fields(density,field,1.d-3,2000,3,gcfstat)
      call Calculate_Mean_Fields(density,newfield)
      mulocal = newfield - field
!      do monomer_type = 1,monomer_types
!        offset(monomer_type) = sum(mulocal(:,:,:,monomer_type))/float(ngrid)
!        mulocal(:,:,:,monomer_type) = mulocal(:,:,:,monomer_type) &
!             - offset(monomer_type)
!      end do
!
!     Local potential in Fourier space
!
      do monomer_type = 1, monomer_types
        dummy = mulocal(:,:,:,monomer_type)
        call Real2Fourier(dummy,dummyfourier)
!        dummyfourier = - laplace*dummyfourier
        mu_q(:,:,:,monomer_type) = dummyfourier
      end do
!
!     Density derivative in Fourier space
!
      ddensity_q = 0.0d0
      do type1 = 1, monomer_types
        do type2 = 1, monomer_types
          ddensity_q(:,:,:,type1) =  ddensity_q(:,:,:,type1) &
            - mobility_qq(:,:,:,type2,type1)*mu_q(:,:,:,type2)
         !  No q^2 because it is already contained in mobility_qq
        end do
      end do
!
!     Density derivative in Real space
!
      do monomer_type = 1, monomer_types
        dummyfourier = ddensity_q(:,:,:,monomer_type)
        call Fourier2Real(dummy,dummyfourier)
        ddensity(:,:,:,monomer_type) = dummy
!        offset(monomer_type) = sum(ddensity(:,:,:,monomer_type))/ngrid
!        ddensity(:,:,:,monomer_type) = ddensity(:,:,:,monomer_type) &
!            + offset(monomer_type)
      end do
!
!     Euler step
!
      dsmall = 1.d-4
      density = density + timestep*ddensity
      if (qdmin) then
        dum=minval(density)
        if (dum .le. dsmall) then
          density = max(density, dsmall)
        end if
      end if
      call Get_Conjugate_Fields(density,field,1.d-3,2000,3,gcfstat)
      if (qdmin .and. (gcfstat.ge.1)) then
        call Calculate_Densities(field,density)
      end if
!
     return
 
   end subroutine DDFT_Timestep_Euler

!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	     T H E    E N D 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end program main_multiblock
