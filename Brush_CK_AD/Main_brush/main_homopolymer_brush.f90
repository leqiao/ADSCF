!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    MAIN PROGRAM (HOMOPOLYMER POLYMER BRUSH MELT)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
program main_homopolymer_brush
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    DECLARATIONS 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  use Global
  use Molecule_Homopolymer
  use Interactions_Virial_Expansion
  use RealSpace
  use mtmod
  use Adaptive
!
  implicit none
!
  type(homopolymer) :: polymer
!
  double precision    ::  phibar    ! dimensionless polymer concentration

  double precision    ::  energy , field_energy ,chain_energy

  integer             ::  itermax_scf    ! maximum number of iterations
  integer             ::  iterations     ! actual number of iterations
  double precision    ::  accuracy_scf   ! desired accuracy in iterations
  double precision    ::  accuracy       ! actual accuracy
  double precision    ::  lambda         ! mixing parameters
  integer             ::  mixing_type    ! type of mixing
!      -------------  1: simple mixing, 2: lambda mixing, 3: anderson mixing ----------

  integer             ::  io_error,x,y,z,iter,i , qread 
  
  double precision, dimension(gridx,gridy,0:gridz,monomer_types) :: density
  double precision, dimension(gridx,gridy,0:gridz,monomer_types) :: newdensity
  double precision, dimension(gridx,gridy,0:gridz,monomer_types) :: field
  double precision, dimension(gridx,gridy,0:gridz,monomer_types) :: newfield, dfield

  double precision, dimension(gridx,gridy,0:gridz,monomer_types) :: saddledensity
  double precision, dimension(gridx,gridy,0:gridz,monomer_types) :: saddlefield


  character*100    :: file_name,output_str0
  integer          :: polymer_label, monomer_type
  integer          :: polymer_ensemble   ! canonical:1, grand canonical:2
  integer          :: every_out

  character(100) :: arg_lz,x1 ! box z size

  real :: start_time, finish_time
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	     INITIALIZATIONS
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



stat_job=0

file_name='input_general'
call read_input_parameters(file_name)


! ------ GENERAL INITIALIZATIONS 
!
  call Initialize

!
! ------ INITIALIZE HOMOPOLYMER COPOLYMER
!
   
    
    polymer_label = 1
    polymer_ensemble  = 1  ! canonical
    monomer_type=1


    call Initialize_Homopolymer(polymer,polymer_label,&
                            monomer_type,polymer_ensemble,&
                            phibar)   



  call Initialize_Interactions


  if (qread.eq.1) then
!
    open (unit=33, file='cnf_densities.in', status='old', action='read', iostat=io_error)
    if(io_error /=0) then                   
      write(*,*) 'Error', io_error , ' while trying to open cnf_densities.in'
    end if
    read (33,*) density
    close (unit=33)

  else if ( qread .eq. -1) then
  do x = 1,gridx
     do y = 1,gridy
        do z = 1,gridz-1
        !SST approximation to the init profile
        density(x,y,z,1) = 1.0/ww*( (3.0*pi/8.0*ww*phibar*sizez)**(2.0/3.0)-(pi/4.0*zr(z))**(2d0) )
        
        where (density .lt. 0.0d0)
            density=0
        endwhere
        end do
     enddo
  end do
end if

density(:,:,0,:)=0.0d0 !dirichlet boundary
density(:,:,gridz,:)=0.0d0

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!

 open(unit=72,file='init_test_densities.dat', status='replace',action='write',iostat=io_error)
do z=0,gridz
   do y=1,gridy
     do x=1,gridx
      write(72, *) z, ' ' , zr(z) , '   ' , density(x,y,z,1)
      
     end do
   end do
end do
 close(unit=72)




call Calculate_Mean_Ext_Fields(density,field)
newdensity = density
newfield   = field
 


  call system("mkdir OUTPUT_DENS")

  call cpu_time(start_time)


! -------- ITERATION LOOP START
  do iter = 1,itermax_scf
    call Calculate_Densities(field,newdensity)
    call Check_Accuracy(density,newdensity,accuracy)
    call Calculate_Energy                                   ! FOR TESTING
    print*, "iter", iter," energy " ,energy ,  " accuracy", accuracy  ! FOR TESTING
    if ( .not.(accuracy>accuracy_scf) .and. (iter > 6) ) then     !double negation to deal with NaNs
      exit
    end if
    call Calculate_Mean_Ext_Fields(newdensity,newfield)
    dfield = newfield - field

    Call Mix_Fields(mixing_type,field,dfield,lambda,iter)
    density=newdensity

 if (((iter/every_out)*every_out).eq. iter) then                         ! FOR TESTING

       write (x1,*) iter
      output_str0='./OUTPUT_DENS/test_densities_iter_'//TRIM(ADJUSTL(x1))//'.out'
      open(unit=72,file=output_str0, status='replace',action='write',iostat=io_error)
        do z=1,gridz
           do y=1,gridy
             do x=1,gridx
              write(72, *) z, ' ' , zr(z) , '   ' , density(x,y,z,1)
              
             end do
           end do
        end do
       close(unit=72)
    
     output_str0='./OUTPUT_DENS/cnf_densities_iter_'//TRIM(ADJUSTL(x1))//'.out'
     open(unit=72,file=output_str0, status='replace',action='write',iostat=io_error)
            write (33,*) density
      close (unit=33)
 end if 

  end do
! -------- ITERATION LOOP END



  iterations = iter
  saddledensity = newdensity
  saddlefield   = newfield

  call Calculate_Energy


  if (.not.(accuracy < accuracy_scf)  ) then  ! DOUBLE NEGATION DEALS WITH NaNs
    print*, 'Caution: SCF saddle point not found within ',iterations,' iteration steps!'
    stat_job=2
  end if


!
  print*, '   Saddle point found after', iterations, 'iterations'
  stat_job=0
  print*, '   Accuracy: ', sqrt(accuracy), ' Energy: ', energy, ' C: ', polymer%rho,  &
          ' mu_scf: ', polymer%mu
!
  call cpu_time(finish_time)
  print*,'finish time',(finish_time-start_time)/60.d0, 'min'

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	   OUTPUT OF SADDLE POINT RESULTS
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!

      if ((iterations .ge. itermax_scf) .or. (accuracy .gt. accuracy_scf) .or. (any(density(:,:,:,1) < 0.0))) then
         stat_job=2 
      end if 


      open(unit=63,file='details.dat', status='replace',action='write',iostat=io_error)
      write(63,"(a,F10.6)") '# 1. Status 2. Accuracy 3. Energy 4. Concentration 5. Mu 6.No Iter 7. Time 8.Mesh type '//&
                '9. LamdhaA 10.Int_ene 11.Ex_ene 12.Polymer ene 13. Delta 14.Qsum 15.Field energy 16.Ds type' 
      write(63,*)   stat_job  , accuracy  , energy  , polymer%rho , polymer%mu , iter , (finish_time-start_time)/60.d0 ,&
                                                       mesh_type , LambdaA ,interaction_energy(saddledensity), &
                                                       external_energy(density) , polymer%energy , delta_rep , polymer%Qsum, &
                                                       field_energy, ds_scheme
      close(unit=63) 



      
      open(unit=72,file='test_densities.dat', status='replace',action='write',iostat=io_error)
        do z=1,gridz
           do y=1,gridy
             do x=1,gridx
              write(72, *) z, ' ' , zr(z) , '   ' , saddledensity(x,y,z,1)
              
             end do
           end do
        end do
       close(unit=72)

      open(unit=72,file='test_field.dat', status='replace',action='write',iostat=io_error)
        do z=1,gridz
           do y=1,gridy
             do x=1,gridx
              write(72, *) z, ' ' , zr(z) , '   ' , saddlefield(x,y,z,1)
              
             end do
           end do
        end do
       close(unit=72)

      open(unit=33, file='cnf_densities.out', status = 'replace', action='write', iostat=io_error)
        write (33,*) saddledensity
      close (unit=33)


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
        read(20,*) phibar
        read(20,*) segments
        read(20,*) ds_scheme
        read(20,*) itermax_scf
        read(20,*) accuracy_scf
        read(20,*) mixing_type 
        read(20,*) lambda
        read(20,*) int_scheme
        read(20,*) mesh_type
        read(20,*) delta_rep 
        read(20,*) gra_pt
        read(20,*) eps 
        read(20,*) eps0
        read(20,*) LambdaA
        read(20,*) delta_z
        read(20,*) wall_scheme
        read(20,*) qread 
        read(20,*) every_out
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
   use RealSpace
    double precision, dimension(gridx,gridy,0:gridz,monomer_types), intent(in)  :: field
    double precision, dimension(gridx,gridy,0:gridz,monomer_types), intent(out) :: density

    call Propagate_Homopolymer(polymer,field)
    density = 0.0d0
    density(:,:,:,1) = polymer%density

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
    double precision              ::  saddle_energy

    call integrator(sum(density*field,dim=4), field_energy)
    chain_energy = polymer%energy - field_energy  

    energy = chain_energy + interaction_energy(density) + external_energy(density)

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
      
     double precision, dimension(gridx,gridy,0:gridz,monomer_types), intent(inout) :: field
     double precision, dimension(gridx,gridy,0:gridz,monomer_types), intent(inout) :: dfield
!
     integer, parameter :: nvar = gridx*gridy*gridz*monomer_types
     double precision, dimension(nvar,0:mixing_dim) :: var,dvar 

     integer :: x,y,z,n,monomer_type

      n=0
      do monomer_type = 1,monomer_types 
        do z = 1,gridz-1
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
        do z = 1,gridz-1
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
    double precision, dimension(gridx,gridy,0:gridz,monomer_types), intent(in) ::    density
    double precision, dimension(gridx,gridy,0:gridz,monomer_types), intent(in) :: newdensity
    double precision, intent(out)       :: accuracy

    accuracy=sum((density-newdensity)*(density-newdensity))

    return
  end subroutine check_accuracy


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!	     T H E    E N D 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end program main_homopolymer_brush
