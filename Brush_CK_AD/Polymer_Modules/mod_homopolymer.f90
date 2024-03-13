!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           ROUTINES FOR HOMOPOLYMERS 
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!

module Molecule_Homopolymer
!
  use Global
  use Adaptive, only: integrator
  use RealSpace
!
  implicit none
!
  private
!

  type homopolymer
    integer          :: label
    integer          :: monomer_type
    double precision :: length
    integer          :: ensemble       ! canonical:1, grand canonical:2
    double precision :: rho            ! total polymer concentration 
    double precision :: mu             ! chemical potential of polymer
    double precision :: energy         ! free energy in SCF field
    integer          :: segm
    double precision :: Qsum
    double precision, dimension(gridx,gridy,0:gridz)    :: distribution
    double precision, dimension(gridx,gridy,0:gridz)    :: density
    double precision, dimension(gridx,gridy,0:gridz)    :: delta_rep_list
    double precision, dimension(:), allocatable         :: dss

  end type homopolymer
!
  public :: homopolymer
  public :: Initialize_Homopolymer
  public :: Propagate_Homopolymer
!
contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!           INITIALIZES HOMOPOLYMER
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!

   function gaus(z) result(otpt)
    
    double precision, intent (in)   :: z   
    double precision                :: sig,aver,otpt  
    
    aver=gra_pt
    sig=aver/5.0d0
   
    if ( (abs(z-aver) .lt. 5.0d0*sig) .and. (z .ge. 1) ) then
        otpt=exp(-(z-aver)**2.0d0/(2.d0*sig**2.0d0)) !Normalization not needed
    else
        otpt=0.0d0
    end if 


  end function


  subroutine Initialize_Homopolymer(polymer,label,monomer_type, &
                                    ensemble,rho)
!
   integer, intent(in)          :: label, ensemble, monomer_type
   double precision, intent(in) ::  rho
   integer                      :: z , i, iAllocStatus
   double precision             :: temp
!
   type(homopolymer), intent(out):: polymer
!
   character*15                  :: string_ensemble
!


   allocate(polymer%dss(0:segments),stat=iAllocStatus)

   if (ensemble.eq.1) then
      string_ensemble = 'canonical'
   else
      print*, 'Invalid Ensemble'
   end if
!
   polymer%label        = label
   polymer%length       = 1.0d0
   polymer%monomer_type = monomer_type
   polymer%ensemble     = ensemble
   polymer%rho          = rho 

   polymer%segm     = segments
   
   polymer%dss(0)=0.0d0
   temp=0.0d0

   !Initialising ds discretization
   do i=1,polymer%segm
       if (ds_scheme .eq. 0) then
             polymer%dss(i)=dble(1.0d0/segments)     
       else if (ds_scheme .eq. 1) then
            polymer%dss(i)=1.0 * 0.5d0 * (DCOS(pi * (dble(i-1)) /  dble(segments)) - DCOS(pi * (dble(i)) /  dble(segments)))
       else if (ds_scheme .eq. 2) then
            polymer%dss(i)=1.0d0*(DCOS(pi*(dble(i-1)) /(2.0d0*dble(segments)) )-DCOS(pi *(dble(i))/(2.0d0*dble(segments))))
       end if
        temp=temp+polymer%dss(i)
   end do

   do i=1,polymer%segm
      polymer%dss(i)=polymer%dss(i)/temp
   end do
   print*, "ds sum ", sum(polymer%dss(:))
!


   
   !Initialising grafting point
   polymer%delta_rep_list=0.0d0
   if (delta_rep .eq. 1) then
       do z=0,gridz
        polymer%delta_rep_list(:,:,z)=gaus(zr(z))
       end do
       if ( count(polymer%delta_rep_list/=0) .eq. 0) then
          stat_job=1
       end if
   end if

   if ( ( count(polymer%delta_rep_list/=0) .le. 1 ) .or. ( delta_rep .eq. 0) ) then
       
       do i=1,gridz
        if (abs(zr(i)-gra_pt) .lt. 0.000001d0 ) then
            polymer%delta_rep_list(:,:,i)=1.0d0 
        end if
       end do
   end if 
   
   if (delta_rep .eq. 2) then 
    polymer%delta_rep_list=1.0d0
   end if
   
   if (delta_rep .ne. 2) then 
       call integrator(polymer%delta_rep_list(:,:,:),temp)   
       polymer%delta_rep_list=polymer%delta_rep_list/temp
       call integrator(polymer%delta_rep_list(:,:,:),temp)
       print*, "Norm const. of delta rep : " , temp
   end if



    !Output test quantities
    open(unit=72,file='test_delta.dat', status='replace',action='write')
        do z=0,gridz
          write(72, *) zr(z) , '   ' , polymer%delta_rep_list(1,1,z)
        end do
    close(unit=72)

    open(unit=72,file='test_ds.dat', status='replace',action='write')
    do i=0,segments
      write(72, *) i , '   ' , polymer%dss(i) 
    end do
    close(unit=72)



   print*, 'Homopolymer ', polymer%label, ' is initialized (', &
        &string_ensemble, ' ensemble).'


  end subroutine Initialize_Homopolymer
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   CALCULATE PROPAGATOR AND DISTRIBUTION FOR A HOMOPOLYMER CHAIN OF LENGTH segments*ds
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Propagate_Homopolymer(polymer,field)
!
    type(homopolymer), intent(inout)  :: polymer
    double precision, dimension(gridx,gridy,0:gridz,monomer_types), intent(in)  :: field

    double precision, dimension(:,:,:,:), allocatable :: qq
    double precision, dimension(:,:,:,:), allocatable :: qqdag
    double precision, dimension(gridx,gridy,0:gridz)  :: prop

    double precision :: total_rho ,temp
!
    integer :: iAllocStatus
    integer :: s,mono ,sd , i


!
    allocate(qq(gridx,gridy,0:gridz,0:polymer%segm),stat=iAllocStatus)
    allocate(qqdag(gridx,gridy,0:gridz,0:polymer%segm),stat=iAllocStatus)

!
    


    mono = polymer%segm
    prop=1.0d0 ! FREE END INITIAL CONDITION BY CONVENTION
    qq(:,:,:,0)=prop
    
    
    do s=1,mono
       call Propagate_Step(prop,field(:,:,:,1),polymer%dss(s))
       qq(:,:,:,s)=prop
    end do
   

    prop=polymer%delta_rep_list
    qqdag(:,:,:,0)=prop
    

    do s = 1,mono
       call Propagate_Step(prop,field(:,:,:,1),polymer%dss(mono+1-s)) !since we are going the other direction
       qqdag(:,:,:,s)=prop
    end do



    call integrator(qqdag(:,:,:,0)*qq(:,:,:,mono),polymer%Qsum)


    
     
    polymer%distribution=polymer%dss(0)*qq(:,:,:,0)*qqdag(:,:,:,mono)
    do s=1,mono-1
       polymer%distribution = polymer%distribution &
          + (polymer%dss(s)+polymer%dss(s+1))*qq(:,:,:,s)*qqdag(:,:,:,mono-s)
    end do
    polymer%distribution = polymer%distribution &
          + polymer%dss(mono)*qq(:,:,:,mono)*qqdag(:,:,:,0)

    polymer%distribution=0.5d0*polymer%distribution



    
    if (polymer%ensemble.eq.1) then
      polymer%density = polymer%rho*volume * polymer%distribution/polymer%Qsum !rho here has meaning of chain number density
      polymer%density(:,:,0)=0.0d0
      polymer%mu = log(polymer%rho*volume/polymer%Qsum)
      polymer%energy = polymer%rho*volume * (polymer%mu - 1.0d0)
    else
      print*, 'Invalid ensemble'
      stop
    end if
!
    deallocate(qq)
    deallocate(qqdag)
!
    return

  end subroutine Propagate_Homopolymer
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module Molecule_Homopolymer
