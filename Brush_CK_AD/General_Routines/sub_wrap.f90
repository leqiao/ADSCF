!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    INITIALIZATIONS
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Initialize

    use Global
    use mtmod
    use RealSpace
    use Adaptive
!
    implicit none
!

!
! ------- INITIALIZATION OF RNG --------
!
    call sgrnd(2007)
!             
!
! ------ OTHER INITIALIZATIONS ---------
!
    call calculate_Volume
    call generate_mesh
!
!--------------------------------------
    return

  end subroutine initialize
