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
!
    implicit none
!
! ------ FFTW INITIALIZATION -----------
!
!    call Initialize_Fourier
!
! -------- TEST OF FOURIER TRANSFORM
!
!   call fourier_test
!
! ------- INITIALIZATION OF RNG --------
!
    call sgrnd(2007)
!  
! -------- LAPLACE OPERATOR ------------
!
!    call Calculate_Laplace           
!
! ------ OTHER INITIALIZATIONS ---------
!
    call calculate_Volume
    call generate_mesh
!
!--------------------------------------
    return

  end subroutine initialize
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!	    CLOSING ROUTINES
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine Finish

    use Global
!
    implicit none
!
! ------ FFTW CLOSING -----------
!

  end subroutine finish
