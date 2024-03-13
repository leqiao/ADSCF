module Adaptive
!
use Global
implicit none
private
public :: generate_mesh
public :: integrator
!
contains

function step(x) result(otpt)
    
    integer, intent (in)   :: x   
    double precision       :: otpt  
    
   
    if ( x .eq. 0 ) then
        otpt=0.5d0
    elseif ( x .gt. 0 ) then
        otpt=1.0d0
    elseif ( x .lt. 0 ) then
        otpt=0.0d0
    end if 

end function

subroutine generate_mesh
!------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------!
integer             :: i
double precision            :: dz_surf, dz_bulk, temp
double precision            :: maax,nrm ,b ,aa,bb,cc ! MG VARIABLES FOR .EQ. 2,3 CASE  

!------------------------------------------------------------------------------------------------------!

if (mesh_type.eq.0) then !uniform
    dz = sizez / dble(gridz)
    
    do i = 0, gridz
        zr(i) = dble(i) * dz(i)
    enddo

    if ((abs(DMOD(gra_pt,dz(1))) .gt. 0.000001d0) .or. ( dz(1) .gt. gra_pt )) then
    write(*,*) 'CANT PUT GRAFTING POINT AT THE POINT SPECIFIED'
    call exit(1)
    end if

elseif (mesh_type.eq.1) then ! nonuniform Cos // Not recommended to use with grafting point as it changes the point
    dz(0)=0.d0
    zr(0)=0.d0
    temp=10**9
    
    do i = 1, gridz
        zr(i) = sizez * 0.5d0 * (1.d0 - DCOS(pi * (dble(i)) /  dble(gridz)))         !symmetric scheme
    end do


elseif (mesh_type.eq.2) then ! dz is constant for points 2*eps0 from each wall and is set to eps 
  dz(0)=0.d0
  zr(0)=0.d0

  if ((abs(DMOD(gra_pt,eps)) .gt. 0.000001d0) .or. ( eps .gt. gra_pt )) then
    write(*,*) 'CANT PUT GRAFTING POINT AT THE POINT SPECIFIED'
    call exit(1)
  end if
  aa=eps

  bb=(4.0d0* eps* eps0 + eps* gridz - 2.0d0* sizez)/(4.0d0* eps0 - gridz)

  print*, aa ,bb , gridz/2, 'aa and bb and gridz/2 '
  if (bb .le. aa) then
    write(*,*) 'WARNING IN MESH GENERATION BB IS SMALLER OR EQUAL TO AA'
    call exit(1)
  end if

  

  do i=1,gridz
    dz(i)=(aa)*step(2*eps0+1-i)&
          +(aa + (bb - aa)/(gridz/2.0d0 - 2.0d0*eps0) * (i - 2*eps0))*step(i - 2*eps0-1)*step(gridz/2 -i) &
          +(bb + (aa - bb)/(gridz/2-2* eps0)*(i - gridz/2))*step(i - gridz/2)*step(gridz-2*eps0 -i) &
          + (aa)*step(i-gridz+2*eps0)
    zr(i)=zr(i-1)+dz(i)
  end do

elseif (mesh_type.eq.3) then ! dz is constant for points 2*eps0 on left wall and then cosine after that

  if ((abs(DMOD(gra_pt,eps)) .gt. 0.000001d0) .or. ( eps .gt. gra_pt )) then
    write(*,*) 'CANT PUT GRAFTING POINT AT THE POINT SPECIFIED'
    call exit(1)
  end if
  aa=eps
    
  temp=0.0d0
  do i=2*eps0,gridz
    temp=temp+DSIN(pi*0.5d0*(DBLE(i)-2.0d0*DBLE(eps0))/(DBLE(gridz)-2.0d0*DBLE(eps0)))
  end do
  bb=(sizez - aa*DBLE(gridz))/temp

                     
  dz(0)=0
  do i=1,gridz
    dz(i)=(aa)*step(2*eps0-i)&
        +step(i-2*eps0)*(aa+bb*DSIN(pi*0.5d0*(DBLE(i)-2.0d0*DBLE(eps0))/(DBLE(gridz)-2.0d0*DBLE(eps0))))
    zr(i)=zr(i-1)+dz(i)
  end do
endif


open(unit=75,file='test_dz.dat', status='replace',action='write')
do i=1,gridz
  write(75, *) i, dz(i) , zr(i)
end do
 close(unit=75)




!------------------------------------------------------------------------------------------------------!
end subroutine generate_mesh
!----------------------------------------------------------------------------------------------------------!


subroutine integrator(integrand, integral)
    double precision, dimension(gridx,gridy,0:gridz), intent(in)  :: integrand
    integer            :: i
    double precision, intent(out)               :: integral
    double precision :: dx, dy
    dx=sizex/dble(gridx)
    dy=sizey/dble(gridy)
    integral=0.d0
    if (int_scheme == 0) then ! rectangle rule
        do i=1, gridz-1
        integral = integral + sum(dx*dy*dz(i)*integrand(:,:,i))
        end do
    else if (int_scheme == 1) then ! trapezoidal rule
        do i = 1, gridz-1
        integral = integral + sum(dx*dy*dz(i)*(integrand(:,:,i-1) + integrand(:,:,i))/2.d0)

        end do
     end if
    
end subroutine integrator

end module Adaptive

