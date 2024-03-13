!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!               SOLVE DIFFUSION EQUATIONS IN REAL SPACE
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module RealSpace
!
use global
implicit none
private
public Propagate_Step, generate_mesh, integrator, calculate_Volume
!
contains
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!          SOLVE DIFFUSION EQUATION
!          CALCULATE PROPAGATOR BY A STEP q(s+ds)
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine Propagate_Step(prop,field)
    real*8, dimension(gridx,gridy,0:gridz), intent(in)       :: field
    real*8, dimension(gridx,gridy,0:gridz), intent(inout)    :: prop
    real*8, dimension(0:gridx+1,0:gridy+1,0:gridz)           :: prop_0 ! temp array for propogator
    real*8, dimension(gridx)                                 :: rhand_inx, rhand_outx
    real*8, dimension(gridy)                                 :: rhand_iny, rhand_outy
    real*8, dimension(0:gridz)                               :: rhand_inz, rhand_outz
    real*8, dimension(gridx)                                 :: ax, cx
    real*8, dimension(gridy)                                 :: ay, cy
    real*8, dimension(0:gridz)                               :: az, cz, bz
    real*8, dimension(gridx,gridy,0:gridz)                   :: coe_x,coe_y,coe_z ! daigonal coefficient containing field 
    real*8                                                   :: bx,by      
    integer                                                  :: x, y, z 
    !!COEFFICIENT IN Z-DIRECTION
    !SETUP DIRICHLIT BOYUNDARY BY MAKEING P=0
    prop_0=0.0d0
    prop_0(1:gridx,1:gridy,:)=prop
    !SETUP PERIODIC BOUNDARY BY ADDING GOHST SITES
    prop_0(0,1:gridy,:)=prop(gridx,:,:)
    prop_0(gridx+1,1:gridy,:)=prop(1,:,:)
    prop_0(1:gridx,0,:)=prop(:,gridy,:)
    prop_0(1:gridx,gridy+1,:)=prop(:,1,:)
    
!  COEFFICIENT ARRAY IN X-DIRECTION
    ax = -ds/dx/dx/2.d0
    bx =  ds/dx/dx
    cx = -ds/dx/dx/2.d0
    do z = 1,gridz-1
       do y = 1,gridy
          do x = 1,gridx
             coe_x(x,y,z) = 1.d0 + bx
          end do
       end do
    end do
!
!  COEFFICIENT ARRAY IN Y-DIRECTION
    ay = -ds/dy/dy/2.d0
    by =  ds/dy/dy
    cy = -ds/dy/dy/2.d0
    do z = 1,gridz-1
       do y = 1,gridy
          do x = 1,gridx
             coe_y(x,y,z) = 1.d0 + by
          end do
       end do
    end do
!
!   COEFFICIENT IN Z-DIRECTION
    az = 0.d0
    cz = 0.d0
    do z = 1,gridz-1
        az(z) = -ds * 1.d0/dz(z)/(dz(z)+dz(z+1))
        bz(z) =  ds * 1.d0/dz(z)/dz(z+1) ! Note dz(0) is not defined, here we approximate h(0)=h(1)
        cz(z) = -ds * 1.d0/dz(z+1)/(dz(z)+dz(z+1))
    end do
    coe_z=1.d0
    do z = 1,gridz-1
        do y = 1,gridy
            do x = 1,gridx
                coe_z(x,y,z) = 1.d0 + bz(z) + ds*field(x,y,z)/2.0
            end do
       end do
    end do
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!               SOLVE IN X-DIRECTION
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   do z = 1,gridz-1
      do y = 1,gridy
         do x = 1,gridx
            rhand_inx(x) = -ax(x)*prop_0(x-1,y,z) - 2.d0*ay(y)*prop_0(x,y-1,z) -2.d0*az(z)* prop_0(x,y,z-1) &
                           -cx(x)*prop_0(x+1,y,z) - 2.d0*cy(y)*prop_0(x,y+1,z) -2.d0*cz(z)* prop_0(x,y,z+1) &
                           +(1.0d0 - bx - by*2.d0 - bz(z)*2.d0 - ds * field(x,y,z)) * prop_0(x,y,z)
         end do
         Call cyclicTriDiag(ax, coe_x(:,y,z), cx, rhand_inx, rhand_outx, gridx)
         do x = 1,gridx
            prop(x,y,z) = rhand_outx(x)
         end do
      end do
   end do
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!               SOLVE IN Y-DIRECTION
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   do z = 1,gridz-1
      do x = 1,gridx
         do y = 1,gridy
            rhand_iny(y) = ay(y) * prop_0(x,y-1,z) + prop(x,y,z) + by*prop_0(x,y,z) &
                + cy(y) *prop_0(x,y+1,z)
         end do
         Call cyclicTriDiag(ay, coe_y(x,:,z), cy, rhand_iny, rhand_outy, gridy)
         do y = 1,gridy
            prop(x,y,z) = rhand_outy(y)
         end do
      end do
    end do
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!             SOLVE IN Z-DIRECTION
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    do x = 1,gridx
       do y = 1,gridy
            rhand_inz(0) = 0.d0
            rhand_inz(gridz) = 0.d0
            do z = 1,gridz-1
                rhand_inz(z) = az(z) *prop_0(x,y,z-1) + prop(x,y,z) &
                            + (bz(z) + ds*field(x,y,z)/2.d0)*prop_0(x,y,z) &
                             + cz(z)*prop_0(x,y,z+1)
            end do
          Call tridiag(az, coe_z(x,y,:), cz, rhand_inz, rhand_outz, gridz+1)
          do z = 0,gridz
             prop(x,y,z) = rhand_outz(z)
          end do
       end do
    end do
   end subroutine Propagate_Step
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!       SOLVE EQUATION WITH PERIODIC BOUNDARIES
!
!         m u = r, where a, b and c are the three main diagonals of matrix  
!         m(n,n), the other terms are 0. r is the right side vector.        
!         m(1,n) = a(1) (or beta from NR); and m(n,1) = c(n) (or alpha from NR)
!
   subroutine cyclicTriDiag(a, b, c, r, x, N)
     integer, intent(in)                  :: N
     double precision, intent(in)         :: a(N), b(N), c(N), r(N)
     double precision, intent(out)        :: x(N)
!
     double precision, dimension(N) :: bb, u, y, q
     integer :: i
     double precision  :: fact, gamma
     if (N <3) then
         x=r(1)/(b(1)+a(1)+c(1))
     else
         gamma     = -b(1)   ! gamma can be arbit; -b1 is recommended, concustruct Matrix B with bb at diaggonals
         bb(1)     =  b(1) - gamma
         bb(N)     =  b(N) - a(1)*c(N)/gamma
         bb(2:N-1) =  b(2:N-1)
         call tridiag(a, bb, c, r, y, N)          ! Solve By = r

         u(1)     = gamma                      ! Set up vector u
         u(N)     = c(N)   !a(1)
         u(2:N-1) = 0.d0
         call tridiag(a, bb, c, u, q, N)           ! Solve Bq = u
! Form v.x/(1 + v.z)
!
         fact = (y(1) + a(1)*y(N)/gamma)/(1.0 + q(1) + a(1)*q(N)/gamma)
         x    = y - fact*q   ! get solution vector z
     end if
!
   end subroutine cyclicTriDiag
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!             SOLVE EQUATION WITH DIRICHLET BOUNDARIES
!  
!       m u = r, where a, b and c are the three main diagonals of matrix  !
!       m(n,n), the other terms are 0. r is the right side vector.        !
   subroutine tridiag(a, b, c, r, u, N)

     integer, intent(in)                  :: N
     double precision, intent(in)         :: a(N), b(N), c(N), r(N)
     double precision, intent(out)        :: u(N)
     double precision, dimension(N)       :: gam
     integer :: j
     double precision  :: bet

     bet = b(1)
     if (bet == 0.0) then
         write(*,*) "TriDiag Error"
         stop
     endif

     u(1) = r(1)/bet

     do j = 2,N
        gam(j) = c(j-1)/bet
        bet    = b(j) - a(j)*gam(j)
        if (bet == 0.0) then
            write(*,*) "TriDiag Error"
            stop
        endif
        u(j) = ( r(j)-a(j)*u(j-1) )/bet
     end do

     do j = N-1, 1, -1
        u(j) = u(j) - gam(j+1)*u(j+1)
     end do

   end subroutine tridiag
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!              INITIALIZE: CALCULATE VOLUME PARAMETERS
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine generate_mesh
!------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------!
integer             :: i, grid_surf

double precision            :: dz_surf, dz_bulk, z_surf
!------------------------------------------------------------------------------------------------------!
dx= sizex/dble(gridx)
dy=sizey/dble(gridy)
if (mesh_type.eq.0) then !uniform
    dz = sizez / dble(gridz)
    dz(0) = 0.d0
    do i = 0, gridz
        zr(i) = dble(i) * dz(i)
    enddo
elseif (mesh_type.eq.1) then ! nonuniform Cos
    dz(0)=0.d0
    zr(0)=0.d0
    do i = 1, gridz
        zr(i) = sizez * 0.5d0 * (1.d0 - dcos(pi * (dble(i)) /  dble(gridz)))         !symmetric scheme
        dz(i) = zr(i) - zr(i-1)
    enddo
elseif (mesh_type.eq.2) then ! nonuniform step function
    print*, 'use discretization from step function'
    dz(0)=0.d0
    zr(0)=0.d0
    grid_surf=20
    z_surf=delta_z
    dz_surf=z_surf/dble(grid_surf)
    dz_bulk=(sizez/2.d0-z_surf)/(dble(gridz)/2.d0-dble(grid_surf)) !grid size in bulk
    do i = 1, gridz
       if (i .le. grid_surf) then
          zr(i) = zr(i-1) + dz_surf
       elseif ((i .gt. z_surf) .and. (i .le. gridz-grid_surf)) then
          zr(i)=zr(i-1)+ dz_bulk
       elseif (i .gt. gridz- grid_surf) then
          zr(i) = zr(i-1) + dz_surf
       end if
        dz(i) = zr(i) - zr(i-1)
    enddo
    print*, 'difference dz',zr(gridz)-sizez 
endif
!------------------------------------------------------------------------------------------------------!
end subroutine generate_mesh

subroutine integrator(integrand, integral)
    double precision, dimension(gridx,gridy,0:gridz), intent(in)  :: integrand
    integer            :: i ! 
    double precision, intent(out)               :: integral
    integral=0.d0
    if (int_scheme == 0) then ! recangle rule
        do i=0, gridz-1
        integral = integral + sum(dx*dy*dz(i)*integrand(:,:,i))
        end do
    else if (int_scheme == 1) then ! trapezoidal rule
        do i = 1, gridz
        integral = integral + sum(dx*dy*dz(i)*(integrand(:,:,i-1) + integrand(:,:,i))/2.d0)
        end do
     end if
end subroutine integrator
!------------------------------------------------------------------------------!
subroutine calculate_Volume
    volume = sizex*sizey*sizez
end subroutine calculate_Volume

end module RealSpace
