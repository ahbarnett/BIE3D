MODULE dspline

IMPLICIT NONE

! Use predictor-corrector schemes 

CONTAINS


  SUBROUTINE point_to_Taylor(x,y,c,x0,m)
!
! Compute the coefficients in Taylor form centered at x0
! of the interpolant of the data y at the nodes x
!
  INTEGER, INTENT(IN) :: m  ! The degree
  DOUBLE PRECISION, DIMENSION(0:m), INTENT(IN) :: x,y
  DOUBLE PRECISION, DIMENSION(0:m), INTENT(OUT) :: c
  DOUBLE PRECISION, INTENT(IN) :: x0
!
  DOUBLE PRECISION, DIMENSION(0:m) :: dx,newton
  INTEGER :: j,k
!
! Compute the Newton polynomial
!
  dx=x-x0
  newton=y
  DO k=1,m
    DO j=0,m-k
      newton(j)=(newton(j+1)-newton(j))/(dx(j+k)-dx(j))
    END DO
  END DO
!
! Now change from Newton to Taylor
!
  c=0.d0
  c(0)=newton(0)
  DO k=1,m
    DO j=k,1,-1
      c(j)=c(j-1)-dx(k)*c(j)
    END DO
    c(0)=newton(k)-dx(k)*c(0)
  END DO
!
  END SUBROUTINE point_to_Taylor

  SUBROUTINE Hermite(cl,xl,cr,xr,ch,x0,m)
!
! Compute the Hermite interpolant of the degree m Taylor polynomials
! centered at xl,xr and express as a degree 2m+1 Taylor polynomial
! centered at x0 
!
  INTEGER, INTENT(IN) :: m
  DOUBLE PRECISION, INTENT(IN) :: xl,xr,x0
  DOUBLE PRECISION, DIMENSION(0:m), INTENT(IN) :: cl,cr
  DOUBLE PRECISION, DIMENSION(0:2*m+1), INTENT(OUT) :: ch 
!
  DOUBLE PRECISION, DIMENSION(0:2*m+1) :: dx,Newton
  INTEGER :: j,k
!
  dx(0:m)=xl-x0
  dx(m+1:2*m+1)=xr-x0
  Newton(0:m)=cl(0)
  Newton(m+1:2*m+1)=cr(0)
!
  DO k=1,2*m+1
    IF (k <= m) THEN
      DO j=0,m-k
        Newton(j)=cl(k)
      END DO
      DO j=m-k+1,m
        Newton(j)=(Newton(j+1)-Newton(j))/(dx(j+k)-dx(j))
      END DO
      DO j=m+1,2*m+1-k
        Newton(j)=cr(k)
      END DO
    ELSE
      DO j=0,2*m+1-k
        Newton(j)=(Newton(j+1)-Newton(j))/(dx(j+k)-dx(j))
      END DO
    END IF
  END DO
!
! Now change from Newton to Taylor
!
  ch=0.d0
  ch(0)=Newton(0)
  DO k=1,2*m+1
    DO j=k,1,-1
      ch(j)=ch(j-1)-dx(k)*ch(j)
    END DO
    ch(0)=Newton(k)-dx(k)*ch(0)
  END DO
!
  END SUBROUTINE Hermite 

!
  SUBROUTINE point_to_dspline(x,y,c,xl,xr,m,n)
!
! Compute the coefficients of the degree 2m+1 dspline on [xl,xr]
!
! The result is in Taylor form centered at (xl+xr)/2 
!
! Assume data for xl starts at x(0) and for xr ends at x(n)
!
  INTEGER, INTENT(IN) :: m,n
  DOUBLE PRECISION, DIMENSION(0:n), INTENT(IN) :: x,y
  DOUBLE PRECISION, DIMENSION(0:2*m+1), INTENT(OUT) :: c
  DOUBLE PRECISION, INTENT(IN) :: xl,xr
!
  DOUBLE PRECISION, DIMENSION(0:m) :: yl,xxl,yr,xxr,cl,cr
  DOUBLE PRECISION :: x0
  INTEGER :: j
!
! Check for enough data
!
  IF (n < m) THEN 
    PRINT *,'Not enough data in point_to_dspline'
    STOP
  END IF
!
  yl=y(0:m)
  xxl=x(0:m)
  yr=y(n-m:n)
  xxr=x(n-m:n)
  x0=(xl+xr)/2.d0 
!
  CALL point_to_Taylor(xxl,yl,cl,xl,m)
  CALL point_to_Taylor(xxr,yr,cr,xr,m)
  CALL Hermite(cl,xl,cr,xr,c,x0,m)    
!
  END SUBROUTINE point_to_dspline 

  SUBROUTINE Interp(tnow,t,dt,uev,upev,u,m,n)
!
  INTEGER, INTENT(IN) :: m,n
  DOUBLE PRECISION, INTENT(IN) :: tnow,t,dt
  DOUBLE PRECISION, DIMENSION(n), INTENT(IN) :: u
  DOUBLE PRECISION, INTENT(OUT) :: uev,upev
!
! tnow - current time
! t    - desired time <= tnow
! dt   - time step
! u    - data u(j)=u at time tnow-(n-j)*dt
! n    - dimension of u Must have t >= tnow-(n-2)*dt
! m    - degree 2m+1 spline   
! uev,upev - dspline interpolant and derivative evaluated at t
!
  INTEGER :: k,nh,jr,nr
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: x,y
  DOUBLE PRECISION, DIMENSION(0:2*m+1) :: c
  DOUBLE PRECISION :: tl,tr,s  
!
! Where is t?
!
  jr=INT((tnow-t)/dt)
  tr=tnow-dt*DBLE(jr)
  tl=tr-dt 
  s=t-.5d0*(tl+tr) 
!
! t is located between n-jr-1 and n-jr 
!
! Are we in a regular case or near the present time?
!
  IF ((n-jr+m/2) <= n) THEN
    nh=m+1
    nr=n-jr+m/2
  ELSE
    nh=m
    nr=n
  END IF
  ALLOCATE(x(0:nh),y(0:nh))
  IF ((nr-nh) < 0) THEN
    y(nh-nr:nh)=u(0:nr)
    y(0:nh-nr-1)=0.d0
  ELSE
    y(0:nh)=u(nr-nh:nr)
  END IF
  DO k=0,nh
    x(k)=tnow-dt*DBLE(n-nr+nh-k)
  END DO 
  CALL point_to_dspline(x,y,c,tl,tr,m,nh)
  uev=c(2*m+1)
  upev=DBLE(2*m+1)*c(2*m+1)
  DO k=2*m,1,-1 
    uev=s*uev+c(k)
    upev=s*upev+DBLE(k)*c(k)
  END DO
  uev=s*uev+c(0)
!
  DEALLOCATE(x,y)
!
  END SUBROUTINE Interp 

  SUBROUTINE Extrap(xcof,m)
!
  INTEGER, INTENT(IN) :: m
  DOUBLE PRECISION, DIMENSION(0:m), INTENT(OUT) :: xcof
!
! Return extrapolation coefficients for degree m
!
  INTEGER :: k,j
!
  DO k=0,m
    xcof(k)=1.d0
    DO j=0,m
      IF (j /= k) THEN
        xcof(k)=xcof(k)*DBLE(m+1-j)/DBLE(k-j)
      END IF
    END DO
  END DO 
!
  END SUBROUTINE Extrap  

END MODULE dspline 
 
      
  
