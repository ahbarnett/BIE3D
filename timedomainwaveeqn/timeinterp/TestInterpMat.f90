PROGRAM TestInterpMat 

USE dspline

IMPLICIT NONE

INTEGER :: m,n,j,k,jmin,jmax,jd,r,ir,nc,kloc 
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: u,p,tinterp,ui,upi 
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: umat,upmat 
DOUBLE PRECISION :: dt,t,tnow,uev,upev,utrue,uptrue,frac,s,erloc,uermax,upermax,tloc   

n=500
tnow=17.137d0
dt=tnow/DBLE(n)

r=2157
ALLOCATE(tinterp(r),ui(r),upi(r))
!
DO m=2,12,2
  ALLOCATE(p(0:m))
  DO j=0,m
    CALL RANDOM_NUMBER(p(j))
    p(j)=-1.d0+2.d0*p(j)
  END DO 
  PRINT *,m
!
! choose a bunch of points in 3 time steps
!
  DO kloc=0,314,314 
!
  tloc=-dt*DBLE(kloc)
!
  DO ir=1,r
    CALL RANDOM_NUMBER(frac)
    tinterp(ir)=tloc-3.d0*dt*frac 
  END DO
  CALL InterpMat(r,tinterp,dt,m,jmax,jmin,umat,upmat)
!
! Create data - here we assume that the polynomial is centered at 0
!
  nc=jmax-jmin+1
  ALLOCATE(u(nc)) 
  DO k=1,nc 
    s=dt*DBLE(jmin+k-1)-tloc
    u(k)=p(m)
    DO jd=m-1,0,-1
      u(k)=s*u(k)+p(jd)
    END DO 
  END DO
!
!  DO ir=1,r
!  DO k=1,nc
!    PRINT *,ir,k,umat(ir,k),upmat(ir,k)
!  END DO
!  END DO 
  CALL dgemv('N',r,nc,1.d0,umat,r,u,1,0.d0,ui,1)
  CALL dgemv('N',r,nc,1.d0,upmat,r,u,1,0.d0,upi,1)
!
  uermax=0.d0
  upermax=0.d0 
  DO ir=1,r  
    s=tinterp(ir)-tloc
    utrue=p(m)
    uptrue=DBLE(m)*p(m)
    DO jd=m-1,1,-1 
      utrue=s*utrue+p(jd)
      uptrue=s*uptrue+DBLE(jd)*p(jd)
    END DO
    utrue=s*utrue+p(0)
    erloc=ABS(ui(ir)-utrue)
    IF (erloc > uermax) THEN
      uermax=erloc
    END IF  
    erloc=ABS(upi(ir)-uptrue)
    IF (erloc > upermax) THEN
      upermax=erloc
    END IF
  END DO
!
  PRINT *,kloc,uermax,upermax   
  DEALLOCATE(u,umat,upmat)
!
  END DO 
!
  DEALLOCATE(p)
END DO

DEALLOCATE(tinterp,ui,upi) 

END PROGRAM TestInterpMat  

  
      
      
