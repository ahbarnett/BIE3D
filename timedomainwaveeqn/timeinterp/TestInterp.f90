PROGRAM TestInterp

USE dspline

IMPLICIT NONE

INTEGER :: m,n,ip,j,k,jmin,jmax,jd
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: u,p,xcof
DOUBLE PRECISION :: dt,t,tnow,uev,upev,utrue,uptrue,frac,s

n=500
tnow=17.137d0
dt=tnow/DBLE(n)
ALLOCATE(u(n))

DO m=2,12,2
  ALLOCATE(p(0:m),xcof(0:m))
  CALL Extrap(xcof,m)
  DO j=0,m
    CALL RANDOM_NUMBER(p(j))
    p(j)=-1.d0+2.d0*p(j)
  END DO 
  DO ip=0,1 
    PRINT *,m,ip
    DO k=(m/2)+1,n-1
      CALL RANDOM_NUMBER(frac) 
      t=dt*(DBLE(k)+frac)
      s=dt*frac
      utrue=p(m)
      uptrue=DBLE(m)*p(m)
      DO jd=m-1,1,-1 
        utrue=s*utrue+p(jd)
        uptrue=s*uptrue+DBLE(jd)*p(jd)
      END DO
      utrue=s*utrue+p(0) 
      u=0.d0
      jmin=MAX(1,k-m)
      jmax=MIN(n,k+(m/2)+2)
      DO j=jmin,jmax 
        s=dt*DBLE(j-k)
        u(j)=p(m)
        DO jd=m-1,0,-1
          u(j)=s*u(j)+p(jd)
        END DO
      END DO 
      IF ((ip==1).AND.(jmax==n)) THEN
        u(n)=0.d0
        DO j=0,m
          u(n)=u(n)+xcof(j)*u(n+j-m-1)
        END DO
      END IF 
      CALL Interp(tnow,t,dt,uev,upev,u,m,n)
      PRINT *,k,ABS(utrue-uev),ABS(uptrue-upev)  
    END DO
  END DO
  DEALLOCATE(p,xcof)
END DO

END PROGRAM TestInterp 

  
      
      
