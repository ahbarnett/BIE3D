PROGRAM volterrasolve

! Tom Hagstrom wrote this to solve the same 2nd-kind Volterra IE as AHB's
! test_volterra.m

USE volterra 

USE SpecFuns

IMPLICIT NONE 

INTEGER :: q,nnodes,n,m,j,it,nsteps,im,jdat  
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: c,u
DOUBLE PRECISION :: ttot,err,size,t,ermax,dt,g,tspan,val  

tspan=1.d0 
ttot=100.d0 

nnodes=128

n=8 
q=4

DO im=1,4 
!
  m=10*im+q+1
  nnodes=16*im 
  ermax=0.d0
  jdat=19+im
  OPEN(UNIT=jdat)
  ALLOCATE(u(0:m),c(0:m))
  CALL BuildCofs(q,nnodes,tspan,Ker,m,dt,c)
  u=0.d0
!
  nsteps=INT(100.d0/dt)
!
  DO it=1,nsteps 
    t=DBLE(it)*dt 
    CALL f(t,val)
    u(m)=val  
    DO j=0,m-1
      u(m)=u(m)-c(j)*u(j)
    END DO 
    u(m)=u(m)/c(m) 
    err=abs(u(m)-sol(t))
    WRITE (jdat,*)t,u(m)
    IF (err > ermax) THEN
      ermax=err
    END IF 
!
! Shuffle
!
    DO j=0,m-1
      u(j)=u(j+1)
    END DO
  END DO
!
  PRINT *,dt,ermax
  DEALLOCATE(c,u) 
  CLOSE(UNIT=jdat)
!
END DO

CONTAINS

DOUBLE PRECISION FUNCTION sol(t)

DOUBLE PRECISION, INTENT(IN) :: t

sol=(1.d0+TANH(2.d0*(t-10.d0)))*SIN(20.d0*t)

END FUNCTION sol

DOUBLE PRECISION FUNCTION Ker(t)

DOUBLE PRECISION, INTENT(IN) :: t

Ker=1.d0

END FUNCTION Ker

SUBROUTINE f(t,val)

DOUBLE PRECISION, INTENT(IN) :: t
DOUBLE PRECISION, INTENT(OUT) :: val 
DOUBLE PRECISION, DIMENSION(0:1024), SAVE :: s,w
LOGICAL :: first_call=.TRUE.
INTEGER :: j

IF (first_call) THEN
  CALL GaussQCofs(s,w,1024)
  s=.5d0*tspan*(s+1.d0)
  w=.5d0*tspan*w
  first_call=.FALSE.
END IF 

val=sol(t)

DO j=0,1024 
  val=val+w(j)*Ker(s(j))*sol(t-s(j))
END DO

END SUBROUTINE f 

END PROGRAM volterrasolve 

