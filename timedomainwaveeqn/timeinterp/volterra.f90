MODULE volterra 

! Subroutines for volterra integral equations using dsplines 

USE dspline

USE SpecFuns 

IMPLICIT NONE 

CONTAINS

  SUBROUTINE buildcofs(q,nnodes,tspan,Ker,m,dt,c)
!
! Build the coefficients of the dspline implicit time stepping
! scheme for the Volterra equation
!
! u(t)+int_0^tspan Ker(s)*u(t-s) ds = f(t)
!
! Assuming u,f=0 t <= -tspan
!
! Time integral using Gaussian quadrature with nnodes 
!
  INTEGER, INTENT(IN) :: q,nnodes,m  ! Choose (m-q-1)*dt = tspan
  DOUBLE PRECISION, INTENT(IN) :: tspan 
  DOUBLE PRECISION, DIMENSION(0:m), INTENT(OUT) :: c
  DOUBLE PRECISION, INTENT(OUT) :: dt
  INTERFACE
    DOUBLE PRECISION FUNCTION Ker(t)
      DOUBLE PRECISION, INTENT(IN) :: t
    END FUNCTION Ker
  END INTERFACE 
!
  DOUBLE PRECISION, DIMENSION(0:nnodes) :: s,w,uwgts
  DOUBLE PRECISION, DIMENSION(nnodes+1) :: z
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: umat,upmat
  INTEGER :: jmax,jmin,j,jr,k 
!
  dt=tspan/DBLE(m-q-1)
  CALL GaussQCofs(s,w,nnodes) 
  s=.5d0*tspan*(s+1.d0)
  w=.5d0*tspan*w 
!
! Compute weights
!
  DO j=0,nnodes
    jr=nnodes-j
    uwgts(jr)=Ker(s(j))*w(j)
    z(jr+1)=-s(j)
  END DO
!
! Compute interpolation matrix
!
  CALL InterpMat(nnodes+1,z,dt,2*q,jmax,jmin,umat,upmat)
!
  c=0.d0
  DO j=1,jmax-jmin+1
    jr=j+m+jmin-jmax-1
    DO k=0,nnodes
      c(jr)=c(jr)+umat(k+1,j)*uwgts(k)
    END DO
  END DO
  c(m)=c(m)+1.d0 
!
!  DO j=0,m
!    PRINT *,j,c(j) 
!  END DO 
  DEALLOCATE(umat,upmat)
!
  END SUBROUTINE BuildCofs
!
END MODULE volterra  


