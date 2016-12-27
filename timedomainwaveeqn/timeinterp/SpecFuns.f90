MODULE SpecFuns

! Module for spectral differentiation and integration 

IMPLICIT NONE

DOUBLE PRECISION, PARAMETER :: pi=3.1415926535897932385d0

CONTAINS 

  SUBROUTINE fft(a,b,nd,id)
!
! Cooley-Tukey FFT as copied from Fornberg Appendix F except
! stride 1 is assumed 
!
! Replaces a+ib with sum_k=0^(nd-1) (a_k+ib_k)*exp(2*pi*i*id*j*k/nd)
!
  INTEGER, INTENT(IN) :: nd
  DOUBLE PRECISION, INTENT(IN) :: id
  DOUBLE PRECISION, DIMENSION(0:nd-1), INTENT(INOUT) :: a,b
  INTEGER :: j,i,k,l,lh,ip
  DOUBLE PRECISION :: tr,ti,s,c,ur,ui
!
  j=0
  DO i=0,nd-2
    IF (i < j) THEN
      tr=a(j)
      a(j)=a(i)
      a(i)=tr
      ti=b(j)
      b(j)=b(i)
      b(i)=ti
    END IF 
    k=nd/2
    DO WHILE (k <= j)
      j=j-k
      k=k/2
    END DO
    j=j+k
  END DO
!
  s=0.d0
  c=-1.d0
  l=1
  DO WHILE (l < nd)
    lh=l
    l=l+l
    ur=1.d0
    ui=0.d0
    DO j=0,lh-1
      DO i=j,nd-1,l
        ip=i+lh
        tr=a(ip)*ur-b(ip)*ui
        ti=a(ip)*ui+b(ip)*ur
        a(ip)=a(i)-tr
        b(ip)=b(i)-ti
        a(i)=a(i)+tr
        b(i)=b(i)+ti
      END DO
      ti=ur*s+ui*c
      ur=ur*c-ui*s
      ui=ti
    END DO
    s=SQRT(.5d0*(1.d0-c))*id
    c=SQRT(.5d0*(1.d0+c))
  END DO
!
  END SUBROUTINE fft

  SUBROUTINE tofour(a,b,nd)
!
! Use complex FFT to transform the real periodic pair a,b
! to ahat bhat - Copied from Fornberg Appendix F
!
  INTEGER, INTENT(IN) :: nd
  DOUBLE PRECISION, DIMENSION(0:nd-1), INTENT(INOUT) :: a,b
  INTEGER :: j
  DOUBLE PRECISION :: scl
!
  scl=1.d0/DBLE(nd)
  CALL fft(a,b,nd,-1.d0)
  DO j=0,nd-1
    a(j)=a(j)*scl
    b(j)=b(j)*scl
  END DO
!
  END SUBROUTINE tofour

  SUBROUTINE diffour(a,b,nd)
!
! Differentiate in Fourier space
! Copied from Fornberg Appendix F
!
  INTEGER, INTENT(IN) :: nd
  DOUBLE PRECISION, DIMENSION(0:nd-1), INTENT(INOUT) :: a,b
  INTEGER :: i,j,nh
  DOUBLE PRECISION :: a1,f
!
  nh=nd/2
  a(0)=0.d0
  b(0)=0.d0
  a(nh)=0.d0
  b(nh)=0.d0
  DO i=1,nh-1
    j=nd-i
    f=DBLE(i)
    a1=a(i)*f
    a(i)=-b(i)*f
    b(i)=a1
    a1=a(j)*f
    a(j)=b(j)*f
    b(j)=-a1
  END DO 
!
  END SUBROUTINE diffour

  SUBROUTINE LaplaceBeltrami(Lw,w,theta,M,N)
!
! Use Fourier PS to approximate the Laplace-Beltrami operator on the
! sphere. Here we use Fornberg's filtered doubly periodic PS approximation.
! See Chapter 6
!
  INTEGER, INTENT(IN) :: M,N
  DOUBLE PRECISION, DIMENSION(M,0:N-1), INTENT(IN) :: w
  DOUBLE PRECISION, DIMENSION(M,0:N-1), INTENT(OUT) :: Lw
  DOUBLE PRECISION, DIMENSION(M), INTENT(IN) :: theta
  DOUBLE PRECISION, DIMENSION(0:N-1) :: ah,bh,ap,bp 
  DOUBLE PRECISION, DIMENSION(M) :: sint,cott
  INTEGER :: j,k,j1,j2,j3,j4,k1,k2 
  LOGICAL, SAVE :: first_call=.TRUE.
!
  IF (first_call) THEN
    DO k=1,m
      sint(k)=SIN(theta(k))
      cott(k)=1.d0/TAN(theta(k))
    END DO
    first_call=.FALSE.
  END IF 
!
  Lw=0.d0
!
! First compute theta derivatives - extend data periodically 
!
  DO j=0,N/2-1 
    j1=2*j
    j2=2*j+1
    j3=j1+(N/2)
    IF (j3 > N-1) THEN
      j3=j3-N
    END IF 
    j4=j2+(N/2)
    IF (j4 > N-1) THEN
      j4=j4-N
    END IF
    DO k=1,M
      ah(k-1)=w(k,j1)
      bh(k-1)=w(k,j2)
      ah(2*M-k)=w(k,j3)
      bh(2*M-k)=w(k,j4)
    END DO 
!
    CALL tofour(ah,bh,N)
    ap=ah
    bp=bh
    CALL diffour(ap,bp,N)
    ah=ap
    bh=bp
    CALL fft(ah,bh,N,1.d0)
    CALL diffour(ap,bp,N)
    CALL fft(ap,bp,N,1.d0)
    DO k=1,M
      Lw(k,j1)=Lw(k,j1)+cott(k)*ah(k-1)+ap(k-1)
      Lw(k,j2)=Lw(k,j2)+cott(k)*bh(k-1)+bp(k-1)
    END DO
  END DO
!
! Now phi derivatives
!
  DO k=1,M/2
    k1=2*k-1
    k2=2*k
    DO j=0,N-1
      ah(j)=w(k1,j)
      bh(j)=w(k2,j)
    END DO
    CALL tofour(ah,bh,N)
!
! Now filter out the high modes according to sin(theta) - 10
! was chosen based on experiments with TestLB.f90 
!
    DO j=1,(N/2)-1
      IF (DBLE(j) > 10.d0*sint(k1)*DBLE(N/2)) THEN
        ah(j)=0.d0
        ah(N-j)=0.d0
      END IF
      IF (DBLE(j) > 10.d0*sint(k2)*DBLE(N/2)) THEN
        bh(j)=0.d0
        bh(N-j)=0.d0
      END IF
    END DO  
!
    CALL diffour(ah,bh,N)
    CALL diffour(ah,bh,N)
    CALL fft(ah,bh,N,1.d0)
    DO j=0,N-1
      Lw(k1,j)=Lw(k1,j)+ah(j)/(sint(k1)*sint(k1))
      Lw(k2,j)=Lw(k2,j)+bh(j)/(sint(k2)*sint(k2))
    END DO 
  END DO
!
  END SUBROUTINE LaplaceBeltrami

  SUBROUTINE Shift(a,b,nd)
!
! Shift periodic arrays a,b to their midpoints 
!
  INTEGER, INTENT(IN) :: nd
  DOUBLE PRECISION, DIMENSION(0:nd-1), INTENT(INOUT) :: a,b
  INTEGER :: i,nh
  DOUBLE PRECISION :: a1,fr,fi
!
  nh=nd/2
  CALL tofour(a,b,nd)
  DO i=0,nd-1
    fr=COS(pi*DBLE(i)/DBLE(nd))
    fi=SIN(pi*DBLE(i)/DBLE(nd))
    IF (i > nh) THEN
      fr=-fr
      fi=-fi
    END IF 
    a1=a(i)*fr-b(i)*fi
    b(i)=b(i)*fr+a(i)*fi
    a(i)=a1
  END DO 
  CALL fft(a,b,nd,1.d0) 
!
  END SUBROUTINE Shift 

  SUBROUTINE DCT2(a,b,nd)
!
! Type II DCT according to Makhoul - IEEE Transactions on Acoustics, Speech
! and Signal Processing 1 1980 - Use Fornberg's real and imaginary trick
! 
  INTEGER, INTENT(IN) :: nd
  DOUBLE PRECISION, DIMENSION(0:nd-1), INTENT(INOUT) :: a,b
  DOUBLE PRECISION, DIMENSION(0:nd-1) :: ah,bh
  INTEGER :: k,j,nh
  DOUBLE PRECISION :: scl,a1,b1,c,s 
!
  scl=1.d0/DBLE(nd)
  nh=nd/2
  DO j=0,nh-1
    ah(j)=a(2*j)
    bh(j)=b(2*j)
    k=nd-j-1
    ah(k)=a(2*(nd-k)-1)
    bh(k)=b(2*(nd-k)-1)
  END DO   
  CALL fft(ah,bh,nd,-1.d0)
  DO j=0,nd-1
    c=COS(pi*DBLE(j)/DBLE(2*nd))
    s=SIN(pi*DBLE(j)/DBLE(2*nd))
    a(j)=(c*ah(j)+s*bh(j))*scl
    b(j)=(c*bh(j)-s*ah(j))*scl
  END DO
  DO j=1,nh
    k=nd-j
    a1=(a(j)-b(k))
    b1=(a(k)+b(j))
    a(k)=(a(k)-b(j))
    b(k)=(a(j)+b(k))
    a(j)=a1
    b(j)=b1
  END DO 
!
  END SUBROUTINE DCT2

  SUBROUTINE IDCT2(a,b,nd)
!
! Inverse DCT according to Makhoul - IEEE Transactions on Acoustics, Speech
! and Signal Processing 1 1980 - Use Fornberg's real and imaginary trick
! 
  INTEGER, INTENT(IN) :: nd
  DOUBLE PRECISION, DIMENSION(0:nd-1), INTENT(INOUT) :: a,b
  DOUBLE PRECISION, DIMENSION(0:nd-1) :: ah,bh
  INTEGER :: k,j,nh
  DOUBLE PRECISION :: a1,b1,c,s 
!
  nh=nd/2
  ah(0)=a(0)
  bh(0)=b(0)
  DO j=1,nd-1
    ah(j)=a(j)+b(nd-j)
    bh(j)=b(j)-a(nd-j)
  END DO   
  DO j=1,nd-1
    c=COS(pi*DBLE(j)/DBLE(2*nd))
    s=SIN(pi*DBLE(j)/DBLE(2*nd))
    a(j)=(c*ah(j)-s*bh(j))/2.d0
    b(j)=(c*bh(j)+s*ah(j))/2.d0 
  END DO
  CALL fft(a,b,nd,1.d0)
  ah=a
  bh=b
  DO j=0,nh-1
    a(2*j)=ah(j)
    b(2*j)=bh(j)
    k=nd-j-1
    a(2*(nd-k)-1)=ah(k)
    b(2*(nd-k)-1)=bh(k)
  END DO   
!
  END SUBROUTINE IDCT2

  SUBROUTINE filter(ah,bh,nd,ifil)
!
! ifil=1 exponential filter
! ifil=2 2/3 rule
!
! assume modes 0, ... ,(nd/2)-1,-(nd/2), ... ,-1
!
  INTEGER, INTENT(IN) :: nd,ifil
  DOUBLE PRECISION, DIMENSION(0:nd-1), INTENT(INOUT) :: ah,bh
!
  INTEGER :: j,k,nh,m
  DOUBLE PRECISION :: alpha,dN 
!
  nh=nd/2
!
  IF (ifil==1) THEN
    m=36
    alpha=36.d0 ! parameters from Hou Acta Numerica 2009
    dN=1.d0/DBLE(nh) 
    DO k=0,nh-1
      ah(k)=EXP(-alpha*(DBLE(k)*dN)**m)*ah(k) 
      bh(k)=EXP(-alpha*(DBLE(k)*dN)**m)*bh(k)
      j=nh-k
      ah(k+nh)=EXP(-alpha*(DBLE(j)*dN)**m)*ah(k+nh) 
      bh(k+nh)=EXP(-alpha*(DBLE(j)*dN)**m)*bh(k+nh)
    END DO
  ELSE IF (ifil==0) THEN
    DO k=0,nh-1
      IF (k > (2*nd)/3) THEN
        ah(k)=0.d0
        bh(k)=0.d0
      END IF
      j=nh-k
      IF (j > (2*nd)/3) THEN
        ah(k+nh)=0.d0
        bh(k+nh)=0.d0
      END IF
    END DO
  END IF
!
  END SUBROUTINE filter 
      
  SUBROUTINE Poisson2d(f,u,nx,ny)
!
! Solve -lap u = f for doubly 2pi-periodic functions
!
! Use fft without the real and imaginary trick 
! for simplicity - harder than necessary 
!
! Assume without checking that f has mean 0
! 
!
  INTEGER, INTENT(IN) :: nx,ny 
  DOUBLE PRECISION, DIMENSION(0:nx-1,0:ny-1), INTENT(IN) :: f
  DOUBLE PRECISION, DIMENSION(0:nx-1,0:ny-1), INTENT(OUT) :: u
!
  DOUBLE PRECISION, DIMENSION(0:nx-1,0:ny-1) :: fhr,fhi
  DOUBLE PRECISION, DIMENSION(0:nx-1) :: ah,bh 
  DOUBLE PRECISION, DIMENSION(0:ny-1) :: ch,dh
  INTEGER :: jx,jy,kx,ky,nhx,nhy
  DOUBLE PRECISION :: den  
!
  nhx=nx/2
  nhy=ny/2
!
  DO jy=0,ny-1
    DO jx=0,nx-1
      ah(jx)=f(jx,jy)
      bh(jx)=0.d0 
    END DO
    CALL tofour(ah,bh,nx)
    DO jx=0,nx-1
      fhr(jx,jy)=ah(jx)
      fhi(jx,jy)=bh(jx)
    END DO 
  END DO
!
  DO jx=0,nx-1
    DO jy=0,ny-1
      ch(jy)=fhr(jx,jy)
      dh(jy)=fhi(jx,jy)
    END DO
    CALL tofour(ch,dh,ny)
    DO jy=0,ny-1
      fhr(jx,jy)=ch(jy) 
      fhi(jx,jy)=dh(jy)
    END DO
  END DO
!
! Solve
!
  DO jx=0,nhx-1
    kx=jx-nhx
    DO jy=0,nhy-1
      ky=jy-nhy 
      den=DBLE(jx**2+jy**2)
      IF (den > 0.d0) THEN
        fhr(jx,jy)=fhr(jx,jy)/den
        fhi(jx,jy)=fhi(jx,jy)/den
      END IF
      den=DBLE(kx**2+jy**2)
      fhr(jx+nhx,jy)=fhr(jx+nhx,jy)/den
      fhi(jx+nhx,jy)=fhi(jx+nhx,jy)/den
      den=DBLE(jx**2+ky**2)
      fhr(jx,jy+nhy)=fhr(jx,jy+nhy)/den
      fhi(jx,jy+nhy)=fhi(jx,jy+nhy)/den
      den=DBLE(kx**2+ky**2)
      fhr(jx+nhx,jy+nhy)=fhr(jx+nhx,jy+nhy)/den
      fhi(jx+nhx,jy+nhy)=fhi(jx+nhx,jy+nhy)/den
    END DO
  END DO 
!
  DO jx=0,nx-1
    DO jy=0,ny-1
      ch(jy)=fhr(jx,jy)
      dh(jy)=fhi(jx,jy)
    END DO
    CALL fft(ch,dh,ny,1.d0)
    DO jy=0,ny-1
      fhr(jx,jy)=ch(jy) 
      fhi(jx,jy)=dh(jy)
    END DO
  END DO
!
  DO jy=0,ny-1
    DO jx=0,nx-1
      ah(jx)=fhr(jx,jy)
      bh(jx)=fhi(jx,jy) 
    END DO
    CALL fft(ah,bh,nx,1.d0)
    DO jx=0,nx-1
      u(jx,jy)=ah(jx)
    END DO 
  END DO
!
  END SUBROUTINE Poisson2d 
 

  SUBROUTINE ChebyInt(ya,yb,a,b,nd)
!
! Use a Chebyshev expansion to compute int_-1^1 f(x)dx 
!                                    = int_0^pi f(cos(g)) sin(g)dg
! Assume a,b data on the Chebyshev nodes x_j=cos((j+1/2)*pi/nd), j=0,..,nd-1
! Use int_0^pi cos(ng)*sin(g) dg = 0, n odd; = -2/(n^2-1) n even 
!
  INTEGER, INTENT(IN) :: nd
  DOUBLE PRECISION, DIMENSION(0:nd-1), INTENT(IN) :: a,b
  DOUBLE PRECISION, INTENT(OUT) :: ya,yb
  DOUBLE PRECISION, DIMENSION(0:nd-1) :: ah,bh
  INTEGER :: k,nh
!
  nh=nd/2
  ah=a
  bh=b
  CALL DCT2(ah,bh,nd)
  ya=0.d0
  yb=0.d0
  DO k=0,nh-1
    ya=ya-2.d0*ah(2*k)/(DBLE(4*k*k)-1.d0)
    yb=yb-2.d0*bh(2*k)/(DBLE(4*k*k)-1.d0)
  END DO 
!
  END SUBROUTINE ChebyInt 

  SUBROUTINE ChebyDiff(a,da,nd)
!
! Differentiate in Chebyshev space using the recursion
!
! T_{n+1}'   T_{n-1}'
! -------  - -------    = 2T_n
!  n+1         n-1 
!
  INTEGER, INTENT(IN) :: nd
  DOUBLE PRECISION, DIMENSION(0:nd-1), INTENT(IN) :: a
  DOUBLE PRECISION, DIMENSION(0:nd-1), INTENT(OUT) :: da
  INTEGER :: k 
!
  da=0.d0
  da(nd-2)=DBLE(2*(nd-1))*a(nd-1)
  DO k=nd-3,1,-1
    da(k)=DBLE(2*(k+1))*a(k+1)+da(k+2)
  END DO 
  da(0)=a(1)+da(2)/2.d0 
!
  END SUBROUTINE ChebyDiff 

  SUBROUTINE ChebyPR(x,w,nd)
!
! Return the rational Chebyshev points x_j=cot((2j+1)*pi/(2*nd))
! and the integration weights 1/sin^2((2j+1)*pi/(2*nd))
!
! Note that now integration will simply be the midpoint rule 
!
  INTEGER, INTENT(IN) :: nd
  DOUBLE PRECISION, DIMENSION(0:nd-1), INTENT(OUT) :: x,w
!
  INTEGER :: j
  DOUBLE PRECISION :: theta
!
  DO j=0,nd-1
    theta=DBLE(2*j+1)*pi/DBLE(2*nd)
    x(j)=1.d0/TAN(theta)
    w(j)=1.d0/(SIN(theta)**2)
  END DO 
!
  END SUBROUTINE ChebyPR 

  SUBROUTINE I3dC(f,intf,ns,nu,nv) 
!
! Integrate in s,u,v - trapezoid in s,u and Chebyshev in v
!
  INTEGER, INTENT(IN) :: ns,nu,nv
  DOUBLE PRECISION, DIMENSION(0:ns-1,0:nu-1,0:nv-1), INTENT(IN) :: f
  DOUBLE PRECISION, INTENT(OUT) :: intf  
  DOUBLE PRECISION, DIMENSION(0:nv-1) :: a,b 
  DOUBLE PRECISION :: hs,hu,ya,yb
  INTEGER :: ks,ku,nuh,kv
!
  nuh=nu/2
  hs=2.d0*pi/DBLE(ns)
  hu=2.d0*pi/DBLE(nu)
  intf=0.d0
  DO ks=0,ns-1
    DO ku=0,nuh-1
      DO kv=0,nv-1
        a(kv)=f(ks,2*ku,kv)
        b(kv)=f(ks,2*ku+1,kv)
      END DO
      CALL ChebyInt(ya,yb,a,b,nv)
      intf=intf+ya+yb
    END DO
  END DO
!
  intf=intf*hs*hu
!
  END SUBROUTINE I3dC

  SUBROUTINE I3dA(f,intf,wq,ns,nu,nv) 
!
! Integrate in s,u,v - trapezoid in s,u and a quadrature rule with
! weights wq in v 
!
  INTEGER, INTENT(IN) :: ns,nu,nv
  DOUBLE PRECISION, DIMENSION(0:ns-1,0:nu-1,0:nv), INTENT(IN) :: f
  DOUBLE PRECISION, INTENT(OUT) :: intf  
  DOUBLE PRECISION, DIMENSION(0:nv), INTENT(IN) :: wq
  DOUBLE PRECISION :: hs,hu
  INTEGER :: ks,ku,kv
!
  hs=2.d0*pi/DBLE(ns)
  hu=2.d0*pi/DBLE(nu)
  intf=0.d0
  DO kv=0,nv
    DO ku=0,nu-1
      DO ks=0,ns-1
        intf=intf+wq(kv)*f(ks,ku,kv)
      END DO
    END DO
  END DO
!
  intf=intf*hs*hu
!
  END SUBROUTINE I3dA

  SUBROUTINE GaussQCofs(x,w,q)
!
! Compute nodes and weights for Gauss-Legendre quadrature
!
  INTEGER, INTENT(IN) :: q
  DOUBLE PRECISION, DIMENSION(0:q), INTENT(OUT) :: x,w
!
  DOUBLE PRECISION, DIMENSION(0:q) :: dg
  DOUBLE PRECISION, DIMENSION(q) :: bg
  DOUBLE PRECISION, DIMENSION(0:q,0:q) :: V
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: wrk
  INTEGER, DIMENSION(:), ALLOCATABLE :: iwrk 
!
  INTEGER :: j,info,lwrk,liwrk  
!
  lwrk=1+4*(q+1)+(q+1)**2
  ALLOCATE(wrk(lwrk))
  liwrk=3+5*(q+1)
  ALLOCATE(iwrk(liwrk))
!
  dg=0.d0
  DO j=1,q
    bg(j)=DBLE(j)/SQRT(4.d0*DBLE(j**2)-1.d0)
  END DO
!
  CALL dstedc('I',q+1,dg,bg,V,q+1,wrk,lwrk,iwrk,liwrk,info)
  IF (info /= 0) THEN
    PRINT *,'Trouble from dstedc in GaussQCofs: info=',info
    STOP
  END IF 
!
  DO j=0,q
    x(j)=dg(j)
    w(j)=2.d0*V(0,j)*V(0,j)
  END DO
!
  END SUBROUTINE GaussQCofs 

  SUBROUTINE LegendreP(p,x,n,m)
!
! Evaluate the standard Legendre polynomial of degrees 0:n
!
  INTEGER, INTENT(IN) :: n,m
  DOUBLE PRECISION, DIMENSION(m), INTENT(IN) :: x
  DOUBLE PRECISION, DIMENSION(m,0:n), INTENT(OUT) :: p
  INTEGER :: k,j
!
  DO j=1,m
    p(j,0)=1.d0
  END DO
!
  IF (n > 0) THEN
!
  DO j=1,m
    p(j,1)=x(j)
  END DO
!
  IF (n > 1) THEN
!
  DO k=2,n
    DO j=1,m
      p(j,k)=(DBLE(2*k-1)*x(j)*p(j,k-1)-DBLE(k-1)*p(j,k-2))/DBLE(k)
    END DO
  END DO
!
  END IF
!
  END IF
!
  END SUBROUTINE LegendreP 

  SUBROUTINE Romberg(f,fnew,RB,errest,nvar,itmax,it)
!
! Build and evaluate a Romberg acceleration table
! assuming mesh halving and an even error function
!
  INTEGER, INTENT(IN) :: nvar,itmax,it
  DOUBLE PRECISION, DIMENSION(nvar), INTENT(IN) :: fnew
  DOUBLE PRECISION, DIMENSION(itmax,itmax,nvar), INTENT(INOUT) :: RB 
  DOUBLE PRECISION, INTENT(OUT) :: errest 
  DOUBLE PRECISION, DIMENSION(nvar), INTENT(OUT) :: f
!
  DOUBLE PRECISION :: err,errtab,fac 
  INTEGER :: q,k
!
  IF (it==1) THEN
    errest=0.d0
    DO k=1,nvar
      RB(1,1,k)=fnew(k)
      f(k)=fnew(k)
      errest=errest+f(k)**2
    END DO
    errest=SQRT(errest)
  ELSE
    IF (it > itmax) THEN
      PRINT *,'Error from Romberg - maximum iterations exceeded'
      RETURN
    END IF 
    errtab=0.d0
    DO k=1,nvar
      RB(it,1,k)=fnew(k)
      errtab=errtab+(RB(it,1,k)-RB(it-1,1,k))**2
      f(k)=fnew(k)
    END DO
    errtab=SQRT(errtab)
    DO q=1,it-1 
      fac=2.d0**(2*q)
      err=0.d0
      DO k=1,nvar
        RB(it,q+1,k)=(fac*RB(it,q,k)-RB(it-1,q,k))/(fac-1.d0)
        IF (q < it-1) THEN
          err=err+(RB(it,q+1,k)-RB(it-1,q+1,k))**2
        END IF
      END DO
      IF (q < it-1) THEN
        err=SQRT(err)
        IF (err < errtab) THEN
          DO k=1,nvar
            f(k)=RB(it,q+1,k)
          END DO
          errtab=err
        END IF 
      END IF
    END DO
    errest=errtab
!    PRINT *,it,errtab
  END IF
!
  END SUBROUTINE Romberg 
      
   
END MODULE SpecFuns

  
