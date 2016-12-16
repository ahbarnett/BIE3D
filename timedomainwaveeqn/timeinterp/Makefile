# Makefile for Hagstrom time-interpolation fortran codes
# Barnett 12/16/16

FC=gfortran
FLIBS= -llapack -lblas
FFLAGS=-fPIC

default: exec

dspline.o: dspline.f90
	$(FC) -c $(FFLAGS) dspline.f90

exec:	TestInterp.f90 TestInterpMat.f90 dspline.o
	$(FC) $(FFLAGS) TestInterp.f90 dspline.o $(FLIBS) -o TestInterp
	./TestInterp
	$(FC) $(FFLAGS) TestInterpMat.f90 dspline.o $(FLIBS) -o TestInterpMat
	./TestInterpMat

mwrap:	gateway.mw
	blah..


clean:
	rm -f *.o *.mod *.out
