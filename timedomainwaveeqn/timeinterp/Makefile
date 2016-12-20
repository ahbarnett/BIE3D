# Makefile for Hagstrom time-interpolation fortran codes
# Barnett 12/16/16
#
# Note the MEX interface is to the f77 version, interpmatnoalloc, which
# uses already-allocated and fixed size incoming umat, upmat arrays

FC=gfortran
FLIBS= -llapack -lblas
FFLAGS=-fPIC -O2

# Change this to location of MATLAB's command line mex executable:
# (For a Mac system, this is something of the form
#  /Applications/MATLAB_R2014a.app/bin/mex )
MEX=mex

# Change this to location of your MWrap executable if you don't want to use
# the shipped version:
MWRAP=mwrap

default: mwrapclean mexfile

dspline.o: dspline.f90
	$(FC) -c $(FFLAGS) dspline.f90

exec:	TestInterp.f90 TestInterpMat.f90 dspline.o
	$(FC) $(FFLAGS) TestInterp.f90 dspline.o $(FLIBS) -o TestInterp
	./TestInterp
	$(FC) $(FFLAGS) TestInterpMat.f90 dspline.o $(FLIBS) -o TestInterpMat
	./TestInterpMat

# to remake the C file that MEX will compile... (not needed by basic user)
gateway.c: gateway.mw Makefile
	$(MWRAP) -c99complex -mex gateway -mb -list gateway.mw
	$(MWRAP) -c99complex -mex gateway -c gateway.c gateway.mw
# hack to link to a F90 module which mwrap seems not to. Replaces 2 lines,
# in particular the f2c line starting #define MWF77_interpmat...
	sed -i 's/interpmatnoalloc_/__dspline_MOD_interpmatnoalloc/' gateway.c

# to do the MEX compilation...
mexfile: gateway.c dspline.o Makefile
	$(MEX) gateway.c dspline.o -largeArrayDims -lgfortran -lm

clean:
	rm -f *.o *.mod *.out *.mex*

# clean the stuff that mwrap makes... (careful; not needed by basic user)
mwrapclean: clean
	rm -f interpmat.m gateway.c
