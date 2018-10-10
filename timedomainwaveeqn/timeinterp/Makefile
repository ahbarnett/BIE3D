# Makefile for Hagstrom time-interpolation fortran drivers & MEX. Unix/MAC only.
# Alex Barnett 12/16/16
#
# Note the MEX interface is to Alex's f77 versions, eg interpmatnoalloc which
# uses already-allocated and fixed size incoming umat, upmat arrays

FC=gfortran
FLIBS= -llapack -lblas          # for driver linking (see MEX libs/flags below)
#FFLAGS=-fPIC -O2
# makes no difference...
FFLAGS=-fPIC -O3 -march=native -funroll-loops

# Change this to location of MATLAB's command line mex executable:
# (For a Mac system, this is something of the form
#  /Applications/MATLAB_R2014a.app/bin/mex )
MEX=mex

# If you want to rebuild the MEX gateway for some reason, set this to location
# of your MWrap executable:
# (after installing mwrap from http://www.cs.cornell.edu/~bindel/sw/mwrap )
MWRAP=mwrap

default: mexfile                # usual: compile a MEX interface (not run mwrap)
#default: mwrapclean mexfile    # Barnett default runs mwrap

# generic fortran modules
%.o: %.f90
	$(FC) -c $(FFLAGS) $<

%.mod: %.f90
	$(FC) -c $(FFLAGS) $<

# needed for some weird f90 thing (actually specfuns.mod is produced)...
volterra.o: SpecFuns.mod

# this rule makes and runs all Tom's f90 drivers...
fdrivers: TestInterp.f90 TestInterpMat.f90 dspline.o volterra.o
	$(FC) $(FFLAGS) TestInterp.f90 dspline.o $(FLIBS) -o TestInterp
	./TestInterp
	$(FC) $(FFLAGS) TestInterpMat.f90 dspline.o $(FLIBS) -o TestInterpMat
	./TestInterpMat
	$(FC) $(FFLAGS) volterrasolve.f90 dspline.o volterra.o SpecFuns.o $(FLIBS) -o volterrasolve
	./volterrasolve

# to remake the C file that MEX will compile... (not needed by basic user)
gateway.c: gateway.mw Makefile
	$(MWRAP) -c99complex -mex gateway -mb -list gateway.mw
	$(MWRAP) -c99complex -mex gateway -c gateway.c gateway.mw
# Hack to link to F90 modules, which mwrap seems not to. Replaces 2 lines,
# in particular the f2c line starting #define MWF77_interpmat...
	sed -i 's/interpmatnoalloc_/__dspline_MOD_interpmatnoalloc/' gateway.c
	sed -i 's/extrapnoalloc_/__dspline_MOD_extrapnoalloc/' gateway.c

# to do the MEX compilation via MATLAB's compiler...
mexfile: gateway.c dspline.o Makefile
	$(MEX) gateway.c dspline.o -largeArrayDims -lgfortran -lm

# remove objects, executables, outputs, but not the mwrap-generated stuff...
clean:
	rm -f *.o *.mod *.out *.mex* fort.*

# clean the stuff that mwrap makes... (careful; not needed by basic user!)
mwrapclean: clean
	rm -f interpmat.m extrap.m gateway.c
