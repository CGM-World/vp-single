FFLAGS=-O 
FC=f95
#FFLAGS=-O -Bstatic
TLIB=/home/matrix/cwc/Programs/tlib/tool.a
#MGOLIBS=-L/usr/X11R6/lib -lmongo -lX11 -lreadline -lhistory -lcurses
#SLATEC=/home/matrix/cwc/Programs/slatec/libslatec.a -llapack
#       Insure that make knows what to do with C preprocessored fortran
.F.o: ; $(FC) $(FFLAGS) -c $*.F
#

default:: rebin

rebin: rebinspec.f rebinspec.h rebinspec.com bindata.f
	$(FC) $(FFLAGS) rebinspec.f -o rebinspec

rebinv0: rebinspec-v0.f rebinspec.h rebinspec.com bindata.f
	$(FC) $(FFLAGS) rebinspec-v0.f -o rebinspec-v0 $(TLIB) 
