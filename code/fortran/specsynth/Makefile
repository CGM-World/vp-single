FFLAGS=-O
#FFLAGS=-O -Bstatic
TLIB=/home/matrix/cwc/Programs/tlib/tool.a
MGOLIBS=-L/usr/X11R6/lib -lmongo -lX11 -lreadline -lhistory -lcurses
SLATEC=/home/matrix/cwc/Programs/slatec/libslatec.a -llapack
#       Insure that make knows what to do with C preprocessored fortran
.F.o: ; f95 $(FFLAGS) -c $*.F
#

default:: specsynth
#

specsynth: fft.f convolve.f getparms.f specsynth.f getatomic.f instrument.f voigt.f noise.f specsynth.h interp.f strappend.f implno.dek const.dek
	f95 -w $(OPTIONS) specsynth.f -o specsynth 
#	f95 -w $(OPTIONS) specsynth.f -o specsynth $(TLIB) 


specsynth2: fft.f convolve.f getparms2.f specsynth2.f getatomic.f instrument2.f voigt.f noise.f specsynth.h 
	f95 -w $(OPTIONS) specsynth2.f -o specsynth2 $(TLIB) 
