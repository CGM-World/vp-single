FFLAGS=-O 
FC=f95
#FFLAGS=-O -Bstatic
#TLIB=/home/matrix/cwc/Programs/tlib/tool.a
#MGOLIBS=-lmongo -lreadline -lhistory -lcurses /usr/lib64/libX11.so.6
#SLATEC=/home/matrix/cwc/Programs/slatec/libslatec.a -llapack
#       Insure that make knows what to do with C preprocessored fortran
.F.o: ; $(FC) $(FFLAGS) -c $*.F
#

default:: sysanal

sysanal: const.dek indata.f zabs.f test.f features.f calc.f parzers.f spectra.f sysanal.com sysanal.f sysanal.h subs.f access.f
	$(FC) -w $(FFLAGS) sysanal.f access.f -o sysanal
