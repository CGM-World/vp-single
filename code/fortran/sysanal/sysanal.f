c.........................................................................
c

      PROGRAM sysanal

c     this code performs analysis on absorption systems for which each
c     ion/transition is stored in a ASCII file over some fixed velocity
c     range.  The file formats are
c
c     wavelength, velocity, flux, sigma1, sigma2, continuum
c
c     where the velocity is a preliminary velocity based upon a visually
c     assigned redshift.  However, this column can contain all zeroes on
c     input. Sigma2 could be a sky spectrum, or backgroud spectrum if so
c     desired.  If the data have been resampled or linearized in a
c     previous step, sigma2 could contain a scaled uncertainty spectrum
c     that accounts for correlations between pixels.
c
c     authored by Chris Churchill (cwc@nmsu.edu)
c
c     requires Lick Mongo for interactive mode.  If you do not want to
c     use interactive mode, comment out the call to routine INTERACT and
c     comment out the statement "include interact.f" at the bottom of
c     this module.
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include            'sysanal.h'

      logical            flag,flag2,error
      integer            i,ikey,itmp,jtmp,ktmp,ptype,action,linei
      integer            iorder,ioni,subregi
      double precision   rtmp,value2
      real*4             xpos,ypos
      character*80       list_file,inp_file,ch_tmp,string1,string2

      include            'sysanal.com'


      error = .false.
      verbose = .false.

c     grab the command line argument; if the comman is "help" or "?", call
c     the helpme routine and stop

      CALL getarg(1,list_file)

      IF ((list_file.eq.'help ').OR.(list_file.eq.'?')) then
       CALL helpme
       STOP
      END IF

      CALL getarg(2,inp_file)

      CALL getarg(3,ch_tmp)

      IF (ch_tmp.eq.'verbose') verbose = .true.

      IF (verbose) then 
       WRITE(6,*) ' ::::::::::::::::::::::::::::::::::::::::::'
       WRITE(6,*) ' '
       WRITE(6,*) '               SYSANAL V3.0'
       WRITE(6,*) '                  (2020) '
       WRITE(6,*) ' '
       WRITE(6,*) '          Christopher W. Churchill    '
       WRITE(6,*) '       (New Mexico State University)     '
       WRITE(6,*) ' '
       WRITE(6,*) ' ::::::::::::::::::::::::::::::::::::::::::'
      ENDIF

c     grab the default settings; if error, bail

      CALL getpars(inp_file,error)
      IF (error) STOP 

c     zero out all of the working arrays

      CALL zeroall

c     read in the ion/transition list, if error, bail

      CALL getlist(list_file,error)
      IF (error) STOP 

c     read in the atomic data, if error, bail

      CALL getatomic(error)
      IF (error) STOP 

c     read in the spectra, if error, bail

      CALL getdata(error)
      IF (error) STOP 

c     read in the masks (NOT FUNCTIONAL)

C      CALL getmasks(error)
C      IF (error) STOP 

c     create species names (this is just a test right now)

      CALL speciesnames

c     communicate inputs

      CALL comm_inputs

c     automated computation of the initial system master regions and
c     subregions; get the initial redshift for the auto regions

      CALL initregs

c     if no master regions are found, then we cannot continue

      IF (nlines.eq.0) then
       IF (verbose) then
        WRITE(6,*) ' NO FEATURES IN THIS SYSTEM.'
        WRITE(6,*) ' terminating'
       ENDIF
       STOP
      END IF

c     compute the final results 

      CALL finalcalc

c     communicate the master regions

      CALL comm_regs

c     write the results

      CALL writefiles

      STOP 

      END

c
c-----------------------------------------------------------------------
c     included modules
      include 'calc.f'       ! routines for calculating EWs, etc.
      include 'features.f'   ! routines for feature finding/manipulation
      include 'indata.f'     ! reads in all input data
      include 'parzers.f'    ! string manpulations
      include 'recipes.f'    ! various required Numerical Recipes (modified)
      include 'spectra.f'    ! routines for computing EW, AOD spectra, etc.
      include 'subs.f'       ! various subroutines
      include 'test.f'       ! makes species name
      include 'zabs.f'       ! computes the systetmic redshift 
c-----------------------------------------------------------------------
