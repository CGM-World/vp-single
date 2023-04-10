c------------------------------------------------------------------------------
c
c
c
	SUBROUTINE        getparams(paramlist,presel)
c
	include           'specsynth.h'
        integer           conflag
        double precision  lamc,rat,presel
        character*80      mf_path,def_file,par_file,paramlist
        include     'const.dek'
c
c
c
c  check current directory for parameter file, if DNE goto error trap, which
c  grabs the default file which lives where the user has set their MFDIR
c  environment variable
c

        convolving = .false.

c
c  did we enter the parameter file on the command line, or
c  if not, use the default
c

        IF (paramlist.eq.' ') then 
          def_file   = 'specsynth.par'
        ELSE
          def_file = paramlist
        END IF

	OPEN(unit=3,file=def_file,err=999,status='old')

 01	read(3,*)                 ! by pass the header
        read(3,*) ndata
        read(3,*) presel
        read(3,*) R_fac
        read(3,*) slit
        read(3,*) conflag
        read(3,*) conwindo
        read(3,*) resfac        ! must be odd integer
	read(3,*) snr
c
	CLOSE(3)

        IF (conflag.eq.1) convolving = .true.
 
        rat = float(ndata)/(2.0d0*presel*R_fac)
        lamc = lambda0(1) * (1.0d0+zabs)
        wave_max = lamc * (1.0d0+rat)   
        wave_min = lamc * (1.0d0-rat)  
        dwave = (wave_max - wave_min)/ float(ndata)

      IF (ndata.gt.maxpix) then
        WRITE(SCREEN,*) 'ERROR(initspectrum): NDATA > MAXPIX'
        STOP
      END IF


      RETURN
c
c
c
c
c

c  error opening local parameter file; check for the default settings
c  get the Unix environment variable MFDIR, which contains the path
c  to where the file lives;  try to open it; on fail, abort 

 999	CALL getenv('MFDIR',mf_path)
        CALL sappend(mf_path,def_file,par_file)
        OPEN(unit=3,file=par_file,err=1999,status='old')
        GOTO 01
c
 1999	write(STDOUT,*) ' NOT FOUND- specsynth.par: parameter listing'
        STOP ' *** SPECSYNTHT(getparams): terminated ***'
c
	END
c

