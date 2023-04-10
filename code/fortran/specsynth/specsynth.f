c..............................................................................
c
c
      PROGRAM    spectra
c
c
c  reads in the minfit.vdat file and computes the model spectra
c  dor plotting purposes
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'specsynth.h' 
      logical            verbose
      integer            m,idum
      double precision   wcen,presel
      character*80       linelist,paramlist,ch_verb,ch_idum



      verbose = .false.

      lines = 0

c     read in the list of ion names, input the ion information and book
c     keeping, input the atomic constants, etc., grab the parameters

      CALL getarg(1,linelist)
      CALL getarg(2,paramlist)  
      CALL getarg(3,ch_idum)  
      CALL getarg(4,ch_verb)  

      READ(ch_idum,'(i6)') idum
      IF (ch_verb.eq.'1') verbose = .true.

      CALL readlines(linelist)             ! load the line list
      CALL getatomic                       ! assigns atomic constants
      CALL getparams(paramlist,presel)     ! load the input deck

      CALL droplines(verbose)              ! removes out of range lines
      CALL mkticks(linelist)               ! make the .ticks file
      CALL initspectrum(presel,verbose)    ! initializes the continuum

      m = ndata
      wcen = 0.5*(wave_max+wave_min)       ! this is an approximation
      CALL instrument(m,wcen,verbose)      ! and will need correcting

      CALL doabslines                      ! line flux radiative transfer
      CALL dolymanlimit                    ! Lyman limit break
                                           
      if (convolving) CALL convolve        ! FFT convolution with ISF

      IF (snr.gt.0.0) CALL addnoise(idum)  ! add gaussian deviates with snr
      CALL output(linelist,verbose)        ! write output file is
                                           ! linelist.out

c     we are done


      STOP

      END

c..............................................................................
c
      SUBROUTINE        readlines(linelist)
c
c     read designated files
c
c==============================================================================

      include           'specsynth.h'

      integer           k
      double precision vline,ckms
      character*80      linelist

      ckms     = 2.99792458e5

c     read the transitions and line lists; on input, columns are log10,
c     convert to 1.e13 units

      lines = 1

      OPEN(unit=1,file=linelist,ERR=999,status='old')

      read(1,*) zabs

      DO 11 while (lines.le.maxlines)

       k = lines
       read(1,*,end=12) ion_name(k),vline,nline(k),bline(k)
       zline(k) = zabs+ (vline/ckms)*(1.0d0+zabs)
      IF (nline(k).ge.0.0d0) then
          nline(k) = 10.0d0**(nline(k)-13.0d0)
       ELSE
          nline(k) = -10.0d0**(abs(nline(k))-13.0d0)
       END IF
       IF (INDEX(ion_name(k),'HIs').ne.0) then
         CALL getHIseries
         lines = lines + 32
       ELSE
         lines = lines + 1
       END IF

 11   CONTINUE

      WRITE(SCREEN,*) 'EOF reached in ionlist before MAXLINES'
      WRITE(SCREEN,*) 'linelist has been truncated...'

 12   CLOSE(unit=1)

      lines = lines - 1

      RETURN


 999  write(SCREEN,*) ' ERROR(readlines): input file DNE'
      STOP ' fatal error in PROGRAM specsynth'

      end

c..............................................................................
c  
      SUBROUTINE        initspectrum(presel,verbose)
c  
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      include           'specsynth.h'
      logical           verbose
      integer           j
      double precision  presel,wzabs,v1,v2


c     old version contant wavelenght
C      DO 11 pixi=1,ndata
C       lambda(pixi) = wave_min + (pixi-1)*dwave
C       wrkflx(pixi) = 1.0d0
C 11   CONTINUE

c     new version, 
      lambda(1) = wave_min
      wrkflx(1) = 1.0d0
      j = 1
      DO WHILE (lambda(j).lt.wave_max)
       dwave       = lambda(j)/(presel*R_fac) 
       lambda(j+1) = lambda(j) + dwave
       wrkflx(j+1) = 1.0d0
       j           = j + 1
      ENDDO

      ndata = j
      wave_max = lambda(ndata)
  
      wzabs = (1.0d0+zabs)*lambda0(1)
      v1    = 2.99792458d5 * (wave_min-wzabs)/wzabs
      v2    = 2.99792458d5 * (wave_max-wzabs)/wzabs

      IF (verbose) then
       WRITE(6,*) 'SPECSYNTH: input parameters are:'
       WRITE(6,'(a,i4)')   ' N pixels                ',ndata
       WRITE(6,'(a,f8.3)') ' Starting wavelength     ',wave_min
       WRITE(6,'(a,f8.3)') ' Ending wavelength       ',wave_max
       WRITE(6,'(a,f8.3)') ' Starting wavelength     ',v1
       WRITE(6,'(a,f8.3)') ' Ending wavelength       ',v2

       IF (convolving) then                                                                                                                                
        WRITE(6,'(a,a)') ' ISF convolution               ','ON'                                                                                            
       ELSE                                                                                                                                                
        WRITE(6,'(a,a)') ' ISF convolution              ','OFF'                                                                                            
       END IF                                                                                                                                              
      ENDIF
   
      RETURN
      END

c..............................................................................
c  
      SUBROUTINE        doabslines
c  
c     compute the model
c  
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      include           'specsynth.h'
      include           'const.dek'
      integer           linei,pixi
      double precision  w0,b1,b2,b3,w,x,y,u,v,tcon
c
c  
c  
c  b1  = column density
c  b2  = rest frame central wavelength
c  b3  = rest frame Doppler width
c  y   = natural broadening term in doppler units
c  x   = wavelength difference in doppler units
c  w0  = rest frame wavelength for transition
c  w   = rest frame wavelength for observed wavelength
c  tau = optical depth
c
c

      DO 05 linei=1,lines
  
        IF (INDEX(ion_name(linei),'LLS').ne.0) GOTO 05


c     set up some stuff


        w0   = lambda0(linei)
        b1   = nline(linei)
        b2   = w0 * (1.0d0+zline(linei))/(1.0d0+zabs)
        b3   = w0 * abs(bline(linei)) /ckms
        y    = con2(linei) / b3
        tcon = con1(linei) * (b1/b3)

c     loop over the data and stuff the workflux


        DO 11 pixi=1,ndata

          w = lambda(pixi)/(1.0+zabs)
          x = (w-b2)/b3 

            CALL voigt(x,y,u,v)
            wrkflx(pixi) = wrkflx(pixi) * exp(-tcon*u)

 11     CONTINUE

 05   CONTINUE


      return

      end
  
c..............................................................................
c  
      SUBROUTINE        dolymanlimit
c  
c     compute the model
c  
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      include           'specsynth.h'
      include           'const.dek'
      integer           linei,pixi
      double precision  w0,w,tau,NHI,zHI

  

      DO 03 linei=1,lines

        IF (INDEX(ion_name(linei),'LLS').eq.0) GOTO 03

        NHI = 1.0e13 * nline(linei)  
        zHI = zline(linei) 
        w0  = lambda0(linei)
        tau = 0.0d0

        DO 11 pixi=1,ndata
          w  = lambda(pixi)/(1.0+zHI)
          IF (w.le.w0) then 
            tau = NHI * 6.3d-18 * (w/w0)**3
            wrkflx(pixi) = wrkflx(pixi) * exp(-tau)
          END IF
          IF ((w.gt.w0).AND.(w.lt.912.645)) then      ! quick fix
            tau = NHI * 6.3d-18 
            wrkflx(pixi) = wrkflx(pixi) * exp(-tau)
          END IF
 11     CONTINUE

 03   CONTINUE


      RETURN

      END
  
c..............................................................................
c
      SUBROUTINE         addnoise(idum)
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      include           'specsynth.h'
      include           'const.dek'
      integer           pixi,idum
      double precision  noise,gasdev,Ic,Isig
      double precision  nrat2,rdnoise


      Ic      = 0.5d0*(snr**2)  ! assumes no read noise
      rdnoise = 2.5 

      nrat2 = 4.0d0 * (rdnoise/snr)**2 
      Ic   = 0.5d0*(snr**2) * (1.0d0 + sqrt(1.0d0+nrat2))

c     add the noise and make the sigma spectrum

       do 05 pixi=1,ndata
        Isig  = (1.0d0/Ic)*sqrt(Ic*wrkflx(pixi) + rdnoise**2)
        noise  = gasdev(idum)*Isig
        wrkflx(pixi) = wrkflx(pixi) + noise 
        sigma(pixi)  = Isig
 05    continue

      return

      end

c..............................................................................
c
      SUBROUTINE         output(linelist,verbose)
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      include           'specsynth.h'
      include           'const.dek'
      logical           verbose
      integer           pixi
      double precision  wave,vel,flux,noise,wzabs
      character*80      out_file,linelist


c     to get the velocity, assume the first line in the list is the
c     primary transition

      wzabs = (1.0d0+zabs)*lambda0(1)

      CALL fappend(linelist,'out',out_file)

      OPEN(unit=1,file=out_file,status='unknown')

       do 05 pixi=1,ndata
        wave  = lambda(pixi)
        vel   = 2.99792458d5 * (wave-wzabs)/wzabs
        flux  = wrkflx(pixi)
        noise = sigma(pixi)
        write(1,100) wave,vel,flux,noise,noise,1.0d0
 05    continue

       CLOSE(unit=1)
 03   continue


      IF (verbose) WRITE(6,'(a,a)') 'OUTPUT FILE is : ',out_file

      return

 100  FORMAT(1x,f10.4,1x,f8.2,1x,1p4e13.5)

      end

c..............................................................................
c
      SUBROUTINE  droplines(verbose)
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      include           'specsynth.h'

      logical           verbose
      integer           i,k,dropped
      double precision  w



c     drop lines that are not in the wavelength range



      IF (verbose) 
     &    WRITE(6,'(a,i5)') ' Lines in input list         ',lines

      k       = 1
      dropped = 0

      DO 21  while (k.le.lines)

        w = lambda0(k)*(1.0d0+zline(k))

        IF ((w.lt.wave_min).OR.(w.gt.wave_max)) then
         lines = lines - 1
         dropped = dropped + 1
         DO 11 i=k,lines
           ion_name(i) = ion_name(i+1)
           lambda0(i)  = lambda0(i+1)
           con1(i)     = con1(i+1)
           con2(i)     = con2(i+1)
           zline(i)    = zline(i+1)
           nline(i)    = nline(i+1)
           bline(i)    = bline(i+1)  
 11      CONTINUE
         
        ELSE

         k = k + 1

        END IF

 21   CONTINUE


      If (verbose) 
     &     WRITE(6,'(a,i5)') ' Lines out of spectral range ',dropped
      IF (verbose) 
     &     WRITE(6,'(a,i5)') ' Lines being processed       ',lines


      RETURN

      END 

c.............................................................................
c..............................................................................
c
      SUBROUTINE  mkticks(linelist)
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      include           'specsynth.h'

      integer           i
      double precision  w
      character*80      out_file,linelist


      CALL fappend(linelist,'ticks',out_file)
      OPEN(unit=8,file=out_file,status='unknown')

      DO 21 i=1,lines
        w = lambda0(i)*(1.0d0+zline(i))
        write(8,800) w,1.25,zline(i),ion_name(i)
 21   CONTINUE

      CLOSE(unit=8)

      RETURN

 800  FORMAT(1x,f10.4,2x,f4.2,2x,f10.6,2x,a10)

      END 

c.............................................................................

      include 'convolve.f'
      include 'fft.f'
      include 'instrument.f'
      include 'getparms.f'
      include 'voigt.f'
      include 'getatomic.f'
      include 'strappend.f'
      include 'noise.f'

      include 'interp.f'

c     eof

