c
      PROGRAM rebinspec

c     This program rebins a spectrum 
c     The output spectrum will be normalized

c     Usaage: rebinspec refspec inspec [printflag]
c     refspec is a reference spectrum that provides the pixelization
c     to which inspec will be rebinned
c     logical printflag =0 or blank will suprress writing to screen

      implicit none

      include            'rebinspec.h'
      logical            pflag,clnflag,hit
      integer            i,j,m,nreg
      double precision   dv1,dv2,lambda0,ckms,oneplusz,dum
      double precision   vreg1(mxpix),vreg2(mxpix)
      character*80       refspec,inspec,outspec,ewregfile
      character*80       header,pfstr,cfstr
      include            'rebinspec.com'
    
      CALL getarg(1,refspec)  ! input reference spectrum
      CALL getarg(2,inspec)   ! spectrum to be rebinned
      CALL getarg(3,cfstr)    ! flag for cleaning the rebinned spectrum
      CALL getarg(4,pfstr)    ! flag for printing progress to screen

c     check the command line

      IF (refspec.eq."") STOP 'command line reference spectrum is blank'
      IF (inspec.eq."") STOP 'command line input spectrum is blank'

c     set the cleaning flag

      clnflag = .false.
      IF (cfstr.eq."1") clnflag = .true.

c     set the print to screen flag

      pflag = .false.
      IF (pfstr.ne."")  pflag = .true.   ! any character will set it high
      IF (pfstr.eq."0") pflag = .false.  ! but if it is a "0", the reset it low

c     central wavelength needed to determine system redshift

      ckms    = 2.99792458e5
      lambda0 = 0.0d0
      nreg    = 1
   
c     read in the template spectrum and get the velocity bin centers, vbin

      OPEN(unit=3,file=refspec,err=993,status='unknown')
      DO m=1,mxpix
       READ(3,*,end=101) wbin(m),vbin(m)
      ENDDO       
      WRITE(6,*) 'file truncated: ',refspec
 101  CLOSE(unit=3)
      nbin = m - 1

c     check that binning of template is linear
c     due to rounding errors I had to apply a NINT trick

      dv1 = 0.010d0*nint(100.0d0*(vbin(2) - vbin(1)))
      dv2 = 0.010d0*nint(100.0d0*(vbin(nbin) - vbin(nbin-1)))
      IF (dv2.ne.dv1) then
       STOP 'dv changes in template spctrum'
      ENDIF

c     make the velocity bin edges based on the template spectrum
     
      DO m=1,nbin
       lbin(m) = vbin(m) - 0.5*dv1
       ubin(m) = vbin(m) + 0.5*dv1
      ENDDO

c     communicate ref spectrum details

      IF (pflag) then
       WRITE(6,*) 'REBINSPEC ...'
       WRITE(6,*) 'template spectrum    : ',refspec
       WRITE(6,601) 'npix,vlo,vhi,dv1,dvN : ',nbin,vbin(1),vbin(nbin),
     &               dv1,dv2
      ENDIF

c
c     open file: the spectrum to be rebinned
c     read in the velocity regions of the absorption and store
c     read in the spectum to be rebinned; normalize if neeed be
c     if clnflag is high, replace flux=1 if outside absorption region
c     any pixel with abs(fluxdec/sigma) >= 3.0 is replace with unity

      CALL fappend(inspec,'ewreg',ewregfile)
      OPEN(unit=4,file=ewregfile,err=994,status='unknown')
      READ(4,*) header
      DO j=1,mxpix
       READ(4,*,end=204) dum,dum,dum,dum,dum,vreg1(j),vreg2(j),lambda0
      ENDDO
      WRITE(6,*) 'file truncated: ',ewregfile
 204  CLOSE(unit=4)
      nreg = j - 1

      OPEN(unit=2,file=inspec,err=992,status='unknown')
      DO i=1,mxpix
       READ(2,*,end=202) wave(i),vel(i),flux(i),sigma(i),
     &                   smsig(i),cont(i) 
       flux(i) = flux(i)/cont(i)
       IF (sigma(i).ne.-1.0) sigma(i) = sigma(i)/cont(i)
       hit = .false.
       DO j=1,nreg  ! cleaning loop followed by cleaning traps
        IF ((vel(i).ge.vreg1(j)).AND.(vel(i).le.vreg2(j))) hit=.true.  ! inside detection region
       ENDDO
       IF (clnflag.AND.(.not.hit)) flux(i) = 1.0d0
       IF ((sigma(i).eq.-1.0).AND.(.not.hit)) flux(i) = 1.0d0  ! bad pixel outside region
       dum = abs(1.0d0-flux(i))/sigma(i)
       IF ((dum.ge.3.0).AND.(.not.hit)) flux(i) = 1.0d0  ! bad pixel outside region
      ENDDO
      WRITE(6,*) 'file truncated: ',inspec
 202  CLOSE(unit=2)
      ndata = i - 1

c     compute the spectrum pixel widths 

      lvel(1) = vel(1) - 0.5d0*(vel(2)-vel(1)) ! reflection
      uvel(1) = vel(1) + 0.5d0*(vel(2)-vel(1))
      DO i=2,ndata-1
       dv1 = vel(i) - vel(i-1)
       dv2 = vel(i+1) - vel(i)
       lvel(i) = vel(i) - 0.50d0*dv1
       uvel(i) = vel(i) + 0.50d0*dv2
      ENDDO
      j = ndata
      lvel(j) = vel(j) - 0.5d0*(vel(j)-vel(j-1))
      uvel(j) = vel(j) + 0.5d0*(vel(j)-vel(j-1)) ! reflection

      IF (pflag) then
       WRITE(6,*) 'target spectrum      : ',inspec
       WRITE(6,601) 'npix,vlo,vhi,dv1,dvN : ',ndata,vel(1),vel(ndata),
     &               vel(2)-vel(1),vel(ndata)-vel(ndata-1)
       IF (clnflag)      WRITE(6,*) '-- cleaning is ON'
       IF (.not.clnflag) WRITE(6,*) '-- cleaning is OFF'
      ENDIF

c     deterine the system redshift and shift the wavelength bins

      oneplusz = (wave(1)/lambda0) / (1.0d0+vel(1)/ckms)
      DO m=1,nbin
       wbin(m) = lambda0*oneplusz*(1.0d0+vbin(m)/ckms)
      ENDDO

c     rebin the spectrum

      CALL bindata
    
c     output the result

      CALL fappend(inspec,'rebin.txt',outspec)
      OPEN(unit=1,file=outspec,err=991,status='unknown')
      DO m=1,nbin
       WRITE(1,300) wbin(m),vbin(m),fbin(m),sbin(m),sbin(m),cbin(m)
      ENDDO
      CLOSE(unit=1)
      IF (pflag) WRITE(6,*) 'output file is ',outspec

c     we are done

      STOP

c     error traps

 994  WRITE(6,*) 'Cannot open : ',ewregfile
      STOP 'rebinspec terminated with error'
 993  WRITE(6,*) 'Cannot open : ',refspec
      STOP 'rebinspec terminated with error'
 992  WRITE(6,*) 'Cannot open : ',inspec
      STOP 'rebinspec terminated with error'
 991  WRITE(6,*) 'Cannot open : ',outspec
      STOP 'rebinspec terminated with error'

 300  FORMAT(1x,f11.4,f9.2,1pe12.4,1pe12.4,1pe12.4,1pe12.4,1pe12.4)
 601  FORMAT(1x,a23,i3,4f9.2)
      END

      include 'bindata.f'

c=======================================================================

      SUBROUTINE fappend(infile,delim,outfile)   


      implicit none
      integer             i,k,lend      
      character*(*)       infile,delim,outfile

      lend = len(infile)
      k    = 0

      DO 09 i=1,lend
        k = i
        IF ((infile(i:i).eq.'.').OR.(infile(i:i).eq.' ')) GOTO 10
 09   CONTINUE

 10   IF (delim.eq.' ') then
       outfile = infile(1:k-1)
      ELSE
       IF (infile(k:k).eq.'.') outfile = infile(1:k)//delim
       IF (infile(k:k).eq.' ') outfile = infile(1:k-1)//"."//delim
      END IF


      RETURN
      END

c=======================================================================
