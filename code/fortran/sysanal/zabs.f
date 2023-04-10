c.........................................................................
c

      SUBROUTINE pullzabs

c
c     this routine is called if the user is inputing the redshift from
c     the "zabs.dat" file, which must exists; if it DNE, or cannot be
c     opened, or the value is not readable, then we halt the program
c     (not pretty, I know- but I can only handle some much of a 100
c     pound Gorilla running my codes)
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      implicit none

      include             'sysanal.h'

      integer              istat,access
      character*80         zabs_file
      include             'sysanal.com'
      
      zabs_file = 'zabs.dat'
      istat    = access(zabs_file,'r')  ! returns 1=DNE, 0="OK"

c     ERROR TRAP - terminal.  if file DNE, we quit

      IF (istat.ne.0) then
        WRITE(6,600) 
        STOP
      END IF

c     Ok to try and read in the redshift

      OPEN(unit=1,file=zabs_file,status='old',err=11)
      READ(1,*,err=12) zbar
      CLOSE(unit=1)

c     everything seems fine

      RETURN

c     ERROR TRAP - teminal, file still not there ????

 11   WRITE(6,601)      
      STOP

c     ERROR TRAP - terminal, trouble reading redshift value

 12   WRITE(6,602)
      STOP

 600  FORMAT(1x,'DEATH(pullzabs):',/,
     &       1x, 'User specified redshift must be in file: zabs.dat',/,
     &       1x, 'This file is not found.  If you want SYSANAL to',/,
     &       1x, 'determine the systemic redshift then change the',/,
     &       1x, 'value of ZUSER=0 in the file: sysanal.inp.')
 601  FORMAT(1x,'DEATH(pullzabs):',/,
     &       1x, 'User specified redshift in file: zabs.dat',/,
     &       1x, 'This file is could not be opened.  Fix and',/,
     &       1x, 'rerun SYSANAL.  Otherwise, if you want SYSANAL to',/,
     &       1x, 'determine the systemic redshift then change the',/,
     &       1x, 'value of ZUSER=0 in the file: sysanal.inp.')
 602  FORMAT(1x,'DEATH(pullzabs):',/,
     &       1x, 'User specified redshift in file: zabs.dat',/,
     &       1x, 'The value of the redshift could not be read. Fix',/,
     &       1x, 'and rerun SYSANAL.  If you want SYSANAL to',/,
     &       1x, 'determine the systemic redshift then change the',/,
     &       1x, 'value of ZUSER=0 in the file: sysanal.inp.')

      END

c.........................................................................
c

      SUBROUTINE getzabs

c
c     Get the system redshift; this uses the technique of finding the
c     wavelength at the optical depth median of the primary ion
c     transition in the input list
c
c     we use several modified Numerical Recipe routines which are in the
c     module "recipes.f", which is included in the main module
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      implicit none

      include             'sysanal.h'

      integer             i,lpix,upix,npix
      external            odepth,halfarea
      double precision    odepth,halfarea,zbrent,fullarea,eps
      double precision    aa,bb,x,y,coeffs

      include             'sysanal.com'

c     an uncommon ocurrance of local common blocks

      COMMON /idata/     npix      
      COMMON /fdata/     x(nmx),y(nmx),coeffs(nmx)
      COMMON /adata/     aa,fullarea


c     set the tolerance for the integration

      eps = toler

c     use the first ion transition in the list

c     get the pixel extremes over all features for this ion transition

      lpix = f_beg(1,1,1)
      upix = f_end(nlines,1,1)
      npix = 0

c     set the integration limits

      aa = wave(lpix,1)
      bb = wave(upix,1)

c     stuff the temporary data for the function (integrand)

      DO 20 i=lpix,upix
       npix = npix + 1
       x(npix)  = wave(i,1)
       y(npix)  = tau(i,1)
 20   CONTINUE

c     compute the spline corefficients for the interative convergence

      CALL spline(x,y,npix,coeffs)

c     compute the full area under the integrand; this value is used in
c     the root solver

      CALL qromb(odepth,aa,bb,eps,fullarea)

c     root solve for the absorber redshift; root solver zbrent returns
c     the wavelength at the median integrated optical depth; we compute
c     the redshift directly upon its success

C      write(6,*) 'CALLING ZBRENT'

      zbar = zbrent(halfarea,aa,bb,eps)/lambda0(1) - 1.0d0

c     return

      RETURN
      END


c.........................................................................
c

      DOUBLE PRECISION FUNCTION odepth(lambda)

c
c     here we use a spline fit to obtain the optical depth at arbitrary
c     wavelength, LAMBDA
c 
c     this function is called by routine qromb iteratively during the
c     root solving by routine zbrent
c                                                                         
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include             'sysanal.h'

      integer             npix
      double precision    taupix,lambda
      double precision    x,y,coeffs

      include             'sysanal.com'

      COMMON /idata/     npix      
      COMMON /fdata/     x(nmx),y(nmx),coeffs(nmx)


c     call the spline at LAMBDA, TAU is returned

      CALL splint(x,y,coeffs,npix,lambda,taupix)

c     set to odepth, the function

      odepth = taupix

c     return

      RETURN
      END


c.........................................................................
c

      DOUBLE PRECISION FUNCTION halfarea(lambda)

c
c     function called by function zbrent to root solve for lambda_zabs
c     set this baby to zero and we are done
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include             'sysanal.h'

      external           odepth
      integer            npix
      double precision   lambda,subarea,eps,odepth
      double precision   x,y,coeffs
      double precision   aa,fullarea

      include             'sysanal.com'

      common /idata/     npix      
      COMMON /fdata/     x(nmx),y(nmx),coeffs(nmx)
      COMMON /adata/     aa,fullarea



c     set the tolerance for the value of HALFAREA from zero

      eps = toler

c     integrate to the get the SUBAREA between AA and LAMBDA; when the
c     SUBAREA is half of the total area, we have our wavelength

      CALL qromb(odepth,aa,lambda,eps,subarea)

c     root solve this function

      halfarea = fullarea - 2.0d0*subarea     

      RETURN
      END

c     eof
