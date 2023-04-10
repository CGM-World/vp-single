c.........................................................................
c

      SUBROUTINE initregs

c
c     Initial computations on the spectra prior to editing
c
c     This routine is called once at the start of the program; 
c
c     (1) creates the EW spectra
c
c     (2) automatigically finds the master regions, 
c
c     (3) automagically finds the subregions for the higher order ion
c     transitions, 
c
c     (4) determines which of the master regions are spurious (have no
c     features in higher order ion transitions), 
c
c     (5) computes the apparent optical depth (AOD) spectra, and 
c
c     (6) determines the initial redshift (prior to editing of the
c     regions and subregions).
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include            'sysanal.h'

      integer            i,linei,sysfeat,subfeat,action
      logical            flag

      include            'sysanal.com'



c     compute the equivalent width spectra

      DO 01 i=1,norders
       CALL ewspec(i)
 01   CONTINUE

c     if the user wants to use the regions from a previous session, we
c     read them in hear, if not, we find them obejctively ourselves

      IF (ewreg) then

       CALL pullregions ! we will terminate within if there are errors

      ELSE

c     first, find the master regions; fill the R_BEG and R_END arrays
c     that store the starting and ending pixels of the system regions

       nlines = sysfeat(flag)

       IF (nlines.eq.0) then
        IF (verbose) then
         WRITE(6,*) ' No features detected in ',order(1)(1:20)
        ENDIF
        RETURN
       END IF

c     then find the regions in higher order ion transitions, if there
c     are any

       IF (norders.gt.1) then 
        DO 02 i=2,norders
         DO 03 linei=1,nlines
          nfind(i,linei) = subfeat(i,linei)
 03      CONTINUE
 02     CONTINUE
       END IF

c     finally, check for spurious master regions, which can be performed
c     only if there are higher order ion transitions

       IF (.not.ewreg.AND.(norders.gt.1).AND.(njtie(1).gt.1)) then

        flag   = .false.
        action = 1
        DO 04 linei=1,nlines
         CALL spurious(linei,action,flag)
 04     CONTINUE
       
        IF (nlines.eq.0) then
         IF (verbose) then
          WRITE(6,*) ' All first pass features in ',order(1)(1:20)
          WRITE(6,*) ' have been removed in routine SPURIOUS'
         END IF
         RETURN
        END IF

       END IF

      END IF ! end block to get the input regions


c     determine the system redshift; this must be called after the
c     aodspec are loaded; the redshift will only be correct on this
c     initialization if the regions a re properly defined.  If we edit
c     them. then we need to recompute the aodspec and the zbaorber
c     redshift

      DO 05 i=1,norders
       CALL aodspec(i)
 05   CONTINUE

c     if the user is fixing the redshift then we pull that from the file
c     "zabs.dat", if not then we find it from the mean optical depth
c     weight of the primary profile

      IF (zfix) then 
       CALL pullzabs   ! we will terminate within if "zabs.dat" not found
       IF (verbose) then
        WRITE(6,*) ' User supplied systemic redshift = ',zbar
       ENDIF
      ELSE
       CALL getzabs
       IF (verbose) then
        WRITE(6,*) ' Initial SYSANAL systemic redshift = ',zbar
       ENDIF
      END IF


c     return

      RETURN
      END


c.........................................................................
c

      INTEGER FUNCTION   sysfeat(flag)

c
c     This routine called once from routine INITREGS
c
c     This function returns the number of significant absorption
c     features, which we call the master regions, and the starting and
c     ending pixels, stored in f_beg and f_end, of each feature.  The
c     method uses the equivalent width spectrum, not the simple flux
c     decrement.  Though it can be used for emission features to, this
c     subroutine is designed to scope out absorption features only
c
c     A feature us found when one pixel has an EW in its pixel that is
c     Nsig beyond the EWSIG uncertainty in that pixel, i.e., EW/EWSIG <
c     -Nsig.  Then the spectrum is scanned backwards to find the blue
c     edge of the feature, which is defined when the ratio EW/EWSIG = -
c     1.  This is the location at about which the flux decrement becomes
c     consistent with the continuum.  Then the spectrum is set back to
c     the initial pixel of the detection and scanned forward for
c     EW/EWSIG = -1.  The pixel locations of the lower and upper bounds
c     are stored in F_BEG and F_END.  Then the pixel number is set to
c     the upper pixel + 1 of the just found feature and the scanning for
c     another Nsig detection is searched for.
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      implicit none

      include            'sysanal.h'

      logical            flag
      integer            ipix,lpix,upix,npix
      integer            iorder,linei
      double precision   Nsig,ratio

      include            'sysanal.com'

C     flag is a dummy for the function call

      sysfeat = 0
      linei    = 0

c     set the ion transistion to the first one in the list; IORDER=1

      iorder   = 1

c     set the significance level

      Nsig     = N_sigma(iorder)

c     scan the full data chunk

      ipix     = 1
      npix     = ndata(iorder)

c     we scan the spectrum; ipix is the pixel meeting the detection
c     criterion; then we scan blueward and then redward looking for the
c     extremities of the feature; we use a while loop because we need to
c     leap frog the feature once we have it defined, thus we do not
c     increment ipix uniformily

      DO 11 WHILE (ipix.lt.npix)

       ratio = ewpix(ipix,iorder)/ewsigpix(ipix,iorder)
       
c     enter scan feature block if we have a significant equivalent
c     width; increment the feature index and then scan for the feature
c     extremities

       IF (ratio.le.-Nsig) then 

         linei = linei + 1
         lpix  = ipix
         upix  = ipix

c     search blueward for the beginning of the feature

         DO 13 WHILE ((ratio.lt.-1.0).AND.(lpix.gt.1))
           lpix  = lpix - 1 
           ratio = ewpix(lpix,iorder)/ewsigpix(lpix,iorder)  
 13      CONTINUE

c     reset the ratio and pixel to the scanning start point

         ratio = ewpix(ipix,iorder)/ewsigpix(ipix,iorder) 

c     search redward for the ending of the feature

         DO 14 WHILE ((ratio.lt.-1.0).AND.(upix.lt.npix))
           upix  = upix + 1 
           ratio = ewpix(upix,iorder)/ewsigpix(upix,iorder)  
 14      CONTINUE

c     front edge effects; warning

         IF (lpix.eq.1) then 
           IF (verbose) then
            WRITE(6,*) ' WARNING: feature',linei,' artificially ',
     &                'terminated start of ',order(iorder)(1:20)
           END IF
         END IF

c     back edge effects; warning

         IF (upix.eq.npix) then 
           IF (verbose) then
            WRITE(6,*) ' WARNING: feature',linei,' artificially ',
     &                'terminated end of ',order(iorder)(1:20)
           ENDIF
         END IF

c     got one- store the pixel indices and set the SL_FLAG high

         nfind(1,linei)        = 1
         f_beg(linei,iorder,1) = lpix
         f_end(linei,iorder,1) = upix
         sf_flag(iorder,linei) = .true. 
         
c     leap frog the feature for continued searching; we increment to the
c     pixel following the upper pixel of the feature and continue
c     scanning the spectrum

         ipix = upix + 1

c     no detection, increment the center pixel and continue scanning

       ELSE  

         ipix = ipix + 1

       END IF

 11   CONTINUE

      sysfeat = linei 

c     return

      RETURN
      END

c.........................................................................
c

      INTEGER FUNCTION   subfeat(iorder,linei)

c
c     This function returns the number of significant absorption
c     subregions within the master region defined by LINEI and the
c     starting and ending pixels, stored in f_beg and f_end, of each
c     feature.  The method uses the equivalent width spectrum, not the
c     simple flux decrement.  Though it can be used for emission
c     features to, this subroutine is designed to scope out absorption
c     features only
c
c     A feature us found when one pixel has an EW in its pixel that is
c     Nsig beyond the EWSIG uncertainty in that pixel, i.e., EW/EWSIG <
c     -Nsig.  Then the spectrum is scanned backwards to find the blue
c     edge of the feature, which is defined when the ratio EW/EWSIG = -
c     1.  This is the location at about which the flux decrement becomes
c     consistent with the continuum.  Then the spectrum is set back to
c     the initial pixel of the detection and scanned forward for
c     EW/EWSIG = -1.  The pixel locations of the lower and upper bounds
c     are stored in F_BEG and F_END.  Then the pixel number is set to
c     the upper pixel + 1 of the just found feature and the scanning for
c     another Nsig detection is searched for.
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      implicit none

      include            'sysanal.h'

      integer            ipix,lpix,upix,pix1,npix,subregi
      integer            iorder,linei
      double precision   Nsig,ratio,v1,v2,vb,ve

      include            'sysanal.com'



      subregi    = 0
      subfeat = 0

c     set the significance level 

      Nsig  = N_sigma(iorder)

c     find the beginning and ending pixels for the master region
c     (obtained from the first ion in the list) the values PIX1 and NPIX
c     are unmodified in this routine

      pix1 = 0
      npix = 0

      vb = vel(f_beg(linei,1,1),1)
      ve = vel(f_end(linei,1,1),1)

      DO 09 ipix=1,ndata(iorder)-1
       v1 = vel(ipix,iorder)
       v2 = vel(ipix+1,iorder)
       IF ((v1.le.vb).AND.(v2.gt.vb)) pix1 = ipix + 1
       IF ((v1.le.ve).AND.(v2.gt.ve)) npix = ipix
 09   CONTINUE

c     there are now two possibilities if we failed to grab PIX1 and/or NPIX

c     condition 1: partial data in this region; find the data edges;
c     first check blue edge, then check red edge

      IF ((pix1.eq.0).AND.(npix.ne.0)) then
         pix1 = 1
C         WRITE(6,*) linei,vb,ve,ion_name(iorder)(1:10),' pix1=1'
      END IF
      IF ((npix.eq.0).AND.(pix1.ne.0)) then
         npix = ndata(iorder)
C         WRITE(6,*) linei,vb,ve,ion_name(iorder)(1:10),' npix=ndata'
      END IF

C     dubug
C
C      WRITE(6,*) linei,vb,ve,ion_name(iorder)(1:10),pix1,npix
C     &           vel(pix1,iorder),vel(npix,iorder)

c     condition 2: data do not extend over range in both blue and red
c     directions or no data in this region at all; which is it?  set
c     accordingly

      IF ((pix1.eq.0).AND.(npix.eq.0)) then
       IF ((vel(1,iorder).gt.vb).AND.(vel(ndata(iorder),iorder).lt.ve)) 
     &   then
        pix1 = 1
        npix = ndata(iorder)
       END IF 

c     no data?; we are too red

       IF (vel(ndata(iorder),iorder).lt.vb) then
        subfeat = 0
        RETURN
       END IF 

c     no data?; we are too blue

       IF (vel(1,iorder).gt.ve) then
        subfeat = 0
        RETURN
       END IF 


      END IF

c     we scan the spectrum over the projected master region and find the
c     subregions for this ion transitions; ipix is the pixel meeting the
c     detection criterion; then we scan blueward and then redward
c     looking for the extremities of the feature; we use a while loop
c     because we need to leap frog the feature once we have it defined,
c     thus we do not increment ipix uniformily

      ipix     = pix1

      DO 11 WHILE (ipix.lt.npix)

       ratio = ewpix(ipix,iorder)/ewsigpix(ipix,iorder)
       
c     enter scan feature block if we have a significant equivalent
c     width; increment the feature index and then scan for the feature
c     extremities

       IF (ratio.le.-Nsig) then 

         subregi = subregi + 1
         lpix    = ipix
         upix    = ipix

c     search blueward for the beginning of the feature; do not scan
c     outside the region (PIX1)

         DO 13 WHILE ((ratio.lt.-1.0).AND.(lpix.gt.pix1))
           lpix  = lpix - 1 
           ratio = ewpix(lpix,iorder)/ewsigpix(lpix,iorder)  
 13      CONTINUE

c     reset the ratio and pixel to the scanning start point

         ratio = ewpix(ipix,iorder)/ewsigpix(ipix,iorder) 

c     search redward for the ending of the feature; do not scan
c     outside the region (NPIX)

         DO 14 WHILE ((ratio.lt.-1.0).AND.(upix.lt.npix))
           upix  = upix + 1 
           ratio = ewpix(upix,iorder)/ewsigpix(upix,iorder)  
 14      CONTINUE

c     front edge effects; warning

C         IF (lpix.eq.pix1) then 
C           WRITE(6,*) ' WARNING: subfeature',subregi,' artificially ',
C     &                '-terminated start of region ',linei,
C     &                '-for ',order(iorder)(1:20)
C         END IF

c     back edge effects; warning

C         IF (upix.eq.npix) then 
C           WRITE(6,*) ' WARNING: subfeature',subregi,' artificially ',
C     &                '-terminated end of region ',linei,
C     &                '-for ',order(iorder)(1:20)
C         END IF

c     got one- store the pixel indices and set the SL_FLAG high

         f_beg(linei,iorder,subregi) = lpix
         f_end(linei,iorder,subregi) = upix
         sf_flag(iorder,linei) = .true.
         
c     leap frog the feature for continued searching; we increment to the
c     pixel following the upper pixel of the feature and continue
c     scanning the spectrum

         ipix = upix + 1

c     no detection, increment the center pixel and continue scanning

       ELSE  

         ipix = ipix + 1

       END IF

 11   CONTINUE

c     if we have no detection for this ion transition withinin this
c     master region then we set the region size to that of the master
c     region and set the SF_FLAG low; this will allow tracking of
c     limits for further computations


      IF (subregi.eq.0) then
       subregi = 1
       f_beg(linei,iorder,subregi) = pix1
       f_end(linei,iorder,subregi) = npix
       sf_flag(iorder,linei) = .false.
      END IF

      subfeat = subregi

c     return

      RETURN
      END

c.........................................................................
c

      SUBROUTINE spurious(k,action,flag)

c     eliminate spurious master regions; this routine must be called
c     before the calculations on the data or there will be trouble
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      implicit none

      include            'sysanal.h'

      integer            i,j,k,action
      logical            flag
      double precision   total,knt

      include            'sysanal.com'


      flag   = .false.
      total  = 0
      knt    = 0

c     check only the second member of the doublet t
c     NOTE: FIX so that we check the sister tied trani

       i = 2

       DO 13 j=1,nfind(i,k)

        total = total + 1
        IF (.not.sf_flag(i,k)) knt = knt + 1

 13    CONTINUE
 

c     if flag is high and action = 1, then delete; if action is low,
c     return flag high

       IF (total.eq.knt) flag = .true.
       IF ((action.eq.1).AND.(flag)) CALL delregion(k)


      RETURN
      END


c.........................................................................
c

      SUBROUTINE pullregions

c     the user has specified that they want to pull the regions from a
c     previous SYSANAL session (or they made their own .ewreg files
c     (complete with a header line).  this routine pulls the regions and
c     subregion from all ions in the input list (default="ions.table").
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include            'sysanal.h'

      logical            error,sfl
      integer            i,iline,isub,iorder,action,ireg
      integer            maxrecord,lastiline,isubmax
      double precision   dum,w1,w2,v1,v2,wav0
      character*80       ew_file,header,chsf

      include            'sysanal.com'


c     set the maximum record length that a file can be, MXLIn is defined
c     in the "sysanal.h" file

      maxrecord = mxlin*mxlin

c     we first pull the primary ion (ORDER) to get the master regions of
c     the system; we use the velocity data and locate the pixel number
c     of the velocity and stuff the arrays F_BEG and F_END we use
c     routine GETPIXEL (see comments therein)

      error  = .false.
      iorder = 1

      CALL fappend(order(iorder),'ewreg',ew_file)
      OPEN(unit=1,file=ew_file,status='old')
      READ(1,*) header
      DO 11 i=1,maxrecord
       READ(1,710,end=12) iline,isub,dum,w1,w2,v1,v2,wav0
       sf_flag(iorder,iline) = .true.
       nfind(iorder,iline) = 1  ! by definition
       CALL getpixel(iline,isub,iorder,3,v1,error)
       CALL getpixel(iline,isub,iorder,4,v2,error)
 11   CONTINUE
 12   CLOSE(unit=1)

c     set the number of master regions

      nlines          = i - 1

c     now we do the remaining ions (ORDERS), but only if there are
c     higher order transitions, otherwise we are done

      IF (norders.eq.1) RETURN

c     first, we zero all significant feature flag (SF_FLAG) and the
c     NFIND book keeping arrays; this is so we can "catch"
c     non-detections

      DO 31 iorder=2,norders
       DO 30 iline=1,nlines
         sf_flag(iorder,iline) = .false.
         nfind(iorder,iline) = 0
 30    CONTINUE
 31   CONTINUE


      DO 21 iorder=2,norders

       lastiline       = 1 ! we need to track the ILINE index
       isubmax         = 1 ! default if no subregions

       CALL fappend(order(iorder),'ewreg',ew_file)
       OPEN(unit=1,file=ew_file,status='old')
       READ(1,*) header

       iline = 1
       isub  = 1

c     this loop reads in the detection regions for the higher order
c     ions, but if there is no entry for a given ILINE, then we need to
c     populate a non-detection region

       DO 19 i=1,maxrecord
        READ(1,710,end=22) iline,isub,dum,w1,w2,v1,v2,wav0,sfl
        sf_flag(iorder,iline) = sfl
        CALL getpixel(iline,isub,iorder,3,v1,error)
        CALL getpixel(iline,isub,iorder,4,v2,error)
        isubmax = max(isubmax,isub)
        IF (iline.ne.lastiline) then
          nfind(iorder,lastiline) = isubmax
          lastiline = iline
          isubmax   = 1
        END IF
 19    CONTINUE ! next record in IORDER

 22    CLOSE(unit=1) 

c     catch the last one in the list
      nfind(iorder,iline) = isub   

c     if there are no features, then we set SF_FLAG=.FALSE. and use the
c     velocities of the master regions (for computing limits), we need
c     to loop over the master regions of ion IORDER=1, using the
c     velocities of those regions to obtain the pixel numbers in the
c     current ion (IORDER); this duplicates similar logic for finding
c     these regions in integer function SUBFEAT

        DO 25 ireg=1,nlines
         IF (.not.sf_flag(iorder,ireg)) then 
          nfind(iorder,ireg) = 1      ! by definition 
          v1 = vel(f_beg(ireg,1,1),1) ! we are using the primary ion
          v2 = vel(f_end(ireg,1,1),1) ! we are using the primary ion
          CALL getpixel(ireg,1,iorder,1,v1,error)
          CALL getpixel(ireg,1,iorder,2,v2,error)
         END IF
 25     CONTINUE

 21   CONTINUE  ! next IORDER (ion)

C DEBUG

      WRITE(6,640) 
      DO 9002 i=1,norders
      WRITE(6,641) order(i),N_sigma(i)
      WRITE(6,642)
       DO 9003 iline=1,nlines
       IF (nfind(i,iline).ne.0) then 
        DO 9004 isub=1,nfind(i,iline)
         chsf = 'limit'
         IF (sf_flag(i,iline)) chsf = 'feature'
         WRITE(6,643) iline,isub,
     &                vel(f_beg(iline,i,isub),i),
     &                vel(f_end(iline,i,isub),i),chsf
 9004    CONTINUE
        ELSE
         WRITE(6,644)
        END IF
 9003  CONTINUE
 9002 CONTINUE
 640  FORMAT(1x,'----- USER SUPPLIED SUBREGIONS -----')
 641  FORMAT(1x,' Ion/Transition: ',a12,' (Nsig ',f4.1,')')
 642  FORMAT(1x,t3,'Reg',t8,'Sub',t18,'vel1',t29,'vel2')
 643  FORMAT(1x,i4,i5,2f11.4,3x,a20)
 644  FORMAT(1x,' No Detection')

 710  FORMAT(1x,i4,i5,f6.2,2x,5f12.4,2x,l2)

      RETURN
      END

c.........................................................................
c

      SUBROUTINE getpixel(iline,isub,iorder,action,vpos,error)

c     this routine finds the pixel of ion IORDER corresponding to the
c     velocity VPOS.  The master region index is ILINE, the subregion
c     index is ISUB.  The ACTION variable dictates the exact function to
c     be carried out by the routine.
c
c     If ACTION=1, then it is assumed that we are on the blue edge and
c     we stuff the pixel number in the array F_BEG(ILINE,IORDER,ISUB)
c
c     If ACTION=2, then it is assumed that we are on the red edge and
c     we stuff the pixel number in the array F_END(ILINE,IORDER,ISUB)
c
c     for ACTION=1 and =2 this routine is used during interactive
c     positioning and automatic finding, and as such scans for the VPOS
c     location between adjacent pixel in case the VPOS value does not
c     precisely equal the velocity value in a given pixel
c
c     ACTION=3, same as ACTION=1 except that it is assumed the VPOS
c     value equals a pixel value and thus we scan for the exact VPOS value 
c
c     ACTION=4, same as ACTION=2 except that it is assumed the VPOS
c     value equals a pixel value and thus we scan for the exact VPOS value 
c
c     for ACTION=3 and =4 this routine is used when the user is reading
c     in the .ewreg files from a previous session and thus the VPOS
c     value will exactly equal the value in one of the pixels
c
c     If no pixel is located the ERROR=.TRUE. and no action is taken, no
c     warning is issued, we simply return
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include            'sysanal.h'

      logical            error
      integer            i,iline,isub,iorder,action,ipix
      double precision   vpos,va,vb,eps

      include            'sysanal.com'

c     For ACTION=3 and =4, the use of EPS is because we had problems
c     with the .EQ. function in the IF statement - probably because of
c     inaccuracies in the storage of the variables VPOS and
c     VEL(PIX,ORDER) in some decimal place on the order of 10**-8 or
c     smaller (this was a pain in the ass work around)

      eps   = 1.0d-3

c     zero the pixel index IPIX and null logical ERROR

      ipix  = 0
      error = .false.

c     if ACTION=1 or =2

      IF ((action.eq.1).OR.(action.eq.2)) then
       DO 11 i=1,ndata(iorder)-1
        IF ((vpos.gt.vel(i,iorder)).AND.
     &      (vpos.le.vel(i+1,iorder))) ipix = i
 11    CONTINUE
      END IF

c     if ACTION=3 or =4

      IF ((action.eq.3).OR.(action.eq.4)) then
       DO 21 i=1,ndata(iorder)-1
        va = vel(i,iorder)-eps
        vb = vel(i,iorder)+eps
        IF ((vpos.gt.va).AND.(vpos.lt.vb)) ipix = i
 21    CONTINUE
      END IF

c     if we failed to locate the pixel at VPOS then set the window
c     within 1 pixel of the the extremes of the data segment, return

      IF (ipix.eq.0) then
        error = .true.
        f_beg(iline,iorder,isub) = 2
        f_end(iline,iorder,isub) = ndata(iorder)-1
        RETURN
      END IF

c     we found IPIX, perform the appropriate ACTION

      IF ((action.eq.1).OR.(action.eq.3)) then ! we are on the blue edge
        f_beg(iline,iorder,isub) = ipix  
        RETURN
      END IF

      IF ((action.eq.2).OR.(action.eq.4)) then ! we are on the red edge
        f_end(iline,iorder,isub) = ipix  
        RETURN
      END IF      

c     if we get here, and we shouldn't then ACTION was set incorrectly
c     or the (IPIX=0) IF block failed. communicate

      WRITE(6,*) ' ERROR(getpix): no action = ',action

      RETURN
      END


c.........................................................................
c

      SUBROUTINE delregion(linei)

c
c     the user has chosen to remove a master region; telescope the book
c     keeping arrays for the removal of the master region, LINEI passes
c     the number of the region to remove
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


      implicit none

      include            'sysanal.h'

      integer            i,j,k,linei

      include            'sysanal.com'


c     if this is the only region, then what it the point?!! we must have
c     monkey at the keyboard! do not allow deltion if this is the only
c     region

      IF (nlines.eq.1) then
       WRITE(6,*) ' *********************************************'
       WRITE(6,*) ' * You cannot delete the one and only region *'
       WRITE(6,*) ' *********************************************'
       RETURN
      END IF

c     if this is the last region in the list, simply reduce the number
c     of regions

      IF (linei.eq.nlines) then
       nlines = nlines -1
       RETURN
      END IF

c     loop over the ions (ORDERS) and the regions (LINES), we are
c     removing the line stored in LINEI, so we start at LINEI+1 and
c     reduce the indices of all higher regions, also carrying along the
c     subregions (if any) and the SF_FLAG values

       DO 09 i=1,norders
        DO 11 k=linei+1,nlines
         DO 13 j=1,nfind(i,k)
          f_beg(k-1,i,j) = f_beg(k,i,j)
          f_end(k-1,i,j) = f_end(k,i,j)
 13             CONTINUE
        nfind(i,k-1)   = nfind(i,k)
        sf_flag(i,k-1) = sf_flag(i,k)
 11          CONTINUE
 09              CONTINUE

c     now decrement the number of regions (NLINES) and return

       nlines = nlines - 1

       RETURN
       END

c     eof


