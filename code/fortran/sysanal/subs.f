c.........................................................................
c

      SUBROUTINE comm_inputs

c
c     communicate intitial inputs
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include            'sysanal.h'

      integer            i

      include            'sysanal.com'


      IF (verbose) WRITE(6,602) 1.0d0/(2.35d0*profile)
      IF (verbose) WRITE(6,610)

      DO 01 i=1,norders
       IF (verbose) then
       WRITE(6,611) ion_name(i),tie(i),wave(1,i),wave(ndata(i),i),
     &              vel(1,i),vel(ndata(i),i),ndata(i)
       ENDIF
 01   CONTINUE

      IF (verbose) WRITE(6,620)
      DO 02 i=1,norders
       IF (verbose) WRITE(6,621) ion_name(i),lambda0(i),fval(i)
 02   CONTINUE



c formats

 602  FORMAT(1x,' Instrumental Resolution: R = ',f8.1)
 610  FORMAT(1x,'---------- PIXEL DATA ----------',/,
     @       t3,'Ion',t13,'Tie',t20,'lam-',t29,'lam+',t38,'vel-',
     @      t47,'vel+',t56,'Npix')
 611  format(1x,a10,i4,f9.2,f9.2,f9.2,f9.2,i8)
 620  FORMAT(1x,'---------- ATOMIC DATA ----------',/,
     &       t3,'Ion',t16,'lambda',t27,'f-val')
 621  FORMAT(1x,a10,f12.4,f8.4)

c     return

      RETURN
      END


c.........................................................................
c

      SUBROUTINE comm_regs

c
c     communicate intitial inputs
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include            'sysanal.h'

      integer            i

      include            'sysanal.com'


      IF (verbose) WRITE(6,630) order(1),N_sigma(1)
      DO 1001 i=1,nlines
       IF (verbose) then
       WRITE(6,631) i,
     &              wave(f_beg(i,1,1),1)/lambda0(1)-1.0d0,
     &              wave(f_end(i,1,1),1)/lambda0(1)-1.0d0,
     &              wave(f_beg(i,1,1),1),wave(f_end(i,1,1),1),
     &              vel(f_beg(i,1,1),1),vel(f_end(i,1,1),1)
       ENDIF

 1001 CONTINUE

 630  FORMAT(1x,'------- SYSTEM REGIONS ---------',/,
     &        1x,' Ion/Transition: ',a12,' (Nsig ',f4.1,')',/,
     &        1x,t3,'Reg',t11,'z-',t22,'z+',t33,'lam-',
     &        t44,'lam+',t55,'vel-',t66,'vel+')
 631  FORMAT(1x,i4,2f11.6,4f11.3)


c     If using the users input regions, then communicate
c     write out the F_BEG and F_END subregions of the weaker transition
c     in the system

C      WRITE(6,640) 
C      DO 9002 i=2,norders
C      WRITE(6,641) order(i),N_sigma(i)
C      WRITE(6,642)
C       DO 9003 linei=1,nlines
C       IF (nfind(i,linei).ne.0) then 
C        DO 9004 subregi=1,nfind(i,linei)
c         WRITE(6,643) linei,subregi,
C     &                wave(f_beg(linei,i,subregi),i)/lambda0(i)-1.0d0,
C     &                wave(f_end(linei,i,subregi),i)/lambda0(i)-1.0d0
C 9004    CONTINUE
C        ELSE
C         WRITE(6,644)
C        END IF
C 9003  CONTINUE
C 9002 CONTINUE

C 640   FORMAT(1x,'----- TRANSITION SUBREGIONS -----')
C 641   FORMAT(1x,' Ion/Transition: ',a12,' (Nsig ',f4.1,')')
C 642   FORMAT(1x,t3,'Reg',t8,'Sub',t17,'z-',t29,'z+')
C 643   FORMAT(1x,i4,i5,f12.6,f12.6)
C 644   FORMAT(1x,' No Detection')


c     return

      RETURN
      END


c.........................................................................
c

      SUBROUTINE writefiles

c
c     write all the files for all ion transitions
c
c     unit 1: the total equivalent widths for each ion transition: there
c     is 1 file for the totals (sysanal.ews); these are science ready
c     numbers
c
c     unit 10: the total AOD optical depths and column densities for
c     each ion transition: there is 1 file for the totals (sysanal.aod);
c     these are science ready numbers
c
c     unit 2: the equivalent widths for each feature (defined by the
c     first feature in the list); there is 1 file per ion transition
c     (*.ews); these are science ready numbers
c
c     unit 3: the apparent optical depth and column density spectra;
c     there is 1 file per ion transition (*.aod); these are science
c     ready numbers; these files are primarily for plotting and further
c     computation to determine the total AOD columns for an ion
c
c     unit 4: the ewlimit spectra; there is 1 file per ion transition
c     (*.ewlim); these are science ready numbers, though there is no real
c     practical use for them other than plotting for inspection
c
c     unit 7: the detection region files: there is 1 file per ion
c     transition (*.ewreg); these are final regions
c
c     unit 8: velocity shifted data file; there is 1 file per ion
c     transition; we backed up these spectra when we read in the data in
c     routine getdata; however that is only done if .orig files don't
c     already exist
c
c     unit 9: the zabs.dat file (Needed for IVPFIT and MINFIT); there is
c     1 file with one line
c
c     unit 11: the ew_regions.dat file (Needed for IVPFIT and MINFIT);
c     there is 1 file with one line; if one exists, then we save it as
c     ew_regions.last and write the new file
c
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit   none
  
      include            'const.dek'
      include            'sysanal.h'

      integer            i,j,k
      integer            istat,access
      double precision   dcu,dcd,lcol,r,w0,w1,w2,w3,v1,v2
      character*80       write_file

      include            'sysanal.com'




c     set the regions to that of the first ion transition in the list

      OPEN(unit=1,file='sysanal.ews',status='unknown')
      OPEN(unit=10,file='sysanal.aod',status='unknown')
      OPEN(unit=12,file='sysanal.vmom',status='unknown')


c     stamp headers for syanal.* files

      WRITE(1,100)     ! sysanal.ews
      WRITE(10,1005)   ! sysanal.aod
      WRITE(12,1000)   ! sysanal.vmom


      DO 11 j=1,norders

c     total measured EWs for each ion transition: 
c     FILE = "sysanal.ews"

       w1 = ewtot(j)/(1.0d0+zbar)
       w2 = ewsigtot(j)/(1.0d0+zbar)
       IF (ewsigtot(j).eq.-1.0) w2 = -1.0d0
       WRITE(1,110) j,ewtot(j),ewsigtot(j),w1,w2,
     &              drtot(j),drsigtot(j),sltot(j),order(j)

c     total measured optical depths and AOD column densities for each
c     ion transition; the errors are not symmetric (because of the
c     natural log), so there is a down error (dd) and an up error (du)
c     for each of the totals: 
c     FILE = "sysanal.aod"

       WRITE(10,1010) j,tautot(j),ddtautot(j),dutautot(j),
     &                coltot(j),ddcoltot(j),ducoltot(j),
     &                order(j)

c     the total velocity moments; i.e., the mean velocity, velocity
c     width, and the asymmetry; we convert the asymmetry to asym/vew,
c     where vew = c*ew/lambda: 
c     FILE = "sysanal.vmom"

       w0 = ewtot(j)/(1.0d0+zbar)
       w1 = vtotasym(j)/(ckms*w0/lambda0(j))
       w2 = sigvtotasym(j)/(ckms*w0/lambda0(j))
       WRITE(12,1010) j,vtotbar(j),sigvtotbar(j),vtotwidth(j),
     &                sigvtotwidth(j),w1,w2,order(j)

c     measured EW per region for each region
c     FILE = order(j)".ews"

       CALL fappend(order(j),'ews',write_file)
       OPEN(unit=2,file=write_file,status='unknown')
       WRITE(2,200)    ! stamp header 
       DO 15 i=1,nlines
        w1 = ew(i,j)/(1.0d0+zbar)
        w2 = ewsig(i,j)/(1.0d0+zbar)
        IF (ewsig(i,j).eq.-1.0) w2 = -1.0d0
        WRITE(2,210) i,wbar(i,j),ew(i,j),ewsig(i,j),w1,w2,
     &               dr(i,j),drsig(i,j),siglevel(i,j),order(j)
 15    CONTINUE
       CLOSE(unit=2)

c     AOD spectra; 1 file per ion transition; 
c     FILE = order(j)".aod"

       CALL fappend(order(j),'aod',write_file)
       OPEN(unit=3,file=write_file,status='unknown')
       WRITE(3,300)    ! stamp header
       DO 17 i=1,ndata(j)
        WRITE(3,310) wave(i,j),vel(i,j),tau(i,j),ddtau(i,j),
     &               dutau(i,j),col(i,j),ddcol(i,j),ducol(i,j),
     &               0,collim(i,j)
 17   CONTINUE
      CLOSE(unit=3)

c     EW per pixel spectra; also provide the ratio, R
c     FILE = order(j)".ewlim"

       CALL fappend(order(j),'ewlim',write_file)
       OPEN(unit=4,file=write_file,status='unknown')
       DO 19 i=1,ndata(j)
        r = ewpix(i,j)/ewsigpix(i,j)  
        WRITE(4,410) wave(i,j),vel(i,j),ewpix(i,j),ewsigpix(i,j),
     &           -ewsigpix(i,j),r,N_sigma(j)*ewsigpix(i,j),
     &           -N_sigma(j)*ewsigpix(i,j)
 19    CONTINUE
       CLOSE(unit=4)

c     EW region files; use the regions defined by first ion transition
c     in the list; r=(1+z) of the region edges

c     REVISION, write no matter what, commented out IF block
c       ONLY WRITE IF user has not supplied the regions OR if there are
c       modifications per the MKMOD flag, which is .TRUE. if the use
c       made modifications FILE = order(j)".ewreg"

C       IF ((.not.ewreg).OR.(mkmod)) then
        CALL fappend(order(j),'ewreg',write_file)
        OPEN(unit=7,file=write_file,status='unknown')
        WRITE(7,700)    ! stamp header
        DO 22 i=1,nlines
C         IF (sf_flag(j,i)) then
          DO 21 k=1,nfind(j,i)
           w1 = wave(f_beg(i,j,k),j)
           w2 = wave(f_end(i,j,k),j)
           v1 = vel(f_beg(i,j,k),j)
           v2 = vel(f_end(i,j,k),j)
           WRITE(7,710) i,k,1.25,w1,w2,v1,v2,lambda0(j),sf_flag(j,i)
 21       CONTINUE
C         END IF
 22     CONTINUE
        CLOSE(unit=7)
C       END IF

c     The spectra themselves; post velocity shifting; note that we have
c     saved the orginal spectra in *.orig files ONLY WRITE IF user has
c     not supplied the redshift (if they have we did not modify the
c     velocity center (i.e., shift the spectra), so there is no change
c     in the files anyhow FILE = order(j)

       IF (.not.zfix) then
       write_file = order(j)
       OPEN(unit=8,file=write_file,status='unknown')
       DO 23 i=1,ndata(j)
        WRITE(8,810) wave(i,j),vel(i,j),flux(i,j),
     &               sigma(i,j),smsig(i,j),cont(i,j)
 23    CONTINUE
       CLOSE(unit=8)
      END IF

c     write the velocity moment data region by region
c     FILE = order(j)".vmom"

       CALL fappend(order(j),'vmom',write_file)
       OPEN(unit=13,file=write_file,status='unknown')
       WRITE(13,1000)
       DO 25 i=1,nlines
        IF (sf_flag(j,i)) then
         w0 = ew(i,j)/(1.0d0+zbar)
         w1 = vasym(i,j)/(ckms*w0/lambda0(j))
         w2 = sigvasym(i,j)/(ckms*w0/lambda0(j))
         WRITE(13,1010) i,vbar(i,j),sigvbar(i,j),vwidth(i,j),
     &                  sigvwidth(i,j),w1,w2,order(j)
        END IF
 25    CONTINUE
       CLOSE(unit=13)

c     next ion transition

 11   CONTINUE

      CLOSE(unit=1)
      CLOSE(unit=10)
      CLOSE(unit=12)

c     we are done looping through the ions

c     zabs.dat file; use the first ion transition in the list
c     ONLY WRITE IF user has not supplied the zbar!

      IF (.not.zfix) then
       OPEN(unit=9,file='zabs.dat',status='unknown')
       WRITE(9,910) zbar,lambda0(1)*(1.0d0+zbar),order(1)
       CLOSE(unit=9)
      END IF

c     REVISION, write no matter what, commented out IF block
c       write the ew_regions.dat file ONLY WRITE IF user has not supplied
c       the regions OR there are modifications to the primary ion per the
c       MKMOD flag, which is .TRUE. if the use made modifications

C       IF ((.not.ewreg).OR.(mkmod)) then

c     stat the file; istat=0 is means it exits

       istat = access('ew_regions.dat','r')  

       IF (istat.eq.0) then
        WRITE(6,*) ' WARNING(writefiles): ew_regions.dat file exists'
        WRITE(6,*) ' Copying current version to ew_regions.last'
        WRITE(6,*) ' before writing...'
        CALL system('/bin/mv ew_regions.dat ew_regions.last')
       END IF

       OPEN(unit=11,file='ew_regions.dat',status='unknown')
       DO 27 i=1,nlines
        w1 = wave(f_beg(i,1,1),1)
        w2 = wave(f_end(i,1,1),1)
        w3 = wbar(i,1)
        WRITE(11,1110) i,1.25,w1,w2,w3,lambda0(1)
 27    CONTINUE
       CLOSE(unit=11)

C      END IF

c     this is a clean return

      RETURN

c     formats

 100  FORMAT(1x,t3,'Idx',t13,'EW_o',t21,'Sig',t29,'EW_r',
     &       t37,'Sig',t43,'EW_rat',t53,'Sig',t61,'SL',
     &       t68,'Ion_Tran')
 110  format(1x,i4,3x,6f8.3,f7.1,4x,a20)
 200  FORMAT(1x,t3,'Idx',t9,'Lambda_o',t21,'EW_o',t29,'Sig',
     &       t37,'EW_r',t45,'Sig',t51,'EW_rat',t61,'Sig',t69,
     &       'SL',t76,'Ion_Tran')
 210  format(1x,i4,3x,7f8.3,f7.1,4x,a20)
 300  FORMAT(1x,t6,'Lambda',t19,'Vel',t30,'Tau',t41,'Sig-',
     &       t52,'Sig+',t63,'logN',t74,'Sig-',t85,'Sig+',t93,
     &       'xlim',t98,'ylim')
 310  format(1x,8f11.4,i6,i5)
 410  format(1x,f11.4,2x,f9.3,3f12.6,f11.2,1x,2f12.6)
 700  FORMAT(1x,t3,'Reg',t8,'Sub',t13,'ypos',t24,'Lambda-',
     &       t36,'Lambda+',t50,'Vel-',t62,'Vel+',t72,'Lambda0',
     &       t82,'Detect')
 710  format(1x,i4,i5,f6.2,2x,5f12.4,2x,l2)
 810  format(1x,6f11.4)
 910  format(1x,f10.7,f12.4,2x,a20)
 1000 FORMAT(1x,t3,'Idx',t14,'Vbar',t23,'Sig',t31,'Width',
     &       t41,'Sig',t50,'Asym',t59,'Sig',t67,'Ion_Tran')
 1005 FORMAT(1x,t3,'Idx',t14,'Tau',t23,'Sig-',t32,'Sig+',
     &       t41,'logN',t50,'Sig-',t59,'Sig+',t67,'Ion_Tran')
 1010 format(1x,i4,3x,6f9.3,4x,a20)
 1110 format(1x,i4,f6.2,4f12.4)

      END

c.........................................................................
c

      SUBROUTINE sortregions

c
c     here we sort the regions just in case they are not in velocity
c     order; keep track of all book keeping (whew!).  We use the
c     straigth insertion method.
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit   none
  
      include            'const.dek'
      include            'sysanal.h'

      integer            i,iorder,j,jidx,k,idx(mxlin),nf(mxlin),n
      logical            sf(mxlin)
      double precision   a(mxlin),fb(mxlin,mxlin),fe(mxlin,mxlin)

      include            'sysanal.com'


c     store the master region array for index sorting

      DO 01 i=1,nlines
        a(i) = f_beg(i,1,1)
 01   CONTINUE  

c     obtain the index order in IDX

      n = nlines
      CALL indexx(n,a,idx)

c     now, loop over each ion transition

      DO 03 iorder=1,norders

c     save the arrays in temporary storage

       DO 05 j=1,nlines
        sf(j) = sf_flag(iorder,j)
        nf(j) = nfind(iorder,j)
        DO 07 k=1,nf(j)
         fb(j,k) = f_beg(j,iorder,k)
         fe(j,k) = f_end(j,iorder,k)
 07     CONTINUE
 05    CONTINUE


c     now move to the IDX location

       DO 15 j=1,nlines
        jidx = idx(j)  
        sf_flag(iorder,j) = sf(jidx)
        nfind(iorder,j)   = nf(jidx) 
        DO 17 k=1,nf(jidx)
         f_beg(j,iorder,k) = fb(jidx,k) 
         f_end(j,iorder,k) = fe(jidx,k)
 17     CONTINUE
 15    CONTINUE

 03   CONTINUE

c     return

      RETURN
      END


c.........................................................................
c

      SUBROUTINE zeroall

c
c     zero all working arrays
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include            'sysanal.h'

      integer            j,i

      include            'sysanal.com'


       mkmod        = .false.

c     intialize all working arrays

c     loop over orders

      DO 61 i=1,mxord

       ndata(i)        = 0
       ewtot(i)        = 0.0d0
       ewsigtot(i)     = 0.0d0
       drtot(i)        = 0.0d0
       drsigtot(i)     = 0.0d0
       sltot(i)        = 0.0d0
       vtotbar(i)      = 0.0d0
       vtotwidth(i)    = 0.0d0
       vtotasym(i)     = 0.0d0
       sigvtotbar(i)   = 0.0d0
       sigvtotwidth(i) = 0.0d0
       sigvtotasym(i)  = 0.0d0
       tautot(i)       = 0.0d0
       dutautot(i)     = 0.0d0
       ddtautot(i)     = 0.0d0
       coltot(i)       = 0.0d0
       ducoltot(i)     = 0.0d0
       ddcoltot(i)     = 0.0d0

c     feature/order arrays

       DO 71 j=1,mxlin

        ew(j,i)        = 0.0d0
        ewsig(j,i)     = 0.0d0
        wbar(j,i)      = 0.0d0
        sigwbar(j,i)   = 0.0d0
        width(j,i)     = 0.0d0
        sigwidth(j,i)  = 0.0d0
        dr(j,i)        = 0.0d0
        drsig(j,i)     = 0.0d0
        siglevel(j,i)  = 0.0d0
        vbar(j,i)      = 0.0d0
        vwidth(j,i)    = 0.0d0
        vasym(j,i)     = 0.0d0
        sigvbar(j,i)   = 0.0d0
        sigvwidth(j,i) = 0.0d0
        sigvasym(j,i)  = 0.0d0

 71    CONTINUE

c     pixel/order arrays

       DO 51 j=1,nmx

        mask(j,i)     = 1      ! this array is NON FUNCTIONAL
        collim(j,i)   = 0
        wave(j,i)     = 0.0d0
        vel(j,i)      = 0.0d0
        flux(j,i)     = 0.0d0
        sigma(j,i)    = 0.0d0
        smsig(j,i)    = 0.0d0
        cont(j,i)     = 0.0d0
        ewpix(j,i)    = 0.0d0
        ewsigpix(j,i) = 0.0d0
        tau(j,i)      = 0.0d0
        dutau(j,i)    = 0.0d0
        ddtau(j,i)    = 0.0d0
        col(j,i)      = 0.0d0
        ducol(j,i)    = 0.0d0
        ddcol(j,i)    = 0.0d0

 51    CONTINUE

 61   CONTINUE

c     return

      RETURN
      END


c.........................................................................
c

      SUBROUTINE helpme

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      WRITE(6,*) ' USAGE: search [$1] [$2] [$3]'
      WRITE(6,*) ' '
      WRITE(6,*) ' [$1], [$2], [$3] = optional command line parameters'
      WRITE(6,*) ' [$1] defaults to ions.table'
      WRITE(6,*) ' [$2] defaults to sysanal.inp'
      WRITE(6,*) ' [$3] defaults to verbose+'
      WRITE(6,*) ' '
      WRITE(6,*) ' if you wish to name your files something other'
      WRITE(6,*) ' other than the default name, place them in the'
      WRITE(6,*) ' current working directory and type their names'
      WRITE(6,*) ' on the command line'

      WRITE(6,*) ' if [$1] = "help" or "?" then this help message'
      WRITE(6,*) ' is printed.'

      RETURN
      END



c     eof


