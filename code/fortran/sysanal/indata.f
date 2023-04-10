c.........................................................................
c

      SUBROUTINE getpars(inp_file,error)

c
c     read in the sysanal.inp file and get the run time parameters; set
c     some COMMON variables
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c


      implicit none

      include            'sysanal.h'

      integer            i,interact,zuser,ewreguser
      integer            istat,access
      logical            error
      double precision   R_isf
      character*80       def_file,inp_file,ch_dum

      include            'sysanal.com'




      zfix     = .false.
      ewreg    = .false.

      def_file = inp_file
      IF (def_file.eq.'') def_file = 'sysanal.inp'
      istat    = access(def_file,'r')  ! returns 1=DNE, 0="OK"

      IF (istat.ne.0) then
       def_file = '/home/matrix/cwc/Programs/Defaults/sysanal.inp'
       istat    = access(def_file,'r') 
       IF (istat.ne.0) GOTO 900
      END IF

      OPEN(unit=1,file=def_file,ERR=900,status='old')
      IF (verbose) WRITE(6,600) def_file

      READ(1,'(a80)') ch_dum  ! skip the header
      READ(1,*,ERR=901) ionfile
      READ(1,*,ERR=907) interact    ! non functional
      READ(1,*,ERR=907) Slev
      READ(1,*,ERR=906) R_isf
      READ(1,*,ERR=908) toler
      READ(1,*,ERR=911) plt_factor
      READ(1,*,ERR=912) zuser
      READ(1,*,ERR=913) ewreguser

      CLOSE(unit=1)

c     if ZUSER=1 then the user is supplying the redshift of the system
c     in the "zabs.dat" file, set ZFIX=.TRUE.

      IF (zuser.eq.1)     zfix = .true.

c     if EWREGUSER=1 then the user is supplying the wavelength regions
c     in the file "ew_regions.dat", set EWREG=.TRUE.

      IF (ewreguser.eq.1) ewreg = .true.

c     set the profile width for the ISF convolution for feature finding

      profile  = 1.0d0 / (2.35d0*R_isf) 

c     return

      RETURN

c     error traps; return error high

 900  WRITE(6,*) ' ERROR(getpars): cannot locate runtime parameter file'
      write(6,*) ' ',def_file(1:30)
      WRITE(6,*) ' cannot continue'
      error= .true.
      RETURN
 901  WRITE(6,*) ' ERROR(getpars): conversion of IONFILE failed'
      WRITE(6,*) ' cannot continue'
      error= .true.
      RETURN
 906  WRITE(6,*) ' ERROR(getpars): conversion of R_ISF failed'
      WRITE(6,*) ' cannot continue'
      error= .true.
      RETURN
 907  WRITE(6,*) ' ERROR(getpars): conversion of SLEV failed'
      WRITE(6,*) ' cannot continue'
      error= .true.
      RETURN
 908  WRITE(6,*) ' ERROR(getpars): conversion of TOLER failed'
      WRITE(6,*) ' cannot continue'
      error= .true.
      RETURN
 911  WRITE(6,*) ' ERROR(getpars): conversion of PLT_FACTOR failed'
      WRITE(6,*) ' cannot continue'
      error= .true.
      RETURN
 912  WRITE(6,*) ' ERROR(getpars): conversion of ZUSER failed'
      WRITE(6,*) ' cannot continue'
      error= .true.
      RETURN
 913  WRITE(6,*) ' ERROR(getpars): conversion of EWREGUSER failed'
      WRITE(6,*) ' cannot continue'
      error= .true.
      RETURN

c     formats

 600  FORMAT(1x,' READING run time parameters: ',a50)

      END

c.........................................................................
c

      SUBROUTINE getlist(list_file,error)

c     read in the ion/transitions (data files names)
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include            'sysanal.h'

      logical            error
      integer            i,j,k
      double precision   dum
      character*80       list_file,ch_dum

      include            'sysanal.com'



      norders = 0
      ntie    = 0

c     if list_file was not specified on the command line, we use the
c     default name "ions.table"

      IF (list_file.eq.' ') list_file = 'ions.table'

      IF (verbose) WRITE(6,500) list_file


      OPEN(unit=1,file=list_file,err=992,status='old')

      DO 51 i=1,mxord
        READ(1,*,end=101) order(i),tie(i),dum,dum,N_sigma(i)
        norders = norders + 1
        ntie    = max(ntie,tie(i))
 51   CONTINUE

c     if you are here, then you truncated the list file.  Warn and
c     continue

      IF (verbose) WRITE(6,600) list_file,order(norders)

 101  CLOSE(unit=1)

c     now set up the transition ties; ntie is the number of unique
c     species (i.e., MgII, FeII, etc); tie(order) is the tie index of
c     the tie for that species, and njtie is the number of ion
c     transitions for the tie index

      DO 43 i=1,ntie
        njtie(i) = 0
       DO 41 j=1,norders
         IF (tie(j).eq.i) njtie(i) = njtie(i) + 1
 41    CONTINUE
 43   CONTINUE

      DO 47 i=1,ntie
       k = 0
       DO 45 j=1,norders
        IF (tie(j).eq.i) THEN
         k = k + 1
         jtie(i,k) = j
        END IF
 45    CONTINUE     
 47   CONTINUE


c     return

      RETURN

c     error trap

 992  WRITE(6,*) ' ERROR(getions): cannot open list file'
      WRITE(6,*) ' ',list_file(1:40)
      error = .true.
      RETURN

c     formats

 500  FORMAT(1x,' READING list of ion transitions: ',a50)

 600  FORMAT(1x,'WARNING(getions): MXORD exceeded while',/,
     &       1x,'- reading list file',a40,/,
     &       1x,'- last successful data file is : ',a12,/,
     &       1x,'- continuing...')

      END

c.........................................................................
c

      SUBROUTINE getatomic(error)

c
c     read in the atomic data file; the name of the file containing the
c     atomic data was read in subroutine getpars
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include             'sysanal.h'

      integer             i,j
      logical             error,flag
      double precision    w,f
      character*80        ion_str

      include             'sysanal.com'



      IF (verbose) WRITE(6,600) ionfile


      DO 09 j=1,norders

        flag = .false.
         
        OPEN(unit=1,file=ionfile,ERR=999,status='old')

 01     DO 11 i=1,mxions
         READ(1,*,END=18) ion_str,w,f
         IF (order(j).eq.ion_str) then
          ion_name(j) = ion_str
          lambda0(j)  = w
          fval(j)     = f
          flag        = .true.
          GOTO 18
         END IF
 11     CONTINUE

c     if you are here you either did not find atomic data for the ion
c     transition or you went to the end of the file without finding it

        IF (.not.flag) then
         WRITE(6,*) ' ERROR(getatomic): cannot locate atomic data'
         WRITE(6,*) ' - for ',order(j)(1:30)
         WRITE(6,*) ' - perhaps the file is too long?'
         WRITE(6,*) ' - cannot continue'
         error = .true.
         RETURN
        END IF

        WRITE(6,*) ' ERROR(getatomic): file is too long- '
        WRITE(6,*) ' atomic data file truncated at ',ion_str(1:10)
        WRITE(6,*) ' - cannot continue'
        RETURN

 18     CLOSE(unit=1)

 09   CONTINUE

c     return

      RETURN

c     error trap

 999  WRITE(6,*) ' ERROR(getatomic): cannot open IONFILE -' 
      WRITE(6,*) ' ions not input; cannot continue'
      error = .true.
      RETURN

c     format

 600  FORMAT(1x,' READING atomic data: ',a50)

      END

c.........................................................................
c

      SUBROUTINE getdata(error)

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit   none
  
      include            'sysanal.h'

      integer            i,j,jj,k,kk,nline,ndum
      integer            access,istat
      parameter          (ndum = 20)
      logical            error
      double precision   dum(ndum)
      character*80       read_file
      character*200      stringline

      include            'sysanal.com'




      IF (verbose) WRITE(6,*) ' READING spectral data ...'

      DO 11 j=1,norders

       read_file = order(j)

       OPEN(unit=1,file=read_file,err=993,status='old')

       DO 15 i=1,nmx
         READ(1,'(a200)',END=16) stringline
         CALL parzeline(nline,ndum,dum,stringline)
         IF (nline.ne.6) GOTO 994
          wave(i,j)  = dum(1)
          vel(i,j)   = dum(2)
          flux(i,j)  = dum(3)
          sigma(i,j) = dum(4)
          smsig(i,j) = dum(5)
          cont(i,j)  = dum(6)
 15    CONTINUE

c     if you get here, then the data was truncated warn and continue

       IF (verbose) WRITE(6,600) order(j),vel(nmx,j)

 16    CLOSE(unit=1)
  
       ndata(j) = i - 1

c     correct zeroed sigma values

      DO i=1,ndata(j)
       IF (sigma(i,j).eq.0.) then
         IF (i.eq.1)        sigma(i,j) = sigma(i+1,j)
         IF (i.eq.ndata(j)) sigma(i,j) = sigma(i-1,j)
         IF ((i.ne.1).AND.(i.ne.ndata(j))) then
           sigma(i,j) = 0.50d0*(sigma(i-1,j)+sigma(i+1,j))
         ENDIF
       ENDIF
      ENDDO

 11   CONTINUE


c     save .orig files if they do not exist; we check if they exist so
c     as to not over write them; if we enter this program more than
c     once, they should be written only the first time

c     stat the file; istat=0 is means it exits

      DO 17 j=1,norders

       CALL fappend(order(j),'orig',read_file)
       istat = access(read_file,'r')  
       IF (istat.ne.0) then
        OPEN(unit=2,file=read_file,err=17,status='new')
        DO 19 i=1,ndata(j)
         WRITE(2,200) wave(i,j),vel(i,j),flux(i,j),
     &                sigma(i,j),smsig(i,j),cont(i,j)
 19     CONTINUE
        CLOSE(unit=2)
       END IF

 17   CONTINUE

c     return

      RETURN

c     error traps

 993  WRITE(6,*) ' ERROR(getdata): data file ',order(j)(1:20)
      WRITE(6,*) ' -not found or cannot be opened'
      WRITE(6,*) ' -cannot continue'
      error = .true.
      RETURN

 994  WRITE(6,*) ' ERROR(getdata): data are not 6 columns'
      WRITE(6,*) ' -cannot continue'
      error = .true.
      RETURN

c     formats

 200  FORMAT(1x,6f11.4)

 600  FORMAT(1x,'WARNING(getdata): NMAX exceeded in ',a12,/,
     &       1x,'-no data redward of velocity ',f10.2,/,
     &       1x,'-continuing...')
      
      END

Cc.........................................................................
Cc NOT FUNCTIONAL (fix someday!)
C
C      SUBROUTINE getmasks(error)
C
Cc
Cc:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Cc
C
C      implicit   none
C  
C      include            'sysanal.h'
C
C      integer            i,j,jj,k,kk,nline,ndum
C      integer            access,istat
C      parameter          (ndum = 20)
C      logical            error
C      double precision   dum(ndum)
C      character*80       read_file
C      character*200      stringline
C
C      include            'sysanal.com'
C
C
C
C
C      WRITE(6,*) ' CKECKING for pixel masks ...'
C
C 01   DO 11 j=1,norders
C
C       CALL fappend(order(j),'mask',read_file)
C
Cc     stat the file; istat=0 is means it exits
C
C       istat = access(read_file,'r')  
C
Cc     if the mask file exists, then read it in
C
C       IF (istat.eq.0) THEN
C
C        OPEN(unit=1,file=read_file,status='old')
C
C        DO 15 i=1,ndata(j)
C         READ(1,'(a200)',END=16) stringline
C         CALL parzeline(nline,ndum,dum,stringline)
C         IF (nline.ne.2) GOTO 994
C         IF (wave(i,j).ne.dum(1)) GOTO 995
C         mask(i,j) = int(dum(2))
C 15     CONTINUE
C
Cc     if you get here, then the data length is not correct
C
C        GOTO 996
C
C 16     CLOSE(unit=1)
C  
C       END IF
C
C 11   CONTINUE
C
Cc     this is a clean return
C
C      RETURN
C
Cc     error traps
C
C 994  WRITE(6,*) ' ERROR(getmasks): mask data are not 2 columns'
C      WRITE(6,*) ' -cannot continue'
C      error = .true.
C      RETURN
C
C 995  WRITE(6,*) ' ERROR(getmasks): wavelength mismatch in mask'
C      WRITE(6,*) ' -data.  bad mask data?'
C      WRITE(6,*) ' -cannot continue'
C      error = .true.
C      RETURN
C
C 996  WRITE(6,*) ' ERROR(getmasks): mask data pixel number mismatch'
C      WRITE(6,*) ' -cannot continue'
C      error = .true.
C      RETURN
C
C
C      END


c     eof


