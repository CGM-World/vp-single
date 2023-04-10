c.........................................................................
c

      SUBROUTINE parzeline(nline,ndum,dum,stringline)

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      logical           error
      integer           ndum,nline,i,j,ii,istart,iend
      double precision  dum(ndum),value2
      character*80      tmpstr
      character*200     stringline


      DO 11 i=1,ndum
        dum(i) = 0.0d0
 11   CONTINUE

      j     = 0
      i     = 1
      nline = 0


 01    IF (stringline(i:i).ne.' ') then 

         istart = i
         error  = .false.

         DO 15 ii=istart,200
          IF (stringline(ii:ii).ne.' ') then
            iend = ii
          ELSE
            GOTO 02
          END IF
 15      CONTINUE

         nline = j
         RETURN

 02      tmpstr = stringline(istart:iend)
         j      = j + 1
         IF (j.gt.ndum) STOP '(parzeline): J>NDUM'
         dum(j) = value2(tmpstr,error)
         IF (error) dum(j) = 0.0d0
         i = iend + 1
         GOTO 01

       ELSE

         i = i + 1
         GOTO 01
 
       END IF
       

      END

c.........................................................................
c

      SUBROUTINE parzstring(inword,outword1,outword2,error)

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      logical             error
      integer             i,k1beg,k1end,k2beg,k2end,lend      
      character*(*)       inword,outword1,outword2


c     find the beginning and ending of each word assuming forms 's1 s2',
c     's1=s2', or 's1 = s2' and allow for leading blanks

      lend  = len(inword)
      k1beg = 0
      k2beg = 0
      k1end = 0
      k1end = 0


c     begin 1st word

      DO 07 i=1,lend
       IF (inword(i:i).ne.' ') then
         k1beg = i
         GOTO 08
       END IF
 07   CONTINUE

c     end first word

 08   DO 09 i=k1beg,lend
       IF ((inword(i:i).eq.' ').OR.(inword(i:i).eq.'=')) then
         k1end = i-1
         GOTO 10
       END IF
 09   CONTINUE

c     begin 2nd word

 10   DO 11 i=k1end+1,lend   
       IF ((inword(i:i).ne.' ').AND.(inword(i:i).ne.'=')) then
         k2beg = i
         GOTO 12
       END IF
 11   CONTINUE

c     end 2nd word

 12   DO 13 i=k2beg,lend
       IF (inword(i:i).eq.' ') then
         k2end = i - 1
         GOTO 14
       END IF
 13   CONTINUE

 14   IF (k1beg*k1end.eq.0) then
        error = .true.
        outword1 = ' '
        outword2 = ' '
        RETURN
      ELSE
        outword1 = inword(k1beg:k1end)
      END IF

      IF (k2beg*k2end.eq.0) then
        outword2 = ' '
      ELSE
        outword2 = inword(k2beg:k2end)
      END IF

      RETURN
      END


c.........................................................................
c

      SUBROUTINE fappend(infile,delim,outfile)   

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

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

c.........................................................................
c

      SUBROUTINE sappend(inword,tie,addword,outword)   

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      implicit none
      integer             i,k,lend      
      character*(*)       inword,addword,outword,tie

      lend = len(inword)
      k    = 0
      DO 09 i=1,lend
        k = i
        IF (inword(i:i).eq.' ') GOTO 10
 09   CONTINUE

 10   IF (tie.eq.' ') then
       outword = inword(1:k-1)//addword
      ELSE
       outword = inword(1:k-1)//tie//addword
      END IF


      RETURN
      END
