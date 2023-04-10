c.........................................................................
c

      SUBROUTINE speciesnames

c
c     construct the species names
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include            'sysanal.h'

      integer            i,j,idx,iorder

      include            'sysanal.com'


c     construct the species names by tie number; we only need to pull
c     the name off of the first member of the tie family (assumes user
c     has set the ties correctly)

       DO 07 i=1,ntie

         iorder = jtie(i,1)

c     loop over the 3rd to 8th characters in the ion transition name and
c     find the index marking the end of theionization level

        DO 09 j=3,8
         IF ((order(iorder)(j:j).ne.'I').AND.
     @       (spec_name(iorder)(j:j).ne.'V').AND.
     @       (spec_name(iorder)(j:j).ne.'X')) THEN
          idx = j - 1
          GOTO 08
         END IF
 09     CONTINUE

c     if we could not locate the position inthe string then warn and set
c     to ???

        WRITE(6,*) ' ERROR(seciesnames): cannot locate name for tie =',i
        spec_name(i) = '???'
        GOTO 07

c     if we jump to here, then set the species name

 08     spec_name(i) = order(iorder)(1:idx)

        IF (verbose) WRITE(6,*) i,spec_name(i)(1:8)

c     next tie value

 07    CONTINUE

c     return

      RETURN
      END

