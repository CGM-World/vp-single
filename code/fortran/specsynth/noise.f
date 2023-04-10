      FUNCTION GASDEV(IDUM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER idum
      SAVE ISET,GSET
      DATA ISET/0/
      IF (ISET.EQ.0) THEN
1       V1=2.*RAN1(IDUM)-1.
        V2=2.*RAN1(IDUM)-1.
        R=V1**2+V2**2
        IF(R.GE.1.)GO TO 1
        FAC=SQRT(-2.*LOG(R)/R)
        GSET=V1*FAC
        GASDEV=V2*FAC
        ISET=1
      ELSE
        GASDEV=GSET
        ISET=0
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION ran1(idum)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     $     NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=3.e-16,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
         idum=max(-idum,1)
         do j=NTAB+8,1,-1
            k=idum/IQ
            idum=IA*(idum-k*IQ)-IR*k
            if (idum.lt.0) idum=idum+IM
            if (j.le.NTAB) iv(j)=idum
         enddo
         iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
