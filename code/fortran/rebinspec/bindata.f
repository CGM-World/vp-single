

      SUBROUTINE bindata

      implicit none
 
      include            'rebinspec.h'
      integer            m,i,j
      double precision   weight,f_flux,f_sig,norm
      double precision   lb,ub,lw,uw
      include            'rebinspec.com'
 

c  loop over the bins one at a time add the flux and sigma that
c  fall in each bin from each order; check order by order running
c  through pixels... very inefficient, but robust

      DO 11 m=1,nbin

        lb       = lbin(m)
        ub       = ubin(m)

        weight = 0.0d0
        f_flux = 0.0d0
        f_sig  = 0.0d0
        norm   = 0.0d0

c     in case norm=0 after scanning

        fbin(m) = 1.0d0
        sbin(m) = -1.0d0
        cbin(m) = 1.0d0

c     loop through the spectrum and rebin it

        DO 31 i=1,ndata

            IF (sigma(i).eq.-1.0) GOTO 31   ! bad pixel

c     shorthand for lower and upper pixel edges

            lw = lvel(i)
            uw = uvel(i)

c     pixel brackets lower edge of bin
            IF ((uw.gt.lb).AND.(lw.le.lb)) THEN
              weight = (uw-lb)/(uw-lw)
              f_flux = f_flux + weight*flux(i)
              f_sig  = f_sig + weight*sigma(i)
              norm   = norm + weight
            END IF
c     pixel embedded within bin
            IF ((lw.gt.lb).AND.(uw.le.ub)) THEN
              weight = 1.0d0
              f_flux = f_flux + weight*flux(i)
              f_sig  = f_sig + weight*sigma(i)
              norm   = norm + weight
            END IF
c     pixel brackets upper edge of bin
            IF ((lw.lt.ub).AND.(uw.ge.ub)) THEN
              weight = (ub-lw)/(uw-lw)
              f_flux = f_flux + weight*flux(i)
              f_sig  = f_sig + weight*sigma(i)
              norm   = norm + weight
            END IF
c     pixel extends beyond both edges of bin
            IF ((lw.lt.lb).AND.(uw.ge.ub)) THEN
              weight = (ub-lb)/(uw-lw)
              f_flux = f_flux + weight*flux(i)
              f_sig  = f_sig + weight*sigma(i)
              norm   = norm + weight
            END IF

 31       CONTINUE

          IF (norm.ne.0.0d0) THEN
            fbin(m) = f_flux/norm
            sbin(m) = f_sig/norm
            cbin(m) = 1.0d0
          END IF

 11   CONTINUE  ! next bin


c  now optimally weight the flux in each bin from each order

C      WRITE(6,*) 'optimally weighting bins...'

C      DO 41 m=1,nbin

C        bigF = 0.0d0
C        bigS = 0.0d0

C        DO 51, j=1,nspec

C          IF (s_pix(j,m).ne.0.0d0) THEN 
C            weight = 1.0d0/s_pix(j,m)**2
C            bigF   = bigF + weight*f_pix(j,m)
C            bigS   = bigS + weight
C          END IF

C 51     CONTINUE

C        IF (bigS.ne.0.0d0) THEN
C          fbin(m) = bigF/bigS
C          sbin(m) = 1.0d0/sqrt(bigS)
C          cbin(m) = 1.0d0
C        END IF

C 41   CONTINUE

      RETURN
      END

