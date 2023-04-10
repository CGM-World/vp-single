
      integer            nmx,mxlin,mxord,mxions,mxfits,nbuff,nmxreg
      parameter          (nmx    =  10000,    ! max number of pixels in a spectrum
     &                    mxlin  =     50,    ! max number of master regions (detections)
     &                    mxord  =     30,    ! max number of transitions to analyze
     &                    mxions =    300)    ! max number of atomic data entries

      integer            M,Jnot
      parameter          (Jnot = 6, M = 2*Jnot + 1)

