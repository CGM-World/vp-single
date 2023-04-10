c     configure the first bin
      dw      = wlow/(resel*R)
      lbin(1) = wlow - 0.5d0*dw

      j = 1
      DO WHILE (lbin(j).lt.whigh)
       dw        = lbin(j)/(resel*R) 
       lbin(j+1) = lbin(j) + dw  
       wbin(j)   = 0.50d0*(lbin(j)+lbin(j+1))
       WRITE(7,*) j,lbin(j),wbin(j),dw        ! visual confirmation
       j = j + 1
      ENDDO

