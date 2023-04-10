

      integer            ndata,nbin

      double precision   wave,vel,flux,sigma,smsig,cont,lvel,uvel
      double precision   wbin,vbin,fbin,sbin,cbin
      double precision   lbin,ubin


      COMMON /iiiblk/ndata,nbin

      COMMON /rrrblk/wave(mxpix),vel(mxpix),flux(mxpix),
     &            sigma(mxpix),smsig(mxpix),cont(mxpix),
     &              lvel(mxpix),uvel(mxpix),vbin(mxpix),
     &              fbin(mxpix),sbin(mxpix),cbin(mxpix),
     &              wbin(mxpix),lbin(mxpix),ubin(mxpix)



