
      logical            sf_flag,int_flag,zfix,ewreg,mkmod,verbose
      integer            ndata,nfind,f_beg,f_end,mask,norders,
     &                   collim,nions,nlines
      integer            ntie,tie,jtie,njtie
      double precision   fval,lambda0,zbar
      double precision   wave,vel,flux,sigma,smsig,cont,
     &                   tau,dutau,ddtau,col,ducol,ddcol,
     &                   ewpix,ewsigpix
      double precision   ew,ewsig,wbar,sigwbar,width,sigwidth,
     &                   dr,drsig,siglevel
      double precision   vbar,vwidth,vasym,vtotbar,vtotwidth,
     &                   vtotasym,sigvbar,sigvwidth,sigvasym,
     &                   sigvtotbar,sigvtotwidth,sigvtotasym
      double precision   ewtot,ewsigtot,drtot,drsigtot,
     &                   sltot,tautot,dutautot,ddtautot,
     &                   coltot,ducoltot,ddcoltot,N_sigma
      double precision   Slev,profile,toler,plt_factor
      character*80       order,ionfile,ion_name,spec_name


      COMMON /lllblk/  int_flag,zfix,ewreg,verbose,
     &                 mkmod,sf_flag(mxord,mxlin)

      COMMON /tieblk/  ntie,tie(mxord),jtie(mxord,mxord),
     &                 njtie(mxord)

      COMMON /iiiblk/  ndata(mxord),nfind(mxord,mxlin),
     &                 f_beg(mxlin,mxord,mxlin),
     &                 f_end(mxlin,mxord,mxlin),
     &                 mask(nmx,mxord),norders,
     &                 collim(nmx,mxord),nions,nlines

      COMMON /rrdblk/ wave(nmx,mxord),vel(nmx,mxord),
     &                flux(nmx,mxord),sigma(nmx,mxord),
     &                smsig(nmx,mxord),cont(nmx,mxord),
     &                ewpix(nmx,mxord),ewsigpix(nmx,mxord),
     &                tau(nmx,mxord),dutau(nmx,mxord),
     &                ddtau(nmx,mxord),col(nmx,mxord),
     &                ducol(nmx,mxord),ddcol(nmx,mxord)

      COMMON /rrlblk/ ew(mxlin,mxord),ewsig(mxlin,mxord),
     &                wbar(mxlin,mxord),sigwbar(mxlin,mxord),
     &                width(mxlin,mxord),sigwidth(mxlin,mxord),
     &                dr(mxlin,mxord),drsig(mxlin,mxord),
     &                siglevel(mxlin,mxord),vbar(mxlin,mxord),
     &                vwidth(mxlin,mxord),vasym(mxlin,mxord),
     &                sigvbar(mxlin,mxord),sigvwidth(mxlin,mxord),
     &                sigvasym(mxlin,mxord)


      COMMON /rroblk/ ewtot(mxord),ewsigtot(mxord),
     &                drtot(mxord),drsigtot(mxord),
     &                coltot(mxord),ducoltot(mxord),
     &                ddcoltot(mxord),tautot(mxord),
     &                dutautot(mxord),ddtautot(mxord),
     &                sltot(mxord),N_sigma(mxord),
     &                vtotbar(mxord),vtotwidth(mxord),
     &                vtotasym(mxord),sigvtotbar(mxord),
     &                sigvtotwidth(mxord),sigvtotasym(mxord),
     &                Slev,profile,toler,plt_factor

      COMMON /miscblk/ fval(mxord),lambda0(mxord),zbar

      COMMON /cccblk/  order(mxord),ionfile,ion_name(mxord),
     &                 spec_name(mxord)

