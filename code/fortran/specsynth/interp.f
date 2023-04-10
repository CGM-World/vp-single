c.. 
c..this file contains routines to interpolate and extrapolate functions,
c..and to find table entries in an array
c..
c..
c..routine polint is a polynomial interpolator/extrapolator
c..routine polins duplicates the above; used to stop recursive programming
c..routine ratint is a rational function interpolator/extrapolator
c..routine spline computes the cubic spline derivatives of a function   
c..routine splint uses the spline derivatives to interpolate a function 
c..routine locate finds a table entry by bisection
c..routine hunt finds a table entry given a guess for a nearby entry
c..routine ssplint does a 3rd order spline interpolation of a function
c..routine polcoe evalutes the coefficients of the interpolating polynomial
c..routine polcof also evalutes the coefficients of the interpolating polynomial
c..routine polin2 does 2d polint interpolations
c..routine bcuint does 2d bicubic interpolations
c..routine bcucof evaluates the coeficients for 2d bicubic interpolations
c..routine splie2 evaluates the coefficients for 2d bicubic spline inter.
c..routine splin2 does 2d bicubic spline interpolation
c.. 
c.. 
c.. 
c..
      subroutine polint(xa,ya,n,x,y,dy)
      include 'implno.dek'
c..
c..given arrays xa and ya of length n and a value x, this routine returns a 
c..value y and an error estimate dy. if p(x) is the polynomial of degree n-1
c..such that ya = p(xa) ya then the returned value is y = p(x) 
c..
c..declare
      integer          n,nmax,ns,i,m
      parameter        (nmax=30)
      double precision xa(1),ya(1),x,y,dy,
     +                 c(nmax),d(nmax),dif,dift,ho,hp,w,den
c..
c..find the index ns of the closest table entry and initialize 
c..the c and d tables
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
c..
c..first guess for y
      y=ya(ns)
c..
c..for each column of the table, loop over the c's and d's and update them
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) stop ' 2 xa entires are the same in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
c..
c..after each column is completed, decide which correction c or d, to add
c..to the accumulating value of y, that is, which path to take in the table
c..by forking up or down. ns is updated as we go to keep track of where we
c..are. the last dy added is the error indicator.
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end
c..
c..
c..
c..
c..
      subroutine polins(xa,ya,n,x,y,dy)
      include 'implno.dek'
c..
c..given arrays xa and ya of length n and a value x, this routine returns a 
c..value y and an error estimate dy. if p(x) is the polynomial of degree n-1
c..such that ya = p(xa) ya then the returned value is y = p(x) 
c..
c..declare
      integer          n,nmax,ns,i,m
      parameter        (nmax=30)
      double precision xa(1),ya(1),x,y,dy,
     +                 c(nmax),d(nmax),dif,dift,ho,hp,w,den
c..
c..find the index ns of the closest table entry and initialize 
c..the c and d tables
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
c..
c..first guess for y
      y=ya(ns)
c..
c..for each column of the table, loop over the c's and d's and update them
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) stop ' 2 xa entires are the same in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
c..
c..after each column is completed, decide which correction c or d, to add
c..to the accumulating value of y, that is, which path to take in the table
c..by forking up or down. ns is updated as we go to keep track of where we
c..are. the last dy added is the error indicator.
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end
c..
c..
c..
c..
c..
      subroutine ratint(xa,ya,n,x,y,dy)
      include 'implno.dek'
c..
c..given arrays xa and ya of length n and a value x, this routine returns a 
c..value y and an error estimate dy. the value returned is that of the
c..diagonal rational function, evaluated at x, which passes through the
c..n points (xa,ya).
c..
c..declare
      integer          n,nmax,ns,i,m
      parameter        (nmax=100)
      double precision xa(1),ya(1),x,y,dy,
     +                 c(nmax),d(nmax),tiny,h,hh,t,w,dd
      parameter        (tiny=1.e-25)
c..
c..initialize , tiny is used to prevent a rare zero over zero condition
      ns=1
      hh=abs(x-xa(1))
      do 11 i=1,n
        h=abs(x-xa(i))
        if (h.eq.0.)then
          y=ya(i)
          dy=0.0
          return
        else if (h.lt.hh) then
          ns=i
          hh=h
        endif
        c(i)=ya(i)
        d(i)=ya(i)+tiny
11    continue
c..
c..first guess for y
      y=ya(ns)
c..
c..for each column of the table, loop over the c's and d's and update them
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          w=c(i+1)-d(i)
          h=xa(i+m)-x
          t=(xa(i)-x)*d(i)/h
          dd=t-c(i+1)
          if(dd.eq.0.) stop 'pole at the value of x in ratint'
          dd=w/dd
          d(i)=c(i+1)*dd
          c(i)=t*dd
12      continue
c..
c..after each column is completed, decide which correction c or d, to add
c..to the accumulating value of y, that is, which path to take in the table
c..by forking up or down. ns is updated as we go to keep track of where we
c..are. the last dy added is the error indicator.
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end
c..
c..
c..
c..
c..
      subroutine spline(x,y,n,y2)   
      include 'implno.dek'
c.. 
c..given arrays x and y of length n containing a tabulated function i.e y=f(x), 
c..with the x monotonically increasing and given values for yp1 and ypn, the
c..first derivative at the pints 1 and n respectively, this routine returns 
c..the array y2 of length n which contains the second derivatives of the    
c..at the tabulated points x.   
c..if yp1 and/or ypn are set larger than 1e30, the routine sets the boundary
c..condtion for a natural spline (one with zero second derivative) at the   
c..boundary.
c..this routine is only called once for any given x and y.  
c.. 
c..declare  
      integer            n,i,k,nmax 
      parameter          (nmax=10000) 
      double precision   x(1),y(1),y2(1),yp1,ypn,u(nmax),sig,p,qn,un
      parameter          (yp1=1.0e33, ypn=1.0e33)
c.. 
c.. 
c..the lower boundary condition is set to be either natural 
      if (yp1 .gt. .99E30) then 
        y2(1)=0.0   
        u(1)=0.0
c..or it has a specified first derivative   
      else  
        y2(1)=-0.5  
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1) 
      end if
c.. 
c..this is the decomposition loop of the tridiagonal algorithm. y2 and u
c..are used for temporary storage of the decomposition factors. 
      do 11 i=2,n-1 
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))   
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) 
     +      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p   
11    continue  
c.. 
c..the upper boundary condition is set to be either natural 
      if (ypn .gt. .99E30) then 
        qn=0.0  
        un=0.0  
c.. 
c..or it has a specified first derivative   
      else  
        qn=0.5  
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1))) 
      end if
c.. 
c..this is the backsubstitution loop of the tridiagonal algorithm   
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)  
      do 12 k=n-1,1,-1  
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue  
c.. 
c..bye bye  
      return
      end   
c.. 
c.. 
c.. 
c.. 
c.. 
c.. 
      subroutine splint(xa,ya,y2a,n,x,y)
      include 'implno.dek'
c.. 
c..given the arrays xa and ya of length n, which tablulate a monotonic function,
c..and given the array y2a, which is the output of spline (above), and given
c..a value of x this routine returns a cubic spline interpolated value of y.
c.. 
c..declare  
      integer             n,klo,khi,k   
      double precision    xa(1),ya(1),y2a(1),x,y,h,a,b  
c.. 
c.. 
c..find the right place in the table by bisection. this is optimal if   
c..the sequential calls to this routine are at random values of x.  
c..if the sequential calls are in order and closely spaced, one might store 
c..the values of klo and khi and test if they remain appropriate on next call   
      klo=1 
      khi=n 
1     if (khi-klo .gt. 1) then  
       k=(khi+klo)/2
       if (xa(k) .gt. x) then   
        khi=k   
       else 
        klo=k   
       end if   
       goto 1   
      end if    
c.. 
c..klo and khi now bracket the input value of x 
c..the xa's must be distinct
      h=xa(khi)-xa(klo) 
      if (h .eq. 0.0d0) then
        write(6,99) khi,xa(khi),klo,xa(klo)
   99   format(1x,'  khi=',i4,'  xa(khi)=',1pe13.5,'  klo=',
     1         i4,'  xa(klo)=',1pe13.5)
        stop 'bad xa input in routine splint' 
      end if
c.. 
c..evaluate the cubic spline
      a=(xa(khi)-x)/h   
      b=(x-xa(klo))/h   
      y=a*ya(klo)+b*ya(khi)+
     +      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0
c.. 
c..bye bye  
      return    
      end   
c..
c..
c..
c..
c..
      subroutine locate(xx,n,x,j)
      include 'implno.dek'
c..
c..given an array xx of length n, and a value of x, this routine returns
c..a value j such that x is between xx(j) and xx(j+1). that is the value
c..in the array that is lower and closest to the value of x.
c..the array xx must be monotonic.
c..j=0 or j=n indicates that x is out of range
c..bisection  is used to find the entry
c..
c..declare
      integer           n,j,jl,ju,jm
      double precision  xx(1),x
c..
c..initialize
      jl=0
      ju=n+1
c..
c..compute a midpoint, and replace either the upper or lower limit
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      go to 10
      endif
c..so this is the table entry
      j=jl
      return
      end
c..
c..
c..
c..
c..
      subroutine hunt(xx,n,x,jlo)
      include 'implno.dek'
c..
c..given an array xx of length n and a value x, this routine returns a value
c..jlo such that x is between xx(jlo) and xx(jlo+1). that is it returns the
c..index of the array whose element is less than and closest to x.
c..the array xx must be monotonic
c..j=0 or j=n indicates that x is out of range
c..on input jlo is a guess for the table entry, and should be used
c..as input for the next hunt. that one assumes the next table entry is
c..relatively close to the old one
c..
c..declare
      logical          ascnd
      integer          n,jlo,jhi,inc,jm
      double precision xx(1),x
c..
c..logical function definition
c..will be true if we are ascending the order of the table; false if descending
      ascnd=xx(n).gt.xx(1)
c..
c..horrid initial guess; go to bisection immediatly
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        go to 3
      endif
c..
c..set the initial hunting increment
      inc=1
c..
c..this is the hunt up section
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
c..
c..if we are the end of the table then we are done hunting
        if(jhi.gt.n)then
          jhi=n+1
c..
c..otherwise replace the lower limit, double the increment and try again
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          go to 1
        endif
c..the value is now bracketed, proceed with the bisection
c..
c..
c..this is the hunt down section
      else
        jhi=jlo
2       jlo=jhi-inc
c..
c..if we are the end of the table then we are done hunting
        if(jlo.lt.1)then
          jlo=0
c..
c..otherwise replace the upper limit, double the increment and try again
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          go to 2
        endif
      endif
c..the value is now bracketed, proceed with the bisection
c..
c..this is the final bisection phase
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      go to 3
      end
c..
c..
c..
c..
c..
      subroutine polcoe(x,y,n,cof)  
      include 'implno.dek'
c..
c..given arrays x and y of length n containing a tabulated function yi = f(xi),
c..this routine returns an array of coefficients cof, also of length n,
c..such that yi = sumj cofj * (xi)**(j-1)
c..this routine is order n**2
c..
c.. 
c..declare
      integer          n,nmax,i,j,k 
      parameter        (nmax=15)
      double precision x(1),y(1),cof(n),s(nmax),phi,ff,b
c..
c..initialze
      do 11 i=1,n   
        s(i)=0. 
        cof(i)=0.   
11    continue  
c..
c..coefficients si of the master polynomial p(x) are found by recurrence
      s(n)=-x(1)
      do 13 i=2,n   
        do 12 j=n+1-i,n-1   
          s(j)=s(j)-x(i)*s(j+1) 
12      continue
        s(n)=s(n)-x(i)  
13    continue  
c..
c..the quantity phi = pi(j.ne.k) (xj - xk) is found as a derivative p(x)
      do 16 j=1,n   
        phi=n   
        do 14 k=n-1,1,-1
          phi=k*s(k+1)+x(j)*phi 
14      continue
        ff=y(j)/phi 
c..
c..coefficients of polynomials in each term of the lagrange formula are
c..found by synthetic division of p(x) by (x-xj). the solution cofk is 
c..accumulated
        b=1.
        do 15 k=n,1,-1  
          cof(k)=cof(k)+b*ff
          b=s(k)+x(j)*b 
15      continue
16    continue  
      return
      end   
c.. 
c.. 
c.. 
c.. 
c.. 
      subroutine polcof(xa,ya,n,cof)
      include 'implno.dek'
c..
c..given arrays x and y of length n containing a tabulated function yi = f(xi),
c..this routine returns an array of coefficients cof, also of length n,
c..such that yi = sumj cofj * (xi)**(j-1)
c..this routine is order n**3 but slightly more stable than the one just given
c..
c.. 
c..declare
      integer          n,nmax,i,j,k 
      parameter        (nmax=15)
      double precision xa(1),ya(1),cof(1),x(nmax),y(nmax),dy,   
     +                 xmin 
c..
c..store a temporary copy
      do 11 j=1,n   
        x(j)=xa(j)  
        y(j)=ya(j)  
11    continue  
c..
c..extrapolate to zero
      do 14 j=1,n   
        call polint(x,y,n+1-j,0.d0,cof(j),dy) 
        xmin=1.e38  
        k=0 
c..
c..and find the remaining xi of smallest absolute value
        do 12 i=1,n+1-j 
          if (abs(x(i)).lt.xmin)then
            xmin=abs(x(i))  
            k=i 
          endif 
c..
c..and keep reducing all the terms
          if(x(i).ne.0.)y(i)=(y(i)-cof(j))/x(i) 
12      continue
c..
c..and eliminate it
        if (k.lt.n+1-j) then
          do 13 i=k+1,n+1-j 
            y(i-1)=y(i) 
            x(i-1)=x(i) 
13        continue  
        endif   
14    continue  
      return
      end   
c.. 
c.. 
c.. 
c.. 
c.. 
      subroutine polin2(x1a,x2a,ya,m,n,x1,x2,y,dy)  
      include 'implno.dek'
c..
c..given arrays x1a (length m) and x2a (length n) of independent variables
c..and an m by n array of funcyion values ya, tabulated at the cartesian
c..grid points defined by x1a and x2a; and values of the desired point
c..x1 and x2; this routine returns an interpolated function value y and
c..an accuracy indication dy (based only on the interpolation in the x1
c..direction).
c..
c..declare
      integer          m,n,nmax,mmax,j,k
      parameter        (nmax=20,mmax=20)
      double precision x1a(1),x2a(1),ya,x1,x2,y,dy, 
     +                 yntmp(nmax),ymtmp(mmax)  
      dimension        ya(m,n)  
c.. 
c..loop over the rows
      do 12 j=1,m   
c..
c..and store the row
        do 11 k=1,n 
          yntmp(k)=ya(j,k)  
11      continue
c..
c..interpolate answer into temporary storage
        call polint(x2a,yntmp,n,x2,ymtmp(j),dy) 
12    continue  
c..
c..do the final interpolation
      call polint(x1a,ymtmp,m,x1,y,dy)  
      return
      end   
c.. 
c.. 
c.. 
c.. 
c.. 
      subroutine bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,  
     +                  ansy,ansy1,ansy2)   
      include 'implno.dek'
c..
c..bicubic interpolation within a square grid. input are y, y1, y2, y12,
c..(as described below in bcucof); x1l and x1u the lower and upper
c..coordinates in the one direction; x2l and x2u likewise in the 2 direction
c..and x1 and x2 the desired point of interpolation. the interpolated
c..function value is returned as ansy with gradients ansy1 and ansy2.
c..
c..this routine calls bcucof below
c..
c.. 
c..declare 
      integer          i        
      double precision y(4),y1(4),y2(4),y12(4),x1l,x1u,x2l,x2u, 
     +                 x1,x2,ansy,ansy1,ansy2,c,d1,d2,t,u   
      dimension        c(4,4)   
c..
      d1=x1u-x1l
      d2=x2u-x2l
      call bcucof(y,y1,y2,y12,d1,d2,c)  
      if(x1u.eq.x1l.or.x2u.eq.x2l) pause 'bad input'    
      t=(x1-x1l)/d1 
      u=(x2-x2l)/d2 
      ansy=0.   
      ansy2=0.  
      ansy1=0.  
      do 11 i=4,1,-1
        ansy=t*ansy+((c(i,4)*u+c(i,3))*u+c(i,2))*u+c(i,1)   
        ansy2=t*ansy2+(3.*c(i,4)*u+2.*c(i,3))*u+c(i,2)  
        ansy1=u*ansy1+(3.*c(4,i)*t+2.*c(3,i))*t+c(2,i)  
11    continue  
      ansy1=ansy1/d1
      ansy2=ansy2/d2
      return
      end   
c.. 
c.. 
c.. 
c.. 
c.. 
      subroutine bcucof(y,y1,y2,y12,d1,d2,c)
      include 'implno.dek'
c.. 
c..given arrays y,y1,y2 and y12 each of length 4, containing the function,
c..gradients and cross derivative at the four corners of a rectangular
c..grid cell (numbered counter-clockwise from the lower left), and given
c..d1 and d2, the lengths of the grid cell in the 1 and 2 directions, this 
c..routine returns the 16 by 16 table c that is used by the routine bcuint 
c..given above for bicubic interpolation
c.. 
c..declare 
      integer          i,j,k,l  
      double precision y(4),y1(4),y2(4),y12(4),d1,d2,c, 
     +                 x(16),wt(16,16),cl(16),d1d2,xx   
      dimension        c(4,4)   
c..
c..god awful table..clean one day 
      data wt/1.,0.,-3.,2.,4*0.,-3.,0.,9.,-6.,2.,0.,-6.,
     +  4.,8*0.,3.,0.,-9.,6.,-2.,0.,6.,-4.,10*0.,9.,-6.,
     +  2*0.,-6.,4.,2*0.,3.,-2.,6*0.,-9.,6.,2*0.,6.,-4.,
     +  4*0.,1.,0.,-3.,2.,-2.,0.,6.,-4.,1.,0.,-3.,2.,8*0.,  
     +  -1.,0.,3.,-2.,1.,0.,-3.,2.,10*0.,-3.,2.,2*0.,3.,
     +  -2.,6*0.,3.,-2.,2*0.,-6.,4.,2*0.,3.,-2.,0.,1.,-2.,  
     +  1.,5*0.,-3.,6.,-3.,0.,2.,-4.,2.,9*0.,3.,-6.,3.,0.,  
     +  -2.,4.,-2.,10*0.,-3.,3.,2*0.,2.,-2.,2*0.,-1.,1.,
     +  6*0.,3.,-3.,2*0.,-2.,2.,5*0.,1.,-2.,1.,0.,-2.,4.,   
     +  -2.,0.,1.,-2.,1.,9*0.,-1.,2.,-1.,0.,1.,-2.,1.,10*0.,
     +  1.,-1.,2*0.,-1.,1.,6*0.,-1.,1.,2*0.,2.,-2.,2*0.,-1.,1./ 
c..
c..initialize and pack a temporary vector x
      d1d2=d1*d2
      do 11 i=1,4   
        x(i)=y(i)   
        x(i+4)=y1(i)*d1 
        x(i+8)=y2(i)*d2 
        x(i+12)=y12(i)*d1d2 
11    continue  
c..
c..matix multiply by the stored table
      do 13 i=1,16  
        xx=0.   
        do 12 k=1,16
          xx=xx+wt(i,k)*x(k)
12      continue
        cl(i)=xx
13    continue  
c..
c..now unpack the result into the output table
      l=0   
      do 15 i=1,4   
        do 14 j=1,4 
          l=l+1 
          c(i,j)=cl(l)  
14      continue
15    continue  
      return
      end   
c.. 
c.. 
c.. 
c.. 
c.. 
      subroutine splie2(x1a,x2a,ya,m,n,y2a) 
      include 'implno.dek'
c..
c..given an m by n tabulated function ya, and tabulated independent variables
c..x1a (length m) and x2a (n values), this routine constructs one dimensional
c..natural cubic splines of the rows of ya and returns the second derivatives
c..in the array y2a
c.. 
c..declare 
      integer          m,n,j,k,nn   
      parameter        (nn=100) 
      double precision x1a(1),x2a(1),ya,y2a,ytmp(nn),y2tmp(nn)  
      dimension        ya(m,n),y2a(m,n) 
c.. 
      do 13 j=1,m   
        do 11 k=1,n 
          ytmp(k)=ya(j,k)   
11      continue
        call spline(x2a,ytmp,n,y2tmp)     
        do 12 k=1,n 
          y2a(j,k)=y2tmp(k) 
12      continue
13    continue  
      return
      end   
c.. 
c.. 
c.. 
c.. 
c.. 
      subroutine splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y) 
      include 'implno.dek'
c..
c..given x1a, x2a, ya, m, and n as described above in the routine splie2 and
c..the y2a produced by that routine and given a desired interpolating point
c..x1 and x2, this routine returns an interpolated function y by bicubic
c..spline interpolation 
c.. 
      integer          m,n,j,k,nn   
      parameter        (nn=100) 
      double precision x1a(m),x2a(1),ya,y2a,x1,x2,y,
     +                 ytmp(nn),y2tmp(nn),yytmp(nn) 
      dimension        ya(m,n),y2a(m,n) 
c.. 
c..perform m evaluations of the row splines constructed by splie2 using
c..the one dimensional spline interpolator splint given above.
      do 12 j=1,m   
        do 11 k=1,n 
          ytmp(k)=ya(j,k)   
          y2tmp(k)=y2a(j,k) 
11      continue
        call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))   
12    continue  
c..
c..get the one dimensional column spline and evaluate it
      call spline(x1a,yytmp,m,y2tmp)
      call splint(x1a,yytmp,y2tmp,m,x1,y)   
      return
      end   
