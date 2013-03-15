      subroutine rotker(idsrkr,x,y,aa,el,sig,iy,ia,nn,nprtkr)
c
c  calculates and outputs rotational kernel and beta,
c  to unit idsrkr
c
c  modified 13/8/87 to standardize output
c
c  modified 21/7/94, adding unit number as argument
c
c  ....................................................................
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      dimension x(1),y(iy,1),aa(ia,1)
      common/csumma/ cs(50)
      common/worksp/ w(4,1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      iw=4
c  return in radial case
      if(el.lt.1.e-6) return
c  set integrands
      ell=el*(el+1)
      fct=1.-2./ell
      do 10 n=1,nn
      w(4,n)=x(n)*x(n)*aa(1,n)*aa(5,n)
      w(3,n)=y(2,n)*y(2,n)/ell
      w(1,n)=w(4,n)*(y(1,n)*y(1,n)+w(3,n))
      w(2,n)=y(1,n)-y(2,n)/ell
   10 w(2,n)=w(4,n)*(w(2,n)*w(2,n)+fct*w(3,n))
c  integrate
      do 20 k=1,2
      k2=k+2
   20 call vinta(x,w(k,1),w(k2,1),nn,iw,iw)
c  set beta and kernel
      beta=w(4,nn)/w(3,nn)
      fct=1./w(4,nn)
      do 30 n=1,nn
   30 w(1,n)=fct*w(2,n)
c  output
      if(istdpr.gt.0) write(istdpr,100) beta
      if(nprtkr.gt.1.and.istdpr.gt.0) then
        nd=max0(1,(nn-1)/(nprtkr-1))
        write(istdpr,110) (n,x(n),w(1,n),n=1,nn,nd)
      end if
c  file result
   40 cs(36)=beta
      write(idsrkr) cs,nn,(x(n),w(1,n),n=1,nn)
      return
  100 format(///' results on rotational splitting'//
     *  ' beta =',f12.7)
  110 format(/' n, x, kernel:'//(i5,0pf10.6,1pe13.5))
      end
c
c
