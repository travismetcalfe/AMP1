      subroutine stxhdr(x,y,xhder,nn,iy,idx)
c
c  sets derivative of X (assumed to be in y(1,.)) wrt x, possibly
c  correcting derivative near centre in case of problems
c
c  Original version: 18/11/03
c
      implicit double precision(a-h, o-z)
      dimension x(1), y(iy,1), xhder(idx,1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      call derive(x,y,xhder,nn,iy,idx,1,1)
c
c  smooth derivatives at the last 20 points (clearly needs
c  refinement!) Note: assumption is that derivative varies
c  linearly with q**(1/3) (i.e., with r).
c
      n1=nn-20
      if(istdpr.gt.0) write(istdpr,
     *  '(/'' Resetting X derivative in stxhdr.''//
     *  '' n, x, orig., reset derivative:''/)')
c
      do n=n1,nn
	xhdern=xhder(1,n1)*10.d0**((x(n)-x(n1))/3.d0)
	if(istdpr.gt.0) write(istdpr,'(i5,f10.5,1p2e13.5)') 
     *    n, x(n), xhder(1,n), xhdern
	xhder(1,n)=xhdern
      end do
      return
      end
