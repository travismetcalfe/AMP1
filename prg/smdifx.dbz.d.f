      subroutine smdifx(x,y,nn,iy,icase)
c
c  smooth hydrogen profile and settling velocity for Xh < xhslim
c  in an attempt to get rid of instability with hydrogen settling
c  near hydrogen exhaustion. Also smooth just outside edge of
c  convective core when Xh < 3*xhslim.
c  xhslim is hardcoded to 0.01.
c  Finally restrict Xh to be no smaller than value in convective core.
c
c  For icase = 1: assume output storage of y
c  For icase = 2: assume internal (tnrkt) storage of y
c
c  Original version: 30/10/00
c
      implicit double precision (a-h, o-z)
c
      include 'engenr.nnz.d.incl'
      dimension x(1), y(iy,1)
      dimension w(7),xxh(nnmax),dxh(nnmax)
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec, idcomp, iccomp
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor,nmxcp,qmxcp,rmxcp,
     *  qmxmin,qmxmax
      common/engcnt/ xhzlm1, xhzlm2, nxhzer, icnocs, nspec
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      data w /0.25d0,0.5d0,1.d0,1.d0,1.d0,0.5d0,0.25d0/
c
c  hard-code limit of Xh
c
      xhslim=0.01
      idiag=0
c
      if(istdpr.gt.0) write(istdpr,'(/a,1pe13.5)')
     *  ' Entering smdifx. xhslim =',xhslim
c
      if(icase.eq.1) then
	ixh=5
	idxh=5+iccomp+idcomp
      else
	ixh=4
	idxh=5+idcomp
      end if
c
      qslim=0
      nslim=nn+1
      nmax=min(nn-3,nxhzer)
      if(nmxcor.gt.0) nmax=min(nmxcor,nmax)
      do 20 n=4,nmax
      q=10.d0**x(n)
      if(y(ixh,n).lt.xhslim.or.
     *  (nmxcor.gt.0.and.q.le.1.1*qmxcor.and.y(ixh,n).lt.3*xhslim)) then
	if(qslim.eq.0) then
	  qslim=q
	  nslim=n
	end if 
	sumx=0
	sumdx=0
	sumw=0
	n1=n-4
	kmax=min(7,nxhzer-n+4)
	do 15 k=1,kmax
	n1=n1+1
	sumx=sumx+w(k)*y(ixh,n1)
	sumdx=sumdx+w(k)*y(idxh,n1)
   15   sumw=sumw+w(k)
c
   18   xxh(n)=sumx/sumw
	dxh(n)=sumdx/sumw
      end if
c 
   20 continue
c
      if(qslim.gt.0) then
c
      if(istdpr.gt.0) write(istdpr,
     *  '(/'' In smdifx smooth X between '',i5,'' and '',i5)')
     *  nslim,nmax
c
	if(idiag.eq.1.and.istdpr.gt.0) 
     *    write(istdpr,*) 'Smoothing Xh, Yh; q, orig., smoothed values'
        do 30 n=nslim,nmax
        q=10.d0**x(n)
        if(q.ge.0.9*qslim) then
  	  fact=(qslim-q)/(0.1*qslim)
        else
  	  fact=1.d0
        end if
        fact1=1.d0-fact
	if(idiag.eq.1.and.istdpr.gt.0) write(istdpr,'(1p6e13.5)') 
     *    q, y(ixh,n), xxh(n), y(idxh,n), dxh(n),
     *    fact
        y(ixh,n)=fact1*y(ixh,n)+fact*xxh(n)
   30   y(idxh,n)=fact1*y(idxh,n)+fact*dxh(n)
c
      end if
c
      if(nxhzer.le.nn) then
	do 40 n=nxhzer,nn
	y(ixh,n)=1.d-10
   40   y(idxh,n)=0.d0
      end if
c
      do 50 n=1,nn
      if(nmxcor.gt.0.and.y(ixh,n).lt.0.99*y(ixh,nn)) then
	if(istdpr.gt.0) write(istdpr,120) 
     *    n, 10.**x(n),y(ixh,n),y(ixh,nn)
	y(ixh,n)=y(ixh,nn)
	y(idxh,n)=0.d0
      end if
   50 continue
      return
  120 format(/' At n, q =',i5,1pe13.5,' X reset from',e13.5,
     *  ' to ',e13.5)
      end
