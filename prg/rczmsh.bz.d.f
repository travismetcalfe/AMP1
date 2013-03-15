      subroutine rczmsh(x,y,yp,nn,iy,nptss,icase,init,istrt2)
c
c  reset mesh in vicinity of borders of convection zones
c  to ensure meshpoint at convection-zone boundary. 
c  Uniform mesh is set over 2*nptss+1 points, with the first and
c  last belonging to the original mesh, and the middle point 
c  at the convection-zone boundary.
c
c  If distance from meshpoint is less than epscz1 do not reset
c  If distance from meshpoint is less than epscz2, reset only 
c  point at convection-zone boundary
c  If distance from meshpoint is less than epscz3, reset only
c  3 points around boundary, regardless of npts
c
c  If istrt2 = 1 reset mesh only at bottom of convective envelope,
c  If istrt2 = 2 in addition reset mesh in convective core, if
c  core size has changed substantially
c
c  epscz1, epscz2, epscz3 are hard-coded below.
c
c  For the time being, composition variables and density (or log f) are
c  set by linear extrapolation around the boundary. Otherwise linear
c  interpolation is used.
c
c  Mesh is reset for both y and yp (generally current and previous
c  time-step solution), assumed to be on the same mesh x.
c
c  icase determines form of the solution:
c
c  icase = 1: output form
c  icase = 2: internal form with settling
c
c  Original version: 13/6/02
c
c  Modified 9/8/03 to include resetting at edge of convective core
c  for istrt2 .ge. 2
c
c  Modified 27/8/03, using instead istrt2 = 2 to reset entire mesh
c  in convective core, if change in core size is too big.
c
      implicit double precision(a-h, o-z)
c
      include 'engenr.nnz.d.incl'
c
      parameter(naztmx = nspcmx+3, iyst = 2*(ivarmx+nspcmx))
c
      dimension x(1), y(iy,1), yp(iy,1)
c
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec, idcomp, iccomp
      common/cxhcnc/ compc(nspcmx)
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,inccor,frcf(6),frcl(6)
      common/convpf/ nff(6),nlf(6),incf,inccrf,frcff(6),frclf(6),ifxcon
      common/cntmsh/ wx,epsr,wr,wpt,wxh,wx3,wdgr,
     *  wmshos,drmsos,wmshcc,dqlmsc,qmscmn,qmscmx,
     *  koint,kvint,iprt,icngr,nnt,iwdgrd,istrtc
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor,nmxcp,qmxcp,rmxcp,
     *  qmxmin,qmxmax
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
c  hardcoded limits for resetting meshpoint
c
      epscz1 = 0.01
      epscz2 = 0.05
      epscz3 = 0.10
c
c  for istrt2 = 2, set initial value of previous size of convective core
c
      if(init.eq.1) qmxcpm = qmxcp
c
      idiag=1
c
      if(nptss.le.0) then
	if(idiag.gt.0) write(istdou,100)
	return
      end if
c
      ncz=nl(1)
      frcz=frcl(1)
      xcz=frcz*x(ncz)+(1.d0-frcz)*x(ncz+1)
      if(frcz.lt.0.5) ncz=ncz+1
c
c  test whether point is already sufficiently close to convection-zone
c  boundary
c
      if(abs(frcz).le.epscz1.or.abs(1.d0-frcz).le.epscz1) then
	if(idiag.gt.0.and.istdpr.gt.0) write(istdpr,110) frcz
	ifxcon=1
        nlf(1)=ncz
        frclf(1)=1.d0
        if(istdpr.gt.0) write(istdpr,140) ncz
	go to 40
      else if(abs(frcz).le.epscz2.or.abs(1.d0-frcz).le.epscz2) then
	if(idiag.gt.0.and.istdpr.gt.0) write(istdpr,111) frcz
	npts = 0
      else if(abs(frcz).le.epscz3.or.abs(1.d0-frcz).le.epscz3) then
	if(idiag.gt.0.and.istdpr.gt.0) write(istdpr,112) frcz
	npts = 1
      else
	npts = nptss
      end if
c
c  set new mesh
c
      call rczms1(x,y,yp,nn,iy,ncz,xcz,npts,icase,istrt2,idiag)
c
c  reset nl(1), frl
c
      nl(1)=ncz
      frcl(1)=1.d0
      ifxcon=1
      nlf(1)=ncz
      frclf(1)=1.d0
      if(istdpr.gt.0) write(istdpr,140) ncz
c
   40 if(istrt2.eq.1) return
c
c  test for presence of convective core
c
      if(qmxcor.le.0) return
      if(idiag.ge.1.and.istdpr.gt.0) write(istdpr,142) 
c
c  test for complete reset of core mesh, with large change in position
c  of core edge
c
      if(abs(log10(qmxcor/(qmxcpm+1.d-10))).gt.1.5*dqlmsc) then
        if(istdpr.gt.0) write(istdpr,145) qmxcor,qmxcpm
        ncz=nmxcor
        xcz=log10(qmxcor)
        call rczms2(x,y,yp,nn,iy,ncz,xcz,npts,icase,istrt2,idiag)
        qmxcpm = qmxcor
        iextrp=10
        call mixcor(x,y,iy,nn,compc,iextrp,0)
      end if
c
c  reset mesh at edge of convective core (possibly) for istrt2 .ge. 3
c
      if(istrt2.eq.2) return
c
      ncz=-1
      xcz=log10(qmxcor)
      do n=1,nn
        if(x(n).le.xcz.and.ncz.eq.-1) ncz=n
      end do
      frcz=(x(ncz-1)-xcz)/(x(ncz-1)-x(ncz))
      if(frcz.lt.0.5) ncz=ncz-1
c
c  test whether point is already sufficiently close to convection-zone
c  boundary
c
      if(abs(frcz).le.epscz1.or.abs(1.d0-frcz).le.epscz1) then
	if(idiag.gt.0.and.istdpr.gt.0) write(istdpr,110) frcz
	ifxcon=1
        nff(inc)=ncz
        frcff(inc)=1.d0
        if(istdpr.gt.0) write(istdpr,150) ncz
	return
      else if(abs(frcz).le.epscz2.or.abs(1.d0-frcz).le.epscz2) then
	if(idiag.gt.0.and.istdpr.gt.0) write(istdpr,111) frcz
	npts = 0
      else if(abs(frcz).le.epscz3.or.abs(1.d0-frcz).le.epscz3) then
	if(idiag.gt.0.and.istdpr.gt.0) write(istdpr,112) frcz
	npts = 1
      else
	npts = nptss
      end if
c
c  set new mesh
c
      call rczms1(x,y,yp,nn,iy,ncz,xcz,npts,icase,istrt2,idiag)
c
c  reset nl(1), frl
c
      nmxcor=ncz
      frmxc=1.d0
      ifxcon=1
      nff(inc)=ncz
      frcff(inc)=1.d0
      if(istdpr.gt.0) write(istdpr,150) ncz
c
      return
  100 format(//' ***** Warning, s/r rczmsh called with no points.'/)
  110 format(//' In s/r rczmsh, convection-zone boundary ',
     *  'already at meshpoint; frcz =',f10.8/)
  111 format(//' In s/r rczmsh, convection-zone boundary ',
     *  'very close to meshpoint; frcz =',f10.8/
     *  ' Reset only one point'/)
  112 format(//' In s/r rczmsh, convection-zone boundary ',
     *  'close to meshpoint; frcz =',f10.8/
     *  ' Reset only three points'/)
  140 format(/' In s/r rczmsh ifxcon has been set to 1, nlf(1) to',i5)
  142 format(//' Resetting mesh at edge of convective core')
  145 format(/' New qmxcor = ',1pe13.5,' far from previous point =',
     *  e13.5/' Reset entire core mesh')
  150 format(/' In s/r rczmsh ifxcon has been set to 1,',
     *  ' nmxcor = nff(inc) to',i5)
      end
      subroutine rczms1(x,y,yp,nn,iy,ncz,xcz,npts,icase,istrt2,idiag)
c
c  carry out the actual resetting of mesh for n = ncz - npts to
c  ncz + npts
c
      implicit double precision (a-h, o-z)
c
      include 'engenr.nnz.d.incl'
c
      parameter(naztmx = nspcmx+3, iyst = 2*(ivarmx+nspcmx))
c
      dimension x(1), y(iy,1), yp(iy,1)
      dimension ist1(100),ist2(100), xnew(nnmax)
c
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec, idcomp, iccomp
      common/work/ ynew(iyst,1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      n1=ncz-npts
      n2=ncz+npts
      dx1=(xcz-x(n1))/npts
      dx2=(x(n2)-xcz)/npts
      if(npts.gt.0) then
        do 20 k=1,npts
        xnew(k)=x(n1)+dx1*(k-1)
   20   xnew(npts+k+1)=xcz+dx2*k
      end if
c
      xnew(npts+1)=xcz
c
c  set storage arrays for interpolation, depending on icase
c
      icomps=idcomp+iccomp
      idcom2=idcomp+idcomp
c
      kst1=3
      kst2=1+iccomp+idcom2
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'iccomp,idcomp,idcom2,kst1,kst2 =',
     *  iccomp,idcomp,idcom2,kst1,kst2
      if(icase.eq.1) then
	ist1(1)=1
	ist1(2)=3
	ist1(3)=4
c
	ist2(1)=2
	ist2(2)=5
	do 25 k=1,iccomp
   25   ist2(2+k)=5+k
	if(idcom2.ge.2) then
	  do 30 j=2,idcom2
   30     ist2(1+iccomp+j)=4+iccomp+j
        end if
c
	ix=5
c
      else
	ist1(1)=1
	ist1(2)=3
	ist1(3)=4+idcomp
c
	ist2(1)=2
	do 35 j=1,idcomp
	ist2(1+j)=3+j
   35   ist2(1+idcomp+j)=4+idcomp+j
	do 40 k=1,iccomp
   40   ist2(1+idcom2+k)=4+idcom2+k
c
	ix=4
c
      end if
c
      if(idiag.ge.1.and.istdpr.gt.0) then
	write(istdpr,120) icase,(ist1(k),k=1,kst1)
	write(istdpr,125) (ist2(k),k=1,kst2)
      end if
c
c  set new variables by interpolation or extrapolation
c
      npts21=2*npts+1
      do 50 k=1,kst1
      i=ist1(k)
      do 50 n=1,npts21
   50 call lir1(xnew(n),x,ynew(i,n),y(i,1),1,iy,nn,n,inter)
c
      do 58 k=1,kst2
      i=ist2(k)
      do 55 n=1,npts+1
   55 call lir1(xnew(n),x,ynew(i,n),y(i,1),1,iy,ncz,n,inter)
      if(npts.gt.0) then
        do 57 n=1,npts
        ns=npts+1+n
   57   call lir1(xnew(ns),x(ncz+1),ynew(i,ns),y(i,ncz+1),1,iy,nn-ncz,
     *     n,inter)
      end if
   58 continue
c
      if(idiag.ge.1.and.istdpr.gt.0) then
	write(istdpr,130)
	do 60 n=1,npts21
	ns=n1+n-1
   60   write(istdpr,135) ns,x(ns),xnew(n),y(2,ns),ynew(2,n)
      end if
c
      kst=kst1+kst2
      do 62 n=1,npts21
      ns=n1+n-1
      do 62 k=1,kst
   62 y(k,ns)=ynew(k,n)
c
c  reset at previous time step. Note: here extrapolation is
c  rather meaningless, since in general the convection zone
c  is at a different point; so linear interpolation is used for
c  all variables.
c
      do 65 i=1,kst
      do 65 n=1,npts21
   65 call lir1(xnew(n),x,ynew(i,n),yp(i,1),1,iy,nn,n,inter)
c
      do 72 n=1,npts21
      ns=n1+n-1
      x(ns)=xnew(n)
      do 72 k=1,kst
   72 yp(k,ns)=ynew(k,n)
c
      return
  120 format(//' In s/r rczmsh, icase =',i2/
     *         ' normal interpol.:   ',20i3)
  125 format(  ' interpol/extrapol.: ',20i3)
  130 format(//' Reset quantities in s/r rczmsh.'/
     *  ' n_orig, x_orig, x_new, Xh_orig, Xh_new:'/)
  135 format(i5,1p2e15.7,0p2f12.7)
      end
      subroutine rczms2(x,y,yp,nn,iy,ncz,xcz,npts,icase,istrt2,idiag)
c
c  Reset mesh in entire core, when edge of convective core is
c  changing rapidly.
c
c     original version: 11/8/03.
c
      implicit double precision (a-h, o-z)
c
      include 'engenr.nnz.d.incl'
c
      parameter(naztmx = nspcmx+3, iyst = 2*(ivarmx+nspcmx))
c
      dimension x(1), y(iy,1), yp(iy,1)
      dimension ist1(100),ist2(100), xnew(nnmax), aks(nnmax)
c
      common/rhcn/ a1,a2,a3,a4,zdummy,nvar,ifwrt
      common/caddvr/ addvar(5,nnmax)
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec, idcomp, iccomp
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,inccor,frcf(6),frcl(6)
      common/convpf/ nff(6),nlf(6),incf,inccrf,frcff(6),frclf(6),ifxcon
      common/cntmsh/ wx,epsr,wr,wpt,wxh,wx3,wdgr,
     *  wmshos,drmsos,wmshcc,dqlmsc,qmscmn,qmscmx,
     *  koint,kvint,iprt,icngr,nnt,iwdgrd,istrtc
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor,nmxcp,qmxcp,rmxcp,
     *  qmxmin,qmxmax
      common/work/ ynew(iyst,1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
c  set limit of reset
c
      qlmxcr=log10(qmxcor)
      qlmxcp=log10(qmxcp)
      qlm=max(qlmxcr,qlmxcp)+3*dqlmsc
      nlim=-1
      do n=1,nn
        if(nlim.eq.-1.and.x(n).le.qlm) nlim=n
      end do
      nlim=nlim-1
c
c  set required derivatives
c
      do n=1,nn
        ynew(1,n)=y(4,n)
	ynew(2,n)=y(5,n)
      end do
c
      do i=1,2
        call derive(x,ynew(i,1),ynew(i+2,1),nn,iyst,iyst,1,1)
      end do
c
c  set integrand and carry out trapezoidal integration
c
      tx=wx/(x(nn)-x(1))
      tx=tx*tx
      aks(nlim)=0
      do n=nlim,nn
        tl=ynew(3,n)/(y(4,1)-y(4,nn))
	tl=tl*tl
	txh=wxh*ynew(4,n)
	txh=txh*txh
c
c  add a term at boundary of mixed region,
c  with some suppression in region where hydrogen abundance 
c  is getting very low
c
c  Also reduce term for very small convective core (so far with
c  hard-coded limits)
c
	dxcvc=x(n)-qlmxcr
        dxcvc=min(10.d0,abs(dxcvc)/dqlmsc)
        trmcvc=((qmxcor+1.d-4)/(qmxcor+5.d-3))*wmshcc*exp(-dxcvc*dxcvc)
        if(y(5,n).le.0.1) trmcvc=trmcvc*(0.01+10.*y(5,n))
c
	sump=sum
        sum=sqrt(tx+tl+txh+trmcvc)
	if(n.gt.nlim) aks(n)=aks(n-1)+0.5d0*(sum+sump)*(x(n-1)-x(n))
c
      end do
c
c  set aks and interpolate to set xnew
c
      fct=(nn-nlim)/aks(nn)
      do n=nlim,nn
        aks(n)=nlim+fct*aks(n)
      end do

c
      init=0
      do n=nlim,nn
        aksn=n
	call lir1(aksn,aks(nlim),xnew(n),x(nlim),1,1,nn-nlim+1,init,
     *    inter)
	init=1
      end do
c
c  interpolate to set ynew
c
      init=0
      do n=nlim,nn
        call lir1(xnew(n),x,ynew(1,n),y,nvar,iy,nn,init,inter)
	init=1
      end do
c
      if(idiag.ge.1.and.istdpr.gt.0) write(istdpr,120)
      do n=nlim,nn
        if(idiag.ge.1.and.istdpr.gt.0) write(istdpr,125) 
     *    n, x(n), xnew(n), y(5,n), ynew(5,n)
        do i=1,nvar
	  y(i,n)=ynew(i,n)
        end do
      end do
c
      init=0
      do n=nlim,nn
        call lir1(xnew(n),x,ynew(1,n),yp,nvar,iy,nn,init,inter)
	init=1
      end do
c
      do n=nlim,nn
        do i=1,nvar
	  yp(i,n)=ynew(i,n)
        end do
      end do
c
c  for possible subsequent call of mixcor also reset addvar
c
      init=0
      do n=nlim,nn
        call lir1(xnew(n),x,ynew(1,n),addvar,5,5,nn,init,inter)
	init=1
      end do
c
      do n=nlim,nn
        do i=1,5
	  addvar(i,n)=ynew(i,n)
        end do
	x(n)=xnew(n)
      end do
c
c  reset nmxcor and frmxc
c
      nmxcor=-1
      xcz=log10(qmxcor)
      do n=1,nn
        if(x(n).le.xcz.and.nmxcor.eq.-1) nmxcor=n
      end do
      frmxc=(x(nmxcor-1)-xcz)/(x(nmxcor-1)-x(nmxcor))
c
      return
  120 format(/' Reset core mesh.'//
     *  ' n, log q(old), log q(new), X(old), X(new)'/)
  125 format(i5,1p4e13.5)
      end
