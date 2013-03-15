      subroutine conmix(x,y,iy,jc,nn,yc,nv,irscvr)
c
c  mixes material between nmn and nmx, making composition homogeneous.
c  abundances must be in y(4+i,n), i=1,nv. iy is first dimension of y.
c
c  Modified 16/3/89, moving storage of xi to common/sooner/.
c
c  Modified 13/8/90, moving storage of xi to common/xnwvar/.
c
c  Modified 14/8/90, extending mixing of first convection zone 
c  to surface, if log q .gt. qlslim (which is hardcoded below)
c
c  Modified 11/8/91, adding input parameter yc containing central
c  abundances (for a convection zone that extends to the centre)
c  Added parameter irscvr to flag for resetting X and Y in cvr.
c  Also removed resetting of azt for convective core.
c
c  Modified 19/7/95, to use general heavy-element abundance
c  as set in common/heavy/. This has to be checked with care in
c  connection with convective mixing.
c
c  Modified 5/1/00, to add more general treatment of convective overshoot
c  from convective envelope and core.
c  With envelope overshoot, envelope mixing is extended to include
c  overshoot region.
c
c  Modified 24/7/00, resetting aztst(3) and axst(3) for convective core.
c  #AI# This needs to be checked for other abundances.
c
c  Modified 3/1/01, to integrate in all cases with respect to q rather
c  than x = log(q) (this was used previously when the centre was
c  not included in the range).
c
c  Modified 22/8/03, to reset nmxcor, qmxcor, etc. only if convective
c  core has not previously been set (as indicated by nmxcor = 0)
c  Also change to mix in region defined by nmxcor, qmxcor, frmxc
c  when convective core has been set previously by call of s/r mixcor.
c
c  Modified 30/9/03, to test for value of yc and possibly use of
c  y(4+i,nn) if yc(i) is unreasonable.
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      implicit double precision (a-h,o-z)
      include 'engenr.bz.d.incl'
c
      parameter(naztmx = nspcmx + 3)
c
      dimension x(1),y(iy,1),yc(1)
      common/noiter/ iter, ntime
      common/anwvar/ data(8), yi(istrmx,1)
      common/xnwvar/ xi(1)
      common/cntmsh/ wx,epsr,wr,wpt,wxh,wx3,wdgr,
     .  wmshos,drmsos,wmshcc,dqlmsc,qmscmn,qmscmx,
     .  intork,intvnt,iprmsh,icngrm,nnt,iwdgrd,istrtc
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,jnc,frcf(6),frcl(6),
     *  xrcf(6),xrcl(6)
      common/covsec/ icsove, icsovc, alpove, alpovc, 
     *  qmxove, qmxovc, rmxove, rmxovc, imxove, imxovc, 
     *  nmxove, nmxovc
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor,nmxcp,qmxcp,rmxcp,
     *  qmxmin,qmxmax
      common/crxcor/ cqc, cqcp, crxmn(nspcmx), crxmnp(nspcmx)
      common/cmxcnt/ ddrmix, inmixc, imixc0, imixc1, imixc2, imixc3
      common/cdiffu/ idiffus, itbdif, rhob, rmxdif, tlfdif, dtldif,
     *  ctdfmx, rctdif, rctdf1, idfxbc, ismdif, nsmdif, idiffc1, idiffc2
      common/compvr/ cvr(icvrmx,1)
      common/cstcvr/ icvh1, icvhe4, icvhe3, icvc12, icvc13, icvn14,
     *  icvo16
      common/engcnt/ xhzlm1, xhzlm2, nxhzer, icnocs, nspec,
     *  iheccs, nspect
      common/heavy/ zatmos, zhc, zh(1)
      common/ksider/ dm1,dm2,aztst(naztmx),axst(naztmx)
      common/ln10/ amm
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
c  hardcoded limit for shifting lower boundary of convective envelope
c
      epsenv=0.05
c
      qlslim=-1.e-5
c
      nmn=nf(jc)
c
c  test for including envelope overshoot region
c  (for now without interpolation fraction)
c
      if(imxove.eq.jc) then
        nmx=nmxove
	f2=0.
      else
        nmx=nl(jc)
        f2=frcl(jc)
      end if
c
      if(istdpr.gt.0) then
        write(istdpr,*)
	write(istdpr,*) 
     *    'Enter conmix with jc =',jc, '  nmn, nmx =', nmn, nmx
        write(istdpr,*) 'Central values of compositions:',
     *    (y(4+i,nn),i=1,nv)
      end if
c
      if(nmn.ge.nmx) return
      f1=frcf(jc)
      g1=1-f1
      g2=1-f2
c
c  test for extending mixing to surface
c
      if(jc.eq.1.and.x(nmn).gt.qlslim) then
	if(istdpr.gt.0) write(istdpr,110) x(nmn)
	nmn=1
	f1=1
	g1=0
	isrfmx=1
      else
	isrfmx=0
      end if
c
c  test for shifting lower limit for convective envelope
c
      istrt2=mod(istrtc/100,10)
      if(jc.eq.1.and.istrt2.ge.1.and.f2.le.epsenv) then
	write(istdou,115) f2, nmx, x(nmx), nmx+1, x(nmx+1)
	nmx=nmx+1
	f2=1.d0
	g2=0.d0
      end if
c
c  test for using nmxcor and frmxc for mixing convective core
c
      if(nmx.eq.nn.and.nmxcor.gt.0) then
	nmn=nmxcor
	f1=frmxc
	g1=1.d0-frmxc
	if(istdpr.gt.0) write(istdpr,
     *    '(/'' For convective core mixing, use nmxcor, frmxc'')')
      end if
c
      if(istdpr.gt.0) then
        write(istdpr,*) 'Set mixing with'
        write(istdpr,*) 
     *    'nmn, q(nmn), f1, g1 =', nmn, 10.d0**x(nmn), f1, g1
        write(istdpr,*) 
     *    'nmx, q(nmx), f2, g2 =', nmx, 10.d0**x(nmx), f2, g2
      end if
c
c  set independent variable and integrands
c
c  central values. Note that as a fudge (added 30/9/03) test for
c  unreasonable value and replace by innermost meshpoint if needed
c  The use of azt to set yc here needs a careful check, in fact.
c
      if(nmx.eq.nn) then
        xi(1)=0
        do i=1,nv
	  if(yc(i).gt.0) then
            yi(i,1)=yc(i)
	  else
            yi(i,1)=y(4+i,nn)
	  end if
	end do
        ns=nmx
      else
        n1=nmx+1
        xx=f2*x(nmx)+g2*x(n1)
        xi(1)=10.d0**xx
        do 13 i=1,nv
   13   yi(i,1)= f2*y(4+i,nmx)+g2*y(4+i,n1)
        ns=nmx
        if(g2.lt.1.e-5) ns=nmx-1
      end if
c
      ntot=ns-nmn+2
      n1=ns+1
      do 20 n=2,ntot
      n1=n1-1
      xi(n)=10.d0**x(n1)
      do 20 i=1,nv
   20 yi(i,n)=y(4+i,n1)
c
      if(g1.ge.1.e-5) then
c final point interpolation
        ntot=ntot+1
        n1=nmn-1
        xx=f1*x(nmn)+g1*x(n1)
        xi(ntot)=10.d0**xx
        do 21 i=1,nv
   21   yi(i,ntot)= f1*y(4+i,nmn)+g1*y(4+i,n1)
      else
c
        do 24 i=1,nv
   24   yi(i,ntot)= y(4+i,nmn)
c
      end if
c
      dm=xi(ntot)-xi(1)
      if(istdpr.gt.0.and.iheccs.gt.0) write(istdpr,'(/a/(i5,1p4e17.9))') 
     *  ' n, x, X, Y, delta:',
     *  (n,xi(n),yi(1,n),yi(3,n),yi(1,n)+yi(3,n)-0.98d0,n=1,10)
c
c  the integration
c  (Replaced by trapezoidal integration, 16/2/02)
c
      do 30 k=1,nv
      k2=k+nv
      yi(k2,1)=0.d0
      do 30 n=2,ntot
   30 yi(k2,n)=yi(k2,n-1)+0.5d0*(yi(k,n)+yi(k,n-1))*(xi(n)-xi(n-1))
c
      if(istdpr.gt.0) then 
	write(istdpr,'(/'' n, qi, X, integral at limits of range'')')
        write(istdpr,'(i5,1p3e15.7)') 
     *   (n,xi(n),yi(1,n),yi(1+nv,n),n=1,10)
        write(istdpr,'(i5,1p3e15.7)') (n,xi(n),yi(1,n),yi(1+nv,n),
     *    n=ntot-10,ntot)
      end if
c  mean values
      do 40 k=1,nv
      abmn=yi(k+nv,ntot)/dm
      xmax=0
      xmin=1
      do 35 n=1,ntot
      xmax=max(xmax,yi(k,n))
   35 xmin=min(xmin,yi(k,n))
      if(xmax-xmin.gt.1.d-12.and.(abmn.lt.xmin.or.abmn.gt.xmax)) then
c
c  error; integrated average is outside range of abundance values
c
	write(istdou,135) k, abmn, xmin, xmax
	if(istdpr.gt.0.and.istdpr.ne.istdou) 
     *    write(istdpr,135) k, abmn, xmin, xmax
	abmn=(xmin+xmax)/2
      end if
      kp4=k+4
      if(k.eq.1.and.istdpr.gt.0) 
     *  write(istdpr,'(/'' n, q, old X, mixed X:'')')
      do 40 n=nmn,nmx
      if(k.eq.1.and.(n.le.nmn+10.or.n.ge.nmx-10).and.istdpr.gt.0) 
     *  write(istdpr,'(i5,1p3e15.7)') n, 10.d0**x(n), y(kp4,n), abmn
      y(kp4,n)=abmn
c
c#ai#  As a temporary fix, store mixed compositions 
c#ai#  of H and He in cvr
c#ai#  This may have to be cleared up with care.
c
      if(k.eq.1.and.irscvr.eq.1) then
        cvr(1,n)=abmn
	cvr(2,n)=1-abmn-zh(n)
      end if
   40 continue
c
      if(icnocs.ge.1) then
	do 42 n=nmn,nmx
   42   call setcno(y(1,n),n)
      end if
c
      if(istdpr.gt.0) then
	write(istdpr,*) 'Mixed composition:',(y(4+k,nmx),k=1,nv)
        write(istdpr,*) 'Mixed delta:',y(5,nmx)+y(7,nmx)-0.98d0
c
        if(isrfmx.eq.1) write(istdpr,140) y(5,1)
      end if
c
c  reset quantities for central boundary composition
c
      if(nmx.eq.nn) then
        aztst(3)=y(5,nn)
	axst(3)=0.d0
        if(ispxx3.eq.1) then
          ishft = 3
        else
          ishft = 2
        end if
        do 70 i=2,nv
        aztst(ishft+i)=y(4+i,nn)
   70   axst(ishft+i) =0.d0
	if(istdpr.gt.0) write(istdpr,155) aztst(3), axst(3)
      end if
c
c  test for resetting qmxcor and cqc
c  (only if nmxcor has not previously been set; modification 22/8/03)
c
      if(nmxcor.le.0.and.nmx.eq.nn) then
	qmxcor=xi(ntot)
	cqc=qmxcor
	nmxcor=nmn
	frmxc=f1
	if(istdpr.gt.0) write(istdpr,150) qmxcor, nmxcor, frmxc
      end if
      return
  110 format(/' log q =',1pe13.5,
     *  ' at top of first convection zone in s/r conmix'/
     *  ' mixing extended to surface')
  115 format(//' In s/r conmix, f2 =',1pe13.5,' is unreasonably small'/
     *  ' Shift lower limit of convective envelope,'/
     *  ' from n =',i5,
     *  ' log(q) =',1pe13.5,' to n =',i5,' log(q) =',e13.5)
  130 format(//' integrand in conmix'//(i4,1p2e13.5))
  135 format(/' ***** Error in s/r conmix for k =',i2,'.'/
     *        '       Integrated average =',1pe15.7,
     *        ' outside range ',2e15.7/
     *        '       Reset to average of range'/)
  140 format(/' Surface X has been reset to',f12.7)
  150 format(/' In s/r conmix, qmxcor and cqc reset to',1pe13.5/
     *        '                nmxcor, frmxc  reset to',i5,0pf10.5)
  155 format(/' In s/r conmix, aztst(3) and axst(3) reset to',1p2e13.5)
      end
