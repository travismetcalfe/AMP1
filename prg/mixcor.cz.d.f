      subroutine mixcor(x,y,iy,nn,compc,iextrp,nmxfrz)
c
c  sets limit of mixed core, by extrapolating 
c  or interpolating to edge of
c  convection zone from stable region.
c  core mass and radius and extrapolation parameters are
c  set in common/cmxcor/
c
c  extrapolation parameters are defined such that for
c  any function a the extrapolated or interpolated
c  value at the boundary is
c
c      amxc=frmxc*a(nmxcor)+(1-frmxc)*a(nmxcor-1)
c
c  For extrapolation frmxc .gt. 1, and nmxcor is the last point
c  outside the mixed region
c  For interpolation frmxc .le. 1, and nmxcor is the first point
c  within the mixed region
c
c  In this revised version, edge is defined as the first point
c  starting suitably far from the centre such that the extrapolated
c  zero is located before the next mesh point
c
c  For iextrp = 1 test may be based on extrapolation 
c  For iextrp = 2, test requires point with instability,
c  but limit is still set from extrapolation from
c  preceding two points (since unstable point is within
c  region that will subsequently be mixed, and hence modified).
c  For iextrp = 3, test requires point with instability, and
c  parameters are set up for interpolation.
c  For iextrp = 4, find boundary near n = nmxfrz, using either
c  interpolation or extrapolation
c  for iextrp = 10, set boundary from first local maximum which
c  exceeds -ddrmix (in common /cmxcnt/)
c
c  Original version: 14/5/92
c
c  Revised definition of edge 22/5/92
c
c  Revised to force frmxc to be never below -1.
c
c  Modified 19/7/95, to use general heavy-element abundance
c  as set in common/heavy/. This has to be checked with care in
c  connection with convective mixing.
c
c  Modified 5/1/00, to add more general treatment of convective overshoot
c  from convective envelope and core.
c  When overshoot from convective core is included, size of mixing
c  region is set simply by overshoot region.
c
c  Modified 5/6/00, to allow no resetting of core parameters if
c  inmixc = -1.
c
c  Modified 15/8/02, resetting test for core convection based 
c  on size of innermost convective core, to avoid problems
c  with extremely deep convective envelopes at 4He ignition.
c
c  Modified 30/9/02, forcing resetting of qmxcor, etc., when frmxc
c  is forced to -1 after statement no 55.
c
      implicit double precision(a-h, o-z)

c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
c  Note: engenr.n.d.incl replaced by engenr.nnz.d.incl, 6/6/02
      include 'engenr.nnz.d.incl'
c
      dimension x(*), y(iy,*),compc(*)
      dimension qcf(6),qcl(6),yn(10),ddads(10)
      common/clshft/ alshft
      common/caddvr/ addvar(5,nnmax)
      common/eqstd/ xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp,gmm1(4)
      common/heavy/ zatmos, zhc, zh(1)
      common/cnvout/ ddrad,rrcnv,aacnv,ddacad,amach,
     *  drr(5), ddr(5), da(5),pturb,ycnv,dyda, dtxh
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,jnc,frcf(6),frcl(6)
      common/covsec/ icsove, icsovc, alpove, alpovc, 
     *  qmxove, qmxovc, rmxove, rmxovc, imxove, imxovc, 
     *  nmxove, nmxovc
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor
      common/cmxcnt/ ddrmix, inmixc, imixc0, imixc1, imixc2, imixc3,
     *  imixc4, imixc5
      common/cxhder/ xhder(2,nnmax)
      common/cdiffu/ idiffus, itbdif, rhob, rmxdif, tlfdif, dtldif,
     *  ctdfmx, rctdif, rctdf1, idfxbc, ismdif, nsmdif, idiffc1
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3,idgmxc
      common/noiter/ iter, ntime, epscon, eamcon, it_force
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      if(istdpr.gt.0) write(istdpr,*) 
     *  'Entering mixcor with idgmxc =',idgmxc,
     *  ' imixcr =',imixc0, imixc1, ' iextrp =', iextrp
c
c  disallow freezing limit in mixcor
c
      ifrz_on = 0
c
      if(iter.le.10) ifrz_mxc=0
c
      if(inmixc.eq.-1) then
        if(istdpr.gt.0) write(istdpr,105) 
        return
      end if
c
c  test for using core composition in previously computed
c  mixed core
c
      if(nmxcor.gt.0) then
        ncrlim=nmxcor+1
      else 
        ncrlim=nn+1
      end if
c
c  set ncrlim large, to avoid resetting X for the time being
c
      ncrlim=nn+1
      if(istdpr.gt.0) write(istdpr,*) ' In mixcor, ncrlim =',ncrlim
c 
c  initialize to zero
c
      qmxcor=0
      rmxcor=0
      frmxc=0
      nmxcor=0
c
c  if no convection zones have been found previously, stop
c
      if(inc.eq.0) return
c
c  test for fully convective star
c
      if(inc.eq.1.and.nl(1).eq.nn.and.x(nf(1)).ge.-1.e-6) then
	write(istdou,102) 
	if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,102)
	nmxcor=1
	qmxcor=1.d0
	frmxc=1.d0
	rmxcor=1.d0
	return
      end if
c
c  test for using region of convective-core overshoot
c  (note that for now we do not set interpolation fraction
c
      if(imxovc.gt.0) then
c
        qmxcor = qmxovc
c
c  #NB# As a rather desparate attempt to fix overshoot, set 
c  slightly larger mixed core
c
        qmxcor=1.00001d0*qmxcor
        if(istdpr.gt.0) write(istdpr,
     *    '(/'' **** Overshoot. qmxcor slightly increased''/)')
        rmxcor = rmxovc/(1.d11*10.d0**y(1,1))
        nmxcor = nmxovc
        frmxc  = 0.
c
        return
c
      end if
c
c  set existing limits of convection zones
c
c  Note: inc = -1 is used to flag for testing, setting large initial
c  try
c
      if(inc.gt.0) then
c
        do 20 i=1,inc
        n1=nf(i)
        frf=frcf(i)
        frl=frcl(i)
c
        if(n1.eq.1) then
          qlcf=x(n1)
        else
          qlcf=frf*x(n1) + (1-frf)*x(n1-1)
        end if
	n2=nl(i)
	if(n2.eq.nn) then
	  qlcl=x(n2)
	else
          qlcl=frl*x(n2) + (1-frl)*x(n2+1)
        end if
        qcl(i) = 10.d0**qlcl
   20   qcf(i) = 10.d0**qlcf
c
        inc1 = inc
c
      else
c
        inc1 = 1
        qlcf = log10(0.5)
        do 21 n=1,nn
        n1=n
        if(x(n).le.qlcf) go to 22
   21   continue
c
   22   nf(1)=n1
        qcf(1)=10.d0**x(n1)
      end if
c
c  test for convection zones in the core
c  (this may well have to be enhanced later)
c
      ic=0
      if(inc.gt.0) then
	qclim=1.5*qcf(inc)
      else
	qclim=0.9
      end if
      do 25 i=1,inc1
      if(qcf(i).lt.qclim) then
        ic=i
        go to 30
      end if
   25 continue
c
c  test that relevant convection zone was found
c
   30 if(ic.eq.0) return
c
c  now start looking for extrapolated limit, assuming that it is
c  within a fraction of the outermost convective core
c  Note: to avoid problems when the core is very small,
c  start at least at q = 1.e-3
c
      nc=nf(ic)
      is=0
c
      if(iextrp.eq.4) then
c
c  test for unreasonable value at freezing
c
        if(nmxfrz.le.1.or.nmxfrz.ge.nn) then
          if(istdpr.gt.0) write(istdpr,110) nmxfrz
          go to 55
        end if
        n1=nmxfrz-1
        qlim=1.5*10.d0**x(nmxfrz)
      else 
        n1=1
        qlim=min(0.9d0,max(1.d-3,1.2*qcf(ic)))
      end if
c
      if(inc.gt.0.and.qcf(1).gt.0.9999) qlim=min(0.99*qcl(1),qlim)
c
c  set limit of search such that point is in stable region
c
      if(istdpr.gt.0) write(istdpr,*) ' Start search for qlim =',qlim
c
      do 35 n=n1,nn
      qx=10.d0**x(n)
      if(qx.lt.qlim) then
	nlim=n
	go to 37
      end if
   35 continue
c
   37 if(imixc0.eq.1.or.addvar(1,nlim).le.0) then
        drad=fdrad(x(nlim),y(1,nlim),zh(nlim),ak,akr,akt,akx)
      else
        pl=addvar(1,nlim)
        drad=fdradp(x(nlim),y(1,nlim),pl,zh(nlim),fl,ak,akr,akt,akx)
      end if
      qlim=10.d0**x(nlim)
      if(idiffc1.eq.0) then
        ddad=drad-dad(1)
      else
	dtxh=fdtxh(x(nlim),y(1,nlim),y(5,nlim),xhder(2,nlim))
        ddad=drad-dad(1)-dtxh
      end if
      if(ddad.gt.0) then
	if(istdpr.gt.0) write(istdpr,112) nlim, qlim, ddad
	nlim=nlim-10
	go to 37
      end if
c
c  When not freezing boundary, shift back a further few points
c  First test for reasonable value
c
      if(nlim.le.4) then
	write(istdou,'(/a,i4)') 
     *    ' ***** In s/r mixcor, unreasonable nlim =', nlim
	stop 'mixcor'
      end if
c
      if(iextrp.ne.4) nlim=nlim-4
      qlim=10.d0**x(nlim)
c
      if(istdpr.gt.0) write(istdpr,*) ' nlim, qlim =',nlim, qlim
c
      iwrite=2
c
      if(iwrite.eq.2.and.istdpr.gt.0) write(istdpr,'(a)')
     *  'in mixcor, n, q, fl, pl, tl, xh, zh, akl, log L*, drad, ddad'
      do 50 n=nlim,nn
      iwritn=1
      qx=10.d0**x(n)
      call store(y(1,n),yn,6)
      if(n.ge.ncrlim.and.compc(1).gt.0) yn(5)=compc(1)
c
c  test for evaluating gradient at fixed log f or log p
c  Note: test that pl is non-zero to (clumsily) avoid problems
c  when calling for initial model
c
      if(imixc0.eq.1.or.addvar(1,n).le.0) then
        drad=fdrad(x(n),yn,zh(n),ak,akr,akt,akx) 
      else
        pl=addvar(1,n)
        drad=fdradp(x(n),yn,pl,zh(n),fl,ak,akr,akt,akx) 
      end if
      if(idiffc1.eq.0) then
        ddad=drad-dad(1)
      else
	dtxh=fdtxh(x(n),yn,yn(5),xhder(2,n))
        ddad=drad-dad(1)-dtxh
      end if
      if(is.lt.10) then
        is=is+1
      else
        do 45 i=1,is-1
   45   ddads(i)=ddads(i+1)
      end if
      ddads(is)=ddad
      if(iwrite.eq.2.and.istdpr.gt.0) write(istdpr,'(i5,1p12e13.5)') 
     *  n, 10.d0**x(n), yn(2),pl,yn(3),yn(5),zh(n),
     *  log10(ak),yn(4),drad,dad(1),dtxh,ddad
c
c  test for output and type of extra/interpolation
c
      if(is.eq.1) then
        if(iextrp.eq.4.and.idgmxc.ge.2.and.iwritn.eq.1.and.istdpr.gt.0) 
     *    write(istdpr,'(a,i5,1p7e13.5)') 
     *      'in mixcor, n, q, fl, tl, xh, akl, log L*, ddad =',
     *      n, 10.d0**x(n), yn(2),yn(3),yn(5),
     *      log10(ak),yn(4),ddad
        iwritn=0
c
      else
c
        if(ddad.le.0.and.ddads(is-1).lt.ddad) then
          frmxc=ddads(is-1)/(ddads(is-1)-ddad)
          if(frmxc.le.5.and.idgmxc.ge.2) iwrite=1
          if(iwrite.eq.1.and.istdpr.gt.0) 
     *      write(istdpr,'(a,i5,1p7e13.5)') 
     *        'in mixcor, n, q, fl, tl, xh, akl, log L*, ddad =',
     *        n, 10.d0**x(n), yn(2),yn(3),yn(5),
     *        log10(ak),yn(4),ddad
          iwritn=0
        end if
c
c  for iextrp = 4, do only two points
c
        if(iextrp.eq.4) then
          if(idgmxc.ge.2.and.iwritn.eq.1.and.istdpr.gt.0) 
     *      write(istdpr,'(a,i5,1p7e13.5)') 
     *        'in mixcor, n, q, fl, tl, xh, akl, log L*, ddad =',
     *        n, 10.d0**x(n), yn(2),yn(3),yn(5),
     *        log10(ak),yn(4),ddad
          nmxcor=nmxfrz
	  if(istdpr.gt.0) write(istdpr,115) nmxcor,
     *      nmxcor, is, ddads(is-1), ddad, frmxc
          frmxc=ddads(is-1)/(ddads(is-1)-ddad)
          qlmxc=frmxc*x(nmxcor)+(1-frmxc)*x(nmxcor-1)
          qmxcor=10.d0**qlmxc
          rlmxc=frmxc*y(1,nmxcor)+(1-frmxc)*y(1,nmxcor-1)
          rmxcor=(10.d0**rlmxc)/10.d0**y(1,1)
          go to 55
c
c  test for proper extrapolation
c
        else if(iextrp.eq.1.and.
     *    ddad.le.0.and.ddads(is-1).lt.ddad) then
          frmxc=ddads(is-1)/(ddads(is-1)-ddad)
          if((imixc5.eq.0.and.frmxc.le.2).or.
     *       (imixc5.eq.1.and.frmxc.le.1.9)) then
            nmxcor=n
	    if(istdpr.gt.0) write(istdpr,120) nmxcor, 
     *          nmxcor, is, ddads(is-1), ddad, frmxc
            qlmxc=frmxc*x(nmxcor)+(1-frmxc)*x(nmxcor-1)
            qmxcor=10.d0**qlmxc
            rlmxc=frmxc*y(1,nmxcor)+(1-frmxc)*y(1,nmxcor-1)
            rmxcor=(10.d0**rlmxc)/10.d0**y(1,1)
            go to 55
          end if
c
c  test for actual zero crossing.
c  In this case, we may go back and extrapolate from previous
c  two points, if iextrp = 2, for consistency with general procedure
c  For iextrp = 3, interpolate
c
c
        else if(
     *    ddad.gt.0.and.ddads(is-1).lt.0.and.is.gt.2) then
          if(idgmxc.ge.2.and.istdpr.gt.0) 
     *      write(istdpr,'(a,i5,1p7e13.5)') 
     *        'in mixcor, n, q, fl, tl, xh, akl, log L*, ddad =',
     *        n, 10.d0**x(n), yn(2),yn(3),yn(5),
     *        log10(ak),yn(4),ddad
          if(iextrp.eq.2.or.iextrp.eq.10) then
            nmxcor=n-1
	    is1=is-1
            frmxc=ddads(is-2)/(ddads(is-2)-ddads(is-1))
          else
            nmxcor=n
	    is1=is
            frmxc=ddads(is-1)/(ddads(is-1)-ddads(is))
          end if
	  if(istdpr.gt.0) write(istdpr,122) iextrp, nmxcor, 
     *      nmxcor, is, is1, ddads(is1-1), ddads(is1), frmxc
          qlmxc=frmxc*x(nmxcor)+(1-frmxc)*x(nmxcor-1)
          qmxcor=10.d0**qlmxc
          rlmxc=frmxc*y(1,nmxcor)+(1-frmxc)*y(1,nmxcor-1)
          rmxcor=(10.d0**rlmxc)/10.d0**y(1,1)
          go to 55
c
c  for iextrp = 10, test on maximum exceeding -ddrmix
c  extrapolate from previous two points
c
        else if(iextrp.eq.10.and.is.gt.2.and.
     *    ddads(is-1).ge.-ddrmix.and. 
     *    ddad.le.ddads(is-1).and.ddads(is-1).gt.ddads(is-2)) then
          if(idgmxc.ge.2.and.iwritn.eq.1.and.istdpr.gt.0) 
     *      write(istdpr,'(a,i5,1p7e13.5)') 
     *        'in mixcor, n,  q, fl, tl, xh, akl, log L*, ddad =',
     *        n, 10.d0**x(n), yn(2),yn(3),yn(5),
     *        log10(ak),yn(4), ddad
          nmxcor=n-1
	  is1=is-1
          frmxc=ddads(is-2)/(ddads(is-2)-ddads(is-1))
	  if(istdpr.gt.0) write(istdpr,125) iextrp, nmxcor, 
     *      nmxcor, is, is1, ddads(is1-1), ddads(is1), frmxc
          qlmxc=frmxc*x(nmxcor)+(1-frmxc)*x(nmxcor-1)
          qmxcor=10.d0**qlmxc
          rlmxc=frmxc*y(1,nmxcor)+(1-frmxc)*y(1,nmxcor-1)
          rmxcor=(10.d0**rlmxc)/10.d0**y(1,1)
          go to 55
c..        else if(qx.le.0.00001) then
        else if(n.ge.nn-2) then
          if(istdpr.gt.0) then
            write(istdpr,'(a,i5,1p7e13.5)') 
     *        'in mixcor, n, q, fl, tl, xh, akl, log L*, ddad =',
     *        n, 10.d0**x(n), yn(2),yn(3),yn(5),
     *        log10(ak),yn(4), ddad
            write(istdpr,130) qcf(ic),is,n,qx,ddad,ddads(is-1)
	  end if
          nmxcor=0
          frmxc=1
          qmxcor=0
          rmxcor=0
          return
        end if
      end if
   50 continue
c
   55 ireset=0
      if(frmxc.lt.-1.) then
        if(istdpr.gt.0) write(istdpr,140) frmxc
        frmxc=-1.d0
	ireset=1
      else if(iter.gt.10.and.imixc5.eq.1.and.ifrz_on.eq.1.and.
     *  (ifrz_mxc.eq.1.or.abs(frmxc-1.d0).le.0.1d0)) then
        if(istdpr.gt.0) write(istdpr,145) frmxc
	frmxc=1.d0
	ireset=1
	ifrz_mxc=1
      end if
c
      if(ireset.eq.1) then
c
c  reset values at boundary
c
        qlmxc=frmxc*x(nmxcor)+(1-frmxc)*x(nmxcor-1)
        qmxcor=10.d0**qlmxc
        rlmxc=frmxc*y(1,nmxcor)+(1-frmxc)*y(1,nmxcor-1)
        rmxcor=(10.d0**rlmxc)/10.d0**y(1,1)
      end if
      if(istdpr.gt.0) write(istdpr,150) 
     *  qcf(ic), qmxcor, rmxcor, nmxcor, frmxc
      return
  102 format(//' ***** In s/r mixcor, fully convective star is assumed')
  105 format(//' In s/r mixcor inmixc = -1. Core parameters not reset'/)
  110 format(/' ***** Warning in mixcor. nmxfrz =',i5,
     *  '  no boundary set')
  112 format(/
     *  ' Warning in s/r mixcor. Starting point in unstable region'/
     *  ' n, q, ddad =',i5,1p2e13.5)
  115 format(/' With frozen starting point, set nmxcor to',i5/
     *        ' nmxcor, is, ddads(is-1), ddad, frmxc =', 
     *        i5, i3, 1p2e13.5, 0pf11.7/)
  120 format(/' With extrapolation, iextrp = 1, set nmxcor to',i5/
     *        ' nmxcor, is, ddads(is-1), ddad, frmxc =',
     *        i5, i3, 1p2e13.5, 0pf11.7/)
  122 format(/' With interpolation or extrapolation, iextrp =',i2,
     *  ' set nmxcor to',i5/
     *        ' nmxcor, is, is1, ddads(is1-1), ddads(is1) frmxc =',
     *        i5,2i3, 1p2e13.5, 0pf11.7/)
  125 format(/' With proper extrapolation, iextrp =',i2,
     *  ' set nmxcor to',i5/
     *        ' nmxcor, is, is1, ddads(is1-1), ddads(is1) frmxc =',
     *        i5,2i3, 1p2e13.5, 0pf11.7/)
  130 format(//' ***** Error in mixcor. No boundary found.'/
     *         '       qcf(ic), is =',1pe13.5,i3/
     *         '       n, qx, ddad, ddadp =',i5,3e13.5)
  140 format(//' ***** Warning in s/r mixcor: frmxc =',1pe13.5,
     *         ' is unreasonably negative.'/
     *         '       frmxc has been reset to -1'/)
  145 format(//' ***** Warning in s/r mixcor: frmxc =',1pe13.5,
     *         ' is close to 1.'/
     *         '       frmxc has been reset to 1'/)
  150 format(/' in mixcor, q at old boundary =',1pe13.5,
     *  ' new q, r/R =',2e13.5/
     *  ' nmxcor, frmxc =',i5,e13.5)
      end
