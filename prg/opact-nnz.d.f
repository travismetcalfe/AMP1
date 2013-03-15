      subroutine opact(rl,tl,xh,zh,ak,akr,akt,akx)
c
c  opacity driving routine
c  ***********************
c
c  Modified to incorporate call of WD/GH routines (for iwdopc .ge. 8)
c
c  controls are in common/opctcl/
c
c  when iwdopc .ne. 1 calls spline interpolation routines. In thic case
c  the value of iwdopc determines the type of interpolation in X:
c 
c  iwdopc = 0: use old scheme with parabolic interpolation.
c    Note: the cluster of points is chosen to be as close as possible
c    to target x. This leads to switch of cluster in the middle of
c    mesh interval, and hence to discontinuous interpolating function.
c    A little unfortunate!
c  iwdopc = 2: use linear interpolation
c  iwdopc = 3: use 4-point Lagrangian interpolation.
c
c  when iwdopc   =  1 calls w. dappen opacity routines with 4 point
c                     lagrangian interpolation in log rho and log t,
c                     and interpolation in x and z. value of z is
c                     taken from zh in common/heavy/ (for consistency
c                     with evolution programme).
c
c  iwdopc .ge. 8: Use WD and GH routines to interpolate in OPAL
c    and Kurucz tables
c
c  when alamop .gt. 0 and iopacm .gt. 0 match separate atmospheric
c  opacity, multiplied by alamop, to interior opacity.
c  when iwdopc = 1 or iopacm = 2 atmospheric
c  opacity is evaluated at heavy element abundance zatmop, from
c  common/opctcl/, otherwise at standard composition.
c
c  for iopacm = 1, iwdopc .ne. 1 makes smooth transition between
c  alamop*kappa below log t = tlmopc and kappa above log t = tlmopc.
c
c  for iopacm = 1 and iwdopc = 1 makes smooth transition between
c  alamop*kappa(z=zatmop) below log t = tlmopc and kappa(z=zh)
c  above log t = tlmopc.
c
c  when iopacm = 2 matches atmospheric opacity, based on auer et al
c  routines called through s/r opacsf to either spline
c  opacity (for iwdopc .ne. 1) or w.d. opacity (for iwdopc = 1),
c  at log t = tlmopc
c
c  when rl .lt. rhlimn when using opacity matching, the pure atmospheric
c  opacity is used. rhlimn is initialized to -9. in data statement.
c
c  when ifdgop = 1 the final opacity is multiplied by the fudge factor
c  10.**fdgopl.
c
c  when ifdgop = 2 the final opacity is multiplied by the fudge function
c  10**(fdgopl*exp(-((tl-tlopfg)/dtlopf)**2)) and the log t derivative
c  is modified correspondingly. this gives a modification of the opacity
c  localized around log t = tlopfg. tlopfg and dtlopf are transmitted
c  in common/opctfg/
c
c  when ifdgop = 3 the final opacity is multiplied by the fudge function
c  10**(fdgopl*exp(-(dtl/dtlopf)**2)) , where dtl = tl - tlopf1 for
c  tl .le. tlopf1, dtl = 0 for tlopf1 .le. tl .le. tlopfg, and
c  dtl = tl - tlopfg for tl .ge. tlopfg; the log t derivative is
c  modified correspondingly. this gives a modification of the opacity
c  between log t = tlopf1 and log t = tlopfg. tlopfg, dtlopf and
c  tlopf1 are transmitted in common/opctfg/
c
c  when ifdgop = 4, read in table of opacity corrections as
c  a function of log T, from unit idopfd, and interpolate in
c  that to obtain opacity modification
c
c  the parameters in common/opctcl/ are initialized in block data
c  program blopac, to the following values:
c
c     data alamop,zatmop,tlmopc,dtmopc,fdgopl,iwdopc,iopacm,ifdgop
c    *    /   0.,  0.02,   5.,    0.2,    0.,   0,     0,     0    /
c
c  the parameters in common/opctfg/ are initialized in block data
c  program blopac, to the following values:
c
c     data tlopfg,dtlopf,tlopf1
c    *    /   0.,  0.15,   0. /
c
c                                -------------------------
c
c  on entry rl = log10(rho), tl = log10(t), xh = x.
c
c  returns  ak = log10(kappa), akr = (dlog kappa/dlog rho)t,x,
c  akt = (dlog kappa/dlog t)rho,x and akx = (dlog kappa/dlog x)t,rho
c  in argument list.
c
c  also returns, in common /opcxdr/, akxa = (dlog kappa/d x) and,
c  when iwdopc = 1, 8 or 9, akz = (dlog kappa/dlog z)
c
c  to initialize tables, s/r opcscf(iwdopc,iopacm) must be called
c  before call of s/r opact. this also initializes iwdopc and
c  iopacm in common/opctcl/. note that for iwdopc .ne. 1, or
c  iopacm = 2, information must be provided in common/opccnt/
c
c  entry opacat(rl,tl,xh,ak,akr,akt,akx) computes pure atmospheric
c  opacity, including factor alamop. opacity computed
c  is determined by value of iopacm, as above.
c
c                         *************************
c
c  modified 15/9/86: save statement added in all routines
c     (needed in f77, apparently)
c
c  modified 15/9/86: make read of tables entirely unformatted
c     (previously composition was read formatted)
c
c  modified 9/8/90, to fix error in extrapolation for inner tables
c
c  modified 14/8/90, to fix up setting of ranges in s/r opcres
c  if desired range is larger than table range.
c
c  modified 30/10/90, making a quick fix in problems in s/r opcext
c  (when trying to extrapolate using point outside tables).
c  **** This badly needs fixing. ****
c
c  modified 10/7/92, to incorporate WD/GH routines
c
c  Modified 19/7/95, including heavy-element abundance as 
c  argument zh in argument list.
c
c  Modified 11/9/96, including option of ifdgop = 4
c
c  Modified 15/9/97, allowing extended iwdopc for use with GH
c  package version 9 (and later)
c
c  Double precision version.
c  ++++++++++++++++++++++++
c 
c  Note: this version still assumes that *opacity tables* are given
c     on single-precision binary form.
c
c  Dated 11/3/90
c
      implicit double precision (a-h,o-z)
      common/opcxdr/ akxa,akz
      common/ln10/ amm
      common/opctcl/ alamop,zatmop,tlmopc,dtmopc,fdgopl,iwdopc,
     *  iopacm,ifdgop
      common/opctfg/ tlopfg,dtlopf,tlopf1,idopfd
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      external blopac, blstio
c
      save
c
      data rhlimn /-9./
c
      go to 5
c
c  pure atmospheric opacity
c
c..      entry opacat(rl,tl,xh,ak,akr,akt,akx)
c..      go to 35
c
c  ** extra diagnostics
c
    5 if(idgopc.eq.1.and.istdpr.gt.0)
     *  write(istdpr,9991) rl,tl,xh,iwdopc,iopacm,ifdgop
 9991 format(/' in s/r opact: rl, tl, xh, iwdopc, iopacm, ifdgop =',
     *  1p3e15.7,3i10)
c
      if((xh.gt.1.or.xh.lt.0).and.istdpr.gt.0) write(istdpr,110) 
     *  xh, rl, tl
c
c  for iwdopc call WD/GH routines, and bypass, for the moment
c  the section with atmospheric matching (this is an awful
c  mess, in any case)
c
      iwdop0=mod(iwdopc,10)
c
      if(iwdop0.eq.8.or.iwdop0.eq.9) then
	xh1=max(xh,1.d-10)
	call opalkz(rl, tl, xh1, zh, ak, akr, akt, akx, ierr)
	akxa = akx/(amm*xh)
	go to 50
      end if
c
c  fudge if iwdopc = -1
c
      if(iwdopc.eq.-1) then
        ak=4-3*(tl-5)*(tl-5)
        akr=0
        akt=6*(5-tl)
        akx=0
        return
      end if
c
c  test for matching
c
      imatch=0
      if(alamop.gt.0.and.iopacm.gt.0) go to 30
c
c  test for spline interpolation
c
   10 if(iwdopc.eq.1) go to 20
c
   15 call opacsi(rl,tl,xh,ak,akr,akt,akx)
      if(imatch.eq.1) go to 40
      go to 50
c
c  call w.d. routine
c
   20 yh=1-xh-zh
      call kappwd(rl,tl,xh,yh,ak,akr,akt,akxa,akz)
      if(imatch) 50,50,40
c
c  matching.
c
c  test for pure atmospheric opacity (when tl .lt. tlmopc - dtmopc or
c  rl .lt. rhlimn)
c
   30 tlmin=tlmopc-dtmopc
      if(tl.le.tlmin.or.rl.lt.rhlimn) go to 35
c  test for pure interior opacity and, if so, keep imatch = 0
      if(tl.lt.tlmopc+dtmopc) imatch=1
c  set interior opacity
      if(iwdopc-1) 15,20,15
c  pure atmosphere
   35 imatch=2
c
c  find atmospheric value
c
   40 if(iopacm.ne.2) go to 45
c  auer et al opacity
c
c  store value of akxa (note that akxa in common/opcxdr/ is zeroed
c  in s/r opacsf)
      akxa0=akxa
c
      call opacsf(rl,tl,xh,ak1,akr1,akt1,akx1)
c  reset akxa
      akxa=akxa0
      akxa1=0
      go to 47
c
c  test for w.d. opacity in atmosphere
c
   45 if(iwdopc.eq.1) go to 46
c  spline opacity. if imatch .ne. 2 opacity has already been
c  calculated.
      if(imatch.eq.2) call opacsi(rl,tl,xh,ak,akr,akt,akx)
      ak1=ak
      akr1=akr
      akt1=akt
      akxa1=akxa
      go to 47
c  w.d. opacity
   46 yh=1-xh-zatmop
      call kappwd(rl,tl,xh,yh,ak1,akr1,akt1,akxa1,akz1)
c  correct with alamop
   47 ak1=ak1+log10(alamop)
      akz1=0
c
c  test for pure atmospheric opacity
c
      if(imatch.eq.1) go to 48
      ak=ak1
      akr=akr1
      akt=akt1
      akxa=akxa1
      go to 50
c  matching function
   48 phi=phism(tl,tlmin,2*dtmopc,dphi)
c  set matched values
      phi1=1-phi
      ak=phi*ak+phi1*ak1
      akr=phi*akr+phi1*akr1
      akt=phi*akt+phi1*akt1+dphi*(ak-ak1)
      akxa=phi*akxa+phi1*akxa1
      akz=phi*akz
c
c  set akx
c
   50 if(iwdopc.eq.1) akx=xh*amm*akxa
c
c  test for inclusion of fudge factor
c
      if(ifdgop.le.0.or.(ifdgop.lt.4.and.fdgopl.eq.0)) then
	return
c
      else if(ifdgop.eq.1.or.dtlopf.eq.0) then
c
c  constant change
c
        ak=ak+fdgopl
        return
c
      else if(ifdgop.eq.2) then
c
c  include gaussian fudge function and derivative
c
        xtl=(tl-tlopfg)/dtlopf
        extl=fdgopl*exp(-xtl*xtl)
        ak=ak+extl
        akt=akt-2*extl*xtl/dtlopf
        return
c
      else if(ifdgop.eq.3) then
c
c  include pseudo-gaussian fudge function and derivative between
c  tlopf1 and tlopfg
c
        xtl=min(tl-tlopf1,0.d0,tlopfg-tl)/dtlopf
        extl=fdgopl*exp(-xtl*xtl)
        ak=ak+extl
        akt=akt-2*extl*xtl/dtlopf
        return
c
      else if(ifdgop.eq.4) then
c
c  include opacity correction interpolated from table on
c  d/s idopfd
c
        call opfdtb(ifdgop, idopfd, tl, extl, dextl)
	ak=ak+extl
	akt=akt+dextl
c
      else
	write(istdou,120) ifdgop
	if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,120) ifdgop
	stop 'Stop in opact, wrong ifdgop'
      end if
c
  110 format(' ****** Warning in s/r opact. xh = ',1pe13.5,
     *  ' at log rho, log T =',0p2f10.5)
  120 format(//' ***** Error in s/r opact. ifdgop = ',i5,
     *  ' not implemented')
      end
      subroutine opcscf(iwdopc,iopacm)
c
c  opacity initialization routine
c  ******************************
c
c  when iwdopc .ne. 1 or .lt. 8 use spline interpolation
c
c  when iwdopc   =  1 use old w. dappen opacity routines
c  (not currently implemented)
c
c  when iwdopc .ge. 8, use WD/GH routines for OPAL and Kurucz
c  opacities
c
c  when iopacm .gt. 0 prepare for matching to atmospheric opacity.
c  for iopacm = 1 use w.d. opacity at z = zatmop.
c  for iopacm = 2 match to auer et al opacities. in this case
c  atmospheric x must be set into xhsopc in common/opccnt/,
c  for initialization of auer et al routines.
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
      implicit double precision (a-h,o-z)
      common/opctcl/ alamop,zatmop,tlmopc,dtmopc,fdgopl,iwdop1,
     *  iopcm1,ifdgop
      common/opccnt/ xhsopc,tsmn,tstr,rhsmn,rhsmx,timx,rhimn,rhimx,
     *  sigstr,inopc,idmopc
      common/xyz/ xop,yop,zop
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      data initwd,initsf /0,0/
c
      if(istdpr.gt.0) write(istdpr,'(//'' Setting up opacity''/)')
c
   10 if(iwdopc.eq.-1) then
c
c  fudge for iwdopc = -1
c
        return
c
      else if(iwdopc.eq.1) then
c
c  initialize old w. d. tables
c
        iwdop1=1
c
c  test whether w.d. tables have already been initialized
c
        if(initwd.ne.1) then
c
          call opinwd
          if(istdpr.gt.0) write(istdpr,110)
c
          initwd=1
        end if
c
      else if(iwdopc.ge.8) then
c
c  initialize new WD/GH tables
c
	iwdop1=iwdopc
	call opingh(inopc,iwdopc)
        if(istdpr.gt.0) write(istdpr,115)
c
      else
c
c  initialize spline tables
c
        iwdop1=iwdopc
        call opinsi
c..      idgopc=idgopp
        if(istdpr.gt.0) write(istdpr,120)
c
      end if
c
c  test for initialization of auer et al routines
c
   30 iopcm1=iopacm
c
      if(iopacm.gt.0.and.istdpr.gt.0) write(istdpr,100) 
     *  iopacm,zatmop,tlmopc,dtmopc
      if(iopacm.ne.2.or.initsf.eq.1) return
      xop=xhsopc
      zop=zatmop
      yop=1-xop-zop
      call iging
      initsf=1
c
      return
  100 format(/' opacity matching type',i2//' zatm =',f10.7,
     *  '  tlmopc, dtmopc =',2f10.5)
  110 format(//' old w. dappen opacity')
  115 format(//' WD/GH OPAL and Kurucz opacities')
  120 format(//' spline interpolation opacity')
      end
      block data blopac
c
c  initialize controlling parameters for opact, in commons /opctcl/ and
c  /opccnt/. also initialize parameters in common/opcmtc/ for
c  s/r opacsf, and set bhh to 1 in common/bfacto/ corresponding to
c  no departure coefficients in s.f. opacity routines.
c
      implicit double precision (a-h,o-z)
      common/opctcl/ alamop,zatmop,tlmopc,dtmopc,fdgopl,iwdopc,
     *  iopacm,ifdgop
      common/opctfg/ tlopfg,dtlopf,tlopf1,idopfd
      common/opccnt/ xhsopc,tsmn,tstr,rhsmn,rhsmx,timx,rhimn,rhimx,
     *  sigstr,inopc,idmopc
      common/opcmtc/ tlm1,dtlm1,fropc,imon
      common/bfacto/ bhh(6)
      data alamop,zatmop,tlmopc,dtmopc,fdgopl,iwdopc,iopacm,ifdgop
     *    /   0.,  0.02,   5.,    0.2,    0.,   0,     0,     0    /
      data tlopfg,dtlopf,tlopf1
     *    /   0.,  0.15,   0. /
      data xhsopc,tsmn,tstr,rhsmn,rhsmx,timx,rhimn,rhimx,sigstr,inopc
     *    /  0.7,  3.5, 6.,  -10.,  0.,  7.3,  -3.,  3.,   -5.,   13/
      data fropc,imon /6.d14,1/
      data bhh /6*1./
      end
c..      double precision function phism(t,t0,dt,dphi)
      function phism(t,t0,dt,dphi)
c
c  calculates phism = phi(t), such that phi = 0 for t .le. t0,
c  phi = 1 for t .ge. t0 + dt, and phi is continuous with
c  continuous first derivative. dphi is derivative of phi with
c  respect to t.
c
      implicit double precision (a-h,o-z)
c
      x=2*(t-t0)/dt-1
      if(abs(x).ge.1) go to 20
      phism=0.75*x*(1-x*x/3)+0.5
      dphi=0.75*(1-x*x)
      return
c
   20 dphi=0
      if(x.le.-1) go to 25
      phism=1
      return
c
   25 phism=0
      return
      end
      subroutine opacsf(rhl,tl,xh,ak,akr,akt,akx)
c  finds opacity at log rho = rhl, log t = tl,
c  from auer et al routines, using s. frandsen interfaces.
c
c  heavy element abundance must be set in zatm in common/opctcl/, for
c  consistency with evolution programme.
c
c  controls in common/opcmtc/ (tlm, dtlm not used here)
c
c  imon = 1: calculate monochromatic opacity at frequency fropc
c  imon .ne. 1: calculate rosseland mean opacity
c
c  defaults (set in block data blopac): imon = 1, fropc = 6.d14
c  (corresponding to monochromatic opacity at 5000 a)
c
c  calls s/r eqstf to iterate on fl to get the given rhl, 
c  as composition, especially z, and hence fl may be different
c  from composition elsewhere in the calculation.
c  assumes that current value of fl is stored in common/eqstfl/,
c  e.g. by previous call of eqstf.
c
c  to save computing time, derivatives are only calculated when
c  change in rhl or tl since last evaluation of the derivatives
c  is greater than epsder (initialized to 0.005 in data statement).
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
      implicit double precision (a-h,o-z)
      logical nosd,notd
      dimension aksfr(5),tld(5),fld(5),eqd(90),eqds(90),xhvs(30)
      common/opcmtc/ tlm,dtlm,fropc,imon
      common/opcxdr/ akxa,akz
      common/opctcl/ alamop,zatmop,tlmopc,dtmopc,fdgopl,iwdopc,
     *  iopacm,ifdgop
      common/eqstd/ xii(4),ane(10),rhh(20),dmm(56)
      common/eqstfl/ flcm
      common/hviond/ xhvmn(10),xhv(125)
      common/xyz/ x,y,z
      common/cscat/ al(64)
      common/cmax/ nmax,iflsfr
      common/eqdpco/ eqdp(3),idpco
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      equivalence(eqd(1),xii(1))
c
      external blopac
c
      save
c
      data rhlp,tlp /-1.d4,-1.d4/
      data rhlkdp,tlkdp /-1.d4,-1.d4/
      data epsder /0.005/
c
c  set composition in common/xyz/
c
      x=xh
      z=zatmop
      y=1-x-z
c  store eqstf data
      call store(eqd,eqds,90)
      call store(xhv,xhvs,30)
c  set equation of state at opacity composition at fixed rhl
      nosd=.true.
      notd=.true.
      nit=0
      fl=flcm
      drhla=abs(rhl-rhlp)
      dtla=abs(tl-tlp)
      if(drhla.lt.0.01.and.dtla.lt.0.01) fl=flp
c  iteration for fl
   15 call eqstf(fl,tl,x,y,z,nosd,notd)
      dfl=(rhl-log10(rhh(1)))/rhh(2)
      if(abs(dfl).lt.1.e-8) go to 18
      if(nit.gt.10) go to 90
      fl=fl+dfl
      go to 15
c
   18 flp=fl
      rhlp=rhl
      tlp=tl
c  call sfr routines
   20 t=10.d0**tl
      aner=rhh(1)*ane(1)
c
      call sfropc(t,aner,opc)
      aksfr(1)=log10(opc)
c
      ak=aksfr(1)
c
c  test for resetting opacity derivatives
c
      if(abs(rhl-rhlkdp).gt.epsder.or.abs(tl-tlkdp).gt.epsder)
     *  go to 20100
c
c  use previous values
c
      akt=aktp
      akr=akrp
      go to 40
c
c  prepare for setting derivatives
c  *******************************
c
20100 dtfl=-0.0001
      do 21 i=2,3
      fld(i)=fl
      tld(i)=tl+dtfl
   21 dtfl=-dtfl
      do 22 i=4,5
      tld(i)=tl
      fld(i)=fl+dtfl
   22 dtfl=-dtfl
c  set flag to keep nmax fixed
      iflsfp=iflsfr
      iflsfr=1
c  fix idpco
      idpcop=idpco
      idpco=0
      if(idpcop.gt.0) idpco=3
c
      do 25 i=2,5
      call eqstf(fld(i),tld(i),x,y,z,nosd,notd)
      t=10.d0**tld(i)
      aner=rhh(1)*ane(1)
      iflcio=0
      if(tl.le.4.3.and.tl.ge.3.85) iflcio=1
      call sfropc(t,aner,opc)
      iflcio=0
   25 aksfr(i)=log10(opc)
c  set derivatives
      akstf=(aksfr(3)-aksfr(2))/(tld(3)-tld(2))
      aksft=(aksfr(5)-aksfr(4))/(fld(5)-fld(4))
c  reset eqd and xhv
      call store(eqds,eqd,90)
      call store(xhvs,xhv,30)
c  reset flag
      iflsfr=iflsfp
c  reset idpco
      idpco=idpcop
c
c  transform derivatives
c
      akt=akstf-aksft*rhh(3)/rhh(2)
      akr=aksft/rhh(2)
c
      aktp=akt
      akrp=akr
      rhlkdp=rhl
      tlkdp=tl
c
   40 akx=0
      akxa=0
      akz=0
c
      return
c
c  diagnostics for convergence failure of rho iteration
c
   90 write(istdou,100) nit,rhl,tl,fl,dfl
      if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,100) 
     *  nit,rhl,tl,fl,dfl
      stop 'Stop in opacsf'
  100 format(//1x,10(1h*),' rho iteration failed to converge in',
     *  ' s/r opacsf after ',i3,' iterations.'/
     *  11x,' rhl,tl =',2f10.5,'  final fl, dfl =',f10.5,1pe13.5)
      end
      subroutine opacm(rhl,tl,fl,xh,ak,akr,akt,ider)
c  finds opacity by matching sfr routines to table interpolation. note:
c  s/r eqstf must be called  before call of opacm.
      implicit double precision (a-h,o-z)
      dimension aksfr(5),tld(5),fld(5),eqd(90),eqds(90),xhvs(30)
      dimension rhld(5)
      common/opcmtc/ tlm,dtlm,fropc,imon
      common/eqstd/ xii(4),ane(10),rhh(20),dmm(56)
      common/hviond/ xhvmn(10),xhv(125)
      common/xyz/ x,y,z
      common/cscat/ al(64)
      common/ciog/ iflcio
      common/cmax/ nmax,iflsfr
      common/eqdpco/ eqdp(3),idpco
      common/bfacto/ bhh(6)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      equivalence(eqd(1),xii(1))
c
      save
c
      iwrdtl=0
c
      akct=-60
      aksfr(1)=-60
      tlmin=tlm-dtlm
      tlmax=tlm+dtlm
c  test for call of opact
      if(tl.lt.tlmin) go to 20
      call opact(rhl,tl,xh,z,ak,akr,akt,akx)
      akct=ak
      al(1)=0
      if(tl.gt.tlmax) go to 50
c  call sfr routines
   20 t=10.d0**tl
      aner=rhh(1)*ane(1)
c
      call sfropc(t,aner,opc)
      aksfr(1)=log10(opc)
      rhld(1)=log10(rhh(1))
c  test for derivatives
      if(ider.ne.1) go to 30
c  store eqstf data
      call store(eqd,eqds,90)
      call store(xhv,xhvs,30)
      dtfl=-0.0001
      do 21 i=2,3
      fld(i)=fl
      tld(i)=tl+dtfl
   21 dtfl=-dtfl
      do 22 i=4,5
      tld(i)=tl
      fld(i)=fl+dtfl
   22 dtfl=-dtfl
c  set flag to keep nmax fixed
      iflsfp=iflsfr
      iflsfr=1
c
      do 25 i=2,5
      call eqstf(fld(i),tld(i),x,y,z,.true.,.true.)
      rhld(i)=log10(rhh(1))
      t=10.d0**tld(i)
      aner=rhh(1)*ane(1)
      iflcio=0
      if(tl.le.4.3.and.tl.ge.3.85) iflcio=1
      call sfropc(t,aner,opc)
      iflcio=0
   25 aksfr(i)=log10(opc)
c  set derivatives
      akstf=(aksfr(3)-aksfr(2))/(tld(3)-tld(2))
      aksft=(aksfr(5)-aksfr(4))/(fld(5)-fld(4))
c  reset eqd and xhv
      call store(eqds,eqd,90)
      call store(xhvs,xhv,30)
c  reset flag
      iflsfr=iflsfp
c
      iwrder=0
      if(iwrder.ne.1) go to 28
      rht=      (rhld(3)-rhld(2))/(tld(3)-tld(2))
      rhf=(rhld(5)-rhld(4))/(fld(5)-fld(4))
      if(istdpr.gt.0) write(istdpr,301) 
     *  tld,fld,rhld,rhf,rht,rhh(2),rhh(3),aksft,akstf
  301 format(/3(5f15.7/),6f15.7)
c  transform derivatives
   28 akstr=akstf-aksft*rhh(3)/rhh(2)
      aksrt=aksft/rhh(2)
c  matching
   30 if(tl.ge.tlmin) go to 40
      ak=aksfr(1)
      if(ider.ne.1) go to 50
      akr=aksrt
      akt=akstr
      go to 50
c
   40 phs=phism(tl,tlmin,2*dtlm,dphs)
      phsi=1-phs
      aktbl=ak
      ak=phs*ak+phsi*aksfr(1)
      al(1)=phsi*al(1)
      if(ider.ne.1) return
      akr=phs*akr+phsi*aksrt
      akt=phs*akt+phsi*akstr+dphs*(aktbl-aksfr(1))
c  output
   50 if(iwrdtl.ne.1) return
      if(tl.gt.4.8.or.tl.lt.3.4) return
      opcct=10.d0**akct
      opcsf=10.d0**aksfr(1)
      opc=10.d0**ak
      t=10.d0**tl
      if(istdpr.gt.0) write(istdpr,100) 
     *  tlm,dtlm,tlmin,tlmax,fropc,imon,t,tl,rhl,
     *  fl,xh,x,y,z,nmax,iflsfr,eqdp,idpco,bhh,opcct,opcsf,opc ,al(1)
      return
  100 format(/' t limits, fropc,imon:',4f10.5,1pe13.5,i4/' t =',0pf10.1,
     *  ' tl,rhl,fl,xh,x,y,z:',0p7f10.5/' nmax,iflsfr:',2i4/' /eqdpco/:'
     *  ,1p3e13.5,i4/' /bfacto/:',6e13.5//' opcct,opcsf,opc:',
     *  1p3e13.5,'  alscat =',e15.7)
      end
      subroutine sfropc(t,ane,opc)
c
c  sets rosseland or monochromatic (at frequency fropc) sfr opacity,
c  depending on whether imon = 0 or imon.ne.0
c
c  Note: constants have not been set consistently
c        Also need to check usage for ihvz = 4.
c
      implicit double precision (a-h,o-z)
      dimension fr(1)
      common/opcmtc/ tlm,dtlm,fropc,imon
      common/eqscnt/ anh0,anhe0,ihvz,idm(3)
      common/hviond/ xhvmn(10),xhv(125)
      common/cscat/ al(64)
c
      save
c
      fxp(x)=exp(max(-87.d0,x))
c  test for resetting of xhv
      go to (2,8,5), ihvz
c  set degrees of ionization from saha equation
    2 fci=4.82907d15*t**1.5/ane
      ri=fci*fxp(-8.866d4/t)
      xhv(6)=ri/(1+ri)
      ri=fci*fxp(-9.458d4/t)
      xhv(8)=ri/(1+ri)
      go to 8
c
    5 xhv(6)=xhv(24)
      xhv(8)=xhv(26)
c
    8 if(imon-1) 10,20,10
c
   10 call kappa(t,ane,opc)
      al(1)=0
      return
c
   20 fr(1)=fropc
      call gingab(1,fr,opc,t,ane)
      return
      end
      subroutine opacsi(rld,tld,xh,ak,rkr,rkt,rkx)
c
      implicit double precision (a-h,o-z)
      logical noder
      dimension akm(3,8),aki(3),zd1(4),zd2(4),zd3(4),zd4(4),zd5(4),
     *  zd6(4),zd7(4),zd8(4)
      common/ln10/ amm
      common/opccnt/ xhs,tsmn,tst,rhsmn,rhsmx,timx,rhimn,rhimx,sig,
     *  inopc,idmopc
      common/opccof/ iopccf,nrs,nts,nri,nti,mc,nrnts,nrpnts,nnns,nrnti,
     *  nrpnti,iopdum,akd(1)
      common/opcxdr/ akxa,akz
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
      common/dgnsrf/ idgsrf
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  initialize counters for points where log rho or log T are
c  out of range
c
      data icnrhs, icnrhi, icnts, icnti /0, 0, 0, 0/
c
      save
c  test if xh = xhs
   10 rl=rld
      tl=tld
      if(abs(xh-xhs).gt.1.e-5.or.tl.gt.tst) go to 30
c
c  test if rl and tl are in proper ranges
c
      if(rl.lt.rhsmn.or.rl.gt.rhsmx) then
	icnrhs=icnrhs+1
	if(icnrhs.le.100.and.istdpr.gt.0) then
          write(istdpr,100) rl,rhsmn,rhsmx
        else if(mod(icnrhs,50).eq.0.and.istdpr.gt.0) then
          write(istdpr,105) rl,rhsmn,rhsmx,icnrhs
        end if
	go to 23
      else if(tl.lt.tsmn) then
	icnts=icnts+1
	if(icnts.le.100.and.istdpr.gt.0) then
          write(istdpr,110) tl,tsmn,tst
        else if(mod(icnts,50).eq.0.and.istdpr.gt.0) then
          write(istdpr,115) tl,tsmn,tst,icnts
        end if
	go to 23
      end if
c
c  standard interpolation
c
      noder=.false.
      ak=surf2(rl,tl,nrs,nts,akd(1),akd(nrs+1),akd(nrpnts+1),
     *  nrs,akd(nrnts+nrpnts+1),sig,noder,rkr,rkt)
      akxa=0
      rkx=0
      return
c
c  linear extrapolation
c
   23 ak=opclne(akd(nrpnts+1),akd(1),akd(nrs+1),rl,tl,nrs,nrs,nts,
     *  rkr,rkt)
      rkx=0
      return
c
c  interior, varying x
c
c  set interpolation coefficients in xh, find m2
c
   30 noder=.false.
      initxi=-1
      call opcxin(akm,ak,akd(nnns+nrpnti+1),xh,3,3,mc,noder,rkx,mi1,mi2,
     *  iext,initxi)
      if(iext.eq.1.and.istdpr.gt.0) then
        write(istdpr,120) xh,(akd(nnns+nrpnti+1+3*(j-1)),j=1,mc)
      end if
c
c  test if rl and tl are in proper ranges
c
      if(rl.lt.rhimn.or.rl.gt.rhimx) then
        irtext=1
	icnrhi=icnrhi+1
	if(icnrhi.le.100.and.istdpr.gt.0) then
          write(istdpr,130) rl,rhimn,rhimx
        else if(mod(icnrhi,50).eq.0.and.istdpr.gt.0) then
          write(istdpr,135) rl,rhimn,rhimx,icnrhi
        end if
      else if(tl.lt.tst.or.tl.gt.timx) then
        irtext=1
	icnti=icnti+1
	if(icnti.le.100.and.istdpr.gt.0) then
          write(istdpr,140) tl,tst,timx
        else if(mod(icnti,50).eq.0.and.istdpr.gt.0) then
          write(istdpr,145) tl,tst,timx,icnti
        end if
      else
        irtext=0
      end if
c
c  find values at different compositions
c
      if(irtext.eq.0) then
c
c  standard interpolation
c
        iak=nnns+nrpnti+3*mc+nrnti*(mi1-2)+1
        izp=iak+nrnti*(mc+2*mi1-4)
        do 35 mm=mi1,mi2
        iak=iak+nrnti
        izp=izp+3*nrnti
        akm(1,mm)=surf2(rl,tl,nri,nti,akd(nnns+1),akd(nnns+nri+1),
     *    akd(iak),nri,akd(izp),sig,noder,akm(2,mm),akm(3,mm))
        if(abs(akm(1,mm)).gt.100) then
          idgsrf=1
          akm(1,mm)=surf2(rl,tl,nri,nti,akd(nnns+1),akd(nnns+nri+1),
     *      akd(iak),nri,akd(izp),sig,noder,akm(2,mm),akm(3,mm))
          idgsrf=0
          write(istdou,150) rl, tl, xh, mm, akm(1,mm)
          if(istdou.ne.istdpr.and.istdpr.gt.0)
     *       write(istdpr,150) rl, tl, xh, mm, akm(1,mm)
          stop 'opact'
        end if
c
c  ** extra diagnostics
c
        if(idgopc.eq.1.and.istdpr.gt.0) 
     *    write(istdpr,34091) mm,akd(nnns+nrpnti+1+3*(mm-1)),
     *    nnns,nrpnti,nri,nrnti,iak,izp,akm(1,mm)
34091   format(' mm,comp(1,mm),nnns,nrpnti,nri,nrnti,iak,izp,akm =',
     *    i4,1pe15.7,6i6,e15.7)
   35   continue
c
      else
c
c  linear extrapolation
c
        iak=nnns+nrpnti+3*mc+nrnti*(mi1-2)+1
        do 47 mm=mi1,mi2
        iak=iak+nrnti
   47   akm(1,mm)=opclne(akd(iak),akd(nnns+1),akd(nnns+nri+1),rl,tl,
     *    nri,nri,nti,akm(2,mm),akm(3,mm))
c
      end if
c
c  interpolation in x
c
      initxi=0
      noder=.false.
c
      do 55 i=1,3
      call opcxin(akm(i,1),aki(i),akd(nnns+nrpnti+1),xh,3,3,mc,noder,
     *  rkx,mi1,mi2,iext,initxi)
   55 noder=.true.
      ak=aki(1)
      rkr=aki(2)
      rkt=aki(3)
      akxa=rkx
      rkx=amm*xh*rkx
      return
  100 format(' env. rl =',f10.3,' outside range:',2f10.3)
  105 format(' env. rl =',f10.3,' outside range:',2f10.3,
     *  '  count no.',i6)
  110 format(' env. tl =',f10.3,' outside range:',2f10.3)
  115 format(' env. tl =',f10.3,' outside range:',2f10.3,
     *  '  count no.',i6)
  120 format(' xh =',f8.5,' outside range. comp =',9f7.3/
     *  (18f7.3))
  130 format(' int. rl =',f10.3,' outside range:',2f10.3)
  135 format(' int. rl =',f10.3,' outside range:',2f10.3,
     *  '  count no.',i6)
  140 format(' int. tl =',f10.3,' outside range:',2f10.3)
  145 format(' int. tl =',f10.3,' outside range:',2f10.3,
     *  '  count no.',i6)
  150 format(//' ***** error in s/r opact, for log rho =',f10.5,
     *  ' log T =',f10.5,'  X =',f10.5/
     *  ' mm =',i5,' resulting log kappa =',1pe13.5/
     *  ' Execution terminated')
      end
      function opclne(akk,rho,t,rl,tl,ia1,
     *  nrm,nsm,dar,dat)
c
      implicit double precision (a-h,o-z)
      dimension  akk(ia1,1),rho(1),t(1)
c
      save
c  find extrapolation points
      do 10 i=1,nrm
      ir2=i
      if(rl.lt.rho(i)) go to 15
   10 continue
   15 if(ir2.eq.1) ir2=2
c
      do 20 j=1,nsm
      jt2=j
      if(tl.lt.t(j)) go to 25
   20 continue
   25 if(jt2.eq.1) jt2=2
      ir1=ir2-1
      jt1=jt2-1
c
      fct=1/(rho(ir2)-rho(ir1))/(t(jt2)-t(jt1))
      ct2=akk(ir1,jt1)*(rho(ir2)-rl)+akk(ir2,jt1)*(rl-rho(ir1))
      ct1=akk(ir1,jt2)*(rho(ir2)-rl)+akk(ir2,jt2)*(rl-rho(ir1))
c
      opclne=fct*(ct2*(t(jt2)-tl)+ct1*(tl-t(jt1)))
c..      write(6,*) ir1, ir2, rl, rho(ir1), rho(ir2)
c..      write(6,*) jt1, jt2, tl, t(jt1), t(jt2)
c..      write(6,*) akk(ir1,jt1), akk(ir2,jt1), akk(ir1,jt2), 
c..     *  akk(ir2,jt2), opclne
      dar=fct*((akk(ir2,jt1)-akk(ir1,jt1))*(t(jt2)-tl)
     *  + (akk(ir2,jt2)-akk(ir1,jt2))*(tl-t(jt1)))
      dat=fct*(ct1-ct2)
      return
      end
      subroutine opcxin(akk,aki,comp,xi,iak,ic,mc,noder,daki,mi1,mi2,
     *  iext,initxi)
c
c  routine for opacity interpolation in x, to value xi.
c
c  the type of interpolation depends on the value of iwdopc in
c  common/opctcl/
c 
c  iwdopc = 0: use old scheme with parabolic interpolation.
c    Note: the cluster of points is chosen to be as close as possible
c    to target x. This leads to switch of cluster in the middle of
c    mesh interval, and hence to discontinuous interpolating function.
c    A little unfortunate!
c  iwdopc = 2: use linear interpolation
c  iwdopc = 3: use 4-point Lagrangian interpolation.
c
c  if initxi = 0 interpolation coefficients are only initialized
c  if the value of xi differs from the previous value.
c
c  if initxi = -1, call only initializes coefficients, without
c  carrying out interpolation.
c
c  if xi is outside range of x, iext is returned as 1, and linear
c  extrapolation is used.
c
c  returns mi1 and mi2 such that the interpolation requires
c  opacities at compositions comp(1,m), m = mi1, ..., mi2
c
c  modified 30/12/1984, to include initxi.
c
      implicit double precision (a-h,o-z)
      logical noder, setder
      dimension akk(iak,mc),comp(ic,mc),x(50),c(4),dc(4)
      common/opctcl/ alamop,zatmop,tlmopc,dtmopc,fdgopl,iwdopc,
     *  iopacm,ifdgop
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      data xp/2./,ids/0/,mi1p,mi2p /-1,-1/
      data x/50*0./
      data icntng /0/
c
      setder=.not.noder
c
c  test for skipping initialization
c
      if(initxi.eq.0.and.
     *  xi.eq.xp.and.x(mc).eq.comp(1,mc).and.(noder.or.ids.eq.1))
     *  then
        mi1=mi1p
        mi2=mi2p
        go to 40
      end if
c
c  find interpolation coefficients
c  comp(1,m) must be monotonic
c
      iext=0
      do 2 m=1,mc
    2 x(m)=comp(1,m)
      xrng=x(mc)-x(1)
      if(xrng.eq.0) go to 50
c
      if(idgopc.eq.1.and.istdpr.gt.0) write(istdpr,105) xi
c
c  test for xi out of range of table
c
      if(xi.lt.0) then
c
c  negative xi, resulting from core hydrogen exhaustion
c
        if(icntng.le.50.and.istdpr.gt.0) write(istdpr,107) xi
	icntng=icntng+1
	xi=1.e-20
      end if
c
      if((x(1)-xi)*xrng.gt.0) then
c  xi to the left of x(1)
        m1=1
        m2=2
        iext=1
      else if((x(mc)-xi)*xrng.lt.0) then
c  xi to the right of x(mc)
        m2=mc
        m1=mc-1
        iext=1
      else
        iext=0
      end if
c
c  test for linear interpolation or extrapolation
c
      if(iext.eq.1) then
c
c  use linear extrapolation
c
        mi1=m1
        mi2=m2
        dxa=x(m2)-x(m1)
        c(1)=(x(m2)-xi)/dxa
        c(2)=(xi-x(m1))/dxa
        if(setder) then
	  ids=1
          dc(2)=1/dxa
          dc(1)=-dc(2)
        end if
c
        go to 35
      end if
c
c  set proper interval for interpolation
c  m2 is set to the index of the first point at or following the
c  target x.
c
      do 10 m=1,mc
      m2=m
      if((xi-x(m))*xrng.le.0) go to 12
   10 continue
c
c  test for type of interpolation, depending on the value of iwdopc
c
   12 if(iwdopc.eq.0) then
c
c  use old scheme, with parabolic interpolation.
c  =============================================
c  test for resetting m2 to use cluster of points closest to target
c  (note that this leads to switch of cluster in the middle of interval
c  and hence to discontinuous interpolating function. Bloody stupid!)
c
        if(m2.gt.2.and.
     *    (m2.eq.mc.or.abs(xi-x(m2)).gt.abs(xi-x(m2-1)))) m2=m2-1
c
        m1=m2-1
        m3=m2+1
        mi1=m1
        mi2=m3
        dxa=x(m2)-x(m1)
        dxb=x(m3)-x(m1)
        dxc=x(m3)-x(m2)
c
        dx1=(xi-x(m1))/dxc
        dx2=(xi-x(m2))/dxb
        dx3=(xi-x(m3))/dxa
c
        c(1)=dx2*dx3
        c(2)=-dx1*dx3
        c(3)=dx1*dx2
        if(setder) then
c
          ids=1
          dc(1)=dx2/dxa+dx3/dxb
          dc(2)=-dx1/dxa-dx3/dxc
          dc(3)=dx1/dxb+dx2/dxc
c
        end if
c
      else if(iwdopc.eq.2) then
c
c  use linear interpolation
c  ========================
c
        mi1=m2-1
        mi2=m2
        dxa=x(m2)-x(mi1)
        c(1)=(x(m2)-xi)/dxa
        c(2)=(xi-x(mi1))/dxa
        if(setder) then
	  ids=1
          dc(2)=1/dxa
          dc(1)=-dc(2)
        end if
c
      else if(iwdopc.eq.3) then
c
c  use 4-point Lagrangian interpolation
c  ====================================
c
c  test for sufficient number of points
c
        if(mc.lt.4) then
          write(istdou,130) mc
          if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,130) mc
          stop 'Stop 1 in opcxin'
        end if
c
c  test for resetting m2 to allow for edge effects
c
        if(m2.eq.2) then
          m2=3
        else if(m2.eq.mc) then
          m2=m2-1
        end if
c
        mi1=m2-2
        mi2=m2+1
c
c  start setting coefficients
c
        do 30 i=1,4
        c(i)=1
        den=1
        xti=x(mi1+i-1)
        do 22 j=1,4
        if(j.ne.i) then
          xtj=x(mi1+j-1)
          if(xtj.eq.xti) go to 50
          c(i)=c(i)*(xi-xtj)
          den=den*(xti-xtj)
        end if
   22   continue
        c(i)=c(i)/den
c
c  test for derivatives
c
        if(setder) then
	  ids=1
          dc(i)=0
          do 26 k=1,4
          if(k.ne.i) then
            prod=1
            do 24 j=1,4
            if(j.ne.i.and.j.ne.k) prod=prod*(xi-x(mi1+j-1))
   24       continue
            dc(i)=dc(i)+prod
          end if
   26     continue
          dc(i)=dc(i)/den
        end if
c
   30   continue
c
      else
        write(istdou,120) iwdopc
        if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,120) iwdopc
        stop 'Stop 2 in opcxin'
      end if
c
c  ** extra output
c
   35 if(idgopc.eq.1) then
	ncoeff=mi2-mi1+1
        if(istdpr.gt.0) then
	  write(istdpr,110) mi1,mi2,(c(i),i=1,ncoeff)
	  if(setder) write(istdpr,115) (dc(i),i=1,ncoeff)
        end if
      end if
c
c  end of initialization
c
   40 xp=xi
      mi1p=mi1
      mi2p=mi2
c
      if(initxi.eq.-1) return
c
c   interpolate to xi
c
      aki=0
      if(setder) daki=0
      i=0
      do 45 m=mi1,mi2
      i=i+1
      aki=aki+c(i)*akk(1,m)
      if(setder) daki=daki+dc(i)*akk(1,m)
   45 continue
c
      return
c
   50 write(istdou,100) (comp(1,m),m=1,mc)
      if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,100) 
     *  (comp(1,m),m=1,mc)
      xi=-1
      return
  100 format(//1x,10(1h*),' comp is not monotonic. comp(1,1-mc)='/
     *  (1p10e13.5))
  105 format('  set X-interpolation coefficients for X =',f12.5)
  107 format(' *** warning in s/r opcxin. X =',f12.7,' reset to 1.e-20')
  110 format(' mi1, mi2,  coefficients:',2i4,1p4e15.7)
  115 format(' derivative coefficients:',8x,1p4e15.7)
  120 format(//' ***** error in s/r opcxin. iwdopc =',i4,
     *  ' is illegal')
  130 format(//' ***** error in s/r opcxin. for iwdopc = 3, more than',
     *  i4,' points are required')
      end
      subroutine opinsi
c
c  initialization routine for spline opacity interpolation
c
c  maximum table sizes set to fit Turck-Chieze version of LAOL tables
c
      parameter(nrhlmx=33,ntlmx=58,ncmpmx=6)
      parameter(idimzp=3*nrhlmx*ntlmx, idimtm=2*nrhlmx+ntlmx)
      implicit double precision (a-h,o-z)
      logical noder
      common/opccnt/ xhs,tsmn,tst,rhsmn,rhsmx,timx,rhimn,rhimx,sig,
     *  inopc,idmopc
      common/work/ akk(nrhlmx,ntlmx,ncmpmx),t(ntlmx),rho(nrhlmx),
     *  nt(nrhlmx),ns(nrhlmx),comp(3,ncmpmx)
      common/sooner/ ak(nrhlmx,ntlmx),zp(idimzp),temp(idimtm)
      common/opccof/ iopccf,nrs,nts,nri,nti,mc,nrnts,nrpnts,nnns,
     *  nrnti,nrpnti,iopdum,akd(1)
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
c..c  apparently vectorization or concurrency causes problems
c..c  with ranges in tables (ir1, ir2, .....).
c..c  thus switch off in this routine.
c..c
c..cvd$r novector
c..cvd$r noconcur
c
c  read tables
c
      ia1=nrhlmx
      ia2=ntlmx
   10 call opcrd(akk,t,rho,comp,nt,ns,nrm,nsm,mc,ia1,ia2,inopc)
c  extrapolate
c..      if(idgopc.eq.1) then
c..        write(6,10191)
c..        do 10091 nr=1,nrm
c..10091   write(6,10192) nr,(akk(nr,j,3),j=1,nsm)
c..      end if
      call opcext(akk,t,rho,comp,nt,ns,ia1,ia2,mc,nrm,nsm)
c..      if(idgopc.eq.1) then
c..        write(6,10193)
c..        do 10093 nr=1,nrm
c..10093   write(6,10192) nr,(akk(nr,j,3),j=1,nsm)
c..      end if
10191 format(/' before call of opcext. nr, akk(1-nsm):')
10192 format(i3,1p7e11.3/(3x,7e11.3))
10193 format(/' after  call of opcext. nr, akk(1-nsm):')

c
c
c  surface opacities
c
c  find range of restricted table
      call opcres(t,rho,nrm,nsm,rhsmn,rhsmx,tsmn,tst,ir1,ir2,ier1,ier2,
     *  it1,it2,iet1,iet2)
c
      if(idgopc.eq.1.and.istdpr.gt.0) write(istdpr,104) 
     *  ier1,ier2,iet1,iet2
c  interpolate in x
      noder=.true.
      iaak=ia1*ia2
c
      initxi=1
      do 20 i=ir1,ir2
      do 20 j=it1,it2
      call opcxin(akk(i,j,1),ak(i,j),comp,xhs,iaak,3,mc,noder,dak,
     *  mi1,mi2,iext,initxi)
c..   20 initxi=0
      initxi=0
   20 continue
c
      if(iext.eq.1.and.istdpr.gt.0) write(istdpr,100) 
     *  xhs,(comp(1,m),m=1,mc)
c
c  find interpolation coefficients
c
      nrmr=ir2-ir1+1
      nsmr=it2-it1+1
      if(idgopc.eq.1.and.istdpr.gt.0) write(istdpr,105) 
     *  nrmr, nsmr, ir1, it1
      call surf1(nrmr,nsmr,rho(ir1),t(it1),ak(ir1,it1),ia1,
     *  zd1,zd2,zd3,zd4,zd5,zd6,zd7,zd8,zp,temp,sig)
      nrntr=nrmr*nsmr
c  store surface values in akd
   30 ir1=ir1+ier1
      ir2=ir2-ier2
      it1=it1+iet1
      it2=it2-iet2
      irs1=ir1
      irs2=ir2
      its1=it1
      its2=it2
      nrs=ir2-ir1+1
      nts=it2-it1+1
      nrnts=nrs*nts
      nrpnts=nrs+nts
      nrntr=nrmr*nsmr
      nnns=4*nrnts+nrpnts
c  test that size of opccof is sufficient
      iopccs=nnns+10
      if(iopccf.lt.iopccs) go to 40
c
      ia=0
      do 32 i=ir1,ir2
      ia=ia+1
   32 akd(ia)=rho(i)
      do 34 j=it1,it2
      ia=ia+1
   34 akd(ia)=t(j)
c
      izzp=ier1-ir1+1+nrmr*(iet1-1)-nrntr
      do 36 j=it1,it2
      izzp=izzp+nrmr
      do 36 i=ir1,ir2
      ia=ia+1
      akd(ia)=ak(i,j)
      izkp=izzp+i
      ib=ia
      do 36 k=1,3
      izkp=izkp+nrntr
      ib=ib+nrnts
   36 akd(ib)=zp(izkp)
c
c  interior opacities
c
c  find range of restricted table
c
   40 call opcres(t,rho,nrm,nsm,rhimn,rhimx,tst,timx,
     *  ir1,ir2,ier1,ier2,it1,it2,iet1,iet2)
c
      if(idgopc.eq.1.and.istdpr.gt.0) write(istdpr,104) 
     *  ier1,ier2,iet1,iet2
c  store rhi and ti
      irr1=ir1+ier1
      irr2=ir2-ier2
      itt1=it1+iet1
      itt2=it2-iet2
c
      nrmr=ir2-ir1+1
      nsmr=it2-it1+1
      nrntr=nrmr*nsmr
      nri=irr2-irr1+1
      nti=itt2-itt1+1
      nrnti=nri*nti
      nrpnti=nri+nti
c  test that size of opccof is sufficient
      iopcct=iopccs+3*mc+4*nrnti*mc+nrpnti
      if(iopccf.lt.iopcct) then
        write(istdou,110) iopccf,iopcct,
     *    irs1,irs2,its1,its2,nrs,nts,nrnts,iopccs,
     *    irr1,irr2,it1,itt2,nri,nti,nrnti,mc
        if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,110) 
     *    iopccf,iopcct,irs1,irs2,its1,its2,nrs,nts,nrnts,iopccs,
     *    irr1,irr2,it1,itt2,nri,nti,nrnti,mc
        stop 'Stop in opinsi'
      end if
c
   41 idopc=iopccf-iopcct
      if(istdpr.gt.0) write(istdpr,115) iopcct,idopc
      ia=ib
      do 42 i=irr1,irr2
      ia=ia+1
   42 akd(ia)=rho(i)
      do 44 j=itt1,itt2
      ia=ia+1
   44 akd(ia)=t(j)
c  store comp
      do 46 m=1,mc
      do 46 i=1,3
      ia=ia+1
   46 akd(ia)=comp(i,m)
c  find interpolation coefficients and store them
      ibmia=nrnti*(mc-3)
   50 do 60 m=1,mc
      if(idgopc.eq.1.and.istdpr.gt.0) write(istdpr,107) 
     *  m, nrmr, nsmr, ir1, it1
      call surf1(nrmr,nsmr,rho(ir1),t(it1),akk(ir1,it1,m),ia1,
     *  zd1,zd2,zd3,zd4,zd5,zd6,zd7,zd8,zp,temp,sig)
      izzp=1-ir1+nrmr*(iet1-1)-nrntr
      ibmia=ibmia+2*nrnti
      do 60 j=itt1,itt2
      izzp=izzp+nrmr
      do 60 i=irr1,irr2
      ia=ia+1
      izkp=izzp+i
      ib=ia+ibmia
      akd(ia)=akk(i,j,m)
      do 60 k=1,3
      ib=ib+nrnti
      izkp=izkp+nrntr
   60 akd(ib)=zp(izkp)
      return
  100 format(///1x,20(1h*),' surface x =',f10.3,
     *  ' is outside range of comp.'/21x,' comp =',10f10.4/
     *  (29x,10f10.4))
  102 format(/' call opcres at surface. ier1, ier2, iet1, iet2 =',
     *  4i4)
  104 format(/' call opcres in interior. ier1, ier2, iet1, iet2 =',
     *  4i4)
  105 format(' call surf1 with nrmr =',i5,'  nsmr =',i5/
     *  '  at ir1 =',i5,'  it1 =',i5)
  107 format(' call surf1 at m =',i3,' with nrmr =',i5,'  nsmr =',i5/
     *  '  at ir1 =',i5,'  it1 =',i5)
  110 format(///1x,130(1h*)//' insufficient store supplied for',
     *  ' common/opccof/.'//i10,' real places supplied.'//
     *  i10,' real places needed.'//
     *  ' in surface part'/
     *  ' ir1 =',i4,' ir2 =',i4,' it1 =',i4,' it2 =',i4/
     *  ' nrs =',i4,' nts =',i4,' nrnts =',i6,'    iopccs =',i6/
     *  ' in interior part'/
     *  ' ir1 =',i4,' ir2 =',i4,' it1 =',i4,' it2 =',i4/
     *  ' nri =',i4,' nti =',i4,' nrnti =',i6,' mc =',i4///
     *  ' execution terminated by s/r opcscf')
  115 format(//i10,' real places needed in common/opccof/.',
     *  i10,' places not used'//)
      end
      subroutine opcrd(akk,tt,rho,comp,nt,ns,nrm,nsm,mc,ia1,ia2,inopc)
c
      implicit double precision (a-h,o-z)
      parameter(ntabmx=500)
      real comps, tts, rhos, akks
      dimension akk(ia1,ia2,1),tt(1),rho(1),comp(3,1),ns(1),nt(1)
      dimension comps(3), tts(ntabmx), akks(ntabmx)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
      common/copcxz/ ixzopc,ixzdum,xzopc(2,20)
c
      save
c
      i=0
c
      if(idgopc.eq.1.and.istdpr.gt.0) 
     *  write(istdpr,*) 'start reading from d/s  ',inopc
c  read opacity
    2 i=i+1
      read(inopc,end=20) nrmrd,(comps(j),j=1,3)
      do 3 j=1,3
      xzopc(1,j)=comps(1)
      xzopc(2,j)=comps(3)
    3 comp(j,i)=comps(j)
c
      nrm=nrmrd
      if(idgopc.eq.1.and.istdpr.gt.0) 
     *   write(istdpr,*) 'read nrm = ',nrm,'  comp =',(comp(j,i),j=1,3)
      read(inopc) nsm,(tts(n),n=1,nsm)
c
c  test for sufficient single-precision storage space
c
      if(nrm.gt.ntabmx.or.nsm.gt.ntabmx) then
	write(istdou,105) nrm, nsm, ntabmx
	if(istdou.ne.istdpr.and.istdpr.gt.0)
     *    write(istdpr,105) nrm, nsm, ntabmx
	stop 'Stop 1 in opcrd'
      end if
c
      do 4 n=1,nsm
    4 tt(n)=tts(n)
c
      if(idgopc.eq.1.and.istdpr.gt.0) write(istdpr,*) 'read nsm =', nsm
      do 8 nr=1,nrm
      do 5 j=1,nsm
    5 akk(nr,j,i)=0
      read(inopc) rhos,nnt,nns,(akks(j),j=1,nnt)
      rho(nr)=rhos
      do 6 j=1,nnt
    6 akk(nr,nns+j,i)=akks(j)
c
      ns(nr)=nns
    8 nt(nr)=nnt
      go to 2
   20 mc=i-1
c
      ixzopc=mc
c
c  test that anything has been read
c
      if(mc.le.0) then
	write(istdou,110)
	if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,110)
	stop 'Stop 2 in opcrd'
      end if
      return
  100 format(i3,3f10.5)
  105 format(//' **** Insufficient single-precision table space in',
     *  ' s/r opcrd'/
     *         '      nrm =',i5,'  nsm =',i5,'  ntabmx =',i5)
  110 format(/'**** No tables read in s/r opcrd'/
     *        '     Execution terminated')
      end
      subroutine opcext(akk,tt,rho,comp,nt,ns,ia1,ia2,mc,nrm,nsm)
c
c  extrapolate opacities outside given table
c
      implicit double precision (a-h,o-z)
      dimension akk(ia1,ia2,mc),comp(3,mc),tt(1),rho(1),nt(1),ns(1),
     *  ake(500)
c
      save
c
   10 nsm1=nsm+1
      nrm1=nrm+1
      do 30 m=1,mc
      akt=log10(0.2*(1+comp(1,m)))
      do 20 n1=1,nsm
      ne=nsm1-n1
      i1=nrm1
   12 i1=i1-1
      if(akk(i1,ne,m).eq.0) go to 12
      if(i1.eq.nrm) go to 20
c
c  test for extrapolation point outside table.
c  If so, arbitrarily set to zero, for consistency with first call
c  of routine. This must be fixed up
c
c  ********************************************************************
c
      if(ne.eq.nsm) then
        da=akk(i1,ne,m)
      else
        da=akk(i1,ne,m)-akk(i1,ne+1,m)
      end if
c
c  ********************************************************************
c
      i1=i1+1
      do 15 i=i1,nrm
c
c  test for extrapolation point outside table.
c  If so, arbitrarily set to zero, for consistency with first call
c  of routine. This must be fixed up
c
c  ********************************************************************
c
      if(ne.eq.nsm) then
        akk(i,ne,m)=da
      else
        akk(i,ne,m)=da+akk(i,ne+1,m)
      end if
   15 continue
c
c  ********************************************************************
c
   20 continue
c
      do 22 n=1,nsm
      i1=0
   21 i1=i1+1
      if(akk(i1,n,m).eq.0) go to 21
   22 ake(n)=akk(i1,n,m)
c
      do 30 n1=1,nrm
      nr=nrm1-n1
      nns=nsm-ns(nr)-nt(nr)
      if(nns.eq.0) go to 30
      a2=rho(nr+1)-rho(nr+2)
      a1=(rho(nr)-rho(nr+2))/a2
      a2=(rho(nr)-rho(nr+1))/a2
      do 28 n2=1,nns
      n=nsm1-n2
      ak1=akk(nr+1,n,m)
      ak2=akk(nr+2,n,m)
      akt1=min(akt,ake(n))
      ak2=max(ak2,(4*ak1-akt1)/3)
      ak2=max(akt1,a1*ak1-a2*ak2)
   28 akk(nr,n,m)=ak2
   30 continue
      return
      end
      subroutine opcres(t,rho,nrm,nsm,rhmn,rhmx,tmn,tmx,
     *  ir1,ir2,ier1,ier2,it1,it2,iet1,iet2)
      implicit double precision (a-h,o-z)
      dimension t(1),rho(1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c..c
c..c  apparently vectorization or concurrency causes problems
c..c  with ranges in tables (ir1, ir2, .....).
c..c  thus switch off in this routine.
c..c
c..cvd$r novector
c..cvd$r noconcur
c..c
c
c  Modified 14/8/90, to fix up setting of ranges if desired range
c  is larger than table range.
c  
c  test for ranges (rhmn, rhmx) and (tmn, tmx) fitting within
c  table ranges
c
      if(rhmn.lt.rho(1)) then
	write(istdou,100) rhmn,rho(1)
	if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,100) 
     *    rhmn,rho(1)
	rhmn=rho(1)
      end if
c
      if(rhmx.gt.rho(nrm)) then
	write(istdou,110) rhmx,rho(nrm)
	if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,110) 
     *    rhmx,rho(nrm)
	rhmx=rho(nrm)
      end if
c
      if(tmn.lt.t(1)) then
	write(istdou,120) tmn,t(1)
	if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,120) tmn,t(1)
	tmn=t(1)
      end if
c
      if(tmx.gt.t(nsm)) then
	write(istdou,130) tmx,t(nsm)
	if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,130) 
     *    tmx,t(nsm)
	tmx=t(nsm)
      end if
c
c  set ranges of indices in tables
c
      ir1=0
      ir2=0
      do 12 i=1,nrm
      if(rho(i).le.rhmn) ir1=i
   12 continue
      do 14 i=1,nrm
      if(rho(i).lt.rhmx) ir2=i+1
   14 continue
c
      it1=0
      it2=0
      do 22 j=1,nsm
      if(t(j).le.tmn) it1=j
   22 continue
      do 24 j=1,nsm
      if(t(j).lt.tmx) it2=j+1
   24 continue
c
      ier1=min0(ir1-1,3)
      ier2=min0(nrm-ir2,3)
      iet1=min0(it1-1,3)
      iet2=min0(nsm-it2,3)
      ir1=ir1-ier1
      ir2=ir2+ier2
      it1=it1-iet1
      it2=it2+iet2
      return
  100 format(//1x,20(1h*),' rhomn =',f10.3,' < rho(1) =',f10.3/
     *                21x,' rhomn has been reset to rho(1)'/)
  110 format(//1x,20(1h*),' rhomx =',f10.3,' > rho(nrm) =',f10.3/
     *                21x,' rhomx has been reset to rho(nrm)'/)
  120 format(//1x,20(1h*),' tmn =',f10.3,' < t(1) =',f10.3/
     *                21x,' tmn has been reset to t(1)'/)
  130 format(//1x,20(1h*),' tmx =',f10.3,' > t(nsm) =',f10.3/
     *                21x,' tmx has been reset to t(nsm)'/)
      end
      subroutine opfdtb(ifdgop, idopfd, tl, extl, dextl)
c
c  set correction to opacity by interpolating in table
c  read in from d/s idopfd. Table is read in first call,
c  or if idopfd is changed
c
c  Original version: 11/9/96
c
      implicit double precision(a-h, o-z)
      parameter(ntbmax=5000)
      dimension tlogtb(ntbmax), dkaptb(2,ntbmax), dkapin(2)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      data idopfp /-1/
      save
c
c  test for reading tables and setting derivative
c
      if(idopfd.ne.idopfp) then
	idopfp=idopfd
	call skpcom(idopfd)
	n=1
   10   read(idopfd,*,end=15) tlogtb(n),dkaptb(1,n)
	n=n+1
	if(n.le.ntbmax) then
	  go to 10
        else
	  if(istdpr.gt.0) write(istdpr,110) ntbmax
        end if
c
   15   ntb=n-1
c
	if(istdpr.gt.0) write(istdpr,*) 'End reading table; ntb =',ntb
	if(ntb.le.0) stop 'No data read in opfdtb'
c
	tltbmn=min(tlogtb(1),tlogtb(ntb))
	tltbmx=max(tlogtb(1),tlogtb(ntb))
c
c  set derivative
c
	call derive(tlogtb,dkaptb(1,1),dkaptb(2,1),ntb,2,2,1,1)
c
      end if
c
c  interpolate and set output
c
      if(tl.le.tltbmn.or.tl.ge.tltbmx) then
        extl=0
        dextl=0
      else	
        call lir(tl,tlogtb,dkapin,dkaptb,2,2,ntb,1,inter)
        extl=dkapin(1) 
        dextl=dkapin(2)
      end if
c
      return
  110 format(/' ***** Warning in s/r opfdtb: ',
     *  ' number of points in table exceeds limit of',i6)
      end
      subroutine kappwd(rl,tl,xh,yh,ak,akr,akt,akxa,akz)
c  dummy subroutine
      implicit double precision (a-h, o-z)
      return
      end
      subroutine opinwd
c  dummy subroutine
      return
      end
      subroutine iging
c  dummy subroutine
      return
      end
      subroutine kappa(t,ane,opc)
c  dummy subroutine
      implicit double precision (a-h, o-z)
      return
      end
      subroutine gingab(idum,fr,opc,t,ane)
c  dummy subroutine
      implicit double precision (a-h, o-z)
      return
      end
      subroutine opacat(rl,tl,xh,ak,akr,akt,akx)
c
c  Dummy subroutine, replacing entry opacat in opact
c
      implicit double precision (a-h, o-z)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'Entry opacat not allowed in present verion'
      stop 'Stop in opacat'
      end
