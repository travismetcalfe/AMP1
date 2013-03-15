      subroutine rhs(x,y,zk,zz,dzdy,ap,aq,f,fd,alam,alamd,h,hd,iz,
     .  ifd,ihd,ialam,n,iter)
c
c  right hand side subroutine for evolution.
c
c  Note: iter = 0 when rhs is called from tnrkt at previous
c  time step, or from mshstr. A consequence of this is that
c  convection zone boundaries are not set.
c  This may cause problems if convective overshoot is
c  used, since this requires that the boundary of the
c  convective envelope is known.
c#ai#  This must be checked.
c
c  On 16/1/2000, option of iter = -1 introduced, to allow setting
c  right-hand sides and convective boundaries, but not derivatives.
c
c  iter1 (in common/noiter/) should be set to the actual iteration number
c  in the routine calling tnrk. When called outside the iteration 
c  (for setting variables for printing, say), iter1 must be set to zero.
c
c  ***************************************************************************
c  
c  Notes on convection-zone boundaries:
c  ===================================
c  
c  Information about the convection-zone boundaries are stored 
c  in commons /convpt/ and /convvr/.
c  In /convpt/ the relevant quantities are:
c    inc: number of convection zones
c    nf, nl: meshpoint number at convection-zone boundaries;
c      nf(i) is the first meshpoint (counted from the surface), and
c      nl(i) is the last, inside or on the edge of the i-th
c      convection zone. 
c    frcf, frcl: linear-interpolation weights to determine quantities
c      at convection-zone boundary.
c      Let nnf = nf(i) and frc = frcf(i), for the i-th zone.
c      Then at the upper edge of the zone, for any variable a(n) the
c      interpolated value is frc*a(nnf)+(1-frc)*a(nnf-1).
c      Similarly, the interpolated value at the lower edge of the zone is
c      frc*a(nnl)+(1-frc)*a(nnl+1), with nnl = nl(i) and frc = frcl(i)
c    inccor: no. of convective core (added 30/12/99)
c  In /convvr/ are stored values at the upper and lower edges of the
c  convection zones, linearly interpolated.
c  At the upper edge:
c    cvvarf(1,i) = m/M
c    cvvarf(2,i) = r/R
c    cvvarf(3,i) = p
c    cvvarf(4,i) = rho
c    cvvarf(5,i) = T
c    cvvarf(6,i) = kappa
c    cvvarf(7,i) = Hp (added 30/12/99)
c    cvvarf(8,i) = nabla_ad (added 12/8/04)
c    cvvarf(9,i) = nabla_rad (added 12/8/04)
c    cvvarf(10,i) = nabla_ac (added 12/8/04)
c    cvvarf(11,i) = nabla_X (added 12/8/04)
c  Similarly, variables at lower edge are stored in cvvarl.
c  Also, the variables at all mesh points are stored in cvvars(.,n)
c  
c  *************************************************************************
c
c
c  modified 4/1/1985 to use log(r/1.d11) and log(l/1.d33) as
c  dependent variables.
c
c  modified 10/3/1986 to straighten up determination of convective
c  zone boundaries  
c
c  21/9/87: implementing modifications from RECKU
c
c  11/5/89: implement possibility of setting entropy terms to zero,
c  if inentr (in common/rhcn/) is zero.
c
c  Modified 5/6/89, adding option for smooth transition in temperature
c     gradient between specified T(tau) relation in atmosphere and
c     diffusion approximation in interior.
c
c  Modified 11/4/91, implementing (fairly roughly) a simplified
c  description of convective overshoot.
c  Note: has to be checked carefully for evolution models.
c
c  Modified 2/8/91, to incorporate detailed treatment of CNO cycle.
c
c  Modified 3/11/94 by MJM to include jcnvos = 3 option for
c  smooth transition in s/r oversh
c
c  Modified 19/4/95, to correct error in connection with mixing
c  at boundary of convective core. Previously, the variables
c  cqc and cqcp were erroneously called cqxc and cqxcp in body
c  of routine (but not in common/crxcor/), and hence were not set.
c
c  Modified 19/7/95, to use general heavy-element abundance
c  as set in common/heavy/. 
c
c  Modified 4/8/96, taking out option of writing data to disk with
c  ifwrt = 1 (but retaining ifwrt in common/rhcn/ for consistency)
c
c  Modified 5/8/96, adding option imixc1 = 1, for setting rhs
c  consistently with treatment of convective-core mixing (we hope)
c
c  Modified 20/8/96, for option of including centrifugal acceleration
c
c  Corrected 5/8/97, by including engenr.nnz.d.incl, to set nspdmx etc
c
c  Modified 18/11/98, adding ptpg to call of mixlng (for treatment 
c  of turbulent pressure)
c
c  Modified 5/1/00, to add more general treatment of convective overshoot
c  from convective envelope and core.
c
c  Modified 19/4/00, to block resetting of composition time derivatives
c  if core was not convective in previous time step (a somewhat 
c  desparate measure)
c
c  Modified 13/7/00, eliminating entropy term in energy generation
c  if it becomes unreasonably large in small convective core
c  (typically as a result of sudden onset of core convection during
c  main-sequence evolution). Parameters etc. in this could do with
c  fine tuning.
c
c  Modified 7/4/01, including derivatives of atmospheric factor in 
c  setting derivatives of f(3) (temperature gradient).
c
c  Modified 5/8/02, including 4He burning
c
c  Modified 8/10/02, adding resetting of zh(n) in region where
c  hydrogen is exhausted.
c
c  Modified 22/6/03, removing thermodynamic and opacity derivatives
c  wrt Y. Still would benefit from further thought.
c
c  Modified 11/9/03, adding z(.) to output on unit 77 when idgn77 .ge. 1
c
c  Modified 18/6/04, increasing cut-off for call of ttau from log T = 4.3
c  to 4.6
c
c  Modified 22/1/08, setting rate of composition change in region of
c  shrinking convective core taking centralization into account, to
c  be consistent with cmpcvc. So far implemented only for imixc6 = 1.
c  
c  Modified 7/3/08, suppressing one-point convection zone (e.g. in connection
c  with local maximum at edge of convective core). 
c  So far implemented only for imixc6 = 1.
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
      implicit double precision (a-h,o-z)
      logical nosd,notd,dtest,skipt,noder,time0,conv,convos,lastmd,
     *  norct,norche
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
c
      parameter(idr1mx = nspcmx+3, idermx = ((nspcmx+3)*(nspcmx+4))/2,
     *  ivrmx1 = ivarmx+1, nspcm2=2*nspcmx, ivrmx4=4*ivarmx,
     *  iyfdmx = 1+ivarmx*(3+2*ivarmx))
c
      dimension y(*),zk(*),f(*),fd(ifd,*),alam(ialam,*),
     .  alamd(ialam,ifd,*),h(*),hd(ihd,*),ap(*),aq(*),
     .  zz(*),dzdy(iz,*),
     .  ddac(6),isig(4),dlmb(3,3),xhe3(4),dwf(35),comp(nspcmx),
     .  yp(ivarmx),zzp(ivarmx),fp(ivarmx),dtst(ivarmx),alamp(1,ivarmx),
     .  cvcpar(10), cvvar(icvvar), cvvarp(icvvar),alamrh(4,krnrmx),
     .  xche3(4), ftprt(nspcmx,nnmax), fthppp(nnmax), fthcnp(nnmax),
     .  y_orig(ivarmx)
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev, syear
      common/ln10/ amm
      common/rhcn/ a1,a2,a3,a4,zdummy,nvar,ifwrt,irhtst,inentr,
     *  ii1,ii2,ii3,icomp
      common/clshft/ alshft
      common/heavy/ zatmos, zhc, zh(1)
      common/comgrt/ isprot, irotcn, a5, riner0, riner, omgrt0,
     *  omgrot(1)
      common/totmss/ am, rs
      common/crhsvr/ qx, radius
      common/bcatms/ tmnbc,tmxbc,sbcfct,flsatm,ntaubc,iopatm,icsrad
      common/cqhopf/ d2qhpf,tgrfct,dtgrfc(5)
      common/cmtime/ age, time0, lastmd
      common/step/ dt
      common/thetac/ theta(ivarmx)
      common/noiter/ iter1, ntime, epspr, eam
      common/cmpder/ ft(idr1mx, nspcmx)
      common/prtvar/ rhl,ak,akt,akp,akx,eps(idr1mx),dr,dac,conv,convos
      common/caddvr/ addvar(5,ngmax)
      common/cnvout/ ddrad,rrcnv,aacnv,ddacad,amach, 
     *  drr(5), ddr(5), da(5),pturb,ycnv,dyda
      common/rnratd/ rnal(10,1)
      common/rnrout/rnzt(10),rnuw(10)
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,jnc,frcf(6),frcl(6),
     *  xrcf(6), xrcl(6), conv_dum(27)
      common/convpf/ nff(6),nlf(6),incf,inccrf,frcff(6),frclf(6),ifxcon
      common/convvr/ cvvarf(icvvar,6), cvvarl(icvvar,6),
     *  cvvars(icvvar,nnmax)
      common/cnvovs/ icnvos, jcnvos, clcovs, cldovs, alphos,
     *   rczl, rczlfx, rcnvos, qlcnos, rhobos, dmrovs, fmrovs,
     *   roscb1, roscb2, ddoscb, aoscb, boscb, coscb, yftlmx, dtlymx
      common/covsec/ icsove, icsovc, alpove, alpovc, 
     *  qmxove, qmxovc, rmxove, rmxovc, imxove, imxovc, 
     *  nmxove, nmxovc
      common/cvcpsc/ cvcpst(20)
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor,nmxcp,qmxcp,rmxcp,
     *  qmxmin,qmxmax,qmxscn,nmxscn
      common/crxcor/ cqc, cqcp, crxmn(nspcmx), crxmnp(nspcmx),
     *  crxxh_sc(nnmax)
      common/cmxcnt/ ddrmix, inmixc, imixc0, imixc1, imixc2, imixc3,
     *  imixc4, imixc5, imixc6
      common/crxstr/ rxstr(nspcm2,1)
      common/cdiffu/ idiffus, itbdif, rhob, rmxdif, tlfdif, dtldif,
     *  ctdfmx, rctdif, rctdf1, idfxbc, ismdif, nsmdif, idiffc1, idiffc2
      common/engfts/ fthpps, fthcns
      common/rnratd/ al(10,krnrmx),norct
      common/rnrhed/ alhe(10,10), norche
      common/he3eql/ xhe3eq,ieqhe3,iequsd
      common/he3fdg/ agesh,ifdhe3,iche30,iche31
      common/compsz/ xzerh, yzer, xzer3, xrz12, xrz13, xrz14, xrz16
      common/compos/ xseth, yset, xset3, xset12, xset13, xset14, xset16
      common/compvr/ cvr(icvrmx, 1)
      common/cstcvr/ icvh1, icvhe4, icvhe3, icvc12, icvc13, icvn14, 
     *  icvo16
      common/cstyst/ iyh1, iyhe3, iyc12, iyc13, iyn14, iyo16, 
     *  iyche4, iycc12
      common/cstazt/ iah1, iahe3, iac12, iac13, ian14, iao16, 
     *  iache4, iacc12
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec
      common/engcnt/ xhzlm1, xhzlm2, nxhzer, icnocs, nspec,
     *  iheccs, nspect
      common/engche/ xmxrhe
      common/enggrv/ epsg(nnmax)
      common/cenghc/ epshec(idermx), fthec(idermx,2)
      common/cepcno/ epscnc, eppcno(4)
      common/degfct/ thte, iscren
      common/eqstcl/ iradp,mode,dtest,skipt,noder
      common/eqprcl/ dmeqpr(12),idiag,iscskp
      common/eqstd/ xi(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     .  dlt(4),gm1,tprh,trhp,rhxp
      common/eqsout/ east(30),xii(30),dne(10),dph(20),he(20),pe(10),
     *  hi(20),pii(10),hh(20),ph(10),hr(20),pr(10),pcoul(10),hcoul(10)
      common/opcdat/ akk(4)
      common/nmbmsh/ nn
      common/cyfstr/ yzfstr(ivrmx4,nnmax)
      common/yprtst/ yprt(ivrmx1,nnmax)
      common/opcxdr/ akxa
      common/cyfdst/ iyfdst,ialamw,yffdst(iyfdmx,nnmax),
     *  yffdsp(iyfdmx,nnmax)
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3,idgmxc,
     *  idgsdt,idgn75,idgn76,idgn77
c
      common/cdgrhb/ kdgrhb
      common/cdgphs/ kdgeos, kdgopc, kdgeng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  dataset designators
c
      common/cdsdsg/ idstrl,idshvz,idszab,idsevl,idssum,idscen,idsnov,
     *  idsoov,idsgng,idsgsm,idstm1,idstm2
c
      save
c
      data isig/2,3,5,6/
      data flpr,tlpr,xhpr,epsskt /-100.,-100.,-100.,1./
      data initrs /0/
c
      data init /0/
c
      ntime_tst=-100
c
c  in present version, initialize flags for physics diagnostics, to
c  allow for return. 
c  Flag initialized to 1: in physics routines, return rather than stop,
c  resetting flag to -1
c
      kdgeos=1
      kdgopc=1
      kdgeng=1
c
c  initialize flag for error in routine
c 
      kdgrhb=1
c
c  initialize surface radius (needs checking)
c
      if(initrs.eq.0) then
        rs=6.9599d10
        initrs=1
      end if
c
c  set internal flag for neglecting entropy term
c
      inenti=inentr
c
c  set storage index for luminosity
c
      ilum=4
c
c  initialize turbulent pressure to zero
c
      pturb=0
      ptpg=0
c
c  initialize flag for convection-zone diagnostics
c
      idgcon=0
c
      qx=10.d0**x
      radius=1.d11*10.d0**y(1)
      amms=amsun*am
      amass=amms*qx
c
c  store original dependent variables
c
      do i=1,ivarmx
	y_orig(i)=y(i)
      end do
c
      if(init.eq.0) then
	xmxrhe=-1.d10
	rczl=1.d20
	init=1
      end if
c
c  when not including convective overshoot, set rcnvos to
c  impossible value
c
c  Note that cqc is used for convective-core mass when setting
c  energy generation in shrinking core, allowing for undercorrection
c
      if(icnvos.ne.1.and.icnvos.ne.2) rcnvos = -1.d10
c
      if(n.eq.1) then
	rs=1.d11*10.d0**y(1)
	hhpcz=1.
        cqc_tfrp=cqc_tfr
      end if
c
c..      if(lastmd.and.n.ge.598) then
c..        idgrhs=1
c..        idgeng=2
c..        idgopc=1
c..        idgeos=1
c..      else
c..        idgrhs=0
c..        idgeng=0
c..        idgopc=0
c..        idgeos=0
c..      end if
c
      if(idgeng.ge.2.and.istdpr.gt.0) then
        write(istdpr,*) 'Entering rhs'
        call dmprcn
      end if
      if(idgeng.eq.-3.and.istdpr.gt.0) write(istdpr,*) 
     *  'On entry to rhs n, thte =',n,thte
c
c  for iter = 0, reset nxhzer to determine cut-off for
c  old and new hydrogen abundance
c
      if((ntime.eq.0.or.iter.eq.0).and.n.eq.1) then
	nxhzer=nn+1
	nxhzrp=nn+1
      end if
c
c  set trivial part of time0 transformation
c
   10 do 11 i=1,4
      if(i.ne.2) then
        zz(i)=y(i)
        dzdy(i,i)=1.d0
      end if
   11 continue
c
c  calculate log(rho),log(kappa) and log(eps)
c
      fl=y(2)
      tl=y(3)
      rlred=y(1)
      allred=y(4)
      xh=y(5)
c
      if(n.eq.1) allrds=allred
c
c  set helium abundance, depending on iheccs. Note that it is assumed
c  that 4He burning is only relevant in regions where hydrogen is exhausted
c
      if(xh.ge.2.e-10.or.iheccs.eq.0) then
        yh=1-xh-zh(n)
	zhh=zh(n)
      else
	yh=y(iyche4)
	zhh=1.d0-xh-yh
	zh(n)=zhh
	if(xh.ge.1.d-9.and.zhh.ge.0.1.and.istdpr.gt.0) write(istdpr,*) 
     *    'Z-error in rhs at n =',n,' zhh =',zhh
      end if
c
      nosd=.false.
      notd=.true.
c
c  test for skip trial
c
      skipt=n.gt.1.and.abs(fl-flpr).lt.epsskt
     *  .and.abs(tl-tlpr).lt.epsskt
     *  .and.abs(xh-xhpr).lt.epsskt
      flpr=fl
      tlpr=tl
      xhpr=xh
c
      idiag=0
      call eqstf(fl,tl,xh,yh,zhh,nosd,notd)
      if(kdgeos.lt.0) then
	kdgrhb=-1
	return
      end if
      if(idgeng.eq.-3.and.istdpr.gt.0) write(istdpr,*) 
     *  'After eqstf n, thte =',n,thte
      if(idgeng.ge.2.and.istdpr.gt.0) then
        write(istdpr,*) 'after call of eqstf'
        call dmprcn
      end if
      pl=log10(pt(1))
c
c  set remaining part of transformation
c
      zz(2)=pl
      dzdy(2,2)=pt(2)
      dzdy(2,3)=pt(3)
c
c  store log p and (d y2/d X)_p,T in addvar (addvar(2,.) is used elsewhere)
c
      addvar(1,n)=pl
      addvar(3,n)=-pt(4)/pt(2)
c
      if(.not.time0) then
c
c  Note that setting of X into zz may be reset below,
c  depending on whether or not this is in region
c  where X should be set to zero, or in convection zone,
c  where X must be decoupled from other variables.
c
c  If 4He burning is considered set derivative of pt wrt Y instead
c  where X is exhausted
c
        zz(5)=xh
        dzdy(5,5)=1.d0
	if(iheccs.eq.0.or.xh.ge.1.e-5) then
          dzdy(2,5)=pt(4)
        else
          dzdy(2,iyche4)=-pt(4)
        end if
c
        if(nvar.ge.6) then
	  do 13 i=6,nvar
          zz(i)=y(i)
  13      dzdy(i,i)=1
        end if
      end if
c
c  store values at previous timestep for testing and computation
c  of gravitational energy release
c
      if(iter.eq.0) then
	do 15 i=1,nvar
   15   yprt(i,n)=y(i)
      end if
c
      p=pt(1)
      t=10.d0**tl
      r=1.d11*1.d1**rlred
      rhl=log10(rho(1))
      call opact(rhl,tl,xh,zhh,ak,rkr,rkt,rkx)
      if(kdgopc.lt.0) then
	kdgrhb=-1
	return
      end if
      if(idgeng.ge.2.and.istdpr.gt.0) then
        write(istdpr,*) 'after call of opact'
        call dmprcn
      end if
c
c  ** extra diagnostics
c
      if(idgrhs.eq.1.and.istdpr.gt.0) write(istdpr,1120) 
     *  n,x,(y(i),i=1,6),rhl,tl,xh,ak
      akf=rkr*rho(2)
      akt=rkt+rkr*rho(3)
      akx=akxa+rkr*rho(4)
c
c  energy generation
c
      ieqheo=ieqhe3
      ift=idr1mx
      comp(1)=y(5)
c
c  set remaining composition variables
c
      if(ifdhe3.eq.1) then
        nosd=.false.
        call he3abd(fl,tl,xh,yh,zhh,agesh,xhe3,anu,nosd)
        if(idgeng.ge.2.and.istdpr.gt.0) then
          write(istdpr,*) 'after call of he3abd'
          call dmprcn
        end if
        if(iter.gt.0) y(6)=xhe3(1)
        comp(2)=xhe3(1)
c..        if(n.eq.1900) write(6,*) '#D# rhs.00 comp(2) =',comp(2)
        notd=.true.
	ixx3 = 2
      else if(ieqhe3.ne.1) then
	comp(2) = y(6)
c..        if(n.eq.1900) write(6,*) '#D# rhs.01 comp(2) =',comp(2)
	ixx3 = 2
      else
	ixx3 = 1
      end if
c
c  At time0, get CNO abundances from common/compvr/, depending on icnocs
c
      if(time0) then
	if(icnocs.eq.1) then
	  comp(ixx3+1) = cvr(icvn14,n)
	else if(icnocs.eq.2) then
	  comp(ixx3+1) = cvr(icvn14,n)
	else if(icnocs.eq.3) then
	  comp(ixx3+1) = cvr(icvc13,n)
	  comp(ixx3+2) = cvr(icvn14,n)
	else if(icnocs.eq.4) then
	  comp(ixx3+1) = cvr(icvc12,n)
	  comp(ixx3+2) = cvr(icvc13,n)
	  comp(ixx3+3) = cvr(icvn14,n)
	end if
c
c  also, possibly set 4He burning abundances
c
	if(iheccs.ne.0) then
	  comp(ixx3+ispcno+1)=yh
	  comp(ixx3+ispcno+2)=cvr(icvc12,n)
	end if
      else 
	if(icnocs.ge.1) call store(y(5+ispxx3),comp(ixx3+1),ispcno)
	if(iheccs.ne.0) then
          call store(y(5+ispxx3+ispcno),comp(ixx3+ispcno+1),2)
	  comp(ixx3+ispcno+1)=yh
        end if
      end if
c
      notd=.true.
      idgeno = idgeng
      if(n.ge.443.and.n.le.446) idgeng=2
c
      xhzlm1_orig=xhzlm1
      if(ntime.eq.ntime_tst.and.nxhzrp.lt.nn.and.n.lt.nxhzrp) 
     *  xhzlm1=-1.d10
      if(ntime.eq.ntime_tst.and.n.ge.460.and.n.le.500.and.istdpr.gt.0)
     *  then
	idgeng=-11
c..	write(istdpr,*) 'n, x, =',n, x
c..	write(istdpr,*) 
c..     *    'fl, tl, comp, yh, zhh =',fl, tl, comp(1), comp(2), yh, zhh
c..	iftprt=1
      else
	idgeng=idgeno
        iftprt=0
      end if
      if(istdpr.le.0) iftprt=0
c..      if(ntime.eq.0.and.n.ge.550) idgeng=2
      call engenr(fl,tl,comp,yh,zhh,eps,ft,ift,notd)
c
      xhzlm1=xhzlm1_orig
c..      if(iter1.eq.0) then 
c..        xr=1.d0
c..        write(82,'(4f11.6,1p11e13.6)') xr, tl, fl, 
c..     *    (comp(i),i=1,5),eps(1),epscnc
c..        write(84,'(3f11.6,1p11e13.6)') xr, tl, fl, 
c..     *    rnuw(1),(rnal(1,i),i=1,10)
c..      end if
      if(kdgeng.lt.0) then
	kdgrhb=-1
	return
      end if
c
c  store rates of change of composition at the previous time step
c
      if(iter.le.0.and..not.time0) then
	do i=1,nspcmx
	  ftprt(i,n)=ft(1,i)
	  fthppp(n)=fthpps
	  fthcnp(n)=fthcns
        end do
      end if
c
c  for testing, set analytical hydrogen abundance in selected range.
c
      if(iter.gt.0.and.ntime.eq.ntime_tst.and.n.lt.nxhzrp) then
	fthppm=(fthpps+fthppp(n))/2.d0
	fthcnm=(fthcns+fthcnp(n))/2.d0
	alphaft=fthppm/fthcnm
	xhtil=yprt(5,n)*exp(fthcnm*dt)/(1+alphaft*yprt(5,n))
	xhnew=xhtil/(1.d0-alphaft*xhtil)
c..	if(n.ge.400.and.n.le.500) write(6,*) 
c..     *    '#D# n, x, fthppm, fthcnm, xhtil, xhnew, y_orig(5)',
c..     *    n, x, fthppm, fthcnm, xhtil, xhnew, y_orig(5)
      end if
c
c  test for resetting hydrogen-burning rate near hydrogen exhaustion
c
      if(ntime.eq.ntime_tst.and.n.lt.nxhzrp.and.iter1.gt.1
c..     *  .and.yprt(5,n).le.0.4) then
     *  .and.y(5).le.xhzlm1) then
	ft1_orig=ft(1,1)
	ft(1,1)=(xhnew-yprt(5,n))/(dt*theta(5))
     *          -(1.d0-theta(5))*ftprt(1,n)/theta(5)
	do i=2,idr1mx
	  ft(i,1)=0.d0
        end do
	y(5)=xhnew
	zz(5)=xhnew
	if(istdpr.gt.0) write(istdpr,
     *    '(/'' At n = '',i4,'' x ='',f10.5,
     *    '' ft(1,1) reset from'',1pe13.5,'' to'',e13.5)')
     *    n, x, ft1_orig, ft(1,1)
c
c  switch off for now. Seems very problematic
c
c..      else if(y(5).lt.xhzlm1.and.iter.gt.0.and.n.le.nxhzrp.and.
c..     *    iter1.gt.0) then
c..	ft1_orig=ft(1,1)
c..	ft(1,1)=(1.d-10-yprt(5,n))/(dt*theta(5))
c..     *          -(1.d0-theta(5))*ftprt(1,n)/theta(5)
c..	do i=2,idr1mx
c..	  ft(i,1)=0.d0
c..        end do
c..	y(5)=1.d-10
c..	if(istdpr.gt.0) write(istdpr,
c..     *    '(/'' At n = '',i4,'' x ='',f10.5,
c..     *    '' ft(1,1) reset from'',1pe13.5,'' to'',e13.5)')
c..     *    n, x, ft1_orig, ft(1,1)
      end if
c
c  test for reset of maximum x at which 4He reactions take place
c
      if(iheccs.ne.0.and..not.norche) then
	if(x.gt.xmxrhe.and.istdpr.gt.0) write(istdpr,*) 
     *    'xmxrhe reset to',x
	xmxrhe=max(x,xmxrhe)
      end if
c
      idgeng = idgeno
      if(idgeng.ge.2.and.istdpr.gt.0) then
        write(istdpr,*) 'after call of engenr'
        call dmprcn
      end if
      if(idgrhs.eq.-1.and.istdpr.gt.0) write(istdpr,1100) 
     *  ' after engenr', (zz(i),i=1,5)
c
      if(idgrhs.eq.1.and.istdpr.gt.0) write(istdpr,1130) (eps(i),i=1,5)
c
c  reset derivatives when using He3 fudge
c
      if(ifdhe3.eq.1.and.xhe3(1).gt.1.e-10) then
        do 18 i=2,4
        eps(i)=eps(i)+eps(5)*xhe3(i)
   18   ft(i,1)=ft(i,1)+ft(5,1)*xhe3(i)
      end if
c
c  test for using constant-coefficient solution for He3 abundance
c
      xhe3or=y(6)
      if(ifdhe3.eq.0.and.iche31.eq.1.and.iter.gt.0) then
        do 19 k=1,4
        alamrh(1,k)=rho(1)*al(1,k)
	do 19 j=2,4
   19   alamrh(j,k)=rho(j)+al(j,k)
	nosd=.false.
	call he3abc(xh,yh,zhh,dt,yprt(6,n),alamrh,xche3,xhe3eq,
     *    anu,4,nosd)
c
	ft(1,2)=(xche3(1)-yprt(6,n))/dt
	do 20 i=2,4
   20   ft(i,2)=(xche3(i)-yprt(6,n))/dt
	ft(5,2)=0
c
      end if
c
      iequsd=ieqhe3
c
      ieqhe3=ieqheo
c
c  store composition variables in common/compvr/ if there have been
c  reactions
c
      cvr(1,n)=xh
      cvr(2,n)=yh
      if(.not.norct) then
        if(ieqhe3.eq.1) then
          cvr(3,n)=xhe3eq
        else
	  cvr(3,n)=comp(2)
        end if
c
        call setcvr(n)
      end if
c
c  possibly store abundances related to 4He burning
c
      if(iheccs.ne.0) then
	if(.not.norche) then
	  cvr(2,n)=y(iyche4)
	  cvr(icvc12,n)=y(iycc12)
	else
c..	  y(iyche4)=yh
c..	  y(iycc12)=cvr(icvc12,n)
        end if
      end if
c
c  test for setting nxhzer, based on values at previous time step
c  for iter = 0. Assume that rx varies linearly with X,
c  so that time dependence is exponential.
c
c  Also, set nxhzrp, for previous time step.
c
      if(iter.eq.0) then
	if(ntime.eq.ntime_tst.and.istdpr.gt.0) write(istdpr,*)
     *    '#D# n, xhzlm1 etc.',
     *    n, xhzlm1, y(5), exp(dt*ft(1,1)/y(5))
	if(y(5)*exp(dt*ft(1,1)/y(5)).le.xhzlm1.and.nxhzer.gt.nn) then
	  nxhzer=n
	  if(istdpr.gt.0) write(istdpr,*) 'nxhzer set to ', nxhzer, 
     *      '  in rhs'
        end if
	if(y(5).le.xhzlm1.and.nxhzrp.gt.nn) then
	  nxhzrp=n
	  if(istdpr.gt.0) write(istdpr,*) 'nxhzrp set to ', nxhzrp, 
     *      '  in rhs'
        end if
      end if
c
c  test for resetting energy generation quantities to zero, where
c  hydrogen abundance is too low
c  Also reset values stored in zz and dzdy.
c
c  Note that zz may be reset later
c
      if(n.ge.nxhzer) then
	zz(5)=1.d-10
        dzdy(2,5)=0.d0
	if(iheccs.eq.0) then
	  do 21 i=1,idr1mx
	  eps(i)=0
	  do 21 k=1,nspcmx
   21     ft(i,k)=0
	else
c
c  test for setting He abundance
c
	  if(x.gt.xmxrhe) then
	    zz(iyche4)=1.d0-zh(n)-zz(5)
	    nspecr=nspect
          else
	    nspecr=nspec
          end if
	  do 22 i=4,nspecr+3
   22     eps(i)=0
	  do 23 i=1,nspecr+3
	  do 23 k=1,nspecr
   23     ft(i,k)=0
        end if
      end if
c
c  test for setting correction factor in temperature gradient for smooth
c  transition from atmospheric T(tau) relation
c
      if(icsrad.ne.1.or.tl.gt.4.6) then
        tgrfct=1
	call zero(dtgrfc,5)
      else
c
        flxrat=7.125771d-15*10.d0**(4.d0*tl+2.d0*rlred-allred)
        tautil=1.3333333d0*(flxrat-1) +tmxbc
        call ttau(tautil,0.d0,0.d0,0.d0,xh,zhh,thopf,qhopf,dqhopf,
     *    dtrs,dtls,dttau)
        tgrfct=1+dqhopf
c
c  set derivatives
c
        d2qfct=1.3333333d0*amm*d2qhpf
	dtgrfc(1)=2*d2qfct
	dtgrfc(2)=0
	dtgrfc(3)=4*d2qfct
	dtgrfc(4)= -d2qfct
	dtgrfc(5)=0
      end if
c
c     *************************************************
c
c  note: to avoid under- and overflow on univac, equations
c  were changed on 18/1/84 to involve log r - 11.
c  this corresponds to change in a1, a2 and a3 in s/r mnevol.
c
c  on 4/1/85 equations were changed to be written in terms of
c  log r - 11 and log l - 33, with corresponding changes to the
c  definition of a1 - a4.
c
c  set the right hand sides
c
      f(1)=a1*10.d0**(x-3.d0*rlred-rhl)
      if(isprot.eq.0) then
        f(2)=-a2*10.d0**(2.d0*x-4.d0*rlred-pl)
      else
c
c  include centrifugal acceleration. Note that a5 contains 
c  factor 2/3
c
        f(2)=-a2*10.d0**(2.d0*x-4.d0*rlred-pl)
     *       +a5*omgrot(n)*omgrot(n)*10.d0**(x-rlred-pl)
      end if
c
      f(3)=-a3*tgrfct*10.d0**(ak+x-4.d0*(tl+rlred))*
     *        ((10.d0**allred)-alshft)
      f(4)=a4*eps(1)*10.d0**(x-allred)
c..      if(ntime.eq.0.and.n.ge.550) write(6,*) 
c..     *  '#D# n, x, fl, tl, comp(1), a4, eps(1), allred',
c..     *  n, x, fl, tl, comp(1), a4, eps(1), allred
      if(idgeng.ge.2.and.istdpr.gt.0) then
        write(istdpr,*) 'after setting f(1) - f(4)'
        call dmprcn
      end if
c
      if(time0.and.iter.gt.0) then
	epsg(n)=0
      end if
c
c  ** extra diagnostics
c
      if(idgrhs.eq.1.and.istdpr.gt.0) write(istdpr,1150) 
     *  a3,ak,x,allred,tl,rlred,f(3)
c
c  test for convection
c
      dr=f(3)/f(2)
      dac=dr
      ddrad=dr-dad(1)
      addvar(5,n)=ddrad
c
c  pressure scale height
c
      hhp=-r*f(1)/f(2)
c
c  do not allow convection in atmosphere (when called from prtsol)
c
      conv=dr.gt.dad(1).and.x.le.0
      ddacad=dr-dad(1)
      ddarad=dr-dad(1)
c
c  skip setting of convective boundaries when iter = 0
c  or in atmosphere (when called from prtsol)
c
      if(iter.eq.0.or.x.gt.0) go to 28
c
c #NOTE#: The following looks most peculiar and is probably best
c         dropped.     27/7/06
c#dr#c
c#dr#c  also skip setting of convective boundaries (except for the 
c#dr#c  first iteration) inside what has otherwise been identified as
c#dr#c  a convective core
c#dr#c
c#dr#      if(iter1.gt.1.and.nmxcor.gt.0.and.
c#dr#     *  (qx.le.qmxcor.or.n.ge.nmxcor)) then
c#dr#c
c#dr#c  test for setting position of core in convection-zone parameters
c#dr#c  Flag for setting proper position of start of convective core
c#dr#c  after call of tnrkt
c#dr#c
c#dr#	if(nl(inc).ne.nn) then
c#dr#	  inc=inc+1
c#dr#	  nf(inc)=-1
c#dr#	  write(istdpr,*) '#D-1# n, inc, nf(inc) =',n,inc,nf(inc)
c#dr#	  nl(inc)=nn
c#dr#	  frcl(inc)=0
c#dr#	end if
c#dr#	go to 28
c#dr#      end if
c
      dcnp=dcn
      call store(cvvar,cvvarp,icvvar)
      amach=0
      if(n.lt.0) go to 28
c
c  initialize flag for convection-zone boundary (1 for start, 2 for end)
c
      iccase=0
c
c  find edges of convection zone
c
c  nnf=nf(inc) is the first meshpoint (counted from the surface), and
c  nnl=nl(inc) is the last inside or on the edge of the inc-th
c  convection zone. further if frc=frcf(inc) (or frcl(inc)) then
c  for any function a=a(n) frc*a(nnf)+(1-frc)*a(nnf-1)
c  (or frc*a(nnl)+(1-frc)*a(nnl+1)) is the value of a, from linear
c  interpolation, at the edge of the convection zone.
c  Variables (m/M, r/R, p, T, rho, kappa) at edges of
c  convection zones are stored in cvvarf(.,inc) and cvvarl(.,inc)
c
      dcn=-ddacad
      cvvar(1)=qx
      cvvar(2)=r/rs
      cvvar(3)=pt(1)
      cvvar(4)=rho(1)
      cvvar(5)=t
      cvvar(6)=10.d0**ak
      cvvar(7)=hhp
      cvvar(8)=dad(1)
      cvvar(9)=dr
      cvvar(10)=dr
      call store(cvvar,cvvars(1,n),icvvar)
c
c  include test on x to avoid resetting in atmosphere
c
      if(n.eq.1.and.x.le.0) then
c  
c  surface point
c  
c  extra variables for output near boundary of convective core
c
	if(inc.gt.0.and.nl(inc).eq.nn.and.idgcon.gt.0) then
	  ncp1=nf(inc)-3
	  ncp2=nf(inc)+3
        else
	  ncp1=0
	  ncp2=0
        end if
c
	if(idgcon.gt.0) then
	  nep1=nl(1)-5
	  nep2=nl(1)+5
        else
	  nep1=0
	  nep2=0
        end if
c
        inc=0
	dcnp=0
        call zero(cvvarp,7)
        inc=0
        dcs=dcn
        n1=n
        frc=1
	if(dcs.le.0) then
c
c  start of convection zone
c
          inc=inc+1
          iccase=1
          nf(inc)=n1
	  if(istdpr.gt.0) write(istdpr,*) 
     *      '1 Set nf n, n1, inc, nf(inc) =',n,n1,inc,nf(inc)
          frcf(inc)=frc
          call store(cvvar,cvvarf(1,inc),icvvar)
	  if(istdpr.gt.0) write(istdpr,*) 'inc, cvvarf(.,inc):',
     *      inc,(cvvarf(i,inc),i=1,icvvar)
        else
	  go to 28
	end if
      end if
c
c  end of test for surface point
c
c  extra output
c
      if(n.ge.ncp1.and.n.le.ncp2.and.istdpr.gt.0) then
	if(n.eq.ncp1) write(istdpr,*) 
     *    ' n, conv, log q, log T, log rho, X, log kappa, ddad'
	write(istdpr,*) n, conv
	write(istdpr,'(i5,6f11.6)') n, x, tl, rhl, xh, ak, ddacad
      end if
      if(n.ge.nep1.and.n.le.nep2.and.istdpr.gt.0) then
	if(n.eq.nep1) write(istdpr,*) 
     *    ' n, conv, log q, log T, log rho, X, log kappa, ddad'
	write(istdpr,*) n, conv
	write(istdpr,'(i5,6f12.7)') n, x, tl, rhl, xh, ak, ddarad
      end if
c  
c  test for change of sign of dcn  
c  
      n1=n 
      if(dcnp.ne.0.and.dcn*dcnp.le.0) then
        frc=dcnp/(dcnp-dcn)
	if(dcn.lt.0) then
c
c  start of convection zone
c
          inc=max(inc+1,1)
          iccase=1
	  if(istdpr.gt.0) then
	    write(istdpr,*) 'change of sign. increment inc to',inc
            write(istdpr,'(3i5,1p10e13.5)')
     *        ntime,iter1,n,10.d0**x,rhl, tl, xh, ak, ddrad,allred, 
     *      rlred,tgrfct
	  end if
          inc=min0(6,inc)
          nf(inc)=n1
	  if(istdpr.gt.0) write(istdpr,*) 
     *      '2 Set nf n, n1, inc, nf(inc) =',n,n1,inc,nf(inc)
          frcf(inc)=frc
c
          do 24 i=1,icvvar
   24     cvvarf(i,inc)=frc*cvvar(i)+(1-frc)*cvvarp(i)
	  if(istdpr.gt.0) write(istdpr,*) 'inc, cvvarf(.,inc):',
     *      inc,(cvvarf(i,inc),i=1,icvvar)
c
c  store parameters for output of quantities relevant to convective
c  core
c
          call store(cvcpar,cvcpst,6)
c
          cvcpst(7)=10.d0**x
          cvcpst(8)=rhl
          cvcpst(9)=tl
          cvcpst(10)=xh
          cvcpst(11)=ak
          cvcpst(12)=ddacad
	  do 25 i=1,6
   25     cvcpst(12+i)=frc*cvcpst(6+i)+(1-frc)*cvcpst(i)
        else if(dcn.gt.0) then
c 
c  end of convection zone
c
	  nl(inc)=n1-1
          iccase=2
	  if(istdpr.gt.0) write(istdpr,*) 
     *      '3 Set nf n, n1, inc, nf(inc) =',n,n1,inc,nf(inc)
          frcl(inc)=1-frc
c
          do 26 i=1,icvvar
   26     cvvarl(i,inc)=frc*cvvar(i)+(1-frc)*cvvarp(i)
	  if(istdpr.gt.0) write(istdpr,*) 'inc, cvvarl(.,inc):',
     *      inc,(cvvarl(i,inc),i=1,icvvar)
	  if(inc.eq.1) then
	    rczl = frcl(inc)*rp+(1-frcl(inc))*r
	    hhpcz = hhp
          end if
	else
c
c  dcn = 0. test on dcnp to decide case
c
          if(dcnp.gt.0) then
c
c  start of convection zone
c
            inc=inc+1
            inc=min0(6,inc)
            nf(inc)=n1
	    if(istdpr.gt.0) write(istdpr,*) 
     *        '4 Set nf n, n1, inc, nf(inc) =',n,n1,inc,nf(inc)
            frcf(inc)=1
            call store(cvvar,cvvarf(1,inc),icvvar)
	    if(istdpr.gt.0) write(istdpr,*) 'inc, cvvarf(.,inc):',inc,
     *        (cvvarf(i,inc),i=1,icvvar)
c
c  store parameters for output of quantities relevant to convective
c  core
c
            call store(cvcpar,cvcpst,5)
c
            cvcpst(7)=10.d0**x
            cvcpst(8)=rhl
            cvcpst(9)=tl
            cvcpst(10)=xh
            cvcpst(11)=ak
            cvcpst(12)=ddarad
          else 
c 
c  end of convection zone
c
	    n1=n+1
	    nl(inc)=n1-1
            iccase=2
	    if(istdpr.gt.0) write(istdpr,*) 
     *        '1 Set nl. n, n1, inc, nl(inc) =',n,n1,inc,nl(inc)
            frcl(inc)=0
            call store(cvvar,cvvarl(1,inc),icvvar)
	    if(istdpr.gt.0) write(istdpr,*) 'inc, cvvarl(.,inc):',inc,
     *        (cvvarl(i,inc),i=1,icvvar)
	    if(inc.eq.1) then
	      rczl = r
	      hhpcz = hhp
            end if
	  end if
        end if
c
c  test for end of 1-point convection zone. If imixc6 = 1, suppress it
c
        if(iccase.eq.2.and.imixc6.eq.1.and.nl(inc).eq.nf(inc)) then
          if(istdpr.gt.0) write(istdpr,'(/'' nf(inc) = nl(inc) = '',i5,
     *      '' Suppress one-point convection zone'')') nf(inc)
          inc=inc-1
        end if
      end if
c  
c  test for convection to centre 
c  
      if(n.eq.nn.and.dcn.lt.0) then  
        frcl(inc)=1
        nl(inc)=nn 
        call store(cvvar,cvvarl(1,inc),icvvar)
      end if
c
c  in convection zone or overshoot region reset f(3)
c
   28 convos = .false.
c
c  for idiffc2 .gt. 0, suppress resetting of temperature gradient
c  in all but the first convection zone
c
      if(conv.and.idiffc2.gt.0.and.inc.gt.1.and.n.ge.nf(2)) then
	if(istdpr.gt.0) 
     *    write(istdpr,*) '#D# Suppress convection at n, x =',n, x
	conv=.false.
      end if
c
      rp=r
c
c  test for new (5/1/00) formulation of overshoot
c
      if(.not.conv) then
c..        if(imxove.gt.0.and.inc.eq.imxove.and.r.ge.rmxove) then
        if(imxove.gt.0.and.r.lt.rs*cvvarf(2,imxove).and.r.ge.rmxove) 
     *    then
	  icnovs = icsove
c..        else if(imxovc.gt.0.and.inc.eq.imxovc-1.and.r.le.rmxovc) then
        else if(imxovc.gt.0.and.r.le.rmxovc) then
	  icnovs = icsovc
        else
	  icnovs = 0
        end if
      end if
c
      if(conv) then
c  store opacity in array for mixlng
        akk(1)=10.d0**ak
        akk(2)=akf
        akk(3)=akt
        akk(4)=akx
        call mixlng(t,hhp,dr,dac,ddac,ptpg)
c
        f(3)=dac*f(2)
        if(idgrhs.eq.1.and.istdpr.gt.0) write(istdpr,1170) f(2),f(3),dac
c
c  Test for overshoot using either old MJM formulation or
c  new overshoot, fitting Matthias Rempel simulations
c  (note: different usage of clcovs and cldovs in this cas.
c  Note also that a fairly arbitrary cut-off at 3*width has been
c  imposed)
c
      else if((icnvos.eq.1.and.dad(1)-dr.le.cldovs.and.
     *    rczl - r.lt.alphos*hhpcz).or.(icnvos.eq.2.and.jcnvos.eq.1
     *  .and.r.lt.1.001*rczl.and.r.ge.rczl-rs*(cldovs+3.*clcovs))) then
	convos=.true.
c  store opacity in array for oversh
        akk(1)=10.d0**ak
        akk(2)=akf
        akk(3)=akt
        akk(4)=akx
	if(icnvos.eq.1) then
	  xovs=0.d0
        else
	  xovs=((r-rczl)/rs+cldovs)/clcovs
	end if
	if(istdpr.gt.0) write(istdpr,*) 'OVS',
     *    n, conv, r/rs, rczl/rs, dad(1), dr, hhp, xovs
        call oversh(t,hhp,dr,dac,ddac,xovs)
c
        f(3)=dac*f(2)
        rcnvos = r
	if(istdpr.gt.0) write(istdpr,*) 
     *    ' Now rcnvos/R, ddacad =',rcnvos/rs, ddacad
        qlcnos = x
c
c  test for new overshoot, with adiabatic temperature gradient
c
      else if(icnovs.eq.1) then
	f(3)=dad(1)*f(2)
	ddacad = 0
	dac = dad(1)
c
      end if
c
c  test for possible subadiabatic MR overshoot formulation
c  Note that option icnvos = 202 was introduced to force
c  usage from prtsol, regardless of the (somewhat opaque)
c  use of iter and iter1
c
      if(icnvos.eq.202.or.
     *  (icnvos.eq.2.and.jcnvos.ge.2.and.iter1.gt.1)) then
	if(n.eq.1.and.icnvos.ne.202) call iovsmr(iter1)
	drovs=(r-rczlfx)/rs+cldovs
	if((drovs.le.6.*dmrovs.or.r.le.rczlfx*1.001).and.
     *    drovs.ge.-clcovs.and.dac.le.dad(1)+1.d-4) then
	  convos=.true.
	  xovs=((r-rczlfx)/rs+cldovs)/dmrovs
	  if(istdpr.gt.0) write(istdpr,*) 'OVS',
     *      n, conv, r/rs, rczlfx/rs, dad(1), dr, hhp, xovs
	  call oversh(t,hhp,dr,dac,ddac,xovs)
	  if(istdpr.gt.0) write(istdpr,*) 'dac =',dac
	  f(3)=dac*f(2)
        end if
      end if
c
c  test for using simple cubic subadiabatic region
c
      if(icnvos.eq.203.or.(icnvos.eq.3.and.iter1.gt.1)) then
	if(n.eq.1.and.icnvos.ne.203) call iovscb(iter1)
	if(n.eq.1.and.icnvos.eq.203.and.istdpr.gt.0) 
     *    write(istdpr,*) '#D# roscb1, roscb2, ddoscb, alpha, beta',
     *    roscb1, roscb2, ddoscb, alpha, beta, aoscb, boscb, coscb
        xovs=(r/rs-roscb1)/(roscb2-roscb1)
	if(xovs.ge.0.and.xovs.le.1) then
	  call oversh(t,hhp,dr,dac,ddac,xovs)
	  if(istdpr.gt.0) write(istdpr,*) 'x, xovs, dac =',
     *      r/rs, xovs, dac
	  f(3)=dac*f(2)
        end if
      end if
c
      cvvar(10)=dac
      cvvars(10,n)=dac
c
c
c  store parameters for later analysis of convective core boundary
c
      cvcpar(1)=10.d0**x
      cvcpar(2)=rhl
      cvcpar(3)=tl
      cvcpar(4)=xh
      cvcpar(5)=ak
      cvcpar(6)=ddacad
c
      if(idgrhs.eq.-1.and.istdpr.gt.0) write(istdpr,1100)
     *  'after testing for convection', (zz(i),i=1,5)
      if(idgeng.ge.2.and.istdpr.gt.0) then
        write(istdpr,*) 'after testing for convection'
        call dmprcn
      end if
c
c  store first variables in yzfstr
c
      do i=1,nvar
	yzfstr(i,n)=y(i)
      end do
      do i=1,4
	yzfstr(nvar+i,n)=zz(i)
	yzfstr(2*nvar+i,n)=f(i)
      end do
c
      if(time0) go to 58
c
c      ---------------------------------------------------
c
c  set quantities for time derivatives
c  Also store in common crxstr for possible later use in s/r cmpcvc
c  Store values at current time step in 1 - icomp, values
c  at previous time step in icomp+1 - 2*icomp
c
      if(iter.eq.0) then
	iftst=icomp
      else
	iftst=0
      end if
      do 32 k=5,nvar
      k1=k-4
      f(k)=ft(1,k1)
   32 rxstr(iftst+k1,n)=ft(1,k1)
c..      if(n.eq.520) then
c..	write(istdpr,*) '#D# n, iter, iftst, nvar, rxstr'
c..	write(istdpr,*) n, iter, iftst, nvar, (rxstr(i,n),i=1,nspcm2)
c..      end if
c
c  test for resetting in mixed core
c  This depends on imixc1, as above, with imixc1 = 0 corresponding to
c  the old, erroneous treatment for iter = 0.
c
      if(iter.eq.0) then
c
c  previous time step
c
        if(((imixc1.eq.0.and.nmxcor.gt.0).or.
     *      (imixc1.gt.0.and.nmxcp.gt.0)).and.qx.le.cqcp) then
          do 34 k=5,nvar
   34     f(k)=crxmnp(k-4)
        end if
      else
c
c  current time step (note that this is independent of
c  imixc1, so far)
c
c  For imixc6 = 1 (may become default later) set fictitious value such
c  that the correct result is obtained with proper centralization
c
c  This was further modified 13/2/08, since previous reset seemed
c  wrong!
c
c
        if(ntime.eq.12.and.qx.le.0.005.and.istdpr.gt.0) write(istdpr,*) 
     *    '#D# n, qx, qmxcor, etc.',n,qx,qmxcor,nmxcor,qmxcp,cqc,cqcp        
        if(inmixc.ne.1.and.(nmxcor.gt.0.or.nmxcp.gt.0)) then
          if(qx.le.qmxcor) then
            if(imixc6.ne.1) then
              do k=5,nvar
                f(k)=crxmn(k-4)
              end do
            else
	      if(qx.le.qmxcp) then
		rxxc=0
                do k=5,nvar
c..                  f(k)=crxmn(k-4)
                  rxxc=(y(k)-yprt(k,n))/dt
		  rxxp=crxmnp(k-4)
                  f(k)=(rxxc-(1.d0-theta(k))*rxxp)/theta(k)
		end do
              else
c
c  set actual local rate of change in composition, in intermediate region
c
		do k=5,nvar
		  k1=k-4
                  rxxc=(y(k)-yprt(k,n))/dt
		  rxxp=rxstr(icomp+k1,n)
                  f(k)=(rxxc-(1.d0-theta(k))*rxxp)/theta(k)
                end do
		if(ntime.eq.20.and.istdpr.gt.0) write(istdpr,*)
     *            '#D# n, qx, qmxcp, qmxcor, rxxc,',
     *            ' rxstr(icomp+k1,n), f(5)',
     *            n, qx, qmxcp, qmxcor, rxxc,rxstr(icomp+k1,n), f(5)
              end if
            end if
c
c  for imixc5 = 1, reset to take into account possible semiconvective region
c
	    if(imixc5.eq.1) f(5)=crxxh_sc(n)
c..          else if(qx.le.cqcp) then
          else if(qx.lt.qmxcp) then
c
c  test for applying undercorrection to cqc
c
            if(iter1.gt.10.and.imixc5.eq.1) then
              if(iter1.le.20) then
                ucyxc=0.7d0
              else if(iter1.le.30) then
                ucyxc=0.4d0
              else 
                ucyxc=0.d0
              end if
	      cqc_tfr=ucyxc*cqc+(1.d0-ucyxc)*cqc_tfrp
            else
              cqc_tfr=cqc
            end if
            if(imixc6.ne.1) then
              tfrct=(cqcp-qx)/(cqcp-cqc_tfr)
            else
              tfrct=(qmxcp-qx)/(qmxcp-qmxcor)
              if(ntime.eq.12.and.istdpr.gt.0) 
     *          write(istdpr,*) '#D# n, qx, tfrct =',
     *          n, qx, tfrct
            end if
c
c  reset ficticious f so that the combined effect, with proper
c  centralization, matches resetting in cmpcvc.
c  For now (22/1/08) apply only when imixc6 = 1
c  Further revised 8/2/08. 
c
            do 37 k=5,nvar
            if(iche31.ne.1.or.k.ne.6) then
	      if(imixc6.ne.1.or.theta(k).eq.0) then
                f(k)=(1-tfrct)*f(k)+tfrct*crxmn(k-4)
              else
                k1=k-4
                tfrct1=1.d0-tfrct
                theta1=1.d0-theta(k)
                f(k)= tfrct*tfrct*crxmn(k1)+
     *            tfrct1*((tfrct1*theta(k)+tfrct)*f(k)+
     *                tfrct1*theta1*rxstr(k1+icomp,n)+
     *                (theta(k)*tfrct-theta1)*crxmnp(k1))/theta(k)
	      end if
	    end if
   37       continue
          end if
        end if
      end if
c
      if(idgrhs.eq.-1.and.istdpr.gt.0) write(istdpr,1100) 
     *  'after setting f(5 - nvar)', (zz(i),i=1,5)
c
      if(irhtst.eq.1.and.iter.eq.0) then
	yprt(7,n)=f(5)
      end if
      bb=a4*10.d0**(x-allred)
      prhm=amm*pt(1)/rho(1)
c
c  entropy terms in energy equation. set to zero if inenti =1
c
      if(inenti.eq.1) then
        do 44 i=1,nvar
   44   alam(1,i)=0
      else
c
        do 46 ia=1,3
        ia1=ia+1
        i=isig(ia)
   46   alam(1,i)=bb*(prhm*pt(ia1)-ht(ia1))
        alam(1,1)=0
        alam(1,4)=0
        if(nvar.ge.6) then
	  do 48 i=6,nvar
   48     alam(1,i)=0
        end if
      end if
c
c  set epsilong
c
c  Note: for the time being assume that backwards differences are
c  used. Hence do not centralize the lambda-s properly.
c  This should be fixed up at some stage.
c
      if(iter.gt.0.and.iter1.gt.0) then
	epsg(n)=0
	do 49 i=1,nvar
   49   epsg(n)=epsg(n)+alam(1,i)*(y(i)-yprt(i,n))/(bb*dt)
c..        if(mod(n-1,100).eq.0) write(87,'(3i5,1p50e13.5)') ntime,iter1,n,
c..     *    10.**x,(yprt(i,n),i=1,nvar),(y(i),i=1,nvar),
c..     *    (alam(1,i),i=1,nvar),epsg(n)
      end if
c
c  test for unreasonably large epsilong in small convective core
c
      epsav=1.d33*((10.d0**allred)-alshft)/(qx*amms)
      if(abs(epsg(n)).gt.0.02*epsav) then
        if(((qx.le.1.0*cqc.or.qx.le.1.0*cqcp).and.xh.ge.0.01)) then
c..     *    .or.epsg(n).lt.-500.*epsav) then
c
	  if(istdpr.gt.0) write(istdpr,'(a,i5,1p3e11.3)') 
     *      ' Excessive epsg',n, qx, epsg(n),epsav
	  inenti=1
	  epsg(n)=0
	  do 50 i=1,nvar
  50      alam(1,i)=0
        end if
      end if
c
      if(idgrhs.eq.-1.and.istdpr.gt.0) write(istdpr,1100) 
     *  'after setting alam', (zz(i),i=1,5)
c
c  set remaining quantities into yzfstr
c
      do i=5,nvar
	yzfstr(nvar+i,n)=zz(i)
	yzfstr(2*nvar+i,n)=f(i)
      end do
      do i=1,nvar
	yzfstr(3*nvar+i,n)=alam(1,i)
      end do
c
c      ---------------------------------------------------
c
c  the derivatives
c
   58 continue
      if(idgeng.ge.2.and.istdpr.gt.0) then
        write(istdpr,*) 'before setting derivatives '
        call dmprcn
      end if
      if(idgrhs.eq.1.and.istdpr.gt.0) then
	write(istdpr,1250) (f(i),i=1,6)
        write(istdpr,1251) (zz(i),i=1,6)
      end if
      if(idgrhs.eq.-1.and.istdpr.gt.0) write(istdpr,1100) 
     *  'before return', (zz(i),i=1,5)
c
      if(iter.le.0.and..not.time0.and.idgn77.gt.0) then
        yffdsp(1,n)=x
        call store(y,yffdsp(2,n),nvar)
        call store(zz,yffdsp(2+nvar,n),nvar)
        call store(f,yffdsp(2+2*nvar,n),nvar)
	if(iftprt.eq.1) write(istdpr,*) 'f(5) =',f(5)
        ii=1+3*nvar
        do j=1,nvar
          do i=1,nvar
            ii=ii+1
            yffdsp(ii,n)=fd(i,j)
          end do
        end do
	if(idgn77.ge.2) then
	  call store(alam,yffdsp(2+nvar*(3+nvar),n),nvar)
          ii=1+nvar*(4+nvar)
          do j=1,nvar
            do i=1,nvar
              ii=ii+1
              yffdsp(ii,n)=alamd(1,i,j)
            end do
          end do
	end if
	if(idgn77.ge.3) then
	  yffdsp(ii+1,n)=10.d0**ak
	  yffdsp(ii+2,n)=rho(1)
	  yffdsp(ii+3,n)=cp(1)
	  yffdsp(ii+4,n)=dad(1)
	  yffdsp(ii+5,n)=dr
c
c  the following two variables could be used for X-derivatives.
c  See diffusive version of routine
c
	  yffdsp(ii+6,n)=0.d0     
	  yffdsp(ii+7,n)=0.d0     
c
c  store nabla_R - nabla_ad as set in mixcor in the last column
c
	  yffdsp(ii+8,n)=addvar(4,n)
	  ii=ii+8
	end if
	iyfdst=ii
        if(n.eq.nn) then
	  iter0=0
	  ialamw=1
	  write(77) iter0,ialamw,iyfdst,nn
          write(77) ((yffdsp(i,k),i=1,iyfdst),k=1,nn)
          if(istdpr.gt.0) write(istdpr,*) 
     *      'RHS quantities written to d/s 77 in rhs, ',
     *      'nn, nvar, iyfdst =', nn, nvar, iyfdst
        end if
      end if
      if(iter.le.0) return
c
c  ***********************************************************************
c
      f1=amm*f(1)
      fd(1,1)=-3.d0*f1
      fd(1,2)=-rho(2)*f1
      fd(1,3)=-rho(3)*f1
      f2=amm*f(2)
      if(isprot.eq.0) then
        fd(2,1)=-4.d0*f2
      else
	fd(2,1)=4.d0*amm*a2*10.d0**(2.d0*x-4.d0*rlred-pl)
     *              -amm*a5*omgrot(n)*omgrot(n)*10.d0**(x-rlred-pl)
      end if
      fd(2,2)=-pt(2)*f2
      fd(2,3)=-pt(3)*f2
      f3=amm*f(3)
      alfact=1.d0/(1.d0-alshft/(10.d0**allred))
      if(conv.or.convos) then
c
c  convection zone or (old) overshoot region
c
        do 60 k=1,3
   60   fd(3,k)=dac*fd(2,k)+ddac(k)*f(2)
	fd(3,4)=     alfact*ddac(4)*f(2)
c
      else if (icnovs.eq.1) then
c
c  new formulation for overshoot
c
	fd(3,1)=dad(1)*fd(2,1)
	fd(3,2)=dad(1)*fd(2,2)+dad(2)*f(2)
	fd(3,3)=dad(1)*fd(2,3)+dad(3)*f(2)
c
      else
c
c  radiative region
c
        fd(3,1)=(dtgrfc(1)/tgrfct-4.d0)*f3
        fd(3,2)=(dtgrfc(2)/tgrfct+akf)*f3
        fd(3,3)=(dtgrfc(3)/tgrfct+akt-4.d0)*f3
        fd(3,4)=alfact*(dtgrfc(4)/tgrfct+1)*f3
c
      end if
c
      f4=amm*f(4)
      fd(4,2)=eps(2)*f4
      fd(4,3)=eps(3)*f4
      fd(4,4)=-f4
c
c  end for time0
c
      if(time0) then
	if(idgn77.gt.0.and.iter.ge.0) then
          yffdst(1,n)=x
	  call store(y,yffdst(2,n),nvar)
	  call store(zz,yffdst(2+nvar,n),nvar)
	  call store(f,yffdst(2+2*nvar,n),nvar)
	  if(iftprt.eq.1) write(istdpr,*) 'f(5) =',f(5)
	  ii=1+3*nvar
          do j=1,nvar
	  do i=1,nvar
	    ii=ii+1
	    yffdst(ii,n)=fd(i,j)
          end do
          end do
	  if(n.eq.nn) then
	    iyfdst=1+nvar*(3+nvar)
	    ialamw=1
	    write(77) iter1,ialamw,iyfdst,nn
            write(77) ((yffdst(i,k),i=1,iyfdst),k=1,nn)
	  end if
	end if
        return
      end if
c
c  number of independent variables for energy generation
c
      iamx=nvar-2
c
c  number of "thermodynamic" variables to take derivatives 
c  with respect to, depending on whether or not X is forced to zero
c
      if(n.lt.nxhzer) then
	ibmx=3
      else
	ibmx=2
      end if
c
c  only set X-derivatives when X is not forced to zero,
c  and outside convection zones
c
      if(conv.or.(inmixc.ne.1.and.nmxcor.gt.0.and.qx.le.qmxcor)) 
     *  then
c
c  reset dzdy(2,5)
c
	dzdy(2,5)=0
	if(iheccs.ne.0) dzdy(2,iyche4) = 0
c
      else 
	if(n.lt.nxhzer) then
          fd(1,5)=-rho(4)*f1
          fd(2,5)=-pt(4)*f2
          if(conv) then
            fd(3,5)=ddac(5)*f(2)+dac*fd(2,5)
          else if(icnovs.eq.1) then
	    fd(3,5)=dad(1)*fd(2,5)+dad(4)*f(2)
          else
            fd(3,5)=akxa*f3
          end if
c
c  switch off thermodynamic and opacity derivatives wrt Y, on the
c  assumption that Y + Z is constant, and so therefore is the electron
c  number per unit mass
c
c..	else if(iheccs.ne.0) then
c..          fd(1,iyche4)=rho(4)*f1
c..          fd(2,iyche4)=pt(4)*f2
c..          if(conv) then
c..            fd(3,iyche4)=-ddac(5)*f(2)+dac*fd(2,iyche4)
c..          else if(icnovs.eq.1) then
c..	    fd(3,iyche4)=dad(1)*fd(2,iyche4)-dad(4)*f(2)
c..          else
c..            fd(3,iyche4)=-akx*f3
c..          end if
	end if
	do 75 i=5,nvar
   75   fd(4,i)=eps(i-1)*f4
c
        do 77 ia=1,iamx
        ia1=ia+1
	if(ia.le.2) then
	  i = ia+1
        else
	  i = ia+2
	end if
        do 77 k=5,nvar
   77   fd(k,i)=ft(ia1,k-4)
      end if
c
c  test for resetting in connection with mixed core
c
      if(inmixc.ne.1.and.(nmxcor.gt.0.or.nmxcp.gt.0)) then
	if(qx.le.qmxcor) then
	  do 78 i=1,nvar
	  do 78 k=5,nvar
	  fd(i,k)=0
   78     fd(k,i)=0
        else if(qx.le.cqcp) then
	  do 79 i=1,nvar
	  do 79 k=5,nvar
          if(iche31.ne.1.or.k.ne.6) then
	    if(imixc6.ne.1.or.theta(k).eq.0) then
              fd(k,i)=(1-tfrct)*fd(k,i)
            else
              fd(k,i)=tfrct1*(tfrct1*theta(k)+tfrct)*fd(k,i)/theta(k)
	    end if
	  end if
	  if(iftprt.eq.1.and.i.eq.6.and.istdpr.gt.0) 
     *      write(istdpr,*) 'Rescaled fd(k,6):',fd(k,6)
   79     continue
        end if
c
c  For imixc6 = 1 reset in region between old and new (growing)
c  limit of core, excluding initial models.
c
c  Switch off again; seems to be a bad idea (13/2/08)
c
c..        if(imixc6.eq.1.and.qx.gt.qmxcp.and.qx.lt.qmxcor.and.
c..     *    theta(5).gt.0.and.ntime.ge.3) then
c..	  fd(5,5)=1.d0/(dt*theta(5))
c..	  write(istdpr,*) '#D# Setting fd(5,5) to',fd(5,5)
c..        end if
      end if
c
c  try setting X-derivative in intermediate region with growing convective core
c  (implemented 23/4/09)
c
      if(imixc6.eq.1.and.nmxcor.gt.0.and.qx.ge.qmxcp.and.qx.le.qmxcor
     *  .and.y(5).ge.0.01) then
        fd(5,5)=f(5)/y(5)
        if(istdpr.gt.0) write(istdpr,*) 
     *    '#D1 n, qx, qmxcp, qmxcor, f(5),y(5),fd(5,5) =',
     *    n, qx, qmxcp, qmxcor, f(5),y(5),fd(5,5)
      end if
c
c  derivatives of alam (set to zero when entropy term is ignored)
c
      if(inenti.ne.1) then
c
        jj=4
        do 84 ia=1,3
        do 84 ib=ia,3
        jj=jj+1
        ddlb=prhm*pt(jj)-ht(jj)
        if(ib.gt.ia) dlmb(ib,ia)=ddlb
   84   dlmb(ia,ib)=ddlb
        prhm=amm*prhm
c
        do 85 ia=1,iamx
        ia1=ia+1
        if(ia.lt.4) then
          i=isig(ia)
          pta=pt(ia1)
          alamd(1,i,4)=-amm*alam(1,i)
          do 83 ib=1,ibmx
          ib1=ib+1
   83     alamd(1,i,isig(ib))=bb*(dlmb(ia,ib)+prhm*(pt(ib1)-rho(ib1))
     .      *pta)
	end if
   85   continue
c
      else
c
        do 87 i=1,nvar
        do 87 j=1,nvar
   87   alamd(1,i,j)=0
c
      end if
c
      if(idgn77.gt.0.and.iter.ge.0) then
        yffdst(1,n)=x
        call store(y,yffdst(2,n),nvar)
        call store(zz,yffdst(2+nvar,n),nvar)
        call store(f,yffdst(2+2*nvar,n),nvar)
	if(iftprt.eq.1) write(istdpr,*) 'f(5) =',f(5)
        ii=1+3*nvar
        do j=1,nvar
          do i=1,nvar
            ii=ii+1
            yffdst(ii,n)=fd(i,j)
          end do
        end do
	if(idgn77.ge.2) then
	  call store(alam,yffdst(2+nvar*(3+nvar),n),nvar)
          ii=1+nvar*(4+nvar)
          do j=1,nvar
            do i=1,nvar
              ii=ii+1
              yffdst(ii,n)=alamd(1,i,j)
            end do
          end do
	end if
	if(idgn77.ge.3) then
	  yffdst(ii+1,n)=10.d0**ak
	  yffdst(ii+2,n)=rho(1)
	  yffdst(ii+3,n)=cp(1)
	  yffdst(ii+4,n)=dad(1)
	  yffdst(ii+5,n)=dr
c
c  the following two variables could be used for X-derivatives.
c  See diffusive version of routine
c
	  yffdst(ii+6,n)=0.d0     
	  yffdst(ii+7,n)=0.d0     
c
c  store nabla_R - nabla_ad as set in mixcor in the last column
c
	  yffdst(ii+8,n)=addvar(4,n)
	  ii=ii+8
	end if
	iyfdst=ii
        if(n.eq.nn) then
	  ialamw=1
	  write(77) iter1,ialamw, iyfdst,nn
          write(77) ((yffdst(i,k),i=1,iyfdst),k=1,nn)
          if(istdpr.gt.0) write(istdpr,*) 
     *      'RHS quantities written to d/s 77 in rhs, ',
     *      'nn, nvar, iyfdst =', nn, nvar, iyfdst
        end if
      end if
c
      return
 1100 format(' In rhs, ',a/(1p5e13.5))
 1120 format(/' s/r rhs at n =',i5,'  x =',1pe15.7/
     *  ' y:',6e15.7/' rhl, tl, xh, ak =',4e15.7)
 1130 format(' eps:',1p5e13.5)
 1150 format(' a3, ak, x, allred, tl, rlred, f(3) =',1p7e14.6)
 1170 format(' f(2),f(3),dac =',1p3e15.7)
 1250 format(' f:',1p6e12.4)
 1251 format(' z:',1p6e12.4)
 1300 format(///
     *  ' test rhs. ntime, iter, n, x, f(1-nvar), dtst(1-nvar):'/)
 1350 format(3i4,1p13e13.5)
 1400 format(' n, x, y(1-6) =',i5,1p7e15.7/'       z(1-6) =',
     *      20x,6e15.7/'       f(1-6) =',20x,6e15.7)
 1500 format(' fd:')
 1501 format(' dzdy:')
 1510 format(1p6e12.4)
      end
      subroutine oversh(t,hhp,dr,dac,ddac,xovs)
c
c  crude model for overshoot, setting gradient as specified
c  function of nabla_ad - nabla_rad
c
c  Original version: 11/4/91
c
c  Modified 3/11/94 by MJM to include jcnvos = 3 option for
c  smooth transition
c
c  Modified 4/8/04, to include fit to Matthias Rempel simulations
c  for icnvos = 2.
c
c  Modified 18/10/07, to include simple cubic subadiabatic layer
c  for icnovs = 3.
c
c  Modified 31/7/08 with generalized Rempel formulation depending on alphos
c
      implicit double precision (a-h, o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
      dimension ddac(*), is(5), dddrad(5)
      dimension xr(nnmax), dgrad(nnmax)
      common/eqstd/ xi(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     .   dlt(4),gm1,tprh,trhp,rhxp
      common/opcdat/ akk(4)
      common/cnvovs/ icnvos, jcnvos, clcovs, cldovs, alphos, 
     *   rczl, rczlfx, rcnvos, qlcnos, rhobos, dmrovs, fmrovs,
     *   roscb1, roscb2, ddoscb, aoscb, boscb, coscb, yftlmx, dtlymx
      common/cnvout/ ddrad,rrcnv,aacnv,ddacad,amach, 
     *  drr(5), ddr(5), da(5),pturb,ycnv,dyda
      common/convvr/ cvvarf(icvvar,6), cvvarl(icvvar,6),
     *  cvvars(icvvar,nnmax)
      common/ln10/ amm
c
      data is/2,3,5,1,4/
c
c  logarithmic derivatives of radiative gradient
c
      do 10 j=1,5
         ddr(j)=0
         dddrad(j)=0
   10    ddac(j)=0
c
      do 15 i=1,3
         j=is(i)
         ddac(j)=dad(i+1)
         dddrad(j)=-dad(i+1)
   15    ddr(j)=pt(i+1)+akk(i+1)
c
      ddr(3)=ddr(3)-4
      ddr(4)=1
c
      do 20 j=1,5
         ddr(j)=amm*dr*ddr(j)
c
c  finally set derivative of rad.grad. - ad. grad.
c
   20    dddrad(j)=ddr(j)+dddrad(j)
c
c  define ddrad=(\nabla_r-\nabla_a), needed in all options
c
      ddrad=dr-dad(1)
c
c
c  Test for MJM or MR formulations or cubic subadiabatic layer
c
      if(icnvos.eq.2.or.icnvos.eq.202) go to 50
      if(icnvos.eq.3.or.icnvos.eq.203) go to 60
c
c  -------------------------------------------------------------------------
c
c  MJM formulation
c
      if(jcnvos.le.1) then
c
c---> Option jcnvos=0 or 1
c
c  sharp overshoot
c
c  set gradients etc.
c
         denom=1+clcovs*(cldovs+ddrad)
         ddacad=ddrad/denom
c
         dac=ddacad+dad(1)
c
         derfct=(1+clcovs*cldovs)/(denom*denom)
         do 30 j=1,5
   30       ddac(j)=ddac(j)+derfct*dddrad(j)
c
      else if (jcnvos.ge.2) then
c
c---> Option jcnvos=2 or 3
c
c  fuzzy overshoot
c
	 yy=-ddrad/cldovs
	 ddacad=-cldovs*yy*yy*(2-yy)
	 dac=ddacad+dad(1)
c
c  derivatives
c
         do 35 j=1,5
   35       ddac(j)=ddac(j)+yy*(4-3*yy)*dddrad(j)
c
      else
c
c---> Option jcnvos=3 or higher
c
c  option for a two zones overshoot layer with
c  continuous first and second derivatives of \nabla
c
	 yy=1-clcovs*(ddrad+cldovs)/cldovs
c
         if (yy.le.0.0) then
c
c  adiabatic section of the overshoot layer
c  for: 0 < (\nabla_a-\nabla_r) < dldovs*(clcovs-1)/clcovs
c
            dac=dad(1)
c
c  the derivatives of \nabla are just the derivatives of \nabla_a
            do 40 j=1,5
 40            ddac(j)=ddac(j)
	 else
c
c  transition from adiabatic stratification to radiative stratification
c  for: dldovs*(clcovs-1)/clcovs < (\nabla_a-\nabla_r) < dldovs
c
            ff1=(1-yy)**3*(1-1/clcovs+(3-2/clcovs)*yy
     *           +3*(2-1/clcovs)*yy**2)
            dac=dr+cldovs*ff1
c
            ff2=(1-yy)**2*(1+2*yy-15*(2*clcovs-1)*yy**2)
            do 45 j=1,5
 45            ddac(j)=ddr(j)+dddrad(j)*ff2
         end if
      end if
c
      return
c
c  ------------------------------------------------------------------
c
c  MR formulation. Note that xovs must be set to (r-r_t)/d
c  in argument list.
c
c  For jcnvos = 2, note that iovsmr must have been called to set
c  MR coefficients.
c
   50 continue
c
c  test for version of MR formulation
c
      if(jcnvos.eq.1) then
c
        phi=0.5d0*(1.d0-tanh(xovs))
        ddacad=phi*ddrad
        dac=ddacad+dad(1)
        do j=1,5
          ddac(j)=phi*ddr(j)
        end do
      else if(jcnvos.eq.2.or.jcnvos.eq.3) then
c
c..        ddacad=-fmrovs*(1.d0-tanh(xovs))
        ddacad=-2.d0*fmrovs/(alphos+exp(2.d0*xovs))
        dac=ddacad+dad(1)
	do j=1,5
	  ddac(j)=0.d0
        end do
      end if
      return
c
c  -------------------------------------------------------------------
c
c  Simple cubic subdiabatic layer (with little physical realism!)
c
   60 ddacad=ddoscb+xovs*(aoscb+xovs*(boscb+coscb*xovs))
      dac=ddacad+dad(1)
      do j=1,5
        ddac(j)=0.d0
      end do
c
      end
      subroutine iovsmr(iter1)
c
c  Initialize coefficients for MR overshoot with possibly
c  smooth subadiabatic gradient
c
c  For larger iter1, do not reset dmrovs but only factor.
c
c  Original version: 12/8/04
c
c  Modified 31/7/08 to include generalized Rempel formulation
c
      implicit double precision(a-h, o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
      dimension xr(nnmax),dgrad(nnmax),cvvrft(icvvar)
      common/totmss/ am, rs
      common/convvr/ cvvarf(icvvar,6), cvvarl(icvvar,6),
     *  cvvars(icvvar,nnmax)
      common/cnvovs/ icnvos, jcnvos, clcovs, cldovs, alphos,
     *   rczl, rczlfx, rcnvos, qlcnos, rhobos, dmrovs, fmrovs,
     *   roscb1, roscb2, ddoscb, aoscb, boscb, coscb, yftlmx, dtlymx
      common/nmbmsh/ nn
      common/noiter/ iter, ntime, eps, eam
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      data alphos_orig, yftlmx, dtlymx /-1.d0, 0.2784d0, 0.4356d0/
c
      save
c
      itermx=100
      ucy=0.6
c
c  test for initializing yftlmx and alphos
c
      if(alphos.ne.alphos_orig) then
        if(jcnvos.eq.3) then
          yftlmx=yftil_max(alphos,dtlymx)
          if(istdpr.gt.0) write(istdpr,'(/
     *      '' For alphos ='',f10.5,'' yftil_max set to '',f10.5)')
     *      alphos, yftlmx
        else
          alphos=1.d0
          if(istdpr.gt.0) write(istdpr,'(/
     *      '' alphos initialized to 1.d0'')')
        end if
        alphos_orig=alphos
      end if
c
c  test for no resetting
c
      if(iter1.ge.itermx.or.(mod(iter1,4).ne.0.and.iter1.ge.4)) then
c..      if(iter1.ge.itermx) then
	if(istdpr.gt.0) then
	  write(istdpr,*) ' No reset in iovsmr'
          write(istdpr,*)' iter1, dmrovs, fmrovs, rczlfx/rs', 
     *      iter1, dmrovs, fmrovs, rczlfx/rs
	end if
	return
      end if
c
      rczlfp=rczlfx
      fmrovp=fmrovs
      rczlfx=rczl
      if(iter1.ge.4) rczlfx=ucy*rczl+(1.d0-ucy)*rczlfp
      rt=rczlfx/rs-cldovs
      rf=rt-clcovs
      if(istdpr.gt.0) then
        write(istdpr,*) 'In iovsmr, set rt, rf =',rt, rf
        write(istdpr,*) 'rczlfx, rs, clcovs, cldovs =',
     *    rczlfx, rs, clcovs, cldovs
      end if
c 
      do n=1,nn
	xr(n)=cvvars(2,n)
	if(xr(n).ge.rf) ifit=n
      end do
c
c  set linear interpolation factors
c
      fct2=(rf-xr(ifit))/(xr(ifit+1)-xr(ifit))
      fct1=1.d0-fct2
      if(istdpr.gt.0) write(istdpr,*) 
     *  'xr(ifit), xr(ifit+1), fct1, fct2 =',
     *  xr(ifit), xr(ifit+1), fct1, fct2
c
      do i=1,icvvar
	cvvrft(i)=fct1*cvvars(i,ifit)+fct2*cvvars(i,ifit+1)
      end do
c
      if(istdpr.gt.0) write(istdpr,*) 'cvvrft =',cvvrft
c
c  test for resetting only factor
c
      if(iter1.ge.itermx-4) then
c..	fmrovs=(cvvrft(8)-cvvrft(9))/
c..     *    (1.d0-tanh((rf-rt)/dmrovs))
	fmrovs=0.5d0*(cvvrft(8)-cvvrft(9))*
     *    (alphos+exp(2.d0*(rf-rt)/dmrovs))
	if(iter1.ge.4) fmrovs=ucy*fmrovs+(1.d0-ucy)*fmrovp
        if(istdpr.gt.0) then
	  write(istdpr,*) 
     *      'Reset only fmrovs. cvvrft(8), cvvrft(9), fmrovs =',
     *      cvvrft(8), cvvrft(9), fmrovs
          write(istdpr,*)' iter1, dmrovs, fmrovs, rczlfx/rs', 
     *      iter1, dmrovs, fmrovs, rczlfx/rs
	end if
	return
      end if
c
      call derive(xr,cvvars(9,1),dgrad,nn,icvvar,1,1,1)
c
      dgrdft=fct1*dgrad(ifit)+fct2*dgrad(ifit+1)
      dgrdf1=(cvvars(9,ifit+1)-cvvars(9,ifit))/(xr(ifit+1)-xr(ifit))
      if(istdpr.gt.0) write(istdpr,*) 'dgrdft, dgrdf1 =',
     *  dgrdft,dgrdf1
      dgrdft=dgrdf1
c
      yf=dgrdft/(cvvrft(8)-cvvrft(9))
c
      if(istdpr.gt.0) write(istdpr,*) 'ifit, dgrdft, etc.',
     *  ifit, dgrdft,cvvrft(8),cvvrft(9),yf
      np1=max(ifit-10,1)
      np2=min(ifit+10,nn)
      if(istdpr.gt.0) write(istdpr,107) 
     *  (n,xr(n),cvvars(9,n),dgrad(n),n=np1,np2)
c
      yftil=(rt-rf)*yf
      if(istdpr.gt.0) write(istdpr,*) ' yf, yftil =', yf, yftil
      if(yf.le.0) then
	write(istdou,105) yf, ifit, dgrdft, cvvrft(8),
     *    cvvrft(9)
	if(istdpr.ne.istdou.and.istdpr.gt.0)
     *    write(istdpr,105) ifit, dgrdft, cvvrft(8),
     *      cvvrft(9)
	if(dgrdft.le.0.and.istdpr.gt.0) then
	  np1=max(ifit-10,1)
	  np2=min(ifit+10,nn)
	  write(istdpr,107) (n,xr(n),cvvars(9,n),dgrad(n),n=np1,np2)
        end if
	stop 'Stop 0 in iovsmr'
      end if
c
c  test for existence of solution
c  If overall iteration is far from converged, attempt to replace by
c  location of maximimum value of yftilde (added 12/7/06)
c
c  Modified 31/7/08 for variable yftlmx. Limit set to be consistent
c  with original case.
c
      if(yftil.ge.yftlmx-0.0034) then
	if(eam.ge.100*eps) then
	  if(istdpr.gt.0) write(istdpr,110) yftil, rt, rf, dgrdft
	  yftil=yftlmx
	  dtil=dtlymx
	else
	  write(istdou,120) yftil, rt, rf, dgrdft
	  if(istdpr.ne.istdou.and.istdpr.gt.0)
     *      write(istdpr,120) yftil, rt, rf, dgrdft
	  stop 'Stop 1 in iovsmr'
        end if
c
      else
c
        call solve_dtil(yftil,dtil,alphos)
c
      end if
c  set final parameters
c
   40 dmrovs=dtil/yf
      fmrovs=alphos*(cvvrft(8)-cvvrft(9))/(2.d0-dmrovs*yf)
      if(iter1.ge.4) fmrovs=ucy*fmrovs+(1.d0-ucy)*fmrovp
      if(istdpr.gt.0) write(istdpr,*)
     *  ' iter1, dmrovs, fmrovs, rczlfx/rs', 
     *  iter1, dmrovs, fmrovs, rczlfx/rs
c
      return
  105 format(/' ***** Error in iovsmr. Negative yf =',1pe13.5/
     *        '       ifit, dgrad, nabla_ad, nabla_rad =',i5,3e13.5)
  107 format(/' n, x, nabla_rad, derivative:'/(i5,0pf10.5,1p2e13.5))
  110 format(/' ***** Warning in MR overshoot formulation. yftil =',
     *        1pe13.5/
     *        '       rt, rf, dgrad(ifit) =',1p3e13.5/
     *        '       Replace by yftilde_max, dtil_max')
  120 format(/' ***** Error in MR overshoot formulation. yftil =',
     *        1pe13.5/
     *        '       rt, rf, dgrad(ifit) =',1p3e13.5)
      end
      subroutine solve_dtil(yftil,dtil,alphos)
c
c  Solves equation for parameter dtil in MR subadiabatic overshoot
c  formulation
c
c  Original version: 12/7/06
c
c  Modified 31/7/08 for generalized Rempel formulation
c
      implicit double precision (a-h, o-z)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  initially solve equation by backsubstitution 
c
      dtil=0.2
      nit=0
c
   10 dtilo=dtil
      xtil=1.d0-dtilo
      dtil=2.d0*yftil/log((1.d0+xtil)/(alphos*(1.d0-xtil)))
      if(istdpr.gt.0) write(istdpr,*) 'dtil =',dtil
      if(abs(dtil-dtilo).lt.1.d-3) then
	go to 30
      else if(nit.ge.50) then
	write(istdou,120) nit, yftil, rt, rf, dgrdft
	if(istdpr.ne.istdou.and.istdpr.gt.0)
     *    write(istdpr,120) nit, yftil, rt, rf, dgrdft
	stop 'Stop 2 in iovsmr'
      else
	nit=nit+1
	go to 10
      end if
c
c continuing Newton iteration
c
   30 ff=dtil*log(2.d0-dtil)-dtil*log(alphos*dtil)-2.d0*yftil
c..      dff=log(2.d0-dtil)+dtil/(2.d0-dtil)-log(dtil)-1.d0
c  Sign of second term changed 31/7/08
      dff=log(2.d0-dtil)-dtil/(2.d0-dtil)-log(alphos*dtil)-1.d0
      dldtil=-ff/dff
      dtil=dtil+dldtil
      if(istdpr.gt.0) write(istdpr,*) 
     *  'dtil, ff, dff, dldtil = ',dtil, ff, dff, dldtil
      if(abs(dldtil).lt.1.d-8) then
	go to 40
      else if(nit.ge.50) then
	write(istdou,120) nit, yftil, rt, rf, dgrdft
	if(istdpr.ne.istdou.and.istdpr.gt.0)
     *    write(istdpr,120) nit, yftil, rt, rf, dgrdft
	stop 'Stop 2 in iovsmr'
      else
	nit=nit+1
	go to 30
      end if

c
   40 continue
c
      return
  120 format(/' ***** No convergence in MR overshoot formulation. ',
     *        ' nit =',i3/
     *        '       yftil, rt, rf, dgrad(ifit) =',1p4e13.5)
      end
      double precision function yftil_max(alpha,dtlymx)
c
c  finds maximum value of yftil for generalized Rempel formulation
c
c  Original version: 31/7/08
c
      implicit double precision(a-h, o-z)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      xx=0.5
      nit=0
c
   10 ff=log((2.d0-xx)/(alpha*xx)) - xx/(2.d0-xx) - 1.d0
      dff=(xx-4.d0)/(2.d0-xx)**2 - 1.d0/xx
      dxx=-ff/dff
      xx=xx+dxx
      nit=nit+1
      if(abs(dxx).lt.1.d-8) then
        go to 20
      else if(nit.ge.20) then
        write(istdou,'(/'' Convergence problems in yftil_max''/
     *    '' alpha, xx, ff, dff, dxx, nit ='',1p5e11.3,i4)')
     *    alpha, xx, ff, dff, dxx, nit
        if(istdpr.ne.istdou.and.istdpr.gt.0)
     *    write(istdpr,'(/'' Convergence problems in yftil_max''/
     *    '' alpha, xx, ff, dff, dxx, nit ='',1p5e11.3,i4)')
     *    alpha, xx, ff, dff, dxx, nit
        yftil_max=0.d0
        dtlymx=0.d0
        return
      end if
      go to 10
c
   20 yftil_max=0.5d0*xx*log((2.d0-xx)/(alpha*xx))
      dtlymx=xx
      return
      end
      subroutine iovscb(iter1)
c
c  Initialize coefficients for cubic subadiabatic gradient
c  smooth subadiabatic gradient
c
c  Original version: 18/10/07
c
      implicit double precision(a-h, o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
      dimension xr(nnmax),dgrad(nnmax),cvvrft(icvvar)
      common/totmss/ am, rs
      common/convvr/ cvvarf(icvvar,6), cvvarl(icvvar,6),
     *  cvvars(icvvar,nnmax)
      common/cnvovs/ icnvos, jcnvos, clcovs, cldovs, alphos,
     *   rczl, rczlfx, rcnvos, qlcnos, rhobos, dmrovs, fmrovs,
     *   roscb1, roscb2, ddoscb, aoscb, boscb, coscb, yftlmx, dtlymx
      common/nmbmsh/ nn
      common/noiter/ iter, ntime, eps, eam
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      roscb1=rczl/rs-clcovs
c 
      do n=1,nn
	xr(n)=cvvars(2,n)
	if(xr(n).ge.roscb1) ifit=n
      end do
c
      roscb1=xr(ifit)
      beta=cvvars(8,ifit)-cvvars(9,ifit)
c
      call derive(xr,cvvars(9,1),dgrad,nn,icvvar,1,1,1)
c
      drmax=3.d0*beta/dgrad(ifit)
      roscb2=rczl/rs+cldovs
      if(roscb2-roscb1.gt.drmax.or.cldovs.le.0) then
        roscb2=roscb1+drmax
	if(istdpr.gt.0) write(istdpr,
     *    '(/'' roscb2 set to maximum value ='',f10.5/)') roscb2
      end if
      alpha=(roscb2-roscb1)*dgrad(ifit)
      aoscb=alpha
      boscb=3.d0*beta-2.d0*alpha
      coscb=alpha-2.d0*beta
      ddoscb=-beta
      if(istdpr.gt.0) write(istdpr,*)
     *  '#D# roscb1, roscb2, ddoscb, alpha, beta',
     *  roscb1, roscb2, ddoscb, alpha, beta, aoscb, boscb, coscb
      return
      end
