      subroutine rhsd(x,yin,zk,zz,dzdy,ap,aq,f,fd,alam,alamd,h,hd,iz,
     .  ifd,ihd,ialam,n,iter)
c
c  right hand side subroutine for evolution.
c
c  Version including diffusion.
c
c  Note: iter = 0 when rhs is called from tnrkt at previous
c  time step, or from mshstr. A consequence of this is that
c  convection zone boundaries are not set.
c  This may cause problems if convective overshoot is
c  used, since this requires that the boundary of the
c  convective envelope is known.
c#ai#  This must be checked.
c
c  iter1 (in common/noiter/) should be set to the actual iteration number
c  in the routine calling tnrk. When called outside the iteration 
c  (for setting variables for printing, say), iter1 must be set to zero.
c
c  Note: for idiffus = 2, derivatives of right-hand side wrt y(5)
c  ( = Z) are not set. This might slightly affect convergence
c  but is unlikely to be significant (and we don't know the
c  derivatives of most quantities wrt to Z anyway!)
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
c  Flag for convection zones (to stabilize iteration, etc; 
c  added 28/10/03):
c
c    icvflg(n) =  0: radiative region
c    icvflg(n) =  1: convective region
c    icvflg(n) = -1: semiconvective region
c  
c  *************************************************************************
c
c  Note on diagnostic output to unit 77:
c  =====================================
c  
c  If idgn77 .ge. 1, results for the last timestep are kept on unit 77.
c  (Note that this creates a rather substantial overhead in execution
c  time and hence should be switched off for general runs.)
c  The results are stored in, and output with, the array yffdst, in the
c  following order:
c  
c  yffdst(1,n):                        x(n)
c  yffdst(2:nvar+1,n):                 y(i,n), i = 1:nvar
c  yffdst(nvar+2:2*nvar+1,n):          z(i,n), i = 1:nvar
c  yffdst(2*nvar+2:3*nvar+1,n):        f(i,n), i = 1:nvar
c  yffdst(3*nvar+2:nvar*(3+nvar)+1,n): fd(i,j,n), i, j = 1:nvar
c  
c  If idgn77 .ge. 2, in addition alam and alamd are stored, 
c  in the following order:
c  
c  yffdst(nvar*(3+nvar)+2:nvar*(3+ialamw+nvar)+1,n): alam(k,i,n),
c                   k = 1:ialamw, i = 1:nvar
c  yffdst(nvar*(3+ialamw+nvar)+2:nvar*(3+ialamw+(ialamw+1)*nvar)+1,n):
c     alamd(k,i,j,n), k=1:ialamw, i, j = 1:nvar
c
c  Here ialamw = 1 + idcomp
c
c  if idgn77 .ge. 3, in addition quantities for diffusion are stored,
c  in yffdst(nvar*(3+ialamw+(ialamw+1)*nvar)+1+2:...,n) (see below,
c  for now).
c  
c  *************************************************************************
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
c  Modified 16/8/92, to include diffusion and settling
c
c  Modified 20/8/92, to store variables at edges of convection zones
c  in common/convvr/
c
c  Modified 14/12/92, to allow diffusion with convective core mixing
c  (this needs to be checked and refined)
c
c  Modified 19/4/95, to correct error in connection with mixing
c  at boundary of convective core. Previously, the variables
c  cqc and cqcp were erroneously called cqxc and cqxcp in body
c  of routine (but not in common/crxcor/), and hence were not set.
c
c  Modified 21/7/95, to use general heavy-element abundance
c  as set in common/heavy/.
c
c  Modified 5/8/95, to allow diffusion of several species
c
c  Modified 20/5/95, setting Z explicitly to y(5) for 
c  idiffus = 2.
c
c  Modified 4/8/96, taking out option of writing data to disk with
c  ifwrt = 1 (but retaining ifwrt in common/rhcn/ for consistency)
c
c  Modified 5/8/96, adding option imixc1 = 1, for setting rhs
c  consistently with treatment of convective-core mixing (we hope)
c
c  Modified 20/8/96, for option of including centrifugal acceleration
c
c  Modified 18/9/96. Correcting error in constant-coefficient
c  evolution of He3. Introducing icphe3 to switch off doubtful
c  fudge (not used in non-diffusive version) in previously
c  convective region. icphe3 is hardcoded at start of routine.
c
c  Corrected 5/8/97, by including engenr.nnz.d.incl, to set nspdmx etc
c
c  Modified 18/11/98, adding ptpg to call of mixlng (for treatment 
c  of turbulent pressure)
c
c  Modified 19/4/00, to block resetting of composition time derivatives
c  if core was not convective in previous time step (a somewhat 
c  desperate measure)
c
c  Modified 13/7/00, eliminating entropy term in energy generation
c  if it becomes unreasonably large in small convective core
c  (typically as a result of sudden onset of core convection during
c  main-sequence evolution). Parameters etc. in this could do with
c  fine tuning.
c
c  Modified 24/10/00, setting large diffusion coefficient for hydrogen
c  when hydrogen abundance is below 2*xhzlm1. This is to avoide convergence
c  problems when hydrogen is nearly exhausted
c
c  Modified 31/10/00, introducing transformation of
c  luminosity, including a shift in the luminosity variable used
c  for computing when the central hydrogen abundance is very small
c  (hard-coded to .le. 1.e-4). The shift is managed by s/r cpydif.
c
c  Modified 13/11/00, resetting Y variable (diffusion velocity)
c  in outer convection zone.
c
c  Modified 27/11/00, switching off resetting of Y, decreasing
c  diffusion coefficient in convection zones (by 1.e10!)
c
c  Modified 7/4/01, including derivatives of atmospheric factor in 
c  setting derivatives of f(3) (temperature gradient).
c
c  Modified 10/7/02, including option for smoothing diffusion coefficient
c  at edge of convection zone, controlled by ismdif and nsmdif
c
c  Modified 13/8/02, starting inclusion of 4He burning
c
c  Modified 3/9/03, to add more general treatment of convective overshoot
c  from convective envelope and core (the procedure introduced earlier
c  in version without diffusion).
c
c  Modified 12/9/03, correcting resetting of reaction-rate term
c  with diffusion, outside contracting convective core.
c
c  Modified 15/9/03, introducing option imixc4 .ge. 1 to treat
c  evolution of composition of diffusing elements in convective core
c  in the case of diffusion, replacing resetting in s/r cmpcvc.
c
c  Modified 22/9/03, suppressing test for excessive gravitational
c  energy release outside convective cores.
c
c  Modified 30/9/03, setting vhx (previously not set) in definition
c  of settling velocity csdifr(i+idcomp,n) to be passed to prtsol.
c
c  Modified 24/10/03, introducing possibility of semiconvection 
c  (using hydrogen-abundance gradient set in common/cxhder/) 
c  flagged by idiffc1 .gt. 0 (idiffc1 is included in input 
c  parameter idiffct, formerly idiffus).
c
c  Modified 28/10/03, correcting derivatives of rhs of equations
c  for gradient of non-hydrogen diffusing elements.
c
c  Modified 18/5/04, including option of idiffc1 to suppress diffusive
c  mixing, for testing.
c
c  Modified 4/8/04, setting diffusion coefficient to convective value in
c  overshoot region (probably needs further refinement)
c
c  Modified 4/12/07, introducing ismdif = 3: this sets the convective
c  diffusion coefficient at the point closest to the lower boundary of
c  the convective envelope and also (somewhat illogically) flags for 
c  setting X-derivates in fd in convective envelope.
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
      implicit double precision (a-h,o-z)
      logical nosd,notd,dtest,skipt,noder,time0,conv,convos,convx,
     *  lastmd,norct,condif,norche,semicn
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
c
      parameter(idr1mx = nspcmx+3, 
     *  ivrmx1 = ivarmx+1, nspcm2=2*nspcmx, ivrmx4=4*ivarmx,
     *  nsdifr=2*nspdmx+2,
     *  iyfdmx = 9+3*nspdmx+ivarmx*(3+nspdmx+ivarmx*(2+nspdmx)))
c
      dimension yin(ifd),zk(*),f(ifd),fd(ifd,*),alam(ialam,*),
     .  alamd(ialam,ifd,*),h(*),hd(ihd,*),ap(*),aq(*),
     .  zz(iz),dzdy(iz,*),
     .  ddac(6),dlmb(3,3),xhe3(4),dwf(35),comp(nspcmx),y(ivarmx),
     .  yp(ivarmx),zzp(ivarmx),fp(ivarmx),dtst(ivarmx),alamp(1,ivarmx),
     .  cvcpar(10), grad(7), clamxy(4), 
     .  cvvar(icvvar), cvvarp(icvvar),alamrh(4,krnrmx),xche3(4), 
     .  icvflg(nnmax)
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev, syear
      common/ln10/ amm
      common/mxlcn/ cc1,cc2,etac,phc,alfa,tprfct,iconcs,imxlng,iturpr
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
      common/cxhcnc/ compc(nspcmx)
      common/thetac/ theta(ivarmx)
      common/noiter/ iter1, ntime
      common/cmpder/ ft(idr1mx, nspcmx)
      common/prtvar/ rhl,ak,akt,akp,akx,eps(idr1mx),dr,dac,conv,convos
      common/caddvr/ addvar(5,ngmax)
      common/cnvout/ ddrad,rrcnv,aacnv,ddacad,amach,
     *  drr(5), ddr(5), da(5),pturb,ycnv,dyda, dtxh
      common/rnrout/rnzt(10),rnuw(10)
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,inccor,frcf(6),frcl(6)
      common/convpf/ nff(6),nlf(6),incf,inccrf,frcff(6),frclf(6),ifxcon
      common/convvr/ cvvarf(icvvar,6), cvvarl(icvvar,6),
     *  cvvars(icvvar,nnmax)
      common/cnvovs/ icnvos, jcnvos, clcovs, cldovs, alphos,
     *   rczl, rczlfx, rcnvos, qlcnos, rhobos, dmrovs, fmrovs,
     *   roscb1, roscb2, ddoscb, aoscb, boscb, coscb
      common/covsec/ icsove, icsovc, alpove, alpovc, 
     *  qmxove, qmxovc, rmxove, rmxovc, imxove, imxovc, 
     *  nmxove, nmxovc
      common/cvcpsc/ cvcpst(20)
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor,nmxcp,qmxcp,rmxcp,
     *  qmxmin,qmxmax
      common/crxcor/ cqc, cqcp, crxmn(nspcmx), crxmnp(nspcmx)
      common/cmxcnt/ ddrmix, inmixc, imixc0, imixc1, imixc2, imixc3,
     *  imixc4
      common/crxstr/ rxstr(nspcm2,1)
      common/cdiffu/ idiffus, itbdif, rhob, rmxdif, tlfdif, dtldif,
     *  ctdfmx, rctdif, rctdf1, idfxbc, ismdif, nsmdif, idiffc1, idiffc2
      common/cdiffv/ cdifr(nspdmx),cvelr(nspdmx,2),cdiftr
      common/cdifsv/ csdifr(nsdifr,1)
      common/cdiffc/ velx(nspdmx,2,7), difx(nspdmx,7), difscx(7)
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
     *  isphec, idcomp, iccomp
      common/cxhder/ xhder(2,nnmax)
      common/engcnt/ xhzlm1, xhzlm2, nxhzer, icnocs, nspec,
     *  iheccs, nspect
      common/enggrv/ epsg(nnmax)
      common/degfct/ thte, iscren
      common/eqstcl/ iradp,mode,dtest,skipt,noder
      common/eqprcl/ dmeqpr(12),idiag,iscskp
      common/eqstd/ xi(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     .  dlt(4),gm1,tprh,trhp,rhxp
      common/opcdat/ akk(5)
      common/nmbmsh/ nn
c
c  NOTE: Storage in yzfstr needs rethinking with diffusion
c
      common/cyfstr/ yzfstr(ivrmx4,nnmax)
c
      common/yprtst/ yprt(ivrmx1,nnmax)
      common/fprtst/ fprt(ivarmx,nnmax)
      common/aprtst/ auxprt(20,nnmax)
      common/opcxdr/ akxa, akz
      common/cdfprt/ prtdvr(10,nnmax)
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
      data flpr,tlpr,xhpr,epsskt /-100.,-100.,-100.,1./
      data initrs /0/
c
c  factor for diffusion coefficient in convection zones
c  (set to dfconfc*amms**2, where amms is mass of star, in g).
c  The previous value was dfconfc = 1.d-10.
c
      data dfconfc /1.d-15/
c
      data init /1/
c..      data rczl /1.d20/
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
      kdgrhb=0
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
c  test for initializing flag for convective regions
c
      if(iter1.eq.1) then
	do i=1,nn
	 icvflg(i)=0
        end do
      end if
c
c  hardcode flag for fudge in previously convective region
c
      icphe3=0
c
c  store yin (as input) in internal file yin
c
      call store(yin,y,ivarmx)
c
      qx=10.d0**x
      radius=1.d11*10.d0**y(1)
c
c  when not including convective overshoot, set rcnvos to
c  impossible value
c
      if(icnvos.ne.1.and.icnvos.ne.2) rcnvos = -1.d10
c
      if(init.eq.0) then
	rczl=1.d20
	init=1
      end if
c
c  set surface radius into rs, make other initializations
c
      if(n.eq.1) then
	if(iter1.le.1) iexcon=0
	izerf3=0
	icqcdg=1
	rs=1.d11*10.d0**y(1)
	hhpcz=1.
c
c  store density at base of convection zone at previous time step,
c  for turbulent diffusion
c
	if(iter.eq.0) then
	  if(icnvos.eq.1.or.icnvos.eq.2) then
	    rhobp=rhobos
          else
	    rhobp=cvvarl(4,1)
	  end if
	  rmxdfp=rs*rmxcor

	  if(istdpr.gt.0) write(istdpr,*) 
     *      'In rhs, set rhobp, rmxdfp to',rhobp, rmxdfp
        end if
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
c  new hydrogen abundance
c
      if((ntime.eq.0.or.iter.eq.0).and.n.eq.1) nxhzer=nn+1
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
c  set storage index for luminosity
c
      ilum=4+idcomp
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
c  for idiffus = 2, store Z from y into zh(n), for use in this
c  routine and later
c
      if(idiffus.eq.2) then
	zh(n)=y(5)
      end if
c
c  calculate log(rho),log(kappa) and log(eps)
c
      fl=y(2)
      tl=y(3)
      rlred=y(1)
c
      irsetx=0
      if(y(4).lt.1.d-10) then
	y(4)=1.d-10
	irsetx=1
      end if
c
      if(iter1.gt.1.and.nmxcor.gt.0.and.qx.le.qmxcor) then
	if(imixc4.eq.0.or.idiffus.eq.0) then
          xh=compc(1)
	  zhh=zhc
	else if(idiffus.eq.1) then
          xh=y(4)
	  zhh=zhc
	else if(idiffus.eq.2) then
          xh=y(4)
	  zhh=zh(n)
	end if
      else
        xh=y(4)
	zhh=zh(n)
      end if
c
      allred=y(ilum)
c
      if(n.eq.1) allrds=allred
c
      yh=1-xh-zhh
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
c     if(n.eq.118) idiag=1
      call eqstf(fl,tl,xh,yh,zhh,nosd,notd)
      if(kdgeos.lt.0) then
	kdgrhb=-1
	return
      end if
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
c  store log p in addvar
c
      addvar(1,n)=pl
c
      if(idgrhs.eq.-1.and.istdpr.gt.0) write(istdpr,1100) 
     *  'just after setting', (zz(i),i=1,5)
c
      if(.not.time0) then
c
c  Note that setting of X into zz may be reset below,
c  depending on whether or not this is in region
c  where X should be set to zero, or in convection zone,
c  where X must be decoupled from other variables.
c
        zz(4)=xh
        dzdy(4,4)=1.d0
        dzdy(2,4)=pt(4)
c
	do 13 i=5,nvar
        zz(i)=y(i)
  13    dzdy(i,i)=1
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
      r=1.d11*10.d0**rlred
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
c
c  set composition variables.
c  Note: in a convective core, this is set directly from compc
c
      if(iter1.gt.1.and.nmxcor.gt.0.and.qx.le.qmxcor) then
	call store(compc,comp,icomp)
c
c  reset X when imixc4 .ge. 1 with diffusion, to use actual hydrogen
c  abundance in the core. Note that this has to be though through
c  with care, particularly when more reacting elements are included
c  with diffusion
c
	if(imixc4.ge.1.and.idiffus.ge.1) comp(1)=y(4)
      else
c
c  general readiative case. Use abundances in y
c  (could also do with a little cleaning-up!)
c
        comp(1)=y(4)
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
          if(iter.ne.0) then
            y(5+2*idcomp)=xhe3(1)
            if(istdpr.gt.0) write(istdpr,*) 
     *        'ifdhe3. set y(',5+2*idcomp,') to',xhe3(1)
          end if
          comp(2)=xhe3(1)
          notd=.true.
	  ixx3 = 2
        else if(ieqhe3.ne.1) then
	  comp(2) = y(5+2*idcomp)
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
        else if(icnocs.ge.1) then
	  call store(y(4+2*idcomp+ispxx3),comp(ixx3+1),ispcno)
        end if
      end if
c
      notd=.true.
      idgeno = idgeng
c
c  as a rather unsatisfactory fudge, use He3 equilibrium in
c  inside convective core in previous model
c
c  *** Note: only implemented if icphe3 = 1 (icphe3 hard-coded
c      at start)
c
      ieqheo=ieqhe3
      if(qx.le.cqcp.and.iter.ne.0.and.icphe3.eq.1) then
	ieqhe3=1
	if(ieqheo.ne.1) comp(2)=comp(3)
	if(icqcdg.eq.1.and.istdpr.gt.0) then
	  write(istdpr,*) 
     *      'Reset He3 in former convective core, ieqhe3 =', ieqheo
	  write(istdpr,*) 'Call engenr with fl, tl, comp =',fl, tl, 
     *      (comp(i),i=1,3)
        end if
      end if
c
      call engenr(fl,tl,comp,yh,zhh,eps,ft,ift,notd)
c..      call zero(ft,idr1mx*nspcmx)
c..      call zero(crxmn,nspcmx)
c..      call zero(crxmnp,nspcmx)
      if(kdgeng.lt.0) then
	kdgrhb=-1
	return
      end if
      if(qx.lt.-0.001.and.istdpr.gt.0) 
     *  write(istdpr,'(a,2f10.5,1p4e11.3)') 
     *  'After engenr fl, tl, comp(1-3), ft(1,1) =',fl, tl, 
     *      (comp(i),i=1,3),ft(1,1)
      if(qx.le.cqcp.and.ieqheo.ne.1.and.iter.ne.0
     *    .and.icphe3.eq.1) then
	ieqhe3=ieqheo
	do 17 j=1,idr1mx
   17   ft(j,2)=0
	ft(6,1)=ft(5,1)
	ft(5,1)=0
	eps(6)=eps(5)
	eps(5)=0
	comp(3)=comp(2)
	comp(2)=xhe3eq
	y(5+2*idcomp)=xhe3eq
	zz(5+2*idcomp)=xhe3eq
	if((icqcdg.eq.1.or.idgeng.ge.1).and.istdpr.gt.0) then
	  write(istdpr,*) 'qx, cqcp =',qx, cqcp
	  write(istdpr,*) 'xhe3eq, eps =',xhe3eq,eps(1)
	  icqcdg=0
        end if
      end if
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
      ihe3=5+2*idcomp
      xhe3or=y(ihe3)
      if(ifdhe3.eq.0.and.iche31.eq.1.and.iter.ne.0) then
        do 19 k=1,4
        alamrh(1,k)=rho(1)*al(1,k)
	do 19 j=2,4
   19   alamrh(j,k)=rho(j)+al(j,k)
	nosd=.false.
	call he3abc(xh,yh,zhh,dt,yprt(ihe3,n),alamrh,xche3,xhe3eq,
     *    anu,4,nosd)
c
	ft(1,2)=(xche3(1)-yprt(ihe3,n))/dt
	do 20 i=2,4
   20   ft(i,2)=(xche3(i)-yprt(ihe3,n))/dt
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
c  test for setting nxhzer, based on values at previous time step
c  for iter = 0. Assume that rx varies linearly with X,
c  so that time dependence is exponential.
c
      if(iter.eq.0.and.nxhzer.gt.nn) then
	if(y(4)*exp(dt*ft(1,1)/y(4)).le.xhzlm1) then
	  nxhzer=n
	  if(istdpr.gt.0) write(istdpr,*) 'nxhzer set to ', nxhzer, 
     *      '  in rhs'
        end if
      end if
c
c  test for resetting energy generation quantities to zero, where
c  hydrogen abundance is too low
c  Also reset values stored in zz and dzdy.
c
      if(n.ge.nxhzer.or.irsetx.eq.1) then
	zz(4)=1.e-10
	zz(5+idcomp)=0.d0
        dzdy(2,4)=0.d0
	do 21 i=1,idr1mx
	eps(i)=0
	do 21 k=1,nspcmx
   21   ft(i,k)=0
      end if
c
c  test for setting correction factor in temperature gradient for smooth
c  transition from atmospheric T(tau) relation
c
      if(icsrad.ne.1.or.tl.gt.4.3) then
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
      f(3)=-a3*tgrfct*10.d0**(ak+x-4.d0*(tl+rlred))*
     *        ((10.d0**allred)-alshft)
c
      f(ilum)=a4*eps(1)*10.d0**(x-allred)
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
c  test for convection, possibly including test for semiconvection
c
      dr=f(3)/f(2)
c..      if(idiffc1.gt.0.and.(iter.eq.0.or.iter1.ge.3)) then
      if(idiffc1.gt.0) then
c
c  hydrogen gradient for semiconvection
c
	if(n.ge.nxhzer) then
	  dxh=0.d0
	else if(iter.eq.0) then
	  dxh=xhder(1,n)
        else
	  dxh=xhder(2,n)
	end if
	dtxh=-trhp*rhxp*dxh/(amm*xh*f(2))
        semicn=dr.gt.dad(1).and.dr.lt.dad(1)+dtxh
c
c  suppress semiconvection for excessive number of iterations
c  where previous iteration did not have convection or semiconvection
c
	if(iter1.ge.10.and.semicn.and.icvflg(n).eq.0) then
	  semicn=.false.
	  if(istdpr.gt.0) write(istdpr,1165) iter1,n
        end if
      else
	dtxh=0.d0
        semicn=.false.
      end if
      dac=dr
      ddrad=dr-dad(1)-dtxh
c
c  pressure scale height
c
      hhp=-r*f(1)/f(2)
c
c  attempt to stabilize boundary of convective region, if
c  excessive number of iterations
c
      conv=dr.gt.dad(1)+dtxh
      if(iter1.gt.10.and.conv) then
        if(n.gt.nl(inc).and.dcnp.lt.0) then
	  drnew=dad(1)+dtxh-1.e-5
	  if(istdpr.gt.0) write(istdpr,1160) iter1,dr,drnew,n
	  dr=drnew
	  if(icnvos.eq.1.or.icnvos.eq.2) convos=.true.
	  conv=.false.
        end if
      end if

c
      drdad=dr-dad(1)
      ddacad=dr-dad(1)-dtxh
      ddarad=dr-dad(1)-dtxh
c
c  reset flag for convection
c
      if(conv) then
	icvflg(n)=1
      else if(semicn) then
	icvflg(n)=-1
      else
	icvflg(n)=0
      end if
c
c  skip setting of convective boundaries when iter = 0
c
      if(iter.eq.0) go to 28
c
c  set differences for test of convective boundary
c
      dcnp=dcn
      dcn=-ddacad
c
      call store(cvvar,cvvarp,icvvar)
      amach=0
      if(n.lt.0) go to 28
c  find edges of convection zone
c
c  nnf=nf(inc) is the first meshpoint (counted from the surface), and
c  nnl=nl(inc) is the last inside or on the edge of the inc-th
c  convection zone. further if frc=frcf(inc) (or frcl(inc)) then
c  for any function a=a(n) frc*a(nnf)+(1-frc)*a(nnf-1)
c  (or frc*a(nnl)+(1-frc)*a(nnl+1)) is the value of a, from linear
c  interpolation, at the edge of the convection zone.
c  Variables (m/M, r/R, p, rho, T, kappa) at edges of
c  convection zones are stored in cvvarf(.,inc) and cvvarl(.,inc)
c
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
      cvvar(11)=dtxh
      call store(cvvar,cvvars(1,n),icvvar)
c
      if(n.eq.1) then
c  
c  surface point
c  
c  extra variables for output near boundary of convective core
c
	if(inc.gt.0.and.nl(inc).eq.nn.and.istdpr.gt.0.and.idgcon.gt.0) 
     *    then
	  ncp1=nf(inc)-3
	  ncp2=nf(inc)+3
	  write(istdpr,*) 
     *      ' n, log q, log T, log rho, X, log kappa, ddad'
        else
	  ncp1=0
	  ncp2=0
        end if
c
	if(idgcon.gt.0) then
	  nep1=nn+1
	  nep2=nep1+20
        else
	  nep1=0
	  nep2=0
        end if
c
        inc_prev=inc
        inc=0
	dcnp=0
	call zero(cvvarp,7)
        inl=0
        dcs=dcn
        n1=n
        frc=1
	if(dcs.le.0) then
c
c  start of convection zone
c
          inc=inc+1
          nf(inc)=n1
          frcf(inc)=frc
	  call store(cvvar,cvvarf(1,inc),icvvar)
	  if(istdpr.gt.0) write(istdpr,*) 'inc, cvvarf(.,inc):',
     *      inc,(cvvarf(i,inc),i=1,icvvar)
        else
	  go to 28
	end if
      end if
c
c  extra output
c
      if(n.ge.ncp1.and.n.le.ncp2.and.istdpr.gt.0) then
	write(istdpr,'(i5,6f11.6)') n, x, tl, rhl, xh, ak, ddacad
      end if
      if(n.ge.nep1.and.n.le.nep2.and.istdpr.gt.0) then
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
          inc=inc+1
	  if(istdpr.gt.0) then
	    write(istdpr,'(/a,i3)') 
     *        'change of sign. increment inc to',inc
            write(istdpr,'(3i5,1p8e13.5)')
     *        ntime,iter1,n,qx,rhl, tl, xh, ak, dtxh, ddrad,allred
	  end if
          inc=min0(6,inc)
          nf(inc)=n1
          frcf(inc)=frc
c
	  do 23 i=1,icvvar
   23     cvvarf(i,inc)=frc*cvvar(i)+(1-frc)*cvvarp(i)
	  if(istdpr.gt.0) then
	    write(istdpr,*) 'cz start. inc, n, dcnp, dcn, frcf =',
     *      inc, n, dcnp, dcn, frc
	    write(istdpr,*) 'x, y:',x,(y(i),i=1,6)
c
            write(istdpr,'(/a,i3,1p12e13.5)') 
     *        'inc, cvvarf(.,inc):',inc,(cvvarf(i,inc),i=1,icvvar)
	  end if
	  if(inc.eq.2.and.istdpr.gt.0) then
	    write(istdpr,*) ' composition:',(comp(i),i=1,3)
	    write(istdpr,*) ' scaled L/M:',10.**(y(ilum)/x)
          end if
c
c  store parameters for output of quantities relevant to convective
c  core
c
          call store(cvcpar,cvcpst,6)
c
          cvcpst(7)=qx
          cvcpst(8)=rhl
          cvcpst(9)=tl
          cvcpst(10)=xh
          cvcpst(11)=ak
          cvcpst(12)=ddacad
	  do 25 i=1,6
   25     cvcpst(12+i)=frc*cvcpst(6+i)+(1-frc)*cvcpst(i)
          if(istdpr.gt.0) write(istdpr,*) 'cvcpst:',(cvcpst(i),i=1,18)
        else if(dcn.gt.0) then
c 
c  end of convection zone
c
	  nl(inc)=n1-1
          frcl(inc)=1-frc
c
	  do 26 i=1,icvvar
   26     cvvarl(i,inc)=frc*cvvar(i)+(1-frc)*cvvarp(i)
	  if(istdpr.gt.0) write(istdpr,*) 'inc, cvvarl(.,inc):',
     *      inc,(cvvarl(i,inc),i=1,icvvar)
c
	  if(istdpr.gt.0) then
	    write(istdpr,*) 'cz end. inc, n, dcnp, dcn, frcl =', 
     *        inc, n, dcnp, dcn, 1-frc
	    write(istdpr,*) 'x, y:',x,(y(i),i=1,6)
	  end if
	  if(inc.eq.1) then
	    rczl = r
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
            cvcpst(7)=qx
            cvcpst(8)=rhl
            cvcpst(9)=tl
            cvcpst(10)=xh
            cvcpst(11)=ak
            cvcpst(12)=ddarad
	    if(istdpr.gt.0) write(istdpr,*) 'cvcpst:',(cvcpst(i),i=1,12)
          else 
c 
c  end of convection zone
c
	    n1=n+1
	    nl(inc)=n1-1
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
      end if
c
c  test for convection to centre 
c  
      if(n.eq.nn.and.dcn.lt.0) then  
        frcl(inc)=1
        nl(inc)=nn 
        call store(cvvar,cvvarl(1,inc),icvvar)
	inccor=inc
	if(istdpr.gt.0) write(istdpr,*) 'Set inccor to',inccor
      end if
c  
c  in convection zone or overshoot region reset f(3)
c
   28 convos = .false.
c
c  for idiffc2 .gt. 0, suppress resetting of temperature gradient
c  in all but the first convection zone
c
      if(conv.and.idiffc2.gt.0.and.inc.gt.1) then
	conv=.false.
	if(istdpr.gt.0) write(istdpr,*) 
     *    '#D# Suppress convection at n, x =', n, x
      end if
c
      drdadp=drdad
c
c  test for new (3/9/03) formulation of overshoot
c
      if(.not.conv) then
        if(imxove.gt.0.and.r.lt.rs*cvvarf(2,imxove).and.r.ge.rmxove) 
     *    then
	  icnovs = icsove
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
	if(idiffus.eq.2) then 
          akk(5)=akz/(amm*y(5))
        else
	  akk(5)=0
        end if
        call mixlng(t,hhp,dr,dac,ddac,ptpg)
c
        f(3)=dac*f(2)
        if(idgrhs.eq.1.and.istdpr.gt.0) write(istdpr,1170) 
     *    f(2),f(3),dac
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
        call oversh(t,hhp,dr,dac,ddac,xovs)
c
        f(3)=dac*f(2)
        rcnvos = r
	rhobos=rho(1)
	if(istdpr.gt.0) write(istdpr,*) 
     *    ' Now rcnvos/R, rhobos =',rcnvos/rs, rhobos
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
c
      if(icnvos.eq.2.and.jcnvos.eq.2.and.iter1.gt.1) then
	if(n.eq.1) call iovsmr(iter1)
	drovs=(r-rczlfx)/rs+cldovs
	if((drovs.le.3.*dmrovs.or.r.le.rczlfx*1.001).and.
     *    drovs.ge.-clcovs.and.(semicn.or.dac.le.dad(1)+1.d-4)) then
	  xovs=((r-rczlfx)/rs+cldovs)/dmrovs
	  if(istdpr.gt.0) write(istdpr,*) 'OVS',
     *      n, conv, r/rs, rczlfx/rs, dad(1), dr, hhp, xovs
	  call oversh(t,hhp,dr,dac,ddac,xovs)
	  f(3)=dac*f(2)
	  convos=.true.
        end if
      end if
c
c  test for using simple cubic subadiabatic region
c
      if(icnvos.eq.3.and.iter1.gt.1) then
	if(n.eq.1) call iovscb(iter1)
        xovs=(r/rs-roscb1)/(roscb2-roscb1)
	if(xovs.ge.0.and.xovs.le.1) then
	  call oversh(t,hhp,dr,dac,ddac,xovs)
	  if(istdpr.gt.0) write(istdpr,*) 'dac =',dac
	  f(3)=dac*f(2)
        end if
      end if
c
      cvvar(10)=dac
      cvvars(10,n)=dac
c
c  store parameters for later analysis of convective core boundary
c
      cvcpar(1)=qx
      cvcpar(2)=rhl
      cvcpar(3)=tl
      cvcpar(4)=xh
      cvcpar(5)=ak
      cvcpar(6)=ddacad
c
c  set flag for "infinite" diffusion coefficients in convection zones
c  For imixc1 .gt. 0 we exclude convective core, where diffusion
c  coefficient is reset below only if core is otherwise treated
c  as mixed. The current criterion is based on assuming that there
c  is an external convection zone, and that any superficial second
c  convection zone does not extend very deeply. This needs more
c  careful thought.
c
c  Replace test variable by convx, fixed when the convection-zone 
c  boundary is fixed to a mesh point with istrt2 .ge. 1.
c  This is flagged by ifxcon set by s/r rczmsh (and initialized to
c  zero in evolmain). This assumes that the actual point is very
c  close to the fixed boundary, and that the call is not at the
c  previous time step. Otherwise the normal criterion is set in 
c  convx also.
c
      if(ifxcon.ge.1.and.iter.ne.0.and.iabs(n-nlf(1)).le.2) then
	convx=n.le.nlf(1)
	if(convx.and.istdpr.gt.0) write(istdpr,
     *    '(/'' Convective mixing forced at n ='',i5,
     *    '' x ='',1pe13.5,'' dr -dad ='',e13.5)') n, x, dr-dad(1)
      else
	convx=conv
      end if
c
      if(imixc1.eq.0.or.imixc4.ge.1) then
	condif=convx
      else
	condif=convx.and.(inc.eq.1.or.qx.ge.0.99)
      end if
c
c  Set condif if boundary of convection envelope is closer to the present
c  point, for ismdif = 3, with some freezing options.
c 
      if(iter1.gt.1.and.ismdif.eq.3.and.inc.eq.1.and.
     *  n.eq.nl(inc)+1) then
c..	write(istdpr,*) '#D# iter1, n, frcl(inc), iexcon =',
c..     *    iter1, n, frcl(inc), iexcon
        if(frcl(inc).lt.0.5d0) then
	  iexcon = 1
	  condif=.true.
	else if(iter1.ge.5.and.iexcon.eq.1.and.frcl(inc).lt.0.6d0) then
	  condif=.true.
        else
	  iexcon=0
	  condif=.false.
          if(istdpr.gt.0) write(istdpr,*) 
     *      '#D# unset condif at n, iter1 =',n, iter1
        end if
	if(condif.and.istdpr.gt.0) 
     *    write(istdpr,*) '#D# set condif at n, iter1 =',n, iter1
      end if
c
      if(idgrhs.eq.-1.and.istdpr.gt.0) write(istdpr,1100) 
     *  'after testing for convection', (zz(i),i=1,5)
      if(idgeng.ge.2.and.istdpr.gt.0) then
        write(istdpr,*) 'after testing for convection'
        call dmprcn
      end if
c
c  set diffusion quantities
c
      grad(1)=dac
      if(conv) then
	do 29 i=2,6
   29   grad(i)=ddac(i-1)
      else
        akf=rkr*rho(2)
        akt=rkt+rkr*rho(3)
        akx=akxa+rkr*rho(4)
	grad(2)=0
	grad(3)=amm*grad(1)*(akf+pt(2))
	grad(4)=amm*grad(1)*(akt+pt(3)-4)
	grad(5)=amm*grad(1)
	grad(6)=amm*grad(1)*(akx+pt(4))
	if(idiffus.ge.2) then 
	  grad(7)=grad(1)*akz/y(5)
        else
	  grad(7)=0
        end if
      end if
      nosd=.false.
      amms=amsun*am
      amass=amms*qx
      if(iter.eq.0) then
	rhob=rhobp
	rmxdif=rmxdfp
      else
	rhob=cvvarl(4,1)
	rmxdif=rs*rmxcor
      end if
      call difcff(t,rho,pt,xh,zhh,r,amass,grad,
     *  idiffus,itbdif,clamxy,velx,difx,nspdmx,nosd)
c
c  possibly set diffusion from semi-convection
c
      if(semicn) then
	call difscn(t,rho,pt,xh,zhh,r,amass,grad,ak,dtxh,difscx,nosd)
      else
	call zero(difscx,6)
      end if
c
c..      do k=1,6
c..	velx(2,2,k)=0
c..      end do
c
c  store variables for later analysis
c
      prtdvr(1,n)=x
      prtdvr(2,n)=y(4)
      prtdvr(3,n)=y(5)
      prtdvr(4,n)=y(5+idcomp)
      prtdvr(5,n)=y(6+idcomp)
      prtdvr(6,n)=difx(1,1)/(amms*amms)
      prtdvr(7,n)=difx(2,1)/(amms*amms)
      prtdvr(8,n)=velx(1,1,1)/amms
      prtdvr(9,n)=velx(2,1,1)/amms+velx(2,2,1)*y(5+idcomp)
c
c  loop over diffusing elements
c
      dmdr=4.d0*pi*rho(1)*r*r
      vhx=amms*y(5+idcomp)
c
      do 30 j=1,idcomp
c
c  reset diffusion velocity in convection zone or close to surface
c  Note that this depends on the value of imixc1, the treatment
c  used before 5/8/96 corresponding to imixc1 = 0
c
c  On 24/6/02, add smoothing near lower boundary of convective envelope,
c  based on value of superadiabatic gradient.
      if(inc.ge.1.and.convx) then
	if(ismdif.eq.1.and.n.eq.nl(1)-nsmdif) then
	  ddrads=ddrad
	  difcon=dfconfc*amms*amms
        else if(ismdif.eq.1.and.n.gt.nl(1)-nsmdif) then
	  diffct=(ddrad/ddrads)**2
          difcon=1.d0/(diffct/(dfconfc*amms*amms)
     *                         +(1.d0-diffct)/difx(j,1))
        else if(idiffc1.eq.1) then
  	  difcon=dmdr*dmdr*alfa*hhp*amach*sqrt(pt(1)*gm1/rho(1))/3.d0
        else
	  difcon=dfconfc*amms*amms
        end if
      else
	difcon=dfconfc*amms*amms
      end if
c
c  As a simple fudge here, suppress any resetting in core
c  if idiffc1 = 9
c
      if(idiffc1.eq.9.and.x.lt.-1.) then
c
c  intentionally empty statement
c
      else if(imixc1.eq.0) then
        if(convx.or.(nmxcor.gt.0.and.qx.le.qmxcor).or.x.gt.-1.e-8) then
	  difx(j,1)=difcon
c..	  call zero(difx(j,2),3)
          do k=2,6
            difx(j,k)=0
          end do
        end if
c
c  new treatment, separating test for convection and core mixing,
c  and fixing treatment of previous time step
c
c  Also try to reset diffusion velocity (Y variable) appropriately
c  in z, to stabilize solution.
c
c  Now switched off again.
c
c  Change test for convection to be based on conv for the first
c  few iterations at the present timestep, even in a convective core.
c
      else
	istzer=0
        if(iter.eq.0) then
	  if(condif.or.x.gt.-1.e-8.or.(nmxcp.gt.0.and.qx.le.cqcp)) then
	    difx(j,1)=difcon
	    do k=2,6
              difx(j,k)=0
            end do
	    if(istzer.eq.1) then
	      velx(j,1,1)=0
	      velx(j,2,1)=0
	      do 29013 k=2,6
	      velx(j,1,k)=0
29013         velx(j,2,k)=0 
	    end if
	  end if
        else
c  Note: in test cqc replaced by qmxcor on 17/9/03
          if(condif.or.(iter1.lt.3.and.conv)
     *      .or.(iter1.ge.3.and.nmxcor.gt.0.and.
     *      qx.le.qmxcor)
     *      .or.x.gt.-1.e-8) then
	    difx(j,1)=difcon
	    do 29015 k=2,6
29015       difx(j,k)=0
	    if(istzer.eq.1) then
	      velx(j,1,1)=0
	      velx(j,2,1)=0
	      do 29017 k=2,6
	      velx(j,1,k)=0
29017         velx(j,2,k)=0 
	    end if
c
          end if
        end if
      end if
c
c  For now, set diffusion coefficient to convection-zone value in
c  overshoot region (assuming that difcon has already been set in
c  the convection zone). This will need elaboration.
c
      if(convos) then 
        difx(j,1)=difcon
        do k=2,6
          difx(j,k)=0
        end do
      end if
c
c  as an attempt at stabilizing convection at boundary of convective
c  envelope, include option of resetting diffusion coeffient for the
c  first few points beneath the convection zone, for ismdif = 2.
c  Certainly needs elaboration (if it works!)
c
      if(ismdif.eq.2.and.n.ge.nl(1).and.n.le.nl(1)+nsmdif) then
        difx(j,1)=difcon
        do k=2,6
          difx(j,k)=0
        end do
      end if
c
c  for ismdif = 4, try to stabilize solution by extending diffusively
c  mixed region below convective envelope
c
      if(ismdif.eq.4.and.iter1.ge.10.and.inc_prev.gt.2.and.
     *  inc.ge.1.and.n.ge.nl(inc).and.n.le.nl(inc)+nsmdif.and.n.lt.nn)
     *  then
        difx(j,1)=difcon
        if(istdpr.gt.0) write(istdpr,*) '#D# Reset difx at n =',n
        do k=2,6
          difx(j,k)=0
        end do
      end if
c
c  when hydrogen content is small, zero settling velocity
c  for hydrogen (rather a desperate measure!)
c
      if(j.eq.1.and.xh.lt.10*xhzlm1) then
	velx(j,1,1)=0
      end if
c
c  possibly include semi-convective contribution
c  Suppress in core for idiffc1 = 9.
c
      if(semicn.and.(idiffc1.ne.9.or.x.ge.-1.)) then
	if(j.eq.1.and.istdpr.gt.0) write(istdpr,*) 
     *    'SEMIC. 2. difx, difscx:', difx(j,1),difscx(1)
	do i=1,7
	  difx(j,i)=difx(j,i)+difscx(i)
	end do
      end if
c
      if(n.gt.1.and.n.lt.nn) then
	if(j.eq.1) then
          f(3+j)=amm*qx*(amms*amms/difx(j,1))*
     *           (y(4+idcomp+j)-(velx(j,1,1)/amms)*y(3+j))
	  if(iter1.eq.1.and.n.ge.400.and.n.le.500.and.mod(n,5).eq.-1
     *      .and.istdpr.gt.0) 
     *      write(istdpr,'(a,i5,1p7e15.7)') 'f(4)',
     *           n,amms,difx(j,1),velx(j,1,1),
     *           amm*qx*(amms*amms/difx(j,1)),
     *           y(4+idcomp+j),(velx(j,1,1)/amms)*y(3+j), f(3+j)
	  if(n.ge.nep1.and.n.le.nep2+10.and.istdpr.gt.0) 
     *      write(istdpr,'(a,i5,1p4e15.7)') 'f(4)',
     *           n,amm*qx*(amms*amms/difx(j,1)),
     *           y(4+idcomp+j),(velx(j,1,1)/amms)*y(3+j), f(3+j)
        else
          f(3+j)=amm*qx*(amms*amms/difx(j,1))*
     *           (y(4+idcomp+j)-
     *           (velx(j,1,1)/amms+velx(j,2,1)*y(5+idcomp))*y(3+j))
	end if
      else
        f(3+j)=0
      end if
c
c  For clarity, as a separate statement, zero hydrogen gradient
c  where X is exhausted
c
      if(n.ge.nxhzer.and.j.eq.1) f(3+j)=0.d0
c
c      ---------------------------------------------------
c
c  set reaction-rate term for diffusing elements 
c  (so far, only hydrogen is assumed to react; this will need fixing
c  up if, for example, diffusion of CNO elements or He3 is included).
c
      if(j.eq.1) then
        f(4+idcomp+j)=-amm*qx*ft(1,1)
      else
        f(4+idcomp+j)=0
      end if
c
c  reset cdifr if convection or other types of turbulence
c  have been involved
c
      cdifr(j)=difx(j,1)/(dmdr*dmdr)
c
c  store diffusion variables for output to gong file
c
      csdifr(j,n)=cdifr(j)
      if(j.eq.1) then
        csdifr(j+idcomp,n)=cvelr(j,1)
      else
        csdifr(j+idcomp,n)=cvelr(j,1)+vhx*cvelr(j,2)
      end if
c
   30 continue
c
      csdifr(2*idcomp+1,n)=cdiftr
c
      if(semicn) then
	csdifr(2*idcomp+2,n)=difscx(1)/(dmdr*dmdr)
      else
	csdifr(2*idcomp+2,n)=0
      end if
c
c  set time derivatives for non-diffusing species (this must be checked
c  if other burning elements than hydrogen are allowed to diffuse)
c
      if(iccomp.gt.0) then
        do 32 k=1,iccomp
   32   f(4+2*idcomp+k)=ft(1,k+1)
      end if
c
c  store in fprt and auxprt
c
      if(iter.ne.0) then 
        do 32090 k=1,nvar
32090   fprt(k,n)=f(k)
        auxprt(1,n)=rho(1)
        auxprt(2,n)=pt(1)
        auxprt(3,n)=dad(1)
        auxprt(4,n)=10.d0**ak
        auxprt(5,n)=dr
        auxprt(6,n)=comp(1)
        auxprt(6,n)=comp(2)
        auxprt(7,n)=comp(3)
        auxprt(8,n)=eps(1)
      end if
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
      do 33 k=1,icomp
   33 rxstr(iftst+k,n)=ft(1,k)
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
c..          f(5+idcomp)=-amm*qx*crxmnp(1)
	  if(iccomp.gt.0) then
	    do 34 k=1,iccomp
   34       f(4+2*idcomp+k)=-amm*qx*crxmnp(k+1)
	  end if
	end if
      else
c
c  current time step (note that this is independent of 
c  imixc1, so far)
c
        if(inmixc.ne.1.and.nmxcor.gt.0) then
	  if(qx.le.qmxcor) then
c..            f(5+idcomp)=-amm*qx*crxmn(1)
	    if(iccomp.gt.0) then
	      do 35 k=1,iccomp
   35         f(4+2*idcomp+k)=crxmn(k+1)
	    end if
	  else if(qx.le.qmxcp) then
	    tfrct=(qmxcp-qx)/(qmxcp-qmxcor)
            f(5+idcomp)=(1-tfrct)*f(5+idcomp)-tfrct*amm*qx*crxmn(1)
c..            f(4)=(1-tfrct)*f(4)-tfrct*amm*qx*crxmn(1)
	    if(iccomp.gt.0) then
	      do 37 k=1,iccomp
              if(iche31.ne.1.or.k.ne.1)
     *          f(4+2*idcomp+k)=
     *              (1-tfrct)*f(4+2*idcomp+k)+tfrct*crxmn(k+1)
   37         continue
	    end if
	  end if
	end if
      end if
c
      if(idgrhs.eq.-1.and.istdpr.gt.0) write(istdpr,1100) 
     *  'after setting f(5 - nvar)', (zz(i),i=1,5)
c
      if(irhtst.eq.1.and.iter.eq.0) then
	yprt(7,n)=f(ilum)
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
        alam(1,1)=0
        alam(1,2)=bb*(prhm*pt(2)-ht(2))
        alam(1,3)=bb*(prhm*pt(3)-ht(3))
        alam(1,4)=bb*(prhm*pt(4)-ht(4))
	do 48 i=5,nvar
   48   alam(1,i)=0
      end if
c
c  set epsilong
c
c  Note: for the time being assume that backwards differences are
c  used. Hence do not centralize the lambda-s properly.
c  This should be fixed up at some stage.
c
      epsg(n)=0
      if(iter.gt.0.and.iter1.gt.0) then
	epsg(n)=0
	do 49 i=1,nvar
   49   epsg(n)=epsg(n)+alam(1,i)*(y(i)-yprt(i,n))/(bb*dt)
      end if
c
c  test for unreasonably large epsilong in small convective core
c  or elsewhere
c
c  Modify test to exclude radiative region, 22/9/03
c  This resetting has to be reconsidered with care.
c
      epsav=1.d33*((10.d0**allred)-alshft)/(qx*amms)
      if((abs(epsg(n)).gt.0.02*epsav.and.
     *  (qx.le.1.1*cqc.or.qx.le.1.1*cqcp).and.xh.ge.0.01)) then
c..     *  .or.abs(epsg(n)).gt.2*epsav) then
c
	if(istdpr.gt.0) write(istdpr,'(a,i5,1p3e11.3)') 
     *    ' Excessive epsg. n, q, epsg, epsav:',
     *    n, qx, epsg(n),epsav
	inenti=1
	epsg(n)=0
	do 50 i=1,nvar
  50   alam(1,i)=0
      end if
c
c  equation for diffusion velocity
c
      do 53 j=1,idcomp
      do 52 i=1,nvar
   52 alam(j+1,i)=0
      alam(j+1,j+3)=amm*qx
c
   53 continue
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
      if(idgrhs.eq.1.and.istdpr.gt.0) write(istdpr,1250) (f(i),i=1,6)
      if(idgrhs.eq.-1.and.istdpr.gt.0) write(istdpr,1100) 
     *  'before return', (zz(i),i=1,5)
c
      if(iter.le.0.and..not.time0.and.idgn77.gt.0) then
        yffdsp(1,n)=x
        call store(y,yffdsp(2,n),nvar)
        call store(zz,yffdsp(2+nvar,n),nvar)
        call store(f,yffdsp(2+2*nvar,n),nvar)
        ii=1+3*nvar
        do j=1,nvar
          do i=1,nvar
            ii=ii+1
            yffdsp(ii,n)=fd(i,j)
          end do
        end do
	if(idgn77.ge.2) then
	  ii=1+nvar*(3+nvar)
	  do j=1,nvar
	    do i=1,idcomp+1
	      ii=ii+1
	      yffdsp(ii,n)=alam(i,j)
	    end do
	  end do
	  do k=1,nvar
            do j=1,nvar
              do i=1,idcomp+1
                ii=ii+1
                yffdsp(ii,n)=alamd(i,j,k)
              end do
            end do
          end do
	end if
	if(idgn77.ge.3) then
	  yffdsp(ii+1,n)=10.d0**ak
	  yffdsp(ii+2,n)=rho(1)
	  yffdsp(ii+3,n)=cp(1)
	  yffdsp(ii+4,n)=dad(1)
	  yffdsp(ii+5,n)=dr
	  yffdsp(ii+6,n)=dxh
	  yffdsp(ii+7,n)=dtxh
	  yffdsp(ii+8,n)=difscx(1)
	  ii=ii+8
	  do i=1,idcomp
	    ii=ii+1
	    yffdsp(ii,n)=difx(i,1)
          end do
	  do i=1,idcomp
	    ii=ii+1
	    yffdsp(ii,n)=velx(i,1,1)
          end do
	  do i=1,idcomp
	    ii=ii+1
	    yffdsp(ii,n)=velx(i,2,1)
          end do
        end if
	iyfdst=ii
        if(n.eq.nn) then
	  iter0=0
	  ialamw=idcomp+1
	  write(77) iter0,ialamw,iyfdst,nn
          write(77) ((yffdsp(i,k),i=1,iyfdst),k=1,nn)
          if(istdpr.gt.0) write(istdpr,*) 
     *      'RHS quantities written to d/s 77 in rhsd, ',
     *      'nn, nvar, iyfdst =', nn, nvar, iyfdst
        end if
      end if
      if(iter.eq.0) return
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
        fd(3,1)=dac*fd(2,1)+ddac(1)*f(2)
        fd(3,2)=dac*fd(2,2)+ddac(2)*f(2)
        fd(3,3)=dac*fd(2,3)+ddac(3)*f(2)
        fd(3,ilum)=  alfact*ddac(4)*f(2)
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
        fd(3,ilum)=alfact*(dtgrfc(4)/tgrfct+1)*f3
c
      end if
c
      f5=amm*f(ilum)
      fd(ilum,2)=eps(2)*f5
      fd(ilum,3)=eps(3)*f5
      fd(ilum,ilum)=-f5
c
c  end for time0
c
      if(time0) then
	if(idgn77.gt.0.and.iter.ge.0) then
          yffdst(1,n)=x
	  call store(y,yffdst(2,n),nvar)
	  call store(zz,yffdst(2+nvar,n),nvar)
	  call store(f,yffdst(2+2*nvar,n),nvar)
	  ii=1+3*nvar
          do j=1,nvar
	    do i=1,nvar
	      ii=ii+1
	      yffdst(ii,n)=fd(i,j)
            end do
          end do
	  if(n.eq.nn) then
	    iyfdst=1+nvar*(2+nvar)
	    ialamw=0
	    write(77) iter1,ialamw,iyfdst,nn
            write(77) ((yffdst(i,k),i=1,iyfdst),k=1,nn)
	  end if
	end if
        return
      end if
c
c  number of independent variables for energy generation
c
      iamx=nvar-3
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
c  Except in convective core, with imixc4 .ge. 1
c  (This is getting rather too messy and deserves a clean-up)
c
c  Note: for now also use ismdif = 3 to flag for keeping X derivatives
c  in convective envelope.
c
c..      if(convx.or.(inmixc.ne.1.and.nmxcor.gt.0.and.qx.le.qmxcor)) then
      if((imixc4.eq.0.or.idiffus.eq.0).and.
     *  (convx.or.(inmixc.ne.1.and.nmxcor.gt.0.and.qx.le.qmxcor)).and.
     *  (ismdif.ne.3.or.inc.gt.1)) then
c
c  reset dzdy(2,4)
c
	dzdy(2,4)=0
c
      else 
	if(idiffus.eq.2) then
	  if(conv) then
            fd(3,5)=ddac(6)*f(2)
	  else
	    fd(3,5)=akz*f3/(amm*y(5))
	  end if
        end if
        if(n.lt.nxhzer) then
          fd(1,4)=-rho(4)*f1
          fd(2,4)=-pt(4)*f2
          if(conv) then
            fd(3,4)=ddac(5)*f(2)+dac*fd(2,4)
          else
            fd(3,4)=akxa*f3
          end if
          fd(ilum,4)=eps(4)*f5
          if(iccomp.gt.0) then
            do 75 k=1,iccomp
   75       fd(ilum,4+2*idcomp+k)=eps(k+4)*f5
c
            do 76 ia=1,iamx
            ia1=ia+1
            if(ia.le.3) then
              i = ia+1
            else
              i = ia+1+2*idcomp
            end if
            do 76 k=1,iccomp
   76       fd(4+2*idcomp+k,i)=ft(ia1,k+1)
          end if
c
c  test for resetting for non-diffusing elements
c  in connection with mixed core
c
          if(inmixc.ne.1.and.nmxcor.gt.0.and.iccomp.gt.0) then
            if(qx.le.qmxcor) then
              do 78 k=1,iccomp
              fd(ilum,4+2*idcomp+k)=0
              do 78 i=1,nvar
   78         fd(4+2*idcomp+k,i)=0
            else if(qx.le.qmxcp) then
              do k=1,iccomp
                if(iche31.ne.1.or.k.ne.1) then
                  do i=1,nvar
                    fd(4+2*idcomp+k,i)=(1-tfrct)*fd(4+2*idcomp+k,i)
                  end do
                end if
              end do
            end if
          end if
        end if
c
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
          i=ia+1
          pta=pt(ia1)
          alamd(1,i,ilum)=-amm*alam(1,i)
          do 83 ib=1,ibmx
          ib1=ib+1
   83     alamd(1,i,ib1)=bb*(dlmb(ia,ib)+prhm*(pt(ib1)-rho(ib1))
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
c  derivatives in diffusion equations
c
      do 88 j=1,idcomp
c
      if(n.gt.1.and.n.lt.nn) then
	j1=4+idcomp+j
	if(j.eq.1) then
          fd(3+j,1)=amm*qx*(amms*amms/difx(j,1))*
     *      (-difx(j,2)/difx(j,1)*(y(j1)-(velx(j,1,1)/amms)*y(3+j)) -
     *        (velx(j,1,2)/amms)*y(3+j))
          fd(3+j,2)=amm*qx*(amms*amms/difx(j,1))*
     *      (-difx(j,3)/difx(j,1)*(y(j1)-(velx(j,1,1)/amms)*y(3+j)) -
     *        (velx(j,1,3)/amms)*y(3+j))
          fd(3+j,3)=amm*qx*(amms*amms/difx(j,1))*
     *      (-difx(j,4)/difx(j,1)*(y(j1)-(velx(j,1,1)/amms)*y(3+j)) -
     *        (velx(j,1,4)/amms)*y(3+j))
          fd(3+j,4)=amm*qx*(amms*amms/difx(j,1))*
     *      (-difx(j,6)/difx(j,1)*(y(j1)-(velx(j,1,1)/amms)*y(3+j)) -
     *        (velx(j,1,6)*y(3+j)+velx(j,1,1))/amms)
          fd(3+j,ilum)=alfact*amm*qx*(amms*amms/difx(j,1))*
     *      (-difx(j,5)/difx(j,1)*(y(j1)-(velx(j,1,1)/amms)*y(3+j)) -
     *        (velx(j,1,5)/amms)*y(3+j))
          fd(3+j,j1)=amm*qx*amms*amms/difx(j,1)
c
c  possibly set Z derivative
c
          if(idiffus.ge.2) fd(3+j,5)=amm*qx*(amms*amms/difx(j,1))*
     *      (-difx(j,7)/difx(j,1)*(y(j1)-(velx(j,1,1)/amms)*y(3+j)))
        else
          dvelj=y(j1)-
     *           (velx(j,1,1)/amms+velx(j,2,1)*y(5+idcomp))*y(3+j)
          fd(3+j,1)=amm*qx*(amms*amms/difx(j,1))*
     *      (-difx(j,2)/difx(j,1)*dvelj -
     *        (velx(j,1,2)/amms+velx(j,2,2)*y(5+idcomp))*y(3+j))
          fd(3+j,2)=amm*qx*(amms*amms/difx(j,1))*
     *      (-difx(j,3)/difx(j,1)*dvelj -
     *        (velx(j,1,3)/amms+velx(j,2,3)*y(5+idcomp))*y(3+j))
          fd(3+j,3)=amm*qx*(amms*amms/difx(j,1))*
     *      (-difx(j,4)/difx(j,1)*dvelj -
     *        (velx(j,1,4)/amms+velx(j,2,4)*y(5+idcomp))*y(3+j))
          fd(3+j,4)=amm*qx*(amms*amms/difx(j,1))*
     *      (-difx(j,6)/difx(j,1)*dvelj -
     *        (velx(j,1,6)/amms+velx(j,2,6)*y(5+idcomp))*y(3+j))
          fd(3+j,ilum)=alfact*amm*qx*(amms*amms/difx(j,1))*
     *      (-difx(j,5)/difx(j,1)*dvelj -
     *        (velx(j,1,5)/amms+velx(j,2,5)*y(5+idcomp))*y(3+j))
	  fd(3+j,3+j)=-amm*qx*(amms*amms/difx(j,1))*
     *           (velx(j,1,1)/amms+velx(j,2,1)*y(5+idcomp))
	  fd(3+j,5+idcomp)=-amm*qx*(amms*amms/difx(j,1))*
     *           velx(j,2,1)*y(3+j)
          fd(3+j,j1)=amm*qx*amms*amms/difx(j,1)
c
c  include Z derivative of difx
c
	  fd(3+j,5)=fd(3+j,5)+amm*qx*(amms*amms/difx(j,1))*
     *      (-difx(j,7)/difx(j,1)*dvelj)
c
        end if
      end if
c
   88 continue
c
c  For clarity, as a separate statement, zero hydrogen gradient
c  where X is exhausted
c
      if(n.ge.nxhzer) then
	do k=1,nvar
	  fd(4,k)=0.d0
        end do
      end if
c
c  derivatives of reaction part of hydrogen equation
c
      fd(5+idcomp,2)=-amm*qx*ft(2,1)
      fd(5+idcomp,3)=-amm*qx*ft(3,1)
      fd(5+idcomp,4)=-amm*qx*ft(4,1)
      if(iccomp.gt.0) then
	do 89 k=1,iccomp
   89   fd(5+idcomp,4+2*idcomp+k)=-amm*qx*ft(k+4,1)
      end if
c
c  test for resetting in connection with mixed core
c  (So far a slightly messy hack; needs revision)
c
      if(convx.or.(inmixc.ne.1.and.nmxcor.gt.0.and.qx.le.qmxcor)) then
	if(inmixc.ne.1.and.nmxcor.gt.0) then
	  if(qx.le.qmxcor) then
	    do i=1,nvar
c..	      fd(5+idcomp,i)=0
	    end do
          else if(qx.le.qmxcp) then
	    do i=1,nvar
c..              fd(5+idcomp,i)=(1-tfrct)*fd(5+idcomp,i)
	    end do
          end if
        end if
c
      end if
c
c  additionally reset Z-derivative near edge of convective core,
c  if wrong sign (we expect Z to decrease outwards!)
c
      iresf5=0
      if(iresf5.gt.0.and.
     *  nmxcor.gt.0.and.n.ge.nmxcor-3.and.f(5).gt.0) then
	if(istdpr.gt.0) write(istdpr,*) 'f(5) reset from ',f(5),
     *    ' to 0 at n =',n
	f(5)=0
	do i=1,nvar
	  fd(5,i)=0
        end do
      end if
c
c  note that alamd(2,.,.) =0
c
      if(mod(n-1,50).eq.-1.and.istdpr.gt.0) then
	write(istdpr,'(/a/i5,1p20e13.5)') 'n, x, y:',
     *    n, x, (y(i),i=1,nvar)
	write(istdpr,'(a/5e13.5)') 'D_1, V_1, D_2, V_21, V_22',
     *    difx(1,1),velx(1,1,1),difx(2,1), velx(2,1,1), velx(2,2,1)
	write(istdpr,'(a/i5,1p20e13.5)') 'n, x, f:',
     *    n, x, (f(i),i=1,nvar)
	write(istdpr,'(a)') 'i, fd:'
	do i=1,nvar
	  write(istdpr,'(i3,1p20e13.5)') i,(fd(i,j),j=1,nvar)
	end do
      end if
c
      if(idgn77.gt.0.and.iter.ge.0) then
        yffdst(1,n)=x
        call store(y,yffdst(2,n),nvar)
        call store(zz,yffdst(2+nvar,n),nvar)
        call store(f,yffdst(2+2*nvar,n),nvar)
        ii=1+3*nvar
        do j=1,nvar
          do i=1,nvar
            ii=ii+1
            yffdst(ii,n)=fd(i,j)
          end do
        end do
	if(idgn77.ge.2) then
	  ii=1+nvar*(3+nvar)
	  do j=1,nvar
	    do i=1,idcomp+1
	      ii=ii+1
	      yffdst(ii,n)=alam(i,j)
	    end do
	  end do
	  do k=1,nvar
            do j=1,nvar
              do i=1,idcomp+1
                ii=ii+1
                yffdst(ii,n)=alamd(i,j,k)
              end do
            end do
          end do
	end if
	if(idgn77.ge.3) then
	  yffdst(ii+1,n)=10.d0**ak
	  yffdst(ii+2,n)=rho(1)
	  yffdst(ii+3,n)=cp(1)
	  yffdst(ii+4,n)=dad(1)
	  yffdst(ii+5,n)=dr
	  yffdst(ii+6,n)=dxh
	  yffdst(ii+7,n)=dtxh
	  yffdst(ii+8,n)=difscx(1)
	  ii=ii+8
	  do i=1,idcomp
	    ii=ii+1
	    yffdst(ii,n)=difx(i,1)
          end do
	  do i=1,idcomp
	    ii=ii+1
	    yffdst(ii,n)=velx(i,1,1)
          end do
	  do i=1,idcomp
	    ii=ii+1
	    yffdst(ii,n)=velx(i,2,1)
          end do
        end if
	iyfdst=ii
        if(n.eq.nn) then
	  ialamw=idcomp+1
	  write(77) iter1,ialamw,iyfdst,nn
          write(77) ((yffdst(i,k),i=1,iyfdst),k=1,nn)
          if(istdpr.gt.0) write(istdpr,*) 
     *      'RHS quantities written to d/s 77 in rhsd, ',
     *      'nn, nvar, iyfdst =', nn, nvar, iyfdst
        end if
      end if
      return
 1100 format(' In rhs, ',a/(1p5e13.5))
 1120 format(/' s/r rhs at n =',i5,'  x =',1pe15.7/
     *  ' y:',6e15.7/' rhl, tl, xh, ak =',4e15.7)
 1130 format(' eps:',1p5e13.5)
 1150 format(' a3, ak, x, allred, tl, rlred, f(3) =',1p7e14.6)
 1160 format(/' ***** Warning in s/r rhs. Excessive no. of iteration =',
     *                i4/
     *        '       Rad. gradient reset from',1pe13.5,' to ',e13.5,
     *        ' at n =',i5)
 1165 format(/' ***** Warning in s/r rhs. Excessive no. of iteration =',
     *                i4/
     *        '       Semiconvection suppressed at n =',i5)
 1170 format(' f(2),f(3),dac =',1p3e15.7)
 1250 format(' f:',1p6e12.4)
 1300 format(///
     *  ' test rhs. ntime, iter, n, x, f(1-nvar), dtst(1-nvar):'/)
 1350 format(3i4,1p13e13.5)
 1400 format(' n, x, y(1-6) =',i5,1p7e15.7/'       z(1-6) =',
     *      20x,6e15.7/'       f(1-6) =',20x,6e15.7)
 1500 format(' fd:')
 1510 format(1p6e12.4)
      end
