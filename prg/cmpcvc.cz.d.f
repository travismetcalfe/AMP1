      subroutine cmpcvc(x,y,in,in1,iy,nn,compc,icomp,iche3,icnocs,
     *  iheccs,dt,it,iprdcr,irsmsh,icrycm)
c
c  handles predictor-corrector integration for composition in
c  convective core.
c
c  in and in1 are storage indices for previous and current solution
c
c  iprdcr = 1 for initial, predictor call
c  iprcdr = 2 for subsequent corrector calls
c
c  Original version: 8/5/92
c
c  Notes on number of composition variables:
c  ========================================
c
c  The variable icomp (passed as argument) is set to number of 
c  reacting species (ispxx3+ispcno). However, with diffusion 
c  mixing in a growing convective core has to include also heavy
c  elements. 
c
c  Thus introduce ifcomp set to the total number of composition 
c  variables (iccomp + idcomp) with diffusion. This is used where
c  nuclear reactions are not involved.
c  Note: it is assumed that the last element, for idiffus .ge. 2,
c  is the heavy-element abundance.
c
c               --------------------------------------
c
c  Modified 19/7/95, to use general heavy-element abundance
c  as set in common/heavy/. This has to be checked with care in
c  connection with convective mixing.
c
c  Modified 3/8/96, fixing additional occurrences of heavy-element
c  abundance (previously unset variable z was used).
c
c  Modified 26/8/96, testing for oscillating core, and
c  fixing boundary, for imixc1 = 1
c
c  Corrected 5/8/97, by including engenr.nnz.d.incl, to set nspdmx etc
c
c  Corrected 18/12/00, in freezing boundary.
c
c  Modified 21/12/00, including resetting of luminosity used for
c  determining edge of core, in s/r mixcor, from the luminosity
c  calculated in s/r rxhmxc.
c
c  Modified 21/12/00 to reset composition due to growing convective
c  core before call of rxhmxc (to obtain consistent luminosity,
c  principally). Also, change to base all tests on qmxcor, qmxcp,
c  as far as possible (NEEDS SOME FURTHER TESTS).
c
c  Modified 28/12/00, including various options for controlling new
c  iterations. These are flagged by the components (imixc0, imixc1, ...)
c  in input parameter imixcr.
c
c  Modified 3/1/01, to mix also heavy elements in case with diffusion,
c  by introducing ifcomp.
c
c  Modified 1/9/03, correcting storage of mixed elements in loop 58
c  to include also heavy elements, with diffusion.
c  NB: However, mixing of heavy elements requires more thought.
c
c  Modified 27/5/05, introducing resetting of log f (or log rho)
c  in case of growing convective core to ensure smooth variation
c  of pressure.
c
c  Modified 15/4/06, to only call resscn if the convective core
c  is bigger than a previous semiconvective region or if the
c  convective core is growing (NEEDS CAREFUL CHECKING!)
c
c  Modified 22/3/07, including contribution from possible diffusive flux
c  in change in convective-core hydrogen abundance.
c
c  Modified 18/1/08, including contribution from growing convective core
c  in crxmn, to be passed back to s/r rhs. On 22/1/08 modified to be
c  used only for imixc6 = 1,
c
c  Modified 22/1/08, resetting y(in1+1,.) throughout core where 
c  composition is reset, when imixc6 = 1.
c
c  *** Note: in old version z was unset before call of he3abd for
c      core He3 abundance.
c
      implicit double precision(a-h,o-z)
      include 'engenr.cz.d.incl'
      logical norct, nosd
c
      parameter(naztmx = nspcmx+3, 
     *  kaxst = 3+naztmx, nbmax = nspcmx+2)
      parameter(nspcm2=2*nspcmx,iymix=ivarmx)
c
      dimension x(*), y(iy,*), compc(*)
      dimension rxmn(nspcmx), drxmn(nspcmx), rxmnp(nspcmx),
     .  drxmnp(nspcmx),
     .  alrhmn(krnrmx), alrhmp(krnrmx), alrhmc(krnrmx), 
     .  dcmpgr(nspcmx),x3new(4),alres(nnmax), 
     .  ycrit(ivarmx),compcp(ivarmx),xrefine(10),
     .  qmxckt(nitmax),qmxnkt(nitmax),nmxckt(nitmax),frmxkt(nitmax),
     .  drmxkt(nitmax)
      common/rhcn/ a1,a2,a3,a4,zdummy,nvar,ifwrt,irhtst,inentc,
     *  ii1,ii2,ii3,icomp1
      common/compvr/ cvr(icvrmx,1)
      common/cstcvr/ icvh1, icvhe4, icvhe3, icvc12, icvc13, icvn14, 
     *  icvo16
      common/cstyst/ iyh1, iyhe3, iyc12, iyc13, iyn14, iyo16, 
     *  iyche4, iycc12
      common/cstazt/ iah1, iahe3, iac12, iac13, ian14, iao16, 
     *  iache4, iacc12
      common/eqstd/ xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp,gmm1(4)
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,jnc,frcf(6),frcl(6)
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor,nmxcp,qmxcp,rmxcp,
     *  qmxmin,qmxmax,qmxscn,qmxscp,nmxscn,nmxscp,icjump,iqcfrz
      common/cmxcit/ qmxcit(nitmax),rmxcit(nitmax),frmxit(nitmax),
     .  nmxcit(nitmax),xmxcit(nitmax),frmrit(nitmax),nmxrit(nitmax),
     .  drmxit(nitmax),qmxnit(nitmax),itcmpc
      common/crxcor/ cqc, cqcp, crxmn(nspcmx), crxmnp(nspcmx),
     *  crxxh_sc(nnmax)
      common/cymixc/ ymix(iymix,nnmax)
      common/cymixp/ ymixp(iymix,nnmax)
      common/ccmpor/ cmporg(nspcmx,nnmax)
      common/cmxcnt/ ddrmix, inmixc, imixc0, imixc1, imixc2, imixc3,
     *  imixc4, imixc5, imixc6, iosc_frz, imix_repeat, icnv_mixcor
      common/crxstr/ rxstr(nspcm2,1)
      common/thetac/ theta(ivarmx)
      common/rnratd/ al(10,krnrmx),norct
      common/heavy/ zatmos, zhc, zh(1)
      common/compos/ xseth, yset, xset3, xset12, xset13, xset14, xset16
      common/caddvr/ addvar(5,ngmax)
      common/ksider/ al01,al21,aztst(naztmx),axst(naztmx)
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec, idcomp, iccomp
      common/cdiffu/ idiffus
      common/cntmsh/ wx,epsr,wr,wpt,wxh,wx3,wdgr,
     *  wmshos,drmsos,wmshcc,dqlmsc,qmscmn,qmscmx,
     *  koint,kvint,iprt,icngr,nnt,iwdgrd,istrtc
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3
      common/cdgphs/ kdgeos, kdgopc, kdgeng
      common/noiter/ iter, ntime, epscon, eamcon, it_force
c
c  dataset designators
c
      common/cdsdsg/ idstrl,idshvz,idszab,idsevl,idssum,idscen,idsnov,
     *  idsoov,idsgng,idsgsm,idstm1,idstm2,iducen,iddgm1
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      data idg_versb, idg_revc6 /0, 0/
c
      save
c
      if(istdpr.gt.0) write(istdpr,'(/'' Entering cmpcvc''/)')
c
c  hardcode imixc6 to 1 as a fairly desparate attempt
c
      imixc6=1
c
c  set flag for reducing to previous version (vers. b)
c
      iset_versb=0
      if(iset_versb.gt.0.and.idg_versb.eq.0) then
	idg_versb=1
	write(istdou,'(/'' Reset to version b.''/
     *   '' Switch off resetting of y(2,.) in predictor step''/
     *   )')
      end if
c
      kdgeos=1
c
c  unset flag to force iteration (may be set below, if convection-zone
c  boundary is fixed to value at previous timestep)
c
      it_force=0
c
c  for now, hardcode imixc5
c
c..      imixc5=1
c
c  set iteration counts for freezing boundaries (fixed at 10 and 12 
c  before 27/8/03)
c
c  changed from 15 and 18 29/5/05
c
      if(imixc5.eq.1.or.imixc6.eq.1) then
        itfrz1=100
        itfrz2=100
      else
        itfrz1=12
        itfrz2=15
      end if
c
c  set full number of composition variables, with diffusion
c
      if(idiffus.gt.0) then
	ifcomp = iccomp+idcomp
      else
	ifcomp = icomp
      end if
c
c  maximum number of iterations for composition in this routine
c  depending in imixc3 (note: before 21/12/00 hardcoded to 0)
c
      if(imixc3.eq.0) then
        nitcmx = 0
      else
        nitcmx = 20
      end if
c
c  set convergence limit
c
      epscmx = 1.e-5
c
      nosd=.true.
c
      icrycm=0
      irsmsh=0
c
c  test for predictor or corrector
c
      if(iprdcr.eq.1) then
c
        if(istdpr.gt.0) 
     *    write(istdpr,'(/'' Start predictor step in cmpcvc''/)')
        if(frmxc.gt.1.and.nmxcor.gt.0) then
          nmxcp=nmxcor+1
        else
          nmxcp=nmxcor
        end if
        qmxcp=qmxcor
        rmxcp=rmxcor
	if(nmxcor.gt.0) then
	  qmxcor_prev=qmxcp
	  dxcrit=x(nmxcor-1)-x(nmxcor)
        else
	  qmxcor_prev=1.d0
	end if
	nn_old=nn
        if(istdpr.gt.0) write(istdpr,*) 'nmxcp, qmxcp =',nmxcp, qmxcp
	nmxscp=nmxscn
	qmxscp=qmxscn
c
	imix_repeat=0
	icnv_mixcor=1
c
c  test for nmxcor = 0. There seems to be no good reason to continue
c  in that case, while the preceding setups might (just) serve a 
c  purpose.
c  Added 3/10/06
c
	if(nmxcor.le.0) then
	  write(istdou,'(/'' nmxcor ='',i5,'' Exit cmpcvc''/)') nmxcor
	  if(istdpr.gt.0.and.istdou.ne.istdpr)
     *      write(istdpr,'(/'' nmxcor ='',i5,'' Exit cmpcvc''/)') nmxcor
c
          return
        end if
c
c  test for very near hydrogen exhaustion in core
c  currently hard-code test to be for X .le. 1.e-4
c  after previous time step, or X .le. 1.e-3 and core
c  mass less than 0.02M. All this, of course, only if
c  there is no 4He burning. Base test on He burning on 
c  crude temperature limit, as in s/r rnrhec.
c  In this case switch off core mixing
c
        if((iheccs.eq.0.or.y(in+2,nn).le.7.778).and.
     *    (y(in+4,nn).lt.1.e-4.or.
     *    y(in+4,nn).lt.1.e-3.and.qmxcp.lt.0.02)) then
          if(istdpr.gt.0) write(istdpr,105) qmxcp,y(in+4,nn)
          nmxcor=-1
          qmxcor=0.
          if(istdpr.gt.0) write(istdpr,'(/'' Exiting cmpcvc''/)')
          return
        else if(iheccs.ne.0.and.y(in+2,nn).gt.7.778.and.
     *    (y(in+iyche4-1,nn).lt.1.e-4.or.
     *    y(in+iyche4-1,nn).lt.1.e-3.and.qmxcp.lt.0.02)) then
          if(istdpr.gt.0) write(istdpr,106) qmxcp,y(in+iyche4-1,nn)
          nmxcor=-1
          qmxcor=0.
          if(istdpr.gt.0) write(istdpr,'(/'' Exiting cmpcvc''/)')
          return
        end if
c
c  set average rx for previous timestep
c
        if(idgeng.eq.-3.and.istdpr.gt.0) write(istdpr,*) 
     *    'before calling rxhmxc thte =',thte
        if(istdpr.gt.0) write(istdpr,*) 
     *    'calling rxhmxc with in, nn, icomp, y(in+icomp,nn) =',
     *    in, nn, icomp, (y(in+3+i,nn),i=1,icomp)
        call rxhmxc(x,y(in,1),iy,nn,y(in+4,nn),qcp,rxmnp,drxmnp,icomp,
     *    alres,iprdcr,icryrx)
	if(icryrx.lt.0) then
	  icrycm=-1
	  go to 90
	end if
        if(idgeng.eq.-3.and.istdpr.gt.0) write(istdpr,*) 
     *    'after calling rxhmxc thte =',thte
c
c  store values at previous time step in common
c
        cqcp=qcp
        call store(rxmnp,crxmnp,icomp)
c
c  set average reaction rates for computing He3 abundance
c
        call almmxc(x,y(in,1),iy,nn,y(in+4,nn),alrhmp,6,1)
c
        if(it.le.0) then
c
c  if this call is before start of iteration loop, set
c  predicted abundances.
c
c  set predicted new core abundance of He3
c
          xh=y(in+4,nn)
          yh=1-xh-zhc
          idghe3=1
          call he3abc(xh,yh,zhc,dt,y(in+5,nn),alrhmp,x3new,x3eq,anu,
     *      1,nosd)
          idghe3=0
c
c  reset predicted new core abundance of hydrogen and CNO elements
c  (note: skip He3)
c
c#ai#  Note: In old version of programme x3new was never set in
c      y at this point. Hence it appears that He3 abundance was not 
c      passed on in convective core. The effects of this needs badly
c      looking into, as does the whole question of He3 in
c      a convective core. For the moment (11/8/91) keep
c      old treatment for compatibility.
c
c      Actually it seems that He3 is set in later loop (50)
c      during iteration loop.
c
	  if(nmxscn.gt.0) then
	    nc1=nmxscn+1
          else if(frmxc.gt.1) then
            nc1=nmxcor+1
          else
            nc1=nmxcor
          end if
	  if(istdpr.gt.0) write(istdpr,*) 
     *      'In cmpcvc, predictor step, set composition with nc1 =',
     *      nc1
	  do i=1,icomp
            if(i.eq.2.and.iche3.eq.1) then
              compc(i)=x3new(1)
            else
              compc(i)=max(y(in+3+i,nn)+dt*rxmnp(i),1.d-10)
	    end if
            if(istdpr.gt.0) write(istdpr,*) 
     *        'predicted compc(i) set to',compc(i)
	  end do
c
          do n=nc1,nn
c
c  prepare for resetting log f (log rho) with new composition
c
	    xhh=y(in1+4,n)
	    yhh=1.d0-xhh-zh(n)
            nosd=.true.
	    fl=y(in1+1,n)
	    tl=y(in1+2,n)
	    call eqstf(fl,tl,xhh,yhh,zh(n),nosd,nosd)
            if(kdgeos.lt.0) then
	      icrycm=-1
	      return
            end if
	    pl=log10(pt(1))
c
            do i=1,icomp
              y(in1+3+i,n)=compc(i)
	    end do
c
c  reset log f (log rho) to keep log p fixed
c
	    xhh=compc(1)
	    yhh=1.d0-xhh-zh(n)
	    call eqstp(pl,tl,xhh,yhh,zh(n),nosd,nosd,fl,nit)
	    if(iset_versb.eq.0) y(in1+1,n)=fl
c
	  end do
c
c  remove possible undesirable effect of extrapolation
c#ai# needs a check; modified to account for 4He burning
c
          do 27 n=1,nmxcor
          if(2.e10.lt.compc(1).and.y(in1+4,n).lt.compc(1)) then
            do 26 i=1,icomp
   26       y(in1+3+i,n)=compc(i)
          end if
   27     continue
c
c#ai#  As a temporary fix, store mixed H and He abundances in cvr
c#ai#  This has to be cleared up with care.
c
c#ai# Also note that here mixed He3 abundance must be stored in cvr
c     later.
c
          do 30 n=nc1,nn
          cvr(1,n)=compc(i)
          cvr(2,n)=1-compc(i)-zhc
c
c  store CNO abundances, preferentially from values set
c  in compc. 
c
	  call setcno(y(in1,n),n)
c
c  test for storing 4He burning abundances
c
	  if(iheccs.ne.0) then
	    cvr(icvhe4,n)=y(in1+iyche4-1,n)
	    cvr(icvc12,n)=y(in1+iycc12-1,n)
          end if
c
   30     continue
c
c..        write(38,'(4i5,1p2e13.5)') 
c..     *    (ntime,iter,-1,n,x(n),y(in+4,n),n=nmxcor-10,nmxcor+10)
c
        end if
c
c  initialize counter for internal storage of core parameters
c  (used when s/r cmpcvc is only called at selected iteration
c  steps)
c
	itcmp=0
	itcmpc=itcmp
c
c  initialize flag for allowing call of resscn
c
	no_resscn=0
c
c  reset log f (log rho) in region of increasing core, to
c  ensure smooth pressure, for both previous and present timestep
c
        call crslgf(x,y(in,1),iy,nn,qmxcor,qmxcp,itcmp,
     *    ncp1,ncp2,plgp1,plgp2,icrycm)
	if(icrycm.lt.0) go to 90
c
        call crslgf(x,y(in1,1),iy,nn,qmxcor,qmxcp,itcmp,
     *    ncp1,ncp2,plgp1,plgp2,icrycm)
	if(icrycm.lt.0) go to 90
c
c  end of predictor step
c
c  --------------------------------------------------------------------
c
      else
c
c  limit for looking for matching core size
c
        dqrmin_0=2.d-2
c
c  corrector step during iteration
c
        if(istdpr.gt.0) 
     *    write(istdpr,'(/'' Start corrector step in cmpcvc''/)')
c
c  For major revision, set flag. Implemented for growing convective
c  core.
c
	if(imixc6.eq.1.and.igrw_core.eq.1) then
	  irevc6=1
	  if(istdpr.gt.0) write(istdpr,'(/
     *      '' ***** Use revised (8/2/08) formulation''/)')
	else
	  irevc6=0
        end if
c
c  Temporarily switch off
c
	if(irevc6.eq.1) then
          if(idg_revc6.eq.0) then
	    write(istdou,*) '*** Switch off mesh refinement'
	    idg_revc6=1
	  end if
	  irevc6=0
        end if
	if(ntime.ge.2) then
          if(idg_revc6.lt.2) then
            write(istdou,*) '*** Switch on new localization'
            write(istdou,*) '    dqrmin_0 =',dqrmin_0
            idg_revc6=2
	  end if
	  irevc6=2
        end if
c..c
c..c  switch off irevc6 (whatever that is!)
c..c
c..        irevc6=0
c..        write(istdpr,*) '#D1 irevc6 set to 0'
c
c  for initial iteration initialize limits for tests for oscillating 
c  core. Also reset iextrp in case it has been set to 9 at a 
c  previous timestep.
c  (So far applied only for imixc6 = 1)
c
        if(itcmp.le.1.and.imixc6.eq.1) then
          itoscl=4
          itoscs=8
	  iosc_frz=0
	  iosc_int=0
	  iextrp=0
	  igrw_core=0
        end if
c
	itcmp=itcmp+1
	itcmpc=itcmp
	if(istdpr.gt.0) write(istdpr,
     *    '(/'' Entering s/r cmpcvc with it, itcmp ='',2i3)') it, itcmp
c
c  store qmxcor_tfr (to be used in setting intermediate points) for
c  previous iteration, for possibly undercorrection
c
        qmxcor_tfrp=qmxcor_tfr
c
c  set qmxcor_itp and nmxcor_itp to test for oscillating boundary
c
	qmxcor_itp=qmxcor
	nmxcor_itp=nmxcor
c
c  store original composition, for resetting outside convective region
c
        do 32 i=1,ifcomp
        do 32 n=1,nn
   32   cmporg(i,n)=y(in1+3+i,n)
c
c  set original integral of abundance
c
        xhintp=0.
        call tstxin(x,y(in1,1),nn,iy,1,xhint,xhintp,xtime,xtime,icry,
     *    'initialize')
        xhintp=xhint
c
c  set model, including composition, for setting of boundary
c  of mixed region
c
c  Also store previous ymix for possible undercorrection
c
        do 34 n=1,nn
	do i=1,4+icomp
	  ymixp(i,n)=ymix(i,n)
        end do
        do 33 i=1,4
   33   ymix(i,n)=y(in1+i-1,n)
c
c  only reset ymix in radiative zone or for first iteration
c  Otherwise it is assumed that the appropriate value of ymix is
c  already in place from previous call of cmpcvc (change 23/11/05)
c
	if(itcmp.eq.1.or.n.lt.nmxcor) then
          do i=1,icomp
	    ymix(4+i,n)=y(in1+3+i,n)
	  end do
        end if
c
c  Suppress this questionable resetting of composition
c  (Might be appropriate for first iteration, at most).
c
c..   34   ymix(4+i,n)=max(y(in+3+i,n)+dt*rxstr(icomp+i,n),1.d-10)
   34   continue
c
c  Set type of setting unstable boundary.
c  Modified 28/12/00 for more clarity, but consistent with
c  older version. Might benefit from a little further thought.
c
c  Add forcing of unstable point (but still with extrapolation) 
c  when there has been a jump in position of core, as flagged
c  by icjump, and imixc5 = 1
c
        if(itcmp.lt.itfrz1) then
c
c  set standard case, depending on imixc2 and imixc5
c  Freeze iextrp = 9 if already set.
c
	  if(imixc5.eq.1.and.icjump.eq.1) then
	    iextrp = 3
          else if(imixc2.eq.0) then
            iextrp = 3
          else if(imixc6.ne.1.or.iextrp.ne.9) then
            iextrp = 1
          end if
c
c  Switch off forcing zero crossing. Not clear that this is
c  important
c
c..c
c..c  #Rev.# Try iextrp = 3 at timestep no. 90.
c..c  
          if(irevc6.eq.2) iextrp=3
c
        else if(iextrp.ne.-1.and.imixc1.eq.1) then
c
c  Test for freezing boundary, or region where boundary is set
c  This is based on behaviour at previous two iterations
c  Choose outermost point, and arrange such that extrapolation
c  is used, from points well removed from boundary
c
c  new type of freeze. Use mean of past two locations and
c  set interpolation weight to unity.
c
c  test whether one or two of the last locations had no convective core
c
c  Only apply when imixc6 = 0
c
	  if(imixc6.ne.1) then
            if(nmxcit(itcmp-1).eq.0.and.nmxcit(itcmp-2).eq.0) then
              nmxcor=-1
            else if(nmxcit(itcmp-1).eq.0.or.nmxcit(itcmp-2).eq.0) then
              nmxcor=max(nmxcit(itcmp-1),nmxcit(itcmp-2))
            else
              nmxcor=0.5*(nmxcit(itcmp-1)+nmxcit(itcmp-2))
            end if
          end if
c
          if(nmxcor.eq.-1) then
            qmxcor=0.
            rmxcor=0.
            frmxc=1.d0
            if(istdpr.gt.0) write(istdpr,102)
            iextrp=-1
          else 
            qmxcor=10.d0**x(nmxcor)
            rmxcor=10.d0**y(in1,nmxcor)/10.d0**y(in1,1)
            frmxc=1.d0
            if(istdpr.gt.0) write(istdpr,103) 
     *        nmxcor, qmxcor, rmxcor, frmxc
            iextrp=-1
          end if
        else if(itcmp.ge.itfrz2) then
c
c  set iextrp = -1, to freeze boundary completely
c
          iextrp=-1
c
c  for it = itfrz2, set frozen interpolation weight
c
          if(itcmp.eq.itfrz2) then
	    frmxc=min(frmxit(itcmp-2),frmxit(itcmp-1))
            if(frmxc.le.0) frmxc=0
            if(istdpr.gt.0) write(istdpr,*) 
     *        'Freeze interpolation weight at frmxc =',frmxc
          else
            if(istdpr.gt.0) write(istdpr,*) 
     *        'Interpolation weight frozen at frmxc =',frmxc
	  end if
        else if(it.ge.itfrz1) then
          if(iextrp.eq.3) then
            if(nmxcit(itcmp-1).lt.nmxcit(itcmp-2)) then
              itfrz=itcmp-1
            else
              itfrz=itcmp-2
            end if
            frmfrz=frmxit(itfrz)
            if(frmfrz.le.1.2) then
              nmxfrz=nmxcit(itfrz)-1
            else
              nmxfrz=nmxcit(itfrz)
            end if
            if(istdpr.gt.0) write(istdpr,*)
     *        'Freeze boundary of mixed region to be near n =',  nmxfrz
            iextrp=4
          end if
        end if
c
c  Force freeze for rapidly growing core, again in desperation
c
        if(itcmp.ge.15.and.qmxcor.ge.2*qmxcp) then
          if(istdpr.gt.0) write(istdpr,*) 
     *      '#D# Force freeze with growing core'
          iextrp=-1
        end if
c
c  test for basing test for mixed region on non-zero
c  gradient difference
c
        if(ddrmix.gt.0.and.iextrp.eq.3) iextrp=10
c
c  when setting semiconvective region, force extrapolation
c  before jump in size of core
c
	if(imixc5.eq.1.and.icjump.eq.0) iextrp=1
c
c  start of iteration loop for composition.
c  currently do simple back substitution, and test for convergence
c  on X(H) only
c
        xcitp=compc(1)
        nitcmp=0
c
c  Entry point for continuing iteration
c
   35   nitcmp=nitcmp+1
c
c  set extent of mixed core
c  Note that we normally set iextrp = 3, and hence use 
c  interpolation with point in the unstable region.
c  However, when freezing region of search, set iextrp = 4
c
        qmxcrp = qmxcrt
        qmxcrt = qmxcor
	frmxcp=frmxc
	nmxcorp=nmxcor
c
c  Here impose new type of varying core, in test case
c
	if(imixc6.eq.1.and.ntime.le.-37.and.itcmp.gt.5) then
          iosc_frz=1
	  if(itcmp.eq.6) then
	    nmx=nmxrit(itcmp-1)
	    dqmsh=10.d0**x(nmx-1)-10.d0**x(nmx)
	    itcmp_frz=itcmp
	    imix_repeat=1
	    icnv_mixcor=0
	  else if(imix_repeat.eq.7) then
c
	    if(istdpr.gt.0.and.icnv_mixcor.ne.1) write(istdpr,
     *        '(/ '' Converged solutions:''/
     *        '' irepeat, nmxcor, qmxcor, frmxc, ddrmax, qmxcor_new,'',
     *        '' qmxcor_new - qmxcor:''/
     *        (2i5, 1pe12.5,0pf8.3,1pe11.3,e12.5,e11.3))')
     *        (i, nmxckt(i), qmxckt(i), frmxkt(i), 
     *         drmxkt(i),qmxnkt(i),qmxnkt(i)-qmxckt(i),
     *         i=1,imix_repeat-1)
	    icnv_mixcor=1
	  else if(mod(itcmp-itcmp_frz,3).eq.0) then
	    qmxckt(imix_repeat)=qmxcor
	    qmxnkt(imix_repeat)=qmxcor_new
            nmxckt(imix_repeat)=nmxcit(itcmp-1)
            frmxkt(imix_repeat)=frmxit(itcmp-1)
            drmxkt(imix_repeat)=drmxit(itcmp-1)
c
	    imix_repeat=imix_repeat+1
c..	    dqmxcor=(-1+2*mod(imix_repeat-1,2))*0.2d0*imix_repeat*dqmsh
	    dqmxcor=-0.1*dqmsh
	    if(istdpr.gt.0) 
     *        write(istdpr,'(//'' Shift qmxcor by'',1pe13.5)') dqmxcor
	    qmxcor=qmxcor+dqmxcor
          end if
        end if
c
	i=-nn
	if(i.gt.0.and.istdpr.gt.0) write(istdpr,*) 
     *    '#D# 1. i =',i,'  y(in1+4,i), ymix(5,i) =', 
     *    y(in1+4,i), ymix(5,i)
c..        if(iextrp.ne.-1.and.(itcmp.eq.1.or.nxcrit.eq.0).and.
c..     *    (imixc6.ne.1.or.(iosc_frz.ne.1.and.(ntime.ne.37.or.
c..     *     itcmp.le.5.or.mod(itcmp,2).eq.0)))) then
        if(iextrp.ne.-1.and.(itcmp.eq.1.or.nxcrit.eq.0).and.
     *    (imixc6.ne.1.or.iosc_frz.ne.1)) then 
          if(istdpr.gt.0) write(istdpr,
     *       '(/'' Call mixcor from cmpcvc, iextrp ='',i3/)') iextrp
	  if(imixc5.ne.1) then
            call mixcor(x,ymix,iymix,nn,compc,iextrp,nmxfrz,itcmp)
	    xmxcor=ymix(5,nn)
          else
            call mixcor(x,y(in1,1),iy,nn,compc,iextrp,nmxfrz,itcmp)
	    xmxcor=y(in1+4,nn)
c
c  possibly apply undercorrection
c
	    if(itcmp.gt.10.and.nmxcorp.eq.nmxcor) then
	      ucy=0.4d0
	      frmxcu=(1-ucy)*frmxcp + ucy*frmxc
	      qmxcoru=10.d0**(frmxcu*x(nmxcor)+(1-frmxcu)*x(nmxcor-1))
	      if(istdpr.gt.0) write(istdpr,'(/
     *          '' ***** Warning in s/r cmpcvc. frmxc reset from '',
     *          f11.7,'' to'',f11.7/
     *          ''       qmxcor reset from'',1pe15.7,'' to'',
     *          e15.7/)') frmxc, frmxcu, qmxcor, qmxcoru
	      qmxcor=qmxcoru
	      frmxc=frmxcu
            end if
	  end if
	else if(imixc6.eq.1.and.iosc_frz.eq.1.and.istdpr.gt.0) then
	  write(istdpr,'(/'' Frozen core location, qmxcor ='',
     *      1pe13.5/)') qmxcor
        end if
c
c  For irevc6 = 2, apply new localization of edge of core 
c
c..	if(irevc6.eq.2.and.ntime.ge.37.and.itcmp.ge.3) then
c..          call loc_qmxcor_n(x,y(in1,1),iy,nn,itcmp,iosc_frz,ilocry,
c..     *      qmxcor_new)
c..	else if(irevc6.eq.2.and.iosc_frz.ne.1) then
c..          call loc_qmxcor(x,y(in1,1),iy,nn,itcmp,iosc_frz,ilocry)
c..        end if
c
        frmxit(itcmp)=frmxc
        nmxcit(itcmp)=nmxcor
        qmxcit(itcmp)=qmxcor
        rmxcit(itcmp)=rmxcor
        xmxcit(itcmp)=xmxcor
c
c  test for potentially growing core
c
	if(imixc6.eq.1.and.qmxcor.gt.qmxcp) igrw_core=1
c
c  For imixc6 = 1, set closest meshpoint and set interpolation 
c  parameters for closest meshpoint
c
	if(imixc6.eq.1) then
	  if(iosc_frz.lt.1) then
	    do n=nmxcor-1,nn
	      q1=10.d0**x(n)
	      if(qmxcor.ge.q1) then
	        q2=10.d0**x(n-1)
	        frq12=(q2-qmxcor)/(q2-q1)
	        frmrit(itcmp)=frq12
	        nmxrit(itcmp)=n
                go to 36
              end if
	    end do
	  else
	    frmrit(itcmp)=frmrit(itcmp-1)
	    nmxrit(itcmp)=nmxrit(itcmp-1)
c
c  set drmxit corresponding to this point (which in many cases should 
c  be local maximum). For iosc_frz .ne.1 this should have been set
c  in loc_qmxcor (as things are now)
c
c  18/5/09: Switch off special treatment at ntime = 37
c
c..	   if(ntime.lt.37) then
	   if(ntime.gt.-37) then
	     nmx=nmxrit(itcmp)
	     drad=fdradp(x(nmx),y(in1,nmx),addvar(1,nmx),zh(nmx),
     *         flmx,akmx,akrmx,aktmx,akxmx)
	     drmxit(itcmp)=drad-dad(1)
           end if
c
          end if
        end if
c
   36   continue
c
	if(istdpr.gt.0) write(istdpr,
     *    '(/ '' it, nmxcor, qmxcor, frmxc, Xc, nmx_res, frmx_res,'',
     *    '' ddrmax, qmxcor_new:''/
     *    (2i5, 1pe12.5,0pf8.3,f12.8,i5,f8.3,1pe11.3,e12.5))')
     *    (i, nmxcit(i), qmxcit(i), frmxit(i), xmxcit(i), 
     *    nmxrit(i), frmrit(i), drmxit(i),qmxnit(i),i=1,itcmp)
c
c  For imixc6 = 1 check for large oscillations in core, to fix near the 
c  outermost choice
c
c  Switch off for irevc6 = 2
c
	if(itcmp.gt.itoscl.and.imixc6.eq.1.and.irevc6.ne.2.and.
     *    nmxcit(itcmp).ge.nmxcit(itcmp-1)+2) then
c
c  step back to test for repeated location of core
c
	  iosc=0
	  do i=itcmp-2,itoscl-2,-1
	    if(nmxcit(i).eq.nmxcit(itcmp).and.iosc.eq.0) iosc=i
          end do
c
c  test whether this is the first case of oscillating core
c
	  if(iosc.gt.0.and.iosc_int.eq.0) then
c
            iosc_int=1
c
c  repeat call of mixcor
c
            if(istdpr.gt.0) write(istdpr,'(/'' Oscillating core.'',
     *        '' nmxcit('',i2,'') = nmxcit('',i2,'') ='',i5/
     *        '' Repeat call of mixcor with iextrp = 9''/)')
     *        iosc, itcmp, nmxcit(iosc)
            iextrp=9
	    iqcfrz=1
            call mixcor(x,ymix,iymix,nn,compc,iextrp,nmxfrz,itcmp)
c
            frmxit(itcmp)=frmxc
            nmxcit(itcmp)=nmxcor
            qmxcit(itcmp)=qmxcor
            rmxcit(itcmp)=rmxcor
c
c  set closest meshpoint and set interpolation 
c  parameters for closest meshpoint
c
            do n=nmxcor-1,nn
              q1=10.d0**x(n)
              if(qmxcor.ge.q1) then
                q2=10.d0**x(n-1)
                frq12=(q2-qmxcor)/(q2-q1)
                frmrit(itcmp)=frq12
                nmxrit(itcmp)=n
                go to 37
              end if
            end do
c
   37       continue
c
            if(istdpr.gt.0) write(istdpr,
     *        '(/ '' it, nmxcor, qmxcor, frmxc, Xc, nmx_res,'',
     *        '' frmx_res, ddrmax:''
     *        /(2i5, 1pe13.5,0pf10.5,f12.8,i5,f10.5,1pe13.5))')
     *        (i, nmxcit(i), qmxcit(i), frmxit(i), xmxcit(i), 
     *        nmxrit(i), frmrit(i), drmxit(i),i=1,itcmp)
c
	  else if(iosc.gt.0) then
c
c  Test for nearly disappearing core
c
	    if(nmxcit(itcmp-1).eq.0) then
	      qmxcor=qmxcit(itcmp)
	      nmxcor=nmxcit(itcmp)
	      if(istdpr.gt.0) write(istdpr,'(/
     *          '' Freeze almost vanishing core'')')
	      iqcfrz=14
c
	    else
c
c  Repeated oscillating core. Use freezing of core at maximum in
c  nabla_R - nabla_ad (as set in addvar)
c  Exclude final point.
c
              iosc_int=iosc_int+1
              iosc_frz=1
	      iqcfrz=13
              n1=nmxcit(itcmp-1)-2
              n2=nmxcit(itcmp)
              nmxddr=0
              ddrmax=-1.d10
              do n=n1,n2
                if(addvar(4,n).gt.ddrmax.and.
     *            addvar(4,n).gt.addvar(4,n+1).and.
     *            addvar(4,n+1).gt.-1.d5) then
                  nmxddr=n
                  ddrmax=addvar(4,n)
                end if
              end do
              if(nmxddr.eq.0) then
                write(istder,'(/'' **** Error in cmpcvc.'',
     *            '' Failed to find maximum in nabla_R - nabla_ad'')')
                write(istdou,'(/'' **** Error in cmpcvc.'',
     *            '' Failed to find maximum in nabla_R - nabla_ad'')')
                if(istdpr.gt.0.and.istdou.ne.istdpr) 
     *            write(istdpr,'(/'' **** Error in cmpcvc.'',
     *            '' Failed to find maximum in nabla_R - nabla_ad'')')
                stop 'cmpcvc'
              end if
              qmxcor=10.d0**x(nmxddr)
	      frmxc=1.d0
              rmxcor=10.d0**y(in1,nmxcor)/10.d0**y(in1,1)
              frmxit(itcmp)=frmxc
              nmxcit(itcmp)=nmxcor
              qmxcit(itcmp)=qmxcor
              rmxcit(itcmp)=rmxcor
              frmrit(itcmp)=frmxc
              nmxrit(itcmp)=nmxcor
c
              if(istdpr.gt.0) write(istdpr,'(/ '' Freeze qmmxcor ='',
     *           1pe13.5,'' at maximum in nabla_R - nabla_ad'')') qmxcor
              nmxcor=nmxddr
c
	    end if
c
            iosc_frz=1
c
c  reset limit for test for small oscillating core
c
            itoscs=max(itoscs,itcmp+1)
            if(istdpr.gt.0) write(istdpr,
     *        '(/ '' it, nmxcor, qmxcor, frmxc, Xc, nmx_res,'',
     *        '' frmx_res, ddrmax:''
     *        /(2i5, 1pe13.5,0pf10.5,f12.8,i5,f10.5,1pe13.5))')
     *        (i, nmxcit(i), qmxcit(i), frmxit(i), xmxcit(i), 
     *        nmxrit(i), frmrit(i), drmxit(i),i=1,itcmp)
          end if
        end if
c
c  For imixc6 = 1 check for small oscillations in core, freeze at
c  average location.
c  For itcmp well above itoscs allow larger step, as a rather 
c  desparate measure.
c
	iset_osc=0
	if(imixc6.eq.1.and.itcmp.gt.itoscs) then
c
c  Switch off for irevc6 = 2
c
          if(irevc6.ne.2.and.((nmxrit(itcmp).eq.nmxrit(itcmp-1)-1).or.
     *      (itcmp.gt.itoscs+4.and.nmxrit(itcmp).le.nmxrit(itcmp-1)-1)))
     *      then
c
c  Step back to test for repeated location of core
c
	    iosc=0
	    do i=itcmp-2,itoscs-2,-1
	      if(nmxrit(i).eq.nmxrit(itcmp).and.iosc.eq.0) iosc=i
            end do
            if(iosc.gt.0) then
	      iset_osc=1
	      iqcfrz=11
              if(istdpr.gt.0) write(istdpr,
     *          '(/'' Small oscillation in qmxcor.'',
     *          '' nmxrit('',i2,'') = nmxrit('',i2,'') ='',i5)')
     *          iosc, itcmp, nmxrit(iosc)
	    end if
          else if(qmxcit(itcmp)/qmxcit(itcmp-1)-1.d0.gt.1.d-5)
     *      then
c
c  No significant change in mesh point number. Check for smaller
c  oscillation in location
c  Step back to test for repeated location of core
c
	    iosc=0
            if(itcmp.le.15) then
              dqrmin=dqrmin_0
            else
              dqrmin=1.4*dqrmin
            end if
            if(istdpr.gt.0) write(istdpr,*) '#D1 test with dqrmin =',
     *        dqrmin
	    dqcmp=max(1.d-2*qmxcit(itcmp),qmxcit(itcmp)-qmxcit(itcmp-1))
	    do i=itcmp-2,itoscs-2,-1
	      if(abs((qmxcit(itcmp)-qmxcit(i))/dqcmp).le.dqrmin.and.
     *          iosc.eq.0) iosc=i
            end do
            if(iosc.gt.0) then
	      iset_osc=1
	      iqcfrz=12
              if(istdpr.gt.0) write(istdpr,
     *          '(/'' Small oscillation in qmxcor.'',
     *          '' qmxcit('',i2,'')= '',1pe13.6,
     *          '' = qmxcit('',i2,'')='',e13.6)')
     *          iosc, qmxcit(iosc), itcmp, qmxcit(itcmp)
	    end if
c
          end if  
c
c  only freeze and reset if matching point is in fact found
c
          if(iset_osc.gt.0) then
c
            qmxcor=0.5d0*(qmxcit(itcmp)+qmxcit(iosc))
c  
            nmxcor=nmxrit(itcmp)
            q1=10.d0**x(nmxcor)
            if(qmxcor.lt.q1) then
              nmxcor=nmxcor+1
              q1=10.d0**x(nmxcor)
            end if
            q2=10.d0**x(nmxcor-1)
            frmxc=(q2-qmxcor)/(q2-q1)
            rmxcor= (frmxc*10.d0**y(in1,nmxcor)+
     *        (1.d0-frmxc)*10.d0**y(in1,nmxcor-1))/10.d0**y(in1,1)
            frmxit(itcmp)=frmxc
            nmxcit(itcmp)=nmxcor
            qmxcit(itcmp)=qmxcor
            rmxcit(itcmp)=rmxcor
            frmrit(itcmp)=frmxc
            nmxrit(itcmp)=nmxcor
c
c  set iosc_frz = 1 to freeze core completely
c
            iosc_frz = 1
c
            if(istdpr.gt.0) write(istdpr,'(/
     *        '' Freeze to average.''
     *        '' qmxcor ='',1pe13.5,'' nmxcor ='',i5,
     *        '' frmxc ='',0pf10.5)') 
     *         qmxcor, nmxcor, frmxc
	    if(nmxrit(itcmp).lt.nmxrit(itcmp-1)-1) then
	      write(istdou,'(/
     *          '' ***** Warning in cmpcvc.'',
     *          '' Freeze with nmxrit(itcmp), nmxrit(itcmp-1) ='',
     *          2i5)') nmxrit(itcmp), nmxrit(itcmp-1)
	      if(istdpr.gt.0.and.istdpr.ne.istdou)
     *          write(istdpr,'(/'' ***** Warning.'',
     *          '' Freeze with nmxrit(itcmp), nmxrit(itcmp-1) ='',
     *          2i5)') nmxrit(itcmp), nmxrit(itcmp-1)
            end if
c
            if(istdpr.gt.0) write(istdpr,
     *        '(/ '' it, nmxcor, qmxcor, frmxc, Xc, nmx_res,'',
     *        '' frmx_res, ddrmax:''
     *        /(2i5, 1pe13.5,0pf10.5,f12.8,i5,f10.5,1pe13.5))')
     *        (i, nmxcit(i), qmxcit(i), frmxit(i), xmxcit(i), 
     *        nmxrit(i), frmrit(i), drmxit(i),i=1,itcmp)
c
	  end if
	end if
c
c  After jump fix mixed core to size at previous time step for the
c  first few iterations for stability(?)
c
	if(istdpr.gt.0) write(istdpr,*) '#D# itcmp, eamcon, epscon =',
     *    itcmp, eamcon, epscon 
	if(icjump.ge.1.and.imixc5.ge.1.and.itcmp.le.5.and.
     *    eamcon.ge.1.5*epscon) then
	  nmxcor=nmxcp
	  qmxcor=qmxcp
	  frmxc=1.d0
	  it_force=1
	  if(istdpr.gt.0) write(istdpr,
     *      '(/'' Warning. Fix nmxcor, '',
     *         ''qmxcor to value at previous time step''/
     *         '' nmxcor, qmxcor ='',i5,1pe13.5/)')
     *      nmxcor, qmxcor
        end if
c
c  test for oscillating mixed core
c
	if(itcmp.eq.1) nxcrit=0
	if(imixc5.ge.1.and.itcmp.ge.10.and.icjump.gt.0
     *    .and.nxcrit.eq.0) then
	  if(nmxcit(itcmp).eq.nmxcit(itcmp-1)-1.and.
     *      nmxcit(itcmp).eq.nmxcit(itcmp-2)) then
	    nxcrit=nmxcit(itcmp)
	    xcrit=y(in1+4,nxcrit)
	    if(istdpr.gt.0) write(istdpr,
     *        '(/'' Oscillating core, n_crit, xcrit  ='', i5,f10.6)') 
     *        nxcrit, xcrit
	    nmxcor=nxcrit
	    frmxc=1.d0
	    qmxcor=10.d0**(x(nxcrit)-1.d-7)
          end if
c..	  nmxcor=nmxcor_itp
c..	  qmxcor=qmxcor_itp
c..	  if(istdpr.gt.0) write(istdpr,
c..     *      '(/'' Warning. Fix nmxcor, qmxcor to'',i5,1pe13.5/)')
c..     *      nmxcor_itp, qmxcor_itp
        end if

c
c  test for presence of mixed core
c
        if(nmxcor.eq.0) then
          call zero(crxmn,icomp)
	  if(istdpr.gt.0) then
            write(istdpr,110)
            write(istdpr,'(/'' Exiting cmpcvc''/)')
	  end if
          return
        end if
c
c  when this is onset of convective core, store central values
c  in compc
c
        if(qmxcp.eq.0) then
          if(istdpr.gt.0) write(istdpr,112) qc
          do i=1,ifcomp
            compc(i)=y(in1+3+i,nn)
	  end do
          qcp=qmxcp
        end if
c
c  carry out secant iteration for actual size of convective core
c
        dqmxp = dqmx
        dqmx = qmxcor - qmxcrt
        if(istdpr.gt.0) write(istdpr,'(a,i3,1p5e13.5)') 
     *  ' nitcmp, qmxcrp, dqmxp, qmxcrt, qmxcor, dqmx',
     *    nitcmp, qmxcrp, dqmxp, qmxcrt, qmxcor, dqmx
c
        if(nitcmp.ge.3.and.imixc3.eq.2) then
          qmxnew = qmxcrt - dqmx*(qmxcrt - qmxcrp)/(dqmx - dqmxp)
c
	  if(qmxnew.gt.0.and.qmxnew.lt.0.9) then
c
c  reset nxmcor, frmxc
c
	    qmxcor=qmxnew
            qmxcrl=log10(qmxcor)
            do n=1,nn
              nmxcor=n
              if(x(n).le.qmxcrl) go to 39
            end do
c
   39       frmxc=(qmxcrl-x(nmxcor-1))/(x(nmxcor)-x(nmxcor-1))
            if(istdpr.gt.0) write(istdpr,*) 
     *        ' Reset qmxcor, nmxcor, frmxc to ',
     *        qmxcor, nmxcor, frmxc
          end if
        end if
c
c  prepare test for adding meshpoints if core has jumped a long
c  distance.
c  For now, only for imixc5 = 1, although it should probably
c  be included more generally.
c  
	if(istdpr.gt.0.and.ntime.eq.22) 
     *    write(istdpr,*) '#D# 2. y2, X =', y(in1+1,483), y(in1+4,483)
      if(nmxcp.gt.0.and.log10(qmxcor/qmxcor_prev).ge.3.d0*dqlmsc.and.
     *  imixc5.ge.1) then
	xrefine(1)=log10(qmxcor)
	ifixnn=0
        call mrefine(x,y,in,in1,iy,nvar,nn,xrefine,1,dxcrit,dqlmsc,
     *    ifixnn,irsmsh)
        qmxcor_prev=qmxcor
      end if
c
c  test for setting composition change from growing mixed region
c
c  To prepare for undercorrection, set value of Xc for previous iteration
c
	xcnuc_prev=xcnuc
c
c  store core abundances from previous timestep
c
	do i=1,ifcomp
	  compcp(i)=y(in+3+i,nn)
        end do
c
	if(imixc5.gt.0) then
c
c  reset hydrogen abundance at previous time step,
c  as average over past or present core 
c  (possible contribution from growing core already included above)
c
	  qmix=min(qmxcp, qmxcor)
          compcp(1)=comp_aver(qmix,y(in+4,1),iy,x,nn)
c
          if(istdpr.gt.0) then
	    write(istdpr,*) '#D# qp, xhp, qxint', qp, xhp, qxint
	    write(istdpr,*) '#D# X_c(prev), Xmean(prev) =',
     *        y(in+4,nn), compcp(1)
          end if
c
        end if
	i=-nn
	if(i.gt.0.and.istdpr.gt.0) write(istdpr,*) 
     *    '#D# 3. i =',i,'  y(in1+4,i), ymix(5,i) =', 
     *    y(in1+4,i), ymix(5,i)
c
        if(qmxcor.gt.qmxcp) then
c
          call cmpgmx(x,y(in,1),y(in1,1),iy,iy,nn,qmxcor,qmxcp,
     *      y(in+4,nn),dcmpgr,ifcomp)
c
c  Reset composition. If qcp = 0, only set average returned in dcmpgr
c  Change to new formulation for imixc6 = 1, 20/2/08
c
          do i=1,ifcomp
	    if(imixc6.eq.1) then
              compc(i)=max(compcp(i)*(qmxcp/qmxcor)+dcmpgr(i),1.d-10)
            else
              compc(i)=max(compcp(i)+dcmpgr(i),1.d-10)
	    end if
            if(istdpr.gt.0) then
              write(istdpr,*) 'Set new compc from previous value'
              write(istdpr,'(''Element'',i2,1p3e12.4)') 
     *          i, compcp(i),dcmpgr(i),compc(i)
	    end if
          end do
c
c  otherwise initialize compc and set dcmpgr to zero, for output 
c  purposes
c
	else
c
c  set value from previous value at innermost mesh point
c
	  do i=1,ifcomp
            compc(i)=max(compcp(i),1.d-10)
            dcmpgr(i)=0.d0
          end do
        end if
c
c  the following statement arises because of confusion between 
c  different measures of core edge.
c
        if(qmxcp.gt.0.and.qcp.eq.0) qcp=qmxcp
c
c  set average rx for this timestep
c
        call rxhmxc(x,y(in1,1),iy,nn,compc,qc,rxmn,drxmn,icomp,alres,
     *    iprdcr,icryrx)
	if(icryrx.lt.0) then
	  icrycm=-1
	  go to 90
        end if
c
c  reset luminosity in ymix from alres
c
c  SWITCH OFF FOR TESTING PURPOSES; NEEDS THOUGHT!!
c
c..        do 42 n=1,nn
c..   42   ymix(4,n)=alres(n)
c
c  store values at current time step in common
c
        cqc=qc
        call store(rxmn,crxmn,icomp)
c
c  with growing convective core, include contribution from enlarged
c  region
c  Note: correcting for centralization now implemented in s/r rhs
c  (from 4/2/08)
c
	if(qmxcor.gt.qmxcp.and.imixc6.eq.1) then
	  do i=1,ifcomp
	    if(istdpr.gt.0) write(istdpr,*) 
     *        '#D# crxmn, dcmpgr term, sum:',
     *        crxmn(i),dcmpgr(i)/dt,
     *        crxmn(i)+dcmpgr(i)/dt
            crxmn(i)=crxmn(i)+dcmpgr(i)/dt
          end do
        end if
c
c  reset corrected new core abundances. Possible
c  contribution from growing core already included.
c  Changed to new formulation for growing core and imixc6 = 1, 20/2/08
c
	icrycm=0
c
        do 45 i=1,icomp
        if(i.ne.2.or.iche3.ne.1) then
          if(qmxcp.gt.0) then
            thetai=theta(4+i)
            phii=1-thetai
	    if(imixc6.ne.1.or.qmxcor.le.qmxcp) then
              compc(i)=max(compc(i)+
     *          dt*(phii*rxmnp(i)+thetai*rxmn(i)),1.d-10)
            else
              compc(i)=max(compc(i)+
     *          dt*(phii*rxmnp(i)*(qmxcp/qmxcor)+thetai*rxmn(i)),1.d-10)
	    end if
            if(istdpr.gt.0) then
              write(istdpr,*) 'Set new compc from previous value'
              write(istdpr,'(''Element'',i2,1p5e12.4)') 
     *          i, compcp(i),dcmpgr(i),
     *          dt*phii*rxmnp(i),dt*thetai*rxmn(i),compc(i)
            end if
          end if
c
c  test for unphysical value
c
	  if(compc(i).gt.1.d0) icrycm=-i
        end if
   45   continue
	if(istdpr.gt.0.and.ntime.eq.22) 
     *    write(istdpr,*) '#D# 3. y2, X =', y(in1+1,483), y(in1+4,483)
c
c  In case of diffusion, add contribution from diffusive flux.
c  So far just for hydrogen.
c  (Added 22/3/07)
c
        if(idiffus.ge.1) then
	  iydif=4+idcomp+iccomp
          thetax=theta(5)
          phix=1-thetai
	  compc_org=compc(1)
c
c  use current timestep
c
          compc(1)=max(compc(1)+dt*y(in1+iydif,nmxcor)/qmxcor,1.d-10)
        end if
c
	i=-nn
	if(i.gt.0.and.istdpr.gt.0) write(istdpr,*) 
     *    '#D# 4. i =',i,'  y(in1+4,i), ymix(5,i) =', 
     *    y(in1+4,i), ymix(5,i)
c
c  Store nuclear Xc (before semiconvective mixing)
c
	xcnuc=compc(1)
c
        if(iche3.eq.1) then
c
c  set average reaction rates for computing He3 abundance
c
          call almmxc(x,y(in1,1),iy,nn,compc,alrhmn,6,1)
c
c  set predicted new core abundance of He3
c
          xh=compc(1)
          yh=1-xh-zhc
          do 47 i=1,6
   47     alrhmc(i)=0.5*(alrhmp(i)+alrhmn(i))
          idghe3=1
          if(istdpr.gt.0) write(istdpr,*) xh, yh, zhc, dt, y(in+5,nn)
          call he3abc(xh,yh,zhc,dt,y(in+5,nn),alrhmc,x3new,x3eq,anu,
     *      1,nosd)
          idghe3=0
          compc(2)=x3new(1)
c
c  reset average rate of change in convective core
c  (this requires further thought; 21/1/96)
c
          rxmn(2)=(x3new(1)-y(in+5,nn))/dt
          crxmn(2)=rxmn(2)
c
        end if
	i=-nn
	if(i.gt.0.and.istdpr.gt.0) write(istdpr,*) 
     *    '#D# 5. i =',i,'  y(in1+4,i), ymix(5,i) =', 
     *    y(in1+4,i), ymix(5,i)
c
        if(frmxc.gt.1) then
          nc1=nmxcor+1
        else
          nc1=nmxcor
        end if
        if(istdpr.gt.0) write(istdpr,115) 
     *    nitcmp,(i, y(in1+3+i,nn),compc(i),i=1,ifcomp)
c
c  test for error in setting new abundance
c
	if(icrycm.lt.0) then
	  write(istdou,120)
	  if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,120)
          if(istdpr.gt.0) write(istdpr,'(/'' Exiting cmpcvc''/)')
	  return
        end if
c
        do 54 n=1,nn
        qx=10.d0**x(n)
c
c  changed from qx.gt.qmxcor 29/3/06
c
        if(qx.ge.qmxcor) then
          do 50 i=1,ifcomp
   50     ymix(4+i,n)=cmporg(i,n)
        else
          do 52 i=1,ifcomp
   52     ymix(4+i,n)=compc(i)
        end if
c
   54   continue
c
c  reset composition in region between old and new convective core
c
        if(qmxcor.ge.qmxcp) then
c
c  growing convective core
c
          if(istdpr.gt.0) write(istdpr,130) qmxcp, qmxcor
c
        else if(nmxcp.lt.nc1.and.imixc6.ne.1) then
c
c  shrinking convective core.
c
c  Note that for imixc6 this should be correctly treated already
c  in rhs and hence no resetting is required here. (The present
c  resetting is inconsistent with the setup in rhs, in any case.)
c
c  reset in intermediate region. Since logics of nmxcor and nc1
c  is a little unclear, make a careful check on whether we
c  are actually in the intermediate region
c
c  test for applying undercorrection to qmxcor
c
          if(itcmp.gt.10.and.imixc5.eq.1) then
            if(itcmp.le.20) then
              ucyxc=0.7d0
            else if(itcmp.le.30) then
              ucyxc=0.4d0
            else 
              ucyxc=0.d0
            end if
            qmxcor_tfr=ucyxc*qmxcor+(1.d0-ucyxc)*qmxcor_tfrp
            if(istdpr.gt.0) write(istdpr,*)
     *        '#D# cmpcvc: qmxcor, qmxcor_tfr', qmxcor, qmxcor_tfr
          else
	    ucyxc=1.d0
            qmxcor_tfr=qmxcor
          end if
c
          if(istdpr.gt.0) write(istdpr,
     *      '(/'' Setting composition in shrinking region'')')
          do 56 n=1,nn
          qx=10.d0**x(n)
	  if(qx.ge.qmxcor.and.qx.le.qmxcp) then
            tfrct=(qmxcp-qx)/(qmxcp-qmxcor_tfr)
	    if(tfrct.gt.1) then
	      if(istdpr.gt.0) write(istdpr,*) 
     *          '#D# tfrct reset from ', tfrct,' to 1'
	      tfrct=1.d0
            end if
            do 55 i=1,icomp
            icp=in +3+i
   55       ymix(4+i,n)=y(icp,n)+dt*(1-tfrct)*rxstr(i,n)+
     *               tfrct*(ymix(4+i,nn)-y(icp,nn))
	    if(istdpr.gt.0) write(istdpr,'(i4,1p6e11.3)') 
     *          n, qx, tfrct, y(in+4,n), dt*rxstr(1,n),
     *          ymix(5,nn)-y(in+4,nn), ymix(5,n)
          end if
   56     continue
        end if
c
c initialize semiconvection quantities to zero
c
        nmxscn=0
        qmxscn=0
	ierr=0
c
c  test for resetting composition in possible semiconvective region
c
	if(imixc5.gt.0.and.qmxcor.ge.1.d-5.and.(qmxcor.gt.qmxcp.or.
     *    (nmxscp.gt.0.and.qmxcor.ge.qmxscp)).and.
     *    no_resscn.eq.0) then
c
	  call resscn(x,ymix,iymix,nn,compc,icomp,itcmp,ierr)
c
c  test for error
c
	  if(ierr.eq.-1) then
	    no_resscn=1
	    if(istdpr.gt.0) write(istdpr,'(/
     *        '' Block further calls of resscn at this timestep''/)')
	  else if(ierr.eq.-2) then
	    icrycm=-2
            if(istdpr.gt.0) write(istdpr,'(/'' Exiting cmpcvc''/)')
            return
          end if
c
	end if
	i=-nn
	if(i.gt.0.and.istdpr.gt.0) write(istdpr,*) 
     *    '#D# 6. i =',i,'  y(in1+4,i), ymix(5,i) =', 
     *    y(in1+4,i), ymix(5,i)
c
c  test for setting zhc from mixed composition
c
	if(idiffus.ge.2) zhc=compc(ifcomp)
	if(istdpr.gt.0.and.ntime.eq.22) 
     *    write(istdpr,*) '#D# 4. y2, X =', y(in1+1,483), y(in1+4,483)
c
c  test reset composition
c
        call tstxin(x,ymix,nn,iymix,1,xhint,xhintp,xtime,xtime,icry,
     *    'after iteration')
c..        write(38,'(4i5,1p2e13.5)') 
c..     *    (ntime,iter,nitcmp,n,x(n),ymix(5,n),n=nmxcor-10,nmxcor+10)
c
c  test for convergence or too many iterations
c
        dxcit=compc(1)-xcitp
        xcitp=compc(1)
        if(nitcmp.lt.nitcmx.and.abs(dxcit).gt.epscmx) go to 35
c
c  end of iteration. Store final composition in y(in1 + .,.)
c
c  test for applying undercorrection to ymix
c
        if(itcmp.gt.10.and.imixc5.eq.1) then
          if(itcmp.le.20) then
            ucyxc=0.7d0
          else if(itcmp.le.30) then
            ucyxc=0.4d0
          else 
            ucyxc=0.d0
          end if
	  do n=1,nn
	    do i=1,4+icomp
	      ymix(i,n)=ucyxc*ymix(i,n)+(1.d0-ucyxc)*ymixp(i,n)
	    end do
	  end do
        else
	  ucyxc=1.d0
        end if
c
	if(istdpr.gt.0) write(istdpr,
     *    '(/'' Restore composition with ucyxc ='',f10.7/)') ucyxc
c
c  test for resetting X at n_crit to force neutral stability
c
	if(nxcrit.gt.0) then
	  do i=1,nvar
	    ycrit(i)=ymix(i,nxcrit)
	  end do
          pl=addvar(1,nxcrit)
          call set_neutralx(x(nxcrit),ycrit,pl,zh(nxcrit),fl)
	  if(istdpr.gt.0) write(istdpr,*) 
     *      '#D# reset at nxcrit. xcrit, ymix(5,nxcrit), reset X =',
     *      xcrit, ymix(5,nxcrit), ycrit(5)
	  ymix(2,nxcrit)=ycrit(2)
	  ymix(5,nxcrit)=ycrit(5)
        end if
	if(istdpr.gt.0.and.ntime.eq.22) 
     *    write(istdpr,*) '#D# 5. y2, X =', y(in1+1,483), y(in1+4,483)
c
	qmxres=max(qmxcor,qmxcp)
	i=-nn
	if(i.gt.0.and.istdpr.gt.0) write(istdpr,*) 
     *    '#D# 7. i =',i,'  y(in1+4,i), ymix(5,i) =', 
     *    y(in1+4,i), ymix(5,i)
	if(istdpr.gt.0.and.ntime.eq.22) 
     *    write(istdpr,*) '#D# 6. y2, X =', y(in1+1,483), y(in1+4,483)
c
	do 58 n=1,nn
	qx=10.d0**x(n)
	if(qx.le.qmxres) then
c
c  for imixc5 = 1 for now (but possibly later extended) 
c  reset log f (log rho) to keep pressure fixed
c  It is assumed that log f has been properly set in ymix from
c  call of resscn
c
	  if(imixc5.eq.1.and.ierr.eq.0) then
            y(in1+1,n)=ymix(2,n)
	  end if
c
c  for imixc6 = 1 set correction to log rho (log f) to keep
c  pressure fixed, in a linearized approximation
c
c  Switch off 25/2/08, 14.00
c
	  if(imixc6.eq.-1.and.ierr.eq.0.and.itcmp.ge.2) then
	    dxh=ymix(5,n)-y(in1+4,n)
	    fl=y(in1+1,n)
	    if(abs(dxh).ge.1.d-10) then
	      tl=ymix(3,n)
	      xh=ymix(5,n)
	      yh=1-xh-zh(n)
	      pl=addvar(1,n)
	      fl2=fl
	      call eqstp(pl,tl,xh,yh,zh(n),nosd,nosd,fl2,nit)
              flnew=fl2
            else
	      flnew=fl
            end if
c..	    if(istdpr.gt.0) 
c..     *        write(istdpr,*) '#D# n, qx, y(in1+1,n),addvar(1,n) etc',
c..     *        n, qx, y(in1+1,n),addvar(1,n),pl,ymix(5,n),y(in1+4,n),dxh,
c..     *        flnew
	    y(in1+1,n)=flnew
          end if
c
c  for imixc5 eq 1 only reset from ymix if there was no error
c  otherwise go back to compc
c
	  if(imixc5.ne.1.or.ierr.eq.0) then
	    do i=1,ifcomp
              y(in1+3+i,n)=ymix(4+i,n)
            end do
          else
	    do i=1,ifcomp
              y(in1+3+i,n)=compc(i)
            end do
          end if
c
c  otherwise reset ymix to be sure to have appropriate value for next
c  iteration
c
	else
	  do i=1,ifcomp
	    ymix(4+i,n)=y(in1+3+i,n)
          end do
	end if
   58   continue

c
c#ai#  As a temporary fix, store mixed compositions in cvr
c#ai#  This has to be cleared up with care.
c
        do 60 n=1,nn
	qx=10.d0**x(n)
        cvr(1,n)=ymix(5,n)
	if(qx.gt.qmxcor) then
          cvr(2,n)=1-ymix(5,n)-zh(n)
        else
          cvr(2,n)=1-ymix(5,n)-zhc
	end if
        if(iche3.eq.1) cvr(3,n)=ymix(6,n)
c
c  store CNO abundances, preferentially from values set
c  in compc. 
c
c  Note: replaced 15/8/02 by call of setcno, based on values in y.
c
	call setcno(y(in1,n),n)
c
c  test for storing 4He burning abundances
c
	if(iheccs.ne.0) then
	  cvr(icvhe4,n)=y(in1+iyche4-1,n)
	  cvr(icvc12,n)=y(in1+iycc12-1,n)
        end if
   60   continue
	do 65 i=1,icvrmx
   65   cvr(i,nn+1)=cvr(i,nn)
	i=-nn
	if(i.gt.0.and.istdpr.gt.0) write(istdpr,*) 
     *    '#D# 8. i =',i,'  y(in1+4,i), ymix(5,i) =', 
     *    y(in1+4,i), ymix(5,i)
c
c  store rate of change of X in mixed region (for use in s/r rhs)
c
	if(imixc5.eq.1) then
	  do n=nmxcor,nn
	    crxxh_sc(n)=(y(in1+4,n)-y(in+4,n))/dt
          end do
        else
	  do n=nmxcor,nn
	    crxxh_sc(n)=crxmnp(1)
          end do
        end if
c
c  test for resetting log f (log rho) in region of increasing core, to
c  ensure smooth pressure.
c
        if(imixc5.eq.1) then
c
c  always reset for imixc5 = 1
c
          call crslgf(x,y(in1,1),iy,nn,qmxcor,qmxcp,itcmp,
     *      ncp1,nn,plgp1,plgp2,icrycm)
        else if(qmxcor.gt.qmxcp.and.iter.le.itfrz1+1) then
c
c  otherwise only reset region of growing convective core
c
          call crslgf(x,y(in1,1),iy,nn,qmxcor,qmxcp,itcmp,
     *      ncp1,ncp2,plgp1,plgp2,icrycm)
        end if
	if(icrycm.lt.0) go to 90
	i=-nn
	if(i.gt.0.and.istdpr.gt.0) write(istdpr,*) 
     *    '#D# 9. i =',i,'  y(in1+4,i), ymix(5,i) =', 
     *    y(in1+4,i), ymix(5,i)
c
c  Locate critical points and possibly add meshpoints
c
	if(irevc6.eq.1.and.itcmp.eq.4) then
c
	  nn_old=nn
          call ddr_ref(x,y,in,in1,iy,nvar,nn,dxcrit,irsmsh)
	  if(irsmsh.eq.1) then
	    if(istdpr.gt.0) write(istdpr,'(/
     *        '' Points added in ddr_ref. New nn ='',i5)') nn
	    frmxit(itcmp)=frmxc
	    nmxcit(itcmp)=nmxcor
	    frmrit(itcmp)=frmxc
	    nmxrit(itcmp)=nmxcor
	    itoscl=itcmp+2
	    itoscs=itcmp+4
            if(istdpr.gt.0) write(istdpr,
     *        '(/ '' it, nmxcor, qmxcor, frmxc, Xc, nmx_res,'',
     *        '' frmx_res, ddrmax:''
     *        /(2i5, 1pe13.5,0pf10.5,f12.8,i5,f10.5,1pe13.5))')
     *        (i, nmxcit(i), qmxcit(i), frmxit(i), xmxcit(i), 
     *        nmxrit(i), frmrit(i), drmxit(i),i=1,itcmp)
c
          end if
        end if
c
c  end corrector step
c
      end if
c
   90 continue
c
      if(istdpr.gt.0) write(istdpr,'(/'' Exiting cmpcvc''/)')
      return
  101 format(/' In s/r cmpcvc, extent of mixed region has been frozen')
  102 format(//' *** Warning. ',
     *         'In s/r cmpcvc, no convective-core boundary ',
     *         ' has been frozen'/)
  103 format(//' *** Warning. In s/r cmpcvc, convective-core boundary',
     *         ' has been frozen at'/
     *         ' nmxcor =',i5,'  qmxcor =',1pe13.5,'  rmxcor =',e13.5,
     *         ' frmxc =',0pf10.4)
  105 format(//' In cmpcvc hydrogen found near exhaustion.',
     *  ' qc =',1pe13.5, ' Xc =',1pe13.5/
     *  ' Mixed core disabled')
  106 format(//' In cmpcvc helium found near exhaustion.',
     *  ' qc =',1pe13.5, ' Yc =',1pe13.5/
     *  ' Mixed core disabled')
  110 format(/' ***** Warning. No mixed core found.',
     *        ' Return from s/r cmpcvc')
  112 format(//' ***** Note. Onset of convective core, qc =',1pe13.5)
  115 format(/' Composition iteration no.',i3//
     *       (' X(',i1,')  set in convective core.',
     *    ' Old value =',1pe13.5,' New value =',e13.5))
  120 format(/' ***** Error in s/r cpmcvc; return')
  130 format(//' **** Warning. Growing convective core.'/
     *         '      Old mass =',1pe13.5,'  New mass =',e13.5)
      end
      subroutine cmpgmx(x,yp,y,iyp,iy,nn,qc,qcp,compc,
     *  dcmpgr,ifcomp)
c
c  sets change in composition resulting from growing mixed 
c  region.
c
      implicit double precision(a-h, o-z)
      include 'engenr.cz.d.incl'
      parameter (nnmaxc=2000, nspcm2=2*nspcmx)
      dimension x(*),yp(iyp,*),y(iy,*),compc(*),dcmpgr(*),
     *  cmp(nspcmx,nnmaxc), qi(nnmaxc), compc_av(nspcmx)
      common/rhcn/ a1,a2,a3,a4,zdummy,nvar,ifwrt,irhtst,inentr,
     *  ii1,ii2,ii3,icomp
      common/step/ dt
      common/cmxcnt/ ddrmix, inmixc, imixc0, imixc1, imixc2, imixc3,
     *  imixc4, imixc5, imixc6, iosc_frz, imix_repeat, icnv_mixcor
      common/crxstr/ rxstr(nspcm2,1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      if(istdpr.gt.0) write(istdpr,*) 'Enter cmpgmx'
c
c  find location of growing region, set composition sum
c
      idiag=1
      if(istdpr.le.0) idiag=0
      n1=0
      iend=0
      qx=1
      if(idiag.gt.0) 
     *  write(istdpr,*) '#D# setting X composition for mixing'
      do 30 n=1,nn
      qxp=qx
      qx=10.d0**x(n)
      nc=n
      if(qx.le.qc.and.n1.eq.0) then
c
c  first point. interpolate.
c  For imixc6 set composition from forward composition change
c
c..        if(istdpr.gt.0) write(istdpr,*) 
c..     *    'first. n, qc, qcp, qx, qxp', n, qc, qcp, qx, qxp
        n1=1
        frct=(qc-qxp)/(qx-qxp)
	if(n.gt.1) then
          do i=1,ifcomp
	    if(imixc6.ne.1) then
              cmp(i,n1)=0.5*(frct*(yp(4+i,n)+y(4+i,n))+
     *            (1-frct)*(yp(4+i,n-1)+y(4+i,n-1)))
	    else
	      i1=i+icomp
              cmp(i,n1)=frct*(yp(4+i,n)+dt*rxstr(i1,n))+
     *          (1-frct)*(yp(4+i,n-1)+dt*rxstr(i1,n-1))
	    end if
	  end do
        else
c
c  A rather desparate attempt at treating the case of a 
c  fully convective star. NEEDS a check.
c
          do i=1,ifcomp
            cmp(i,n1)=0.5*(yp(4+i,n)+y(4+i,n))
	  end do
	end if
        qi(n1)=qc
      end if
      if(qx.le.qcp) then
c
c  last point. interpolate, depending in imixc6.
c
        n1=n1+1
c..        if(istdpr.gt.0) write(istdpr,*) 
c..     *    'last. n, qc, qcp, qx, qxp', n, qc, qcp, qx, qxp
        frct=(qcp-qxp)/(qx-qxp)
        do i=1,ifcomp
	  if(imixc6.ne.1) then
            cmp(i,n1)=0.5*(frct*(yp(4+i,n)+y(4+i,n))+
     *         (1-frct)*(yp(4+i,n-1)+y(4+i,n-1)))
	  else
            cmp(i,n1)=frct*yp(4+i,n)+(1-frct)*yp(4+i,n-1)
	  end if
	end do
        qi(n1)=qcp
        go to 35
      else if(qx.lt.qc) then
c
c  set intermediate point. For imixc6 = 1 set mixed composition 
c  from forward composition evolution, assuming linear change in
c  the core boundary with time.
c
	n1=n1+1
        if(imixc6.ne.1) then
          do i=1,ifcomp
            cmp(i,n1)=0.5*(yp(4+i,n)+y(4+i,n))
	  end do
        else
	  tfrct=dt*(qx-qcp)/(qc-qcp)
          do i=1,ifcomp
            cmp(i,n1)=yp(4+i,n)+tfrct*rxstr(i+icomp,n)
	  end do
	end if
	if(idiag.gt.0) write(istdpr,'(i5,1p4e13.5)') n1, qx,
     *    yp(5,n), y(5,n), cmp(1,n1)
        qi(n1)=qx
      end if
   30 continue
c
   35 continue
c
      if(istdpr.gt.0) then 
        write(istdpr,110) qcp, qc, (compc(i),i=1,ifcomp)
        do 37 n=1,n1
   37   write(istdpr,115) qi(n),(cmp(i,n),i=1,ifcomp)
      end if
c
c  set integral
c  For now, as reference use average between previous and present
c  central abundance. This needs further development.
c
      do 40 i=1,ifcomp
      compc_av(i)=0.5d0*(y(4+i,nn)+yp(4+i,nn))
   40 dcmpgr(i)=0
c
      if(n1.eq.1) then
        if(istdpr.gt.0) write(istdpr,130) qcp, qc
      else
        do n=2,n1
          ww=0.5*(qi(n-1)-qi(n))
	  if(imixc6.ne.1) then
            do i=1,ifcomp
c
c  Replace compc by compc_av 21/1/08
c
              dcmpgr(i)=dcmpgr(i)+(qi(n-1)-qi(n))*
     *          (0.5*(cmp(i,n-1)+cmp(i,n))-compc_av(i))
            end do
	  else
c
c  Change to new formulation, without average subtraction, 20/2/08
c
            do i=1,ifcomp
              dcmpgr(i)=dcmpgr(i)+(qi(n-1)-qi(n))*
     *          0.5*(cmp(i,n-1)+cmp(i,n))
	    end do
          end if
	end do
      end if
c
      do 60 i=1,ifcomp
   60 dcmpgr(i)=dcmpgr(i)/qc
c
      if(istdpr.gt.0) write(istdpr,140) (i,dcmpgr(i),i=1,ifcomp)
      return
  110 format(/' In cmpgmx qcp, qc =',1p2e13.5/
     *  '       compc:',10e13.5)
  112 format(/' q, composition:')
  115 format(1p10e13.5)
  130 format(//' **** Error in cmpgmx. Only one point'/
     *         '      qcp, qc =',1p3e13.5)
  140 format(//' Change in composition from growing core:'/
     *  ' i, delta X(i):'/(i5,1pe13.5))
      end
      subroutine tstmxc(x,y,iy,nn,icomp,icry)
c
c  scan through size of mixed core (or perhaps iterate) 
c  to test for extent of mixing
c
c  Original version: 18/8/03
c
      implicit double precision(a-h,o-z)
      include 'engenr.cz.d.incl'
      logical norct, nosd
c
      parameter(naztmx = nspcmx+3, 
     *  kaxst = 3+naztmx, nbmax = nspcmx+2)
      parameter(nspcm2=2*nspcmx,iymix=ivarmx)
c
      dimension x(*), y(iy,*)
      dimension qmxst(8),qmxstp(8), compc(10)
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor,nmxcp,qmxcp,rmxcp,
     *  qmxmin,qmxmax,qmxscn,qmxscp,nmxscn,nmxscp,icjump,iqcfrz
      common/noiter/ iter, ntime, epscon, eamcon, it_force
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      equivalence(qmxst(1),qmxcor)
c
      icry=0
c
      if(istdpr.gt.0) write(istdpr,
     *  '(// '' Test mixing and instability'')')
c
      call store(qmxst,qmxstp,8)
c
      nqmx=40
      dqmn=1.5d-2/nqmx
      if(qmxcor.le.0) then
        dqmx=dqmn
      else
        dqmx=4.0*qmxcor/nqmx
      end if
      dqmx=max(dqmx,10.d0**x(nn-3))
      if(dqmx.lt.dqmn) then
	i0=3
	dqmx=(1.5d-2-1.5*qmxcor)/(nqmx-3)
	dqmx0=0.5*qmxcor
	qmx0=1.5*qmxcor
      else
	i0=0
	qmx0=dqmx
      end if
      if(mod(ntime,5).eq.0) then
	ifile=1
      else
	ifile=0
      end if
c
c  reset for test during iteration
c
      i0=0
      dqmx=0.01*dqmn
      qmx0=10.**x(nn-2)
      ifile=1
c
      init=1
      do i=1,nqmx
        if(i0.gt.0.and.i.le.3) then
	  qmx=dqmx0*i
        else
	  qmx=qmx0+dqmx*(i-i0)
	end if
        call resmxc(x,y,iy,nn,compc,icomp,qmx,ddradc,ifile,init,icry)
	if(icry.lt.0) return
	init=0
  	if(istdpr.gt.0) write(istdpr,'(/a,1p2e13.5)') 
     *    'In tstmxc, qmx, ddradc =', qmx,ddradc
      end do
c
c  restore cmxcor
c
      call store(qmxstp,qmxst,8)
      return
      end
      subroutine resmxc(x,y,iy,nn,compc,icomp,qmx,ddradc,ifile,init,
     *  icry)
c
c  evaluates properties of core mixed out to q = qmx, resetting
c  also luminosity but assuming that temperature and log f structure
c  are unchanged.
c
c  Original version: 18/8/03
c
      implicit double precision(a-h,o-z)
      include 'engenr.cz.d.incl'
      logical norct, nosd, time0, lastmd
c
      parameter(naztmx = nspcmx+3, 
     *  kaxst = 3+naztmx, nbmax = nspcmx+2)
      parameter(nspcm2=2*nspcmx,iymix=ivarmx)
c
      dimension x(*), y(iy,*), compc(*)
      dimension rxmn(nspcmx), drxmn(nspcmx), rxmnp(nspcmx),yn(10),
     .  drxmnp(nspcmx),
     .  alrhmn(krnrmx), alrhmp(krnrmx), alrhmc(krnrmx), 
     .  dcmpgr(nspcmx),x3new(4),alres(nnmax)
      common/rhcn/ a1,a2,a3,a4,zdummy,nvar,ifwrt,irhtst,inentc,
     *  ii1,ii2,ii3,icomp1
      common/cmtime/ age, time0, lastmd
      common/compvr/ cvr(icvrmx,1)
      common/cstcvr/ icvh1, icvhe4, icvhe3, icvc12, icvc13, icvn14, 
     *  icvo16
      common/anwvar/ data(8), yi(istrmx,1)
      common/xnwvar/ xi(1)
      common/caddvr/ addvar(5,ngmax)
      common/eqstd/ xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp,gmm1(4)
      common/cstyst/ iyh1, iyhe3, iyc12, iyc13, iyn14, iyo16, 
     *  iyche4, iycc12
      common/cstazt/ iah1, iahe3, iac12, iac13, ian14, iao16, 
     *  iache4, iacc12
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,jnc,frcf(6),frcl(6)
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor,nmxcp,qmxcp,rmxcp,
     *  qmxmin,qmxmax,qmxscn,qmxscp,nmxscn,nmxscp,icjump,iqcfrz
      common/cmxcit/ qmxcit(nitmax),rmxcit(nitmax),frmxit(nitmax),
     .  nmxcit(nitmax),xmxcit(nitmax),frmrit(nitmax),nmxrit(nitmax),
     .  drmxit(nitmax),qmxnit(nitmax),itcmpc
      common/crxcor/ cqc, cqcp, crxmn(nspcmx), crxmnp(nspcmx),
     *  crxxh_sc(nnmax)
      common/cymixc/ ymix(iymix,nnmax)
      common/cmxcnt/ ddrmix, inmixc, imixc0, imixc1, imixc2, imixc3,
     *  imixc4, imixc5, imixc6, iosc_frz, imix_repeat, icnv_mixcor
      common/crxstr/ rxstr(nspcm2,1)
      common/thetac/ theta(ivarmx)
      common/rnratd/ al(10,krnrmx),norct
      common/heavy/ zatmos, zhc, zh(1)
      common/compos/ xseth, yset, xset3, xset12, xset13, xset14, xset16
      common/ksider/ al01,al21,aztst(naztmx),axst(naztmx)
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec, idcomp, iccomp
      common/cdiffu/ idiffus
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3
      common/noiter/ iter, ntime, epscon, eamcon, it_force
      common/work/ ydrad(20,nnmax)
c
c  dataset designators
c
      common/cdsdsg/ idstrl,idshvz,idszab,idsevl,idssum,idscen,idsnov,
     *  idsoov,idsgng,idsgsm,idstm1,idstm2,iducen,iddgm1
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      icry=0
c
c  test for output to file of drad quantities for original model
c
      if(ifile.gt.0.and.init.eq.1) then
c
        do n=1,nn
          pl=addvar(1,n)
          drad=fdradp(x(n),y(1,n),pl,zh(n),fl,ak,akr,akt,akx)
	  call store(y(1,n),ydrad(1,n),nvar)
	  ydrad(nvar+1,n)=pl
	  ydrad(nvar+2,n)=zh(n)
	  ydrad(nvar+3,n)=ak
	  ydrad(nvar+4,n)=drad-dad(1)
        end do
	icount=0
      end if
c
c  reset core parameters to match qmx (note that previous values
c  should have been stored for safekeeping before call of
c  this routine)
c
      if(istdpr.gt.0) write(istdpr,*) 
     *  'Entering resmxc with qmx, icomp =',qmx, icomp
      qlmx=log10(qmx)
      nmxcor=-1
      do n=1,nn
        if(x(n).le.qlmx.and.nmxcor.eq.-1) nmxcor=n
      end do
      frmxc=(x(nmxcor-1)-qlmx)/(x(nmxcor-1)-x(nmxcor))
      if(istdpr.gt.0) write(istdpr,*) 
     *  'In resmxc, nmxcor, frmxc =',nmxcor,frmxc
c
c  set mixed composition (assume original composition is the
c  same at centre and innermost meshpoint)
c
      xi(1)=0.d0
      do i=1,icomp
        yi(i,1)=y(4+i,nn)
      end do
      ni=1
      do n=nn,nmxcor,-1
        ni=ni+1
        xi(ni)=10.d0**x(n)
        do i=1,icomp
          yi(i,ni)=y(4+i,n)
        end do
      end do
      if(frmxc.lt.1) then
        ni=ni+1
	frmxc1=1.d0-frmxc
	xi(ni)=qmx
        do i=1,icomp
          yi(i,ni)=frmxc*y(4+i,nmxcor)+frmxc1*y(4+i,nmxcor-1)
        end do
      end if
      do i=1,icomp
        compc(i)=0
	do n=2,ni
          compc(i)=compc(i)+0.5d0*(yi(i,n)+yi(i,n-1))*(xi(n)-xi(n-1))
        end do
	compc(i)=compc(i)/qmx
      end do
      if(istdpr.gt.0) write(istdpr,*) 
     *  'qmx, ni, compc(1)',qmx, ni, compc(1)

c
c  set average rx and reset luminosity for this mixed region
c
      call rxhmxc(x,y,iy,nn,compc,qc,rxmn,drxmn,icomp,alres,iprdcr,
     *  icryrx)
      if(icryrx.lt.0) then
	icry=-1
	return
      end if
c
c  set radiative gradient at border of mixed region, using
c  extrapolated composition from radiative side
c
      frmxc1=1.d0-frmxc
      do i=1,3
        yn(i)=frmxc*y(i,nmxcor)+frmxc1*y(i,nmxcor-1)
      end do
      yn(4)=frmxc*alres(nmxcor)+frmxc1*alres(nmxcor-1)
      frext=(qlmx-x(nmxcor-2))/(x(nmxcor-1)-x(nmxcor-2))
      frext1=1.d0-frext
      do i=1,icomp
        yn(4+i)=frext*y(4+i,nmxcor-1)+frext1*y(4+i,nmxcor-2)
      end do
      pl=frmxc*addvar(1,nmxcor)+frmxc1*addvar(1,nmxcor-1)
      zhh=frext*zh(nmxcor-1)+frext1*zh(nmxcor-2)
c
      drad=fdradp(qlmx,yn,pl,zhh,ak,fl,akr,akt,akx)
      ddradc=drad-dad(1)
c
c  test for output to file of drad quantities for remixed model
c
      if(ifile.gt.0) then
c
	icount=icount+1
        do n=1,nn
          pl=addvar(1,n)
	  call store(y(1,n),yn,3)
	  yn(4)=alres(n)
	  if(n.ge.nmxcor) then
	    call store(compc,yn(5),icomp)
          else
	    call store(y(5,n),yn(5),icomp)
	  end if
          drad=fdradp(x(n),yn,pl,zh(n),fl,ak,akr,akt,akx)
	  call store(yn,ydrad(1,n),nvar)
	  ydrad(nvar+1,n)=pl
	  ydrad(nvar+2,n)=zh(n)
	  ydrad(nvar+3,n)=ak
	  ydrad(nvar+4,n)=drad-dad(1)
        end do
      end if
      return
      end
      subroutine crslgf(x,y,iy,nn,qc,qcp,itcmp,ncp1,ncp2,plgp1,plgp2,
     *  icrycm)
c
c  Reset log f (log rho) near region of growing convective core
c  to ensure essentially smooth pressure
c
c  If itcmp.ge.0 base region of resetting on previous and current size
c  of convective core.
c
c  If itcmp = -1, reset over interval ncp1, ncp2 (for use in predictor
c  step)
c
c  If itcmp = -2, use same interval, but reset only log p; store reset
c  log p in addvar(1,.); here plgp1 and plgp2 must be log(p) at ends
c  of interval.
c
c  If itcmp = -3, use same interval, with plgp1 and plgp2 being 
c  log(p) at ends c  of interval. Do full reset of log f
c
c  In any case, ncp1, ncp2, plgp1, plgp2 return the values at the ends
c  of interval.
c
c  Original version: 27/5/05
c
c  Modified 14/11/05 to allow resetting of pressure by integration
c  of hydrostatic equilibrium
c
c  Modified 16/11/05, to allow for flexible use.
c
      implicit double precision(a-h, o-z)
      logical nosd, notd, time0, lastmd
      include 'engenr.cz.d.incl'
      dimension x(*),y(iy,*)
      dimension plg(nnmax), plgr(nnmax)
      common/ln10/ amm
      common/rhcn/ a1,a2,a3,a4,zdummy,nvar,ifwrt,irhtst,inentc,
     *  ii1,ii2,ii3,icomp1
      common/cmtime/ age, time0, lastmd
      common/eqstd/ xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp,gmm1(4)
      common/heavy/ zatmos, zhc, zh(1)
      common/caddvr/ addvar(5,ngmax)
      common/cdgphs/ kdgeos, kdgopc, kdgeng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      kdgeos=1
c
      idiag=0
      if(istdpr.le.0) idiag=0
c
      if(idiag.gt.0) write(istdpr,'(/'' Entering s/r crslgf''/)')
c
c  set irint = 1 for resetting through integration of hydrostatic
c  equilibrium and using maximal region (added 14/11/05)
c
      irint=1
c
      if(idiag.gt.0) 
     *  write(istdpr,*) 'Enter crslgf, itcmp, nc11, nc22  =',
     *  itcmp, nc11, nc22
c
c  limiting departure
c
      epslim=1.d-9
c
      if(itcmp.ge.0) then
c
c  find location of growing region
c
        qx=1
        nc1=0
        nc2=0
        do n=1,nn
          qxp=qx
          qx=10.d0**x(n)
          if(qx.le.qc.and.nc1.eq.0) nc1=n-1
c
          if(qx.le.qcp.and.nc2.eq.0) nc2=n
        end do
c
        if(itcmp.le.1.or.irint.ne.1) then
          nc11=nc1-3
          nc22=nc2+3
        else
          nc11=min(nc11,nc1-3)
          nc22=max(nc22,nc2+3)
        end if
	if(ncp2.eq.nn) nc22=nn
      else
	nc11=ncp1
	nc22=ncp2
      end if
c
      if(itcmp.ge.-1) then
c
c  step through critical region setting log(pressure) with current log f
c
        nosd=.true.
        notd=.true.
        do n=nc11,nc22
          fl=y(2,n)
          tl=y(3,n)
          xh=y(5,n)
          yh=1-xh-zh(n)
          call eqstf(fl,tl,xh,yh,zh(n),nosd,notd)
          if(kdgeos.lt.0) then
	    icrycm=-1
	    return
          end if
          plg(n)=log10(pt(1))
        end do
      else
	plg(nc11)=plgp1
	plg(nc22)=plgp2
      end if
c
c  reset pressure in critical region
c
      if(irint.ne.1) then
c
c  compare with linearly interpolated pressure 
c
        dplg=(plg(nc22)-plg(nc11))/(x(nc22)-x(nc11))
	do n=nc11,nc22
	  plgr(n)=plg(nc11)+(x(n)-x(nc11))*dplg
        end do
c
      else
c
c  integrate hydrostatic equilibrium at given log r
c
	plgrr=plg(nc11)
	do n=nc11,nc22
	  ftr=-a2*10.d0**(2.d0*x(n)-4.d0*y(1,n))
	  fztr=ftr*10.d0**(-plgrr)
	  if(n.gt.nc11) then
c
c  iteration for new log p
c
	    nit=0
c
   10       nit=nit+1
	    dx=x(n)-x(n-1)
	    phi=plgrr-plgrp-0.5d0*dx*(fztrp+fztr)
	    dphi=1.d0+0.5d0*dx*amm*fztr
	    dplgrr=-phi/dphi
	    plgrr=plgrr+dplgrr
	    fztr=ftr*10.d0**(-plgrr)
	    if(nit.ge.10) then
	      write(istdou,*) 
     *          'Iteration for log(p) in crslgf unconverged. dplrr =',
     *          dplrr
	      if(idiag.gt.0.and.istdpr.ne.istdou) write(istdpr,*) 
     *          'Iteration for log(p) in crslgf unconverged. dplrr =',
     *          dplrr
	      stop 'In s/r crslgf'
	    else if(abs(dplgrr).ge.1.d-8) then
	      go to 10
            end if
          end if
	  plgr(n)=plgrr
	  fztrp=fztr
	  plgrp=plgrr
        end do
c
c  rescale for continuity at both ends (except if last point is at centre)
c
        if(nc22.lt.nn) then
  	  plgfct=(plg(nc22)-plg(nc11))/(plgr(nc22)-plgr(nc11))
        else
  	  plgfct=1.d0
	end if
	do n=nc11,nc22
	  plgrr=plgr(n)
	  plgr(n)=plgr(nc11)+plgfct*(plgr(n)-plgr(nc11))
c..     *      n, x(n), plgrr, plgr(n)
	end do
      end if
c
c  test for just resetting log p in addvar(1.,)
c
      if(itcmp.eq.-2) then
	do n=nc11,nc22
	  if(istdpr.gt.0) write(istdpr,*) '#D# n, x, old, new log p',
     *      n, x(n), addvar(1,n), plgr(n)
	  addvar(1,n)=plgr(n)
        end do
        if(idiag.gt.0) write(istdpr,'(/'' Exiting s/r crslgf''/)')
	return
c
      else if(itcmp.eq.-3) then
c
c  reset log f (log rho) to match reset pressure
c
	do n=nc11, nc22
	  tl=y(3,n)
	  xh=y(5,n)
	  yh=1-xh-zh(n)
	  pl=plgr(n)
	  call eqstp(pl,tl,xh,yh,zh(n),nosd,nosd,fl,nit)
          y(2,n)=fl
	end do
        if(idiag.gt.0) write(istdpr,'(/'' Exiting s/r crslgf''/)')
	return
      else
	ncp1=nc11
	if(ncp2.ne.nn) ncp2=nc22
	plgp1=plgr(nc11)
	plgp2=plgr(nc22)
      end if
c
c  compare with reset pressure and possibly reset log(f)
c
      if(idiag.gt.0) 
     *  write(istdpr,'(//'' n, x, log f, log p, (log p)_i in crslgf'')')
      do n=nc11, nc22
	if(idiag.gt.0)
     *    write(istdpr,'(i5,4f12.7)') n, x(n), y(2,n),plg(n),plgr(n)
	if(abs(plg(n)-plgr(n)).gt.epslim) then
c
c  iterate to determine new log f
c
	  tl=y(3,n)
	  xh=y(5,n)
	  yh=1-xh-zh(n)
	  fl=y(2,n)
	  nit=0
c
   30     call eqstf(fl,tl,xh,yh,zh(n),nosd,notd)
          if(kdgeos.lt.0) then
	    icrycm=-1
	    return
          end if
	  dfl=(plgr(n)-log10(pt(1)))/pt(2)
	  fl=fl+dfl
	  nit=nit+1
	  if(nit.ge.10) then
	    write(istdou,*) 'Iteration in crslgf unconverged. dfl =',
     *        dfl
	    if(idiag.gt.0.and.istdpr.ne.istdou) 
     *        write(istdpr,*) 'Iteration in crslgf unconverged. dfl =',
     *          dfl
	    fl=y(2,n)
	  else if(abs(dfl).ge.1.d-8) then
	    go to 30
          end if
c
	  y(2,n)=fl
	  plg(n)=log10(pt(1))
	  if(idiag.gt.0) write(istdpr,'(i5,4f12.7,a)') 
     *      n, x(n), y(2,n),plg(n),plgr(n),' reset'
        end if
      end do
      if(idiag.gt.0) write(istdpr,'(/'' Exiting s/r crslgf''/)')
      return
      end
      subroutine resscn(x,y,iy,nn,compc,icomp,itcmp,ierr)
c
c  Resets conposition in growing convective core, including 
c  semiconvective region (with composition set such that 
c  nabla_rad = nabla_ad). Resets log(f) (or log rho) to keep
c  fixed pressure, as given in addvar(1,.)
c
c  Note that qmxcor is taken (for now) from common/cmxcor/
c
c  Original version: 8/11/05
c
      implicit double precision(a-h,o-z)
      include 'engenr.cz.d.incl'
      logical norct, nosd, time0, lastmd
c
      parameter(naztmx = nspcmx+3, 
     *  kaxst = 3+naztmx, nbmax = nspcmx+2)
      parameter(nspcm2=2*nspcmx,iymix=ivarmx)
c
      dimension x(*), y(iy,*), compc(*)
      dimension rxmn(nspcmx), drxmn(nspcmx), rxmnp(nspcmx),yn(10),
     .  drxmnp(nspcmx),
     .  alrhmn(krnrmx), alrhmp(krnrmx), alrhmc(krnrmx), 
     .  dcmpgr(nspcmx),x3new(4),alres(nnmax)
      common/ln10/ amm
      common/rhcn/ a1,a2,a3,a4,zdummy,nvar,ifwrt,irhtst,inentc,
     *  ii1,ii2,ii3,icomp1
      common/cmtime/ age, time0, lastmd
      common/compvr/ cvr(icvrmx,1)
      common/cstcvr/ icvh1, icvhe4, icvhe3, icvc12, icvc13, icvn14, 
     *  icvo16
      common/anwvar/ data(8), yi(istrmx,1)
      common/xnwvar/ xi(1)
      common/caddvr/ addvar(5,ngmax)
      common/eqstd/ xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp,gmm1(4)
      common/cstyst/ iyh1, iyhe3, iyc12, iyc13, iyn14, iyo16, 
     *  iyche4, iycc12
      common/cstazt/ iah1, iahe3, iac12, iac13, ian14, iao16, 
     *  iache4, iacc12
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,jnc,frcf(6),frcl(6)
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor,nmxcp,qmxcp,rmxcp,
     *  qmxmin,qmxmax,qmxscn,qmxscp,nmxscn,nmxscp,icjump,iqcfrz
      common/cmxcit/ qmxcit(nitmax),rmxcit(nitmax),frmxit(nitmax),
     .  nmxcit(nitmax),xmxcit(nitmax),frmrit(nitmax),nmxrit(nitmax),
     .  drmxit(nitmax),qmxnit(nitmax),itcmpc
      common/crxcor/ cqc, cqcp, crxmn(nspcmx), crxmnp(nspcmx),
     *  crxxh_sc(nnmax)
      common/cymixc/ ymix(iymix,nnmax)
      common/cmxcnt/ ddrmix, inmixc, imixc0, imixc1, imixc2, imixc3,
     *  imixc4, imixc5, imixc6, iosc_frz, imix_repeat, icnv_mixcor
      common/crxstr/ rxstr(nspcm2,1)
      common/thetac/ theta(ivarmx)
      common/rnratd/ al(10,krnrmx),norct
      common/heavy/ zatmos, zhc, zh(1)
      common/compos/ xseth, yset, xset3, xset12, xset13, xset14, xset16
      common/ksider/ al01,al21,aztst(naztmx),axst(naztmx)
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec, idcomp, iccomp
      common/cdiffu/ idiffus
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3
      common/noiter/ iter, ntime, epscon, eamcon, it_force
      common/work/ ydrad(20,nnmax)
      common/cdgphs/ kdgeos, kdgopc, kdgeng
c
c  dataset designators
c
      common/cdsdsg/ idstrl,idshvz,idszab,idsevl,idssum,idscen,idsnov,
     *  idsoov,idsgng,idsgsm,idstm1,idstm2,iducen,iddgm1
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      kdgeos=1
c
      ierr=0
c
      if(istdpr.gt.0) write(istdpr,'(/'' Entering s/r resscn''/)')
c
c  integrate for initial hydrogen mass and reset boundary of core
c
      qxhcor=0.d0
      qp=0
      xhp=y(5,nn)
      do n=nn,1,-1
	q=10.d0**x(n)
	if(q.gt.qmxcor) go to 15
	nmxcrl=n
	xh=y(5,n)
	qxhcor=qxhcor+0.5d0*(q-qp)*(xh+xhp)
	qp=q
	xhp=xh
      end do
c
   15 continue
c
c  add contribution from remaining bit, extrapolating from interior
c
      q1=10.d0**x(nmxcrl)
      q2=10.d0**x(nmxcrl+1)
      fctxtr=(qmxcor-q2)/(q1-q2)
      xhextr=fctxtr*y(5,nmxcrl)+(1.d0-fctxtr)*y(5,nmxcrl+1)
      qxhcor=qxhcor+0.5d0*(xhextr+y(5,nmxcrl))*(qmxcor-q1)
      if(istdpr.gt.0) write(istdpr,*) '#D# q1, q2, fctxtr, xhextr', 
     *  q1, q2, fctxtr, xhextr
c
        if(istdpr.gt.0) 
     *    write(istdpr,*) '#D# q, xh(nn), q*xh(nn), qxhcor',
     *    qp, y(5,nn), qp*y(5,nn), qxhcor
c
c  step through core from the edge, testing for resetting composition
c  extrapolate from radiative region to get composition at boundary
c
      irs_xcor=1
      istop_scn=0
      qxhscn=0.d0
      q1=10.d0**x(nmxcrl-1)
      q2=10.d0**x(nmxcrl-2)
      fctxtr=(qmxcor-q2)/(q1-q2)
      xhextr=fctxtr*y(5,nmxcrl-1)+(1.d0-fctxtr)*y(5,nmxcrl-2)
      xhp=xhextr
      qp=qmxcor
      xhcor=y(5,nn)
c
      if(itcmp.eq.1) xhcor0=xhcor
      do n=nmxcrl-3, nn
        q=10.d0**x(n)
	pl=addvar(1,n)
	if(n.lt.nmxcrl) then
	  xhh=y(5,n)
        else
	  xhh=xhcor
	end if
	do i=1,nvar
	  yn(i)=y(i,n)
        end do
	yn(5)=xhh
        drad=fdradp(x(n),yn,pl,zh(n),fl,ak,akr,akt,akx)
        ddrad=drad-dad(1)
	if(n.le.nmxcrl+20.and.istdpr.gt.0) write(istdpr,'(a,i5,6f11.7)') 
     *    '#D# n, x(n), fl, tl, xhh, yn(4), ddrad', 
     *    n, x(n), yn(2),yn(3), xhh, yn(4), ddrad
c
	if(ddrad.lt.0.and.n.ge.nmxcrl.and.q.lt.qmxcor) then
c
	  it=0
c
c  start iterating to set new X
c
   20     it = it + 1
	  rhopx = (rho(4) - pt(4)*rho(2)/pt(2))*amm*xhh
	  akpx = akx + akr*rhopx
	  if(istdpr.gt.0) write(istdpr,'(a,i3,5f11.7)') 
     *      '#D# it, xhh, kappa, akpx, yn(4), ddrad', 
     *      it, xhh, ak, akpx, yn(4), ddrad
	  dxhh=-ddrad*xhh/(drad*akpx)
	  if(abs(dxhh).lt.1.d-7) then
          else if(it.gt.10) then
	    write(istdou,'(/
     *        '' ***** Warning. Excessive number of iterations'',
     *        '' in s/r resscn''/
     *        ''       n, log(q), log(f), log(T), X ='',i5, 1pe13.5,
     *        0p3f10.5)') n, x(n), yn(2), yn(3), xhh
	    if(istdpr.gt.0.and.istdou.ne.istdpr) write(istdpr,'(/
     *        '' ***** Warning. Excessive number of iterations'',
     *        '' in s/r resscn''/
     *        ''       n, log(q), log(f), log(T), X ='',i5, 1pe13.5,
     *        0p3f10.5)') n, x(n), yn(2), yn(3), xhh
          else
	    xhh=xhh + dxhh
	    yn(5)=xhh
            drad=fdradp(x(n),yn,pl,zh(n),fl,ak,akr,akt,akx)
            ddrad=drad-dad(1)
	    go to 20
          end if
c
c  test for (probably unphysical) growing X
c
c  Temporarily switch off blocking of resetting of X_c
c  (doubtful effect).
c
	  if(xhh.ge.xhp.and.n.ge.nmxcor+2) then
	    if(istop_scn.eq.0.and.istdpr.gt.0) write(istdpr,
     *        '('' X increasing with depth''/)')
c..	    irs_xcor=0
            istop_scn=1
          end if
c
c  set corresponding core composition and test for stop of
c  resetting
c
	  if(istdpr.gt.0) write(istdpr,*) '#D# qp, q, xhh, xhp, qxhscn',
     *      qp, q, xhh, xhp, qxhscn
	  qxhscn=qxhscn+0.5d0*(qp-q)*(xhh+xhp)
	  nmxscn=n
	  qp=q
	  xhp=xhh
	  if(irs_xcor.eq.1) then
	    xhcor=(qxhcor-qxhscn)/q
	  else
	    if(istdpr.gt.0) write(istdpr,
     *        '(/'' Resetting of xhcor blocked''/)')
	  end if
c
c  test for unphysical core abundance
c
          if(xhcor.le.0.d0) then
	    ierr=-1
	    if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,
     *        '(/'' ***** Error in s/r resscn. X_c = '',1pe13.5,
     *        '' negative''/)') xhcor
	    write(istdou,'(/'' ***** Error in s/r resscn. X_c = '',
     *        1pe13.5,'' negative''/)') xhcor
	    return
          end if
	  if(istdpr.gt.0) write(istdpr,'(a,i5,1p6e14.6)') 
     *      '#D# n, x, q, qxhcor, qxhscn, X, X_c:',
     *      n, x(n), q, qxhcor, qxhscn, xhh, xhcor
c
c  test out stop on increasing X, for now (seems to give problems)
c
c..	  if(xhh.le.xhcor.or.istop_scn.eq.1) then
	  if(xhh.le.xhcor) then
	    nmxscn=n-1
	    if(istdpr.gt.0) write(istdpr,*) 
     *        '#D# ntime, eamcon', ntime, eamcon
	    if((ntime.eq.50.or.ntime.eq.268.or.ntime.eq.269).and.
     *        eamcon.le.1.d-7) then
	      if(istdpr.gt.0) write(istdpr,*) '#D# calling test_kappa'
	      call test_kappa(x(n),yn,pl,zh(n))
	    end if
	    go to 30
          end if
	  y(5,n)=xhh
	  y(2,n)=fl
        else if(n.eq.nmxscn+1) then
	  if(istdpr.gt.0) write(istdpr,*) 
     *      '#D# ntime, eamcon', ntime, eamcon
	  if((ntime.eq.50.or.ntime.eq.268.or.ntime.eq.269).and.
     *      eamcon.le.1.d-7) then
	    if(istdpr.gt.0) write(istdpr,*) '#D# calling test_kappa'
	    call test_kappa(x(n),yn,pl,zh(n))
c
c  add test to skip rest of core if unstable (added 16/4/06)
c
	  end if
          go to 30
        end if
      end do
c
   30 continue
c
c  Note: nmxscn is last point in semiconvective region
c
      if(istdpr.gt.0) write(istdpr,*)
     *  '#D# End of scan; nmxscn, qxhscn =', nmxscn, qxhscn
      if(nmxscn.gt.0) then
        qmxscn=10.d0**x(nmxscn)
      else
        qmxscn=0.d0
      end if
c
c  test for presence of semiconvective region
c
      if(qxhscn.eq.0) then
	if(istdpr.gt.0) write(istdpr,'(/
     *    '' No semiconvective region found in s/r resscn''/)')
	return
      end if
c
c  reset X and log(f) in the inner homogeneous part of the core
c  For late iterations, try undercorrection for xhcor
c
      if(istdpr.gt.0) write(istdpr,*) 
     *  '#D# Core resetting. itcmp, xhcor, xhcor0 =',
     *    itcmp, xhcor, xhcor0
      if(itcmp.ge.10) xhcor=xhcor0+0.5*(xhcor-xhcor0)
      xhcor0=xhcor
c
      nosd=.true.
      do n=nmxscn+1,nn
	xhh=y(5,n)
	yhh=1.d0-xhh-zh(n)
	fl=y(2,n)
	tl=y(3,n)
        call eqstf(fl,tl,xhh,yhh,zh(n),nosd,nosd)
        if(kdgeos.lt.0) then
	  ierr=-2
	  return
        end if
	pl=log10(pt(1))
	pl=addvar(1,n)
	y(5,n)=xhcor
	yhcor=1.d0-xhcor-zh(n)
c
c  reset log f (log rho) to keep log p fixed
c
	call eqstp(pl,tl,xhcor,yhcor,zh(n),nosd,nosd,fl,nit)
	y(2,n)=fl
      end do
      compc(1)=xhcor
      return
      end
      subroutine set_neutralx(x,y,pl,zh,fl)
c
c  sets value of X in y(5) such that the point is convectively neutral,
c  at fixed log(p) = pl. Also resets fl (log f or log rho) to the required
c  value
c
c  Original version: 27/1/06
c
      implicit double precision(a-h, o-z)
      dimension y(*)
      common/eqstd/ xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp,gmm1(4)
      common/ln10/ amm
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      xhh=y(5)
      fl=y(2)
      drad=fdradp(x,y,pl,zh,fl,ak,akr,akt,akx)
      ddrad=drad-dad(1)
      it=0
c
c  start iterating to set new X
c
   20 it = it + 1
      rhopx = (rho(4) - pt(4)*rho(2)/pt(2))*amm*xhh
      akpx = akx + akr*rhopx
      if(istdpr.gt.0) write(istdpr,'(a,i3,5f11.7)') 
     *  '#D# it, xhh, kappa, akpx, y(4), ddrad', 
     *  it, xhh, ak, akpx, y(4), ddrad
      dxhh=-ddrad*xhh/(drad*akpx)
      if(abs(dxhh).lt.1.d-7) then
      else if(it.gt.10) then
        if(istdpr.gt.0) write(istdou,'(/
     *    '' ***** Warning. Excessive number of iterations'',
     *    '' in s/r set_neutralx''/
     *    ''       n, log(q), log(f), log(T), X ='',i5, 1pe13.5,
     *    0p3f10.5)') n, x, y(2), y(3), xhh
	if(istdpr.gt.0.and.istdou.ne.istdpr) write(istdpr,'(/
     *    '' ***** Warning. Excessive number of iterations'',
     *    '' in s/r resscn''/
     *    ''       n, log(q), log(f), log(T), X ='',i5, 1pe13.5,
     *    0p3f10.5)') n, x, y(2), y(3), xhh
      else
        xhh=xhh + dxhh
        y(5)=xhh
        drad=fdradp(x,y,pl,zh,fl,ak,akr,akt,akx)
        ddrad=drad-dad(1)
        go to 20
      end if
      y(2)=fl
      return
      end
      double precision function comp_aver(qmix,xh,ixh,x,nn)
c
c  sets average abundance in xh out to qmix
c  
c  Original version: 17/3/06
c
      implicit double precision (a-h, o-z)
      dimension xh(ixh,*),x(*)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      if(istdpr.gt.0) write(istdpr,*)
     *  '#D# Entering comp_aver with qmix =',qmix
      xint=0.d0
      xhp=xh(1,nn)
      qp=0
      istop=0
      do n=nn,1,-1
        q=10.d0**x(n)
        if(q.ge.qmix) then
c
c  interpolate to last point
c
	  xhh=xhp+(xh(1,n)-xhp)*(qmix-qp)/(q-qp)
	  q=qmix
	  istop=1
c
	else
c
          xhh=xh(1,n)
        end if
        xint=xint+0.5d0*(q-qp)*(xhh+xhp)
        xhp=xhh
        qp=q
	if(istop.eq.1) go to 40
      end do
c
   40 comp_aver=xint/qmix
      return
      end
      subroutine test_kappa(x,yn,pl,zh)
c
c  run scan of opacity in log T, log p and X
c
c  Original version: 10/4/06
c
      implicit double precision(a-h, o-z)
      dimension yn(*)
      common/eqstd/ xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp,gmm1(4)
      common/ln10/ amm
      common/noiter/ iter, ntime, epscon, eamcon, it_force
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      data init/0/
c
c  steps in log T, log p and X
c
      nstp_t=5
      nstp_p=5
      nstp_xh=5
      d_logt=3.d-5
      d_logp=8.d-4
      d_xh=0.002
      tl0=yn(3)
      xh0=yn(5)
c
      if(init.eq.0) then
        init=1
c..        open(98,file='ttt.kappa',status='unknown')
c..        write(98,*) 2*nstp_t+1, 2*nstp_p+1, 2*nstp_xh+1
      end if
      pl0=pl
      do i=-nstp_t,nstp_t
	yn(3)=tl0+d_logt*i
        do j=-nstp_p,nstp_p
	  pl=pl0+d_logp*j
          do k=-nstp_xh,nstp_xh
	    yn(5)=xh0+d_xh*k
            drad=fdradp(x,yn,pl,zh,fl,ak,akr,akt,akx)
c..	    write(98,'(2i5,1p8e15.7)')
c..     *        ntime, iter, yn(3),pl,yn(5),rho(1),ak,akr,akt,akx
          end do
        end do
      end do
      return
      end
      subroutine mrefine(x,y,in,in1,iy,nvar,nn,xrefine,krefine,dxcrit,
     *  dxwdth,ifixnn,irsmsh)
c
c  adds meshpoints in the region corresponding to a growing 
c  convective core.
c
c  This is implemented as adding meshpoints near cricital points
c  defined by xrefine(k), k = 1,krefine.
c  This could, for example, be qmxcor.
c  This is to be used if the core jumps a substantial distance
c  as part of the evolution
c  dxcrit specifies the spacing in x at a critical point.
c  The variation at other points is determined roughly in accordance
c  with the mesh-setting procedure.
c
c  If ifixnn = 1, reset mesh to keep the number of points unchanged.
c
c  Original version 8/6/06
c
c  Modified 8/2/08, to allow remeshing near several critical points
c  Now also resets nmxcor, frmxc
c
      implicit double precision(a-h,o-z)
      include 'engenr.cz.d.incl'
      dimension x(*), y(iy,*),xrefine(*)
      dimension xst(nnmax), x_new(nnmax), xi(nnmax), yst(ivarmx,nnmax,2)
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor,nmxcp,qmxcp,rmxcp,
     *  qmxmin,qmxmax,qmxscn,qmxscp,nmxscn,nmxscp,icjump,iqcfrz
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      idiag=0
      if(istdpr.le.0) idiag=0
c
      if(istdpr.gt.0) write(istdpr,'(/
     *  '' Entering s/r mrefine with dxcrit, dxwdth ='',1p2e13.5/)') 
     *  dxcrit, dxwdth
c
      irsmsh=0
c
      n0=1
      n1=0
      do n=1,nn
	dxmax=1.d10
	if(n.gt.1) then
	  do k=1,krefine
	    xrf=xrefine(k)
	    dxcvc=min(abs(x(n)-xrf),abs(x(n-1)-xrf))/dxwdth
	    if(dxcvc.le.10.d0) then
	      dxmax=min(dxmax,dxcrit*exp(dxcvc*dxcvc))
	    end if
	  end do
	end if
	n_add=0
	if(n.gt.1.and.x(n-1)-x(n).gt.dxmax) then
c
c  add extra points between n-1 and n
c
	  if(istdpr.gt.0) write(istdpr,*) 
     *      '#D# n-1, n1, x, q, X', n-1, n1, x(n-1), 10.d0**x(n-1),
     *       y(in1+4,n-1)
	  n_add=(x(n-1)-x(n))/dxmax
	  dx=(x(n-1)-x(n))/(n_add+1)
	  do k=1,n_add
	    n1=n1+1
	    xst(n1)=x(n-1)-k*dx
c
c  if not fixing number of points, interpolate here
c
	    if(ifixnn.ne.1) then
	      call lir1(xst(n1),x,yst(1,n1,1),y(in1,1),nvar,iy,nn,n0,
     *          inter)
	      call lir1(xst(n1),x,yst(1,n1,2),y(in ,1),nvar,iy,nn,n0,
     *          inter)
	      n0=n0+1
	    end if
          end do
        end if
c
c  store model point
c
        n1=n1+1
	xst(n1)=x(n)
	do i=1,nvar
	  yst(i,n1,1)=y(in1-1+i,n)
	  yst(i,n1,2)=y(in -1+i,n)
        end do
	if(idiag.gt.0.and.(n_add.gt.0.or.n1.gt.n))
     *    write(istdpr,*) 
     *      '#D# n  , n1, x, q, X', n, n1, xst(n1), 10.d0**xst(n1),
     *       yst(5,n1,1)
c
      end do
      nn_new=n1
      if(nn_new.eq.nn) then
	if(istdpr.gt.0) write(istdpr,'(/
     *    '' No points added in mrefine''/)')
      else
	if(istdpr.gt.0) write(istdpr,'(/'' Points added in mrefine.''
     *    '' Old nn ='',i5,'' new nn ='',i5/)') nn, nn_new
c
c  if ifixnn = 1, reset new mesh to keep constant number of
c  meshpoints
c
	if(ifixnn.ne.1) then
	  do n=1,nn_new
	    x_new(n)=xst(n)
          end do
	else
	  if(istdpr.gt.0) write(istdpr,'(/
     *      '' Move to mesh with original number of points'')')
	  do n=1,nn_new
	    xi(n)=(n-1.d0)/(nn_new-1.d0)
          end do
	  do n=1,nn
	    xx=(n-1.d0)/(nn-1.d0)
	    call lir1(xx,xi,x_new(n),xst,1,1,nn_new,n,inter)
          end do
        end if
	  
c
c  move variables
c
	call move_vars(nn, nn_new, x, x_new)
c
c  test for fixed number of meshpoints
c
	if(ifixnn.eq.1) then
	  call move_var(nn, nn, x, x_new, y(in1,1),nvar,iy,yst)
	  call move_var(nn, nn, x, x_new, y(in ,1),nvar,iy,yst)
	  do n=1,nn
	    x(n)=x_new(n)
          end do
	else
	  do n=1,nn_new
	    x(n)=x_new(n)
	    do i=1,nvar
	      y(in1-1+i,n)=yst(i,n,1)
	      y(in -1+i,n)=yst(i,n,2)
	    end do
	  end do
	  nn=nn_new
c
	end if
c
c  reset nmxcor (for interpolation)
c
        qlmxcr=log10(qmxcor)
        qlmxcp=log10(qmxcp)
	nmxcor_new=0
	nmxcp_new=0
	do n=1,nn
	  if(x(n).le.qlmxcr.and.nmxcor_new.eq.0) then
	    q1=10.d0**x(n-1)
	    q2=10.d0**x(n)
	    nmxcor_new=n
	    frmxc=(q1-qmxcor)/(q1-q2)
          end if
	  if(x(n).le.qlmxcp.and.nmxcp_new.eq.0) nmxcp_new=n
        end do
	nmxcor=nmxcor_new
	nmxcp=nmxcp_new
	if(istdpr.gt.0) write(istdpr,'(/'' Model on new mesh:''/
     *    '' n, x, y(1-4):''/(i5,1p6e12.4))')
     *    (n,x(n),(y(in1-1+i,n),i=1,5),n=1,nn)
c
c  set flag for reset mesh
c
	irsmsh=1
      end if
      return
      end
      subroutine move_vars(nn, nn_new, x, x_new)
c
c  Moves various arrays originally on mesh x(1-nn) to 
c  mesh x_new(1-nn_new)
c
c  Original version 9/6/06
c
      implicit double precision(a-h,o-z)
      include 'engenr.cz.d.incl'
c
c  note that istore must be at least as big as icvvar, ivarmx1 and icvrmx
c
      parameter(ivrmx1 = ivarmx+1, istore=ivrmx1, iymix=ivarmx,
     *  nspcm2=2*nspcmx)
      dimension x(*), x_new(*)
      dimension store(istore,nnmax)
      common/rhcn/ a1,a2,a3,a4,zdummy,nvar,ifwrt,irhtst,inentc,
     *  ii1,ii2,ii3,icomp
      common/cymixc/ ymix(iymix,nnmax)
      common/cymixp/ ymixp(iymix,nnmax)
      common/ccmpor/ cmporg(nspcmx,nnmax)
      common/caddvr/ addvar(5,ngmax)
      common/convvr/ cvvarf(icvvar,6), cvvarl(icvvar,6),
     *  cvvars(icvvar,nnmax)
      common/crxcor/ cqc, cqcp, crxmn(nspcmx), crxmnp(nspcmx),
     *  crxxh_sc(nnmax)
      common/heavy/ zatmos, zhc, zh(1)
      common/compvr/ cvr(icvrmx,1)
      common/enggrv/ epsg(nnmax)
      common/crxstr/ rxstr(nspcm2,1)
      common/yprtst/ yprt(ivrmx1,nnmax)
      common/c_tstcon/ yptst(ivarmx,nnmax)
      common/cxhder/ xhder(2,nnmax)
      common/comgrt/ isprot, irotcn, a5, riner0, riner, omgrt0,
     *  omgrot(nnmax)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      call move_var(nn, nn_new, x, x_new, zh, 1, 1, store)
      call move_var(nn, nn_new, x, x_new, epsg, 1, 1, store)
      call move_var(nn, nn_new, x, x_new, crxxh_sc, 1, 1, store)
      call move_var(nn, nn_new, x, x_new, omgrot, 1, 1, store)
      call move_var(nn, nn_new, x, x_new, addvar, 5, 5, store)
      call move_var(nn, nn_new, x, x_new, xhder, 2, 2, store)
      call move_var(nn, nn_new, x, x_new, cvvars, icvvar, icvvar, 
     *  store)
      call move_var(nn, nn_new, x, x_new, yprt, ivrmx1, ivrmx1, 
     *  store)
      call move_var(nn, nn_new, x, x_new, yptst, ivarmx, ivarmx, 
     *  store)
      call move_var(nn, nn_new, x, x_new, cvr, icvrmx, icvrmx, 
     *  store)
      call move_var(nn, nn_new, x, x_new, ymix, iymix, iymix, 
     *  store)
      call move_var(nn, nn_new, x, x_new, ymixp, iymix, iymix, 
     *  store)
      call move_var(nn, nn_new, x, x_new, cmporg, nspcmx, nspcmx, 
     *  store)
      call move_var(nn, nn_new, x, x_new, rxstr, icomp, nspcm2, store)
      call move_var(nn, nn_new, x, x_new, rxstr(icomp+1,1), icomp,
     *  nspcm2, store)
      return
      end
      subroutine move_var(nn, nn_new, x, x_new, a, ii, ia, w)
c
c  interpolates array a(1-ii,1-nn) to new mesh x_new(1-nn_new)
c  w is a scratch array of size at least w(ia,nn)
c
c  Original version: 9/6/06
c
      implicit double precision(a-h,o-z)
      dimension x(*), x_new(*), a(ia,*), w(ia,*)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      do i=1,ii
        do n=1,nn
	  w(i,n)=a(i,n)
        end do
      end do
      do n=1,nn_new
	call lir1(x_new(n),x,a(1,n),w(1,1),ii,ia,nn,n,inter)
      end do
      return
      end
      subroutine ddr_ref(x,y,in,in1,iy,nvar,nn,dxcrit,irsmsh)
c
c  Analyse behaviour of ddr = nabla_R - nabla_ad, as set in
c  addvar(5,.), and refine mesh at critical points
c
c  Original version: 8/2/08
c
      implicit double precision (a-h, o-z)
      include 'engenr.cz.d.incl'
c
      parameter(naztmx = nspcmx+3, 
     *  kaxst = 3+naztmx, nbmax = nspcmx+2)
      parameter(nspcm2=2*nspcmx,iymix=ivarmx)
c
      dimension x(*), y(iy,*)
      dimension xrefine(10)
      common/caddvr/ addvar(5,ngmax)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      if(istdpr.gt.0) write(istdpr,'(/'' Enter ddr_ref''/)')
c
c  locate zero crossing in nabla_R - nabla_ad, starting from centre
c
      do n=nn,1,-1
	if(addvar(5,n).le.0) go to 15
      end do
c
   15 if(istdpr.gt.0) write(istdpr,'(/'' zero crossing found at n ='',
     *  i5,'' q ='',1pe13.5)') n, 10.d0**x(n)
c
      nz=n
c
c  locate next local maximum in nabla_R - nabla_ad near edge 
c  of core.
c  Limit test to region just beyond mixed core.
c
      nmax=0
      nmix=nz
      xhc=y(in1+4,nn)
      do n=nz,1,-1
	if(addvar(5,n-1).le.addvar(5,n).and.
     *     addvar(5,n+1).le.addvar(5,n)) then
	  nmax=n
	  go to 25
        end if
	if(nmix.eq.nz.and.y(in1+4,n)-xhc.le.1.d-3) nmix=n
	if(n.ge.nmix+10) go to 25
      end do
c
   25 if(istdpr.gt.0.and.nmax.gt.0)
     *  write(istdpr,'(/'' local maximum found at n ='',
     *  i5,'' q ='',1pe13.5/
     *  '' nabla_r - nabla_ad ='',e13.5)') n, 10.d0**x(n), addvar(5,n)
c
c  prepare for mesh refinement
c
      xrefine(1)=x(nz)
      if(addvar(5,nmax).ge.-1.d-4) then
	krefine=2
      else
c
c  test for extrapolation to zero
c
	if(addvar(5,nmax+1)/(addvar(5,nmax+1)-addvar(5,nmax)).le.4)
     *    then
	  krefine=2
        else
	  krefine=1
        end if
      end if
      if(krefine.eq.2) xrefine(2)=x(nmax)
c
      if(istdpr.gt.0) write(istdpr,'(/'' Refine near x ='',2f10.6)')
     *  (xrefine(k),k=1,krefine)
c
c  set limited range of mesh increase in vicinity of cricital points
c 
      dxwdth=x(nz-1)-x(nz+1)
      if(krefine.eq.2) dxwdth=max(dxwdth,x(nmax-1)-x(nmax+1))
c
      ifixnn=1
      call mrefine(x,y,in,in1,iy,nvar,nn,xrefine,krefine,dxcrit,
     *  dxwdth,ifixnn,irsmsh)
c
      if(istdpr.gt.0) write(istdpr,'(/'' Exit ddr_ref''/)')
      return
      end
      subroutine loc_qmxcor(x,y,iy,nn,itcmp,iosc_frz,icry)
c
c  locates boundary of convective core through suirable extrapolation
c  based on nabla_R - nabla_ad as stored in addvar(4,.)
c
c  If iosc_frz .gt. 0 only sets maximum of nabla_R - nabla_ad into
c  drmxit
c
c  Original version: 25/2/08
c
      implicit double precision(a-h, o-z)

c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
c  Note: engenr.n.d.incl replaced by engenr.nnz.d.incl, 6/6/02
c  Note: engenr.nnz.d.incl replaced by engenr.cz.d.incl, 8/1/07
      include 'engenr.cz.d.incl'
c
      dimension x(*), y(iy,*)
      common/caddvr/ addvar(5,ngmax)
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor,nmxcp,qmxcp,rmxcp,
     *  qmxmin,qmxmax,qmxscn,qmxscp,nmxscn,nmxscp,icjump,iqcfrz
      common/cmxcit/ qmxcit(nitmax),rmxcit(nitmax),frmxit(nitmax),
     .  nmxcit(nitmax),xmxcit(nitmax),frmrit(nitmax),nmxrit(nitmax),
     .  drmxit(nitmax),qmxnit(nitmax),itcmpc
      common/noiter/ iter, ntime, epscon, eamcon, it_force
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      idiag=1
      if(istdpr.le.0) idiag=0
c
      icry=0
      drmxit(itcmp)=-1.d10
c
c  find range in n where addvar(4,.) is set
c
      n1=0
      n2=0
      do n=1,nn
	if(n1.eq.0.and.addvar(4,n).ge.-1.d5) n1=n
	if(addvar(4,n).ge.-1.d5) n2=n
      end do
c
c  test for error
c
      if(n1.eq.0.or.n2.eq.0.or.n2-n1.le.2) then
	write(istdou,'(//'' ***** Error in loc_qmxcor. n1, n2 ='',2i5)')
     *    n1, n2
	if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdou,'(//
     *   '' ***** Error in loc_qmxcor. n1, n2 ='',2i5)') n1, n2
        icry=-1
	return
      end if
c
c  search for the first local maximum in nabla_R - nabla_ad
c
      nmax=0
      do n=n1+1,n2-1
	if(addvar(4,n).ge.addvar(4,n-1).and.
     *    addvar(4,n).ge.addvar(4,n+1).and.nmax.eq.0) nmax=n
      end do
c
c  test for no maximum
c
      if(nmax.eq.0) then
	if(istdpr.gt.0) write(istdpr,'(/
     *    '' No maximum found in loc_qmxcor. Keep existing value''/)')
	return
      else 
	drmxit(itcmp)=addvar(4,nmax)
	if(addvar(4,nmax).ge.0.d0) then
c
c  as a test, try to determine a better qmxcor by interpolating
c  Note that this ddrmax corresponds to the _previous_ qmxcor,
c  in an appropriate (and possibly relevant) sense
c
	   if(itcmp.ge.4) then
c
c  find latest negative value
c
	     qmxcor_tst=0.d0
	     do i=itcmp-1,2,-1
	       if(istdpr.gt.0) write(istdpr,*) 
     *           '#D# itcmp, i, drmxit(i),qmxcit(i)',
     *           itcmp, i, drmxit(i),qmxcit(i)
	       if(drmxit(i).lt.0.and.drmxit(i).gt.-1.d5) then 
c..		 qmxcor_tst=
c..     *            (qmxcit(i-1)*drmxit(itcmp)-qmxcit(itcmp-1)*drmxit(i))/
c..     *             (drmxit(itcmp)-drmxit(i))
c
c  try setting qmxcor_test with undercorrection factor ucor
c
		 ucor=0.8d0
		 qmxcor_tst=qmxcit(i-1)+ucor*
     *            (qmxcit(i-1)-qmxcit(itcmp-1))*drmxit(i)/
     *             (drmxit(itcmp)-drmxit(i))
		 go to 25
               end if
             end do
	  end if
   25     if(istdpr.gt.0) write(istdpr,'(/
     *      '' Maximum found in loc_qmxcor is positive.'',
     *      '' Keep existing value''/
     *      '' qmxcor_tst ='',1pe13.5/)') qmxcor_tst
c..	  if(ntime.eq.41) then
	  if(ntime.ge.0.and.iosc_frz.eq.0) then
	    if(istdpr.gt.0) write(istdpr,'(//
     *       '' Use qmxcor_tst after all''/)')
            qmxcor=qmxcor_tst
	  end if
	  return
        end if
      end if
c
      if(iosc_frz.gt.0) return
c
c  possible diagnostic output in vicinity of maximum
c
      if(idiag.gt.0) write(istdpr,'(/'' Values near local maximum''/
     *  '' n, x, q, nabla_r - nabla_ad:''/(i5,0pf10.5,1p2e13.5))')
     *  (n, x(n), 10.d0**x(n), addvar(4,n),n=nmax-4,nmax+4)
c
c  Set qmxcor from extrapolation from maximum
c
      nmxcor=nmax
      frmxc=addvar(4,nmax-1)/(addvar(4,nmax-1)-addvar(4,nmax))
      qmxcor_new=frmxc*10.d0**x(nmax)+(1-frmxc)*10.d0**x(nmax-1)
      if(istdpr.gt.0) write(istdpr,'(/
     *  '' In loc_qmxcor, qmxcor reset from'',1pe13.5,'' to'',e13.5/
     *  '' nxmcor, frmxc ='',i5,0pf10.5)')
     *  qmxcor, qmxcor_new, nmxcor, frmxc
      qmxcor=qmxcor_new
      return
      end
      subroutine loc_qmxcor_n(x,y,iy,nn,itcmp,iosc_frz,icry,qmxcor_new)
c
c  locates boundary of convective core through suirable extrapolation
c  based on nabla_R - nabla_ad as stored in addvar(4,.)
c
c  If iosc_frz .gt. 0 only sets maximum of nabla_R - nabla_ad into
c  drmxit
c
c  Original version: 25/2/08
c
      implicit double precision(a-h, o-z)

c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
c  Note: engenr.n.d.incl replaced by engenr.nnz.d.incl, 6/6/02
c  Note: engenr.nnz.d.incl replaced by engenr.cz.d.incl, 8/1/07
      include 'engenr.cz.d.incl'
c
      dimension x(*), y(iy,*)
      dimension ddrad(nnmax)
      common/eqstd/ xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp,gmm1(4)
      common/heavy/ zatmos, zhc, zh(1)
      common/caddvr/ addvar(5,ngmax)
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor,nmxcp,qmxcp,rmxcp,
     *  qmxmin,qmxmax,qmxscn,qmxscp,nmxscn,nmxscp,icjump,iqcfrz
      common/cmxcit/ qmxcit(nitmax),rmxcit(nitmax),frmxit(nitmax),
     .  nmxcit(nitmax),xmxcit(nitmax),frmrit(nitmax),nmxrit(nitmax),
     .  drmxit(nitmax),qmxnit(nitmax),itcmpc
      common/noiter/ iter, ntime, epscon, eamcon, it_force
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      idiag=1
      if(istdpr.le.0) idiag=0
c
      icry=0
      drmxit(itcmp)=-1.d10
c
c  find range in n for setting nabla_R - nabla_ad to search for
c  maximum and extrapolate
c
      n1=nn+1
      n2=0
      do it=1,itcmp-1
	nmx=nmxcit(it)
	if(nmx.gt.0) then
	  n1=min(n1,nmxcit(it))
	  n2=max(n2,nmxrit(it))
	end if
      end do
      n1=max(1,n1-5)
      n2=min(nn,n2+5)
c
      do n=n1,n2
	drad=fdrad(x(n),y(1,n),zh(n),ak,akr,akt,akx)
	ddrad(n)=drad-dad(1)
      end do
c
c  test for error
c
      if(n1.eq.0.or.n2.eq.0.or.n2-n1.le.2) then
	write(istdou,'(//'' ***** Error in loc_qmxcor_n. n1, n2 ='',
     *    2i5)') n1, n2
	if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdou,'(//
     *   '' ***** Error in loc_qmxcor_n. n1, n2 ='',2i5)') n1, n2
        icry=-1
	return
      end if
c
c  search for the first local maximum in nabla_R - nabla_ad
c
      nmax=0
      do n=n1+1,n2-1
	if(ddrad(n).ge.ddrad(n-1).and.
     *    ddrad(n).ge.ddrad(n+1).and.nmax.eq.0) nmax=n
      end do
c
c  test for no maximum
c
      if(nmax.eq.0) then
	if(istdpr.gt.0) write(istdpr,'(/
     *    '' No maximum found in loc_qmxcor. Keep existing value''/)')
	return
      else 
	drmxit(itcmp)=ddrad(nmax)
	if(ddrad(nmax).ge.0.d0) then
c
c  as a test set qmxcor_new to maximum
c
	  qmxcor_new=10.d0**x(nmax)
	  qmxnit(itcmp)=qmxcor_new
c
   25     if(istdpr.gt.0) write(istdpr,'(/
     *      '' Maximum found in loc_qmxcor is positive,'',1pe13.5,
     *      '' Keep existing value''/
     *      '' qmxcor_new ='',e13.5/)') ddrad(nmax), qmxcor_new
	  if(ntime.ge.0.and.iosc_frz.eq.0) then
	    if(istdpr.gt.0) write(istdpr,'(//
     *       '' Use qmxcor_new after all''/)')
            qmxcor=qmxcor_new
	  end if
	  return
        end if
      end if
c
c  possible diagnostic output in vicinity of maximum
c
      if(idiag.gt.0) write(istdpr,'(/'' Values near local maximum''/
     *  '' n, x, q, nabla_r - nabla_ad:''/(i5,0pf10.5,1p2e13.5))')
     *  (n, x(n), 10.d0**x(n), ddrad(n),n=nmax-4,nmax+4)
c
c  Set qmxcor from extrapolation from maximum
c
      nmxcor=nmax
      frmxc=ddrad(nmax-1)/(ddrad(nmax-1)-ddrad(nmax))
      qmxcor_new=frmxc*10.d0**x(nmax)+(1-frmxc)*10.d0**x(nmax-1)
      qmxnit(itcmp)=qmxcor_new
      if(iosc_frz.eq.0) then
        if(istdpr.gt.0) write(istdpr,'(/
     *    '' In loc_qmxcor, qmxcor reset from'',1pe13.5,'' to'',e13.5/
     *    '' nxmcor, frmxc ='',i5,0pf10.5)')
     *    qmxcor, qmxcor_new, nmxcor, frmxc
	qmxcor=qmxcor_new
      else
        if(istdpr.gt.0) write(istdpr,'(/
     *    '' In loc_qmxcor, qmxcor_new ='',1pe13.5,''. No reset''/
     *    '' nxmcor, frmxc ='',i5,0pf10.5)')
     *    qmxcor_new, nmxcor, frmxc
      end if
      return
      end
