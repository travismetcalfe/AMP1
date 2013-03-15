      subroutine cmpcvc(x,y,in,in1,iy,nn,compc,icomp,iche3,icnocs,
     *  iheccs,dt,it,iprdcr,icrycm)
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
c  Modified 15/9/03, introducing option imixc4 .ge. 1 to suppress
c  resetting of composition of diffusing elements in this routine 
c  in the case of diffusion.
c
c  *** Note: in old version z was unset before call of he3abd for
c      core He3 abundance.
c
      implicit double precision(a-h,o-z)
      include 'engenr.nnz.d.incl'
      logical norct, nosd
c
      parameter(naztmx = nspcmx+3, 
     *  kaxst = 3+naztmx, nbmax = nspcmx+2)
      parameter(nitmax=50,nspcm2=2*nspcmx,iymix=ivarmx)
c
      dimension x(1), y(iy,1), compc(1)
      dimension rxmn(nspcmx), drxmn(nspcmx), rxmnp(nspcmx),
     .  drxmnp(nspcmx),
     .  alrhmn(krnrmx), alrhmp(krnrmx), alrhmc(krnrmx), 
     .  qmxcit(nitmax),rmxcit(nitmax),frmxit(nitmax),nmxcit(nitmax),
     *  cmporg(nspcmx,nnmax),
     .  dcmpgr(nspcmx),x3new(4),alres(nnmax),cmpqcp(nspcmx)
      common/compvr/ cvr(icvrmx,1)
      common/cstcvr/ icvh1, icvhe4, icvhe3, icvc12, icvc13, icvn14, 
     *  icvo16
      common/cstyst/ iyh1, iyhe3, iyc12, iyc13, iyn14, iyo16, 
     *  iyche4, iycc12
      common/cstazt/ iah1, iahe3, iac12, iac13, ian14, iao16, 
     *  iache4, iacc12
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,jnc,frcf(6),frcl(6)
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor,nmxcp,qmxcp,rmxcp,
     *  qmxmin,qmxmax
      common/crxcor/ cqc, cqcp, crxmn(nspcmx), crxmnp(nspcmx)
      common/cymixc/ ymix(iymix,nnmax)
      common/cmxcnt/ ddrmix, inmixc, imixc0, imixc1, imixc2, imixc3,
     *  imixc4
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
      common/noiter/ iter, ntime
c
c  dataset designators
c
      common/cdsdsg/ idstrl,idshvz,idszab,idsevl,idssum,idscen,idsnov,
     *  idsoov,idsgng,idsgsm,idstm1,idstm2,iducen
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      if(istdpr.gt.0) write(istdpr,'(/a)') 'Entering cmpcvc'
c
c  set full number of composition variables, with diffusion
c
      if(idiffus.gt.0) then
	ifcomp = iccomp+idcomp
      else
	ifcomp = icomp
      end if
c
      if(istdpr.gt.0) write(istdpr,*) 
     *  'In cmpcvc, ifcomp, imixc3 =',ifcomp,imixc3
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
c
      if(istdpr.gt.0) write(istdpr,*) 
     *  'testing for predictor with iprdcr =',iprdcr
c
c  test for predictor or corrector
c
      if(iprdcr.eq.1) then
	if(istdpr.gt.0) write(istdpr,'(/a)') 'Start predictor section'
        if(frmxc.gt.1.and.nmxcor.gt.0) then
          nmxcp=nmxcor+1
        else
          nmxcp=nmxcor
        end if
        qmxcp=qmxcor
        rmxcp=rmxcor
        if(istdpr.gt.0) write(istdpr,*) 
     *    'nmxcp, qmxcp =',nmxcp, qmxcp
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
          return
        else if(iheccs.ne.0.and.y(in+2,nn).gt.7.778.and.
     *    (y(in+iyche4-1,nn).lt.1.e-4.or.
     *    y(in+iyche4-1,nn).lt.1.e-3.and.qmxcp.lt.0.02)) then
          if(istdpr.gt.0) write(istdpr,106) qmxcp,y(in+iyche4-1,nn)
          nmxcor=-1
          qmxcor=0.
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
     *    alres)
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
          if(frmxc.gt.1) then
            nc1=nmxcor+1
          else
            nc1=nmxcor
          end if
	  if(istdpr.gt.0) write(istdpr,*) 
     *      'In cmpcvc, predictor step, set composition with nc1 =',
     *      nc1
          do 25 i=1,icomp
          if(i.eq.2.and.iche3.eq.1) then
            compc(i)=x3new(1)
          else
            compc(i)=max(y(in+3+i,nn)+dt*rxmnp(i),1.d-10)
            if(istdpr.gt.0) write(istdpr,*) 
     *        'predicted compc(i) set to',compc(i)
          end if
   25     continue
c
c  Store predicted composition in y(in1 + .,.)
c  unless for diffusing elements, in case with diffusion, 
c  and imixc4 .ge. 1
c
	if(idiffus.lt.1.or.imixc4.eq.0) then
	  icomp1=1
	  icomp2=ifcomp
        else
	  icomp1=2
	  icomp2=1+iccomp
        end if
	if(icomp2.ge.icomp1) then
	  do n=nc1,nn
	    do i=icomp1,icomp2
              y(in1+3+i,n)=compc(i)
            end do
          end do
        end if
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
c  end of predictor step
c
      else
c
c  corrector step during iteration
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
c..        write(38,'(4i5,1p2e13.5)') 
c..     *    (ntime,iter,0,n,x(n),y(in1+4,n),n=nmxcor-10,nmxcor+10)
        xhintp=xhint
c..      do n=400,430
c..        write(6,'(i5,1p10e13.5)') n,(rxstr(i,n),i=1,nspcm2)
c..      end do
c
c  set model, including composition, for setting of boundary
c  of mixed region
c
c..     write(6,*) 'n, q, Xp, change, new X:'
        do 34 n=1,nn
        do 33 i=1,4
   33   ymix(i,n)=y(in1+i-1,n)
c..     if(iabs(n-nmxcp).le.10.or.iabs(n-nmxcor).le.5) 
c..     *    write(6,'(i5,1p4e13.5)')
c..     *    n, 10.d0**x(n), y(in+4,n), dt*rxstr(icomp+1,n),
c..     *    y(in+4,n)+dt*rxstr(icomp+1,n)
        do 34 i=1,icomp
   34   ymix(4+i,n)=max(y(in+3+i,n)+dt*rxstr(icomp+i,n),1.d-10)
c
c  Set type of setting unstable boundary.
c  Modified 28/12/00 for more clarity, but consistent with
c  older version. Might benefit from a little further thought.
c
        if(it.lt.10) then
c
c  set standard case, depending on imixc2
c
          if(imixc2.eq.0) then
            iextrp = 3
          else
            iextrp = 1
          end if
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
          if(nmxcit(it-1).eq.0.and.nmxcit(it-2).eq.0) then
            nmxcor=-1
          else if(nmxcit(it-1).eq.0.or.nmxcit(it-2).eq.0) then
            nmxcor=max(nmxcit(it-1),nmxcit(it-2))
          else
            nmxcor=0.5*(nmxcit(it-1)+nmxcit(it-2))
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
        else if(it.ge.12) then
c
c  set iextrp = -1, to freeze boundary completely
c
          iextrp=-1
c
c  for it = 12, set frozen interpolation weight
c
          if(it.eq.12) then
	    frmxc=min(frmxit(it-2),frmxit(it-1))
            if(frmxc.le.0) frmxc=0
            if(istdpr.gt.0) write(istdpr,*) 
     *        'Freeze interpolation weight at frmxc =',frmxc
          else
            if(istdpr.gt.0) write(istdpr,*) 
     *        'Interpolation weight frozen at frmxc =',frmxc
	  end if
        else if(it.ge.10) then
          if(iextrp.eq.3) then
            if(nmxcit(it-1).lt.nmxcit(it-2)) then
              itfrz=it-1
            else
              itfrz=it-2
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
c  test for basing test for mixed region on non-zero
c  gradient difference
c
        if(ddrmix.gt.0.and.iextrp.eq.3) iextrp=10
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
c
        if(iextrp.ne.-1) then
          if(istdpr.gt.0) 
     *      write(istdpr,*) 'call mixcor from cmpcvc, iextrp =',iextrp
          call mixcor(x,ymix,iymix,nn,compc,iextrp,nmxfrz)
          frmxit(it)=frmxc
          nmxcit(it)=nmxcor
          qmxcit(it)=qmxcor
          rmxcit(it)=rmxcor
        end if
c
c  test for presence of mixed core
c
        if(nmxcor.eq.0) then
          if(istdpr.gt.0) write(istdpr,110)
          return
        end if
c
c  when this is onset of convective core, store central values
c  in compc
c
        if(qmxcp.eq.0) then
          if(istdpr.gt.0) write(istdpr,112) qc
          do 37 i=1,ifcomp
   37     compc(i)=y(in1+3+i,nn)
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
            do 38 n=1,nn
            nmxcor=n
            if(x(n).le.qmxcrl) go to 39
   38       continue
c
   39       frmxc=(qmxcrl-x(nmxcor-1))/(x(nmxcor)-x(nmxcor-1))
            if(istdpr.gt.0) write(istdpr,*) 
     *        ' Reset qmxcor, nmxcor, frmxc to ',
     *          qmxcor, nmxcor, frmxc
          end if
        end if
c
c  test for setting composition change from growing mixed region
c
        if(qmxcor.gt.qmxcp) then
c
          call cmpgmx(x,y(in,1),y(in1,1),iy,iy,nn,qmxcor,qmxcp,
     *      y(in+4,nn),dcmpgr,ifcomp)
c
c  Reset composition. If qcp = 0, only set average returned in dcmpgr
c
          do 40 i=1,ifcomp
          if(istdpr.gt.0) write(istdpr,*) 
     *      'Set new compc from previous value'
          compc(i)=max(y(in+3+i,nn)+dcmpgr(i),1.d-10)
          if(istdpr.gt.0) write(istdpr,'(''Element'',i2,1p3e12.4)') 
     *      i, y(in+3+i,nn),dcmpgr(i),compc(i)
   40     continue
c
c  otherwise initialize compc and set dcmpgr to zero, for output 
c  purposes
c
	else
	  do 41 i=1,ifcomp
          compc(i)=max(y(in+3+i,nn),1.d-10)
   41     dcmpgr(i)=0.d0
c
        end if
c
c  the following statement arises because of confusion between 
c  different measures of core edge.
c
        if(qmxcp.gt.0.and.qcp.eq.0) qcp=qmxcp
c
c  set average rx for this timestep
c
        call rxhmxc(x,y(in1,1),iy,nn,compc,qc,rxmn,drxmn,icomp,alres)
c
c  reset luminosity in ymix from alres
c
        do 42 n=1,nn
   42   ymix(4,n)=alres(n)
c
c  store values at current time step in common
c
        cqc=qc
        call store(rxmn,crxmn,icomp)
c
c  reset corrected new core abundances. Possible
c  contribution from growing core already included.
c
	icrycm=0
c
        do 45 i=1,icomp
        if(i.ne.2.or.iche3.ne.1) then
          if(qmxcp.gt.0) then
            if(istdpr.gt.0) write(istdpr,*) 
     *        'Set new compc from previous value'
            thetai=theta(4+i)
            phii=1-thetai
            compc(i)=max(compc(i)+
     *        dt*(phii*rxmnp(i)+thetai*rxmn(i)),1.d-10)
            if(istdpr.gt.0) write(istdpr,'(''Element'',i2,1p5e12.4)') 
     *        i, y(in+3+i,nn),dcmpgr(i),
     *        dt*phii*rxmnp(i),dt*thetai*rxmn(i),compc(i)
          end if
c
c  test for unphysical value
c
	  if(compc(i).gt.1.d0) icrycm=-i
        end if
   45   continue
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
c..       write(6,*) 
c..     *      'Before call of he3abc in cmpcvc, xh, yh, zhc, dt, y5:'
          if(istdpr.gt.0) write(istdpr,*) 
     *      xh, yh, zhc, dt, y(in+5,nn)
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
	  return
        end if
c
        do 54 n=1,nn
        qx=10.d0**x(n)
        if(qx.gt.qmxcor) then
          do 50 i=1,ifcomp
   50     ymix(4+i,n)=cmporg(i,n)
        else
c..	  write(6,*) 'At q =',qx,' reset X to', compc(1)
          do 52 i=1,ifcomp
   52     ymix(4+i,n)=compc(i)
        end if
c
   54   continue
c
c  reset composition in region between old and new convective core
c
        if(qmxcor.ge.qmxcp) then
          if(istdpr.gt.0) write(istdpr,130) qmxcp, qmxcor
        else if(nmxcp.lt.nc1) then
c
c  reset in intermediate region. Since logics of nmxcor and nc1
c  is a little unclear, make a careful check on whether we
c  are actually in the intermediate region
c
          if(istdpr.gt.0) write(istdpr,*) 'n, qx, tfrct, ',
     *      'dt*(1-tfrct)*rxstr(2,n),tfrct*(ymix(2)-ypc(2))'
          do 56 n=1,nn
          qx=10.d0**x(n)
	  if(qx.gt.qmxcp) then
	    do i=1,icomp
	      cmpqcp(i)=y(in+3+i,n)
            end do
	  end if
	  if(qx.ge.qmxcor.and.qx.le.qmxcp) then
            tfrct=(qmxcp-qx)/(qmxcp-qmxcor)
            if(istdpr.gt.0) write(istdpr,'(i5,1p4e13.5)') 
     *        n, qx, tfrct, dt*(1-tfrct)*rxstr(2,n), 
     *        tfrct*(ymix(6,nn)-y(in+5,nn))
c
c  include test for unreasonable deviation from linear interpolation
c
            do i=1,icomp
              icp=in +3+i
              ymix(4+i,n)=y(icp,n)+dt*(1-tfrct)*rxstr(i,n)+
     *               tfrct*(ymix(4+i,nn)-y(icp,nn))
	      yintp=(1-tfrct)*cmpqcp(i)+tfrct*compc(i)
	      if(ymix(4+i,n).gt.1.5*yintp.or.
     *           ymix(4+i,n).lt.0.5*yintp) then
		if(istdpr.gt.0) write(istdpr,*) 
     *            'Error: n, i, ymix, yintp:',
     *       	  n, i, ymix(4+i,n), yintp
		ymix(4+i,n)=yintp
	      end if
	    end do
          end if
   56     continue
        end if
c
c  test for setting zhc from mixed composition
c
	if(idiffus.ge.2) zhc=compc(ifcomp)
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
c  unless for diffusing elements, in case with diffusion, 
c  and imixc4 .ge. 1
c
	if(idiffus.lt.1.or.imixc4.eq.0) then
	  qmxres=max(qmxcor,qmxcp)
c..	  write(6,*) 'Reset composition for q .le. ',qmxres
	  do 58 n=1,nn
	  qx=10.d0**x(n)
	  if(qx.le.qmxres) then
	    do 57 i=1,ifcomp
   57       y(in1+3+i,n)=ymix(4+i,n)
	  end if
   58     continue
c
c  otherwise reset composition for non-diffusing elements,
c  and do not reset composition of diffusing elements
c  but restore original values in ymix 
c  (the latter for setting of cvr, and in need of careful thought)
c
	else
	  do n=1,nn
	    ymix(5,n)=y(in1+4,n)
	    if(idcomp.gt.1) then
	      do i=2,idcomp
	        ymix(4+iccomp+i,n)=y(in1+3+iccomp+i,n)
	      end do
	    end if
	    if(iccomp.gt.0) then
	      do i=1,iccomp
	        y(in1+4+i,n)=ymix(5+i,n)
	      end do
	    end if
	  end do
        end if
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
      end if
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
      include 'engenr.nnz.d.incl'
      parameter (nnmaxc=2000)
      dimension x(1),yp(iyp,1),y(iy,1),compc(1),dcmpgr(1),
     *  cmp(nspcmx,nnmaxc), qi(nnmaxc)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'Enter cmpgmx'
      call flush(6)
c
c  find location of growing region, set composition sum
c
      n1=0
      iend=0
      qx=1
      do 30 n=1,nn
      qxp=qx
      qx=10.d0**x(n)
      nc=n
      if(qx.le.qc.and.n1.eq.0) then
c
c  first point. interpolate
c
        if(istdpr.gt.0) write(istdpr,*) 
     *    'first. n, qc, qcp, qx, qxp', n, qc, qcp, qx, qxp
        call flush(6)
        n1=1
        frct=(qc-qxp)/(qx-qxp)
	if(n.gt.1) then
          do i=1,ifcomp
            cmp(i,n1)=0.5*(frct*(yp(4+i,n)+y(4+i,n))+
     *            (1-frct)*(yp(4+i,n-1)+y(4+i,n-1)))
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
c  last point. interpolate
c
        n1=n1+1
        if(istdpr.gt.0) write(istdpr,*) 
     *    'last. n, qc, qcp, qx, qxp', n, qc, qcp, qx, qxp
        call flush(6)
        frct=(qcp-qxp)/(qx-qxp)
        do 26 i=1,ifcomp
   26   cmp(i,n1)=0.5*(frct*(yp(4+i,n)+y(4+i,n))+
     *          (1-frct)*(yp(4+i,n-1)+y(4+i,n-1)))
        qi(n1)=qcp
        go to 35
      else if(qx.lt.qc) then
        n1=n1+1
        do 28 i=1,ifcomp
   28   cmp(i,n1)=0.5*(yp(4+i,n)+y(4+i,n))
        qi(n1)=qx
c..     write(6,*) n, n1, qx
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
c
      do 40 i=1,ifcomp
   40 dcmpgr(i)=0
c
      if(n1.eq.1) then
        if(istdpr.gt.0) write(istdpr,130) qcp, qc
      else
        do 45 n=2,n1
        ww=0.5*(qi(n-1)-qi(n))
        do 45 i=1,ifcomp
   45   dcmpgr(i)=dcmpgr(i)+(qi(n-1)-qi(n))*
     *    (0.5*(cmp(i,n-1)+cmp(i,n))-compc(i))
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
