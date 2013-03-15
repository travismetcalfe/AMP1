      subroutine bciter(az,gh,dgh,acy,eam,iconv,nmax,iter,nb)
c
c  Iterates to find values and second derivatives of p, t and x
c  at the centre. furthermore sets expansion coefficients for m and l
c  and the quantities needed in s/r bcs, i. e. values of m and l at
c  the first mesh point and their derivatives
c
c  note on modifications 20/8/84, 31/7/91:
c
c  variables have been rescaled. now az(1-5) = (1.e-17*p,1.e-7*t,
c  x,x3,1.e-11*r). azt(1-4) are central values of these variables,
c  and ax(1-4) are half their second derivatives with respect to
c  1.e-11*r.
c
c  Treatment of He3: determined by variable ihe3bc which is
c  transferred in common/bccn/. If ihe3bc = 1, He3 is either assumed 
c  to be in equilibrium or given by time-dependent fudge. 
c  all calculations done so far (31/7/1991) have used this option.
c  To simplify organization we introduce a second array
c  az1(1-5), defined by az1(i)=az(i), i = 1, ..., nb, az1(nb1)=az(5).
c  nb transmits number of variables included
c  in boundary condition calculation. 
c
c  Calculation of azt(4) and ax(4) has been fixed up.
c
c  modified 4/1/1985 to set boundary conditions in terms of
c  m/1.d33 and l/1.d33.
c
c  modified 27/3/1986, to check for zero central energy generation in  
c  setting gh. 
c
c  modified 1/3/88, setting central composition derivatives to zero
c  when X .le. 2.e-3 at the innermost meshpoint. This should be
c  thought through more carefully
c  
c  modified 15/3/88, setting consistent numerical constants
c
c  modified 11/5/89 to ignore entropy term in energy equation for
c  inentr=1. this is implemented in s/r bccoef by setting sds, sdd 
c  and sdt to zero
c
c  modified 12/6/89, in treatment of abundance derivatives for
c  prescribed abundance profile. Now assumes that central values
c  are given in the correct azt(kp), and sets the second derivatives
c  into ax(kp). Previously assumed ax(kp) to be given, and set
c  azt(kp)
c  Note: problems for convective core, at "do 49" and "do 198"
c
c  modified 12/8/90, setting al01 and al21 
c  (for storage in common/ksider/) to include terms in time 
c  derivatives. Note that this has no effect on solution 
c  if theta(4) = 1, as it has been in all calculations so far,
c  but does affect setting the central mesh point when stretching a
c  trial model.
c
c  Corrected 5/8/97, by including engenr.nnz.d.incl, to set nspdmx etc
c
c  Modified 14/7/00, eliminating entropy term in energy generation
c  if it becomes unreasonably large in convective core
c  (typically as a result of sudden onset of core convection during
c  main-sequence evolution). Parameters etc. in this could do with
c  fine tuning.
c
c  ****************************************************************
c
c  Implementation problems:
c
c  #P# System: sundog.ucar.edu. 27/7/04: the line
c  if(idgbcs.ge.33333) write(6,*) rhc
c  was added to suppress problems with setting of dxt(1,1) on 
c  sundog.ucar.edu. Similar problems have not been encountered on
c  tuc47.
c
c  ****************************************************************
c
c  Major modification initiated 26/7/91, to incorporate details
c  of CNO cycle.
c
c  26/7/91: change storage of derivatives to compressed form
c  31/7/91 - : Change dimensioning of arrays to general case.
c
c  ****************************************************************
c
c  Modified 24/8/96, setting to zero derivative of X for
c  idfxbc = 1. Note that this is a fairly crude fudge,
c  to be replaced by a proper treatment of the central
c  boundary conditions with diffusion and settling.
c
c  Modified 31/10/00, introducing transformation of
c  luminosity, including a shift in the luminosity variable used
c  for computing when the central hydrogen abundance is very small
c  (hard-coded to .le. 1.e-4). The shift is managed by s/r cpydif.
c  Note that this only has an effect in the version of the programme
c  with diffusion. Otherwise (at least for now) alshft is set to 0.
c
c  Modified 10/9/02, changing setting of thtp from theta(kp) to
c  theta(kp+2), to be consistent with use of theta in other parts
c  of the calculation.
c
c  Modified 20/10/04, introducing flag for testing for errors in physics
c  routines and return with flag.
c
c  Modified 5/12/07, setting proper centralization of energy equation in
c  case with diffusion. (Centralization of composition equations still
c  needs thought.)
c
c  Modified 7/3/08, introducing icncbc .ge. 1 to reset composition derivatives
c  in convective core.
c
c  Old storage in  common/ksider/ (up to and including xt, this
c  corresponds to new case with nspcmx = 2, idthm = 3):
c  
c  1:             al01
c  2:             al21
c  3 - 7:         azt(1-5)
c  8 - 12:        ax(1-5)
c  13:            pt
c  14:            tt
c  15 - 17:       th(1-3)
c  18 - 20:       xt(1-3)
c  21 - 35:       xxh(1-5,1-3)
c  36:            alt
c  37:            altt
c  38 - 40:       alh(1-3)
c  41 - 43:       alhh(1-3)
c  44 - 47:       albb(1-3)
c  48 - 77:       dztz(1-5,1-6)
c  78 - 92:       dax(1-5,1-5)
c  93 - 97:       dpt(1-5)
c  98 - 102:      dtt(1-5)
c  103 - 111:     dth(1-3,1-3)
c  112 - 120:     dxt(1-3,1-3)
c  121 - 196:     dxh(1-5,1-5,1-3)
c  197 - 201:     dlt(1-5)
c  202 - 206:     dltt(1-5)
c  207 - 215:     dlh(1-3,1-3)
c  216 - 224:     dlhh(1-3,1-3)
c  225 - 233:     dlbb(1-3,1-3)
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
      implicit double precision (a-h,o-z)
      logical time0,nosd,notd,conv,dtest,skipt,noder,norct,
     *  dgbcs2, dgbcs3
      character cxnuc*4
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.nnz.d.incl'
c
c  set additional dimensions. Note that idthm corresponds to
c  number of independent thermodynamic variables
c
      parameter(naztmx = nspcmx+3, nbmax = nspcmx+2, 
     *  iwrk1 = nbmax, iwrk2 = 2*naztmx)
c
      dimension az(*),gh(*),dgh(2,*)
      dimension az1(naztmx),aztp(nbmax),axp(nbmax),
     .  wrk(iwrk1, iwrk2), cxnuc(nspcmx),
     .  ddrc(nbmax),ec(2,2),dec(nbmax,2,2),coef(nbcprv),coefp(nbcprv)
      common/bcprev/ al0p,al2p,aztstp(naztmx),axstp(naztmx),
     *  ptp,ttp,thp(idthm),xtp(idthm),xxhp(nbmax,nspcmx),
     *  altp,alttp,alhp(idthm),alhhp(idthm),albbp(idthm)
      common/ksider/ al01,al21,aztst(naztmx),axst(naztmx),pt,tt,
     *  th(idthm),xt(idthm),xxh(nbmax,nspcmx),alt,altt,
     *  alh(idthm),alhh(idthm),albb(idthm),
     *  dztz(nbmax,naztmx),dax(nbmax,nbmax),dpt(nbmax),dtt(nbmax),
     *  dth(idthm,idthm),dxt(idthm,idthm),
     .  dxh(nbmax,nbmax,nspcmx),dlt(nbmax),dltt(nbmax),
     *  dlh(idthm,idthm),dlhh(idthm,idthm),dlbb(idthm,idthm)
      common/caztax/ azt(nbmax),ax(nbmax)
      common/rhcn/ crhcn(5),icrhcn(3),inentr
      common/clshft/ alshft
      common/bccn/ b1,b2,b3,b4,nnb,iveb,icncbc,ihe3bc
      common/cmtime/ age, time0
      common/nmbmsh/ nn
      common/excf/ am0,al0,am2,al2,
     *  dm0(nbmax),dl0(nbmax),dm2(nbmax),dl2(nbmax)
      common/thetac/ theta(1)
      common/thetcf/ thetdf(1)
      common/step/ dt
      common/heavy/ zatmos, zhc, zh(1)
      common/logf/ fl
      common/he3int/ xhe3(4)
      common/rbcder/ rhc1,rhds(3),rhdd(3,3),sds(3),sdd(3,3),sdt(3,3,3)
      common/eqstd/dummy(54),p(20),cp(4),dad(4),dummy1(8)
      common/eqstcl/ iradp,mode,dtest,skipt,noder
      common/ln10/ amm
      common/rnratd/ al(10,krnrmx),norct
      common/engcnt/ xhzlm1, xhzlm2, nxhzer, icnocs, nspec,
     *  iheccs, nspect
      common/he3eql/ xhe3eq,ieqhe3,iequsd
      common/he3fdg/ agesh,ifdhe3
      common/compos/ xseth, yset, xset3, xset12, xset13, xset14, xset16
      common/compvr/ cvr(icvrmx, 1)
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec, idcomp, iccomp
      common/cstcvr/ icvh1, icvhe4, icvhe3, icvc12, icvc13, icvn14, 
     *  icvo16
      common/cstyst/ iyh1, iyhe3, iyc12, iyc13, iyn14, iyo16, 
     *  iyche4, iycc12
      common/cstazt/ iah1, iahe3, iac12, iac13, ian14, iao16, 
     *  iache4, iacc12
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor
      common/rnrout/ rnzt(10),rnuw(10)
      common/frecno/ rcno
      common/cntval/ rcnoc,uwc,dradc
      common/crxcor/ cqc, cqcp, crxmn(nspcmx), crxmnp(nspcmx)
      common/cdiffu/ idiffus, itbdif, rhob, rmxdif, tlfdif, dtldif,
     *  ctdfmx, rctdif, rctdf1,idfxbc
      common/noiter/ iterrh, ntime
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3
c
      common/cdgrhb/ kdgrhb
      common/cdgphs/ kdgeos, kdgopc, kdgeng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
c
      equivalence (am0,ec(1,1)),(dm0(1),dec(1,1,1)),(al01,coef(1)),
     .  (al0p,coefp(1))
c
      dgbcs2=idgbcs.ge.2.and.istdpr.gt.0
      dgbcs3=idgbcs.ge.3.and.istdpr.gt.0
c
      if(istdpr.gt.0) then
	write(istdpr,'(/'' Entering bciter''/)')
        if(iter.eq.0) call strbcs(1)
      end if
c
c  set storage parameters and nb, depending on ihe3bc and
c  icnocs
c
      call engcse(ieqhe3, icnocs, iheccs, nspec, nspxx3, nspcno,
     *  nsphec, nspect, irnrat, cxnuc)
c
c  initialize internal flag for entropy terms in energy equation
c
      inenti=inentr
c
      if(ihe3bc.eq.0) then
	nbh = 4 + nspcno
	nb  = 4 + nspcno + nsphec
      else
	nbh = 3 + nspcno
	nb  = 3 + nspcno + nsphec
      end if
      if(dgbcs2) write(istdpr,*) 'ihe3bc, nspect', ihe3bc, nspect
      if(dgbcs2) write(istdpr,*) 'nb set to',nb,' in bciter'
c
      nb1=nb+1
c
c  set original value of ieqhe3 and test for equilibrium he3
c  at centre
c
      ieqho=ieqhe3
      if(ihe3bc.eq.1.and.ifdhe3.ne.1) ieqhe3=1
   10 if(.not.time0.and.iter.eq.0) then
c
c  set values at previous time level
c
        do 12 i=1,nbcprv
   12   coefp(i)=coef(i)
c
	call strazt(aztstp,aztp,nb,ihe3bc,1)
	call strazt(axstp,axp,nb,ihe3bc,1)
      end if
c
c  set az1
c
      call storaz(az,az1,nb,ihe3bc,1)
      if(dgbcs2) then
	write(istdpr,'(a,1p10e13.5)') 
     *    'After storaz, az :',(az(i),i=1,nb+2)
	write(istdpr,'(a,1p10e13.5)') 
     *    'After storaz, az1:',(az1(i),i=1,nb+2)
      end if
c
      call zero(dztz,nbmax*naztmx)
c
      rhat=az1(nb1)
      rhsq=rhat**2
c
      if(time0) then
        jmax=2
c
c  for time0 azt(3 - nb) are assumed to be already set in
c  aztst in common/ksider/. hence we may set ax(3 - nb) directly.
c  note: this was modified 12/6/89. previously azt(kp) was set, from
c  assumed knowledge of ax(kp)
c
c  set azt from aztst in common/ksider/
c
        call strazt(aztst,azt,nb,ihe3bc,1)
	if(dgbcs2) write(istdpr,*) 'After time0 storage aztst =',
     *    aztst,'  azt =',azt
c
        do 26 kp=3,nb
c..      azt(kp)=az1(kp)-rhsq*ax(kp)
        ax(kp)=(az1(kp)-azt(kp))/rhsq
c
c  as a hack, set ax(kp) to 0 for time 0. Should be o.k., in general
c
        ax(kp)=0.d0
        if(dgbcs2) then
          write(istdpr,*) 'setting ax(kp). ax(kp),az1(kp),azt(kp),rhsq:'
          write(istdpr,*) ax(kp),az1(kp),azt(kp),rhsq
        end if
c  zero the corresponding dax-es
        do 26 i=1,nbmax
   26   dax(kp,i)=0
c
      else
c
        jmax=nb
c
      end if
c
      iconv=1
c
c  initial values of azt
c
      do 28 i=1,jmax
   28 azt(i)=az1(i)
c
c  store azt in aztst in common/ksider/
c
      call strazt(azt,aztst,nb,ihe3bc,2)
      if(dgbcs2) write(istdpr,*) 'After general storage aztst =',
     *  aztst,'  azt =',azt
c
   30 if(idgbcs.ge.1.and.istdpr.gt.0) write(istdpr,481)
c
      skipt=.false.
c
      do 90 nit=1,nmax
c
      if(dgbcs3) write(istdpr,*) 'Before call 1 of bccoef'
      call bccoef(azt,aztst,zhc,pt,dpt,tt,dtt,th,dth,xt,dxt,xxh,dxh,
     *  idthm,nbmax,nb,inenti)
c
      if(kdgrhb.lt.0) return
c
      if(dgbcs3) write(istdpr,*) '2. dxt(1,1):', dxt(1,1)
c
      skipt=.true.
c
c  calculate the ax's and dax's
c
c  ax(1)
c
   32 ax(1)=pt
      do 35 k=1,nb
   35 dax(1,k)=dpt(k)
c
c ax(2)
c
   40 ax2=tt
      if(.not.time0) then
c
	if(idiffus.gt.0) then
	  ilum=4+idcomp
	else
	  ilum=4
	end if
        phi4=(1-theta(ilum))/theta(ilum)
c
        ax2=ax2+phi4*(ttp-axp(2))
        do 43 j=1,3
        ax2=ax2+(th(j)+phi4*thp(j))*(azt(j)-aztp(j))/dt
   43   continue
      end if
c
      ax(2)=ax2
      dradc=ax2*azt(1)/(azt(2)*ax(1))
c  test for convection
      axc2=azt(2)*dad(1)*ax(1)/azt(1)
      if(dgbcs2) write(istdpr,*) 'axr2, axc2:',ax(2),axc2
      conv=ax2.lt.axc2
      if(conv) go to 48
c  radiative case
      do 47 k=1,nb
      daxk=dtt(k)
      if(.not.time0.and.k.le.3) then
        do 46 j=1,3
   46   daxk=daxk+dth(j,k)*(azt(j)-aztp(j))/dt
        daxk=daxk+(th(k)+phi4*thp(k))/dt
      end if
c
   47 dax(2,k)=daxk
c
      go to 50
c
c  convective case
c
   48 ax(2)=axc2
c
c  derivatives, neglecting derivatives of dad
c
      dax(2,1)=axc2*(dpt(1)/pt-1/azt(1))
      dax(2,2)=axc2*(1/azt(2)+dpt(2)/pt)
      dax(2,3)=axc2*dpt(3)/pt
      dax(2,4)=0
c
c  here we should ideally set composition derivatives to zero.
c  however, for testing prescribed X(q), only write warning, and
c  do not reset.
c
c  Modified 7/3/08, to reset when icncbc .ge. 1
c
      if(icncbc.lt.1) then
        if(istdpr.gt.0) write(istdpr,*) 
     *    '***** Warning. Composition derivatives not',
     *    ' reset in convective core'
c
      else
c
c  set composition derivatives to zero
c
        do 49 kp=3,nb
        ax(kp)=0
        if(time0) azt(kp)=az(kp)
        do 49 m=1,nb
   49   dax(kp,m)=0
        if(istdpr.gt.0) write(istdpr,*) 
     *    'Composition derivatives set to zero in convective core'
c
      end if
c
      go to 70
c
c    -------------------------------------
c
c  start setting composition quantities
c
c  when finding expansion for a trial solution , or for time0, ax(3)
c  and ax(4) are already set
c
   50 if(iter.eq.-1.or.time0) go to 70
c
c  test for resetting abundances to very small values, possibly
c  with zero derivatives, under various circumstances (still needs
c  some thought). This also depends on whether He burning is included.
c  Logics: derivatives from nbc1 to nbc2 are set to zero.
c  Derivatives from nbc2+1 to nb are calculated (if there are any).
c  nbc2 = nb is used to flag for no setting of derivatives.
c  Values and derivatives from 3 to nbh are involved in hydrogen burning.
c
      nbc1=3
      nbc2=2
c
c  Test for very small previous value of hydrogen abundance
c
      if(aztp(3).lt.1.e-8) then
        az1(3)=1.e-10
	nbc1=3
	nbc2=nbh
	if(istdpr.gt.0) write(istdpr,*) 
     *    'In bciter, reset central hydrogen abundance to 1.e-10'
      end if
c
c  Test for very small previous value of helium abundance
c
      if(iheccs.ne.0.and.aztp(nb-1).lt.1.e-8) then
        az1(nb-1)=1.e-10
	nbc1=3
	nbc2=nb
	if(istdpr.gt.0) write(istdpr,*) 
     *    'In bciter, reset central helium abundance to 1.e-10'
      end if
c
c  as a half-baked fudge, set derivatives to zero 
c  when hydrogen abundance at innermost point is too small. 
c#ai# This really requires more thought.
c
      if(az1(3).le.2.e-3) then
	nbc1=3
	nbc2=nbh
	if(istdpr.gt.0) write(istdpr,*) 
     *    'In bciter, set central hydrogen derivatives to zero'
      end if
c
c  as a half-baked fudge, set derivatives to zero 
c  when helium abundance at innermost point is too small. 
c#ai# This really requires more thought.
c
      if(iheccs.ne.0.and.az1(nb-1).le.2.e-3) then
	nbc1=3
	nbc2=nb
	if(istdpr.gt.0) write(istdpr,*) 
     *    'In bciter, set central helium derivatives to zero'
      end if
c
c  Zero derivatives when iterating with mixed core.
c
      if(nmxcor.gt.0) then
	nbc1=3
	nbc2=nb
      end if
c
c  zero selected derivatives
c
      if(nbc1.le.nbc2) then
        do 51 kp=nbc1,nbc2
        ax(kp)=0
        azt(kp)=az1(kp)
        do 51 m=1,nb
   51   dax(kp,m)=0
c
c  test if this is all
c
	if(nbc2.ge.nb) go to 70
      end if
c
c  set remaining values of ax and derivatives, from nbc2+1 to nb
c
      if(idgbcs.ge.1.and.istdpr.gt.0) 
     *  write(istdpr,'('' Setting ax(i) from'',i2,
     *    '' to'',i2)') nbc2+1,nb
c
c  calculate (delta t)*(-2/3 dln rhoc/dt) and its derivatives
c
      drc=0
      drcp=0
      call zero(ddrc,nbmax)
      do 52 i=1,3
      dzt=azt(i)-aztp(i)
      drc=drc+xt(i)*dzt
      drcp=drcp+xtp(i)*dzt
      if(dgbcs3) write(istdpr,*) 
     *  'i, azt(i), aztp(i), dzt, xt(i), drc',
     *  i, azt(i), aztp(i), dzt, xt(i), drc
      do 52 j=1,3
      ddrc(j)=ddrc(j)+dxt(i,j)*dzt
      if(dgbcs3) write(istdpr,*) 
     *  'j, dxt(i,j), dxt(i,j)*dzt, ddrc(j)',
     *  j, dxt(i,j), dxt(i,j)*dzt, ddrc(j)
   52 continue
      do 53 j=1,3
      ddrc(j)=ddrc(j)+xt(j)
      if(dgbcs3) write(istdpr,*) 'j, xt(j), ddrc(j)',
     *  j, xt(j), ddrc(j)
   53 continue
c
c  set up equations for ax and dax
      irpt=-1
   54 irpt=irpt+1
      nbm2=nb-nbc2
      do 61 kp=nbc2+1,nb
      thtp=theta(kp+2)
      thtu=1-thtp
      if(dgbcs3) write(istdpr,'(/a,i2,1p2e13.5)') 
     *  'kp, theta(kp), drc', kp, thtp, drc
      k1=kp-nbc2
      k2=kp-2
      do 55 l=1,nbm2
   55 wrk(k1,l)=xxh(l+nbc2,k2)*dt*thtp
      wrk(k1,k1)=1+wrk(k1,k1)+drc*thtp
      sum=0
      if(dgbcs3) write(istdpr,'(/a)') 'm, pres., prev., sum'
      do 57 m=1,nb
      sum=sum+thtu*xxhp(m,k2)*axp(m)
      if(m.le.2) sum=sum+thtp*xxh(m,k2)*ax(m)
      if(dgbcs2) then
	if(m.le.2) then
	  write(istdpr,'(i2,1p7e11.3)') 
     *      m, xxh(m,k2), ax(m), xxhp(m,k2), axp(m),
     *      thtp*xxh(m,k2)*ax(m), thtu*xxhp(m,k2)*axp(m), sum
        else
	  write(istdpr,'(i2,1p7e11.3)') 
     *      m, 0.d0, 0.d0, xxhp(m,k2), axp(m),
     *      0.d0, thtu*xxhp(m,k2)*axp(m), sum
	end if
      end if
   57 continue
      wrk(k1,nb-1)=axp(kp)*(1-drcp*thtu)-dt*sum
      if(dgbcs3) write(istdpr,'(a,1p3e13.5)')
     *  'axp(kp)*(1-drcp*thtu), dt*sum, wrk(k1,nb-1)',
     *  axp(kp)*(1-drcp*thtu), dt*sum, wrk(k1,nb-1)
      do 61 m=1,nb
      sum=0
      do 59 l=1,nb
      if(l.le.2) then
	sum=sum+xxh(l,k2)*dax(l,m)
	if(dgbcs3) write(istdpr,*) 
     *    'kp, m, l, xxh(l,k2)*dax(l,m), sum',
     *    kp, m, l, xxh(l,k2)*dax(l,m), sum
      end if
      sum=sum+dxh(l,m,k1)*ax(l)
      if(dgbcs3) write(istdpr,*) 'kp, m, dxh(l,m,k1)*ax(l), sum',
     *  kp, m, dxh(l,m,k1)*ax(l), sum
   59 continue
      sum=-dt*sum*thtp
      if(m.le.3) sum=sum-xtp(m)*axp(kp)*thtu
      wrk(k1,nb-1+m)=sum-ddrc(m)*ax(kp)*thtp
      if(dgbcs3) then
	write(istdpr,*) 'm, kp, ddrc(m),ax(kp),thtp',
     *    m, kp, ddrc(m),ax(kp),thtp
        write(istdpr,*) 'wrk,k1,nb-1+m,sum,ddrc(m)*ax(kp)*thtp',
     *    wrk(k1,nb-1+m),k1,nb-1+m,sum,ddrc(m)*ax(kp)*thtp
      end if
   61 continue
c  solve for ax and dax
      if(dgbcs2) then
	write(istdpr,'(/a,i3)') ' Call leq with nbm2 =',nbm2
	write(istdpr,*) ' Coefficients in first call of leq:'
	do j=1,nbm2
	  write(istdpr,'(i3,1p10e13.5)') j,(wrk(j,k),k=1,nbm2)
        end do
	write(istdpr,*) ' RHS in first call of leq:'
	do j=1,nbm2
	  write(istdpr,'(i3,1p10e13.5)') j,(wrk(j,nbm2+k),k=1,nb+1)
        end do
      end if
      call leq(wrk,wrk(1,nb-1),nbm2,nb+1,iwrk1,iwrk1,det)
      if(det.eq.0) then
        if(kdgrhb.gt.0) then
          kdgrhb = -1
          write(istder,*) 'Error in solution for ax'
          return
        else
          stop 'Error in solution for ax'
        end if
      end if
c
      if(dgbcs2) write(istdpr,'(/a,1p10e13.5)') ' ax set to',
     *  (wrk(i,nb-1),i=1,nbm2)
c
      do 68 kp=nbc2+1,nb
      k1=kp-nbc2
      k2=kp-2
c
c#ai# as a slightly doubtful fudge, reset ax to zero in 
c     case of no reactions (determined from xxh(2,.), i.e., the
c     temperature derivative)
c
      if(xxh(2,k2).eq.0) then
	if(idgbcs.ge.1.and.istdpr.gt.0) 
     *    write(istdpr,'(/'' no reactions for kp ='',i2,
     *      '' zero ax'')') kp
        ax(kp)=0
        do 64 m=1,nb
   64   dax(kp,m)=0
      else
        ax(kp)=wrk(k1,nb-1)
        do 66 m=1,nb
   66   dax(kp,m)=wrk(k1,nb-1+m)
      end if
c
   68 continue
c  repeat with corrected ax
      if(irpt.eq.0) go to 54
c
c  With hydrogen settling, and idfxbc = 1, zero gradient
c  at centre, at fairly rough fudge.
c
   70 if(idiffus.ge.1.and.idfxbc.eq.1) then
	ax(3)=0
	azt(3)=az1(3)
	if(istdpr.gt.0) write(istdpr,*) 
     *    'Zero hydrogen derivative for diffusion in bciter'
      end if
c  
c  set up and solve equations for the corrections to azt and
c  the derivatives dztz(i,j) = (d azt(i))/(d az1(j)).
c
      call zero(wrk,iwrk1*iwrk2)
      nrhs=nb+2
      if(iter.eq.-1) nrhs=1
      if(dgbcs2) then
	write(istdpr,'(/a)') ' Terms in rhs(i):'
	do j=1,jmax 
	  write(istdpr,'(i3,1p4e15.7)') j, az1(j),azt(j),ax(j),rhsq
	end do
      end if
      do 76 j=1,jmax
      do 72 k=1,jmax
   72 wrk(j,k)=rhsq*dax(j,k)
      wrk(j,j)=1.d0+wrk(j,j)
      wrk(j,nb1)=az1(j)-azt(j)-ax(j)*rhsq
      if(iter.eq.-1) go to 76
      wrk(j,j+nb1)=1.d0
      wrk(j,2*nb1)=-2.d0*rhat*ax(j)
   76 continue
c
      if(dgbcs2) then
	write(istdpr,*) ' Coefficients in second call of leq:'
	do j=1,jmax
	  write(istdpr,'(i3,1p10e13.5)') j,(wrk(j,k),k=1,jmax)
        end do
	write(istdpr,*) ' RHS in second call of leq:'
	do j=1,jmax
	  write(istdpr,'(i3,1p10e13.5)') j,(wrk(j,nb+k),k=1,nrhs)
        end do
      end if
c
      call leq(wrk,wrk(1,nb1),jmax,nrhs,iwrk1,iwrk1,det)
      if(det.eq.0) then
        if(kdgrhb.gt.0) then
          kdgrhb = -1
          write(istder,*) 'Error in iteration for azt'
          return
        else
          stop 'Error in iteration for azt'
        end if
      end if
c
      if(dgbcs2) write(istdpr,'(/a,1p15e13.5)')
     *  'Corrections to zeta:',(wrk(j,nb1),j=1,jmax)
c
c  correct azt, set dztz and determine relative change
c
      sum1=0.d0
      sum2=0.d0
   80 do 82 i=1,jmax
      dzt=wrk(i,nb1)
      sum1=sum1+ abs(dzt)
      dzt=dzt+azt(i)
      sum2=sum2+ abs(dzt)
      azt(i)=dzt
      if(iter.ne.-1) then
        do 81 ib=1,nb1
   81   dztz(i,ib)=wrk(i,ib+nb1)
      end if
   82 continue
c
      eam=sum1/sum2
c
c  store new azt and ax in aztst and axst in common/ksider/
c
      call strazt(azt,aztst,nb,ihe3bc,2)
      call strazt(ax,axst,nb,ihe3bc,2)
c
c  diagnostic output
c
      if(idgbcs.ge.1.and.istdpr.gt.0) then
        do 83 j=1,nb
   83   wrk(1,j)=1-(azt(j)+ax(j)*rhsq)/az1(j)
        write(istdpr,485) nit,eam,(azt(j),j=1,nb),(wrk(1,j),j=1,nb)
      end if
c
c  test for catastrophic failure of iteration
c
      if(azt(1).le.0.or.azt(2).le.0) then
	write(istdou,487) (azt(j),j=1,nb)
	if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,487) 
     *    (azt(j),j=1,nb)
        if(kdgrhb.gt.0) then
          kdgrhb = -1
          write(istder,*) 'Disaster in bciter iteration'
          return
        else
	  stop 'Disaster in bciter iteration'
        end if
      end if
c
c  iterate for fl
c
      pl=log10(1.d17*azt(1))
      tl=log10(1.d7*azt(2))
      y=1-azt(3)-zhc
      nosd=.true.
      notd=.true.
      nitfl=0
c
   84 dpl=pl-log10(p(1))
      if( abs(dpl).lt.acy) go to 86
      if(nitfl.gt.10) go to 85
c
      dfl=dpl/p(2)
      fl=fl+dfl
      nitfl=nitfl+1
      if(dgbcs2) write(istdpr,*) 'nitfl,dpl,fl,tl,azt(3),y,zhc',
     *  nitfl,dpl,fl,tl,azt(3),y,zhc
      call eqstf(fl,tl,azt(3),y,zhc,nosd,notd)
      if(kdgeos.lt.0) then
        kdgrhb=-1
        return
      end if
      go to 84
c
c  diagnostics for unconverged fl-iteration
c
   85 write(istdou,493) pl,tl,fl,azt(3),dpl
      if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,493) 
     *  pl,tl,fl,azt(3),dpl
c
c  test for stopping
c
      if(abs(dpl).gt.1.e-3) then
        if(kdgrhb.gt.0) then
          kdgrhb = -1
          write(istder,*) 'Failed iteration in bciter'
          return
        else
          stop 'Failed iteration in bciter'
        end if
      end if
c
   86 continue
c
      if(eam.lt.acy) go to 95
   90 continue
c  end iteration for azt and ax
c  ****************************
      iconv=0
c  iteration unconverged,write diagnostics
      do 92 j=1,nb
   92 wrk(1,j)=1-(azt(j)+ax(j)*rhsq)/az1(j)
      write(istdou,495) nmax,eam,rhsq,(az1(j),j=1,nb)
      write(istdou,496) (azt(j),j=1,nb)
      write(istdou,497) (ax(j),j=1,nb)
      write(istdou,498) (wrk(1,j),j=1,nb)
c
      if(istdpr.ne.istdou.and.istdpr.gt.0) then
        write(istdpr,495) nmax,eam,rhsq,(az1(j),j=1,nb)
        write(istdpr,496) (azt(j),j=1,nb)
        write(istdpr,497) (ax(j),j=1,nb)
        write(istdpr,498) (wrk(1,j),j=1,nb)
      end if
c
   95 continue
c
c  store composition variables in common/compvr/ if there have been
c  reactions
c
      if(.not.norct) then
        cvr(1,nn+1)=aztst(3)
        cvr(2,nn+1)=1 - aztst(3) - zhc
        cvr(3,nn+1)=aztst(4)
        call setcvr(nn+1)
      end if
c  calculate coefficients for m and l
  100 call expcff(azt,ax,dax,am0,dm0,am2,dm2,alt,dl0,alh,dlh,altt,
     .  dl2,alhh,dlhh,albb,dlbb,idthm,nbmax,nb)
c
      if(kdgrhb.lt.0) return
c
c  restore ieqhe3 to original value
c
      ieqhe3=ieqho
c
      al0=alt
      al2=altt
      al01=al0
      al21=al2
c
      if(iter.eq.-1) go to 205
c
      call store(dl0,dlt,nbmax)
      call store(dl2,dltt,nbmax)
c
      if(time0) go to 150
c
  111 al0=al0+phi4*(altp-al0p)
      al2=al2+phi4*(alttp-al2p)
c
c  test for skipping entropy term
c
      if(inenti.ne.1) then
c
c  set entropy contributions separately, to allow test for
c  unreasonable size
c
        alen0=0.d0
        alen2=0.d0
        do 120 k=1,3
        dztk=(azt(k)-aztp(k))/dt
        alen0=alen0+(phi4*alhp(k)+alh(k))*dztk
  120   alen2=alen2+(phi4*alhhp(k)+alhh(k))*dztk
     .    +(phi4*albbp(k)+albb(k))*(ax(k)-axp(k))/dt
c
c  test for size
c
        if(abs(alen0).gt.0.02*al0.and.(cqc.gt.0.or.cqcp.gt.0)
     *     .and.azt(3).ge.0.01) then
          if(istdpr.gt.0) write(istdpr,510) alen0,al0
	  inenti=1
	  go to 30
        else
	  al0=al0+alen0
	  al2=al2+alen2
        end if
      end if
c
c  reset al01 and al21 for storage in common/ksider/
c
      al01=al0
      al21=al2
c
c  reset derivatives if entropy terms are included
c  (basic values are set by call of expcff above)
c
      if(inenti.ne.1) then
        do 130 j=1,3
        sum=0.d0
        sum1=0.d0
        do 125 k=1,3
        dztk=azt(k)-aztp(k)
        sum=sum+dlh(k,j)*dztk
  125   sum1=sum1+dlhh(k,j)*dztk+dlbb(k,j)*(ax(k)-axp(k))
     .    +(phi4*albbp(k)+albb(k))*dax(k,j)
        dl0(j)=dl0(j)+(sum+phi4*alhp(j)+alh(j))/dt
  130   dl2(j)=dl2(j)+(sum1+alhh(j)+phi4*alhhp(j))/dt
c
      end if
c
c  reset dztz to (d azt(i))/(d eta(j)), where eta(1) = log p,
c  eta(2) = log t, eta(3) = x, ..., eta(nb+1) = log r.
c
  150 do 160 ia=1,nb1
      if(ia.gt.2.and.ia.lt.nb1) go to 160
      za=amm*az1(ia)
      do 155 j=1,nb
  155 dztz(j,ia)=za*dztz(j,ia)
  160 continue
c
c  set quantities for boundary conditions
c
      rhat3=rhat**3
      algr=3.d0*log10(rhat)
      if(dgbcs2) then
	do 170 k=1,2
  170   write(istdpr,16091) k, (dec(j,2,k),j=1,nb)
      end if
16091 format(/' dec(j,2,',i1,'):'/(1p5e15.7))
c
      if(dgbcs2) write(istdpr,*) 'In bciter, alshft =',alshft
      do 190 ial=1,2
      ghh= abs(ec(ial,1)+rhsq*ec(ial,2))
c
      if(dgbcs3) write(istdpr,*) '10. dxt(1,1):', dxt(1,1)
c  
c  test for zero energy generation (the value adopted may have to be
c  checked)
c  
      if(ghh.le.0) ghh=1.d-10  
c
      if(ial.eq.1.or.alshft.eq.0) then
        gh(ial)=algr+log10(ghh)
	alfact=1.d0
      else
	gh(ial)=log10(alshft+rhat3*ghh)
	alfact=rhat3*ghh/(alshft+rhat3*ghh)
      end if
c
      ghh=amm*ghh
      do 180 ia=1,nb1
      sum=0.d0
      do 175 j=1,nb
  175 sum=sum+(dec(j,ial,1)+rhsq*dec(j,ial,2))*dztz(j,ia)
  180 dgh(ial,ia)=alfact*sum/ghh
      dgh(ial,nb1)=dgh(ial,nb1)
     *             +alfact*(3.d0+2.d0*rhsq*amm*ec(ial,2)/ghh)
c
      if(dgbcs2.and.ial.eq.2) then
	write(istdpr,18091) (dgh(2,k),k=1,nb1)
18091   format(//' dgh(2,.) before transformation:'/
     *    (1p5e15.7))
      end if
c
c  transform to derivatives wrt (log f, log t, x)
c
      do 182 k=2,3
  182 dgh(ial,k)=dgh(ial,k)+dgh(ial,1)*p(k+1)
      dgh(ial,1)=dgh(ial,1)*p(2)
c
      if(dgbcs2.and.ial.eq.2) then
	write(istdpr,18091) (dgh(2,k),k=1,nb1)
      end if
  190 continue
      do 197 i=1,nb
      do 196 k=2,3
  196 dztz(i,k)=dztz(i,k)+dztz(i,1)*p(k+1)
  197 dztz(i,1)=dztz(i,1)*p(2)
c
c  diagnostic output?
c
      if(idgbcs.ge.1.and.istdpr.gt.0) then
        write(istdpr,550) (azt(i),i=1,nb)
        write(istdpr,551) (ax(i),i=1,nb)
        write(istdpr,552) ec,gh(1),gh(2)
      end if
c
c  set values at centre for prtsol
c
      rcnoc=rcno
      uwc=rnuw(1)
c
c  central he3
c
      if(conv) then
c
c  convective core
c
c  here we should ideally set composition derivatives to zero.
c  however, for testing prescribed X(q), only write warning, and
c  do not reset.
c
c  Modified 7/3/08, to reset when icncbc .ge. 1
c

        if(icncbc.lt.1) then
          if(istdpr.gt.0) write(istdpr,*) 
     *      '***** Warning. Composition derivatives not',
     *      ' reset in convective core'
c
        else
          do 198 k=3,nb
          azt(k)=az(k)
  198     ax(k)=0
          if(istdpr.gt.0) write(istdpr,*) 
     *      'Composition derivatives set to zero in convective core'
        end if
        if(ihe3bc.eq.1) then
          aztst(4)=xhe3eq
          axst(4)=0
	end if
c
      else if((time0.and.agesh.gt.0).or.ihe3bc.eq.1) then
c
c  radiative core. test if he3 is included in solution
c
c  central value
c
        aztst(4)=xhe3eq
c
c  second derivatives from forcing second-order expansion
c
        axst(4)=(az(4)-aztst(4))/rhsq
	if(idgbcs.ge.1.and.istdpr.gt.0) 
     *    write(istdpr,610) az(4), aztst(4), axst(4)
c
      end if
c
c  test for initial, static model
c
  205 if(time0.or.iter.eq.-1) then
c
c  set complete set of coefficients
c
         time0=.false.
         ieqho=ieqhe3
         if(ihe3bc.eq.1.and.ifdhe3.ne.1) ieqhe3=1
	 if(dgbcs3) write(istdpr,*) 'Before call 2 of bccoef'
         if(dgbcs3) write(istdpr,*) '10x. dxt(1,1):', dxt(1,1)
         call bccoef(azt,aztst,zhc,pt,dpt,tt,dtt,th,dth,xt,dxt,xxh,dxh,
     *     idthm,nbmax,nb,inenti)
c
         if(kdgrhb.lt.0) return
c
         if(dgbcs3) write(istdpr,*) '10a. dxt(1,1):', dxt(1,1)
         call expcff(azt,ax,dax,amb0,dm0,amb2,dm2,alt,dlt,alh,dlh,altt,
     .     dltt,alhh,dlhh,albb,dlbb,idthm,nbmax,nb)
c
         if(kdgrhb.lt.0) return
c
         if(dgbcs3) write(istdpr,*) '10b. dxt(1,1):', dxt(1,1)
         time0=.true.
         ieqhe3=ieqho
      end if
c
      if(dgbcs3) write(istdpr,*) '11. dxt(1,1):', dxt(1,1)
c
c  output coefficients in ksider if idgbcs .ge. 1
c
      if(istdpr.gt.0) then
        if(idgbcs.ge.1) write(istdpr,560) (coef(i),i=1,nbcprv)
        if(dgbcs2) call strbcs(2)
      end if
c
      return
  481 format(//' output from iteration in bciter: nit,eam,zeta(j)',
     .  ',relative deviation.'/)
  485 format(i5,1pe13.5,2x,20e13.5)
  487 format(//
     *   ' ***** Error in s/r bciter: negative pressure or temperature'/
     *   '       zeta:',1p15e13.5)
  493 format(//' ***** iteration for fl unconverged in s/r bciter'/
     *  7x,'pl, tl, fl, xh =',4f12.6,'  final dpl =',1pe13.5)
  495 format(///1x,10(1h*),' iteration in bciter failed to converge',
     .  ' in',i5,' iterations.'/11x,' final mean relative correction =',
     .  1pe13.5///14x,'z:',5e15.7)
  496 format(/11x,'zeta:',1p5e15.7)
  497 format(/12x,'xsi:',1p5e15.7)
  498 format(/' rel. departure:',1p5e15.7)
  510 format(//' ***** Warning. ',
     *  'Excessive entropy contribution in bciter'/
     *         '       Entropy, energy terms in al0:',1p2e13.5/
     *         '       Restart iteration without entropy term'/)
  550 format(//' dump from bciter.'/' azt:',1p5e13.5/(5x,e13.5))
  551 format('  ax:', 1p5e13.5/(5x,e13.5))
  552 format('  ec:',1p4e13.5/'  gh:',2e13.5/)
  560 format(//' coefficients in /ksider/:'/(1p10e13.5))
  610 format(//' setting X3 derivative. az(4), aztst(4), axst(4) =',
     *  1p3e13.5)
      end
      subroutine bccoef(azt,aztst,zhc,pt,dpt,tt,dtt,th,dth,xt,dxt,xxh,
     *  dxh,idth,idxh,nb,inenti)
c
c  calculates coefficients needed for the second-order expansion
c  at the centre
c
c  modified 20/8/1984 to rescale variables in order to avoid
c  underflow problems on the recku univac.
c  now azt(1-4) = (1.e-17*p, 1.e-7*t, x, x3).
c  expansion coefficients returned are those relevant for
c  these variables.
c
c  modified 29/12/1984 to include resetting of derivatives of
c  ft (in loop 17) when using he3 fudge.
c
c  modified 11/5/89 to ignore entropy term in energy equation for
c  inentr=1. this is implemented by setting sds, sdd and sdt to zero
c
c  Modified 31/7/91, for generalizing storage to allow full treatment
c  of CNO cycle.
c
c  Modified 14/7/00, to pass internal value of inenti in argument
c  list, to contrast from global value inentr set in common.
c
c  Note: aztst is now passed in argument list as well as azt.
c  This is used for calling engenr in case where He3 abundance
c  fudge (for ifdhe3 = 1) is used. Also ensures that He3
c  abundance is set correctly in aztst.
c
      implicit double precision (a-h,o-z)
      logical time0,nosd,notd,norche,dgbcs2, dgbcs3
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.nnz.d.incl'
c
      parameter(idermx=((nspcmx+3)*(nspcmx+4))/2,nbmax=nspcmx+2)
c
      dimension azt(*),aztst(*),dpt(*),dtt(*),th(*),dth(idth,*),
     *  xt(*),dxt(idth,3),xxh(idxh,3),dxh(idxh,idxh,3)
      dimension eps(idermx),ft(idermx,nspcmx),eps1(idermx),rho1(14)
      common/ebcder/ epsc,ftc(nspcmx),ftds(nbmax,nspcmx),
     *  ftdd(nbmax,nbmax,nspcmx),eds(nbmax),edd(nbmax,nbmax)
      common/rbcder/ rhc1,rhds(3),rhdd(3,3),sds(3),sdd(3,3),sdt(3,3,3)
      common/logf/ fl
      common/eqstd/ dummy(14),rho(20),ht(20),p(20),dm1(8),dlt(4),
     *  dummy1(4)
      common/eqsout/ eaeq(30),xiieq(30),dneeq(10),dpheq(20),
     *  heeq(20),peeq(10),
     *  hieq(20),pieq(10),hheq(20),pheq(10),hreq(20),preq(10)
      common/cmtime/ age, time0
      common/ln10/ amm
      common/rhcn/ crhcn(5),icrhcn(3),inentr
      common/bccn/ b1,b2,b3,b4,nnb,iveb,icncbc,ihe3bc
      common/opcder/ akc,akd(3)
      common/opcxdr/ akxa
      common/engcnt/ xhzlm1, xhzlm2, nxhzer, icnocs, nspec,
     *  iheccs, nspect
      common/engfdg/ epsfdg, qepsf1, qepsf2, ifdep1
      common/he3fdg/ agesh,ifdhe3,iche30, iche31
      common/he3int/ xhe3(4)
      common/he3eql/ xhe3eq,ieqhe3,iequsd
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor
      common/cxhcnc/ compc(nspcmx)
      common/rnrhed/ alhe(10,10), norche
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3
c
c  common giving fundamental constants
c
      common /fconst/ pibase, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev
c
c  common giving derived constants, and constants used for
c  equation of state
c
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev1,exhm
c
      common/cdgrhb/ kdgrhb
      common/cdgphs/ kdgeos, kdgopc, kdgeng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
c  eps1,rhc1 and /opcder/ are for the later call of s/r expcff
c
      dgbcs2=idgbcs.ge.2.and.istdpr.gt.0
      dgbcs3=idgbcs.ge.3.and.istdpr.gt.0
c
      if(dgbcs2) write(istdpr,*) 
     *  ' Entering bccoef with idth, idxh, nb =', idth, idxh, nb
c
      pc=1.d17*azt(1)
      tc=1.d7*azt(2)
      xc=azt(3)
      tcl=log10(tc)
      yc=1-xc-zhc
      nosd=.false.
      notd=.false.
      ift=idermx
      call eqstf(fl,tcl,xc,yc,zhc,nosd,notd)
      if(kdgeos.lt.0) then
        kdgrhb=-1
        return
      end if
      rhc=rho(1)
      rhcl=log10(rhc)
      call opact(rhcl,tcl,xc,zhc,ak,rkr,rkt,rkx)
      if(kdgopc.lt.0) then
        kdgrhb=-1
        return
      end if
      akc=1.d1**ak
c
c  derivatives of kappa wrt (1.e-17*p,1.e-7*t,x)
c
      al=rho(2)/p(2)
      akd(1)=1.d17*akc/pc*al*rkr
      akd(2)=1.d7*akc/tc*(rkt-dlt(1)*rkr)
      akd(3)=akc*amm*(akxa+rkr*(rho(4)-al*p(4)))
c
c  energy generation
c
      if(ifdhe3.eq.1) then
        nosd=.false.
        call he3abd(fl,tcl,xc,yc,zhc,agesh,xhe3,anu,nosd)
        aztst(4)=xhe3(1)
c
	if(dgbcs2) write(istdpr,'(a,1p5e13.5)') 
     *    'Calling engenr(1) with fl, tcl, aztst(3-5) =',
     *    fl, tcl, aztst(3), aztst(4), aztst(5)
        call engenr(fl,tcl,aztst(3),yc,zhc,eps,ft,ift,nosd)
        if(kdgeng.lt.0) then
          kdgrhb=-1
          return
        end if
c
c  test for switching off hydrogen burning
c
	if(ifdep1.eq.10) call zero(ft,idermx*nspcmx)
c
c  reset derivatives
        do 17 i=2,4
        eps(i)=eps(i)+eps(5)*xhe3(i)
   17   ft(i,1)=ft(i,1)+ft(5,1)*xhe3(i)
c
c  test for convective core, He3 not in equilibrium in solution
c  (this still required fixing up with care; at the moment
c   use ihe3bc = 1 to avoid including He3 in solution in core).
c
      else if(ihe3bc.eq.1.and.iche30.eq.1.and.nmxcor.gt.0) then
	ieqheo=ieqhe3
	ieqhe3=0
        aztst(4)=compc(2)
	if(istdpr.gt.0) write(istdpr,*) 
     *    'In bccoef aztst(4) set to',compc(2),
     *    ' in convective core'
   	if(dgbcs2) write(istdpr,'(a,1p5e13.5)') 
     *    'Calling engenr(2) with fl, tcl, aztst(3-5) =',
     *    fl, tcl, aztst(3), aztst(4), aztst(5)
        call engenr(fl,tcl,aztst(3),yc,zhc,eps,ft,ift,nosd)
        if(kdgeng.lt.0) then
          kdgrhb=-1
          return
        end if
c
c  test for switching off hydrogen burning
c
	if(ifdep1.eq.10) call zero(ft,idermx*nspcmx)
	ieqhe3=ieqheo
c
      else
        nosd=.false.
   	if(dgbcs2) write(istdpr,'(a,1p5e13.5)') 
     *    'Calling engenr(3) with fl, tcl, azt(3-5) =',
     *    fl, tcl, azt(3), azt(4), azt(5)
        call engenr(fl,tcl,azt(3),yc,zhc,eps,ft,ift,nosd)
        if(kdgeng.lt.0) then
          kdgrhb=-1
          return
        end if
c
c  test for switching off hydrogen burning
c
	if(ifdep1.eq.10) call zero(ft,idermx*nspcmx)
c
      end if
c
c  change to derivatives of nonlogarithmic quantities
c
c  set number of variables for derivatives, depending on
c  nspect set by engenr
c
      iveb = nspect + 2
c
      if(dgbcs2) write(istdpr,*) 'rho =',(rho(i),i=1,10)
      call chdr1(rho,rho1,3)
      call chdr1(eps,eps1,iveb)
c  change to derivatives wrt (p,t,x)
      notd=.true.
      if(dgbcs2) write(istdpr,*) 'rho1 =',(rho1(i),i=1,10)
      call trsder(tc,rho1,rhc1,rhds,rhdd,sdt,3,3,3,notd,1)
      call trsder(tc,eps1,epsc,eds,edd,sdt,nbmax,3,iveb,notd,2)
      if(time0) go to 50
      do 36 k=3,nb
      k1=k-2
   36 call trsder(tc,ft(1,k1),ftc(k1),ftds(1,k1),ftdd(1,1,k1),
     .  sdt,nbmax,3,iveb,notd,2)
      if(dgbcs2) then
	nbprt=min0(nb,3)
	write(istdpr,24091) 'rhds',(rhds(j),j=1,nbprt)
	write(istdpr,24093) 'rhdd:'
	do 24041 k=1,nbprt
24041   write(istdpr,24095) (rhdd(k,j),j=1,nbprt)
	do k1=1,nb-2
	  write(istdpr,24092) k1,ftc(k1),k1,(ftds(j,k1),j=1,nb)
	end do
	write(istdpr,24093) 'ftdd(.,.,1):'
	do 24042 k=1,nb
24042   write(istdpr,24095) (ftdd(k,j,1),j=1,nb)
      end if
      notd=.false.
c
c  set coefficients for the energy equation
c  for inenti=1 ignore entropy terms. this is implemented by setting
c  sds, sdd and sdt to zero
c
      if(inenti.ne.1) then
        call trsder(tc,ht,hht,sds,sdd,sdt,3,3,3,notd,2)
	if(dgbcs2) then
	  write(istdpr,24093) 'Preliminary sdd:'
	  do 24049 k=1,nbprt
24049     write(istdpr,24095) (sdd(k,j),j=1,nbprt)
	end if
        sds(1)=sds(1)-1.d17/rhc
        do 42 j=1,3
        sdd(1,j)=sdd(1,j)+1.d17*rhds(j)/rhc**2
        do 42 k=1,3
   42   sdt(1,j,k)=sdt(1,j,k)
     *            +1.d17*(rhdd(j,k)-2.d0*rhds(j)*rhds(k)/rhc)/rhc**2
c
      else
        call zero(sds,3)
        call zero(sdd,9)
        call zero(sdt,27)
      end if
c
      if(dgbcs2) then
	write(istdpr,24091) 'sds',(sds(j),j=1,nbprt)
	write(istdpr,24093) 'Final sdd:'
	do 24051 k=1,nbprt
24051   write(istdpr,24095) (sdd(k,j),j=1,nbprt)
c
24091 format(/a/(1p5e15.7))
24092 format(/'ftc(',i2,') =',1pe15.7,' ftds(.,',i2,'):'/(1p5e15.7))
24093 format(/a)
24095 format(1p5e15.7)
      end if
c
c  set pt and its derivatives
c
c  note: with consistent value of g = 6.6732e-8, numerical
c  constant should be 1.39763e-2.
c
c  note: the factor 1.d5 comes from measuring p in units of 1.d17 and
c  r in units of 1.d11
c
c++   50 cpt=1.3965d-2
   50 cpt=0.666666666667d5*pibase*cgrav
      pt=-cpt*rhc**2
      do 52 i=1,nb
      dppt=0
      if(i.gt.3) go to 52
      dppt=2.d0*pt*rhds(i)/rhc
   52 dpt(i)=dppt
c
c  set th,tt,xxh,xt and their derivatives
c
c  central conductivity times 1.d15
c
c++      ccndc=5.51234d17
      ccndc=1.d15/(8*clight*car)
      cndc=ccndc*akc*rhc**2/tc**3
      tt=-cndc*epsc
      do 53 i=1,nb
      dttt=eds(i)
      if(i.le.3) dttt=dttt+epsc*(akd(i)/akc+2*rhds(i)/rhc)
   53 dtt(i)=-cndc*dttt
      dtt(2)=dtt(2)-3.d7*tt/tc
c
      if(time0) return
c
      do 56 i=1,3
      th(i)=cndc*sds(i)
      xt(i)=-0.6666666667*rhds(i)/rhc
      do 54 k=1,3
c
c  #AI# Note that the following statement is pure desparation
c       to make the code work on sundog
c
      if(idgbcs.ge.33333.and.istdpr.gt.0) write(istdpr,*) rhc
      dth(i,k)=cndc*(sds(i)*(akd(k)/akc+2.d0*rhds(k)/rhc)+sdd(i,k))
      dxt(i,k)=0.6666666667d0*(rhds(i)*rhds(k)/rhc-rhdd(i,k))/rhc
   54 continue
      dth(i,2)=dth(i,2)-3.d7*th(i)/tc
   56 continue
      do 58 k=3,nb
      k1=k-2
      do 58 i=1,nb
      xxh(i,k1)=-ftds(i,k1)
      do 58 j=1,nb
   58 dxh(i,j,k1)=-ftdd(i,j,k1)
      return
      end
      subroutine expcff(azt,ax,dax,amt,dmt,amtt,dmtt,alt,dlt,alh,dlh,
     .  altt,dltt,alhh,dlhh,albb,dlbb,idth,idax,nb)
c
c  calculates coefficients in the equations for the expansion coeffi-
c  cients for m and l. the thermodynamic quantities needed are trans-
c
c  consistent with modifications caused by rescaling of variables,
c  20/8/1984.
c
      implicit double precision (a-h,o-z)
      logical time0, dgbcs2, dgbcs3
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.nnz.d.incl'
c
      parameter(idermx=((nspcmx+3)*(nspcmx+4))/2,nbmax=nspcmx+2)
c
      dimension azt(*),ax(*),dax(idax,*),dmt(*),dmtt(*),dlt(*),
     .  dlh(3,*),dltt(*),dlhh(idth,*),alh(*),alhh(*),albb(*),
     *  dlbb(idth,*)
      common/ebcder/ epsc,ftc(nspcmx),ftds(nbmax,nspcmx),
     *  ftdd(nbmax,nbmax,nspcmx),eds(nbmax),edd(nbmax,nbmax)
      common/rbcder/ rhc,rhds(3),rhdd(3,3),sds(3),sdd(3,3),sdt(3,3,3)
      common/opcder/ akc,akd(3)
      common/cmtime/ age, time0
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3
c
c  common giving fundamental constants
c
      common /fconst/ pibase, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
c++      data pi43,pi45 /4.18879021d0,2.51327412d0/
c
      dgbcs2=idgbcs.ge.2.and.istdpr.gt.0
      dgbcs3=idgbcs.ge.3.and.istdpr.gt.0
c
      pi43 = (4./3.)*pibase
      pi45 = (4./5.)*pibase
c
c  rho2,eps2,s2
c
      s2=0.d0
      eps2=0.d0
      rh2=0.d0
      if(dgbcs2) then
        write(istdpr,*) 'Setting eps2, rh2 with nb =',nb
        write(istdpr,'(a,1p10e13.5)') 'ax:  ',(ax(i),i=1,nb)
        write(istdpr,'(a,1p10e13.5)') 'eds: ',(eds(i),i=1,nb)
        write(istdpr,'(a,1p10e13.5)') 'rhds:',(rhds(i),i=1,3)
      end if
      do 10 i=1,nb
      axi=ax(i)
      eps2=eps2+axi*eds(i)
      if(i.gt.3) go to 10
      rh2=rh2+axi*rhds(i)
      if(time0) go to 10
      s2=s2+axi*sds(i)
   10 continue
      if(dgbcs2)
     *  write(istdpr,'(a,1p4e13.5)') 'rhc, epsc, rh2, eps2', 
     *     rhc, epsc, rh2, eps2
c
c  coefficients for m
c
   12 amt=pi43*rhc
      amtt=pi45*rh2
      do 20 j=1,3
      dmt(j)=pi43*rhds(j)
      sum=0.d0
      do 15 i=1,3
   15 sum=sum+rhdd(i,j)*ax(i)+rhds(i)*dax(i,j)
   20 dmtt(j)=pi45*sum
      if(nb.gt.3) then
        do 21 i=4,nb
        dmt(i)=0
   21   dmtt(i)=0
      end if
c
c  coefficients for l not multiplying time derivatives
c
      alt=pi43*rhc*epsc
      altt=pi45*(epsc*rh2+eps2*rhc)
      do 24 j=1,nb
      dd=rhc*eds(j)
      if(j.le.3) dd=dd+epsc*rhds(j)
   24 dlt(j)=dd
c
      if(dgbcs2) then
	write(istdpr,24091) (eds(j),j=1,nb)
	write(istdpr,24093)
	do 24051 k=1,nb
24051   write(istdpr,24095) (edd(k,j),j=1,nb)
c
24091 format(/' eds:'/(1p5e15.7))
24093 format(/' edd:')
24095 format(1p5e15.7)
      end if
c
      do 27 j=1,nb
      sum=0
      do 25 k=1,nb
      dd=rhc*edd(k,j)
      if(j.le.3.and.k.le.3) dd=dd+epsc*rhdd(k,j)
   25 sum=sum+dd*ax(k)+dlt(k)*dax(k,j)
      dd=sum+eds(j)*rh2
      if(j.le.3) dd=dd+eps2*rhds(j)
   27 dltt(j)=pi45*dd
      do 30 j=1,nb
   30 dlt(j)=pi43*dlt(j)
c
      if(time0) return
c
c  coefficients for l multiplying time derivatives
c
   40 do 50 k=1,3
      alh(k)=-pi43*rhc*sds(k)
      sk2=0.d0
      do 42 l=1,3
   42 sk2=sk2+sdd(k,l)*ax(l)
      alhh(k)=-pi45*(rh2*sds(k)+rhc*sk2-0.666666667*s2*rhds(k))
      albb(k)=-pi45*rhc*sds(k)
      do 50 j=1,3
      ssdd=rhds(j)*sds(k)+rhc*sdd(k,j)
      dlh(k,j)=-pi43*ssdd
      sum=0.d0
      do 44 l=1,3
   44 sum=sum+(rhdd(l,j)*sds(k)+rhc*sdt(k,l,j)-0.6666667*sdd(l,j)
     .  *rhds(k))*ax(l)+(rhds(l)*sds(k)+rhc*sdd(k,l)-0.6666667*
     .  sds(l)*rhds(k))*dax(l,j)
      dlhh(k,j)=-pi45*(sum+sdd(k,j)*rh2+rhds(j)*sk2-0.666666666667*
     .  s2*rhdd(k,j))
   50 dlbb(k,j)=-pi45*ssdd
      return
      end
      subroutine strazt(azt1,azt2,nb,ihe3bc,icase)
c
c  Given azt1, set azt2, depending on icase, and on ihe3bc
c
c  icase = 1: azt1 is aztst as stored in common/ksider/,
c             azt2 is format internal to s/r bciter
c  icase = 2: azt1 is aztst as in format internal to s/r bciter,
c             azt2 is as stored in common/ksider/.
c
c  Original version: 31/7/91
c
c  Modified 9/8/02, allowing for possible inclusion of 4He burning
c
      implicit double precision (a-h, o-z)
      dimension azt1(*), azt2(*)
c
      if(icase.eq.1) then
        do 15 i=1,nb
        if(ihe3bc.eq.0.or.i.le.3) then
	  i1=i
	else
	  i1=i+1
        end if
   15   azt2(i) = azt1(i1)
c
      else
c
	do 35 i=1,nb
        if(ihe3bc.eq.0.or.i.le.3) then
	  i1=i
	else
	  i1=i+1
        end if
   35   azt2(i1) = azt1(i)
c
      end if
      return
      end
      subroutine storaz(az1,az2,nb,ihe3bc,icase)
c
c  Given az1, set az2, depending on icase, and on ihe3bc
c
c  icase = 1: az1 is az as passed in argument list of s/r bciter
c             from s/r evolbcs, az2 is format internal to s/r bciter
c  icase = 2: az1 is az as in format internal to s/r bciter,
c             az2 is as passed in argument list.
c
c  Original version: 31/7/91
c
c  Modified 9/8/02, allowing for possible inclusion of 4He burning
c
      implicit double precision (a-h, o-z)
      dimension az1(*), az2(*)
c
      if(ihe3bc.eq.0) then 
	ilim=4
	ishft=1
      else
	ilim=3
	ishft=2
      end if
c
      if(icase.eq.1) then
	do 15 i=1,nb
	if(i.le.ilim) then
	  i1=i
        else
	  i1=i+ishft
        end if
   15   az2(i)=az1(i1)
	az2(nb+1)=az1(5)
c
      else
c
	do 20 i=1,nb
	if(i.le.ilim) then
	  i1=i
        else
	  i1=i+ishft
        end if
   20   az2(i1)=az1(i)
	az2(5)=az1(nb+1)
c
      end if
      return
      end
      subroutine chdr1(a,a1,iv)
c
c  changes from first and second derivatives of log(a) to
c  derivatives of a. The number of independent variables must
c  be given in iv.
c
      implicit double precision (a-h,o-z)
      dimension a(*),a1(*)
      common/ln10/ amm
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      ii=iv+1
      aa=amm*a(1)
      a1(1)=a(1)
      do 10 i=1,iv
      i1=i+1
      ai=a(i1)
      do 5 j=i,iv
      ii=ii+1
    5 a1(ii)=aa*(a(ii)+amm*ai*a(j+1))
   10 a1(i1)=aa*ai
      return
      end
      double precision function summa(a,b)
      implicit double precision (a-h,o-z)
      c=a+b
      if( abs(c).lt.1.d-12* abs(a)) c=0.d0
      summa=c
      return
      end
      subroutine trsder(t,a1,a2,da2,dda2,ddda2,id1,id2,iv,
     .  notd,inb)
c
c  transforms from derivatives wrt t(i) = (log f, log t, x, x3) to
c  derivatives wrt tau(i) = (1.e-17*p, 1.e-7*t, x, x3)
c
c  on input a1(1, .... ) contain the function value and first,
c  second and possibly third derivatives, on standard form.
c
c  t inputs the temperature.
c
c  returns function value in a2, first derivatives in da2(i),
c  second derivatives in dda2(i,j) and possibly (if notd .ne. true)
c  third derivatives in ddda2(i,j,k).
c
c  modified 20/8/1984 to use scaled variables.
c
      implicit double precision (a-h,o-z)
      logical notd, dgbcs2, dgbcs3
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.nnz.d.incl'
      parameter (idim = nspcmx+2)
c
      dimension a1(*),da2(*),dda2(id1,*),ddda2(id2,id2,*)
      dimension db(idim),ddb(idim,idim),dddb(idim,idim,idim),wrk(3,3),
     *  tnv(idim,idim),dp(4),ddp(4,4),dddp(4,4,4),st(idim),st1(idim)
      common/eqstd/ dum(54),pt(20),dum1(16)
      common/ln10/ amm
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3
c
      common/cdgrhb/ kdgrhb
      common/cdgphs/ kdgeos, kdgopc, kdgeng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      dgbcs2=idgbcs.ge.2.and.istdpr.gt.0
      dgbcs3=idgbcs.ge.3.and.istdpr.gt.0
c
c  test for initialization of derivative matrix
c
      if(inb.gt.1) go to 20
c
c  find d ti/d tauj
c
      dp(4)=0
      call zero(ddp,16)
      call zero(dddp,64)
      call filder(pt,dp,ddp,dddp,4,4,3,.false.)
      do 1 i=1,3
    1 dp(i)=1.d-17*amm*pt(1)*pt(i+1)
      do 3 i=1,3
      do 3 j=1,3
      ppij=ddp(i,j)
      do 2 k=1,3
    2 ddda2(i,j,k)=ppij*pt(k+1)
    3 ddp(i,j)=1.d-17*amm*pt(1)*(ppij+amm*pt(i+1)*pt(j+1))
      do 4 i=1,3
      do 4 j=1,3
      do 4 k=1,3
    4 dddp(i,j,k)=1.d-17*amm*pt(1)*(dddp(i,j,k)+amm*(ddda2(i,j,k)+
     .  ddda2(i,k,j)+ddda2(j,k,i)+amm*pt(i+1)*
     .  pt(j+1)*pt(k+1)))
c
      if(dgbcs3) then
	write(istdpr,*) 'iv =',iv
	write(istdpr,'(a,1p3e13.5)') 'dp:', (dp(i),i=1,3)
	write(istdpr,'(a/(1p3e13.5))') 'ddp:', ((ddp(i,j),i=1,3),j=1,3)
      end if
c
c  now dp, ddp and dddp contain first, second and third
c  derivatives of 1.e-17
c
    5 do 6 j=1,3
    6 wrk(1,j)=dp(j)
      do 10 j=1,3
      do 10 i=2,3
   10 wrk(i,j)=0
      wrk(2,2)=1.d-7*amm*t
      wrk(3,3)=1
      do 15 i=1,idim
      do 12 j=1,idim
   12 tnv(i,j)=0
   15 tnv(i,i)=1
      call leq(wrk,tnv,3,3,3,idim,det)
      if(det.eq.0) then
        if(kdgrhb.gt.0) then
          kdgrhb = -1
          write(istder,*) 'Stop in s/r trsder'
          return
        else
          stop 'Stop in s/r trsder'
        end if
      end if
c
      if(dgbcs3) then
	write(istdpr,'(a/(1p3e13.5))') 'tnv:', ((tnv(i,j),i=1,3),j=1,3)
      end if
c
c  now  dti/d tauj is in tnv(i,j)
c
c     -------------------------------------
c
c  start computing derivatives of a
c
c  file derivatives to db,....
c
   20 call filder(a1,db,ddb,dddb,idim,idim,iv,notd)
c
      a2=a1(1)
c  find first derivatives
      do 25 j=1,iv
      sum=0
      do 23 i=1,iv
   23 sum=sum+db(i)*tnv(i,j)
   25 da2(j)=sum
  140 format(1p4e15.5)
c
c  equation for second derivatives
c
      tscl=1.d-7*t
      tt1=amm*tscl
      tt2=amm*tt1
      tt3=amm*tt2
c
      do 30 i=1,iv
      do 30 j=1,iv
   30 ddb(i,j)=ddb(i,j)-da2(1)*ddp(i,j)
      ddb(2,2)=ddb(2,2)-tt2*da2(2)
c
c  solve for second derivatives
c
      do 40 k=1,iv
      do 40 l=k,iv
      sum=0
      do 35 i=1,iv
      do 35 j=1,iv
   35 sum=sum+ddb(i,j)*tnv(i,k)*tnv(j,l)
   40 dda2(k,l)=sum
c
      call symm(dda2,iv,2,id1)
      if(notd) return
c
c  equations for third derivatives
c
      do 45 i=1,iv
      do 45 j=1,iv
      t1=ddp(i,j)
      t2=0
      if(i.eq.2.and.j.eq.2) t2=tt2
      do 42 m=1,iv
   42 st1(m)=dda2(1,m)*t1+dda2(2,m)*t2
      do 43 k=1,iv
   43 st(k)=st1(1)*dp(k)
      st(2)=st(2)+st1(2)*tt1
      do 44 k=3,iv
   44 st(k)=st(k)+st1(k)
      do 45 k=1,iv
   45 ddda2(i,j,k)=st(k)
c
      do 50 i=1,iv
      do 50 j=1,iv
      do 50 k=1,iv
   50 dddb(i,j,k)=dddb(i,j,k)-ddda2(i,j,k)-ddda2(i,k,j)-ddda2(j,k,i)
     .  -da2(1)*dddp(i,j,k)
      dddb(2,2,2)=dddb(2,2,2)-tt3*da2(2)
c
c  solve for third derivatives
c
      do 60 l=1,iv
      do 60 m=l,iv
      do 60 n=m,iv
      sum=0
      do 55 i=1,iv
      do 55 j=1,iv
      do 55 k=1,iv
   55 sum=sum+tnv(i,l)*tnv(j,m)*tnv(k,n)*dddb(i,j,k)
   60 ddda2(l,m,n)=sum
  160 format(//(1p4e15.5))
c
      call symm(ddda2,iv,3,id2)
      return
      end
      subroutine filder(a,db,ddb,dddb,id1,id2,iv,notd)
      implicit double precision (a-h,o-z)
      logical notd
      dimension a(*),db(*),ddb(id1,*),dddb(id2,id2,*)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      ii=iv+1
      jj=iv*(iv+1)
      jj=jj/2+ii
      do 10 i=1,iv
      db(i)=a(i+1)
      do 10 j=i,iv
      ii=ii+1
      ddb(i,j)=a(ii)
      if(notd) go to 10
      do 8 k=j,iv
      jj=jj+1
    8 dddb(i,j,k)=a(jj)
   10 continue
      call symm(ddb,iv,2,id1)
      if(notd) return
      call symm(dddb,iv,3,id2)
      return
      end
      subroutine suberr
c  dumps commons from s/r eqstf
      implicit double precision (a-h,o-z)
      common/eqstd/ c1(90)
      common/eqsout/ c2(210)
      common/dmuder/ c3(10)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      write(istdpr,105)
      write(istdpr,110) c1
      write(istdpr,115) c2
      write(istdpr,120) c3
      return
  105 format(///' output from suberr:'/1x,20(1h*))
  110 format(/' common/eqstd/:'/1p4e13.5/(10e13.5))
  115 format(/' common/eqsout/:'/(1p10e13.5))
  120 format(/' common/dmuder/:'/1p10e13.5)
      end
      subroutine dmpbcs
c  dumps commons involved in bcs at centre
      implicit double precision (a-h,o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
c
c  set additional dimensions. Note that idthm corresponds to
c  number of independent thermodynamic variables
c
      parameter(naztmx = nspcmx+3, nbmax = nspcmx+2)
c
      dimension bc(nbcprv)
      common/ksider/ al01,al21,aztst(naztmx),axst(naztmx),pt,tt,
     *  th(idthm),xt(idthm),xxh(nbmax,nspcmx),alt,altt,
     *  alh(idthm),alhh(idthm),albb(idthm),
     *  dztz(nbmax,naztmx),dax(nbmax,nbmax),dpt(nbmax),dtt(nbmax),
     *  dth(idthm,idthm),dxt(idthm,idthm),
     .  dxh(nbmax,nbmax,nspcmx),dlt(nbmax),dltt(nbmax),
     *  dlh(idthm,idthm),dlhh(idthm,idthm),dlbb(idthm,idthm)
      common/excf/ ex(4),
     *  dm0(nbmax),dl0(nbmax),dm2(nbmax),dl2(nbmax)
      common/ebcder/ epsc,ftc(nspcmx),ftds(nbmax,nspcmx),
     *  ftdd(nbmax,nbmax,nspcmx),eds(nbmax),edd(nbmax,nbmax)
      common/rbcder/ rho,rhds(3),rhdd(3,3),sds(3),sdd(3,3),sdt(3,3,3)
      common/he3int/ xhe3(4)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      equivalence(bc(1),al01)
c
      write(istdpr,100) (bc(i),i=1,nbcprv)
      write(istdpr,110) ex
      write(istdpr,120) epsc,rho,xhe3
      return
  100 format(//' output from dmpbcs.'//' ksider:'/(1p10e13.5))
  110 format(//' excf:',1p4e13.5)
  120 format(//' eps =',1pe13.5,'  rho =',e13.5//' xhe3:',1p4e13.5)
      end
      subroutine symm(a,nd,nr,ndt)
c  symm completes the completely symmetric array a from those elements
c  a(i,j,k,l,...) for which i.le.j.le.k.le.l.le.... .
c   nd: range of i,j,k,.....
c   nr: rank of a (i. e. number of indices)
c   ndt: dimension of a (must be the same for all indices)
      implicit double precision (a-h,o-z)
      logical ord,ordnr2,ordit
      dimension a(*)
      dimension i(30),ord(30),iam(30)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
    5 if(nd.eq.1.or.nr.eq.1) return
      nr1=nr-1
      nr2=nr-2
      nd1=nd-1
      ndt1=ndt-1
      do 10 k=1,nr1
      ord(k)=.true.
   10 i(k)=1
      ndr=ndt**nr
      ndr1=ndr/ndt
      do 15 k=2,nd
   15 iam(k)=ndr
      iam(1)=ndr1
      nf=ndr1
      nfc=ndt
      ngc0=ndr1/ndt
      ng=1+ngc0
      ngc=ngc0
      it=nr1
      iit=2
      ordit=.true.
      inr2=1
      ordnr2=.true.
   25 im=iit-1
   26 do 30 l=1,im
      if(l.eq.1) go to 28
      l1=l-1
      iaml=iam(l1)/ndt
      nf=nf+iaml
      ng=ng+ndr1
      iam(l1)=iaml
   28 a(ng)=a(1+nf)
   30 continue
      im1=im-1
      ng=ng-im1*ndr1
      if(im.eq.nd) go to 38
      nfc=nfc/ndt
c  reset nf because of change of im to 1
      im2=im-2
      if(im2) 38,34,32
   32 do 33 k=1,im2
   33 nf=nf+k*(iam(k+1)-iam(k))
   34 nf=nf-im1*iam(im1)
      do 35 k=1,im1
   35 iam(k)=iam(k)*ndt
   38 if(iit.lt.nd) go to 40
      if(it.le.1) return
      it=it-1
      iit=i(it)
      nfc=nfc*ndt
      ngc=ngc/ndt
      go to 38
c  reset nf because of change of nd's to 1
   40 go to (44),nfc
      do 42 k=1,nd1
   42 iam(k)=nfc*iam(k)
   44 iaml=iam(iit)/ndt
      nf=nfc*nf+iaml-ndr*nd1*(nfc-1)/ndt1
      iam(iit)=iaml
      nfc=ndt
c  reset ng
      ng=ng-nd1*(ndr1-ngc*ndt)/ndt1+ngc
c  reset i and ord
      iit=iit+1
      i(it)=iit
      if(nr2) 25,25,45
   45 if(it.eq.nr1) go to 50
      ngc=ngc0
      it1=it+1
      do 46 k=it1,nr1
      i(k)=1
   46 ord(k)=.false.
      im=nd
      go to(48),it
      if(ord(it-1)) ord(it)=iit.ge.i(it-1)
   48 it=nr1
      iit=1
      inr2=i(nr2)
      ordnr2=ord(nr2)
      ordit=.false.
   50 if(ordnr2) ordit=iit.ge.inr2
      if(ordit) go to 25
      go to 26
      end
      subroutine strbcs(icase)
c  Outputs storage of variables in common/ksider/ (for icase = 1)
c  or named variables (for icase = 2), or named variables and
c  derivatives (for icase = 3; NB: a huge set)
c
c  Original version: 18/11/03
c
      implicit double precision (a-h,o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
c
c  set additional dimensions. Note that idthm corresponds to
c  number of independent thermodynamic variables
c
      parameter(naztmx = nspcmx+3, nbmax = nspcmx+2)
c
      dimension bc(nbcprv)
      common/ksider/ al01,al21,aztst(naztmx),axst(naztmx),pt,tt,
     *  th(idthm),xt(idthm),xxh(nbmax,nspcmx),alt,altt,
     *  alh(idthm),alhh(idthm),albb(idthm),
     *  dztz(nbmax,naztmx),dax(nbmax,nbmax),dpt(nbmax),dtt(nbmax),
     *  dth(idthm,idthm),dxt(idthm,idthm),
     .  dxh(nbmax,nbmax,nspcmx),dlt(nbmax),dltt(nbmax),
     *  dlh(idthm,idthm),dlhh(idthm,idthm),dlbb(idthm,idthm)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      equivalence(bc(1),al01)
c
      write(istdpr,'(//'' Storage in common/ksider/ bc(.)''/
     *                 '' *******************************''/)')
      i1=1
      i2=2
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)')
     *  i1,i2, 'al01,al21'
      i1=i2+1
      i2=i2+naztmx
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)') i1,i2, 'aztst'
      i1=i2+1
      i2=i2+naztmx
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)') i1,i2, 'axst'
      i1=i2+1
      i2=i1+1
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)') i1,i2, 'pt,tt'
      i1=i2+1
      i2=i2+idthm
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)') i1,i2, 'th'
      i1=i2+1
      i2=i2+idthm
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)') i1,i2, 'xt'
      i1=i2+1
      i2=i2+nbmax*nspcmx
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)') i1,i2, 'xxh'
      i1=i2+1
      i2=i1+1
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)') 
     *  i1,i2, 'alt,altt'
      i1=i2+1
      i2=i2+idthm
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)') i1,i2, 'alh'
      i1=i2+1
      i2=i2+idthm
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)') i1,i2, 'alhh'
      i1=i2+1
      i2=i2+idthm
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)') i1,i2, 'albb'
      i1=i2+1
      i2=i2+nbmax*naztmx
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)') i1,i2, 'dztz'
      i1=i2+1
      i2=i2+nbmax*nbmax
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)') i1,i2, 'dax'
      i1=i2+1
      i2=i2+nbmax
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)') i1,i2, 'dpt'
      i1=i2+1
      i2=i2+nbmax
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)') i1,i2, 'dtt'
      i1=i2+1
      i2=i2+idthm*idthm
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)') i1,i2, 'dth'
      i1=i2+1
      i2=i2+idthm*idthm
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)') i1,i2, 'dxt'
      i1=i2+1
      i2=i2+nbmax*nbmax*nspcmx
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)') i1,i2, 'dxh'
      i1=i2+1
      i2=i2+nbmax
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)') i1,i2, 'dlt'
      i1=i2+1
      i2=i2+nbmax
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)') i1,i2, 'dltt'
      i1=i2+1
      i2=i2+idthm*idthm
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)') i1,i2, 'dlh'
      i1=i2+1
      i2=i2+idthm*idthm
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)') i1,i2, 'dlhh'
      i1=i2+1
      i2=i2+idthm*idthm
      write(istdpr,'( '' bc('',i3,'' - '',i3,''): '',a)') i1,i2, 'dlbb'
c
      if(icase.ge.2) then
	write(istdpr,'(/'' Values of variables:''/)')
        write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'al01,al21',
     *    al01,al21 
        write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'aztst',
     *    (aztst(i),i=1,naztmx)
        write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'axst',
     *    (axst(i),i=1,naztmx)
        write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'pt,tt',
     *    pt,tt 
        write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'th',
     *    (th(i),i=1,idthm)
        write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'xt',
     *    (xt(i),i=1,idthm)
        write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'xxh',
     *    ((xxh(i,j),i=1,nbmax),j=1,nspcmx)
        write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'alt,altt',
     *    alt,altt 
        write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'alh',
     *    (alh(i),i=1,idthm)
        write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'alhh',
     *    (alhh(i),i=1,idthm)
        write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'albb',
     *    (albb(i),i=1,idthm)
c
	if(icase.ge.3) then
          write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'dztz',
     *      ((dztz(i,j),i=1,nbmax),j=1,naztmx)
          write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'dax',
     *      ((dax(i,j),i=1,nbmax),j=1,nbmax)
          write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'dpt',
     *      (dpt(i),i=1,nbmax)
          write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'dtt',
     *      (dtt(i),i=1,nbmax)
          write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'dth',
     *      ((dth(i,j),i=1,idthm),j=1,idthm)
          write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'dxt',
     *      ((dxt(i,j),i=1,idthm),j=1,idthm)
          write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'dxh',
     *      (((dxh(i,j,k),i=1,nbmax),j=1,nbmax),k=1,nspcmx)
          write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'dlt',
     *      (dlt(i),i=1,nbmax)
          write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'dltt',
     *      (dltt(i),i=1,nbmax)
          write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'dlh',
     *      ((dlh(i,j),i=1,idthm),j=1,idthm)
          write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'dlhh',
     *      ((dlhh(i,j),i=1,idthm),j=1,idthm)
          write(istdpr,'(1x,a,'':''/(1p5e13.5))') 'dlbb',
     *      ((dlbb(i,j),i=1,idthm),j=1,idthm)
        end if
      end if
c
      return
  100 format(//' output from dmpbcs.'//' ksider:'/(1p10e13.5))
  110 format(//' excf:',1p4e13.5)
  120 format(//' eps =',1pe13.5,'  rho =',e13.5//' xhe3:',1p4e13.5)
      end
