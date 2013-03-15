      subroutine engenr(fl, tl,x,y,z,eps,ft,ift,nosd)
c
c  calculates energy generation rate and rate of change of H and 
c  other elements considered. These may include He3 and some or
c  all of the elements in the CNO cycle. Also, the triple-alpha and
c  the 4He + 12C processes may be included.
c
c  Treatment of He3 is controlled by ieqhe3 (in common/he3eql/).
c  If ieqhe3=1 equilibrium abundance
c  of He3 is assumed everywhere. If ieqhe3 = 2 equilibrium
c  abundance is used in calculating eps, but not for ft.
c
c  Treatment of CNO cycle depends on the parameter icnocs
c  in common /engcnt/.
c  If icnocs = 0, use old (i.e., before July 1991) treatment,
c  assuming CN cycle to be in equilibrium and neglecting ON branch.
c  If icnocs = 1, treat CN cycle in equilibrium as in old
c  treatment, assuming all CN to be in N14, but include conversion
c  between N14 and O16.
c  If icnocs = 2, treat CN cycle in equilibrium, taking into account
c  the distribution amongst C12, C13, and C14 but include conversion
c  between N14 and O16.
c  If icnocs = 3, treat evolution of C12, C13 and C14, but neglect
c  conversion between N14 and O16.
c  If icnocs = 4, follow evolution of C12, C13, N14, and O16.
c
c  Note: in treatment of ON cycle, equilibrium between O16
c  and O17 is assumed and the cycle going through O18 is neglected.
c
c  Inclusion of helium burning is controlled by iheccs. If iheccs = 1,
c  the triple-alpha and the 4He + 12C reactions are included, through
c  call of s/r enghec.
c
c  Arguments: t is temperature (in k). Array x contains abundances
c             of 'active' elements by mass (see below).
c             z is abundance of heavy elements. 
c             If nosd is .false. second derivatives are
c             calculated.
c
c  Meaning of elements in array x of abundances depends on ieqhe3
c  and icnocs.
c  ieqhe3 = 1,      icnocs = 0: x(1) = X
c  ieqhe3 = 0 or 2, icnocs = 0: x(1) = X, x(2) = X(He3)
c  ieqhe3 = 1,      icnocs = 1: x(1) = X, x(2) = X(N14)
c  ieqhe3 = 0 or 2, icnocs = 1: x(1) = X, x(2) = X(He3), x(3) = X(N14)
c  ieqhe3 = 1,      icnocs = 2: x(1) = X, x(2) = X(N14)
c  ieqhe3 = 0 or 2, icnocs = 2: x(1) = X, x(2) = X(He3), x(3) = X(N14)
c  ieqhe3 = 1,      icnocs = 3: x(1) = X, x(2) = X(C13), x(3) = x(N14)
c  ieqhe3 = 0 or 2, icnocs = 3: x(1) = X, x(2) = X(He3), x(3) = X(C13),
c                               x(4) = X(N14)
c  ieqhe3 = 1,      icnocs = 4: x(1) = X, x(2) = X(C12), x(3) = X(C13), 
c                               x(4) = x(N14)
c  ieqhe3 = 0 or 2, icnocs = 4: x(1) = X, x(2) = X(He3), x(3) = X(C12),
c                               x(4) = X(C13), x(5) = X(N14)
c
c  Screening is controlled by iscren = iscnuc+ 10*iscbe7.
c  If iscren = 0, no electron screening is included.
c  Otherwise, iscnuc controls screening of nuclear reactions:
c
c  iscnuc = 1: use weak screening in original formulation, with
c     thte = thtec (as given in input)
c  iscnuc = 2: use weak screening, with consistent thte
c  iscnuc = 3: use intermediate screening, as done by Bahcall.
c
c  iscbe7 controls screening of electron capture in Be7:
c
c  iscbe7 = 0: compute screening as Bahcall and Moeller 
c     (ApJ, vol. 155, 511)
c  iscbe7 = 1: use approximate screening factor = 1.2, as
c     suggested by Bahcall
c  
c
c  Returns: eps(1) is energy generation rate per unit mass (cgs
c           units) and eps(2-   ) are first and possibly second
c           derivatives of log eps with respect to log f, log T and
c           x(k).
c           ft(1,k) is rate of change of x(k) (in 1/sec) and
c           ft(2...,k) are first and possibly second derivatives
c           of ft(1,k).
c
c
c  Also sets abundances into common/compos/ for those elements
c  that have been treated. The remaining abundances are set to -1.
c
c  Note that s/r eqstf must be called before call of engenr.
c
c  updated january 1984 to avoid over- and underflow on
c  recku univac. this should eventually be made more consistent,
c  in the parameters employed. also treatment of rates of change
c  of abundances should be fixed up, along the lines of the
c  calculation of eps from a q eff.
c
c  arguments changed 18/1/84. old call was
c
c      call engenr(t,x,eps,ft,ift,nosd)
c
c  heavy element abundance was transmitted in
c
c      common /heavy/ z
c
c  modified on 29/12/1984, to use expansion in equilibrium he3
c  abundance for small ca1 (cf. statement 37 ff.)
c
c  modified 15/9/86: save statement added in all routines
c     (needed in f77, apparently)
c
c  modified 23/10/87, to set reaction rates less than about exp(-60)
c (or 1.e-26) to zero. See s/r rnrate
c
c  modified 12/8/90, making sure that energy generation is set to zero
c  when there are no reactions, or no hydrogen
c
c  ********************************************************************
c
c  Major update started 17/7/91, to allow for more elements,
c  complete treatment of CNO cycle.
c
c  17/7/91 - 22/7/91: Generalize storage in s/r engenr.
c
c  22/7/91 - 25/7/91: Extend number of reactions handled in s/r rnrate,
c  to include reactions in CNO cycle.
c
c  26/7/91: Continue generalizing storage in s/r engenr. Using
c  compressed storage for derivatives (currently only implemented
c  for nspec = 1 or 2). Separating CNO quantities in special
c  subroutine engcno.
c
c  29/7/91: Set up energy generation and rates of change of composition
c  for the various CNO cases in s/r engcno.
c
c  19/2/93: Include various options for screening, selected by
c  iscren, and implemented in new separate subroutine elscrn
c
c  19/2/93: Change calling sequence, from
c     call engenr(t,x,z,eps,ft,ift,nosd)
c  to
c     call engenr(fl, tl,x,y,z,eps,ft,ift,nosd)
c
c  ********************************************************************
c
c  Modified 17/3/94, fixing setting of qe when there is no He3
c
c  Modified 12/7/96, including JOP additions to energy generation
c  (principally to start pre-main-sequence evolution)
c
c  ********************************************************************
c
c  Corrected 5/8/97, by including engenr.nnz.d.incl, to set nspdmx etc
c
c  Modified 1/8/02, including helium burning.
c
c  Modified 7/8/05, including NACRE reaction rates (M. Bazot)
c
      implicit double precision (a-h,o-z)
      character cxnuc*4
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
      parameter (idermx = ((nspcmx+3)*(nspcmx+4))/2, idalt = 10)
c
      logical nosd,norct,calrnr,secder,nohrct,norche
      dimension x(*),eps(*),ft(ift,*)
      dimension st1(idermx),st2(idermx),alt(idalt,krnrmx),
     *  ddr1(5),ddr2(5),drx(4),dftx(4),rhe32(idermx),rhe33(idermx),
     *  cxnuc(nspcmx) 
c
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
      common/nmbmsh/ nn,nc1
      common/noiter/ iter1,ntime
      common/consts/ av,ah,ahe,az
      common/engout/ a1(idermx),a2(idermx),q34(idermx),qe(idermx),
     *  rt1(idermx),rt2(idermx),epp(idermx),epscno(idermx),
     *  fth1(idermx),fth2(idermx),fthcno(idermx)
      common/eqstd/ xii(4),ane(10),rho(10)
      common/rnratd/ al(10,krnrmx),norct
      common/rnrhed/ alhe(10,10),norche
      common/rcncns/ a(knucst),b(knucst),s(knucst,isnuc),
     *  zz(knucst),zz12(2,knucst),zsmh,acno,awght(knucwg),q(knucq), 
     *  knucpm, kqvals, kawght
      common/he3eql/ xhe3eq,ieqhe3,iequsd
      common/compos/ xseth, yset, xset3, xset12, xset13, xset14, xset16
      common/frecno/ rcno, rhec
      common/engcnt/ xhzlm1, xhzlm2, nxhzer, icnocs, nspec, 
     *  iheccs, nspect
      common/engfdg/ epsfdg, qepsf1, qepsf2, ifdeps
      common/rnrcnt/ irnrat, krnrat, irn12, irn13, irn14, 
     *  irn15c, irn15o, irn16, irn17
      common/cenghc/ epshec(idermx), fthec(idermx,2)
      common/ln10/ amm,amm2
      common/cnofrc/ fcno, xtlcno
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      external blengr, blstio
c..      data idrprt /1/
c
      save
c
c  statement function defining storage of second derivatives
c
      jjsder(i, j) = (i*(2*nspect+5-i))/2 +j +1
c
      calrnr=.true.
      go to 10
c
      entry engen1(fl,tl,x,y,z,eps,ft,ift,nosd)
c
c  entry engen1 assumes that s/r rnrate has already been called.
c
      calrnr=.false.
c
   10 continue
c
      t=10.d0**tl
c
      if(idgeng.gt.0.and.istdpr.gt.0) then
        write(istdpr,*) 
     *    'Enter engenr with fl, tl, xx:',fl, tl, (x(i),i=1,4)
        write(istdpr,*) 'icnocs, iheccs =', icnocs, iheccs
      end if
c
c  for simplicity, set ahe3 from array awght
c
      ahe3 = awght(1)
c
c  set nspec, nspect and irnrat, depending on ieqhe3, icnocs and iheccs
c
      call engcse(ieqhe3, icnocs, iheccs, nspec, nspxx3, nspcno, 
     *  nsphec, nspect, irnrat, cxnuc)
c
      secder = .not.nosd
c
c  number of derivatives
c
      ider1=nspect+3
      if(idgeng.eq.-2.and.istdpr.gt.0) 
     *  write(istdpr,*) 'nspec, nspect, ider1',
     *  nspec, nspect2, ider1
      ider2=((nspect+3)*(nspect+4))/2
      if(secder) then
        ider=ider2
      else
        ider=ider1
      end if
c
c  Set index for derivatives going only over H and possibly He3
c  Also set storage index for rates of change of CNO abundances
c
      if(ieqhe3.eq.0.or.ieqhe3.eq.2) then
	iderh = 5
	nspech = 2
      else
	iderh = 4
	nspech = 1
      end if
c
c  storage indices for derivatives
c
      nspcm2 = nspect - 2
      nspcmh = nspect - iderh + 3
      jjxx =   jjsder(3, 3)
      jjxx3 =  jjsder(3, 4)
      jjx3x3 = jjsder(4, 4)
c
      y=1-x(1)-z
c
c  set hydrogen and helium abundances in common/compos/
c
      xseth = x(1)
      yset = y
      if(ieqhe3.ne.1) then
	xset3 = x(2)
      else
        xset3 = -1.
      end if
c
c  initialize CNO variables to -1
c  actual values are set in s/r engcno
c
      xset12 = -1.
      xset13 = -1.
      xset14 = -1.
      xset16 = -1.
c
c  test for setting artificial increment of energy generation rate
c
      if(ifdeps.ge.1) then
        tl=log10(t)
        call epsjop(epsj,tl)
      else
	epsj=0
      end if
c
c  if no hydrogen zero contribution from PP and CNO 
c  and go to calculation of possible 4He contribution
c  (Note: `zero' rate reduced from 1.d-5 to 1.d-10, 2/8/02)
c
      if(x(1).le.xhzlm2) then
        epp(1)=1.d-10
	call zero(epp(2),ider-1)
	call zero(epscno(1),ider)
	do 11 k=1,nspec
	do 11 i=1,ider
   11   ft(i,k)=0
	nohrct=.true.
	if(idgeng.gt.0.and.istdpr.gt.0) 
     *    write(istdpr,*) 'No hydrogen reactions'
	go to 54
      end if
c
c  calculation of reaction rates and their derivatives
c
      if(calrnr) call rnrate(fl,tl,x(1),y,z,nosd)
c
      if(norct) then
	norche=.true.
	eps(1)=0
	go to 190
      end if
c
      nohrct=.false.
c
c  set alt to be derivatives of al (rather than log(al))
c
      do 15 k=1,krnrat
      ii=4
      alk=max(al(1,k)*amm,1.d-30)
      alt(1,k)=al(1,k)
      do 15 i=2,4
      if(secder) then
        do 12 j=i,4
        ii=ii+1
   12   alt(ii,k)=alk*(amm*al(i,k)*al(j,k)+al(ii,k))
      end if
   15 alt(i,k)=alk*al(i,k)
c
c  find qeff
c
      xh=x(1)
      xhe3=0
      if(ieqhe3.ne.1) xhe3=x(2)
c
c  test for 'infinite' alpha 2
c
      if(al(1,4).eq.0.or.al(1,6).eq.0) then
c
c  set q34 tilde constant
c
        q34(1)=q(3)+q(4)
        do 17 i=2,ider
        a2(i)=0
   17   q34(i)=0
        ifxq34=1
c
      else
c
c  set alpha 2 and derivatives
c
        aa2=ah*al(1,6)/(rho(1)*xh*al(1,4))
        ifxq34=0
        a2(1)=aa2
        do 20 i=2,4
   20   a2(i)=aa2*amm*(al(i,6)-al(i,4)-rho(i))
        a2(5)=0
        a2(4)=a2(4)-aa2/xh
c
        ca2=1+aa2
        q34(1)=q(3)+(aa2*q(4)+q(5))/ca2
c
c  test for zero derivative
c
        if(ca2.ge.1.d17) then
          dq=0
          do 22 i=2,ider
   22     q34(i)=0
        else
c
          dq=(q(4)-q(5))/(ca2*ca2)
          do 25 i=2,iderh
   25     q34(i)=dq*a2(i)
        end if
c
      end if
c
c  test for he3 in equilibrium
c
      if(ieqhe3.eq.0) then
c
c  he3 not in equilibrium
c
        drt1=ah/(ahe3*xh)
        drt1=drt1*drt1*al(1,2)/al(1,1)
        rrt1=drt1*xhe3*xhe3
        rt1(1)=rrt1
        drt2=2*ah*ah*y*al(1,3)/(ahe3*ahe*xh*xh*al(1,1))
        rrt2=drt2*xhe3
        rt2(1)=rrt2
        do 30 i=2,4
        ddr1(i)=amm*(al(i,2)-al(i,1))
   30   ddr2(i)=amm*(al(i,3)-al(i,1))
        ddr1(4)=ddr1(4)-2/xh
        ddr2(4)=ddr2(4)-2/xh-1/y
        ddr1(5)=0
        ddr2(5)=0
        do 31 i=2,4
        rt1(i)=rrt1*ddr1(i)
   31   rt2(i)=rrt2*ddr2(i)
        rt1(5)=2*drt1*xhe3
        rt2(5)=drt2
        qee=q(1)+rrt1*q(2)+rrt2*q34(1)
        qe(1)=qee
        amqe=amm*qee
        do 32 i=2,iderh
   32   qe(i)=(q(2)*rt1(i)+q34(1)*rt2(i)+rrt2*q34(i))/amqe
        call zero(a1,ider)
c
      end if
c
c  set equilibrium he3 abundance
c
      if(al(1,2).eq.0.or.al(1,3).eq.0) then
c
c  zero abundance and derivatives
c
        xhe3eq=0
        ca1=1
        aaa=1
        call zero(a1,ider)
c
        if(ieqhe3.eq.1) then
          fa1=1
c
c  These statements corrected on 17/3/94
c
c..          qee=q(1)+q34(1)
c..          qe(1)=qee
c..          dqe=2*q34(1)-q(2)
c..          amqe=amm*qee
c..          do 36 i=2,iderh
c..   36     qe(i)=q34(i)/amqe
c
          qee=q(1)
          qe(1)=qee
          dqe=0
          amqe=amm*qee
          do 36 i=2,iderh
   36     qe(i)=0
        end if
c
      else
c
   37   ca1=xh*ahe/(y*ah*al(1,3))
        ca1=ca1*al(1,1)*al(1,2)*ca1
c
c  test for expansion
c
        if(ca1.le.1.e-4) then
          aa1=ca1*(1-0.5d0*ca1)
          aaa=1+aa1
        else
          aaa= sqrt(1+2*ca1)
          aa1=aaa-1
        end if
c
        xhe3eq=ahe3*y*al(1,3)*aa1/(ahe*al(1,2)*2)
	if(idgeng.ge.1.and.istdpr.gt.0) 
     *    write(istdpr,*) 'Setting xhe3eq.',
     *      ' ahe3, y, al(1,3), aa1, ahe, al(1,2), xhe3eq =',
     *        ahe3, y, al(1,3), aa1, ahe, al(1,2), xhe3eq
c
c  store in common/compos/
c
	xset3 = xhe3eq
c
c  test for using equilibrium he3 abundance everywhere
c
        if(ieqhe3.gt.0) then
          a1(1)=aa1/4
          ca1=ca1/(4*aaa)
          do 38 i=2,4
   38     a1(i)=ca1*amm*(al(i,1)+al(i,2)-2*al(i,3))
          a1(4)=a1(4)+2*ca1*(1/xh+1/y)
          a1(5)=0
        end if
c
c  test for using equilibrium he3 abundance for eps
c
        if(ieqhe3.eq.1) then
          fa1=1+2*a1(1)
          qee=q(1)+(a1(1)*q(2)+q34(1))/fa1
          qe(1)=qee
          dqe=2*q34(1)-q(2)
          amqe=amm*qee
          do 40 i=2,iderh
   40     qe(i)=(q34(i)-dqe*a1(i)/fa1)/fa1/amqe
c
c  Note: the following statement was previously
c
c..          call zero(rt1,30)
c
c  which makes little sense. The following replacement seems 
c  plausible, but must be checked against detailed notes
c
          call zero(rt1,ider)
c
        end if
c
      end if
c
c  final pp energy generation rate
c
   50 call zero(epp,ider)
c
      cpp=xh/ah
      rhh=rho(1)
      epp(1)=cpp*cpp*av*rhh*al(1,1)*qe(1)/2
      do 52 i=2,iderh
      depp=qe(i)
      if(i.lt.5) depp=depp+rho(i)+al(i,1)
   52 epp(i)=depp
      epp(4)=epp(4)+2/(amm*xh)
c
c  reset from MeV to ergs
c
      epp(1)=1.d6*ergev*epp(1)
c
c  set contribution to energy generation and rate of change of
c  abundances from CNO cycle
c
      call engcno(t,x,z,alt,epscno,fthcno,ft,idalt,ift,
     *  ider1, ider2, ider, iderh, nspech, nspcmh, secder)
c
c
c  Entry point in case of no hydrogen burning
c
c  set contribution to energy generation, and rate of change of 4He and
c  12C abundances
c
   54 call enghec(fl,t,x,z,epshec,fthec,idermx,secder)
c
c  total energy generation rate
c
      epss=epp(1)+epscno(1)+epshec(1)
c
      eps(1)=epss
c
      if(epss.gt.1.d-10) then
        rcno=epscno(1)/epss
        rhec=epshec(1)/epss
        do 56 i=2,ider1
   56   eps(i)=(epp(1)*epp(i)+epscno(1)*epscno(i)+epshec(1)*epshec(i))/
     *         epss
c
      else
        rcno=0.d0
        rhec=0.d0
	call zero(eps(2),ider-1)
      end if
c
c  possibly add contribution from artificial increase 
c  in energy generation rate and scale derivatives
c
      if(ifdeps.ge.1) then
        eps(1)=eps(1)+epsj
        if(epsj.ne.0) then
          derfct=epss/(epss+epsj)
          do 57 i=2,ider1
   57     eps(i)=derfct*eps(i)
        end if
      end if
c
c  If no hydrogen burning, skip setting rates of change and
c  go to treatment of 4He burning
c
      if(nohrct) go to 86
c
c  rates of change of H and He3
c
   60 call zero(fth1,ider)
c
      cc1=xh/ah
      ccy=y/ahe
      che3=xhe3/ahe3
      cc3=che3*ccy
      cc1=cc1*cc1
      cc2=che3*che3
c
      if(ieqhe3.ne.1) then
c
c  he3 not in equilibrium
c
c  set auxiliary arrays
c
        if(xhe3.gt.0) then
          xhei3=1.d-28/xhe3
        else
          xhei3=1.d37
        end if
c
        if(alt(1,2).le.xhei3) then
          do 70050 i=1,ider
70050     rhe32(i)=0
        else
c
          do 70080 i=1,ider
70080     rhe32(i)=che3*alt(i,2)
        end if
c
        if(alt(1,3).le.xhei3) then
          do 70150 i=1,ider
70150     rhe33(i)=0
        else
c
          do 70180 i=1,ider
70180     rhe33(i)=che3*alt(i,3)
c
        end if
c
c  now set reaction rates
c
        do 72 i=1,4
   72   fth1(i)=-1.5d0*cc1*alt(i,1)+che3*rhe32(i)-ccy*rhe33(i)
        fth1(4)=fth1(4)+rhe33(1)/ahe-3*alt(1,1)*xh/(ah*ah)
        fth1(5)=(2*rhe32(1)-alt(1,3)*y/ahe)/ahe3
        do 74 i=1,4
   74   fth2(i)=cc1*alt(i,1)/2-che3*rhe32(i)-ccy*rhe33(i)
        fth2(4)=fth2(4)+xh*alt(1,1)/(ah*ah)+rhe33(1)/ahe
        fth2(5)=-(alt(1,3)*y/ahe+2*rhe32(1))/ahe3
c
c  find ft(i,2)
c
        rha3=rhh*ahe3
        ft(1,2)=fth2(1)*rha3
        do 76 i=2,iderh
        dft2=fth2(i)
        if(i.lt.5) dft2=dft2+amm*fth2(1)*rho(i)
   76   ft(i,2)=dft2*rha3
c
c  test for diagnostic output
c
        if(idgeng.gt.0.and.istdpr.gt.0) then
          write(istdpr,76091) t,rhh,xh,z,xhe3,(alt(1,i),i=1,krnrat)
          write(istdpr,76092) cc1,cc2,cc3,ccy,che3,
     *      rhe32(1),rhe33(1),fth2(1),ft(1,2)
76091     format(/' diagnostic output from s/r engenr.'/
     *    ' t, rho, x, z, xhe3 =',1p5e15.7/ ' al:',20e15.7)
76092     format(' cc1, cc2, cc3, ccy, che3 =',6e15.7/
     *    ' rhe32, rhe33 =',2e15.7,'  fth2, ft(1,2) =',2e15.7)
c
          write(istdpr,'(a,1p4e15.7)')
     *     'In engenr: epp, epscno, epshec, eps',
     *      epp(1),epscno(1),epshec(1),eps(1)
        end if
c
      else
c
c  he3 in equilibrium
c
   80   frx1=1+1/fa1
        crx1=cc1*al(1,1)
        fth1(1)=-frx1*crx1
        do 81 i=2,4
   81   drx(i)=amm*al(i,1)
        drx(4)=drx(4)+2/xh
        do 82 i=2,4
        dftx(i)=2*a1(i)/(fa1*fa1)-frx1*drx(i)
   82   fth1(i)=crx1*dftx(i)
c
c  test for diagnostic output
c
        if(idgeng.gt.0.and.istdpr.gt.0) then
          write(istdpr,76091) t,rhh,xh,z,xhe3,(alt(1,i),i=1,krnrat)
          write(istdpr,82091) cc1,cc2,cc3,ccy,che3,fth1(1)
82091     format(' cc1, cc2, cc3, ccy, che3 =',6e15.7/
     *    '  fth1 =',1e15.7)
c
          write(istdpr,'(a,1p3e15.7)')
     *     'In engenr: epp, epscno, eps',epp(1),epscno(1),eps(1)
        end if
c
      end if
c
c  set ft(i,1), including term from CNO
c
      rha1=rhh*ah
      ft(1,1)=rha1*fth1(1)+fthcno(1)
      do 85 i=2,ider1
      dft1=fth1(i)
      if(i.lt.5) dft1=dft1+amm*fth1(1)*rho(i)
   85 ft(i,1)=rha1*dft1+fthcno(i)
c
c  rates of change of 4He and 12C. Note that second derivatives are
c  not set.
c#ai# As a rather desparate attempt, set 4He rate from H rate when
c  there are no 4He reactions. This should be cleaned up.
c
   86 if(iheccs.ne.0) then
c
        do 90 k=1,2
        do 87 i=1,ider1
	if(k.eq.1.and.norche) then
          ft(i,nspec+k)=-ft(i,1)
        else
          ft(i,nspec+k)=fthec(i,k)
	end if
   87   continue
	if(secder) then
	  do 88 i=ider1+1,ider2
   88     ft(i,nspec+k)=0
	end if
c
   90   continue
c
      end if
c
      if(nosd.or.nohrct) return
c
c  ***************************************************
c
c  second derivatives
c
c  test for fixed q34 tilde
c
      if(ifxq34.ne.1.and.dq.ne.0) then
c
c  set second derivatives of q34 tilde
c
        do 92 i=2,iderh
   92   st1(i)=a2(i)/aa2
        ii=4
        jj=ider1
        do 96 i=2,iderh
        do 94 j=i,iderh
        jj=jj+1
        if(i.lt.5.and.j.lt.5) then
          ii=ii+1
          a2(jj)=aa2*(st1(i)*st1(j)+amm*(al(ii,6)-al(ii,4)-rho(ii)))
          q34(jj)=dq*(a2(jj)-2*a2(i)*a2(j)/ca2)
        else
          a2(jj)=0
          q34(jj)=0
        end if
   94   continue
   96   jj = jj + nspcmh
c
        da2=aa2/(xh*xh)
        a2(jjxx)=a2(jjxx)+da2
        q34(jjxx)=q34(jjxx)+dq*da2
c
      end if
c
c  test for he3 in equilibrium
c
      if(ieqhe3.eq.0) then
c
c  He3 not in equilibrium
c
        ii=4
        jj=ider1
        do 104 i=2,iderh
        do 102 j=i,iderh
        jj=jj+1
        dr1=ddr1(i)*ddr1(j)
        dr2=ddr2(i)*ddr2(j)
        if(i.lt.5.and.j.lt.5) then
          ii=ii+1
          dr1=dr1+amm*(al(ii,2)-al(ii,1))
          dr2=dr2+amm*(al(ii,3)-al(ii,1))
        end if
        rt1(jj)=rrt1*dr1
        rt2(jj)=rrt2*dr2
        if(i.eq.5) then
          rt1(jj)=rt1(jj)+2*ddr1(j)*xhe3*drt1
          rt2(jj)=rt2(jj)+ddr2(j)*drt2
        end if
        if(j.eq.5) then
          rt1(jj)=rt1(jj)+2*ddr1(i)*xhe3*drt1
          rt2(jj)=rt2(jj)+ddr2(i)*drt2
        end if
  102   continue
  104   jj = jj + nspcmh
c
        rt1(jjxx)=rt1(jjxx)+2*rrt1/xh/xh
        rt1(jjx3x3)=rt1(jjx3x3)+2*drt1
        rt2(jjxx)=rt2(jjxx)+(2/xh/xh-1/y/y)*rrt2
c
c..	if(idrprt.eq.1) write(6,*) 'rt1',rt1
c..	if(idrprt.eq.1) write(6,*) 'rt2',rt2
c
        do 106 i=2,iderh
        st1(i)=qe(i)*amqe
        if(abs(st1(i)).lt.1.e-15) st1(i)=0
  106   continue
c
        jj=ider1
        do 108 i=2,iderh
        do 107 j=i,iderh
        jj=jj+1
  107   qe(jj)=(-st1(i)*st1(j)/qee+q(2)*rt1(jj)+q34(1)*rt2(jj)
     .    +q34(i)*rt2(j)+q34(j)*rt2(i)+rrt2*q34(jj))/amqe
  108   jj = jj + nspcmh
c..	if(idrprt.eq.1) write(6,*) 'qe',qe
c
      else
c
c  He3 in equilibrium
c
        do 112 i=2,4
  112   st1(i)=a1(i)/ca1
        ii=4
        jj=ider1
        do 115 i=2,iderh
        do 114 j=i,iderh
        jj=jj+1
        if(i.lt.5.and.j.lt.5) then
          ii=ii+1
          a1(jj)=ca1*(st1(i)*st1(j)*(1-4*ca1/aaa)+
     .      amm*(al(ii,1)+al(ii,2)-2*al(ii,3)))
        else
          a1(jj)=0
        end if
  114   continue
  115   jj = jj + nspcmh
c
        a1(jjxx)=a1(jjxx)-2*ca1*(1/xh/xh-1/y/y)
c
        if(ieqhe3.eq.1) then
          do 117 i=2,iderh
  117     st1(i)=amqe*qe(i)
          jj=ider1
          do 119 i=2,iderh
          do 118 j=i,iderh
          jj=jj+1
  118     qe(jj)=(-st1(i)*st1(j)/qee+(q34(jj)-
     .       (2*(q34(j)*a1(i)+q34(i)*a1(j))
     .      +dqe*(a1(jj)-4*a1(i)*a1(j)/fa1))/fa1)/fa1)/amqe
  119     jj = jj + nspcmh
        end if
c
      end if
c
c  now set derivatives of energy generation rate
c
      ii=4
      jj=ider1
      do 125 i=2,iderh
      do 124 j=i,iderh
      jj=jj+1
      depp=qe(jj)
      if(i.lt.5.and.j.lt.5) then
        ii=ii+1
        depp=depp+rho(ii)+al(ii,1)
      end if
  124 epp(jj)=depp
  125 jj = jj + nspcmh
c
      epp(jjxx)=epp(jjxx)-2/(amm*xh*xh)
c
      do 127 i=2,ider1
  127 st1(i)=eps(i)*epss
c
c..      if(idrprt.eq.1) write(6,*) 'st1',st1
c..      if(idrprt.eq.1) write(6,*) 'epp',epp
c..      if(idrprt.eq.1) write(6,*) 'epscno',epscno
c
      jj=ider1
      do 128 i=2,ider1
      do 128 j=i,ider1
      jj=jj+1
  128 eps(jj)=(-amm*st1(i)*st1(j)/epss+epp(1)*(amm*epp(i)*epp(j)
     .  +epp(jj))+epscno(1)*(amm*epscno(i)*epscno(j)+epscno(jj)))/epss
      if(idgeng.ge.1.and.istdpr.gt.0) then
        write(istdpr,*) 'epp:',(epp(i),i=1,ider)
        write(istdpr,*) 'epscno:',(epscno(i),i=1,ider)
        write(istdpr,*) 'eps:',(eps(i),i=1,ider)
      end if
c
c..      idrprt = 0
c
c  rates of change of H and He3
c
      call zero(st1,ider)
      call zero(st2,ider)
c
      if(ieqhe3.ne.1) then
c
c  He3 not in equilibrium
c
        do 132 i=2,4
        st1(i)=xh*alt(i,1)/(ah*ah)+rhe33(i)/ahe3
  132   st2(i)=-(2*rhe32(i)+y*alt(i,3)/ahe)/ahe3
        st1(5)=0
        st2(5)=0
        ii=4
        jj=ider1
        do 135 i=2,iderh
        do 134 j=i,iderh
        jj=jj+1
        dft=0
        dft2=0
        if(i.lt.5.and.j.lt.5) then
          ii=ii+1
          dft=cc1*alt(ii,1)/2-che3*rhe32(ii)-ccy*rhe33(ii)
          dft2=fth2(1)*(rho(ii)+amm*rho(i)*rho(j))
        end if
        if(i.eq.4) dft=dft+st1(j)
        if(i.eq.5) dft=dft+st2(j)
        if(j.eq.4) dft=dft+st1(i)
        if(j.eq.5) dft=dft+st2(i)
        if(i.lt.5) dft2=dft2+rho(i)*fth2(j)
        if(j.lt.5) dft2=dft2+rho(j)*fth2(i)
        fth2(jj)=dft
  134   ft(jj,2)=rha3*(amm*dft2+fth2(jj))
  135   jj = jj + nspcmh
c
        dft=alt(1,1)/(ah*ah)
        fth2(jjxx)=fth2(jjxx)+dft
        ft(jjxx,2)=ft(jjxx,2)+rha3*dft
        dft=alt(1,3)/(ahe*ahe3)
        fth2(jjxx3)=fth2(jjxx3)+dft
        ft(jjxx3,2)=ft(jjxx3,2)+rha3*dft
        dft=-2*alt(1,2)/(ahe3*ahe3)
        fth2(jjx3x3)=fth2(jjx3x3)+dft
        ft(jjx3x3,2)=ft(jjx3x3,2)+rha3*dft
c
        do 137 i=2,4
        st1(i)=-3*xh*alt(i,1)/(ah*ah)+rhe33(i)/ahe
  137   st2(i)=(2*rhe32(i)-y*alt(i,3)/ahe)/ahe3
c
        ii=4
        jj=ider1
        do 139 i=2,iderh
        do 138 j=i,iderh
        jj=jj+1
        dft=0
        if(i.lt.5.and.j.lt.5) then
          ii=ii+1
          dft=-1.5d0*cc1*alt(ii,1)+che3*rhe32(ii)-ccy*rhe33(ii)
        end if
  138   fth1(jj)=dft
  139   jj = jj + nspcmh
c
        fth1(jjxx)=fth1(jjxx)-3*alt(1,1)/(ah*ah)
        fth1(jjxx3)=fth1(jjxx3)+alt(1,3)/(ahe*ahe3)
        fth1(jjx3x3)=fth1(jjx3x3)+2*alt(1,2)/(ahe3*ahe3)
c
      else
c
c  He3 in equilibrium
c
        ii=4
        jj=ider1
        do 144 i=2,iderh
        do 143 j=i,iderh
        jj=jj+1
        dft=0
        if(i.eq.5.or.j.eq.5) go to 144
        if(i.lt.5.and.j.lt.5) then
          ii=ii+1
          dft=crx1*(dftx(i)*drx(j)+2*(a1(jj)+a1(j)*drx(i)-4*a1(i)*a1(j)
     .      /fa1)/(fa1*fa1)-frx1*amm*al(ii,1))
        end if
  143   fth1(jj)=dft
  144   jj = jj + nspcmh
c
        fth1(jjxx)=fth1(jjxx)+2*crx1*frx1/(xh*xh)
        do 146 i=2,ider1
        st1(i)=0
  146   st2(i)=0
c
      end if
c
      ii=4
      jj=ider1
      do 155 i=2,ider1
      do 155 j=i,ider1
      jj=jj+1
      dft=fth1(jj)
      dft2=0
      if(i.lt.5.and.j.lt.5) then
        ii=ii+1
        dft2=fth1(1)*(rho(ii)+amm*rho(i)*rho(j))
      end if
      if(i.eq.4) dft=dft+st1(j)
      if(j.eq.4) dft=dft+st1(i)
      if(i.eq.5) dft=dft+st2(j)
      if(j.eq.5) dft=dft+st2(i)
      if(i.lt.5) dft2=dft2+rho(i)*fth1(j)
      if(j.lt.5) dft2=dft2+rho(j)*fth1(i)
      fth1(jj)=dft
  155 ft(jj,1)=rha1*(amm*dft2+fth1(jj))+fthcno(jj)
      return
c
c  no reactions. Only include possible artificial increment
c
  190 eps(1)=epsj
c
      do 192 i=2,ider
  192 eps(i)=0
      do 195 i=1,ider
      do 195 j=1,nspect
  195 ft(i,j)=0
      xhe3eq=0
c
c  set flag for no reactions (in case this is caused by no hydrogen)
c
      norct = .true.
c
      return
      end
      subroutine rnrate(fl,tl,x,y,z,nosd)
c
c  Sets reaction rates and their derivatives wrt (log f, log t, x)
c  into al. rate for reaction no i and its derivatives are set into
c  al(k,i), k = 1, ..., 10, where k = 1 gives the reaction rate and
c  k = 2, ..., 10 give the derivatives of log10(al(1,i)), 
c  stored in the usual way.
c
c  NB. in common with Fowler, Caughlan & Zimmermann (1975) the
c      rate for nuclear reations is defined such that al(1,i) is 
c      set to av*lambda(i,j) where av is Avogadro's number and
c      lambda(i,j) is the number of reactions per pair of particle
c      (i,j) per second. In the present implementation this holds
c      for i .noteq. 6. al(1,6) is ne*lambda(Be7,electron), where
c      ne is the electron density per unit volume, and again
c      lambda(Be7, electron) is the reaction rate per pair of particles.
c
c  Data for calculation is set in set in s/r srncns.
c
c  Note that equation of state
c  s/r must have been called before the call of rnrate.
c
c  The reactions that are considered depend on the value of 
c  the parameter irnrat passed in common /rnrcnt/.
c  In all cases the first 6 reactions are stored in
c
c   al(1,1): H1(H1,e+ nu)H2
c   al(1,2): He3(He3,2H1)He4
c   al(1,3): He3(He4,gamma)Be7
c   al(1,4): Be7(H1,gamma)Be8
c   al(1,5): N14(H1,gamma)O15
c   al(1,6): Be7(e-,nu)Li7
c
c  The remaining reactions depend in irnrat. 
c
c  For irnrat = 1, only those 6 reactions are computed.
c
c  For irnrat = 2:
c
c   al(1,7): N15(H1,He4)C12
c   al(1,8): N15(H1,gamma)O16
c   al(1,9): O16(H1,gamma)F17
c  al(1,10): O17(H1,He4)N14
c
c  For irnrat = 3:
c
c   al(1,7): C12(H1,gamma)N13
c   al(1,8): C13(H1,gamma)N14
c
c  For irnrat = 4:
c
c   al(1,7): C12(H1,gamma)N13
c   al(1,8): C13(H1,gamma)N14
c   al(1,9): N15(H1,He4)C12
c  al(1,10): N15(H1,gamma)O16
c  al(1,11): O16(H1,gamma)F17
c
c  For irnrat = 5:
c
c   al(1,7): C12(H1,gamma)N13
c   al(1,8): C13(H1,gamma)N14
c   al(1,9): N15(H1,He4)C12
c  al(1,10): N15(H1,gamma)O16
c  al(1,11): O16(H1,gamma)F17
c  al(1,12): O17(H1,He4)N14
c
c  The total number of reactions set is returned in krnrat.
c  Also, the routine sets indices irn12, irn13, irn14, 
c  irn15c, irn15o, irn16, irn17 such that in all cases 
c  the reactions that have been set are stored as
c
c   al(1,irn12):  C12(H1,gamma)N13
c   al(1,irn13):  C13(H1,gamma)N14
c   al(1,irn14):  N14(H1,gamma)O15
c   al(1,irn15c): N15(H1,He4)C12
c   al(1,irn15o): N15(H1,gamma)O16
c   al(1,irn16):  O16(H1,gamma)F17
c   al(1,irn17):  O17(H1,He4)N14
c
c  Indices for reactions that have not been set are returned as 0.
c
c  ******************************************************************
c
c  updated january 1984 to correct for under- and overflow on
c  the recku univac.
c
c  Modified 25/8/87, so that now uw(1) in common/rnrout/ returns
c  Uw. Previously it returned Uw/ln(10).
c
c  Modified 23/10/87, to use limits +- 60 in fxp, and in test for
c  zero reaction rate
c
c  Modified 23/7/91, allowing calculation of reactions in CNO cycle
c
c  19/2/93: Change calling sequence, from
c     call rnrate(t,x,y,z,nosd)
c  to
c     call rnrate(fl,tl,x,y,z,nosd)
c
c  Modified 29/5/03, replacing cut-off of 60 in exponential by 
c  parameter coexp, set in parameter statement (currently to -120)
c
c  Modified 07/08/05 call to cnrnacre for NACRE reaction rates
c  inclusion of common/cversn/ to distinguish the case ivreng=8
c
c  Note: there is a discrepancy for the second derivatives
c        involving log f of the electron screening for be7 + e
c        between the low-temperature approximation and the full
c        expression. this must be checked sometime. 18/1/84.
c
      implicit double precision (a-h,o-z)
      logical nosd,norct,norct1,secder
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
      parameter (krnrm1 = krnrmx + 1, iptdat = krnrm1 - 12)
c
c  cut-off parameter for exponential
c
      parameter(coexp = 120.d0)
c
      dimension iptrnr(krnrm1, 4)
      dimension st1(10),st2(10),tpw(5),xhi(2),anel(10),
     *  thte(10), zzscr(krnrmx)
      common/cversn/ ivrbcs,ivrrhs,ivreos,ivropc,ivreng
      common/rnrout/ zt(10),uw(10),rs(10),ee(20),alr(10),
     *  dd(10),sr(10),f(10),alc(10)
      common/rnratd/ al(10,krnrmx),norct
      common/rnrcnt/ irnrat, krnrat, irn12, irn13, irn14, 
     *  irn15c, irn15o, irn16, irn17
      common/rcncns/ a(knucst),b(knucst),s(knucst,isnuc),
     *  zz(knucst),zz12(2,knucst),zsmh,acno,awght(knucwg),q(knucq), 
     *  knucpm, kqvals, kawght
      common/eqstd/ xii(4),ane(10),rho(10)
      common/degfct/ thtec,iscren
      common/consts/ av,ah,ahe,az
      common/ln10/ amm,amm2
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      external blengr
c
c  set pointing indices to data in common /rcncns/
c
c                  1  2  3  4  5  6  7  8   9  10  11  12
      data iptrnr /1, 2, 3, 4, 5, 0, 0, 0,  0,  0,  0,  0, iptdat*0,
     *             1, 2, 3, 4, 5, 0, 8, 9, 10, 11, 10, 11, iptdat*0,
     *             1, 2, 3, 4, 5, 0, 6, 7,  0,  0,  0,  0, iptdat*0,
     *             1, 2, 3, 4, 5, 0, 6, 7,  8,  9, 10, 11, iptdat*0/
c
      save
c
      fxp(aa)= exp(min(coexp,max(aa,-coexp)))
c
      t=10.d0**tl
c
c  set storage indices, depending on irnrat
c
      irn14 = 5
c
      if(irnrat.eq.1) then
	krnrat = 6
	irn12 = 0
	irn13 = 0
	irn15c = 0
	irn15o = 0
	irn16 = 0
	irn17 = 0
      else if (irnrat.eq.2) then
	krnrat = 10
	irn12 = 0
	irn13 = 0
	irn15c = 7
	irn15o = 8
	irn16 = 9
	irn17 = 10
      else if (irnrat.eq.3) then
	krnrat = 8
	irn12 = 7
	irn13 = 8
	irn15c = 0
	irn15o = 0
	irn16 = 0
	irn17 = 0
      else if (irnrat.eq.4) then
	krnrat = 11
	irn12 = 7
	irn13 = 8
	irn15c = 9
	irn15o = 10
	irn16 = 11
	irn17 = 0
      else if (irnrat.eq.5) then
	krnrat = 12
	irn12 = 7
	irn13 = 8
	irn15c = 9
	irn15o = 10
	irn16 = 11
	irn17 = 12
      else
	write(istdou,205) irnrat
	if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,205) irnrat
	stop 'rnrate'
      end if
c
c  set derivatives of log ne into anel
c
      secder=.not.nosd
      ii=4
      amne=amm*ane(1)
      do 3 i=2,4
      anel(i)=ane(i)/amne
      if(secder) then
        do 2 j=i,4
        ii=ii+1
    2   anel(ii)=(ane(ii)-ane(i)*(ane(j)/ane(1)))/amne
      end if
    3 continue
c
c  reactions between heavy nuclei
c
      t9=t/1.d9
      t13=t9**0.3333333333333333d0
      tpw(1)=t13
      do 5 l=2,5
    5 tpw(l)=t13*tpw(l-1)
      t62= sqrt(t)/1000
      t623=t62*t62*t62
      rh=rho(1)
c  test for reactions
      if(t13.le.0.1) then
        norct=.true.
        uw(1)=0
        return
      end if
c
      norct=.false.
c
      if(idgeng.ge.2.and.istdpr.gt.0) then
        write(istdpr,210) t, x, y, z, rh
      end if
c
c  initialize internal test for no reactions
c
      norct1=.true.
c
c  range of derivatives
c
      if(secder) then
        ider=10
      else
	ider=4
      end if
c
c  electron screening
c
      call elscrn(fl, tl, x, y, z, iscren, zt, uw, zzscr, thte, nosd)
c
c  ********************************************
c
c
c  calculate nuclear rates
c
c Not NACRE reaction rates
c
      if (ivreng.ne.8) then
c
        do 50 k=1,krnrat
c
c  skip reaction Be7(e-,nu)Li7
c
        if(k.eq.6) go to 50
c
c  set index for reaction parameters
c
        kpar = iptrnr(k,irnrat)
c
        if(idgeng.ge.2.and.istdpr.gt.0) then
          write(istdpr,221) k, kpar, a(kpar), b(kpar), zzscr(k), t13
        end if
        exx=b(kpar)/t13
c
c  test for no reactions
c
        if(a(kpar).le.0) stop 'rnrate'
        aexx=log(a(kpar))-exx
        if(aexx.le.-coexp) then
c
c  zero reaction rates and logarithmic derivatives
c
c  Loop corrected 7/8/05 by M. Bazot, from al(1,k) = 0
c
          do 36 i=1,ider
   36     al(i,k)=0
          go to 50
c
        end if
c
c  set reaction rates
c
        norct1=.false.
c
        sum1=1
        sum2=0
        sum3=0
        do 38 l=1,5
        t1=s(kpar,l)*tpw(l)
        sum1=sum1+t1
        sum2=sum2+l*t1
   38   sum3=sum3+l*l*t1
        sum2=sum2/3
        sum3=sum3/9
c
        if(idgeng.ge.2.and.istdpr.gt.0) then
          write(istdpr,222) k, aexx, sum1, sum2, sum3, tpw(2)
        end if
c
        alt=fxp(aexx)*sum1/tpw(2)
        dalt=exx/3-0.6666666666666667d0+sum2/sum1
        if(secder) then
          ddalt=-amm*(exx/9+(sum2*sum2/sum1-sum3)/sum1)
        end if
c
        al(1,k)= exp(zzscr(k)*uw(1))*alt
        do 45 i=2,ider
   45   al(i,k)=zzscr(k)*uw(i)
        al(3,k)=al(3,k)+dalt
        if(secder) then
          al(8,k)=al(8,k)+ddalt
        end if
c
        if(idgeng.ge.2.and.istdpr.gt.0) then
          write(istdpr,230) k, al(1,k)
        end if
   50   continue
      else
c
c Nacre reaction rates
c
c Test
c
        if(idgeng.ne.0.and.istdpr.gt.0) 
     *     write(istdpr,*) 'Call s/r cnrnacre'
c
c End test
c
        call cnrnacre(tl,ider,norct1,secder,zzscr) 
c
      end if
c
c Test
c
      if(idgeng.ne.0.and.istdpr.gt.0) then
	write(istdpr,*) 'Reaction rates in rnrate'
        do j=1,krnrat
          write(istdpr,*) al(1,j)
        enddo
        write(istdpr,*) 'norct1 = ',norct1
c
c End test
c
      end if
c
c  test whether there were actually any reactions. otherwise
c  skip Be7+e- reaction.
c
      if(norct1) then
        norct=.true.
        uw(1)=0
        return
      end if
c
c  the Be7+e- reaction, From Bahcall and Moeller (ApJ, vol. 155, 511)
c
      iscbe7=mod(iscren/10,10)
c
      if(iscren.eq.0) then
	f(1)=1
	go to 65
      else if(iscbe7.eq.1) then
	f(1)=1.2
	go to 65
      end if
c
c  full Bahcall and Moeller screening expression
c
      rs(1)=t62* sqrt(2.845d0/rh)/zt(1)
      r=rs(1)
      amrs=amm*r
      do 52 i=2,4
   52 rs(i)=amrs*(-rho(i)/2-zt(i))
      rs(3)=rs(3)+amrs/2
      sgr=-0.431d0+r*(2.091d0+r*(0.401d0*r-1.481d0))
      cr2=-0.6064d0+r*(4.859d0+r*(1.907d0*r-5.283d0))
      dsgr=2.091d0+r*(1.203d0*r-2.962d0)
      dcr2=4.859d0+r*(5.721d0*r-10.566d0)
      xhi(1)=-7.35d5
      xhi(2)=2.515d6
      ix=0
      do 58 i=1,11,10
      ix=ix+1
      xx=xhi(ix)
      ee(i)=fxp(xx*sgr/t)
      fct=ee(i)*xx/t
      do 55 k=1,3
   55 ee(i+k)=fct*dsgr*rs(k+1)
   58 ee(i+2)=ee(i+2)-amm*fct*sgr
c
c  test for very large lr
c
      if(ee(11).lt.1.d17) go to 59
c
c  set approximate screening and first derivatives
c
      iapscr=1
      allr=0.246d0*rh*(ane(1)/av)/t623
      fscr=8.822d6/(allr*t)
      f(1)=1+fscr*cr2
c
      fscr=fscr/(amm*f(1))
      do 58100 i=2,4
58100 f(i)=fscr*(dcr2*rs(i)-cr2*amm*(rho(i)+anel(i)))
      f(3)=f(3)+0.5d0*amm*fscr*cr2
      go to 67
c
c  set full expression for screening
c
   59 iapscr=0
c
      alr(1)=0.246d0*rh*(ane(1)/av)*ee(11)/t623
      allr=alr(1)
      amlr=amm*allr
      do 60 i=2,4
   60 alr(i)=amlr*(rho(i)+anel(i)+ee(i+10)/ee(11)/amm)
      alr(3)=alr(3)-1.5d0*amlr
      ddd=1+allr*(1+allr*ee(1)/4)
      dd(1)=ddd
      alr2=allr*allr
      alr1=1+allr*ee(1)/2
      do 62 i=2,4
   62 dd(i)=alr(i)*alr1+alr2*ee(i)/4
      fsr=1+0.435d0*allr*ee(1)
      ssr=cr2*fsr/ddd
      sr(1)=ssr
      dlcr=dcr2/cr2
      fsr=0.435d0/fsr
      do 63 i=2,4
   63 sr(i)=ssr*(dlcr*rs(i)+fsr*(alr(i)*ee(1)+allr*ee(i))
     .  -dd(i)/ddd)
      fcf=5.07d6/t
      fcf1=fcf*ssr*ee(11)
      f(1)=1+fcf1
      amff=amm*f(1)
      do 64 i=2,4
   64 f(i)=fcf*(ee(10+i)*ssr+ee(11)*sr(i))/amff
      f(3)=f(3)-amm*fcf1/amff
      go to 67
c
c  no screening
c
   65 do 66 i=2,10
   66 f(i)=0
c
c  set lambda c
c
   67 flc=1+4*(t9-0.016d0)
      allc=4.62d-9*rh*ane(1)/av/t62*flc
      alc(1)=allc
      do 68 i=2,4
   68 alc(i)=rho(i)+anel(i)
      alc(3)=alc(3)+4*t9/flc-0.5d0
      al(1,6)=allc*f(1)
      do 70 i=2,4
   70 al(i,6)=alc(i)+f(i)
c
      if(idgeng.ge.2.and.istdpr.gt.0) then
        write(istdpr,240) al(1,6)
      end if
c
c
      if(nosd) return
c
c  second derivatives
c
      if(iscren.eq.0.or.iscbe7.eq.1) go to 97
      ddsgr=-2.962d0+2.406d0*r
      ddcr2=-10.566d0+11.442d0*r
      do 72 i=2,4
   72 st1(i)=rs(i)/amrs
      ii=4
      do 74 i=2,4
      do 74 j=i,4
      ii=ii+1
   74 rs(ii)=amrs*(amm*st1(i)*st1(j)-zt(ii)-rho(ii)/2)
c
c  test for low-temperature approximation
c
      if(iapscr.ne.1) go to 75
      ii=4
      do 74100 i=2,4
      do 74100 j=i,4
      ii=ii+1
      f(ii)=f(i)*f(j)
      if(i.eq.3) f(ii)=f(ii)+f(j)
      if(j.eq.3) f(ii)=f(ii)+f(i)
c
c
c  note that the inclusion of the effect of lr tilde is not
c  strictly correct, but this is unlikely to have a significant
c  effect.
c
      f(ii)=-amm*f(ii)+fscr*(ddcr2*rs(i)*rs(j)+dcr2*rs(ii))-
     *  anel(ii)-rho(ii)
c
74100 continue
      f(8)=f(8)-fscr*amm2*cr2
      go to 97
c
c  set second derivative of full expression
c
   75 do 76 i=2,4
   76 st1(i)=amm*dsgr*rs(i)
      ix=0
      do 84 k=1,11,10
      k1=k-1
      ix=ix+1
      fc2=xhi(ix)/t
      fcx=ee(k)*fc2
      fc2=fc2*fcx
      do 78 i=2,4
   78 st2(i)=ee(k1+i)/fcx
      ii=4
      do 82 i=2,4
      do 82 j=i,4
      ii=ii+1
      dde=ddsgr*rs(i)*rs(j)+dsgr*rs(ii)
      if(i.eq.3) dde=dde-st1(j)
      if(j.eq.3) dde=dde-st1(i)
   82 ee(k1+ii)=fc2*st2(i)*st2(j)+fcx*dde
   84 ee(k+7)=ee(k+7)+amm*amm*fcx*sgr
      do 86 i=2,4
   86 st1(i)=alr(i)/amlr
      ii=4
      do 88 i=2,4
      do 88 j=i,4
      ii=ii+1
      alr(ii)=amlr*(amm*st1(i)*st1(j)+rho(ii)+anel(ii)+(ee(10+ii)-
     .  ee(10+i)*ee(10+j)/ee(11))/ee(11)/amm)
   88 dd(ii)=alr1*alr(ii)+ee(1)*alr(i)*alr(j)/2+allr*(alr(i)*ee(j)
     .  +alr(j)*ee(i))/2+alr2*ee(ii)/4
      do 90 i=2,4
      st1(i)=sr(i)/ssr
   90 st2(i)=ee(1)*alr(i)+ee(i)*alr(1)
      cdcr=(ddcr2-dcr2*dcr2/cr2)/cr2
      ii=4
      do 92 i=2,4
      do 92 j=i,4
      ii=ii+1
   92 sr(ii)=ssr*(st1(i)*st1(j)+cdcr*rs(i)*rs(j)+dlcr*rs(ii)
     .  +fsr*(alr(ii)*ee(1)+alr(i)*ee(j)+alr(j)*ee(i)+allr*ee(ii)
     .  -fsr*st2(i)*st2(j))+(dd(i)*(dd(j)/ddd)-dd(ii))/ddd)
      do 94 i=2,4
   94 st1(i)=f(i)*amff
      ii=4
      do 96 i=2,4
      do 96 j=i,4
      ii=ii+1
      ddf=-st1(i)*st1(j)/f(1)+fcf*(ee(10+ii)*ssr+ee(10+i)*sr(j)
     .  +ee(10+j)*sr(i)+ee(11)*sr(ii))
      if(j.eq.3) ddf=ddf-amm*st1(i)
      if(i.eq.3) ddf=ddf-amm*fcf*(ee(11)*sr(j)+ee(10+j)*ssr)
   96 f(ii)=ddf/amff
c
   97 ii=4
      do 98 i=2,4
      do 98 j=i,4
      ii=ii+1
   98 alc(ii)=rho(ii)+anel(ii)
      alc(8)=alc(8)+3.744d0*amm*t9/(flc*flc)
      ii=4
      do 100 i=2,4
      do 100 j=i,4
      ii=ii+1
  100 al(ii,6)=alc(ii)+f(ii)
c
      f(1)=log10(f(1))
      alc(1)=log10(alc(1))
      return
  205 format(//' **** Error in s/r rnrate. irnrat = ',i4,
     *  ' illegal')
  210 format(' diagnostics from s/r rnrate. T, X, Y, Z, rho =',
     *  1p5e13.5)
  220 format(' electron screening with zeta =',1pe13.5,
     *  '  Uw = ',e13.5)
  221 format(' start setting lambda for k =',i3,'  kpar =',i3/
     *  ' a(k), b(k), zzscr(k), t13 =',1p4e13.5)
  222 format(' k =',i3,'  aexx, sum1, sum2, sum3, tpw(2) =',
     *  1p5e13.5)
  230 format(' reaction rate no',i3,' set to',1pe13.5)
  240 format(' Be7 + e- reaction rate =',1pe13.5)
      end
      subroutine dmprcn
c
c  for testing, dump constants in common/rcncns/
c
      implicit double precision(a-h, o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.nnz.d.incl'
c
      common/rcncns/ a(knucst),b(knucst),s(knucst,isnuc),
     *  zz(knucst),zz12(2,knucst),zsmh,acno,awght(knucwg),q(knucq), 
     *  knucpm, kqvals, kawght
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      write(istdpr,100) (a(i),i=1,knucpm)
      write(istdpr,110) (b(i),i=1,knucpm)
      do 15 k=1,knucst
   15 write(istdpr,120) k,(s(k,i),i=1,isnuc)
      write(istdpr,130) (zz(i),i=1,knucpm)
      write(istdpr,135) ((zz12(i,j),j=1,2),i=1,knucpm)
      write(istdpr,140) zsmh, acno
      write(istdpr,150) (awght(i),i=1,knucpm)
      write(istdpr,160) (q(i),i=1,knucpm)
      write(istdpr,170) knucpm, kqvals, kawght
      return
  100 format(' Constants in common/rcncns/'/
     *  ' a ='/(1p5e13.5))
  110 format(' b ='/(1p5e13.5))
  120 format(' s(',i2,') ='/(1p5e13.5))
  130 format(' zz ='/(1p5e13.5))
  135 format(' zz:'/(1p2e13.5))
  140 format(' zsmh, acno =',1p2e13.5)
  150 format(' awght ='/(1p5e13.5))
  160 format(' q ='/(1p5e13.5))
  170 format(' knucpm, kqvals, kawght =',3i5)
      end
      subroutine engcse(ieqhe3, icnocs, iheccs, nspec, nspxx3, nspcno,
     *  nsphec, nspect, irnrat,cxnuc)
c
c  sets number nspec of elements considered in hydrogen burning
c  and case number irnrat for
c  s/r rnrate, depending on flags ieqhe3 and icnocs for
c  treatment of He3 and CNO cycle. Also sets character array
c  cxnuc(k), k = 1, ..., nspec, with names of elements.
c  Finally sets variables ixc12, ..., ixo17 to point to the
c  respective elements in the array x(k). If the element is
c  not included, the corresponding ix is set to 0.
c
c  If iheccs ne 0 also include helium burning. In this case the total 
c  of composition variables is set to nspect = nspec + 2 (including
c  also, as separately treated elements, 4He and 12C 
c  NOTE: THIS MAY NEED LATER CHANGE
c  Otherwise, nspect is returned equal to nspec.
c
c  Original version: 26/7/91.
c
c  Modified 1/8/02, allowing for helium burning.
c
      implicit double precision (a-h, o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.nnz.d.incl'
c
      character*(*) cxnuc
      dimension cxnuc(*)
      dimension ixx(nspcmx)
      common/cengcs/ ixc12, ixc13, ixn14, ixn15, ixo16, ixo17
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      equivalence(ixc12, ixx(1))
c
c  initialize element pointers to zero
c
c..      call zero(ixx, nspcmx)
      do 10 i=1,nspcmx
   10 ixx(i)=0
c
      if(ieqhe3.eq.0.or.ieqhe3.eq.2) then
	nspxx3 = 2
	cxnuc(1)='H1'
	cxnuc(2)='He3'
      else if(ieqhe3.eq.1) then
	cxnuc(1)='H1'
	nspxx3 = 1
      else 
	write(istdou,110) ieqhe3
	if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,110) ieqhe3
	stop 'engcse'
      end if
c
      if(icnocs.eq.0) then
	nspcno = 0
	irnrat = 1
      else if(icnocs.eq.1) then
	nspcno = 1
	irnrat = 2
	cxnuc(nspxx3+1)='N14'
	ixn14 = nspxx3+1
      else if(icnocs.eq.2) then
	nspcno = 1
	irnrat = 4
	cxnuc(nspxx3+1)='N14'
	ixn14 = nspxx3+1
      else if(icnocs.eq.3) then
	nspcno = 2
	irnrat = 3
	cxnuc(nspxx3+1)='C13'
	cxnuc(nspxx3+2)='N14'
	ixc13 = nspxx3+1
	ixn14 = nspxx3+2
      else if(icnocs.eq.4) then
	nspcno = 3
	irnrat = 4
	cxnuc(nspxx3+1)='C12'
	cxnuc(nspxx3+2)='C13'
	cxnuc(nspxx3+3)='N14'
	ixc12 = nspxx3+1
	ixc13 = nspxx3+2
	ixn14 = nspxx3+3
      else
	write(istdou,120) icnocs
	if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,120) icnocs
	stop 'engcse'
      end if
c
      nspec = nspcno + nspxx3
c
c  total number of elements, possibly including contribution from 
c  helium burning
c
      if(iheccs.ne.0) then
	nspect = nspec + 2
	nsphec = 2
      else
	nspect = nspec
	nsphec = 0
      end if
      if(idgeng.ne.0.and.istdpr.gt.0) 
     *  write(istdpr,*) 'nspcno, nspxx3, nspec, nspect',
     *  nspcno,nspxx3,nspec,nspect
c
      return
  110 format(//' ***** Error in s/r engcse. ieqhe3 =',i5,
     *  '  not implemented.')
  120 format(//' ***** Error in s/r engcse. icnocs =',i5,
     *  '  not implemented.')
      end
      subroutine epsjop(epsj,tl)
c
c  sets artificial increase in energy generation rate.
c  Controlled by parameters in common/engfdg/
c  ifdeps = 1: Increase epsilon by epsfdg for qepsf1 le (m/M) le qepsf2
c  ifdeps = 2: Increase epsilon by epsfdg*(T/Tc)
c              for qepsf1 le (m/M) le qepsf2
c
c  Original version: 12/7/96 (based on earlier routine written with
c  JOP)
c
      implicit double precision(a-h,o-z)
      include 'engenr.nnz.d.incl'
      parameter(idr1mx = nspcmx+3, naztmx = nspcmx+3, 
     *  nbmax = nspcmx + 2)
c..      character*80 file
      common/crhsvr/ qx
c..      common/cofile/ nfiles, idsfil(20), file(20),iopen(20)
      common/ksider/ al0,al2,aztst(naztmx),axst(naztmx)
      common/engfdg/ epsfdg, qepsf1, qepsf2, ifdeps
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
      if(ifdeps.eq.1) then
        if(qx.ge.qepsf1.and.qx.le.qepsf2) then
	  epsj=epsfdg
        else
	  epsj=0
        end if
      else
	tc=1.d7*aztst(2)
	t=10.d0**tl
	if(tc.ge.t) then 
	  tfact=t/tc
        else
	  tfact=1
        end if
        if(qx.ge.qepsf1.and.qx.le.qepsf2) then
	  epsj=epsfdg*tfact
        else
	  epsj=0
        end if
      end if
      return
      end
