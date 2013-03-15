      subroutine mshstr(x,y,in1,in,iy,nnin,nn,rhs,resrc)
c
c  routine for resetting mesh in evolution calculation
c
c  note: nnin is number of mesh points in input model, and
c  nn is desired number of mesh points in stretched model.
c
c  modified 20/8/1984 to take into account rescaling of variables
c  in central boundary conditions.
c
c  modified 6/9/1984 to correct for underflow in second derivative
c  of dgrad.
c
c  modified 29/12/1984 to use expansion at central meshpoint
c  for x and x3 in all cases
c
c  modified 4/1/1985 to use log(r/1.d11) and log(l/1.d33) as
c  dependent variables.
c
c  modified 7/2/1986 to allow choice between old and new
c  weighting of term in second derivative of dgrad,
c  as determined by iwdgrd 
c  
c  21/9/87: implementing modifications from RECKU
c
c  22/1/88: add output to file ttt.mshstr.n.out when iprt .ne. 0
c    (changed to /scratch/jcd/evolprg/ttt.mshstr.out on 30/5/92)
c
c  29/2/88: check for size of X in determining central radius
c
c  21/3/88: add nnin as argument.
c
c  30/8/88: zero term in X3 near transition point, where X3
c  becones non-zero, to avoid effect of transient in derivative.
c
c  1/3/89:  test for whether expansion coefficients of m are set
c  (they are not in stretching of trial solution). If not, use
c  approximate expressions based on ideal gas and homogeneous
c  composition.
c
c  11/8/90: Introduce factor to reduce effect of rX on mesh when
c  X is very small. For the time being size of factor is determined
c  by the constant xmshcs hardcoded below.
c
c  13/8/90: Introduce option iwdgrd = 3 to switch off contribution from
c  second derivative of superadiabatic gradient.
c
c  16/8/90: Introduce test for very large increment between consecutive
c  meshpoints. If so, use linear interpolation, and repeat stretch
c  on new mesh.
c
c  11/11/90: Introduce term in X in setting am2 for trial model
c
c  11/4/91: Introduce option for increasing density of points
c     at bottom of convective envelope, to resolve overshooting
c
c  2/8/91: Modified for consistency with new energy generation etc.
c     involving CNO consistently.
c
c  10/8/91: Modified to reset nn1 if stretch is repeated.
c
c  7/5/92: Modified to reset boundary of possible convective core,
c     and to use linear interpolation in the vicinity of core
c     boundary
c
c  8/5/92: Modified to include option of extra points at convective
c     core boundary.
c
c  4/6/92: Modified to do stretch at boundary of mixed core in
c     terms of log q.
c
c  8/6/92: Modified to include term in gradient in hydrogen abundance
c     in stretch. Also reorganize storage, update documentation.
c
c  7/12/92: Modified setting of luminosity at central mesh point after
c     resetting rc. If al0 is not set, and for extrapolated model
c     now scale luminosity by mass ratio. Previously, no change.
c
c  10/12/92: modified to set X3 at new mesh point by extrapolation,
c     when central derivative is zero.
c
c  18/3/94: modified to take out dependence on X and rX when wxh .lt. 0
c
c  2/8/96: Modified to restrict increase in innermost mesh point
c  to be at most 10 per cent in radius.
c
c  Modified 19/7/95, to use general heavy-element abundance
c  as set in common/heavy/.
c
c  Modified 4/8/96, taking out option of writing data to disk with
c  ifwrt = 1 (but retaining ifwrt in common/rhcn/ for consistency)
c
c  Corrected 5/8/97, by including engenr.nnz.d.incl, to set nspdmx etc
c
c  Modified 24/9/98, including a test for non-monotonicity and
c  correction of mesh, if needed, in s/r testmono.
c
c  Modified 16/11/99, blocking too rapid decrease in central radius,
c  and resetting luminosity setting in core, when hydrogen is exhausted
c
c  Modified 15/5/00, including istrtc to flag for resetting
c  of composition variables in cvr.
c
c  Modified 10/8/02, starting implementation of 4He burning
c
c  Modified 15/8/02, correcting setting of limits for special treatment 
c  at edge of convective core from
c       aksc1=aks(nmxcor)-2
c       aksc2=aks(nmxcor)+2
c  (incorrect) to
c       aksc1=aks(nmxcor-2)
c       aksc2=aks(nmxcor+2)
c
c  Also introduced option istrt1 = 2, to interpolate/extrapolate 
c  linearly on either side of convective-core boundary.
c
c  Modified 29/8/02, allowing effects of several discontinuities in
c  composition, etc., and modifying the definition of discontinuities
c  in s/r resxdr
c
c  Modified 30/8/02, reducing effects of very small convective cores.
c
c  Modified 10/10/02, correcting testmono so that it actually resets
c  mesh in case of non-monotonicity
c
c  Modified 16/10/02 in s/r resxdr, changing definition of discontinuity
c  to be more sensitive, and introducing rescaling of reset derivative
c  to ensure same trapezoidal integral. Also forcing use of trapezoidal
c  integration (kvint = 1) in setting aks.
c
c  Modified 18/10/02, including epsg in resetting
c
c  Modified 12/5/03, using linear interpolation near base of convective
c  envelope when istrt2 .ge.1 and with diffusion.
c
c  Modified 5/6/03, forcing continuation even with excessive no. of 
c  repeats, by hard-coded flag ignrpt.
c
c  Modified 26/6/03, resetting zh(.) before a possible repeat.
c  Also, include zh in the basic interpolation, to ensure that it
c  is always consistent.
c
c  Modified 21/8/03, changing treatment of discontinuties in
c  s/r resxdr.
c
c  Modified 22/8/03, enforcing linear interpolation across largest
c  jump in aks if the jump exceeds 3.
c
c  Modified 20/10/03, correcting setting of composition variables
c  in convective core (previously, Y_j was erroneously set to central
c  value in case with diffusion).
c
c  Modified 15/11/03, correcting resetting of extrapolated solution
c  by mass scaling
c
c  Modified 18/10/05, to use negative wmshos for skewed setting of
c  mesh points at edge of convective envelope (or overshoot region)
c
c  Modified 19/10/05, to force linear interpolation wherever the
c  interval in aks exceeds dakslin, if i_dakslin is 1. So far
c  i_dakslin is hardcoded in the routine, and dakslin is set to 6.
c
c  Modified 15/3/06, correcting setting the position of the lower
c  boundary of convection zones (qll(i))
c
c  Modified 6/8/10 initializing qmxmin and qmxmax for time0
c
c
c  Notes on variables etc in stretch:
c  ---------------------------------
c  
c  The overall strategy is to stretch on the basis of derivatives of
c  various quantities. Currently the following variables are used
c  in the routine (and stored in the array trm, etc):
c  
c   Element    variable        weight
c  
c     1           r             wr
c     2         log p           wpt
c     3         log T           wpt
c     4         log L            1
c     5           Xh            wxh
c     6           rXh            1?
c     7           X3            wx3
c     8          ddad           wdgr
c     9        grad(ddad)       wdgr
c  
c  For further details of the precise form of the corresponding terms,
c  see source code and CD82. In addition extra contributions may
c  be added at, e.g., the boundary of a mixed region or the base of 
c  the convection zone.
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
      implicit double precision (a-h,o-z)
      logical time0,resrc,conv,convc,nosd,notd,theccs
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
c
      parameter(idr1mx = nspcmx + 3, 
     *  iyint = 2*ivarmx+1,
     *  nbmax = nspcmx+2, naztmx = nspcmx + 3,
     *  itrmmx=15)
c
      dimension x(*),y(iy,*)
      dimension zk(1),ap(1),aq(1),
     *  trmin(itrmmx),trmax(itrmmx),rg(itrmmx),trm(itrmmx),
     *  weight(itrmmx),z(itrmmx),
     *  f(ivarmx),fd(ivarmx,1),alam(1,ivarmx),alamd(1,ivarmx,1),
     *  h(1),hd(1,1),zrhs(ivarmx),dzdy(ivarmx,ivarmx),
     .  yint(iyint),azt1(5),ax1(5),qlf(6),qll(6),dgnout(30,nnmax),
     *  nrxdrc(20),xrxdrc(20),ncdisc(20),xcdisc(20),aksms(2,20),
     *  xspcpt(100)
c
      common/work/ yw(iywstr,1)
      common/sooner/ aks(1)
      common/ln10/ amm
      common/cmtime/ age, time0
      common/step/ dt
      common/noiter/ iter, ntime
      common/heavy/ zatmos, zhc, zh(1)
      common/he3fdg/ agesh,ifdhe3
      common/he3eql/ xhe3eq,ieqhe3,iequsd
      common/engcnt/ xhzlm1, xhzlm2, nxhzer, icnocs, nspec,
     *  iheccs, nspect
      common/engche/ xmxrhe
      common/compvr/ cvr(icvrmx,1)
      common/cstcvr/ icvh1, icvhe4, icvhe3, icvc12, icvc13, icvn14, 
     *  icvo16
      common/cstyst/ iyh1, iyhe3, iyc12, iyc13, iyn14, iyo16, 
     *  iyche4, iycc12
      common/cstazt/ iah1, iahe3, iac12, iac13, ian14, iao16, 
     *  iache4, iacc12
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec, idcomp, iccomp
      common/eqstd/ dumeq(54),pt(20),cp(4),dad(4)
      common/cmpder/ ft(idr1mx,nspcmx)
      common/excf/ am0,al0,am2,al2
      common/rhcn/ a1,a2,a3,a4,zdummy,nvar,ifwrt
      common/enggrv/ epsg(nnmax)
      common/clshft/ alshft
      common/bccn/ b1
      common/ksider/ dm1,dm2,aztst(naztmx),axst(naztmx)
      common/caztax/ azt(nbmax),ax(nbmax)
      common/cntmsh/ wx,epsr,wr,wpt,wxh,wx3,wdgr,
     *  wmshos,drmsos,wmshcc,dqlmsc,qmscmn,qmscmx,
     *  koint,kvint,iprt,icngr,nnt,iwdgrd,istrtc
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,jnc,frcf(6),frcl(6)
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor,nmxcp,qmxcp,rmxcp,
     *  qmxmin,qmxmax,qmxscn,qmxscp,nmxscn,nmxscp,icjump
      common/cdiffu/ idiffus, itbdif, rhob, rmxdif, tlfdif, dtldif,
     *  ctdfmx, rctdif, rctdf1, idfxbc, ismdif, nsmdif
      common/cnvout/ ddrad,rr,a,ddacad
      common/cnvovs/ icnvos, jcnvos, clcovs, cldovs, alphos,
     *  rczl, rczlfx, rcnvos, qlcnos
      common/prtvar/ rhl,ak,akt,akp,akx,eps(idr1mx),dr,dac,conv
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3,idgmxc,
     *  idgsdt,idgn75,idgn76,idgn77
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      data iopnfl /0/
      data iopn75 /0/
c
c  constant for mesh contribution from rX
c
      data xmshcs / 0.05/
c
c  number of special points
c
      data nspcpt /0/
c
c  initial limits for linear interpolation in connection with
c  convective core
c
      data qmxmax, qmxmin /-1.d0, -1.d0/
c
c  limit for increment between meshpoints
c#ai#  Note: this has been set unreasonably high to make the
c  programme go, with a growing convective core. Has to
c  be fixed, somehow.
c
      data xsimax /20.d0/
c
      external rhs
c
      save
c
      if(istdpr.gt.0) write(istdpr,'(/'' Enter mshstr'')')
c
c  Warning: hardcode ignrpt to ignore failure after several repeats
c
      ignrpt=1
      write(6,*) '#D# (1) xsimax =', xsimax
c
c  hard-coded range of running-mean smoothing and Gaussian replacement
c
      nsmth=21
      deltax=1.e-3
c
c  hardcode reset to trapezoidal integration
c  (Implemented explicitly instead of call to vintk below)
c
      kvint=1
c
c  hardcode forcing linear interpolation in all large intervals of
c  aks
c
      i_dakslin=0
      dakslin=6.
c
c  internal printing flag (might be reset in case of errors)
c
      iprtlc=iprt
c
c  test for re-initializing qmxmin and qmxmax
c
      if(time0) then
        qmxmin=-1.d0
        qmxmax=-1.d0
        if(istdpr.gt.0) write(istdpr,
     *    '(/'' Initialize qmxmin, qmxmax to -1.d0''/)')
      end if
c
      istrt0=mod(istrtc,10)
      istrt1=mod(istrtc/10,10)
      istrt2=mod(istrtc/100,10)
      istrt3=mod(istrtc/1000,10)
      istrt4=mod(istrtc/10000,10)
      istrt5=mod(istrtc/100000,10)
c
      if(istdpr.gt.0) then
        write(istdpr,*) 
     *    'wx,epsr,wr,wpt,wxh,wx3,wdgr,wmshos,drmsos,wmshcc,dqlmsc,',
     .    'koint,kvint,iprt,',
     .    'icngr,nnt,iwdgrd,nvar,in1,in,ivarmx'
        write(istdpr,*) 
     *    wx,epsr,wr,wpt,wxh,wx3,wdgr,wmshos,drmsos,wmshcc,dqlmsc,
     .    koint,kvint,iprt,
     .    icngr,nnt,iwdgrd,nvar,in1,in,ivarmx
      end if
c
c  store internal version of rczl, before starting 
c  calls of rhs
c
      rczli=rczl
      if(istdpr.gt.0) write(istdpr,*) 'Setting rczli to',rczli
c
c  set weight for hydrogen abundance to zero, when wxh .lt. 0
c
      if(wxh.ge.0) then
	wxhc=wxh
	frxfac=1
      else
	wxhc=0
	frxfac=0
	if(istdpr.gt.0) write(istdpr,*) 
     *    'Zero weights for hydrogen abundance'
      end if
c
c  test for setting qmxmin and qmxmax from input
c
      if(qmscmn.gt.0) qmxmin=qmscmn
      if(qmscmx.gt.0) qmxmax=qmscmx
c
c  set flag for convectively mixed core
c
      if(nmxcor.gt.0) then
	convc=.true.
      else
	convc=.false.
      end if
c
c  set flag for 4He burning (note hard-coded temperature limit)
c
      theccs=iheccs.ne.0.and.y(in1+2,nnin).gt.7.78
      if(istdpr.gt.0) then
	write(istdpr,*) 
     *    'In mshstr, iheccs, nnin, y(in1+2,nnin), theccs =',
     *    iheccs, nnin, y(in1+2,nnin), theccs
c
        write(istdpr,102) koint,kvint
      end if
c
      iyw=iywstr
c
c  set number of numerical derivatives
c
      if(theccs) then
	knder = 7
      else 
	knder = 5
      end if
c
      in11=in1-1
c
c  set input number of mesh points
c
      nn1=nnin
c
      ifdho=ifdhe3
      ifdhe3=0
      if(icngr.eq.1) then
c  congruent mesh with nn points
        do 3 n=1,nn1
    3   yw(2,n)=n-1
	kointc=koint
        go to 56
      end if
      write(6,*) '#D# (2) xsimax =', xsimax
c
c  test for redetermination of rc
    1 if(.not.resrc) go to 9
c
c  determine central r
c
    2 call store(aztst,azt1,3)
      call store(axst,ax1,3)
c
c  exclude contribution from X when X .lt. 0.001
c  or for initial model
c
      if(time0.or.aztst(3).lt.0.001) then
        azt1(3)=1
        ax1(3)=0
      end if
c
c  test if am0 and am2 are set
c
      if(am0.eq.0.or.am2.eq.0) then
c
c  set rhoc from xi(1). Note: here we use G = 6.67232e-8,
c  and assume the ideal gas law
c
        rhoc=sqrt(-ax1(1)/(1.397447e-2))
        am0 =4.188790*rhoc
        am2 =0.6*am0*(ax1(1)/azt1(1)-ax1(2)/azt1(2)
     *       -5*ax1(3)/(3+5*azt1(3)-zhc))
        if(istdpr.gt.0) write(istdpr,103) rhoc,am0,am2
      end if
c
      azt1(4)=am0
      ax1(4)=am2
c
c  test for excluding luminosity term
c
      if(al0.eq.0.or.istrt4.ge.1) then 
	if(istdpr.gt.0) then
          if(al0.eq.0) write(istdpr,104)
          if(istrt4.ge.1) write(istdpr,
     *      '(//'' In s/r mshstr, istrt4 ='',i3,
     *       '' Exclude luminosity term''/)') istdt4
        end if
	azt1(5)=1
	ax1(5)=0
      else
        azt1(5)=al0
        ax1(5)=al2
      end if
c
      sum=0
      do 5 i=1,5
      ttrm= abs(ax1(i)/azt1(i))
      trm(i)=ttrm
    5 sum=sum+ttrm
      if(istdpr.gt.0) write(istdpr,105) (trm(i),i=1,5),sum
      rhc=epsr/ sqrt(sum)
      rc=1.d11*rhc
c
c  test that size of central region is not increasing too fast
c
      rcold=1.d11*10.**y(in1,nn)
      if(rc.gt.1.1*rcold) then
        if(istdpr.gt.0) write(istdpr,106) rc, rcold
        rc=1.1*rcold
        rhc=1.d-11*rc
      end if
c
c  do not allow central radius to be reduced by more than 10 per cent
c  when hydrogen is exhausted in the core (a little artificial,
c  admittedly)
c
      rhcp=10.d0**y(in1,nn1)
      if(y(in1+4,nn1).le.1.e-5.and.rhc.lt.0.9*rhcp) then
	if(istdpr.gt.0) write(istdpr,107)
	rhc=0.9*rhcp
      end if
c
      rc=1.d11*rhc
c
c  remove points interior to new central meshpoint
c
      rhcl=log10(rhc)
    6 nn1=nn1-1
      if(y(in1,nn1).le.rhcl) go to 6
      nn1=nn1+1
      write(6,*) '#D# (3) xsimax =', xsimax
c
c  reset x and y at central meshpoint, using expansion
c  Ensure that log q is still monotonic (this may be a problem when
c  stretching trial model).
c
      rhc2=rhc*rhc
      rhc3=rhc*rhc2
      xnew=log10(rhc3*(am0+rhc2*am2))-b1
      if(xnew.lt.x(nn1-1)) then
	delqlc=xnew-x(nn1)
        x(nn1)=xnew
      else
	delqlc=0
	if(istdpr.gt.0) write(istdpr,110)
      end if
      if(istdpr.gt.0) write(istdpr,109) rc, xnew, nn1
c
c  set coefficients for extrapolation, if needed
c
      rh12=10.d0**(2*y(in1,nn1))
      rh22=10.d0**(2*y(in1,nn1-1))
      fcct1=(rh22-rhc2)/(rh22-rh12)
      fcct2=1-fcct1
c
      y(in1,nn1)=rhcl
      tl=7+log10(aztst(2)+rhc2*axst(2))
c
c  test for zero al0. If so, reset by mass scaling
c
      if(al0.le.0.or.y(in1+4,nn1).le.1.e-5) then
	if(istdpr.gt.0) write(istdpr,111) al0, y(in1+4,nn1)
        y(in1+3,nn1)=
     *    log10((10.d0**y(in1+3,nn1)-alshft)*10.d0**delqlc+alshft)
      else
        y(in1+3,nn1)=log10(rhc3*(al0+rhc2*al2)+alshft)
      end if
c
c  as a temporary fix, set composition variables by extrapolation
c  if the central derivatives are set to zero.
c  This needs attention. Also, should be fixed for CNO abundances
c
      if(axst(3).ne.0) then
        xh=aztst(3)+rhc2*axst(3)
      else
        xh=fcct1*y(in1+4,nn1)+fcct2*y(in1+4,nn1-1)
      endif
      if(ispxx3.eq.2) then
        if(axst(4).ne.0) then
          xhe3=aztst(4)+rhc2*axst(4)
        else
          xhe3=fcct1*y(in1+5,nn1)+fcct2*y(in1+5,nn1-1)
        end if
      end if
c
    7 pl=17+log10(aztst(1)+rhc2*axst(1))
      fl=y(in1+1,nn1)
c
c  iterate to calculate log f
c
      yh=1-xh-zhc
      nosd=.true.
      notd=.true.
      nit=0
    8 call eqstf(fl,tl,xh,yh,zhc,nosd,notd)
      nit=nit+1
      dfl=(pl-log10(pt(1)))/pt(2)
      fl=fl+dfl
      if( abs(dfl).lt.1.e-7) go to 8100
c
      if(nit.lt.20) go to 8
      if(istdpr.gt.0) write(istdpr,150) pl,fl,tl,dfl
c
 8100 y(in1+1,nn1)=fl
      y(in1+2,nn1)=tl
      y(in1+4,nn1)=xh
      if(ispxx3.eq.2) y(in1+5,nn1)=xhe3
c
c  set CNO and 4He abundances, possibly again using extrapolation
c
      if(ispcno+isphec.gt.0) then
	do 8105 i = 1, ispcno+isphec
	if(axst(4+i).ne.0) then
	  y(in1+ispxx3+3+i,nn1)=aztst(4+i)+rhc2*axst(4+i)
        else
	  y(in1+ispxx3+3+i,nn1)= fcct1*y(in1+ispxx3+3+i,nn1)+
     *                           fcct2*y(in1+ispxx3+3+i,nn1-1)
	end if
 8105   continue
      end if
      if(istdpr.gt.0) then
        write(istdpr,*) 'ispxx3, ispcno, isphec', ispxx3, ispcno, isphec
        write(istdpr,*) 'aztst:',aztst
        write(istdpr,*) 'axst:',axst
        write(istdpr,*) 'azt:',azt
        write(istdpr,*) 'ax:',ax
        write(istdpr,*) 'Now solution at innermost meshpoint is'
        write(istdpr,*) (y(in1+i-1,nn1),i=1,nvar)
      end if
c
c  test for resetting innermost Y, to keep zhc unchanged.
c#ai# Note: the treatment of zhc (and zh, in general, in the present
c     routine, needs more thought
c
      if(iheccs.ne.0) then
	ycnew=1.d0-zhc-y(in1+4,nn1)
	if(abs(ycnew-y(in1+iyche4-1,nn1)).gt.1.e-10.and.istdpr.gt.0) 
     *    write(istdpr, '(/'' ***** Warning in s/r mshstr: ''/
     *       ''       Y(nn1) changed by'',
     *    1pe17.10,'' to '',e17.9,'' to keep Z fixed''/)')
     *    ycnew-y(in1+iyche4-1,nn1),ycnew
        y(in1+iyche4-1,nn1)=ycnew
      end if
      write(6,*) '#D# (4) xsimax =', xsimax
c
c  Reset for extrapolated solution. For the time being, only
c  reset radius and luminosity by mass scaling
c
      if(.not.time0) then
        y(in,nn1)= y(in,nn1) + delqlc/3
        y(in+3,nn1)=
     *    log10((10.d0**y(in+3,nn1)-alshft)*10.d0**delqlc+alshft)
      end if
c
c  Start section to reset mesh
c  ***************************
c
c  Possibly reset location of special points
c
    9 if(istrt3.ge.1.and.ntime.ge.1) then
          call spcpts(nspcpt, xspcpt, x, nn, qmxcor, qmxscn)
      end if
      write(6,*) '#D# (5) xsimax =', xsimax
c
c  set g(1-3), ranges and quantities that must be differentiated
c  numerically
c
      iz=ivarmx
      ifd=1
      iter=0
      jvar=4+knder
c 
      kointc=koint
      krpeat=0
c
c  start first loop through model
c
   10 do 30 n=1,nn1
c
c  recompute quantities for stretch
c
      idgrho = idgrhs
      call rhs(x(n),y(in1,n),zk,zrhs,dzdy,ap,aq,f,fd,alam,alamd,h,hd,
     .  iz,ifd,ifd,ifd,n,iter)
c
      idgrhs=idgrho
c
      call store(f,yw(1,n),3)
      call store(zrhs,z,5)
c
      z(1)=1.d1**z(1)
      yw(1,n)=z(1)*amm*yw(1,n)
      yw(4+knder,n)=z(4)
      yw(5+knder,n)=z(5)
      z(6)=1.d20*ft(1,1)
      yw(6+knder,n)=z(6)
      if(cvr(3,n).gt.-1.e-5) then
        z(7)=log10(1.e-5+cvr(3,n))
      else
        write(istdou,115) n, x(n), cvr(3,n)
        if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,115) 
     *    n, x(n), cvr(3,n)
        stop 'Stop 1 mshstr'
      end if
      if(ispxx3.eq.1) z(7)=1.e-3
      if(cvr(3,n).le.1.e-8) ntrsx3=n
      yw(7+knder,n)=z(7)
      dgrad=ddacad
      z(8)=dgrad
c
c  as temporary fudge, switch off term in dgrad beneath base of
c  convective envelope, when istrt2 .ge. 1 (to avoid sharp feature
c  at base of convection zone)
c
      if(istrt2.ge.1.and.n.ge.nl(1)-5) z(8)=0
      yw(8+knder,n)=z(8)
c
c  test for inclusion of 4He burning quantities
c
      if(theccs) then
	if(z(5).gt.2.e-10) then
	  z(9)=1-zh(n)
	  z(10)=0
        else
	  z(9)=zrhs(iyche4)
	  z(10)=1.d20*ft(1,iyche4-4)
        end if
	yw(9+knder,n)=z(9)
	yw(10+knder,n)=z(10)
c
      end if
c
      if(n.le.2) go to 22
c
c  second derivative of dgrad
c
c  two formulations have been used.
c  1) scale by factor 1.e-11. this was used for original model 1
c     calculation, up to about 1982.
c  2) scale by min(abs(x(n)),1.e-9). this was probably introduced
c     for computation of models at other masses, around 1982.
c
c  these may be selected by using iwdgrd = 1 or 2 respectively 
c
c  In addition introduce option iwdgrd = 3, to switch off second
c  derivative entirely.
c  
      n1=n-1
      n2=n-2
      dx2=x(n)-x(n2)
      if(dx2.eq.0) then
        write(istdou,116) n, x(n)
        if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,116) n, x(n)
        stop 'Stop 2 mshstr'
      else
        drx1=(x(n)-x(n1))/dx2
        drx3=(x(n1)-x(n2))/dx2
        if(drx1.eq.0.or.drx3.eq.0) then
          write(istdou,118) n, x(n)
          if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,118) n, x(n)
          stop 'Stop 3 in mshstr'
        end if
      end if
c
c  set weight for second derivative of dgrad
c
      if(iwdgrd.eq.1) then
        fcdgrd=1.e-11
      else if(iwdgrd.eq.2) then
	fcdgrd=min(abs(x(n)),1.d-9)
      else
	fcdgrd=0
      end if
c
      yw(20,n1)=0.5*fcdgrd*dgrad*
     *    (drx1*yw(13,n2)-yw(13,n1)+drx3*yw(13,n))/(dx2*dx2*drx1*drx3)
c
   22 if(n.eq.1.or.n.eq.nn1) yw(20,n)=0
      do 30 i=1,jvar
      yi=z(i)
      if(n.eq.1) then
        trmax(i)=yi
        trmin(i)=yi
      else
        trmax(i)=max(trmax(i),yi)
        trmin(i)=min(trmin(i),yi)
      end if
   30 continue
      write(6,*) '#D# (6) xsimax =', xsimax
c
c  differentiate numerically
c
      do 31 k=4,3+knder
   31 call derive(x,yw(k+knder,1),yw(k,1),nn1,iyw,iyw,1,1)
c
c  replace X derivative by Gaussian at discontinuities, set up 
c  combined discontinuity array
c
      nndisc=0
c
      call smngdl(x,yw(4+knder,1),yw(4,1),nn1,iyw,iyw)
      call resxdr(x,yw(5+knder,1),yw(5,1),nn1,iyw,iyw,deltax,xrxdrc,
     *  nrxdrc,nnrxdr,'X')
      call comdsc(ncdisc,xcdisc,nndisc,nrxdrc,xrxdrc,nnrxdr,deltax)
      call resxdr(x,yw(6+knder,1),yw(6,1),nn1,iyw,iyw,deltax,xrxdrc,
     *  nrxdrc,nnrxdr,'r_X')
      call comdsc(ncdisc,xcdisc,nndisc,nrxdrc,xrxdrc,nnrxdr,deltax)
c
c  to smooth possible corners in Xh-derivative, perform
c  running mean
c
      if(ismth.eq.1) then
        call rnmean(yw(5,1),yw(5+knder,1),nn1,iyw,iyw,nsmth)
        do 32 n=1,nn1
   32   yw(5,n)=yw(5+knder,n)
      end if
c
c  set reduction factor for effect of rX
c
      trxmax=0
      nrxmax=nn
      do 33 n=1,nn
      if(abs(yw(6,n)).gt.trxmax) then
	trxmax=abs(yw(6,n))
	nrxmax=n
      end if
   33 continue
c
c  with He burning, replace Y derivative by Gaussian at discontinuity
c
      if(theccs) then
        call resxdr(x,yw(9+knder,1),yw(9,1),nn1,iyw,iyw,deltax,xrxdrc,
     *    nrxdrc,nnrxdr,'Y')
        call comdsc(ncdisc,xcdisc,nndisc,nrxdrc,xrxdrc,nnrxdr,deltax)
        call resxdr(x,yw(10+knder,1),yw(10,1),nn1,iyw,iyw,deltax,xrxdrc,
     *    nrxdrc,nnrxdr,'r_Y')
        call comdsc(ncdisc,xcdisc,nndisc,nrxdrc,xrxdrc,nnrxdr,deltax)
c
c  to smooth possible corners in Yh-derivative, perform
c  running mean
c
	if(ismth.eq.1) then
          call rnmean(yw(9,1),yw(9+knder,1),nn1,iyw,iyw,nsmth)
          do 34 n=1,nn1
   34     yw(9,n)=yw(9+knder,n)
        end if
c
c  set reduction factor for effect of rY
c
        trxmax=0
        nrxmax=nn
        do 35 n=1,nn
        if(abs(yw(10,n)).gt.trxmax) then
	  trxmax=abs(yw(10,n))
	  nrxmax=n
        end if
   35   continue
c
        if(istdpr.gt.0) write(istdpr,*) 'nrxmax, trxmax, Y(nrxmax) =',
     *    nrxmax,trxmax,y(in11+iyche4,nrxmax)
        fry=
     *    frxfac*y(in11+iyche4,nrxmax)/(xmshcs+y(in11+iyche4,nrxmax))
c
c  reduce factor in case of very low reaction rates
c
	dyest=abs(dt*(trmax(10)-trmin(10))*1.d-20)
	fryest=dyest/(dyest+1.d-8)
	if(fryest.le.0.99d0.and.istdpr.gt.0) 
     *    write(istdpr,'(/a,1pe13.5)') ' rY term reduced by factor',
     *    fryest
	fry=fry*fryest
c
      end if
c
c  to avoid problems due to large step at innermost meshpoint,
c  after rc has been reset, replace the last derivatives by 
c  constant value
c
      nn12=nn1-1
      do 36 k=4,3+knder
      do 36 n=nn12,nn1
   36 yw(k,n)=yw(k,nn1-2)
c
c  set ranges
c
      do 38 j=1,jvar
   38 rg(j)=max(trmax(j)-trmin(j),1.d-30)
c
c  weights
c
      weight(1)=wr/rg(1)
      weight(2)=wpt/rg(2)
      weight(3)=wpt/rg(3)
      weight(4)=1.d0/rg(4)
      weight(5)=wxhc
      weight(6)=frx/rg(6)
      weight(7)=wx3/rg(7)
      weight(8)=wdgr/rg(8)
      weight(jvar)=weight(8)
      if(theccs) then
	weight(9)=wxhc
	weight(10)=fry/rg(10)
      end if
c
      tx=wx/(x(nn1)-x(1))
      tx=tx*tx
c
      if(iprtlc.gt.0.and.istdpr.gt.0) then
	write(istdpr,120) (i,rg(i),i=1,jvar)
	write(istdpr,122) (i,weight(i),i=1,jvar)
      end if
c
c  set lower limit of X3 term
c
      ntrsx3=ntrsx3+3
c
c  set stretching integrand
c
      rs=1.d11*10.d0**y(in1,1)
      xbconp=xbconv
      if(rcnvos.gt.0.5*rs.and.rcnvos.lt.0.99*rs) then
        xbconv=rcnvos/rs
      else
        xbconv=rczli/rs
      end if
c
c  test for extrapolating to position in next model
c  Certainly needs refinement.
c
      if(wmshos.lt.0) then
	xbconx=2*xbconv-xbconp
	if(abs(xbconx-xbconv).lt.drmsos) then
	  if(istdpr.gt.0) write(istdpr,
     *      '(/'' xbconp, xbconv, xbconx'',1p3e13.5/
     *      '' reset xbconv to xbconx''/)') xbconp, xbconv, xbconx
	  xbconv=xbconx
        end if
      end if
      if(istdpr.gt.0) write(istdpr,*) 
     *  'rs, rczli, rcnvos/rs,rczli/rs,xbconv =',
     *  rs, rczli, rcnvos/rs,rczli/rs,xbconv
c
c  test for setting up for convective core
c
      if(convc) then
c
c  The following two lines have given problems, because apparently
c  frmxc and nmxcor are not consistent with qmxcor. Needs fixing.
c  for now, replace. Although retain for testing.
c
	qlctst=frmxc*x(nmxcor)+(1-frmxc)*x(nmxcor-1)
	qqctst=10.**qlctst
	if(istdpr.gt.0) then
	  write(istdpr,*) 'In mshstr, nmxcor, frmxc =', nmxcor, frmxc
	  write(istdpr,*) '  qqctst, qmxcor =', qqctst, qmxcor
        end if
c
	qqc=qmxcor
	qlc=log10(qqc)
	if(nmxscn.gt.0) then
	  qlmxsc=log10(qmxscn)
        else
	  qlmxsc=-100.d0
	end if
	if(istdpr.gt.0) write(istdpr,*) 'qqc =',qqc
	qmxmax=max(qmxmax,qmxcor)
	if(qmxmin.gt.0) then
	  qmxmin=min(qmxmin,qmxcor)
        else
	  qmxmin=qmxcor
        end if
	if(istdpr.gt.0) write(istdpr,*) 
     *    'In mshstr, qmxmin, qmxmax =',qmxmin,qmxmax
	qlmxcr=log10(qmxcor)
      end if
c
c  set boundaries of convective regions
c
      if(inc.gt.0) then
        do 40 i=1,inc
        kf=nf(i)
        kl=nl(i)
        frf=frcf(i)
        frl=frcl(i)
c
        if(kf.eq.1) then
          qlf(i)=x(kf)
        else
          qlf(i)=frf*x(kf) + (1-frf)*x(kf-1)
        end if
c
        if(kl.eq.nn) then
          qll(i)=-60.
        else
c
c  corrected 15/3/06. Previous, erroneous, expression was
c..          qll(i)=frl*x(kl) + (1-frl)*x(kl-1)
c
          qll(i)=frl*x(kl) + (1-frl)*x(kl+1)
        end if
   40   continue
      end if
c
      if(iprtlc.gt.0.and.istdpr.gt.0) write(istdpr,112)
      write(6,*) '#D# (7) xsimax =', xsimax
c
      ismth=0
c
      if(icngr.ne.1.and.idgn75.gt.0) then
	if(iopn75.eq.0) then
          open(75,file='ttt.mshstr',status='unknown',form='unformatted')
	  if(istdpr.gt.0) write(istdpr,*) 
     *      'Opening ttt.mshstr for idgn75 =',idgn75
	end if
	ndgn=nvar+jvar+5
	write(75) ntime,ndgn,nn1,nvar
	if(istdpr.gt.0) write(istdpr,*) ' idgn75 =',idgn75,
     *    ' ntime,nvar+jvar+5,nn1,jvar,nvar',ntime,ndgn,nn1,jvar,nvar
      end if
c
      do 45 n=1,nn1
      xr=1.d1**(y(in1,n)-y(in1,1))
      qq=10.**x(n)
      sum=tx
c
c  add a term centred on overshoot radius
c  For wmshos .lt. 0, skew it towards radiative interior
c
      dxovs=(xr-xbconv)/drmsos
      dxovs=max(-10.d0,min(10.d0,dxovs))
      if(wmshos.gt.0) then
        trmovs=wmshos*exp(-dxovs*dxovs)
      else
	dxmlim=0.5
	wmshao=abs(wmshos)
	if(dxovs.le.dxmlin) then
	  dxmsq=dxovs*dxovs
        else if(dxovs.le.sqrt(10.d0*dxmlim)) then
	  dxmsq=dxovs**4/(dxmlim*dxmlim)
	else
	  dxmsq=100.d0
        end if
	trmovs=wmshao*exp(-dxmsq)
      end if
      sum=sum+trmovs
c
c  add a term at boundary of mixed region,
c  with some suppression in region where hydrogen abundance 
c  (or, in the case of 4He burning, the helium abundance)
c  is getting very low
c
c  Also reduce term for very small convective core (so far with
c  hard-coded limits)
c
      if(convc) then
	dxcvc=x(n)-qlmxcr
        dxcvc=min(10.d0,abs(dxcvc)/dqlmsc)
        trmcvc=((qmxcor+1.d-4)/(qmxcor+5.d-3))*wmshcc*exp(-dxcvc*dxcvc)
        if(((iheccs.eq.0.or.y(in1+2,n).le.7.7).and.y(in1+4,n).le.0.1)
     *    .and.istrt5.eq.0) then
          trmcvc=trmcvc*(0.01+10.*y(in1+4,n))
        else if(theccs.and.y(in11+iyche4,n).le.0.1.and.istrt5.eq.0) then
          trmcvc=trmcvc*(0.01+10.*y(in11+iyche4,n))
        end if
        sum=sum+trmcvc
      end if
c
      xh=abs(y(in11+5,n))
      yw(jvar,n)=yw(20,n)
      do 42 j=1,jvar
      ttrm=max(abs(yw(j,n)),1.d-36*rg(j))
      if(j.eq.7.and.n.le.ntrsx3) ttrm=0
      if((j.eq.5.or.j.eq.6).and.xh.le.1.e-2) then 
	redfct=(100.*xh)**2
	ttrm=redfct*ttrm
      end if
      if(theccs.and.(j.eq.9.or.j.eq.10)
     *  .and.y(in11+iyche4,n).le.1.e-2) then 
	redfct=(100.*y(in11+iyche4,n))**2
	ttrm=redfct*ttrm
      end if
      ttrm=weight(j)*ttrm
      ttrm=ttrm*ttrm
      trm(j)=ttrm
   42 sum=sum+ttrm
      yw(1,n)= sqrt(sum)
      if(iprtlc.gt.0.and.istdpr.gt.0) then
        if(mod(n-1,iprtlc).eq.0) 
     *    write(istdpr,100) n,x(n),sum,tx,(trm(j),j=1,jvar)
      end if
c
c  test for output of last case to file
c
      if(icngr.ne.1.and.idgn75.gt.0) then
	dgnout(1,n)=x(n)
	dgnout(2,n)=sum
	dgnout(3,n)=tx
	do j=1,jvar
	  dgnout(3+j,n)=trm(j)
        end do
	dgnout(4+jvar,n)=trmovs
	dgnout(5+jvar,n)=trmcvc
      end if
c
c  test for smoothing, if integrand varies too rapidly
c
      if(n.gt.1.and.ismth.eq.0.and.
     *  (yw(1,n).gt.2.*yw(1,n-1).or.yw(1,n-1).gt.2.*yw(1,n))) then
	ismth=-1
	if(istdpr.gt.0) write(istdpr,132) n, yw(1,n-1), yw(1,n)
      end if
   45 continue
c
c  make output
c
      if(icngr.ne.1.and.idgn75.gt.0) then
	jvar5=jvar+5
        write(75) ((dgnout(j,n),j=1,jvar5),(y(in11+j,n),j=1,nvar),
     *    n=1,nn1)
      end if
c
c  possibly smooth integrand
c
      if(ismth.eq.1) then
	if(istdpr.gt.0) write(istdpr,133)
        call rnmean(yw(1,1),yw(2,1),nn1,iyw,iyw,nsmth)
	do 47 n=1,nn1
   47   yw(1,n)=yw(2,n)
      end if
c
c  integrate to get ksi.
c
      yw(2,1)=0
c..      call vintk(x,yw,yw(2,1),1,nn1,iyw,iyw,kvint)
      if(istdpr.gt.0) write(istdpr,
     *    '(/'' Use trapezoidal rule for mshstr integration'')')
c
      do n=2,nn1
        yw(3,n)=0.5d0*(x(n)-x(n-1))*(yw(1,n-1)+yw(1,n))
      end do
c
c  finally, do integration
c
      do n=2,nn1
        yw(2,n)=yw(2,n-1)+yw(3,n)
      end do
c
c  entry point for start of interpolation
c
   56 fct= float(nn-1)/yw(2,nn1)
      ismth=0
      daks=0
      write(6,*) '#D# (8) xsimax =', xsimax
      do 57 n=1,nn1
      yw(1,n)=fct*yw(1,n)
      aks(n)=1+fct*yw(2,n)
c
      if(n.gt.1) then
        daksp=daks
	daks=aks(n)-aks(n-1)
      end if
c
c  set maximum value of change in aks
c
      if(n.eq.1) then
        daksmx=0
        naksmx=1
      else if(aks(n)-aks(n-1).gt.daksmx) then
        daksmx=aks(n)-aks(n-1)
        naksmx=n
      end if
c
c  test for excessively rapid change, and hence smoothing
c
      if(n.gt.2.and.ismth.eq.0.and.
     *  (daks.gt.2.*daksp.or.daksp.gt.2.*daks)) then
	ismth=-1
	if(istdpr.gt.0) write(istdpr,142) n, aks(n-2),aks(n-1),aks(n)
      end if
   57 continue
c
c  for i_dakslin = 1, find all intervals of excessive change in aks
c
      if(i_dakslin.ge.1) then
	i_daks=0
	n2=0
	do n=1,nn
	  if(aks(n)-aks(n-1).gt.dakslin) then
	    n1=max(1,n-2)
	    if(n1.gt.n2) then
	      i_daks=i_daks+1
	      aksms(1,i_daks)=aks(n1)
	    end if
	    n2=min(nn,n+1)
	    aksms(2,i_daks)=aks(n2)
          end if
        end do
	if(i_daks.gt.0.and.istdpr.gt.0) write(istdpr,
     *    '(/'' Intervals of excessive change in aks''/
     *       '' i, aks1, aks2:''/(i3,2f10.3))')
     *    (i,aksms(1,i),aksms(2,i),i=1,i_daks)
      end if
c
      if(icngr.ne.1.and.idgn75.gt.0) then
	write(75) (aks(n),n=1,nn1)
	if(idgn75.eq.1) then
	  close(75)
	  iopn75=0
        else
	  iopn75=1
        end if
      end if
c
c  possibly smooth integral
c
      if(ismth.eq.1) then
	if(istdpr.gt.0) write(istdpr,143)
        call rnmean(aks,yw(2,1),nn1,1,iyw,nsmth)
	do 59 n=1,nn1
   59   aks(n)=yw(2,n)
      end if
c
c  test for using linear interpolation and repeat pass
c
      if(daksmx.ge.xsimax) then
        irpeat=1
        kointc=1
        write(istdou,160) daksmx, naksmx, aks(naksmx)
        if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,160) 
     *    daksmx, naksmx, 
     *    aks(naksmx)
        if(iprtlc.eq.0.and.istdpr.gt.0) then 
	  write(istdpr,165) (aks(n),n=1,nn1)
	  iprtlc=1
	end if
        krpeat=krpeat+1
        if(krpeat.ge.3) then
	  if(ignrpt.ne.1) then
            write(istdou,162)
            if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,162)
            stop 'Stop 4 in mshstr'
          else
            write(istdou,163)
            if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,163)
	    irpeat=0
          end if
        end if
c
      else
        irpeat=0
      end if
c
c  set limits for using linear extrapolation on largest jump in
c  aks, if too big
c  NOTE: this should really be extended to all large intervals in aks
c
      if(daksmx.ge.3.) then
	aksm1=aks(naksmx-1)
	aksm2=aks(naksmx)
	if(istdpr.gt.0) write(istdpr,
     *    '(/'' Force linear interpolation for aks between'',2f10.3)')
     *    aksm1, aksm2
      else
	aksm1=1
	aksm2=-1
      end if
c
      
c
      if(iprtlc.ne.0.and.istdpr.gt.0) write(istdpr,165) (aks(n),n=1,nn1)
c
      dks=(aks(nn1)-aks(1))/(nn-1)
c  store old x and y in yw
      if(time0) then
        do 62 n=1,nn1
        call store(y(in1,n),yw(1,n),nvar)
   62   continue
        kint=nvar
      else
        do 64 n=1,nn1
        call store(y(in1,n),yw(1,n),nvar)
   64   call store(y(in,n),yw(nvar+1,n),nvar)
        kint=2*nvar
      end if
c
      ix=kint+2
      do 65 n=1,nn1
      yw(ix-1,n)=epsg(n)
   65 yw(ix,n)=x(n)
c
c  test for setting up for convective core
c
      if(convc) then
        aksc1=aks(nmxcor-2)
        akscc=aks(nmxcor)
        aksc2=aks(nmxcor+2)
      else
        nfc=-1
        aksc1=-1
        akscc=-1
        aksc2=-1
      end if
c
c  test for setting up for bottom of convective envelope 
c  (NEEDS refining, probably)
c
      if(istrt2.ge.1.and.idiffus.ge.1) then
        aksce1=aks(nl(1)-3)
        aksce2=aks(nl(1)+3)
      else
        aksce1=-1
        aksce2=-1
      end if
c
c  interpolate to uniform ksi-mesh
c
      init=1
c
c  initialize counter for discontinuities
c
      idisc=1
c
      if(istdpr.gt.0) write(istdpr,*) 'In mshstr, zhc =',zhc
c
      j_daks=0
      do 70 n=2,nn
      aksn=aks(1)+dks*(n-1)
      if(n.eq.nn) aksn=aks(nn1)
c
c  possibly reset interval for linear interpolation
c
      if(i_dakslin.ge.1.and.j_daks.le.i_daks) then
	if(j_daks.eq.0.or.
     *     (aksn.gt.aksms2.and.j_daks.lt.i_daks)) then
	  j_daks=j_daks+1
	  aksms1=aksms(1,j_daks)
	  aksms2=aksms(2,j_daks)
	end if
      else
	aksms1=1
	aksms2=-1
      end if
      
c
c  test for using linear interpolation near edge of convective core
c
      if(istrt1.eq.2.and.icngr.ne.1.and.
     *  aksn.ge.aksc1.and.aksn.lt.akscc) then
        init=1
        call lir1(aksn,aks,yint,yw,ix,iyw,nmxcor,init,inter)
        if(iprtlc.ne.0.and.istdpr.gt.0)
     *    write(istdpr,*) 'Linear interpolation at q =',10.d0**yint(ix)
c
c  test for using linear interpolation near bottom of convective envelope
c
      else if(istrt2.ge.1.and.icngr.ne.1.and.idiffus.ge.1.and.
     *  aksn.ge.aksce1.and.aksn.lt.aksce2) then
        init=1
        call lir1(aksn,aks,yint,yw,ix,iyw,nn1,init,inter)
        if(iprtlc.ne.0.and.istdpr.gt.0)
     *    write(istdpr,*) 'Linear interpolation at q =',10.d0**yint(ix)
      else if(istrt1.eq.2.and.icngr.ne.1.and.
     *  aksn.ge.akscc.and.aksn.le.aksc2) then
        init=1
        call lir1(aksn,aks(nmxcor),yint,yw(1,nmxcor),ix,iyw,
     *    nn1-nmxcor+1,init,inter)
        if(iprtlc.ne.0.and.istdpr.gt.0)
     *    write(istdpr,*) 'Linear interpolation at q =',10.d0**yint(ix)
      else if(aksn.ge.aksc1.and.aksn.le.aksc2.and.icngr.ne.1) then
        init=1
        call lir1(aksn,aks,yint,yw,ix,iyw,nn1,init,inter)
        if(iprtlc.ne.0.and.istdpr.gt.0)
     *    write(istdpr,*) 'Linear interpolation at q =',10.d0**yint(ix)
      else if(icngr.ne.1.and.((aksn.ge.aksm1.and.aksn.le.aksm2).or.
     *        (aksn.ge.aksms1.and.aksn.le.aksms2))) then
        init=1
        call lir1(aksn,aks,yint,yw,ix,iyw,nn1,init,inter)
        if(iprtlc.ne.0.and.istdpr.gt.0)
     *    write(istdpr,*) 'Linear interpolation at q =',10.d0**yint(ix)
      else
        call lirk(aksn,aks,yint,yw,ix,iyw,nn1,init,inter,kointc)
        init=2
	if(icngr.ne.1) then
	  qx=10.d0**yint(ix)
          if(abs(qx-qmxmax).lt.0.002.or.abs(qx-qmxmin).lt.0.002) then
	    init=1
          else if(idisc.le.nndisc.and.
     *      abs(yint(ix)-xcdisc(idisc)).lt.3*deltax) then
            init=1
          else if(idisc.le.nndisc.and.
     *      yint(ix).le.xcdisc(idisc)-3*deltax) then
	    idisc=idisc+1
	  end if
	  if(init.eq.1) then
            call lir1(aksn,aks,yint,yw,ix,iyw,nn1,init,inter)
            if(iprtlc.ne.0.and.istdpr.gt.0)
     *        write(istdpr,*) 'Linear interpolation at q =',
     *          10.d0**yint(ix)
	  end if
        end if
      end if
c
c  test for errors in interpolation
c
      if(inter.eq.-1) stop 'Stop 5 in mshstret'
c
      if(time0) then
        call store(yint,y(in1,n),nvar)
      else
        call store(yint,y(in1,n),nvar)
        call store(yint(nvar+1),y(in,n),nvar)
      end if
c
c  test for setting CNO and possibly 4He abundances in cvr
c
      if(istrt0.ge.1.and.istrt1.ge.1) then
	call setcno(y(in1,n),n)
        if(iheccs.ne.0) then
          cvr(icvhe4,n)=y(in11+iyche4,n)
          cvr(icvc12,n)=y(in11+iycc12,n)
        end if
      end if
      epsg(n)=yint(ix-1)
      x(n)=yint(ix)
   70 continue
c
c  as a crude fix, reset unreaonable values of composition
c
c #ai# Needs fixing in case of diffusion
c
      if(istdpr.gt.0) write(istdpr,*) 'In mshstr, xmxrhe =',xmxrhe
      if(idiffus.eq.0) then
c..	write(istdpr,'(/a/(i5,1p5e15.7))') 'n, aks, X, Y, Xp, Yp',
c..     *   (n,aks(n),yw(5,n),yw(iyche4,n),
c..     *    yw(nvar+5,n),yw(nvar+iyche4,n),n=500,nn1)
	init=1
	do 72 n=1,nn
	if(y(in1+4,n).lt.xhzlm1) then
	  if(iheccs.ne.0.and.x(n).gt.xmxrhe) then
            aksn=aks(1)+dks*(n-1)
            call lir1(aksn,aks,yint,yw(iyche4,1),1,iyw,
     *        nn1,init,inter)
	    yh1=y(in1+iyche4-1,n)+y(in1+4,n)-1.d-10
	    if(abs(yint(1)-yh1).gt.1.d-8.and.istdpr.gt.0) then
	      write(istdpr,*) 
     *          ' Error in curr.: interpolated, reset Y =',yint(1), yh1
	    end if
	    y(in1+iyche4-1,n)=yh1
	    if(y(in1+iyche4-1,n).lt.0.9.and.istdpr.gt.0) then
	      write(istdpr,*) 'Error in curr.: init, aksn, Y =',
     *         init, aksn, y(in1+iyche4-1,n)
	    end if
	    init=2
	  end if
	  y(in1+4,n)=1.d-10
	end if
        do 71 i=6,nvar
        if(y(in1+i-1,n).lt.1.e-20) y(in1+i-1,n)=1.d-10
   71   continue
   72   continue
      end if
c
c  possibly reset at new (extrapolated) timestep
c
      if(idiffus.eq.0.and..not.time0) then
	init=1
	do 75 n=1,nn
	if(y(in+4,n).lt.xhzlm1) then
	  if(iheccs.ne.0.and.x(n).gt.xmxrhe) then
            aksn=aks(1)+dks*(n-1)
            call lir1(aksn,aks,yint,yw(nvar+iyche4,1),1,iyw,
     *        nn1,init,inter)
	    yh1=y(in+iyche4-1,n)+y(in+4,n)-1.d-10
	    if(abs(yint(1)-yh1).gt.1.d-8.and.istdpr.gt.0) then
	      write(istdpr,*) ' Error in new: interpolated, reset Y =',
     *          yint(1), yh1
	    end if
	    y(in+iyche4-1,n)=yh1
	    if(y(in+iyche4-1,n).lt.0.9.and.istdpr.gt.0) then
	      write(istdpr,*) 'Error in new: init, aksn, Y =',
     *          init, aksn, y(in+iyche4-1,n)
	    end if
	    init=2
	  end if
	  y(in+4,n)=1.d-10
	end if
        do 74 i=6,nvar
        if(y(in+i-1,n).lt.1.e-20) y(in+i-1,n)=1.d-10
   74   continue
   75   continue
      end if
c
c  test for resetting luminosity at new timestep by linear interpolation 
c  where it is close to the shift value
c
      if(alshft.gt.0) then
	init=1
	do 77 n=1,nn
	if(y(in+3,n).le.log10(3*alshft)) then
	  aksn=aks(1)+dks*(n-1)
	  call lir1(aksn,aks,yint,yw(4+nvar,1),1,iyw,nn1,init,inter)
	  init=2
	  y(in+3,n)=yint(1)
        end if
   77   continue
      end if
c
c  test for non-monotonic mesh, and, if so, reset by
c  interpolation in x
c
      call testmono(x,y,in1,in,iy,nn,nn1,yw,iywstr,nvar,ix,time0)
c
c  for istrt2 .ge. 1, and with diffusion, reset near base of convective
c  envelope.
c  NEEDS further refinement.
c  Changed 2/9/03 to be applied in all cases
c
      if(icngr.ne.1) then 
	if(istrt2.ge.1) then
          call rescbz(x,y,in1,in,iy,nn,nn1,yw,iywstr,nvar,ix,qll,time0)
        else if(istrt3.ge.1.and.ntime.ge.1) then
c
c  use location of special points, set before mesh stretch with call
c  of spcpts
c
	  idiag=1
	  nczp=-1
	  do kspcpt=1,nspcpt
            ncz=-1
            xcz=xspcpt(kspcpt)
            do n=1,nn
              if(x(n).le.xcz.and.ncz.eq.-1) ncz=n
            end do
            frcz=(x(ncz-1)-xcz)/(x(ncz-1)-x(ncz))
            if(frcz.lt.0.5) ncz=ncz-1
	    write(istdpr,*) '#D# kspcpt, xcz, ncz, frcz',
     *        kspcpt, xcz, ncz, frcz
c
c  reset mesh if point is sufficiently removed from previous point
c  (note that it is assumed that xspcpt is stored in order of 
c  decreasing q).
c
	    if(nczp.eq.-1.or.ncz.gt.nczp+3) then
	      npts=2
	      call rczms1(x,y(in1,1),y(in,1),nn,iy,ncz,xcz,npts,1,
     *          istrt3,idiag)
	    end if
	    nczp=ncz
	  end do
        end if
      end if
c
c  test for resetting boundary of convective regions and 
c  mixed core
c
      if(inc.gt.0) then
        initmc=0
	initsn=0
	icf=1
	icl=1
        do 82 n=1,nn
c
c  test for convective region
c
        if(icf.le.inc.and.x(n).le.qlf(icf)) then
	  nf(icf)=n
	  if(n.eq.1) then
	    frcf(icf)=1
          else
	    if(istdpr.gt.0) then
	      write(istdpr,*) 
     *          'Setting frcf at n =',n, ' for icf =',icf
	      write(istdpr,*) 'qlf(icf), x(n-1),x(n) =',
     *          qlf(icf), x(n-1),x(n) 
	    end if
	    frcf(icf)=(qlf(icf)-x(n-1))/(x(n)-x(n-1))
	  end if
	  icf=icf+1
        end if
c
        if(icl.le.inc.and.x(n).lt.qll(icl)) then
	  nl(icl)=n-1
          if(istdpr.gt.0) then
	    write(istdpr,*) 'Setting frcl at n =',n, ' for icl =',icl
	    write(istdpr,*) 'qll(icl), x(n-1),x(n) =',
     *        qll(icl), x(n-1),x(n) 
	  end if
	  frcl(icl)=(qll(icl)-x(n))/(x(n-1)-x(n))
	  icl=icl+1
        end if
c
c  test for mixed core
c

        if(convc) then
          if(x(n).lt.qlc) then
c
c  test for initialization
c
            if(initmc.eq.0) then
              initmc=1
              nmxcor=n-1
              if(istdpr.gt.0) then
	        write(istdpr,*) 'Setting frmxc at n =',n
	        write(istdpr,*) 'qlc, x(n-2),x(n-1) =',
     *            qlc, x(n-2),x(n-1) 
	      end if
              frmxc=(qlc-x(n-2))/(x(n-1)-x(n-2))
            end if
c
c  test for resetting nmxscn
c
	    if(x(n).lt.qlmxsc.and.nmxscn.gt.0.and.initsn.eq.0) then
	      nmxscn=n-1
	      initsn=1
	      if(istdpr.gt.0) write(istdpr,'(/
     *          '' Reset nmxscn to '',i5/)') nmxscn
	    end if
c
c  reset composition to central value in fully mixed region
c  (possibly excluding semiconvective region)
c  Treatment with diffusion deserves further thought, particularly
c  when core mixing is handled diffusively.
c
	    if(nmxscn.le.0.or.x(n).lt.qlmxsc) then
	      if(idiffus.le.0) then
	        kmax=nvar
              else
	        kmax=4+iccomp+idcomp
              end if
c
              do 81 k=5,kmax
   81         y(in11+k,n)=y(in11+k,nn)
            end if
          end if
	end if
   82   continue
c
	if(qll(inc).le.-60) then
	  nl(inc)=nn
	  frcl(inc)=1
	end if
        if(istdpr.gt.0) then
	  write(istdpr,'(/'' i, nf, frcf, nl, frcl reset to:''//
     *     (i4,i6,f10.6,i6,f10.5))') 
     *     (i,nf(i),frcf(i),nl(i),frcl(i),i=1,inc)
          if(convc) write(istdpr,'(//a,i6,f10.5)') 
     *      'nmxcor, frmxc reset to ',nmxcor, frmxc
        end if
      end if
c
      if(iprtlc.gt.0.and.istdpr.gt.0) then
        write(istdpr,128)
        rs=1.d1**y(in1,1)
        do 90 n=1,nn,iprtlc
        rr=1.d1**y(in1,n)/rs
   90   write(istdpr,130) n,x(n),rr,(y(in11+j,n),j=1,nvar)
c
      end if
c
c  in case of 4He burning, reset zh(n) for consistency
c
      if(iheccs.ne.0) then
	do 95 n=1,nn
	if(y(in1+4,n).lt.xhzlm1.and.x(n).le.xmxrhe) then
	  zh(n)=1.d0-y(in11+5,n)-y(in11+iyche4,n)
        end if
   95   continue
      end if
c
c  test for repeating
c
      if(irpeat.eq.1) then
        if(koint.eq.1) then
          if(istdpr.gt.0) write(istdpr,172)
        else
          if(istdpr.gt.0) write(istdpr,174)
          kointc=koint
c
c  reset nn1, so that mesh is set over entire (new) interval
c
          nn1 = nn
c
          go to 10
        end if
      end if
c
      ifdhe3=ifdho
      if(icngr.eq.1.and.istdpr.gt.0) then
	write(istdpr,*) '#D# Reset mesh. n, x, X'
	write(istdpr,'(i5,1p2e13.5)') (n, x(n), y(in1+4,n),n=nn-100,nn)
      end if
      return
  100 format(i4,1p21e11.3)
  101 format(//' output from s/r mshstr on file ttt.mshstr.n.out')
  102 format(///' in mshstr koint=',i2,
     .  ' and kvint =',i2//)
  103 format(//' set m0, m2. rhoc =',1pe13.5,'  m0, m2 =',2e13.5)
  104 format(//' ***** error in s/r mshstr. al0 = 0.'/
     *         '       term in al0 in setting of rc is suppressed')
  105 format(//' terms in sum to determine rc:',1p6e13.4)
  106 format(//' Warning in s/r mshstr:'/
     *         ' new rc =',1pe13.5,
     *  ' exceeds old value =',e13.5,' by more than 10 per cent'/
     *  ' Use rc = 1.1*(old value)')
  107 format(//' *** Warning in mshstr: new rhc restricted to ',
     *  '0.9 times previous value')
  109 format(//' new value of r     at innermost meshpoint:',1pe13.5/
     *         ' new value of log q at innermost meshpoint:',1pe13.5/
     *         ' nn1=',i4)
  110 format(//' ***** error in s/r mshstr.',
     *  ' reset log q at centre nonmonotonic'/
     *         '       use old value.')
  111 format(//' *** Warning in s/r mshstr. al0 = ',1pe13.5,
     *         ' Xc =',e13.5/
     *         '       y(4,nn1) reset by mass scaling')
  112 format(//' terms in stretching integrand'/
     *  ' n,x(n),sum,tx,(trm(j),j=1,jvar)')
  115 format(//' ***** error in s/r mshstr. At n =',
     *  i5,'   log q =',1pe13.5,'  X3 =',e13.5)
  116 format(//' ***** error in s/r mshstr. dx2 = 0 at n =',
     *  i5,'   log q =',1pe13.5)
  118 format(//' ***** error in s/r mshstr. drx1 or drx3 = 0 at n =',
     *  i5,'   log q =',1pe13.5)
  120 format(//(' rg(',i2,') =',1pe13.5))
  122 format(//(' weight(',i2,') =',1pe13.5))
  125 format('# terms in stretching integrand.'/'#'/
     *  '# Timestep no. ',i4/'#'/
     *  '#   n       q       r/R     comb. int.  cst. term    ',
     *  'ovs.    mix. core     r         log p      log T      ',
     *  'log L        Xh        rX    log(1e-5+X3)   ddacad   ',
     *  'd2(ddacad)'/'#')
  128 format(//' solution on new mesh'///' n,x,r/rs,y(1-6):'/)
  130 format(i4,1p15e13.5)
  132 format(//' **** Warning in s/r mshstr.',
     *  ' Excessively rapid change in integrand at n =', i5/
     *  '      yw(1,n-1) =',1pe13.5,'  yw(1,n) =',1pe13.5)
  133 format(/'      Apply running-mean smoothing to integrand')
  135 format(i5,1pe12.4,0pf10.7,1p15e11.3)
  142 format(//' **** Warning in s/r mshstr.',
     *  ' Excessively rapid change in integral at n =', i5/
     *  '      aks(n-2), aks(n-1), aks(n) =',3f12.5)
  143 format(/'      Apply running-mean smoothing to integral')
  150 format(//' **** iteration for fl not converged in s/r mshstr'//
     *  ' pl, tl, fl, dfl =',3f12.7,1pe13.5)
  160 format(/' ***** error in mshstr. Max. change in xsi =',f8.3,
     *  ' at n =',i6,'  is too large.'/
     *        '       This occurs at xsi =',f10.3/
     *        '       Use linear interpolation initially and repeat')
  162 format(//'      Excessive number of repeats. Give up.')
  163 format(//
     *  '      Excessive number of repeats. Continue nonetheless.')
  165 format(//' aks:'/(10f13.5))
  172 format(/' No need to repeat stretch, since koint = 1')
  174 format(/' Go back and repeat stretch on new mesh')
      end
      subroutine resxdr(x,y,dydx,nn,iy,idy,deltax,xc,nc,nnc,label)
c
c  replaces derivative in dydx by Gaussian of width deltax
c  in the vicinity of discontinuities
c
c  xc(i) and nc(i), i = 1, ..., nnc, return central value of x and
c  mesh number of discontinuities. If no discontinuities are found,
c  nnc is returned as zero.
c
c  Original version: 5/8/92
c
c  Modified 20/5/93, scaling derivative in test by x
c  (to avoid erroneously detecting discontinuities near the surface)
c
c  Modified 19/11/98, ensuring that discontinuity is spread over
c  at least 5 points
c
c  Modified 29/8/02, to allow for several discontinuities and modifying
c  condition for discontinuity.
c
c  Modified 21/8/03, setting combined old and new criterion for 
c  discontinuity, and modifying the resetting of the derivatives
c  to ensure a more continuous and consistent results. This is
c  flagged by the (hardcoded) parameter setting inewcr = 2
c
      implicit double precision (a-h, o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
      character*(*) label
      logical ldisc
      dimension x(*), y(iy,*), dydx(idy,*), xc(*), nc(*)
      dimension dxx(nnmax), dydxnw(nnmax), sdydx(nnmax), 
     *  dydxng(nnmax), sdydxg(nnmax)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      fxp(xx)=exp(max(xx,-40.d0))
c
c  hardcode flag for new discontinuity criterion
c  inewcr = 0 for old criterion
c  As a hack, apply both criteria when inewcr = 2 (20/8/03)
c
      inewcr=2
c
c  find range of function, for testing for discontinuity
c
      ymin=y(1,1)
      ymax=y(1,1)
      do 5 n=2,nn
      ymin=min(ymin,y(1,n))
      ymax=max(ymax,y(1,n))
    5 dxx(n)=abs(x(n)*(y(1,n)-y(1,n-1))/(x(n)-x(n-1)))
c
c  step through, to locate discontinuities
c  (crudely hardcoding criterion for discontinuity)
c
      n1=0
      n2=0
      discrt=100*(ymax-ymin)
      ic=0
      do 50 n=3,nn-2
      if(inewcr.eq.0) then
        ldisc=dxx(n).gt.discrt
      else if(inewcr.eq.1) then
	ldisc=dxx(n)/(dxx(n-2)+dxx(n+2)+ymax-ymin).gt.1.d0
      else
	ldisc=dxx(n)/(dxx(n-2)+dxx(n+2)+ymax-ymin).gt.1.d0.or.
     *        dxx(n).gt.discrt
      end if
c
      if(ldisc) then
	if(n1.eq.0) n1=max0(n-3,1)
	n2=min0(n+3,nn)
      else if(n1.gt.0) then
c
c  end of present discontinuity. Test for resetting derivative on
c  interval (n1,n2)
c
        if(abs(x(n2)-x(n1)).lt.500*deltax) then
c
c  as a desperate fudge, try to spread discontinuity over at least
c  7 points (but the whole thing is a bloody mess!)
c
          if(n2-n1.lt.7) then 
	    nadd=max0(1,(7-n2+n1)/2)
	    n1=n1-nadd
	    n2=n2+nadd
          end if
c
          if(istdpr.gt.0) write(istdpr,110) 
     *      label,n1,x(n1), y(1,n1), n2, x(n2), y(1,n2)
c
c  set gaussian centred on this interval
c
	  delta1=max(deltax,0.1d0*abs(x(n1)-x(n2)))
	  ic=ic+1
c..          xxc=0.5*(x(n2)+x(n1))
          xxc=x((n2+n1)/2)
	  xc(ic)=xxc
	  nc(ic)=0.5*(n1+n2)
          ag=(y(1,n1)-y(1,n2))/(1.77245*delta1)
          if(istdpr.gt.0) write(istdpr,*) 'xc, ag, deltax =',
     *      xxc, ag, delta1
c
          nr1=max0(1,n1-10)
          nr2=min0(nn,n2+10)
          if(istdpr.gt.0) write(istdpr,112)
	  ns1=-1
c
	  ns2=-1
          do 25 j=1,nn
          ddydx=dydx(1,j)
          if(j.gt.n1.and.j.lt.n2) then
            fct=(x(j)-x(n1))/(x(n2)-x(n1))
            ddydx=(1-fct)*dydx(1,n1)+fct*dydx(1,n2)
          end if
          xx=(x(j)-xxc)/delta1
	  dydxng(j)=fxp(-xx*xx)
	  if(inewcr.ne.2) then 
            dydxnw(j)=ddydx+ag*dydxng(j)
          else
            dydxnw(j)=ddydx
	  end if
          if(j.eq.1) then
            sdydx(j)=y(1,1)
            sdydxg(j)=0.
          else
            sdydx(j)=sdydx(j-1)
     *        +0.5d0*(dydxnw(j)+dydxnw(j-1))*(x(j)-x(j-1))
            sdydxg(j)=sdydxg(j-1)
     *        +0.5d0*(dydxng(j)+dydxng(j-1))*(x(j)-x(j-1))
          end if
          if(j.ge.nr1.and.j.le.nr2.and.istdpr.gt.0)
     *      write(istdpr,120) j, x(j), y(1,j),dydx(1,j),dydxnw(j),
     *         sdydx(j)
c
c  set ns1 and ns2 for rescaling with inewcr = 2
c
          if(ns1.eq.-1.and.x(j).le.xxc+3*delta1) ns1=j-1
          if(ns2.eq.-1.and.x(j).le.xxc-3*delta1) ns2=j
c
   25     continue
c
          ns1=max0(1,ns1)
	  ns2=min0(nn,ns2)
c
c  set scaling factor such that integral of new derivative matches
c  jump in y(1,.)
c
c  For inewcr = 2, this is applied at points +- 3*delta1 (set above)
c  from centre of discontinuity
c
          if(inewcr.ne.2) then
	    ns1=n1+1
	    ns2=n2-1
c
	    if(abs(sdydx(ns2)-sdydx(ns1)).gt.1.d-8*(ymax-ymin)) then
              scfct=(y(1,nn)-y(1,1)
     *            -sdydx(nn)+sdydx(ns2)-sdydx(ns1)+sdydx(1))/
     *            (sdydx(ns2)-sdydx(ns1))
            else
	      scfct=1.d0
            end if
c
            if(istdpr.gt.0) write(istdpr,114) scfct
c
	  else
	    if(istdpr.gt.0) write(istdpr,*) 
     *        'Rescale over range in x',x(ns1), x(ns2)
	    agn=(y(1,ns2)-y(1,ns1)-sdydx(ns2)+sdydx(ns1))/
     *         (sdydxg(ns2)-sdydxg(ns1))
            if(istdpr.gt.0) write(istdpr,*) 'ag reset from ',ag,
     *        '  to', agn
	    ag=agn
          end if
c
c  do final resetting
c
          if(istdpr.gt.0) write(istdpr,115)
          do 30 j=1,nn
          if(j.lt.ns1.or.j.gt.ns2) then
            ddydx=dydxnw(j)
          else if(inewcr.ne.2) then
            ddydx=scfct*dydxnw(j)
	  else
	    ddydx=dydxnw(j)+ag*dydxng(j)
          end if
          if(j.ge.nr1.and.j.le.nr2.and.istdpr.gt.0)
     *      write(istdpr,120) j, x(j), y(1,j),dydx(1,j),ddydx
   30     dydx(1,j)=ddydx
        end if
        n1=0
        n2=0
      end if
   50 continue
      nnc=ic
      return
  110 format(//' In s/r resxdr, discontinuity found. in ',a/
     *  ' n, x, y:'/
     *  (i5,1p2e13.5,i5,1p2e13.5,'   #resxdr'))
  112 format(/' n, x, y, old, new dydx, int(dydx):'/)
  114 format(/' derivative rescaled by factor ',1pe13.5)
  115 format(/' n, x, y, old, new dydx:'/)
  120 format(i5,1p5e13.5)
      end
      subroutine testmono(x,y,in1,in,iy,nn,nn1,yw,iywst,nvar,ix,time0)
c
c  tests monotonicity of mesh set up in s/r mshstr and
c  resets model in case of problems
c
c  Original version: 24/9/98
c
c  Corrected (and finally brought to have an effect) 10/10/02
c
c  Modified 11/6/03, requiring two decreasing points to mark
c  end of non-monotonous region
c
      implicit double precision (a-h, o-z)
      logical time0
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
      parameter(iyint = 2*ivarmx+1)
c
      dimension x(*), y(iy, *), yw(iywst,*)
      dimension yint(iyint)
      common/noiter/ iter, ntime
      common/sooner/ xw(1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      if(istdpr.gt.0) then 
        idiag=1
      else
        idiag=0
      end if
c
      do 10 n=1,nn1
   10 xw(n)=yw(ix,n)
c  
      kint=ix-1
c
c  scan through model to test for non-monotonicity
c
      xp=x(1)
      np=1
      n=1
      n1=0
c
   20 n=n+1
      if(x(n).ge.x(n-1).and.n1.eq.0) then
c
c  non-monotonic region starts
c
        if(idiag.gt.0) write(istdpr,105) n,x(n-1),x(n)
        n1=n-1
      end if
c
c  test whether this is the end of a non-monotonic region
c  Note: redefined 11/6/03 to require two consecutive points
c  going down
c
      if(n1.gt.0.and.(n.eq.nn.or.
     *  (x(n).lt.x(n-1).and.x(n+1).lt.x(n)))) then
c
c  step back to find previous monotonuous point
c
   25   n1=n1-1
        if(x(n1).le.x(n)) go to 25
c
c  non-monotonous region localized. Reset model
c
        n1=max(1,n1-1)
        n2=min(nn,n+1)
        dx=(x(n2) - x(n1))/(n2-n1)
        if(idiag.gt.0) write(istdpr,110)
        do 30 nr=n1,n2
        xr=x(n1)+dx*(nr-n1)
        call lir1(xr,xw,yint,yw,kint,iywst,nn1,nr-n1+1,inter)
        if(idiag.gt.0) write(istdpr,120) 
     *    nr, x(nr), xr, y(in1+2,nr),yint(3)
        if(time0) then
          call store(yint,y(in1,nr),nvar)
        else
          call store(yint,y(in1,nr),nvar)
          call store(yint(nvar+1),y(in,nr),nvar)
        end if
   30   x(nr)=xr
c
        n1=0
        n=n2
c
      end if
c
      if(n.lt.nn) go to 20
c
      return
  105 format(/' Start of non-monotonic region in s/r testmono at n = ',
     *        i5/
     *        ' x(n), x(n-1) = ',1p2e17.9)
  110 format(//' Resetting mesh because of non-monotonicity.'/
     *         ' n, old x, new x, old, new log(T):'/)
  120 format(i5,1p2e15.7,0p2f12.7)
      end
      subroutine rescbz(x,y,in1,in,iy,nn,nn1,yw,iywst,nvar,ix,qll,time0)
c
c  Reset mesh near base of convective envelope, through linear 
c  interpolation in x. In particular, ensure point at base of
c  convective envelope.
c
c  Original version: 12/5/03
c
c
      implicit double precision (a-h, o-z)
      logical time0
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
      parameter(iyint = 2*ivarmx+1)
c
      dimension x(*), y(iy, *), yw(iywst,*),qll(*)
      dimension yint(iyint)
      common/noiter/ iter, ntime
      common/sooner/ xw(1)
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,jnc,frcf(6),frcl(6)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      if(istdpr.gt.0) then 
        idiag=1
      else
        idiag=0
      end if
c
c  find new location of convection-zone base
c
      ncbz=0
      dxcbz=1.d10
      do n=1,nn1
	if(abs(x(n)-qll(1)).lt.dxcbz) then
	  ncbz=n
	  dxcbz=abs(x(n)-qll(1))
	end if
      end do
c
      do 10 n=1,nn1
   10 xw(n)=yw(ix,n)
c
      kint=ix-1
c
      n1=ncbz-4
      n2=ncbz+4
      dx1=(qll(1) - x(n1))/(ncbz-n1)
      dx2=(x(n2) - qll(1))/(n2-ncbz)
      if(idiag.gt.0) write(istdpr,110)
      do 30 n=n1,n2
      if(n.le.ncbz) then
        xr=x(n1)+dx1*(n-n1)
      else
        xr=x(ncbz)+dx2*(n-ncbz)
      end if
      call lir1(xr,xw,yint,yw,kint,iywst,nn1,n-n1+1,inter)
      if(idiag.gt.0) write(istdpr,120) n, x(n), xr, y(in1+2,n),yint(3)
      if(time0) then
        call store(yint,y(in1,n),nvar)
      else
        call store(yint,y(in1,n),nvar)
        call store(yint(nvar+1),y(in,n),nvar)
      end if
   30 x(n)=xr
c
      return
  110 format(//' Resetting mesh near base of convective envelope'/
     *         ' n, old x, new x, old, new log(T):'/)
  120 format(i5,1p2e15.7,0p2f12.7)
      end
      subroutine smngdl(x,al,dl,nn,ial,idl)
c
c  smooth negative luminosity derivative with parabolic weight,
c  resetting integral to be unchanged
c
c  Original version: 18/10/02
c
      implicit double precision(a-h,o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
      dimension x(*),al(ial,*),dl(idl,*)
      dimension dl1(nnmax),w(101)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      idiag=1
      if(istdpr.le.0) idiag=0
c
      nsmth=20
      ntrunc=7
      nsmth2=2*nsmth+1
c
      sum=0
      do 10 k=1,nsmth2
      w(k)=(k-1)*(nsmth2-k)
   10 sum=sum+w(k)
c
      do 12 k=1,nsmth2
   12 w(k)=w(k)/sum
c
      do 14 n=1,nn
   14 dl1(n)=dl(1,n)
c
c  step through, looking for negative values
c
      n1=0
      do 50 n=1,nn
      if(dl(1,n).lt.-0.1d0.and.n1.eq.0) then
	n1=max(1,n-1)
      else if(dl(1,n).ge.0.1) then
c 
c  now end of negative region has been found
c
c  test for resetting
c
	if(n1.gt.0.and.n.ge.n1+nsmth) then
c
c  Reset derivative
c
          ns1=max(1,n1)
          ns2=min(nn,n)
          if(idiag.gt.0) write(istdpr,'(/a,2i5)')
     *      ' In s/r smngdl, reset derivative over range ',ns1, ns2
          sdl=0
          do l=ns1,ns2
            dl1(l)=0
            l1=l-nsmth-1
            do k=1,nsmth2
              l1=l1+1
	      if(l1.lt.ns1.or.l1.gt.ns2) then
		ww=0
              else if(l1.le.ns1+ntrunc.or.l1.ge.ns2-ntrunc) then
		ww=ww*min(((l1-ns1)/ntrunc)**2,((l1-ns2)/ntrunc)**2)
	      else
		ww=1.d0
	      end if
              dl1(l)=dl1(l)+ww*w(k)*dl(1,l1)
            end do
            if(idiag.gt.0) write(istdpr,*) 
     *        'l, dl, dl1, al', l, dl(1,l),dl1(l), al(1,l)
            if(l.gt.ns1) sdl=sdl+0.5d0*(dl1(l)+dl1(l-1))*(x(l)-x(l-1))
          end do
c
c  rescale to preserve total change over this region
c
          sfact=(al(1,ns2)-al(1,ns1))/sdl
          if(idiag.gt.0) write(istdpr,'(/a,1pe13.5)')
     *      ' Rescale reset derivative by factor',sfact
          do l=ns1,ns2
            dl1(l)=sfact*dl1(l)
          end do
c
	end if
        n1=0
      end if
c
   50 continue
c
c  reset derivative
c  
      do n=1,nn
        dl(1,n)=dl1(n)
      end do
      return
      end
      subroutine comdsc(ncdisc,xcdisc,nndisc,npdisc,xpdisc,nnpdsc,
     *   deltax)
c
c  Adds discontinuities in npdisc(i), xpdisc(i), i = 1,nnpdsc,
c  to combined set of discontinuities in ncdisc, xcdisc.
c  A discontinuity is only added if it is located more than 3*deltax
c  from any other discontinuity
c
c  NOTE: for now it is assumed that xpdisc and xcdisc decrease with
c  increasing i
c
c  Original version: 29/8/02
c
      implicit double precision(a-h, o-z)
      dimension ncdisc(*), xcdisc(*), npdisc(*), xpdisc(*)
      dimension is(2,100), ncnew(100), xcnew(100)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      idiag=1
      if(istdpr.le.0) idiag=0
c
c  test for any new discontinuities
c
      if(nnpdsc.le.0) return
c
c  test for previous discontinuities
c
      if(nndisc.le.0) then
        nndisc=nnpdsc
        do 10 n=1,nnpdsc
        ncdisc(n)=npdisc(n)
   10   xcdisc(n)=xpdisc(n)
c
        return
c
      else
c
c  go through and merge the two sets
c
        ic=1
        ip=1
        xc=xcdisc(1)
        xp=xpdisc(1)
        ii=1
c
   15   continue
        if(ip.le.nnpdsc.and.abs(xp-xc).lt.3*deltax) then
          xp=xpdisc(min(nnpdsc,ip+1))
          ip=ip+1
        else if((ip.le.nnpdsc.and.xp.gt.xc).or.ic.gt.nndisc) then
          is(1,ii)=2
          is(2,ii)=ip
          xp=xpdisc(min(nnpdsc,ip+1))
          ip=ip+1
          ii=ii+1
        else
          is(1,ii)=1
          is(2,ii)=ic
          xc=xcdisc(min(nndisc,ic+1))
          ic=ic+1
          ii=ii+1
        end if
        if(ip.le.nnpdsc.or.ic.le.nndisc) go to 15
c
c  now sorting array is set. Store in intermediate arrays
c
        ii=ii-1
        if(idiag.gt.0) write(istdpr,110)
        do 20 i=1,ii
        if(is(1,i).eq.1) then
          i1=is(2,i)
          xcnew(i)=xcdisc(i1)
          ncnew(i)=ncdisc(i1)
          if(idiag.gt.0) write(istdpr,115) i, ncnew(i), xcnew(i), 'old'
        else
          i1=is(2,i)
          xcnew(i)=xpdisc(i1)
          ncnew(i)=npdisc(i1)
          if(idiag.gt.0) write(istdpr,115) i, ncnew(i), xcnew(i), 'new'
        end if
   20   continue
c
c  restore in xcdisc, ncdisc
c
        do 25 i=1,ii
        xcdisc(i)=xcnew(i)
   25   ncdisc(i)=ncnew(i)
c
        nndisc=ii
c
      end if
      return
  110 format(//' In s/r comdsc. i, ndisc, xdisc:'/)
  115 format(i3,i5,1pe13.5,3x,a)
      end 
      subroutine spcpts(nspcpt, xspcpt, x, nn, qmxcor, qmxscn)
c
c  resets information in xspcpt(k), k = 1, ..., nspcpt, on
c  points where a meshpoint should be placed. For now this
c  is the location of the present and possibly past convective
c  cores.
c  #AI# Desciption to be updated
c  The points should be placed in order of decreasing q.
c  qmxcor (the mass of the convective core) is returned as
c  innermost point.
c
      implicit double precision (a-h, o-z)
      dimension xspcpt(*), x(*)
      dimension xf(6), xl(6)
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,jnc,frcf(6),frcl(6)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      data qbenv_min, qmxcor_max, qmxscn_max, init 
     *       /1.d0, 0.d0, 0.d0, 0/
      save
c
c  Initialize nspcpt to 0
c
      nspcpt=0
c
c  Set log(mass) at edges of convection zones
c
      do i=1,inc
	nnf=nf(i)
	xf(i)=frcf(i)*x(nnf)+(1.d0-frcf(i))*x(nnf-1)
	nnl=nl(i)
	xl(i)=frcl(i)*x(nnl)+(1.d0-frcl(i))*x(nnl+1)
	write(istdpr,*) '#D# i, frcl(i), nnl, x(nnl), x(nnl+1)',
     *    i, frcl(i), nnl, x(nnl), x(nnl+1)
      end do
      write(istdpr,*) '#D# nl:',(nl(i),i=1,inc)
      write(istdpr,*) '#D# frcl:',(frcl(i),i=1,inc)
      write(istdpr,*) '#D# xl:',(xl(i),i=1,inc)
c
c  test for fully convective star
c
      if(inc.eq.1.and.nl(inc).eq.nn) then
	if(istdpr.gt.0) write(istdpr,'(//
     *    '' Fully convective star in s/r spcpts.''/
     *    '' No setting of special points.''/)')
	return
      end if
c
c  set relevant bottom of convective envelope
c  (tests need refinement, most likely, and consistency with other
c  parts of the code)
c
      if(nl(inc).eq.nn) then
	ibenv=0
	do i=inc-1,1,-1
	  if(ibenv.eq.0.and.nl(i).lt.nf(i+1)-3) ibenv=i
	end do
	xbenv=xl(ibenv)
      else
	xbenv=xl(inc)
      end if
      write(istdpr,*) '#D# in spcpts, xbenv, xmxcor, xmxscn =',
     *  xbenv, log10(qmxcor), log10(qmxscn)
      write(istdpr,*) '#D# xbenv_min, xmxcor_max, xmxscn_max =',
     *  xbenv_min, log10(qmxcor_max), log10(qmxscn_max)
c
c  test for new extreme
c
      nspcpt=1
      xspcpt(nspcpt)=xbenv
      write(istdpr,*) '#D1# set nspcpt, log(q(nspcpt)) to',
     *  nspcpt, xspcpt(nspcpt)
      if(xbenv.gt.xbenv_min) then
	nspcpt=2
	xspcpt(nspcpt)=xbenv_min
        write(istdpr,*) '#D2# set nspcpt, log(q(nspcpt)) to',
     *    nspcpt, xspcpt(nspcpt)
      else
	xbenv_min=xbenv
      end if
c
c  test for edge of convective core
c
      if(qmxcor_max.gt.0.or.(qmxcor.gt.0.and.init.eq.0)) then
	init=1
	if(qmxcor.ge.qmxcor_max) then
	  nspcpt=nspcpt+1
	  xspcpt(nspcpt)=log10(qmxcor)
          write(istdpr,*) '#D3# set nspcpt, log(q(nspcpt)) to',
     *      nspcpt, xspcpt(nspcpt)
	  qmxcor_max=qmxcor
	else
	  nspcpt=nspcpt+1
	  xspcpt(nspcpt)=log10(qmxcor_max)
          write(istdpr,*) '#D4# set nspcpt, log(q(nspcpt)) to',
     *      nspcpt, xspcpt(nspcpt)
	  if(qmxcor.gt.0) then
	    nspcpt=nspcpt+1
	    xspcpt(nspcpt)=log10(qmxcor)
            write(istdpr,*) '#D5# set nspcpt, log(q(nspcpt)) to',
     *        nspcpt, xspcpt(nspcpt)
          end if
        end if
      end if
c
c  test for maximum extent of semiconvection
c
      if(qmxscn_max.gt.qmxcor) then
        nspcpt=nspcpt+1
	xspcpt(nspcpt)=log10(qmxscn_max)
        write(istdpr,*) '#D6# set nspcpt, log(q(nspcpt)) to',
     *    nspcpt, xspcpt(nspcpt)
      end if
c
c  test for semiconvection border
c
      if(qmxscn.gt.0) then
	if(qmxscn.gt.qmxscn_max) qmxscn_max=qmxscn
        nspcpt=nspcpt+1
	xspcpt(nspcpt)=log10(qmxscn)
        write(istdpr,*) '#D7# set nspcpt, log(q(nspcpt)) to',
     *    nspcpt, xspcpt(nspcpt)
      end if
c
      if(istdpr.gt.0) write(istdpr,
     *  '(/'' Special mesh points set in s/r spcpts.''/
     *  '' q, log(q):''/(1p2e15.7))')
     *  (10.d0**xspcpt(k), xspcpt(k), k = 1, nspcpt)
      return
      end
