      subroutine prtsol(x,y,yp,zk,ap,aq,rhs,nn,iy,intr,ifprt,istosc,
     *  icasex,iform,agey,datmod,ndtmod,ndsmth)
c
c  Note: diagnostic output from calculation of A added 23/2/90
c
c  Output from convection zone for Svend T. P. added 10/4/91
c
c  prints solution with step intr. if istosc = 1 sets quantities for
c  pulsation calculation on file (dsn 1)
c  when ifprt .lt. 0, sets quantities at centre and in atmosphere by
c  calling bcs. for ifprt .ne. -2 this includes resetting quantities in
c  common/ksider/.
c  if ifprt .ne. -1 writes summary on d/s 10
c
c  Modified 20/8/84 to include scaled variables in central boundary
c  conditions
c
c  Modified 4/1/85 to use as dependent variables log(r/1.d11) and
c  log(l/1.d33).
c
c  Modified 25/8/87, so that now uw(1) in common/rnrout/ returns
c  Uw. This is also the quantity printed at the centre.
c  Previously it returned Uw/ln(10).
c
c  21/9/87: implementing modifications from RECKU
c
c  09/11/87: correcting initial setting of dcz.
c     modifying value of iter in call of s/r rhs.
c     note: this has to be thought through more carefully.
c
c  18/12/87: include setting of GONG model file
c
c  18/1/87: add Ne and rX to GONG model file
c
c  Modified 15/3/88 to use consistent numerical constants
c
c  Modified 5/3/90, to increase size of arrays for atmospheric solution
c  from 51 to 201
c
c  Modified 7/3/90, adding parameter icasex to choose optimal setting
c  of r/R in adiabatic pulsation model. icasex = 1 corresponds
c  to using optimal setting. Otherwise old setting is used.
c
c  Modified 2/8/91, to incorporate detailed treatment of CNO cycle.
c
c  Modified 25/7/92, adding reset of adiabatic variable A_4 near
c  centre, to correct for problems near hydrogen exhaustion in
c  models with radiative cores. 
c#ai#  This needs to be fixed properly
c
c  Modified 2/4/93, adding enthalpy and gravitational energy generation
c  rate to GONG output
c
c  Modified 3/4/93, adding integrated gravitational and thermal
c  energy, and gravitational luminosity, to datgng
c
c  Modified 25/7/93, adding convective velocity to GONG output
c
c  Modified 19/7/95, to use general heavy-element abundance
c  as set in common/heavy/. 
c
c  Modified 22/8/95, to include output of diffusion quantities
c
c  Modified 4/8/96, taking out option of writing data to disk with
c  ifwrt = 1 (but retaining ifwrt in common/rhcn/ for consistency)
c
c  Modified 15/9/97, allowing extended iwdopc for use with GH
c  package version 9 (and later)
c
c  Modified 15/10/97, fixing setting of nabla_rad - nabla_ad
c  (i.e., ddrad) in atmosphere.
c
c  Modified 16/2/98, setting adiabatic-oscillation (amdl) variables
c  consistently with rotation. Note that here aa(6,.) has been added,
c  as g/(g tilde), and data(8) has been set to 20.
c
c  Modified 31/3/98 to correct (I hope) storage of Y_H in 
c  GONG model
c
c  Modified 5/1/00, to add more general treatment of convective overshoot
c  from convective envelope and core.
c
c  Modified 23/3/01, to add diffusion coefficients etc to gong file 
c  output.
c
c  Modified 13/6/02, correcting error in setting of q = m/M at lower 
c  boundary of convection zone, for inclusion in datgng
c
c  Modified 6/6/03, resetting Ledoux stability quantity A in regions
c  with 4He burning, for iheccs .gt. 1. Currently the resetting is
c  hard-coded to take place for log(T) .gt. 7.5, outside convective
c  regions.
c
c  Modified 16/9/03, to set CNO variables to yw by using extcno
c  (previously in gong output CNO abundances were set to zero)
c
c  Modified 30/9/03, correcting common/cmpstr/ by including isphec
c  such that idcomp and iccomp are passed correctly
c
c  Modified 21/7/04, adding unformatted output of csum to iducen
c
c  Modified 1/11/04, adding semiconvection diffusion coefficient to
c  gong output.
c
c  Modified 14/5/05, adding call of mixcor to set convective-core
c  properties when called for reading in model and setting quantities.
c
c  Modified 2/11/05, adding output of Gamma_1 derivatives to iddgm1,
c  if iddgm1 gt 0. For the time being use yw(65 - 70) to store the
c  required output quantities.
c
c  ****************************************************************
c
c  Implementation problems:
c
c  #P# System: sundog.ucar.edu. 28/7/04: the line
c  if(.false.) write(6,*) 'Start setting diff. quantities', ...
c  was added to suppress problems with output of diffusion quantities
c  in gong file on sundog.ucar.edu. Similar problems have not been 
c  encountered on tuc47.
c
c  ****************************************************************
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
      implicit double precision (a-h,o-z)
      logical newfil
      logical conv,time0,nosd,notd,timep,dtest,skipt,noder,nscfil
      character cxcno*50, cxcnoc*30
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.cz.d.incl'
      parameter(idr1mx = nspcmx+3, naztmx = nspcmx+3, 
     *  nbmax = nspcmx + 2)
c
      dimension x(*),y(iy,*),yp(iy,*),zk(*),ap(*),aq(*)
      dimension f(ivarmx),fd(ivarmx,ivarmx),zz(ivarmx),
     *  dzdy(ivarmx,ivarmx),alam(1,ivarmx),alamd(1,ivarmx,ivarmx),
     .  h(1),hd(1),data(5),bcc(nbcprv),bccp(nbcprv),
     *  datmod(1),ndtmod(1),datgng(ndtgng),
     *  epscen(idr1mx),ftcen(idr1mx,nspcmx),xcno(nspcmx),cxcno(4),
     *  cxcnoc(4),dgmm1(3),trder(3,3)
      dimension yy(ivarmx),aztus(4)
      common/mxlcn/ c1,c2,etac,phc,alfa,tprfct,iconcs,imxlng,iturpr
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,jnc,frcf(6),frcl(6)
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor
      common/rhcn/ a1,a2,a3,a4,zdummy,nvar,ifwrt
      common/ln10/ amm
      common/clshft/ alshft
      common/heavy/ zatmos, zhc, zh(1)
      common/totmss/ amass, rscmm
      common/logf/ flc
      common/consts/ av
      common/cntval/ rcnoc,uwc,dradc
      common/prtvar/ rhl,akl,akt,akp,akx,eps(idr1mx),dr,dac,conv
      common/cnvovs/ icnvos, jcnvos, clcovs, cldovs, alphos,
     *   rczl, rczlfx, rcnvos, qlcnos, rhobos, dmrovs, fmrovs
      common/cnvout/ ddrad,rr,a,ddacad,amach,cnvdum(15),pturb
      common/cvcpsc/ cvcpst(20)
      common/comgrt/ isprot, irotcn, a5, riner0, riner, omgrt0,
     *  omgrot(1)
      common/comgrp/ isprtp, irotcp, omgrtp(1)
      common/he3eql/ xhe3eq,ieqhe3,iequs
      common/he3fdg/ agesh,ifdhe3
      common/bccn/ b1,b2,b3,b4,nnb,iveb,icncbc,ihe3bc
      common/compvr/ cvr(icvrmx,1)
      common/cstcvr/ icvh1, icvhe4, icvhe3, icvc12, icvc13, icvn14, 
     *  icvo16
      common/cstyst/ iyh1, iyhe3, iyc12, iyc13, iyn14, iyo16, 
     *  iyche4, iycc12
      common/cstazt/ iah1, iahe3, iac12, iac13, ian14, iao16, 
     *  iache4, iacc12
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec, idcomp, iccomp
      common/engcnt/ xhzlm1, xhzlm2, nxhzer, icnocs, nspec,
     *  iheccs, nspect
      common/enggrv/ epsg(nnmax)
      common/ksider/ al0,al2,aztst(naztmx),axst(naztmx)
      common/caztax/ azt(nbmax),ax(nbmax)
      common/ebcder/ epsc
      common/rbcder/ rhoc,rhds(3)
      common/catmos/ natm,matm,vatm(5,201),patm(5,201)
      common/work/ pl(nnmax),yw(iywstr,1)
      common/sooner/ xw(1)
      common/eqstd/ xi(4),ane(10),rho(20),ht(20),pt(20),cp(4),
     .  dad(4),dlt(4),gm1,tprh,trhp,rhxp,gmm1(4)
      common /eqstcl/ iradp,mode,dtest,skipt,noder
      common/rnrout/ rnzt(10),rnuw(10)
      common/cmtime/ age, time0
      common/noiter/ iter, ntime
      common/step/ dt
      common/cxhcnc/ compc(nspcmx)
      common/cdiffu/ idiffus, itbdif, rhob, rmxdif, tlfdif, dtldif,
     *  ctdfmx, rctdif, rctdf1
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev, syear
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  dataset designators
c
      common/cdsdsg/ idstrl,idshvz,idszab,idsevl,idssum,idscen,idsnov,
     *  idsoov,idsgng,idsgsm,idstm1,idstm2,iducen,iddgm1
c
c  initialize headers for CNO abundances
c
      data cxcno /'X(N14), X(O16)',
     *            'X(C12), X(C13), X(N14), X(O16)',
     *            'X(C12), X(C13), X(N14)',
     *            'X(C12), X(C13), X(N14), X(O16)'/
      data cxcnoc /'X14,c, X16,c',
     *             'X12,c, X13,c, X14,c, X16,c',
     *             'X12,c, X13,c, X14,c',
     *             'X12,c, X13,c, X14,c, X16,c'/
c
      external rhs
c
      save
c
      equivalence (bcc(1),al0)
c
      data nmodos /0/
c
      if(istdpr.gt.0) then
	write(istdpr,'(/'' Enter prtsol''/)')
c
        write(istdpr,*) 'nrdtmd, nidtmd, ndtgng, iywstr:'
        write(istdpr,*)  nrdtmd, nidtmd, ndtgng, iywstr
      end if
c
c  zero xcno
c
      call zero(xcno,nspcmx)
c
c  zero datgng
c
      call zero(datgng,ndtgng)
c
      iz=ivarmx
      iyw=iywstr
c
c  surface radius
c
      rhats=1.d1**y(1,1)
      rs=1.d11*rhats
c
c  test for call of bcs, to set up quantities at centre and in
c  atmosphere
c
      if(ifprt.lt.0) then
c  store previous values for later resetting
        if(ifprt.eq.-2) call store(bcc,bccp,nbcprv)
        iterp=iter
        iter=2
        timep=time0
        time0=.true.
c
        matm=0
c
        call bcs(x(1),x(nn),y(1,1),y(1,nn),zk,ap,aq,f,
     .    dzdy,iz,iz,nn)
c
c  reset to previous values
c
        if(ifprt.eq.-2) call store(bccp,bcc,nbcprv)
        iter=iterp
        time0=timep
c
      end if
c
      pc=1.d17*aztst(1)
c
c  differentiate x wrt r/rs numerically
c
      n1=1
      if(matm.gt.0) n1=natm
      n0=n1-1
      do 7 n=1,nn
      n0=n0+1
      xw(n0)=(1.d1**y(1,n))/rhats
    7 yw(2,n0)=y(5,n)
c
      call derive(xw(n1),yw(2,n1),yw(1,n1),nn,iyw,iyw,1,1)
c
c  equation of state at centre
c
      nosd=.true.
      notd=.true.
      skipt=.false.
      fl=flc
      tl=log10(1.d7*aztst(2))
      xh=aztst(3)
      yh=1-xh-zhc
    8 call eqstf(fl,tl,xh,yh,zhc,nosd,notd)
c
      rhoc=rho(1)
      rhl=log10(rhoc)
      call opact(rhl,tl,xh,zhc,aklc,akrc,aktc,akxc)
      akc=10.d0**aklc
      pc=pt(1)
      dadc=dad(1)
c
c  call engenr at centre 
c  Note: details of setting of composition (by using either 
c  the external array aztst in common/ksider/ or the internal
c  array azt in common/caztax/) are still a little uncertain
c  and ought to be checked. The effect of (or need for) this
c  call should also be tested.
c
      tc=1.d7*aztst(2)
      tcl=log10(tc)
      nosd=.true.
      if(ifdhe3.eq.1.or.(ieqhe3.eq.0.and.ihe3bc.eq.1)) then
        call engenr(fl,tcl,aztst(3),yh,zhc,epscen,ftcen,nspcmx,nosd)
      else
        call engenr(fl,tcl,azt(3),yh,zhc,epscen,ftcen,nspcmx,nosd)
      end if
      if(idgeng.ge.2.and.istdpr.gt.0) then
        write(istdpr,*) ' after calling engenr in prtsol'
        call dmprcn
      end if
c
c  cdp changed 18/1/84 to correct for changed definition of
c  a1 and a2 in s/r mnevol. this also causes change in
c  calculation of dpdx below.
c
   10 cdp=1.d-11*rs*a2/a1
      caf=rhoc/pc
      cbv=rhoc/pc
c
      intrdo=intr
c
      if(istosc.ne.1) go to 11
c
c  values at centre for oscillations
c
      nn1=nn+1
      nnosc=nn1
      if(matm.gt.0)  nnosc=nn+natm
c
      call zero(yw(1,nnosc),iyw)
c
      yw(1,nnosc)=pt(1)
      yw(3,nnosc)=rho(1)
      yw(5,nnosc)=gm1
      yw(6,nnosc)=1.d7*aztst(2)
      yw(7,nnosc)=tprh
      yw(8,nnosc)=trhp
      yw(9,nnosc)=epsc
      yw(10,nnosc)=xh
      yw(11,nnosc)=yh
      yw(12,nnosc)=zhc
      yw(13,nnosc)=aztst(4)
      yw(14,nnosc)=flc
      yw(16,nnosc)=2.d17*axst(1)
      yw(17,nnosc)=-60
      yw(18,nnosc)=dad(1)
      yw(22,nnosc)=rhxp
      yw(31,nnosc)=cp(1)
      yw(32,nnosc)=dlt(1)
      yw(34,nnosc)=akc
      yw(35,nnosc)=-60
      yw(36,nnosc)=tl
      yw(37,nnosc)=-60
      yw(38,nnosc)=ane(1)
      yw(39,nnosc)=ftcen(1,1)
c
c  yw(41 - 45,nnosc) set to CNO abundances below
c
      yw(46,nnosc)=dradc-dad(1)
      yw(47,nnosc)=ht(1)
c
c  set central angular velocity to value at innermost non-zero point
c
      if(isprot.ne.0) yw(59,nnosc)=omgrot(nn)
c
      intrdo=1
c
c  test for setting Gamma_1 derivatives
c
      if(iddgm1.gt.0) then
	call trfder(fl, tl, xh, yh, zhc, gmm1(2),dgmm1,trder,3,0)
	do i=1,3
	  yw(64+i,nnosc)=dgmm1(i)
	  yw(67+i,nnosc)=trder(2,i)
        end do
      end if
c
c  Test for (re)setting convection zone boundaries.
c  This is flagged by inc (in common/convpt/) .le. 0.
c  Should probably only be used when setting oscillation variables
c  for a model read in.
c
c  Note that when convection zone boundaries are reset, we possibly
c  also need to reset boundary of convective core by calling
c  s/r conmxc
c
   11 if(inc.eq.-1) then
        itrrhs=1
      else
        itrrhs=0
      end if
c  depth of surface convection zone
      if(inc.gt.0.and.x(nf(1)).ge.-1.d-5) then
        n=nl(1)
        frc=frcl(1)
        n1=n
        if(matm.gt.0) n1=n+natm-1
        dcz=1-frc*xw(n1)-(1-frc)*xw(n1+1)
      else
	dcz=0
      end if
c
      if(istdpr.gt.0) then
        write(istdpr,115)
        if(icnocs.ge.1) write(istdpr,116) cxcno(icnocs)
        write(istdpr,117)
        if(idgeng.ge.2) then
          write(istdpr,*) 
     *      ' before starting loop through atmosphere in prtsol'
          call dmprcn
        end if
      end if
c
c  test whether atmosphere is included  in printout
c
      n=1-intrdo
      n1=n
      inatm=0
      if(matm.gt.0) then
c  include atmosphere
        inatm=1
        call store(y,yy,nvar)
        xhe3 = cvr(3,1)
        if(icnocs.ge.1) call extcno(xcno,1)
c
        if(istdpr.gt.0) write(istdpr,118)
        fcaxx=1./(amsun*amass*amm)
        ratm0=rs+vatm(5,natm)
c  zero yw(1,.), set xw, yw(2,.) and CNO abundances in atmosphere
        natm1=natm-1
        do 11005 k=1,natm1
        xw(k)=(ratm0-vatm(5,k))/rs
        yw(1,k)=0
        yw(2,k)=yy(5)
        if(icnocs.gt.0.or.iheccs.gt.0) call store(xcno,yw(42,k),4)
c
11005   continue
c
      end if
c
c  step through  atmosphere or interior
c
11010 n=n+intrdo
      n1=n1+intrdo
c
      if(inatm.ne.0) then
c
c  point  in atmosphere
c
c  test for final point
c
        if(n.ge.natm) then
          inatm=0
          n1=natm
          n=1
          if(istdpr.gt.0) write(istdpr,119)
        else
c  set  up for atmospheric point
          xx=fcaxx*(vatm(4,natm)-vatm(4,n))
          yy(1)=log10((ratm0-vatm(5,n))/1.d11)
          yy(2)=vatm(2,n)
          yy(3)=vatm(3,n)
	  zhh=zatmos
	  nrhs=1
          go to 12
        end if
      end if
c  set up for interior point
      if(n.gt.nn) go to 40
      nrhs=n
      xx=x(n)
      call store(y(1,n),yy,nvar)
      xhe3 = cvr(3,n)
      zhh=zh(n)
      if(icnocs.ge.1) call extcno(xcno,n)
c
c  call right hand  side subroutine
c
   12 timep=time0
      time0=.false.
      iterp=iter
      iter=0
c
      if(idgeng.ge.2.and.istdpr.gt.0) then
        write(istdpr,*) ' before entering rhs'
        call dmprcn
      end if
c
c  as a fairly crude hack, force calling of MR subadiabatic overshoot
c  if relevant
c
      icnvop=icnvos
      if(icnvos.eq.2.and.jcnvos.eq.2.and.inatm.eq.0) then
	if(nrhs.eq.1.and.istdpr.gt.0) 
     *    write(istdpr,*) 'Force MR subadiabatic overshoot in prtsol'
	icnvos=202
      end if
      call rhs(xx,yy,zk,zz,dzdy,ap,aq,f,fd,alam,alamd,h,hd,
     .  iz,iz,1,1,nrhs,itrrhs)
      icnvos=icnvop
c
      if(idgeng.ge.2.and.istdpr.gt.0) then
        write(istdpr,*) ' after calling rhs'
        call dmprcn
      end if
c
      time0=timep
      iter=iterp
c
c  in atmosphere reset actual and super-adiabatic gradient
c  and zero gravitational energy
c
      if(inatm.ne.0) then
        dac=patm(5,n)
	dr=dac
        ddacad=dac-dad(1)
	ddrad=ddacad
        conv=ddacad.gt.0
c
	epsgrv=0
c
c  in case of convection in atmosphere, we have to think carefully
c  about appropriate output, but later
c
        if(conv) then
          ff=1
	  amach=0
        end if
c
c  else set epsgrv from array
c
      else
	epsgrv=epsg(n)
      end if
c
      ff=dac/dr
c
c  dpdx changed 18/1/84 to correct for changed definition of
c  a1 and a2 in s/r mnevol.
c  changed again 4/1/85 in connection with redefinition of yy(1)
c
      dpdx=-cdp*rho(1)*1.d1**(xx-2*yy(1))
      pp=1.d1**zz(2)
      if(abs(yw(1,n1)).lt.1.e-10.or.n.lt.20) yw(1,n1)=0
c
c  brunt-vaisala frequency and acoustic frequency
c  (Note: with turbulent pressure brunt-vaisala frequency is
c   note set correctly)
c
      bv=cbv*dpdx*(yw(1,n1)*rhxp/yy(5)
     .  -dpdx*dlt(1)*ddacad/pp)/rho(1)
      af=caf*pp*gm1/(rho(1)*xw(n1)*xw(n1))
c
      amuel=ane(1)/av
      akk=10.d0**akl
      rl=11+yy(1)
      alum=1.d33*(10.d0**yy(4)-alshft)
      if(alum.gt.0) then
	all=log10(alum)
      else
	all=-60.d0
      end if
c
c  test for convective output
c
      if(conv) then
        vcon=sqrt(pp*gm1/rho(1))*amach
      else
        amach=0
        vcon=0
        ddacad=dac-dad(1)
      end if
c
      if(mod(n-1,intr).eq.0.and.istdpr.gt.0) then
c
        write(istdpr,125) n,xx,rl,zz(2),yy(3),rhl,yy(2),
     .    all,yy(5),xhe3,xhe3eq,dr,dac,
     .    xw(n1),xi,amuel,akk,eps(1),bv,af,dlt(1),dad(1)
	if(icnocs.ge.1) then
	  write(istdpr,126) (xcno(i),i=1,icvcno)
        end if
        if(conv) then
          write(istdpr,128) ddacad,ff,amach
        end if
      end if
   30 continue
      if(istosc.ne.1) go to 11010
c
c  store data in yw for oscillation file
c  note: derivatives are set wrt to r/1.d11.
c  rs is surface radius in cm
c
      xh=yy(5)
      t=1.d1**yy(3)
      rr=rs*xw(n1)
      dlmr=1.d11/(rr*f(1))
      dxdr=1.d11*yw(1,n1)/rs
      yw(1,n1)=pp
      yw(2,n1)=pp*f(2)*dlmr
      yw(3,n1)=rho(1)
c
c  set convective stability parameter A.
c  Note: with turbulent pressure this must be reset later,
c  to include correction for gradient of turbulent pressure
c
      yw(33,n1)=dlt(1)*f(2)*dlmr*ddacad-rhxp*dxdr/xh
      yw(4,n1)=rho(1)*(yw(2,n1)/(gm1*pp)-yw(33,n1))
      yw(33,n1)=1.d-11*rr*yw(33,n1)
c
c  store composition variables
c
      yw(41,n1) = cvr(2,n)
      if(icnocs.gt.0.or.iheccs.gt.0) call store(xcno,yw(42,n1),4)
c
      yw(5,n1)=gm1
      yw(6,n1)=t
      yw(7,n1)=tprh
      yw(8,n1)=trhp
      yw(9,n1)=eps(1)
      yw(10,n1)=xh
      yw(11,n1)=1-xh-zhh
      yw(12,n1)=zhh
      yw(13,n1)=xhe3
      yw(14,n1)=yy(2)
      yw(15,n1)=t*dac*f(2)*dlmr
      yw(16,n1)=yw(2,n1)*(dlmr+yw(4,n1)/rho(1)-2.d11/rr)
      yw(17,n1)=xx
      yw(18,n1)=dad(1)
      yw(19,n1)=ddacad
      yw(20,n1)=ff
      yw(21,n1)=dxdr
      yw(22,n1)=rhxp
      yw(23,n1)=dac
      yw(24,n1)=0
      yw(25,n1)=1.d1**yy(4)-alshft
      yw(26,n1)=1.d-11*rr
      yw(27,n1)=amach
      yw(28,n1)=rnuw(1)
      yw(29,n1)=eps(1)
      yw(30,n1)=xhe3eq
      yw(31,n1)=cp(1)
      yw(32,n1)=dlt(1)
c
c  note that yw(33,n1) was set above.
c
      yw(34,n1)=akk
c
c  store evolution variables appropriately 
c
      yw(35,n1)=yy(1)
      yw(36,n1)=yy(3)
      yw(37,n1)=yy(4)
c
      yw(38,n1)=ane(1)
      yw(39,n1)=f(5)
c
c  set optimal r/R into yw(40,.)
c
      yw(40,n1)=10.d0**(yy(1) - y(1,1))
c
c  yw(41 - 45,nnosc) set to CNO abundances below
c
      yw(46,n1)=ddrad
      yw(47,n1)=ht(1)
      yw(48,n1)=epsgrv
      yw(49,n1)=vcon
      yw(50,n1)=pturb
c
c  set log10(p) and log10(pg), for later differentiation
c  Also log(rho), for case with 4He burning
c
      pl(n1)=zz(2)
      yw(57,n1)=log10(pp-pturb)
      yw(60,n1)=log(rho(1))
c
c  test for setting angular velocity in yw(59,.)
c
      if(isprot.ne.0) then
	if(inatm.ne.0) then
	  yw(59,n1)=omgrot(1)
        else
	  yw(59,n1)=omgrot(n)
        end if
      end if
c
c  test for setting Gamma_1 derivatives
c
      if(iddgm1.gt.0) then
	fl=yy(2)
	tl=yy(3)
	xhh=yy(5)
	yhh=1.d0-xhh-zhh
c
	call trfder(fl, tl, xhh, yhh, zhh, gmm1(2),dgmm1,trder,3,0)
	do i=1,3
	  yw(64+i,n1)=dgmm1(i)/gm1
	  yw(67+i,n1)=amm*trder(2,i)
        end do
      end if
c
      go to 11010
c
   40 continue
c
c  test for updating yw(33,.) to include effects of turbulent
c  pressure
c
      if(iturpr.gt.0) then
	call derive(pl,yw(57,1),yw(58,1),nn,iyw,iyw,1,1)
	do 41 n=1,nnosc
   41   yw(33,n)=yw(33,n)-(1.d0-yw(58,n))*yw(26,n)*yw(2,n)*yw(7,n)/
     *           (yw(1,n)*yw(8,n))
      end if
c
c  test for updating yw(33,.) in case of 4He burning
c  (note that we do not have the proper thermodynamic derivative
c  of rho wrt Y, so need to differentiate rho wrt r instead)
c
      if(iheccs.gt.0) then
	call derive(xw,yw(60,1),yw(61,1),nn,iyw,iyw,1,1)
	if(istdpr.gt.0) write(istdpr,
     *    '(//'' n, log(q), log(T), der(log rho), old A, new A:'')')
	do n=1,nn
	  if(yw(36,n).ge.7.5.and.yw(46,n).lt.0) then
	    dlrho=1.d11*yw(26,n)*yw(61,n)/rs
	    anew=yw(2,n)*yw(26,n)/(yw(1,n)*yw(5,n))-dlrho
	    if(mod(n,5).eq.0.and.istdpr.gt.0) write(istdpr,
     *        '(i5,2f10.5,1p3e13.5)')
     *        n, x(n), yw(36,n), dlrho, yw(33,n), anew
	    yw(33,n)=anew
          end if
        end do
      end if
c
c  test for setting mixed-core boundary
c
      if(ifprt.lt.0.and.inc.gt.0.and.nl(inc).eq.nn) then
        iextrp=10
        iterp=iter
        iter=1
        call mixcor(x,y(1,1),iy,nn,compc,iextrp,0)
        iter=iterp
      end if
c
c  set (for the time being) central epsgrv to be value at
c  preceding point
c
      yw(48,nnosc)=yw(48,nnosc-1)
c
c  central values
c
      if(ifprt.ne.-1) then
c
c  write global values to file
c
        if(newfil(idssum).and.nscfil(idssum)) then
          write(idssum,130)
	  if(icnocs.ge.1) write(idssum,132) cxcnoc(icnocs)
	  write(idssum,133)
        end if
        if(newfil(idscen).and.idscen.gt.0) then
	  if(icnocs.lt.1) then
            write(idscen,135)
          else
            lxcnoc=length(cxcnoc(icnocs))
            write(idscen,136) cxcnoc(icnocs)(1:lxcnoc)
	  end if
        end if
c
        te=1.d1**y(3,1)
        als=1.d33*(1.d1**y(4,1)-alshft)
        r1=(1.d1**y(1,nn))/rhats
        amc=1.d1**x(nn)
        alc=1.d33*1.d1**y(4,nn)
c
c  set unscaled central values
c
        aztus(1)=1.d17*aztst(1)
        aztus(2)=1.d7*aztst(2)
        aztus(3)=aztst(3)
        aztus(4)=aztst(4)
c
c  when convection zone boundaries have been reset, redetermine dcz,
c
        if(itrrhs.eq.1) then
	  if(inc.gt.0.and.x(nf(1)).gt.-1.d-5) then
            n=nl(1)
            frc=frcl(1)
            n1=n
            if(matm.gt.0) n1=n+natm-1
            dcz=1-frc*xw(n1)-(1-frc)*xw(n1+1)
c
	    if(nl(inc).eq.nn) then
	      call conmxc(x,y,iy,inc,nn,2,0)
            end if
	  else
	    dcz=0
          end if
        end if
c
	if(nl(inc).eq.nn) then
c
c  set mass and radius in convective core, if any
c
	  n=nf(inc)
          n1=n
          if(matm.gt.0) n1=n+natm-1
	  frc=frcf(inc)
          qlcc=frc*x(n) + (1-frc)*x(n-1)
	  qcc = 10.d0**qlcc
	  xcc = frc*xw(n1) + (1-frc)*xw(n1-1)
	else
	  qcc = 0
	  xcc = 0
        end if
c
c  output quantities relevant to convective-core boundary
c
	datgng(18) = qcc
	datgng(19) = xcc
c
c  mass and radius in possible core-mixed region
c
	datgng(29)=rmxcor
	datgng(30)=qmxcor
c
c  set quantities at convective region boundaries into datgng
c
	if(inc.gt.0) then
	  do 50 i=1,inc
	  kf=nf(i)
	  kl=nl(i)
          if(matm.gt.0) then
	    n1=kf+natm-1
	    n2=kl+natm-1
          else
            n1=kf
            n2=kl
	  end if
	  frf=frcf(i)
	  frl=frcl(i)
c
	  if(kf.eq.1) then
	    qlcf=x(kf)
	    xcf =xw(n1)
	  else
            qlcf=frf*x(kf) + (1-frf)*x(kf-1)
	    xcf = frf*xw(n1) + (1-frf)*xw(n1-1)
	  end if
	  qcf = 10.d0**qlcf
c
	  if(kl.eq.nn) then
	    qlcl=x(kl)
	    xcl =xw(n2)
	  else
            qlcl=frl*x(kl) + (1-frl)*x(kl+1)
	    xcl = frl*xw(n2) + (1-frl)*xw(n2+1)
            if(istdpr.gt.0) then
	      write(istdpr,*) 
     *          'Set lower convective boundary with kl, n2 =', kl,n2
	      write(istdpr,*) 
     *          'frl, 1-frl, x(kl), x(kl+1), xw(n2), xw(n2+1) ='
	      write(istdpr,*)  
     *          frl, 1-frl, x(kl), x(kl+1), xw(n2), xw(n2+1)
	      write(istdpr,*) 'qlcl, xcl', qlcl, xcl
	    end if
	  end if
	  qcl = 10.d0**qlcl
c
	  i4=4*(i-1)
	  datgng(31+i4)=xcf
	  datgng(32+i4)=qcf
	  datgng(33+i4)=xcl
   50     datgng(34+i4)=qcl
c
	end if
c
c  store surface composition into datgng
c
	datgng(55)=y(5,1)
	datgng(56)=zh(1)
c
        if(nscfil(idssum)) then
          write(idssum,140) agey,rs,te,als,dcz,(aztus(i),i=1,4),rhoc,
     .      epsc,rcnoc,uwc,r1,amc,alc,cvr(2,nn+1)
	  if(icnocs.ge.1) then
            call extcno(xcno,nn+1)
	    write(idssum,142) (xcno(i),i=1,icvcno)
          end if
        end if
c
c  set first part of datgng by writing to and reading 
c  from scratch file
c
        rewind 91
        write(91) agey,rs,te,als,dcz,(aztus(i),i=1,4),rhoc,epsc,akc
        rewind 91
        read(91) (datgng(i),i=1,12)
c
c  file central values and store in common /csum_param/
c
	if(idscen.gt.0) then
          if(icnocs.lt.1) then
            write(idscen,150) amass, (datgng(i),i=1,12), 
     *        dradc, dadc,qmxcor,rmxcor,cvr(2,nn+1),y(5,1),zh(1)
            if(iducen.ge.0) write(iducen) amass, (datgng(i),i=1,12), 
     *        dradc, dadc,qmxcor,rmxcor,cvr(2,nn+1),y(5,1),zh(1)
          else
            call extcno(xcno,nn+1)
            write(idscen,150) amass, (datgng(i),i=1,12), dradc,
     *        dadc,qmxcor,rmxcor,cvr(2,nn+1),(xcno(k),k=1,icvcno),
     *        y(5,1),zh(1)
            if(iducen.ge.0) write(iducen) amass, (datgng(i),i=1,12), 
     *        dradc, dadc,qmxcor,rmxcor,cvr(2,nn+1),
     *        (xcno(k),k=1,icvcno),y(5,1),zh(1)
          end if
	  call flush(idscen)
	end if
        call store_csum(ntime, icnocs, amass, datgng, 12, 
     *      dradc, dadc,qmxcor,rmxcor,cvr(2,nn+1),xcno,icvcno,
     *      y(5,1),zh(1))
c
      end if
c
c  test for setting oscillation variables
c
      if(istosc.ne.1) then
        if(istdpr.gt.0) write(istdpr,'(/'' Exit prtsol''/)')
	return
      end if
c
c  write data on disc for oscillations
c
c  old oscillation variables
c
      nmodos=nmodos+1
      n1=nnosc-1
      yw(19,nnosc)=yw(19,n1)
      yw(20,nnosc)=yw(20,n1)
      yw(23,nnosc)=yw(23,n1)
c
c  store central abundances in yw and datgng
c
      yw(41,nnosc)=cvr(2,nn+1)
      datgng(21)=cvr(2,nn+1)
      if(icnocs.gt.0.or.iheccs.gt.0) then
	call store(xcno,yw(42,nnosc),4)
	call store(xcno,datgng(22),4)
      end if
c
c  store second derivatives at centre in yw(24,.)
c
      yw(24,nnosc)=2.d17*axst(1)
      yw(24,nnosc-1)=2.*(rhds(1)*axst(1)+rhds(2)*axst(2)
     *                  +rhds(3)*axst(3))
      yw(24,nnosc-2)=2.d7*axst(2)
      yw(24,nnosc-3)=2.*axst(3)
      yw(24,nnosc-4)=2.*axst(4)
c
      n1=1
      if(matm.gt.0) n1=natm
c
      data(1)=amass*amsun
      data(2)=yw(25,n1)/3.9
      data(3)=yw(6,n1)
      data(4)=rs/1.d11
      data(5)=agey
c
      nn2=nnosc+1
      if(nscfil(idsoov)) then
        write(idsoov) nmodos,nnosc,data,((yw(k,nn2-n),k=1,26),
     *    n=1,nnosc)
c
        if(istdpr.gt.0)
     *    write(istdpr,170) nmodos,data(5),data(2),data(3),data(4)
      end if
c
c  set new variables
c
c  test for smoothing at atmospheric-interior matching
c
      nsm3=0
      nsm1=0
      nsm2=0
      if(ndsmth.gt.0.and.matm.gt.0) then
        nsm1=nn-ndsmth
        nsm2=nn+ndsmth
      end if
c
      call newvar(nmodos,nnosc,data,yw,iyw,0,2,intr,nsm1,nsm2,
     *  nsm3,icasex,idsnov)
c
c  set integrated energies and luminosities
c
      call seteng(x,y,yw,datmod,nnosc,iy,iyw,utot,grvtot,algtot)
c
      datgng(26)=utot
      datgng(27)=grvtot
      datgng(28)=algtot
c
c  GONG model variables
c
      call stgong(nmodos,nnosc,data,datmod,ndtmod,datgng,
     *  x,y,yw,iy,iyw,iform,idsgng)
c
c  Derivatives of Gamma_1
c
      if(iddgm1.gt.0) call stdgm1(nnosc,yw,iyw,iddgm1)
c
      if(istdpr.gt.0) write(istdpr,'(/'' Exit prtsol''/)')
      return
  115 format(//' Detailed solution.'//
     *         ' ******************'/
     *  ' n,log m, log r, log p, log T, log rho, log f, log L,',
     .  ' X, X3, X3eq, rad. grad.T, actual grad.T'/
     .  5x,' r/Rs, xi(1-4), 1/mue, opacity, eps, Br.-Vai. fr, ac. fr.,',
     .  ' delta, dad')
  116 format(18x,a)
  117 format(17x,' In convection zone: grad.T - ad. grad.T,',
     .  ' radiative/total flux, Mach number:'/)
  118 format(/' atmospheric solution:'/)
  119 format(/' interior solution:'/)
  125 format(/i5,1pe12.4,0p7f10.5,1p2e11.4,e10.3,0pf8.5/
     .  5x,f12.8,0p4f10.6,1p3e10.3,2e11.3,0pf10.5,f8.5)
  126 format(17x,7f10.5)
  128 format(17x,1p3e10.3)
  130 format(//////' rs (cm), te (k), ls (cgs),',
     .  ' (depth of outer c. z.)/rs,',
     .  ' pc (cgs), tc (k), Xc, X3c, rhoc (cgs), epsc (cgs)'/
     .  ' (ecno/etot)c, uwc,',
     .  ' values at central meshpoint: r/rs, m/mtot, l (cgs), Yc')
  132 format(1x,a)
  133 format(/)
  135 format('# Summary of evolution calculation.'/'#'/
     *  '# M/Msun, age, R, Teff, L, Dcz, pc, Tc, Xc, X3c,',
     *  ' rhoc, epsilonc, kappac, dradc, dadc, qcc, xcc, Yc',
     *  ' Xs, Zs:'/'#')
  136 format('# Summary of evolution calculation.'/'#'/
     *  '# M/Msun, age, R, Teff, L, Dcz, pc, Tc, Xc, X3c,',
     *  ' rhoc, epsilonc, kappac, dradc, dadc, qcc, xcc, Yc, ',a,
     *  ', Xs, Zs:'/'#')
  140 format(/' age = ',1pe13.5,' years.'/1pe12.4,0pf12.2,
     .  1p8e12.4/6e12.4)
  142 format(0p10f12.6)
  150 format(f8.3,1p2e13.5,0pf10.2,1pe13.5,0pf10.5,1p7e13.5,
     *  0p2f10.5,1p8e13.5)
  170 format(///' model no ',i4,' has been written on disc for',
     .  ' oscillations'/' age =',1pe13.5,' years, alum =',
     .  0pf10.5,' te =',f10.2,' vr =',f10.5/)
      end
      subroutine seteng(x,y,yw,datmod,nnosc,iy,iyw,utot,grvtot,algtot)
c
c  Compute integrated thermal and gravitational energy, and 
c  gravitational luminosity.
c  Assumes that relevant variables have been set in array yw.
c  Uses yw(51-56,.) as work area, as well as common sooner
c
c  Original version: 3/4/93
c
      implicit double precision (a-h, o-z)
      dimension x(*), y(iy,*), yw(iyw,*), datmod(*)
      common/sooner/ xql(1)
      common/totmss/ am, rscmm
      common/ln10/ amm
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev, syear
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  set integrands
c
      amass=am*amsun
      grvfct=cgrav*amass*1.d-11
      nn1=nnosc-1
      do 20 n=1,nn1
      xql(n)=yw(17,n)
      yw(54,n)=amm*10.d0**yw(17,n)
      yw(51,n)=yw(54,n)*(yw(47,n)-yw(1,n)/yw(3,n))
      yw(52,n)=-grvfct*yw(54,n)*10.d0**(xql(n)-yw(35,n))
      yw(53,n)=yw(54,n)*yw(48,n)
   20 continue
c
c  integrate (with respect to mass fraction)
c
      do 30 k=51,53
      k1=k+3
   30 call vinta(xql,yw(k,1),yw(k1,1),nn1,iyw,iyw)
c
c  include approximate contributions from innermost region
c
      yw(54,nnosc)=yw(54,nn1)-yw(51,nn1)/amm
      yw(55,nnosc)=yw(55,nn1)-0.75*yw(52,nn1)/amm
      yw(56,nnosc)=yw(56,nn1)-yw(53,nn1)/amm
c
c  reset to integral from centre to surface
c
      do 35 k=54,56
      do 35 n=1,nnosc
   35 yw(k,n)=yw(k,n)-yw(k,nnosc)
c
c  set final surface values, including the mass
c
      utot=amass*yw(54,1)
      grvtot=amass*yw(55,1)
      algtot=amass*yw(56,1)
c
      if(istdpr.gt.0) write(istdpr,120) 
     *  utot, grvtot, algtot, algtot/datmod(25)
      return
  120 format(//' Total thermal       energy =',1pe13.5,'  ergs'/
     *         ' Total gravitational energy =',  e13.5,'  ergs'/
     *         ' Gravitational luminosity   =',  e13.5,'  ergs/sec =',
     *         e13.5,' Ltot')
      end
      subroutine newvar(nmod,nn,data,aa,iaa,iread,idir,intr,nsm1,nsm2,
     *  nsm3,icasex,idsout)
c  set variables for new adiabatic pulsation programme
c
c  for iread = 1, read controls from namelist /execnv/, and read
c  old oscillation variables from d/s 2.
c  otherwise take old oscillation variables from arguments data and
c  aa.
c  ***** NOTE: This option disabled 7/7/05
c
c  idir = 1: assume that aa(.,1) corresponds to centre and aa(.,nn) to
c            surface (standard case).
c  idir = 2: assume that aa(.,nn) corresponds to centre and aa(.,1) to
c            surface.
c
c  output new model on d/s 3 when iread = 1, on d/s 15 otherwise.
c
c  may reset v rho by calling s/r rsamdl depending on
c  values of nsm1, nsm2, nsm3.
c
c  Modified 16/3/89, assigning storage for x and aan in commons.
c
c  Modified 7/3/90, adding parameter icasex to choose optimal setting
c  of r/R in adiabatic pulsation model. icasex = 1 corresponds
c  to using optimal setting. Otherwise old setting is used.
c  Also use datn(8) as case number, setting it to 1 if icasex = 1.
c
c  Modified 7/7/05, disabling option to read in old oscillation variables.
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      implicit double precision (a-h,o-z)
      include 'engenr.cz.d.incl'
c
      character tail*5
      logical nscfil
      dimension aa(iaa,610),data(5)
      common/xnwvar/ x(1)
      common/anwvar/ datn(8), aan(istrmx,1)
      common/nmbmsh/ nn_evol, nn_adi, ivar_adi
      common/comgrt/ isprot, irotcn, a5, riner0, riner, omgrt0,
     *  omgrot(1)
      common/comgrp/ isprtp, irotcp, omgrtp(1)
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev, syear
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c$nl      namelist/execnv/ age,ifrmat,iwrbas
c
      data tail /'    @'/
c
c
c  test for erroneous call with iread = 1
c
      if(iread.eq.1) then
c
        write(isdtou,90)
        if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,90)
        stop 'prtsol'
      end if
c
      iwrbas=0
c
c  test for direction
c
c
c  test for direction
c
      if(idir.eq.1) then
        ns=nn
        nc=1
        id=1
      else
        ns=1
        nc=nn
        id=-1
      end if
c
c  set new variables
c
      rs=data(4)
      rsi=1./rs
      ami=1./data(1)
      pc=aa(1,nc)
      rhc=aa(3,nc)
      if(aa(24,nc).eq.0) go to 5
      pdd=aa(24,nc)
      rhdd=aa(24,nc+id)
      go to 7
    5 pdd=aa(16,nn)
      r=aa(26,nc+id)
      rhdd=2*(aa(3,nc+id)-rhc)/(r*r)
    7 pdd=-rs*rs*pdd/(pc*aa(5,nc))
      rhdd=-rs*rs*rhdd/rhc
      cu=4*pi*1.d33*ami*rs**3
c
c  step through model
c
      n1=nc
      do 10 n=2,nn
      n1=n1+id
c
      r=aa(26,n1)
c
c  set r/R, depending on icasex
c
      if(icasex.eq.1) then
	xx=aa(40,n1)
      else
        xx=r*rsi
      end if
c
      rho=aa(3,n1)
c..      cq=1.49853d18
      cq=1.d11/cgrav
      q =10.d0**aa(17,n1)
      qt=-cq*r*r*ami*aa(2,n1)/rho
      if(isprot.eq.0) then
	q=qt
      end if
      qtx=qt/(xx*xx*xx)
      qx =q/(xx*xx*xx)
      vg=-aa(2,n1)*r/(aa(1,n1)*aa(5,n1))
      aan(4,n)=aa(33,n1)
      aan(5,n)=cu*rho/qx
      if(isprot.ne.0.and.xx.gt.0) then
	aan(6,n)=q/qt
      else
	aan(6,n)=1.d0
      end if
      x(n)=xx
      aan(1,n)=qtx
      aan(2,n)=vg
      aan(3,n)=aa(5,n1)
c 
c  test for setting angular velocity in omgrtp
c
      if(isprot.ne.0) omgrtp(n)=aa(59,n1)
   10 continue
c  central values
      x(1)=0
      aan(1,1)=cu*aa(3,nc)/3
      aan(2,1)=0
      aan(3,1)=aa(5,nc)
      aan(4,1)=0
      aan(5,1)=3
      if(isprot.ne.0) then
	aan(6,1)=1.d0
	omgrtp(1)=aa(59,nc)
	isprtp=isprot
	irotcp=irotcn
      else
	isprtp=0
      end if
c  other variables
      datn(1)=data(1)
      datn(2)=rs*1.d11
      datn(3)=pc
      datn(4)=rhc
      datn(5)=pdd
      datn(6)=rhdd
      datn(7)=-1
c
c  set case number in datn(8), depending on icasex, and whether or
c  not rotation is included.
c
      if(icasex.eq.1) then
        datn(8)=1
      else
        datn(8)=0
      end if
c
      if(isprot.ne.0) then
	datn(8)=20+datn(8)
	ivarn=6
      else
	ivarn=5
      end if
c
c  test for smoothing
c
      if((nsm1.gt.0.and.nsm2.gt.0.and.nsm1.lt.nsm2).or.(nsm3.gt.0
     *  .and.nsm3.lt.nn)) call rsamdl(x,aan,datn,nsm1,nsm2,nsm3,
     *  nn,istrmx,aan(6,1),istrmx)
c
c  reset A_4 near centre
c
      call rseta4(x,aan,nn,datn,istrmx)
c
c  output
c
      if(istdpr.gt.0) then
        write(istdpr,130) aa(1,ns),aa(10,ns)
        write(istdpr,100) (datn(i),i=1,6)
        do 20 n=1,nn,intr
        qx=aan(1,n)
        vg=max(aan(2,n),1.d-10)
        bv=qx*aan(4,n)
        ac=qx/vg
   20   write(istdpr,110) n,x(n),(aan(i,n),i=1,5),bv,ac
      end if
c
      if(nscfil(idsout))
     *  write(idsout) nmod,nn,datn,(x(n),(aan(i,n),i=1,ivarn),n=1,nn)
      nn_adi=nn
      ivar_adi=ivarn
      return
   90 format(//
     *  ' ***** Error in s/r prtsol: iread = 1 no longer allowed')
  100 format(///' new variables:'//'  m =',1pe13.5,' rs   =',e13.5/
     .  ' pc =',e13.5,' rhoc =',e13.5,' p2 =',e13.5,' rho2 =', 
     .  e13.5///' n,x,qx,vg,vr,a,u,b.v.,acoust.:'//)
  110 format(i5,f10.6,1p5e13.5,5x,2e13.5)
  130 format(//' surface p =',1pe13.5,'    x =',0pf10.5)
      end
      subroutine rseta4(x,aa,nn,data,iaa)   
c  
c  Reset A4 near centre, to correct for problems near
c  end of hydrogen burning
c  Resetting is only applied if A4/x**2 is non-monotonic
c  near centre
c
c  Orignial version: 25/7/92
c  
c  Modified 8/6/03, suppressing resetting in convective core. In this
c  case, also reset data(6) to ensure consistency.
c
      implicit double precision (a-h, o-z)
      dimension x(*),aa(iaa,*),data(*)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  reference values
c
      axc=data(6)-data(5)
      nref=5
      axref=aa(4,nref)/x(nref)**2
c
      ireset=0
      do 10 n=2,nref-1
      ax=aa(4,n)/x(n)**2
      if((axc-ax)*(ax-axref).lt.0) ireset=1
   10 continue
c
      if(ireset.eq.1) then
c
c  suppress resetting in case of convective core
c
	if(axref.lt.0) then
	  write(istdou,105) 
	  if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,105)
	  datan6=data(5)+axref
	  write(istdou,107) data(6), datan6
	  if(istdou.ne.istdpr.and.istdpr.gt.0) 
     *      write(istdpr,107) data(6), datan6
	  data(6)=datan6
	  return
        end if
c
	if(istdpr.gt.0) write(istdpr,110)
	axcoef=(axref-axc)/x(nref)**2
	do 20 n=2,nref-1
	x2=x(n)**2
	ax=aa(4,n)/x2
	axnew=axc+axcoef*x2
	if(istdpr.gt.0) write(istdpr,115) n, x(n), ax, axnew
   20   aa(4,n)=axnew*x2
c
      end if
c  
      return   
  105 format(//
     *  ' Resetting of A4 near centre suppressed in convective core')
  107 format(/ 'Reset data(6). Old, new values =',1p2e13.5)
  110 format(//' Reset A4 near centre. ',
     *   ' n, x, old, new values of A4/x**2:'/)
  115 format(i5,1p3e13.5)
      end  
      subroutine stgong(nmod,nn,data,datmod,ndtmod,datgng,
     *  x,y,yw,iy,iyw,iform,idsgng)
c
c  sets extensive set of model variables, suffient to set the
c  GONG model set, and containing also all evolution variables.
c  note that variables is stored in same order as evolution
c  variables, i.e. with surface at point 1.
c
c  Here nn refers to the total number of points, including the
c  atmosphere. natm gives the number of points in the atmosphere
c
c  output format:
c
c  write(idsgng) (cdata(i),i=1,4),nmod,nn,nvar,(datmod(i),i=1,31),
c    (datgng(i),i=1,30),(bc(i),i=1,54),
c    ((yvar(i,n),i=1,nvar),n=1,nn)
c
c  here cdata is character*80, and contains
c  cdata(1): date and time
c  cdata(2): trial model file name
c  cdata(3): output model file name
c  cdata(4): opacity file name
c
c  datmod is as defined in evolmain (see programme notes).
c
c  datgng contains the following variables:
c  datgng( 1): age
c  datgng( 2): Rs
c  datgng( 3): Teff
c  datgng( 4): Ls
c  datgng( 5): (depth of outermost C.Z.)/Rs
c  datgng( 6): central p
c  datgng( 7): central T
c  datgng( 8): central X
c  datgng( 9): central X3
c  datgng(10): rhoc
c  datgng(11): epsc
c  datgng(12): kappac
c  datgng(13): d2rho/dr2
c
c  the remainder of datgng is zero so far.
c
c  the array yvar contains variables at each mesh point. 
c  currently the number is nvar = 20:
c
c  yvar( 1,n): log q 
c  yvar( 2,n): log(r/1.d11) 
c  yvar( 3,n): log f 
c  yvar( 4,n): log T 
c  yvar( 5,n): log (L/1.d33) 
c  yvar( 6,n): X 
c  yvar( 7,n): X3 
c  yvar( 8,n): T 
c  yvar( 9,n): p 
c  yvar(10,n): rho 
c  yvar(11,n): L
c  yvar(12,n): kappa 
c  yvar(13,n): epsilon 
c  yvar(14,n): GAMMA1 
c  yvar(15,n): ad. grad. 
c  yvar(16,n): delta 
c  yvar(17,n): cp 
c  yvar(18,n): 1/(GAMMA1) (d ln p)/(d ln r) - (d ln rho)/(d ln r) 
c  yvar(19,n): Ne
c  yvar(20,n): rX (rate of change of X).
c
c  original version: 18/12/87
c
c  Modified 20/8/96, including potential output of rotation profile.
c  Note: room has been left in yvar for diffusion of all relevant
c  elements, when turbulent pressure or angular velocity is included.
c  This is wasteful in space, but simplifies organization
c
      implicit double precision (a-h,o-z)
      character*80 cdata, file, timest, loctim
      character*256 copcvr
      character*24 notava
      logical nscfil
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
      parameter( nsdifr=2*nspdmx+2)
      dimension data(5),datmod(*),ndtmod(*),
     *  datgng(*),x(*),y(iy,*),yw(iyw,*)
      dimension cdata(4)
      common/cofile/ nfiles, idsfil(99), file(99), iopen(99)
      common/opccnt/ xhsopc,tsmn,tstr,rhsmn,rhsmx,timx,rhimn,
     .  rhimx,sigstr,inopc,idmopc
      common/opctcl/ alamop,zatmop,tlmopc,dtmopc,fdgopl,iwdopc,
     *  iopcm1,ifdgop
      common/mxlcn/ c1,c2,etac,phc,alfa,tprfct,iconcs,imxlng,iturpr
      common/covsec/ icsove, icsovc, alpove, alpovc, 
     *  qmxove, qmxovc, rmxove, rmxovc, imxove, imxovc, 
     *  nmxove, nmxovc
      common/comgrt/ isprot, irotcn, a5, riner0, riner, omgrt0,
     *  omgrot(1)
      common/opcver/ copcvr
      common/cgngvr/ yvar(igvrmx,nnmax)
      common/compvr/ cvr(icvrmx,1)
      common/cstcvr/ icvh1, icvhe4, icvhe3, icvc12, icvc13, icvn14, 
     *  icvo16
      common/cstyst/ iyh1, iyhe3, iyc12, iyc13, iyn14, iyo16, 
     *  iyche4, iycc12
      common/cstazt/ iah1, iahe3, iac12, iac13, ian14, iao16, 
     *  iache4, iacc12
      common/ksider/ bc(nbcprv)
      common/catmos/ natm,matm,vatm(5,201),patm(5,201)
      common/cdiffu/ idiffus, itbdif, rhob, rmxdif, tlfdif, dtldif,
     *  ctdfmx, rctdif, rctdf1
      common/cdifsv/ csdifr(nsdifr,1)
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec, idcomp, iccomp
c
c  dataset designators
c
      common/cdsdsg/ idstrl,idshvz,idszab,idsevl,idssum,idscen,idsnov,
     *  idsoov,idsgn1,idsgsm,idstm1,idstm2,iducen,iddgm1
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      data notava /'not available'/
c
c  initialize number of variables (may be reset later)
c
      if(isprot.gt.0) then
	nvarg=45
      else if(iturpr.gt.0) then
	nvarg=44
      else
	nvarg=30
      end if
c
c  set date and time to cdata
c
      timest=loctim()
c
c.. This needs to be removed, when loctim is fixed
c
      timest='Dummy_date_and_time'
c
      write(cdata(1),120) ' date and time     : ',timest
c
c  set file names to cdata
c
      call stfile(idstrl,nfin)
      write(cdata(2),120) ' trial model file  : ',file(nfin)
      call stfile(idsevl,nfin)
      if(nfin.gt.0) then
        write(cdata(3),120) ' output model file : ',file(nfin)
      else
        write(cdata(3),120) ' output model file : ',notava
      end if
      if(iwdopc.lt.8) then
        call stfile(inopc,nfin)
        write(cdata(4),120) ' opacity file      : ',file(nfin)
      else
        write(cdata(4),'(a80)') copcvr
      end if
      if(istdpr.gt.0) write(istdpr,'(//(a/))') ' cdata:',cdata
c
c  set second derivatives in datgng
c
      nn1=nn+1
      do 20 i=1,5
   20 datgng(12+i)=1.d-22*yw(24,nn1-i)
c
c  now we are largely ready for output
c
      do 25 n=1,nn
      yvar( 1,n)=yw(17,n)
      yvar( 2,n)=yw(35,n)
      yvar( 3,n)=yw(14,n)
      yvar( 4,n)=yw(36,n)
      yvar( 5,n)=yw(37,n)
      yvar( 6,n)=yw(10,n)
      yvar( 7,n)=yw(13,n)
      yvar( 8,n)=yw(6,n)
      yvar( 9,n)=yw(1,n)
      yvar(10,n)=yw(3,n)
      yvar(11,n)=1.d33*yw(25,n)
      yvar(12,n)=yw(34,n)
      yvar(13,n)=yw(9,n)
      yvar(14,n)=yw(5,n)
      yvar(15,n)=yw(18,n)
      yvar(16,n)=yw(32,n)
      yvar(17,n)=yw(31,n)
      yvar(18,n)=yw(33,n)
      yvar(19,n)=yw(38,n)
      yvar(20,n)=yw(39,n)
      yvar(21,n)=yw(41,n)
      yvar(22,n)=yw(42,n)
      yvar(23,n)=yw(43,n)
      yvar(24,n)=yw(44,n)
      yvar(25,n)=yw(45,n)
      yvar(26,n)=yw(46,n)
      yvar(27,n)=yw(19,n)
      yvar(28,n)=yw(47,n)
      yvar(29,n)=yw(48,n)
   25 yvar(30,n)=yw(49,n)
c
c  set final variables, depending on diffusion case
c
      if(istdpr.gt.0) write(istdpr,*) 'idiffus, idcomp, iccomp =', 
     *  idiffus, idcomp, iccomp
      if(idiffus.eq.1) then
	nvarg=max0(nvarg,31)
	iyh=4+idcomp+iccomp
	do 27 n=1,nn
	if(n.lt.natm.or.n.eq.nn) then
          yvar(31,n)=0
        else
          yvar(31,n)=y(iyh,n-natm+1)
        end if
   27   continue
      else if(idiffus.eq.2) then
	nvarg=max0(nvarg,29+2*idcomp)
	iyh=5+idcomp+iccomp
	do 28 n=1,nn
	if(n.lt.natm.or.n.eq.nn) then
          yvar(31,n)=0
        else
          yvar(31,n)=y(iyh,n-natm+1)
        end if
   28   continue
c
	do 29 j=2,idcomp
	ixi=4+iccomp+j
	iyi=4+idcomp+iccomp+j
	jxi=30+j
	jyi=29+idcomp+j
	do 29 n=1,nn
	if(n.lt.natm) then
          yvar(jxi,n)=y(ixi,1)
          yvar(jyi,n)=0
        else if(n.lt.nn) then
	  n1=n-natm+1
          yvar(jxi,n)=y(ixi,n1)
          yvar(jyi,n)=y(iyi,n1)
        else
          yvar(jxi,n)=y(ixi,n1)
          yvar(jyi,n)=0
        end if
   29   continue
      end if
c
c  test for including diffusion coefficients etc. 
c  (added 23/3/01)
c
      if(idiffus.gt.0) then
c
c..        if(.false.) write(istdpr,*) 'Start setting diff. quantities',
c..     *    ' csdifr(1,300) = ',csdifr(1,300)
c
	do 33 n=1,nn
        if(n.lt.natm) then
	  do 31 i=45,47+2*idcomp
   31     yvar(i,n)=0
	else
c
c  turbulent diffusion coefficient
c
	  n1=n-natm+1
	  yvar(46,n)=csdifr(2*idcomp+1,n1)
c
c  diffusion coefficients and settling velocities
c
          do 32 i=1,idcomp
	  yvar(46+i,n)=csdifr(i,n1)
   32     yvar(46+idcomp+i,n)=csdifr(idcomp+i,n1)
c
c  semiconvection diffusion coefficient
c
	  yvar(47+2*idcomp,n)=csdifr(2*idcomp+2,n1)
        end if
c
   33   continue
c
        nvarg=max0(nvarg,47+2*idcomp)
	if(istdpr.gt.0) write(istdpr,*) 'nvarg set to ',nvarg
      end if
c
c  test for including angular velocity
c
      if(isprot.gt.0) then
	do 35 n=1,nn
	if(n.lt.natm) then
          yvar(45,n)=omgrot(1)
	else if(n.lt.nn) then
          yvar(45,n)=omgrot(n-natm+1)
	else 
          yvar(45,n)=omgrot(nn-natm)
        end if
   35   continue
	nvarg=max0(45,nvarg)
      end if
c
c  test for including turbulent pressure
c  In this case, set max(pturb/p) into datgng(57)
c
      if(iturpr.gt.0) then
	ptpmax=0
	do 40 n=1,nn
	yvar(44,n)=yw(50,n)
   40   ptpmax=max(ptpmax,yvar(44,n)/yvar(9,n))
	datgng(57)=ptpmax
	nvarg=max0(44,nvarg)
      end if
c
c  test for new version of overshoot. In that case, include limits
c  of overshoot regions
c
      if(imxove.gt.0) then
        datgng(61)=qmxove
	datgng(62)=rmxove
      end if
      if(imxovc.gt.0) then
	datgng(63)=qmxovc
	datgng(64)=rmxovc
      end if
c
      if(nscfil(idsgng)) then
        write(idsgng) (cdata(i),i=1,4),
     *    nmod,iform,nn,nrdtmd,nidtmd,ndtgng,nvarg,
     *    nbcprv,(datmod(i),i=1,nrdtmd),(ndtmod(i),i=1,nidtmd),
     *    (datgng(i),i=1,ndtgng),(bc(i),i=1,nbcprv),
     *    ((yvar(i,n),i=1,nvarg),n=1,nn)
        call flush(idsgng)
      end if
c
c  output GONG summary
c
      nnsum=1
      n=nn
      if(nscfil(idsgsm)) then
        write(idsgsm) (cdata(i),i=1,4),
     *    nmod,iform,nnsum,nrdtmd,nidtmd,ndtgng,nvarg,
     *    nbcprv,(datmod(i),i=1,nrdtmd),(ndtmod(i),i=1,nidtmd),
     *    (datgng(i),i=1,ndtgng),(bc(i),i=1,nbcprv),
     *    (yvar(i,n),i=1,nvarg)
        call flush(idsgsm)
      end if
c
      return
  120 format(a20,a60)
      end
      subroutine store_csum(ntime, icnocs, amass, datgng, istdtg, 
     *      dradc, dadc,qmxcor,rmxcor,yhc,xcno,icvcno, xhs, zhs)
c
c  store values for csum in common /csum_param/ or /csum_indiv/ 
c  for access from calling programme
c
c  If isetos = 0 storage is assumed to be for an evolution sequence
c  and is made in /csum_param/.
c  ntime is assumed to give the time step, starting from ntime = 0
c  this is used for the storage, for evolution calculation.
c  If isetos .gt. 0 (setting oscillation variables) storage is made
c  in /csum_indiv/ and nstep_ind is used as a counting index.
c
      implicit double precision (a-h, o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.cz.d.incl'
c
      dimension datgng(*), xcno(*)
      dimension csum_loc(icsum_max)
      common /csum_param/ icsum, nstep, csum_st(icsum_max, nstep_max)
      common /csum_indiv/ icsum_ind, nstep_ind, 
     *  csum_ind(icsum_max, nstep_max)
      common/cevlio/ isetos, iastr
c
      if(isetos.eq.0) then
        nt=ntime+1
      else
	nstep_ind=nstep_ind+1
	nt=nstep_ind
      end if
c
      csum_loc(1) = amass
      do i = 1,istdtg
        csum_loc(i+1) = datgng(i)
      end do
      i1 = 2+istdtg
      csum_loc(i1) = dradc
      csum_loc(i1+1) = dadc
      csum_loc(i1+2) = qmxcor
      csum_loc(i1+3) = rmxcor
      csum_loc(i1+4) = yhc
      if(icnocs.ge.1) then
	i2 = i1+4
	do i=1,icvcno
	  i2 = i2+1
          csum_loc(i2) = xcno(i)
	end do
	i1 = i2+1
      else
	i1 = i1+5
      end if
      csum_loc(i1) = xhs
      csum_loc(i1+1) = zhs
c
      if(isetos.eq.0) then
        icsum = i1+1
	nstep = nt
	do i=1,icsum
	  csum_st(i,nt)=csum_loc(i)
        end do
      else
        icsum_ind = i1+1
	do i=1,icsum_ind
	  csum_ind(i,nt)=csum_loc(i)
        end do
      end if
      return
      end
      subroutine stdgm1(nnosc,yw,iyw,iddgm1)
c
c  Output Gamma_1 derivatives as stored in yw
c
c  Original version: 2/11/05
c
      implicit double precision(a-h, o-z)
      dimension yw(iyw,*)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      write(iddgm1,110)
      do n=1,nnosc
        write(iddgm1,120) n,yw(40,n),yw(6,n),yw(1,n),yw(3,n),
     *    yw(11,n),yw(5,n),(yw(i,n),i=65,70)
      end do
      return
c
  110 format('# n, r/R, T, p, rho, Y,',
     *  ' Gamma_1, (d ln Gamma_1/d ln p), (d ln Gamma_1/d ln rho),',
     *  ' (d ln Gamma_1/d Y), (d ln T/d ln p), (d ln T/d ln rho),',
     *  ' (d ln T/d Y):'/'#')
  120 format(i5,f12.7,1p15e13.5)
      end
