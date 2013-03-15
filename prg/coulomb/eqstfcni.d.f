      subroutine eqstf(fl,tl,x,y,z,nosd,notd)
c
c  This version of the EFF routine includes
c  effects of Coulomb interaction in the Debye-Hueckel approximation,
c  using the F4 from the MHD package.
c
c  Modifications initiated 19/4/90.
c
c  jc-d equation of state updated to be compatible with  dog equation
c  of state routines. notice in particular that common/eqscnt/ has
c  been changed, by the inclusion of anhe0 after anh0. in addition  a
c  dummy subroutine seteqs has been added. furthermore the dimensions
c  in common /potetc/ have been changed, to allow for only up to
c  10 heavy elements.
c
c  controls:
c  ********
c
c  in argument list:
c
c  logical nosd and notd are switches for no calculation of second and
c  third derivatives, respectively.
c
c  in common/eqscnt/:
c
c  ihvz determines treatment of heavy elements.
c  ihvz = 0: assume heavy elements to be fully ionized everywhere
c       = 1: complete treatment of ionization of c and o; first level
c            of ionization of fe. reset abundances to fit results of
c            full treatment of all 10 elements considered.
c       = 2: first level of ionization of all 10 elements.
c       = 3: complete treatment of c, n and o; first level of
c            ionization of remaining elements.
c       = 4: complete treatment of ionization of all 10 elements
c            included.
c
c  iprrad .ne. 0: include pressure and enthalpy of radiation (in
c                 diffusion approximation)
c  iprrad  =   0: do not include pressure and enthalpy of radiation.
c
c  ihmin   =   1: include h-
c         .ne. 1: do not include h-
c
c  these are initialized  in block data bleqst to ihvz = 1,
c  iprrad = 1, ihmin = 0.
c
c  in common/eqdpco/, relating to fudge departure coefficients to
c  mimick nlte effects:
c
c  idpco = 0: do not include departure coefficients
c        = 1: determine departure coefficient for h (returned in
c             bdcoh) such that n(h+)/n(h neutr.) = frhi.
c        = 2: treat h as for idpco = 1. take departure coefficient
c             for he and heavy elements to be bdcoz.
c        = 3: take departure coefficient for h to be bdcoh, and
c             departure coefficient for he and heavy elements to be
c             bdcoz.
c
c  frhi, bdcoh and bdcoz are initialized to 1., and idpco to 0.
c
c  addition 5/1/84: flag iomfll in common/hvomcl/.
c    if iomfll = 0 statistical weight omega for fully ionized state 
c    of heavy elements is set to 15. 
c    this corresponds to the situation before 5/1/84. the
c    origin of this error is unclear.
c    if iomfll .ne. 0 omega for fully ionized state is set to 1,
c    as it should be.
c
c  Controls in common /ccoulm/ 
c    epsdmu: convergence criterion for iteration for Coulomb
c    effect in Saha equation. (Default: 1.e-12)
c    icoulm: determines how Coulomb terms are included:
c      icoulm = -1: do not include Coulomb effects (should correspond
c         to basic eqstf).
c      icoulm = 0: include Coulomb effects fully.
c      icoulm = 1: include Coulomb effects in Saha equations, but not
c         explicit contributions to pressure and enthalpy.
c      icoulm = 2: include explicit contributions to pressure and
c         enthalpy from Coulomb effects, but not effect
c         in Saha equations.
c      icoulm = 10 - 12: as icoulm = 0 - 2, but suppressing 
c         tau correction
c      (Default: 0)
c    iclmit: determines type of Coulomb iteration
c      iclmit = 0: backsubstitution
c      iclmit = 1: Newton iteration
c      (Default: 1)
c    iclmsd; if iclmsd = 1, set include second derivatives of
c      Coulomb terms
c      (Default: 0)
c         
c
c  modification 15/9/86: save statement added in all subroutines.
c     (needed in f77, apparently)
c
c  modification 3/1/88: dimension of dxt1 increased to 20.
c
c  Modification 6/8/88: second derivative of He ionization fractions
c     corrected. Also in commom/eqsout/ is now stored ea before rescaling,
c     for test of derivatives.
c
c  Modification 17/5/90: include possibility of including detailed 
c     ionization of heavy elements, for ihvz = 4. Note that this 
c     involves changing size of arrays in commons /hvredu/ and /xhv/
c
c  Modified 22/5/90 to include all levels of ionization of Ar and Fe.
c     Previously only the first 15 levels were included, the remainder
c     being forced to be unionized
c
c  Modification 4/6/90: include average degree of ionization for
c     each element in array xhvmn(10) at the beginning of
c     common/hviond/
c
c  Modified 4/6/90, adding array gmm1(4) in common/ eqstd/ containing
c     Gamma1 and derivatives.
c
c  Modified 4/6/90, to include numerical second derivatives of
c     Coulomb terms.
c
c  Modified 5/11/95, including option of suppressing tau correction
c
c  Modified 21/6/03: Include common/cdgphs/ kdgeos, kdgopc, kdgeng.
c  If, on input, kdgeos .gt. 0, do not stop on fatal error, but
c  return with kdgeos = -1.
c
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 9/3/90
c
c            *********************************
c
c  results returned in common/eqstd/. in addition value of fl
c  is set into common/eqstfl/.
c
      parameter(nspe = 6, npar = 3, npar2=npar+2)
      implicit double precision (a-h,o-z)
      logical tstl,nosd,notd,cmplio, secder, thrder
      dimension phi(30),hst(10),ex(3),ea(30),xi(30),dxt(4),dxt1(20),
     .  dxt2(4),ddmu(11),
     .  anuh(10),ueh(10),anuhr(23),uehr(23),
     *  ehm(10),aneder(10), sn(nspe), dmucpr(npar), 
     *  eahat(30), dmucc(npar,4), dmucp(npar,4,6), 
     *  pcoulc(4), pcoulp(4,6), hcoulc(4), hcoulp(4,6)
      common/eqscnt/ anh0,anhe0,ihvz,iprrad,ihmin,igndeg
      common/ln10/ amm,amm2,amm3
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm
      common/dmuder/ dmu,dmuf,dmut,dmux,dmuff,dmuft,dmufx,dmutt,dmutx,
     .  dmuxx,idmu
c
c  note: second dimension of dmuc must be at least as great as npar + 2
c
      common/dmucdr/ dmuc(npar,10)
      common/df4der/ df4(npar), d2f4i(npar,npar), d2f4f(npar,npar),
     .          d2f4x(npar),d2f4r(npar),d2f4t(npar),
     .          dp4i(npar),dh4i(npar), p4, h4, dp4dr, dp4dt, dp4dx, 
     .          dh4dr, dh4dt, dh4dx
      common/eqsout/ east(30),xii(30),dne(10),dph(20),he(20),pe(10),
     *  hi(20),pi(10),hh(20),ph(10),hr(20),pr(10), pcoul(10), hcoul(10)
      common/eqstd/ xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp,gmm1(4)
      common/eqstfl/ flcm
      common/eqdpco/ frhi,bdcoh,bdcoz,idpco
      common/xhminc/ xhm(10)
c
c  controls for Coulomb effect.
c  epsdmu is accuracy requirement in backsubstitution 
c  icoulm is used for selectively switching off parts of Coulomb effect
c  and/or suppressing tau correction
c
      common/ccoulm/ epsdmu, icoulm, iclmit, iclmsd, nitdmu, epssdr
      common/cnttau/ notau
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
      common/cversn/ ivrbcs,ivrrhs,ivreos,ivropc,ivreng
c
      common/cdgphs/ kdgeos, kdgopc, kdgeng
c
      external bleqst, blstio
c
      equivalence(ex(1),exh),(ddmu(1),dmu)
c
      data flprev, tlprev, xprev /-1.d10, -1.d10, -1.d10/
      data kdgeos, kdgopc, kdgeng /0, 0, 0/
c
      save
c
      fxp(a)=exp(min(85.d0,max(a,-85.d0)))
      nuu(l,n)=((l-1)*(8-l))/2+n-l+5
c
c  *****************************************************************
c
c  set equation of state version number, depending on icoulm
c
      if(icoulm.eq.-1) then
        icoul1=icoulm
        notau=0
      else
        icoul1=mod(icoulm,10)
        notau=icoulm/10
      end if
c
      if(icoulm.eq.0) then
        ivreos=1
      else if(icoulm.eq.2) then
        ivreos=2
      else if(icoulm.eq.10) then
        ivreos=3
      else if(icoulm.eq.12) then
        ivreos=4
      end if
c
c  *****************************************************************
c
c  for changing internal logics, set logical variables for setting 
c  second and third derivatives
c
      secder = .not. nosd
      thrder = .not. notd
c
c  set y
      y=1-x-z
c  store fl in common
      flcm=fl
c  number of derivatives
      if(thrder) then
        ider=20
      else if(secder) then
        ider=10
      else 
        ider=4
      end if
      iders=min0(ider,10)
c
      iclsdr = 0
c
c  entry point for setting numerical second derivatives of 
c  Coulomb terms
c  
    5 f=1.d1**fl
      t=1.d1**tl
c
      call phder(fl,tl,phi,hst,nosd,notd)
c
c
c  psi and derivatives
      wf= sqrt(1+f)
      psi=2*wf+log(f/((1+wf)*(1+wf)))
      psif=amm*wf
      psff=amm2*f/(2*wf)
c
      zf=x+2*y+anhe0*z
      zf2=zf*zf
      zf3=zf*zf2
c
      ak0=ckh*zf2
c
c  set k*t, in ev and ergs
c
      tkev=ck1*t
      tkergs=ck2*t
c  a0**3*ne/(k*t)
      aa=caa*phi(1)/(tkev*zf3)
c
c  delta mu and derivatives
c
      bmu=tkev+20*ak0
      dmu=aa*bmu
c
      ref=phi(2)
      ret=phi(3)
      dmuf=dmu*ref*amm
      dmut=aa*(bmu*ret-20*ak0)*amm
      dmux=aa*(3*tkev+20*ak0)/zf
c
      if(secder) then
        dmuff=(dmuf*ref+dmu*phi(4)*amm)*amm-psff
        dmuft=(dmut*ref+dmu*phi(5)*amm)*amm
        dmufx=dmux*ref*amm
        dmutt=aa*(20*ak0*(1-2*ret)+bmu*(phi(6)+ret*ret))*amm2
        dmutx=aa*(ret*(3*tkev+20*ak0)-20*ak0)*amm/zf
        dmuxx=aa*(12*tkev+40*ak0)/zf2
      end if
c
      dmu=dmu-psi
      dmuf=dmuf-psif
      idmu=1
c
c  ***************************************
c
c  ionization of heavy elements
c
      if(ihvz.le.0) then
        anuh(1)=anh0
        call zero(anuh(2),9)
        call zero(ueh,10)
c
      else
c
        call hviona(fl,tl,x,y,z,nosd,notd,anuh,ueh,anuhr,uehr)
c
      end if
c
c  test for complete ionization of h and he
c
      cmplio=(dmu-54.4/tkev).gt.19
c  test for use of departure coefficient
      if(idpco .ge.1) cmplio=cmplio.and.frhi.gt.1.e6
      if(cmplio) then
c
c  complete ionization
c
        xi(1) =1
        xi(11)=0
        xi(21)=1
        do 10 i=2,22,10
        jj=i+8
        do 10 j=i,jj
   10   xi(j)=0
c
c  to avoid problems with logics in Coulomb term treatment,
c  initialize iclder 
c
        iclder = 0
c
        go to 30
c
      end if
c
c  A note on logics of Coulomb iteration: iclder = 0 flags for
c  initial passes iterating for Coulomb correction, and iclder = 1
c  flags for final pass setting derivatives
c  Similarly, when setting second derivatives numerically, iclsdr 
c  flags for modifications of variables
c
c  e h, e he, e he+ and derivatives
c
c  test for initializing Coulomb terms or using previous values
c  Note: for icoulm = 2 switch off effect of Coulomb term on ionization
c  For icoulm = -1 switch off effect of Coulomb term entirely
c  Initialization also depends on whether Newton iteration or
c  backsubstitution is used
c
      if(icoul1.eq.2.or.icoul1.eq.-1) then
        call zero(dmuc,npar*10)
        call zero(dmucpr,npar)
      else 
c
c  initialize for Coulomb iteration. Switch off second derivatives
c  for iteration passes.
c
        iclder=0
        secder=.false.
        call zero(dmuc,npar*10)
        if(abs(fl-flprev).gt.0.1.or.abs(tl-tlprev).gt.0.1
     *    .or.abs(x-xprev).gt.0.1) then
          call zero(dmucpr,npar)
        end if
      end if
c
      flprev=fl
      tlprev=tl
      xprev=x
      nitdmu=0
c
c  entry point for setting Coulomb derivatives
c
   15 k=-10
      do 25 ia=1,3
      k=k+10
      ext=ex(ia)/tkev
      eea=dmu-ext+dmuc(ia,1)
      if(eea+30.le.0) then
c
c  no ionization
c
        do 16 i=1,10
   16   ea(k+i)=0.e0
        go to 25
c 
      end if
c
      eea=fxp(eea)
c
c  set statistical weights
c
      if(ia.ne.2) then
        eea=eea/2
      else
        eea=2*eea
      end if
c
c  test for departure coefficient for h or he
c
      if(ia.eq.1.and.idpco.gt.0) then
        if(idpco.le.2) bdcoh=frhi/eea
        eea=bdcoh*eea
      else if(ia.ge.2.and.idpco.ge.2) then
        eea=bdcoz*eea
      end if
c
      ea(k+1)=eea
c  first derivatives
      ea(k+2)=(dmuf+dmuc(ia,2))*eea
      ea(k+3)=(amm*ext+dmut+dmuc(ia,3))*eea
      ea(k+4)=(dmux+dmuc(ia,4))*eea
      if(secder) then
c  second derivatives, with Coulomb contribution
        ea(k+5)=(dmuff+dmuc(ia,5))*eea+(dmuf+dmuc(ia,2))*ea(k+2)
        ea(k+6)=(dmuft+dmuc(ia,6))*eea+(dmuf+dmuc(ia,2))*ea(k+3)
        ea(k+7)=(dmufx+dmuc(ia,7))*eea+(dmux+dmuc(ia,4))*ea(k+2)
        ea(k+8)=(dmutt+dmuc(ia,8)-amm2*ext)*eea
     *         +(amm*ext+dmut+dmuc(ia,3))*ea(k+3)
        ea(k+9)=(dmutx+dmuc(ia,9))*eea+(dmux+dmuc(ia,4))*ea(k+3)
        ea(k+10)=(dmuxx+dmuc(ia,10))*eea+(dmux+dmuc(ia,4))*ea(k+4)
      end if
c
   25 continue
c
c  test for setting up for Coulomb iteration
c
      call store(ea,eahat,30)
c
c  in iteration pass, include dmuc from previous case in ea
c  unless ignoring Coulomb effect on ionization
c
      if(iclder.ne.1.and.icoul1.ne.-1.and.icoul1.ne.2) then
        do 26 i=1,3
        is=1+10*(i-1)
   26   ea(is)=eahat(is)*exp(dmucpr(i))
      end if
c
c  entry point for continuing Coulomb iteration by backsubstitution
c  For subsequent iterations, include Coulomb term
c
   27 if(iclder.ne.1.and.nitdmu.gt.0) then
        do 28 i=1,3
        is=1+10*(i-1)
   28   ea(is)=eahat(is)*exp(dmuc(i,1))
      end if
c
c  entry point for continuing Coulomb iteration by Newton iteration
c  for diagnostic or iteration purposes, store current ea
c
   29 call store(ea,east,30)
c
c  set degrees of ionization of H and He
c
      call hheion(ea(1), ea(11), ea(21), xi(1), xi(11), xi(21), secder)
c
c  inclusion of h-
c
   30 if(ihmin.eq.1) then
        call hmnion(tl, ea(1), ehm, xhm, secder)
      else
        call zero(xhm,10)
      end if
c
c  ne and derivatives
c
c  combine he fractions for an in xi(20+.), and for ionization
c  energy in xi(10+.)
c
      exhc=exhe+exhep
c
      if(secder) then
        imx=10
      else
        imx=4
      end if
      do 35 i=1,21,10
   35 call store(xi(i),xii(i),imx)
      xia=xi(1)
      imx=imx+20
c
      do 36 i=21,imx
      i10=i-10
      xio=exhe*xi(i10)+exhc*xi(i)
      xi(i)=2*xi(i)+xi(i10)
   36 xi(i10)=xio
c
      ider2=min0(ider,10)
      if(ihmin.eq.1) then
c
c  in the case with H- there is a risk of negative (unphysical)
c  Ne. Write error message and reset to values without 
c  combine h and h-
c
        dxt1(1)=xi(1)-xhm(1)
        ane(1)=av*(dxt1(1)*x/ah+xi(21)*y/ahe+z*anuh(1))
        if(ane(1).le.0) then
          if(istdpr.gt.0) write(istdpr,1010) xi(1), xhm(1), anuh(1)
          do 37 l=1,ider2
   37     dxt1(l)=xi(l)
          ane(1)=av*(dxt1(1)*x/ah+xi(21)*y/ahe+z*anuh(1))
        else
          do 38 l=2,ider2
   38     dxt1(l)=xi(l)-xhm(l)
        end if
      else
c
c  no H-
c
        do 39 l=1,ider2
   39   dxt1(l)=xi(l)
c
        ane(1)=av*(dxt1(1)*x/ah+xi(21)*y/ahe+z*anuh(1))
c
      end if
c
c  terms needed for x-derivatives
c
      do 40 l=1,4
   40 dxt(l)=av*(dxt1(l)/ah-xi(20+l)/ahe)
c
      ii=4
      do 44 l=1,3
      l1=l+1
      anel=(dxt1(l1)*x/ah+xi(l+21)*y/ahe+z*anuh(l1))*av
      tstl=l.eq.3
      if(tstl) anel=anel+dxt(1)
      ane(l1)=anel
      if(secder) then
        do 42 m=l,3
        m1=m+1
        ii=ii+1
        anelm=(dxt1(ii)*x/ah+xi(20+ii)*y/ahe+z*anuh(ii))*av
        if(tstl) anelm=anelm+dxt(m1)
        if(m.eq.3) anelm=anelm+dxt(l1)
   42   ane(ii)=anelm
      end if
   44 continue
c
c  test that Ne is not zero
c
      if(ane(1).le.0) then
	write(istdou,1110) fl, tl, x, z
	if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,1110) 
     *    fl, tl, x, z
	if(kdgeos.gt.0) then
	  kdgeos=-1
	  return
        else
	  stop 'eqstf'
        end if
      end if
c
c  as auxiliary quantities, to avoid over- and underflow, set
c  derivatives of ne divided by ne
c
      do 46 i=2,iders
   46 aneder(i)=ane(i)/ane(1)
c
c  the density and derivatives (note that rho(2) = dlog rho/dlog f,
c  and so on).
c
      anee=ane(1)
      rho(1)=phi(1)*(crho/anee)
      ii=4
      jj=10
      do 50 l=1,3
      l1=l+1
      rhol=-aneder(l1)/amm
      tstl=l.le.2
      if(tstl) rhol=rhol+phi(l1)
      rho(l1)=rhol
      if(secder) then
        do 48 m=l,3
        ii=ii+1
        m1=m+1
        lm=l+m
        rholm=(aneder(l1)*aneder(m1)-aneder(ii))/amm
        if(tstl.and.m.le.2) rholm=rholm+amm*phi(2+lm)
        if(thrder) then
          do 47 n=m,3
          jj=jj+1
          rhd=-2*aneder(l1)*aneder(m1)*aneder(n+1)/amm
          if(l.lt.3.and.m.lt.3.and.n.lt.3) rhd=rhd+amm2*phi(4+lm+n)
   47     rho(jj)=rhd
        end if
   48   rho(ii)=rholm
      end if
   50 continue
c
c  test for skipping Coulomb terms entirely
c  Also skip this bit for final pass when setting second derivatives
c
      if(icoul1.eq.-1.or.iclsdr.eq.7) go to 60
c
c  set Coulomb terms. Prepare for calling WD f4 routine
c  Note: we currently skip this in the pass for setting derivatives,
c  unless for full ionization, where Coulomb part has not been passed
c  previously
c  If not setting derivatives in f4der is implemented in future,
c  this may require further thought
c
      if(iclder.ne.1.or.cmplio.or.icoul1.eq.2) then
c
        anx=x*av/ah
        any=y*av/ahe
c
        sn(1)=(1.d0-xii(1))*anx
        sn(2)=xii(1)*anx
        sn(3)=(1.d0-xii(11)-xii(21))*any
        sn(4)=xii(11)*any
        sn(5)=xii(21)*any
c
c  to avoid rounding error problems in subtractions in sn(1) and sn(3),
c  reset them from the ea-s, or to zero for full ionization
c
        if(cmplio) then
          sn(1)=0
          sn(3)=0
        else
          sn(1)=anx/(1+east(1))
          sn(3)=any/(1+east(11)*(1+east(21)))
        end if
c
c  number of free electrons. For consistency include only those
c  coming from H and He
c
        sn(6)=xii(1)*anx+(xii(11)+2*xii(21))*any
c
        rhl=log10(rho(1))
        if(idgeos.ge.3.and.istdpr.gt.0) write(istdpr,*) 
     *    'calling f4der at',rhl,tl,sn
c
        call f4der(rhl,tl,sn,nspe,f4,e4,p4,h4,d2f4r2,d2f4t2,
     .             d2f4rt,df4,d2f4i,d2f4f,d2f4x,d2f4r,d2f4t,
     .             dp4i,dp4dr,dp4dt,dp4dx,dh4i,dh4dr,dh4dt,dh4dx,npar)
c
c  test for error in f4der
c
	if(kdgeos.lt.0) return
c
        if(icoul1.ne.2) then
          dmuc(1,1)=-df4(1)/tkergs
          dmuc(2,1)=-df4(2)/tkergs
          dmuc(3,1)=-df4(3)/tkergs
c
c  test for continuing Coulomb iteration
c
          ddmuc1=dmuc(1,1)-dmucpr(1)
          ddmuc2=dmuc(2,1)-dmucpr(2)
          ddmuc3=dmuc(3,1)-dmucpr(3)
          if(idgeos.ge.2.and.istdpr.gt.0) 
     *      write(istdpr,*) ' ddmuc1-3:', ddmuc1, ddmuc2, ddmuc3
        else
          ddmuc1=0
          ddmuc2=0
          ddmuc3=0
        end if
c
      end if
c
c
c  test for convergence
c
      if(.not.cmplio.and.icoul1.ne.2.and.iclder.eq.0.and.
     *   (abs(ddmuc1).ge.epsdmu.or.
     *    abs(ddmuc2).ge.epsdmu.or.
     *    abs(ddmuc3).ge.epsdmu)) then
c
c  test for failure to converge
c
        if(nitdmu.gt.20) then
          write(istdou,1100) fl, tl, x, z,
     *      ddmuc1, ddmuc2, ddmuc3
          if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,1100) 
     *      fl, tl, x, z, ddmuc1, ddmuc2, ddmuc3
	  if(kdgeos.gt.0) then
	    kdgeos=-1
	    return
          else
	    stop 'eqstf'
          end if
        end if
c
c
        nitdmu=nitdmu+1
c
c  store previous value
c
        do 51 i=1,npar
   51   dmucpr(i)=dmuc(i,1)
c
c  test for backsubstitution or Newton iteration
c
        if(iclmit.eq.0) then
          go to 27
        else
c
c  store derivatives in dmuc
c
          do 52 i=1,npar
          dmuc(i,2)=-d2f4r(i)/tkergs
          do 52 j=1,npar
   52     dmuc(i,j+2)=-d2f4f(i,j)/tkergs
c
          if(idgeos.ge.4.and.istdpr.gt.0) 
     *      write(istdpr,*) 'xi(1-3):',xii(1), xii(11), xii(21)
          call clmnew(east,eahat,dmuc,ane(1),x,y,ea,npar,nitdmu)
          go to 29
        end if
c
c
c  test for needing last pass of ionization section, to set
c  derivatives
c
      else if(.not.cmplio.and.icoul1.ne.2.and.iclder.eq.0) then
c
c  set derivatives of delta mu c
c
        do 54 i=1,npar
        do 53 j=2,4
        dmuc(i,j)=0
        do 53 k=1,npar
   53   dmuc(i,j)=dmuc(i,j)-d2f4i(i,k)*xii(10*(k-1)+j)
c
c  add explicit contributions from derivatives of
c  dmuc wrt log rho, log T and X at fixed xi
c  Note: last term in dmuc(i,3) corrects for division by kT.
c
        dmuc(i,2)=dmuc(i,2)-d2f4r(i)*rho(2)
        dmuc(i,3)=dmuc(i,3)-d2f4r(i)*rho(3)-d2f4t(i)+amm*df4(i)
        dmuc(i,4)=dmuc(i,4)-d2f4r(i)*rho(4)-d2f4x(i)
c
        do 54 j=2,4
   54   dmuc(i,j)=dmuc(i,j)/tkergs
c
        iclder=1
c
c  restore secder, unless a later pass will be made to compute
c  second derivatives
c
        secder=.not.nosd
c
        go to 15
c
      else 
c
c  otherwise store pressure and enthalpy from Coulomb effect,
c  transforming derivatives
c
        pcoul(1)=p4
	hcoul(1)=h4
c
c
	call zero(pcoul(2),9)
	call zero(hcoul(2),9)
c
	do 56 j=2,4
        do 56 k=1,npar
        pcoul(j)=pcoul(j)+dp4i(k)*xii(10*(k-1)+j)
   56   hcoul(j)=hcoul(j)+dh4i(k)*xii(10*(k-1)+j)
c
c  add explicit contributions from derivatives of
c  dmuc wrt log rho, log T and X at fixed xi
c
        pcoul(2)=pcoul(2)+dp4dr*rho(2)
        pcoul(3)=pcoul(3)+dp4dr*rho(3)+dp4dt
        pcoul(4)=pcoul(4)+dp4dr*rho(4)+dp4dx
        hcoul(2)=hcoul(2)+dh4dr*rho(2)
        hcoul(3)=hcoul(3)+dh4dr*rho(3)+dh4dt
        hcoul(4)=hcoul(4)+dh4dr*rho(4)+dh4dx
      end if
c
c  test for setting second derivatives of Coulomb terms numerically
c
      if(iclmsd.eq.1.and..not.nosd) then
c
	if(iclsdr.eq.0) then
	  flc=fl
	  tlc=tl
	  xc=x
	  if(icoul1.ne.2) call store(dmuc,dmucc,4*npar)
	  call store(pcoul,pcoulc,4)
	  call store(hcoul,hcoulc,4)
c
	  fl=flc-epssdr
	  iclsdr=1
	  go to 5
c
	else
	  if(icoul1.ne.2) call store(dmuc,dmucp(1,1,iclsdr),4*npar)
	  call store(pcoul,pcoulp(1,iclsdr),4)
	  call store(hcoul,hcoulp(1,iclsdr),4)
c
	  if(iclsdr.le.5) then
	    if(iclsdr.eq.1) then
	      fl=flc+epssdr
	    else if(iclsdr.eq.2) then
	      fl=flc
	      tl=tlc-epssdr
	    else if(iclsdr.eq.3) then
	      tl=tlc+epssdr
	    else if(iclsdr.eq.4) then
	      tl=tlc
	      x=xc-epssdr
	      y=1-x-z
	    else if(iclsdr.eq.5) then
	      x=xc+epssdr
	      y=1-x-z
	    end if
c
	    iclsdr=iclsdr+1
	    go to 5
c
	  else
c
	    x=xc
	    y=1-x-z
c
            if(icoul1.ne.2) call store(dmucc,dmuc,4*npar)
            call store(pcoulc,pcoul,4)
            call store(hcoulc,hcoul,4)
c
c  set second derivatives
c
	    epsdi2=1.d0/(2*epssdr)
            pcoul(5) =(pcoulp(2,2)-pcoulp(2,1))*epsdi2
            hcoul(5) =(hcoulp(2,2)-hcoulp(2,1))*epsdi2
            pcoul(6) =(pcoulp(3,2)-pcoulp(3,1))*epsdi2
            hcoul(6) =(hcoulp(3,2)-hcoulp(3,1))*epsdi2
            pcoul(7) =(pcoulp(4,2)-pcoulp(4,1))*epsdi2
            hcoul(7) =(hcoulp(4,2)-hcoulp(4,1))*epsdi2
            pcoul(8) =(pcoulp(3,4)-pcoulp(3,3))*epsdi2
            hcoul(8) =(hcoulp(3,4)-hcoulp(3,3))*epsdi2
            pcoul(9) =(pcoulp(3,6)-pcoulp(3,5))*epsdi2
            hcoul(9) =(hcoulp(3,6)-hcoulp(3,5))*epsdi2
            pcoul(10)=(pcoulp(4,6)-pcoulp(4,5))*epsdi2
            hcoul(10)=(hcoulp(4,6)-hcoulp(4,5))*epsdi2
c
            if(icoul1.ne.2) then
	      do 58 i=1,npar
              dmuc(i,5) =(dmucp(i,2,2)-dmucp(i,2,1))*epsdi2
              dmuc(i,6) =(dmucp(i,3,2)-dmucp(i,3,1))*epsdi2
              dmuc(i,7) =(dmucp(i,4,2)-dmucp(i,4,1))*epsdi2
              dmuc(i,8) =(dmucp(i,3,4)-dmucp(i,3,3))*epsdi2
              dmuc(i,9) =(dmucp(i,3,6)-dmucp(i,3,5))*epsdi2
   58         dmuc(i,10)=(dmucp(i,4,6)-dmucp(i,4,5))*epsdi2
c
c  make final pass of ionization section to set second derivatives
c  correctly
c
	      iclder=1
	      secder=.true.
	      iclsdr=7
              go to 15
            end if
c
	  end if
c
	end if
c
      end if

c
c  ****************************************************
c
c  end skipping Coulomb terms
c
   60 continue
c
c  start setting total pressure and enthalpy
c
c  as initialization, zero arrays
c
      call zero(pt,20)
      call zero(ht,20)
c
c
c  delta p, delta h and derivatives.
c
c  test for complete ionization
      if(cmplio) go to 80
c
c  delta ne**2
c
c  note: on 3/1/84 definition of dne was changed to
c        (delta ne **2)/av**2 to avoid overflows in calculation
c        of delta p and delta h on univac (with exponents limited
c        to 38). a corresponding change was made in the definition
c        of ca03 in s/r setcns.
c
c
c  reset xi(1) and xi(21) (note that xi(2-10) and xi(22-30) are still
c  consistent).
      xi(1)=-1/(east(1)+1)
      xi(21)=-(2+east(11))/(1+east(11)*(1+east(21)))
c
      do 65 l=1,4
      dxt1(l)=xi(l+20)/ahe-xi(l)/ah
      dxtl=-xi(l+20)*y/ahe-xi(l)*x/ah-anuh(l)*z
      if(l.eq.1) dxtl=dxtl+anh0*z
      dxtl2=ane(l)
      if(l.le.3) go to 64
      dxtl=dxtl+dxt1(1)
      dxtl2=dxtl2+avda
   64 dxt(l)=dxtl
   65 dxt2(l)=dxtl2/av
c
c  ne0/av:
c
      ann=ane(1)/av+(x/ah+2*y/ahe+z*anh0)
      xtt=dxt(1)
      xttiav=xtt/av
      dne(1)=ann*xtt
c
      ii=4
      do 70 l=1,3
      l1=l+1
      dne(l1)=dxt2(l1)*xtt+ann*dxt(l1)
      if(nosd) go to 70
      tstl=l.eq.3
      do 68 m=l,3
      m1=m+1
      ii=ii+1
      dnlm=-xi(20+ii)*y/ahe-xi(ii)*x/ah
      if(tstl) dnlm=dnlm+dxt1(m1)
      if(m.eq.3) dnlm=dnlm+dxt1(l1)
   68 dne(ii)=ane(ii)*xttiav+dxt2(l1)*dxt(m1)+dxt2(m1)*dxt(l1)
     .  +ann*dnlm
   70 continue
c
c  quantities common to delta p and delta h (dph(15-20) is used as
c  intermediate storage).
c
c  Note: here the constant for delta p and delta H is set to C1 = ca03/2
c
      a03=ca03/zf3/2
c
      call zero(dph(15),6)
      dxt1(1)=0
      c1=amm*tkev
      dxt1(2)=c1
      dph(18)=amm2*tkev
      dph(19)=3*c1/zf
      dnee=dne(1)
      a03=a03*rho(1)*rho(1)
      nu=2
      nb=20
      k1=1
   75 do 77 l=2,4
   77 dxt(l-1)=dne(l)+nu*amm*dnee*rho(l)
c
      dxt1(3)=(3*tkev+nb*ak0)/zf
      bmu=tkev+nb*ak0
      dph(k1)=a03*bmu*dnee
      ii=k1+3
      jj=4
      do 79 l=1,3
      l1=l+1
      dph(k1+l)=a03*(bmu*dxt(l)+dnee*dxt1(l))
      if(nosd) go to 79
      do 78 m=l,3
      ii=ii+1
      jj=jj+1
      m1=m+1
      dphlm=bmu*(dne(jj)+nu*amm*(dne(l1)*rho(m1)+dne(m1)*rho(l1)
     .  +dnee*(rho(jj)+nu*amm*rho(m1)*rho(l1))))+dxt(l)*dxt1(m)
     .  +dxt(m)*dxt1(l)+dnee*dph(10+jj)
      if(m.eq.3.and.l.eq.3) dphlm=dphlm+dnee*(12*tkev+2*nb*ak0)/zf2
   78 dph(ii)=a03*dphlm
   79 continue
      if(nu.eq.1) go to 90
      nu=1
      nb=40
      a03=a03/rho(1)
      k1=11
      go to 75
c  complete ionization
   80 call zero(dne,10)
      call zero(dph,20)
      go to 100
   90 continue
      do 95 i=1,iders
      pt(i)=pt(i)+dph(i)
   95 ht(i)=ht(i)+dph(i+10)
c
  100 continue
c
c  *********************************************
c
c  electron pressure and enthalpy
c
      ii=4
      jj=10
      pee=cpe*phi(11)
      pt(1)=pt(1)+pee
      hsst=hst(1)
      hee=che*anee*hsst
      ht(1)=ht(1)+hee
      pe(1)=pee
      he(1)=hee
      do 110 l=1,3
      l1=l+1
      hll=0
      hel=hsst*ane(l1)
      pel=0
      tstl=l.eq.3
      if(tstl) go to 102
      pel=amm*pee*phi(11+l)
      hll=hst(l1)
      hel=hel+anee*hll
      pt(l1)=pt(l1)+pel
  102 pe(l1)=pel
      hel=che*hel
      ht(l1)=ht(l1)+hel
      he(l1)=hel
      if(nosd) go to 110
      do 108 m=l,3
      ii=ii+1
      m1=m+1
      pelm=0
      helm=ane(ii)*hsst+hll*ane(m1)
      tstl=tstl.or.m.eq.3
      if(tstl) go to 104
      lm=l+m
      pelm=amm2*pee*(phi(11+m)*phi(11+l)+phi(12+lm))
      pt(ii)=pt(ii)+pelm
      helm=helm+ane(l1)*hst(m1)+anee*hst(2+lm)
  104 pe(ii)=pelm
      helm=che*helm
      ht(ii)=ht(ii)+helm
      he(ii)=helm
      if(notd) go to 108
      do 106 n=m,3
      jj=jj+1
      helm=0
      if(tstl) go to 106
      helm=ane(n+1)*hst(2+lm)
      if(n.eq.3) go to 105
      helm=helm+ane(l1)*hst(2+m+n)+ane(m1)*hst(2+l+n)
     .  +anee*hst(4+lm+n)
      pt(jj)=pt(jj)+amm3*pee*(phi(11+l)*phi(11+m)*phi(11+n)
     .  +phi(12+lm)*phi(11+n)+phi(12+l+n)*phi(11+m)
     .  +phi(12+m+n)*phi(11+l)+phi(14+l+m+n))
  105 helm=che*helm
      ht(jj)=ht(jj)+helm
  106 he(jj)=helm
  108 continue
  110 continue
c
c  *********************************************
c
c  ionization enthalpy
c
c  (pi is introduced for consistency)
      call zero(pi,10)
      call zero(hi(11),10)
c
      xi(1)=xia
      averg=av*ergev
c  combine h and h-
      do 111 l=1,ider
  111 dxt1(l)=exh*xi(l)-exhm*xhm(l)
      do 112 l=1,4
  112 dxt(l)=dxt1(l)/ah-xi(10+l)/ahe
c
      hi1=averg*(dxt1(1)*x/ah+xi(11)*y/ahe+z*ueh(1))
      hi(1)=hi1
      ht(1)=ht(1)+hi1
      ii=4
      do 115 l=2,4
      dhi=dxt1(l)*x/ah+xi(l+10)*y/ahe+z*ueh(l)
      tstl=l.eq.4
      if(tstl) dhi=dhi+dxt(1)
      dhi=averg*dhi
      ht(l)=ht(l)+dhi
      hi(l)=dhi
      if(nosd) go to 115
      do 114 m=l,4
      ii=ii+1
      dhi=dxt1(ii)*x/ah+xi(10+ii)*y/ahe+z*ueh(ii)
      if(tstl) dhi=dhi+dxt(m)
      if(m.eq.4) dhi=dhi+dxt(l)
      dhi=averg*dhi
      ht(ii)=ht(ii)+dhi
  114 hi(ii)=dhi
  115 continue
c
c  *********************************************
c
c  pressure and enthalpy of heavy particles
c
  120 anh=(x/ah+y/ahe+z/az)*av
      rhh=rho(1)
      phh=anh*rhh*tkergs
      ph(1)=phh
      ph(2)=amm*phh*rho(2)
      drht=1+rho(3)
      ph(3)=amm*phh*drht
      ph(4)=amm*phh*rho(4)+rhh*tkergs*avd1
      do 121 i=1,4
  121 pt(i)=pt(i)+ph(i)
      if(nosd) go to 125
      ph(10)=amm*(phh*rho(10)+rho(4)*(ph(4)+rhh*tkergs*avd1))
      do 122 k=1,3
      k1=k+1
      ph(k+4)=amm*(ph(k1)*rho(2)+phh*rho(k+4))
      if(k.gt.1) ph(k+6)=amm*(ph(k1)*drht+phh*rho(k+6))
  122 continue
      do 123 i=5,10
  123 pt(i)=pt(i)+ph(i)
c
      if(notd) go to 125
      do 124 k=1,3
      k1=k+1
      do 124 l=k,3
      kl=nuu(k,l)
      pt(6+kl)=pt(6+kl)+amm*(ph(kl)*rho(2)+ph(k1)*rho(l+4)+
     .  ph(l+1)*rho(k+4)+phh*rho(kl+6))
      if(k.gt.1) pt(9+kl)=pt(9+kl)+amm*(ph(kl)*drht+ph(k1)*
     .  rho(6+l)+ph(l+1)*rho(6+k)+phh*rho(9+kl))
  124 continue
      pt(20)=pt(20)+amm*(ph(10)*rho(4)+2*ph(4)*rho(10)
     .  +phh*rho(20)+rhh*tkergs*avd1*(amm*rho(4)*rho(4)+rho(10)))
c
  125 hhh=2.5*anh*tkergs
      hh(1)=hhh
      hh(2)=0
      hh(3)=amm*hhh
      hh(4)=2.5*tkergs*avd1
      hh(5)=0
      hh(6)=0
      hh(7)=0
      hh(8)=amm2*hhh
      hh(9)=amm*hh(4)
      hh(10)=0
      call zero(hh(11),10)
      hh(17)=amm*hh(8)
      hh(18)=amm*hh(9)
      ht(1)=ht(1)+hhh
      ht(3)=ht(3)+amm*hhh
      hh4=2.5*tkergs*avd1
      ht(4)=ht(4)+hh4
      if(nosd) go to 130
      ht(8)=ht(8)+amm2*hhh
      ht(9)=ht(9)+amm*hh4
      if(notd) go to 130
      ht(17)=ht(17)+amm3*hhh
      ht(18)=ht(18)+amm2*hh4
c
c  *********************************************
c
c  pressure and enthalpy of radiation (included if iprrad.ne.0)
c
  130 call zero(pr,10)
      call zero(hr,20)
      if(iprrad.eq.0) go to 145
      t4=t*t*t*t
      prr=car*t4/3
      pr(1)=prr
      call zero(pr(2),9)
      pr(3)=4*amm*pr(1)
      pr(8)=4*amm*pr(3)
      pt(1)=pt(1)+prr
      pt(3)=pt(3)+4*amm*prr
      if(nosd) go to 135
      pt(8)=pt(8)+16*amm2*prr
      if(notd) go to 135
      pt(17)=pt(17)+64*amm3*prr
c
  135 hrr=4*pr(1)/rhh
      hr(1)=hrr
      do 136 i=2,4
  136 dxt(i)=-rho(i)
      dxt(3)=4+dxt(3)
      ii=4
      jj=10
      do 140 l=1,3
      l1=l+1
      hr(l1)=amm*hrr*dxt(l1)
      if(nosd) go to 140
      do 138 m=l,3
      ii=ii+1
      m1=m+1
      hr(ii)=amm*(hr(l1)*dxt(m1)-hrr*rho(ii))
      if(notd) go to 138
      do 137 n=m,3
      jj=jj+1
      ln=nuu(l,n)
      mn=nuu(m,n)
  137 hr(jj)=amm*(hr(ii)*dxt(n+1)-amm*hrr*dxt(m1)*rho(ln)-hr(l1)*rho(mn)
     .  -hrr*rho(jj))
  138 continue
  140 continue
      do 142 i=1,ider
  142 ht(i)=ht(i)+hr(i)
c
c  pressure and enthalpy from Coulomb effect
c
  145 if(icoul1.ne.1.and.icoul1.ne.-1) then
        do 150 i=1,ider
        pt(i)=pt(i)+pcoul(i)
  150   ht(i)=ht(i)+hcoul(i)
      end if
c
c  *********************************************
c
c  change to derivatives of log p
c
      ptt=pt(1)
c
c  first divide all derivatives by pt, to avoid over- and
c  underflow
c
      do 15010 i=2,ider
15010 pt(i)=pt(i)/ptt
c
      if(notd) go to 155
      jj=10
      do 152 l=1,3
      do 152 m=l,3
      do 152 n=m,3
      lm=nuu(l,m)
      ln=nuu(l,n)
      mn=nuu(m,n)
      jj=jj+1
  152 pt(jj)=(pt(jj)+(-(pt(lm)*pt(n+1)+pt(ln)*pt(m+1)+pt(mn)*pt(l+1))
     .  +2*pt(m+1)*pt(l+1)*pt(n+1)))/amm
  155 ii=4
      do 158 l=2,4
      ptl=pt(l)/amm
      if(nosd) go to 158
      do 156 m=l,4
      ii=ii+1
  156 pt(ii)=(pt(ii)-pt(l)*pt(m))/amm
  158 pt(l)=ptl
c     if(.not.notd.and.istdpr.gt.0) write(istdpr,15801) pt
15801 format(/' pt:',1p10e12.4/4x,10e12.4)
c
c  cp and dad
c
  160 pf=pt(2)
      hf=ht(2)
      dxtt=ht(3)*pf-hf*pt(3)
      lt=6
      do 165 l=1,3
      lf=l+4
      if(l.gt.1) lt=6+l
  165 dxt(l+1)=ht(lt)*pf+ht(3)*pt(lf)-ht(lf)*pt(3)-hf*pt(lt)
c
      fcpp=1./(amm*t*pf)
      cpp=dxtt*fcpp
      cp(1)=cpp
      if(nosd) go to 173
      do 170 l=2,4
      dcp=dxt(l)*fcpp-cpp*pt(l+3)/pf
      if(l.eq.3) dcp=dcp-cpp*amm
  170 cp(l)=dcp
c
  173 prh=amm*ptt/rhh
      anum=pf*prh-hf
      dad(1)=anum/dxtt
      if(nosd) go to 177
      do 175 l=2,4
      lf=l+3
  175 dad(l)=(prh*(amm*(pt(l)-rho(l))*pf+pt(lf))-ht(lf)-dad(1)*dxt(l)
     .  )/dxtt
c  further derivatives
  177 rhf=rho(2)
      dxtt=pt(3)*rhf-pf*rho(3)
c
c  Gamma1 and derivatives
c
      dxt(1)=rhf-dad(1)*dxtt
      gm1=pf/dxt(1)
      gmm1(1)=gm1
      if(secder) then
        dxt(2)=rho(5)-dad(2)*dxtt
     *    -dad(1)*(pt(6)*rho(2)+pt(3)*rho(5)-pt(5)*rho(3)-pt(2)*rho(6))
        dxt(3)=rho(6)-dad(3)*dxtt
     *    -dad(1)*(pt(8)*rho(2)+pt(3)*rho(6)-pt(6)*rho(3)-pt(2)*rho(8))
        dxt(4)=rho(7)-dad(4)*dxtt
     *    -dad(1)*(pt(9)*rho(2)+pt(3)*rho(7)-pt(7)*rho(3)-pt(2)*rho(9))
        gmm1(2)=(pt(5)-gm1*dxt(2))/dxt(1)
        gmm1(3)=(pt(6)-gm1*dxt(3))/dxt(1)
        gmm1(4)=(pt(7)-gm1*dxt(4))/dxt(1)
      else
        call zero(gmm1(2),3)
      end if
c
      tprh=rhf/dxtt
      trhp=-pf/dxtt
      rhxp=amm*x*(rho(4)-rhf*pt(4)/pf)
c  delta and derivatives
      delta=-1/trhp
      dlt(1)=delta
      lt=6
      if(nosd) go to 190
      do 180 l=2,4
      lf=l+3
      if(l.gt.2) lt=l+5
  180 dlt(l)=(pt(lt)*rhf+pt(3)*rho(lf)-pt(lf)*rho(3)-pf*rho(lt))/pf
     .  -delta*pt(lf)/pf
  190 xii1(1)=xii(1)
      xii1(2)=xii(11)
      xii1(3)=xii(21)
      xii1(4)=anuh(1)/anh0
      return
 1010 format(/' **** error in eqstf. With H-, Ne .le. 0.'/
     *         '      x(H+) =',1pe11.3,'  x(H-) =',e11.3,
     *         ' x(heavies) =',e11.3/
     *         '      Ignore H-')
 1100 format(/' ***** error in s/r eqstf.',
     *  ' Iteration for Coulomb term failed to converge.'/
     *  7x,' log f, log T, X, Z =',4f10.5//
     *  ' last changes in delta mu_c(1-3):'/ 1p3e13.5)
 1110 format(//' ***** Error in eqstf. Ne = 0'/
     *         '       log f, log T, X, Z =',4f12.5/
     *         '       Execution terminated.')
      end
      subroutine phder(fl,tl,phi,hst,nosd,notd)
c
c  computes quantities for Eggleton, Faulker & Flannery approximation
c  to partial degeneracy. On input fl and tl are log(F) and log(T).
c  Returns phi(1-30) and hst(1-10). Here phi(1) is defined such that
c  the density is rho = phi(1)*(crho/ane), where ane is the number
c  of electrons per unit mass, phi(11) is defined such that the
c  electron pressure is pe = cpe*phi(11), and hst(1) such that
c  the electron enthalpy per unit mass is He = che*ane*hst(1).
c  The constants crho, cpe and che are given in common/consts/.
c  phi(21) is related to hst. First, second and third derivatives
c  with respect to log(f) and log(T)
c  of log(phi(1)), log(phi(11)), log(phi(21)) and hst(1) are
c  given in phi(2 - 10), phi(12 - 20), phi(22 - 30) and hst(2 - 10),
c  respectively.
c
      implicit double precision (a-h,o-z)
      logical nosd,notd
      dimension phi(30)
      dimension sig(10),tmn(10),ff(4),gg(4),hst(10)
      common/phdsms/ s0,sf,sg,sff,sfg,sgg,sfff,sffg,sfgg,sggg,
     .  cfg,tf,tg,tff,tfg,tgg,tfff,tffg,tfgg,tggg
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car
      common/eqphcs/ ic,c(48)
      common/ln10/ amm
c  logical equivalence
      equivalence(s0,sig(1)),(cfg,tmn(1))
c     ***********************************
c
      save
c
c  number of sums
c
      imax=10
      if(notd) imax=6
      if(nosd) imax=3
c  f,g
      f=1.d1**fl
      t=1.d1**tl
      ts=ct*t
      wf= sqrt(1+f)
      g=ts*wf
c
c  1/(1+f), 1/(1+g) etc.
c
      vf=1/(1+f)
      vg=1+g
      fdf=vg*g
      fdf=vf*fdf* sqrt(fdf)
      vg=1/vg
      vfg=vf*vg
      ug=g*vg
      uf=f*vf
c
      ug2=ug*vg
      ug3=vg*ug2
c  powers of f and g
      ff(1)=f
      gg(1)=1
      do 10 i=2,ic
      ff(i)=f*ff(i-1)
      gg(i)=g*gg(i-1)
   10 fdf=fdf*vfg
c
c  test on size of f and g
c  icase is set as flag for size, as follows
c   icase = 1: general case
c   icase = 2: f small
c   icase = 3: g small
c   icase = 4: f and g small
c
      icase=1
      if(.not.notd) go to 12
      if(f.lt.1.e-4) icase=icase+1
      if(g.lt.1.e-4) icase=icase+2
c
   12 ic2=ic*ic
c
c  calculate phi* and derivatives
c
      l=1
      anu=1.5e0
      mio=ic
      an32=2.5e0-ic
      kk=-10
c
      l0=1
      do 50 k=1,3
      kk=kk+10
      if(k-2) 18,16,17
c  reset fdf for k=2 and 3
   16 anu=2.5e0
      fdf=g*fdf
      go to 18
   17 mio=mio+1
      fdf=fdf*vf
   18 do 19 i=1,imax
   19 sig(i)=0.e0
      annu=anu-1
c
c  the summation
c
      l=l0
c
c  select case, based on size of f and g, as determined by icase
c
      go to (20,25,30,35), icase
c  the general case
   20 do 23 in=1,ic
      annu=annu+1
      do 23 im=1,ic
c  phimn*(f**(m+1))*(g**n)
      cfg=c(l)*ff(im)*gg(in)
c
      tg=annu*cfg
      tf=im*cfg
      if(nosd) go to 21
c  second derivatives
      tgg=annu*tg
      tfg=im*tg
      tff=im*tf
      if(notd) go to 21
c  third derivatives
      tggg=annu*tgg
      tfgg=im*tgg
      tffg=im*tfg
      tfff=im*tff
c  summing
   21 do 22 i=1,imax
   22 sig(i)=sig(i)+tmn(i)
   23 l=l+1
c  the summation is finished
      if(nosd) go to 45
c
c  the sigma tilde (cf (22.2)) are stored in the corresponding sigma.
c  this is o.k. provided that we go backwards.
c
      s02=s0*s0
      sg2=sg*sg
      sf2=sf*sf
      if(notd) go to 24
c  third derivatives
      s03=s02*s0
      sfff=(sfff*s02-sf*(3*sff*s0-2*sf2))/s03
      sffg=(sffg*s02-sff*sg*s0-2*sf*(sfg*s0-sg*sf))/s03
      sfgg=(sfgg*s02-sgg*sf*s0-2*sg*(sfg*s0-sg*sf))/s03
      sggg=(sggg*s02-sg*(3*sgg*s0-2*sg2))/s03
c  second derivatives
   24 sff=(sff*s0-sf2)/s02
      sfg=(sfg*s0-sf*sg)/s02
      sgg=(sgg*s0-sg2)/s02
      go to 45
c  f is small
   25 do 28 in=1,ic
      annu=annu+1
      do 27 im=1,2
      cfg=c(l)*ff(im)*gg(in)
      sig(im)=sig(im)+cfg
      tg=annu*cfg
      sg=sg+tg
      if(nosd) go to 27
      sgg=sgg+annu*tg
      sfg=sfg+im*tg
   27 l=l+1
   28 l=l+2
c  the summation is finished. set precursors for sigma tilde
      sff=s0*sf
      s0=s0+sf
      sf=s0+sf
      if(nosd) go to 45
      sgg=sgg*s0-sg*sg
      go to 40
c  g is small
   30 ig=1
      do 33 in=1,2
      annu=annu+1
      do 32 im=1,4
      cfg=c(l)*ff(im)*gg(in)
      sig(ig)=sig(ig)+cfg
      tf=im*cfg
      sf=sf+tf
      if(nosd) go to 32
      sff=sff+im*tf
      sfg=sfg+annu*tf
   32 l=l+1
   33 ig=3
c  the summation is finished. set precursors for sigma tilde.
      sgg=s0*sg
      s0=s0+sg
      sg=anu*s0+sg
      if(nosd) go to 45
      sff=sff*s0-sf*sf
      go to 40
c  both f ang g are small
   35 ig=3
c  in this case we must also zero sfg
      sfg=0.e0
      do 38 in=1,2
      annu=annu+1
      do 37 im=1,2
      cfg=c(l)*ff(im)*gg(in)
      sig(im)=sig(im)+cfg
      sig(ig)=sig(ig)+cfg
   37 l=l+1
      ig=5
   38 l=l+2
c  the summation is finished. set precursors for sigma tilde.
      sff=s0*sf
      s0=s0+sf
      sf=s0+sf
      sgg=sg*sfg
      sg=anu*s0+sfg
      if(nosd) go to 45
c  set final values of the sigma tilde.
   40 s02=s0*s0
      sff=sff/s02
      sgg=sgg/s02
      if(f*g.lt.1.00001e-8) go to 42
      sfg=(sfg*s0-sf*sg)/s02
      go to 45
c  approximate expression for sfg (may need fixing up, if f = o(1)
c  or g = o(1))
   42 sfg=f*g*(c(l0+5)-c(l0+1)*c(l0+4)/c(l0))/c(l0)
c
c  phi* and first derivatives
c
   45 phi(kk+1)=fdf*s0
      pht=an32*ug+sg/s0
      phi(kk+3)=pht
      phi(kk+2)=(pht/2-mio)*uf+sf/s0
      if(nosd) go to 50
c
c  second derivatives of phi*.
c
      phtt=an32*ug2+sgg
      phi(kk+6)=phtt
      phi(kk+5)=phtt*uf/2+sfg
      phi(kk+4)=sff+uf*(sfg+vf*(pht/2-mio+f*phtt/4))
c
      if(notd) go to 50
c  third derivatives
      phttt=an32*ug3*(1-g)+sggg
      phi(kk+10)=phttt
      phi(kk+9)=sfgg+uf*phttt/2
      phfft=sffg+uf*(sfgg+vf*(phtt+f*phttt/2)/2)
      phi(kk+8)=phfft
      phi(kk+7)=sfff+uf*(sffg+phfft/2+vf*(1.5*sfg+f*sfgg/4
     .  +vf*((1-f)*(pht/2-mio)+f*phtt/2)))
   50 l0=l0+ic2
c
c  h* and its derivatives (pp 23-25)
c
      do 55 i=2,imax
      ik=20+i
   55 phi(ik)=phi(ik)-phi(i)
c
      hs=phi(21)/phi(1)
      wft1=2*g
      hst(1)=hs+wft1
c
      hf=phi(22)
      ht=phi(23)
      wft2=ts*f/wf
      hst(2)=hs*hf+wft2
      hst(3)=hs*ht+wft1
      if(nosd) go to 58
c  second derivatives
      hff=phi(24)
      hft=phi(25)
      htt=phi(26)
      wft3=uf*(1+f/2)*ts/wf
      hst(4)=hs*(hf*hf+hff)+wft3
      hst(5)=hs*(hf*ht+hft)+wft2
      hst(6)=hs*(ht*ht+htt)+wft1
      if(notd) go to 58
c  third derivatives
      hst(7)=hs*(hf*(hf*hf+3*hff)+phi(27))+uf*vf*(1+f*(2+f)/4)*ts/wf
      hst(8)=hs*(hf*(hf*ht+2*hft)+ht*hff+phi(28))+wft3
      hst(9)=hs*(ht*(ht*hf+2*hft)+hf*htt+phi(29))+wft2
      hst(10)=hs*(ht*(ht*ht+3*htt)+phi(30))+wft1
c  change to derivatives wrt log10 f and log10 t
   58 fct=amm
      do 60 i=2,imax
      if(i.eq.4.or.i.eq.7) fct=fct*amm
   60 hst(i)=hst(i)*fct
      return
      end
      subroutine hheion(eah, eahe, eahep, xih, xihe, xihep, secder)
c
c  given ratios between succesive states of ionization, and
c  derivatives, for H, He and He+ in eah(1-10), eahe(1-10),
c  and eahep(1-10), sets the corresponding degrees of ionization,
c  and derivatives, into xih(1-10), xihe(1-10) and xihep(1-10)
c
      implicit double precision (a-h,o-z)
      logical secder
      dimension eah(10), eahe(10), eahep(10), xih(10), xihe(10),
     *  xihep(10), eahs(10), eahes(10), eahepc(10)
c
      if(secder) then
        iders=10
      else
        iders=4
      end if
c
c  set combined ratio for he+
c
      eea=eahep(1)
      if(eea.gt.0) then
        eea1=eahe(1)
        eahepc(1)=eea*eea1
        ii=4
        if(secder) then
          do 15 l=2,4
          do 15 m=l,4
          ii=ii+1
   15     eahepc(ii)=eea*eahe(ii)+eahep(l)*eahe(m)+eahe(l)*eahep(m)+
     .      eea1*eahep(ii)
        end if
        do 20 i=2,4
   20   eahepc(i)=eahep(i)*eea1+eea*eahe(i)
c
      else
c
        do 22 i=1,iders
   22   eahepc(i)=0
c
      end if
c
c  set x h+, x he+, x he++ and derivatives
c
      dnm=1+eah(1)
      xih(1)=eah(1)/dnm
c
c  hydrogen fraction
c
c  to avoid over- and underflow, divide derivatives
c  of ea by dnm (modification 3/1/84)
c
      do 25 i=2,iders
   25 eahs(i)=eah(i)/dnm
c
      ii=4
      do 30 l=1,3
      l1=l+1
      eal=eahs(l1)
      xih(l1)=eal/dnm
      if(secder) then
        do 28 m=l,3
        m1=m+1
        ii=ii+1
   28   xih(ii)=(eahs(ii)-2*eal*eahs(m1))/dnm
      end if
   30 continue
c
c  helium fractions
c
      dnm=1+eahe(1)+eahepc(1)
c
c  to avoid over- and underflow, divide derivatives
c  of ea by dnm (modification 3/1/84)
c
      do 35 i=2,iders
      eahes(i)=eahe(i)/dnm
   35 eahepc(i)=eahepc(i)/dnm
c
      ii=4
      eeahe=eahe(1)
      eeahep=eahepc(1)
      xihe(1)=eeahe/dnm
      xihep(1)=eeahep/dnm
      do 40 l=2,4
      ealhe=eahes(l)
      ealhep=eahepc(l)
      anmhe=(1+eeahep)*ealhe-eeahe*ealhep
      anmhep=(1+eeahe)*ealhep-eeahep*ealhe
      xihe(l)=anmhe/dnm
      xihep(l)=anmhep/dnm
c
c  second derivatives
c
      if(secder) then
        do 38 m=l,4
        ii=ii+1
        eamhe=eahes(m)
        eamhep=eahepc(m)
c
c  for xi .lt. 1.e-10 zero second derivatives
c
        if(xihe(1).le.1.e-10) then
          xihe(ii)=0
        else
          eamhe=eahes(m)
          eamhep=eahepc(m)
          xihe(ii)=ealhe*eamhep-ealhep*eamhe+
     *      ((1+eeahep)*eahes(ii)-eeahe*eahepc(ii)
     *      -2*anmhe*(eamhe+eamhep))/dnm
        end if
c
        if(xihep(1).le.1.e-10) then
          xihep(ii)=0
        else
          xihep(ii)=ealhep*eamhe-ealhe*eamhep+
     *      ((1+eeahe)*eahepc(ii)-eeahep*eahes(ii)
     *      -2*anmhep*(eamhep+eamhe))/dnm
        end if
   38   continue
      end if
   40 continue
c
      return
      end
      subroutine hmnion(tl, eah, ehm, xhm, secder)
c
c  given ratio between succesive states of ionization, and
c  derivatives, for H in eah(1-10)
c  sets fraction of H-, and derivatives, into xhm(1-10)
c
c  Note: assumes that fraction of H- is always very small.
c
      implicit double precision (a-h, o-z)
      logical secder
      dimension eah(10), ehm(10),xhm(10), eahs(10)
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm
      common/ln10/ amm
      common/dmuder/ dmu,dmuf,dmut,dmux,dmuff,dmuft,dmufx,dmutt,dmutx,
     .  dmuxx,idmu
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      if(idgeos.ge.3.and.istdpr.gt.0) write(istdpr,*) 'Entering hmnion'
      t=10.d0**tl
      tkev=ck1*t
c
      call zero(xhm,10)
c
      ext=exhm/tkev
      eea=ext-dmu
c
c  test for no h-
c
      if(eea.le.-100) then
        if(idgeos.ge.3.and.istdpr.gt.0) write(istdpr,*) 'No H-'
        return
      end if
      eea=0.5*exp(eea)
      dnm=1+eah(1)
      xhm(1)=eea/dnm
c
c  to avoid over- and underflow, divide derivatives
c  of ea by dnm (modification 3/1/84)
c
      if(secder) then
        iders=10
      else
       iders=4
      end if
c
      do 10 i=2,iders
   10 eahs(i)=eah(i)/dnm
c
      ehm(1)=eea
      ehm(2)=-dmuf*eea
      ehm(3)=-(amm*ext+dmut)*eea
      ehm(4)=-dmux*eea
      if(secder) then
        ehm(5)=-dmuff*eea-dmuf*ehm(2)
        ehm(6)=-dmuft*eea-dmuf*ehm(3)
        ehm(7)=-dmufx*eea-dmux*ehm(2)
        ehm(8)=(amm2*ext-dmutt)*eea-(amm*ext+dmut)*ehm(3)
        ehm(9)=-dmutx*eea-dmux*ehm(3)
        ehm(10)=-dmuxx*eea-dmux*ehm(4)
      end if
c
c  derivatives of xh-
c
      ii=4
      do 20 l=1,3
      l1=l+1
      ehml=ehm(l1)
      eal=eahs(l1)
      xhm(l1)=(ehml-eea*eal)/dnm
      if(secder) then
        do 15 m=l,3
        m1=m+1
        ii=ii+1
   15   xhm(ii)=(ehm(ii)-(ehml*eahs(m1)+ehm(m1)*eal+eea*eahs(ii))
     *    +2*eea*eal*eahs(m1))/dnm
      end if
   20 continue
      if(idgeos.ge.3.and.istdpr.gt.0) write(istdpr,*) 'x(H-) =',xhm(1)
      return
      end
      subroutine clmnew(east,eahat,dmuc,ane,x,y,ea,npar,nitdmu)
c
c  finds Newton iteration correction to ionization fractions
c  Input: east: current value of ionization fractions
c         eahat: Coulomb-independent part of ionization fractions
c         dmuc: Coulomb corrections (dmuc(i,1) contains correction,
c               dmuc(i,2) derivative wrt log rho, dmuc(i, 2+j)
c               derivative wrt ionization fraction no. j.
c         ane: electron number density (per unit mass)
c         x and y: hydrogen and helium abundances
c  Returns updated ionization fractions in ea.
c
c  Original version: 10/5/90
c
      implicit double precision (a-h, o-z)
      parameter(nparw=3, nparw1=nparw+1)
      dimension east(30), eahat(30), ea(30), dmuc(npar,1),
     *  rhld(nparw),w(nparw,nparw1)
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm
      common/ln10/ amm
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  set derivatives of log rho wrt ionization fractions
c
      anxfct=x*av/(amm*ane*ah)
      anyfct=y*av/(amm*ane*ahe)
c
      denom=1+east(1)
      denom=denom*denom
      rhld(1)=-anxfct/denom
c
      denom=1+east(11)*(1+east(21))
      denom=denom*denom
      rhld(2)=-anyfct*(1+2*east(21))/denom
      rhld(3)=-anyfct*east(11)*(2+east(11))/denom
c
c  set up linear equations
c
      do 20 i=1,npar
      is=1+10*(i-1)
      eea=eahat(is)*exp(dmuc(i,1))
      w(i,nparw1)=eea-east(is)
c
      do 15 j=1,npar
   15 w(i,j)=-eea*(dmuc(i,j+2)+dmuc(i,2)*rhld(j))
c
   20 w(i,i)=1+w(i,i)
c
      if(idgeos.ge.4.and.istdpr.gt.0) then
        write(istdpr,110)
        do 22 i=1,npar
   22   write(istdpr,120) (w(i,j),j=1,npar),w(i,nparw1)
      end if
c
c  solve linear equations
c
      call leq(w,w(1,nparw1),npar,1,nparw,nparw,err)   
c
c  as pure ad hoc fudge, halve correction if there are convergence
c  problems
c
      if(nitdmu.gt.20) then
        do 24 i=1,npar
   24   w(i,nparw1)=0.5*w(i,nparw1)
      end if
c
c  set updated ea
c
      call store(east,ea,30)
      do 30 i=1,npar
      is=1+10*(i-1)
   30 ea(is)=ea(is)+w(i,nparw1)
c
      if(idgeos.ge.4.and.istdpr.gt.0) then
        write(istdpr,*) 'Corr. to ea', (w(i,nparw1),i=1,npar)
      end if
c
      return
  110 format(' equations in clmnew:')
  120 format(1p5e13.5)
      end
      subroutine hviona(fl,tl,x,y,z,nosd,notd,anu,ue,anur,uer)
c
c  Calculate ionization of heavy elements.
c
c  Modification 17/5/90: include possibility of including full ionization
c     of heavy elements, for ihvz = 4. Note that this involves changing
c     size of arrays in commons /hvredu/ and /xhv/
c
c  Modification 4/6/90: include average degree of ionization for
c     each element in array xhvmn(10) at the beginning of
c     common/hviond/
c
      implicit double precision (a-h,o-z)
      logical nosd,notd
      character name*5
      dimension anu(1),ue(1),anur(1),uer(1)
      dimension eir(10),dr(10),hr(10),gr(10),dmup(10),
     *  xi(29),phi(30),hst(30)
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev
      common/potetc/ chi(125),am(10),iz(10)
      common/hvname/ name(10)
      common/hvabnd/ ab(10),iab
      common /hvredu/ irab,jfirst,ir(10),chir(125),izr(10),abr(10),
     *  amr(10)
      common /hvcntl/ icount,iwrite,dptst0,dptst1
      common/hviond/ xhvmn(10),xhv(125)
      common/eqscnt/ anz0,anze0,ihvz,iprrad,ihmins,igndeg
      common/dmuder/ dmup,idmu
      common/eqdpco/ frhi,bdcoh,bdcoz,idpco
      common/ln10/ amm,amm2
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
      save
c
c  test for restricting set of elements
c
      if(icount.gt.0) go to 10
      icount=1
      if(iwrite.eq.1.and.istdpr.gt.0) then
        write(istdpr,200)
        write(istdpr,205)
        lj=0
        do 3 k=1,10
        izj=iz(k)
        do 2 i=1,izj
        lj=lj+1
    2   xi(i)=chi(lj)
        im=min0(izj,13)
    3   write(istdpr,210) name(k),ab(k),izj,(xi(i),i=1,im)
      end if
c
c  test for restricting set of elements or levels included.
c  for ihvz = 4 keep everything, but for consistency store 
c  full set of information in restricted arrays
c
      if(ihvz.eq.1) then
c
c  reset to restricted set of variables
c
        irab=3
        ir(1)=1
        ir(2)=3
        ir(3)=10
c
c  reset abundances, to get reasonable fit to full case
c
        abr(1)=1.02*ab(1)
        abr(2)=1.35*ab(3)
        abr(3)=3.71*ab(10)
c
c  element with lowest ionization potential (here fe)
c
        jfirst=3
c
c  number of elements treated fully
c
        jiz1=2
c
      else
c
c  use full set of elements
c
        irab=10
        do 76 i=1,iab
        abr(i)=ab(i)
   76   ir(i)=i
c
c  lowest potential is of na
c
        jfirst=5
c
c  number of elements treated fully
c
        if(ihvz.eq.2) then
          jiz1=0
        else if(ihvz.eq.3) then
          jiz1=3
        else if(ihvz.eq.4) then
          jiz1=10
        else
          write(istdou,110) ihvz
          if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,110) ihvz
          stop 'in hvionz'
        end if
      end if
c
c  set possibly reduced set of ionization data
c
   79 j=1
      lj0=0
      ljr=0
      do 5 i=1,10
c
      izi=iz(i)
      if(i.eq.ir(j)) then
c
c  test for inclusion of all levels
c
        if(j.le.jiz1) then
          izrj=izi
        else
          izrj=1
        end if
        izr(j)=izrj
        amr(j)=am(i)
c
c  restore ionization potentials
c
        lj=lj0
        do 4 k=1,izrj
        lj=lj+1
        ljr=ljr+1
    4   chir(ljr)=chi(lj)
        j=j+1
      end if
c
    5 lj0=lj0+izi
c
c  reset anz0 and anze0
c
c  Note: here we also should have reset az, but for consistency with
c  earlier calculations leave it unchanged at this point.
c
      anz0=0
      anze0=0
      do 77 j=1,irab
      anz0=anz0+abr(j)*izr(j)/amr(j)
   77 anze0=anze0+abr(j)*izr(j)
c
      if(iwrite.eq.1.and.istdpr.gt.0) then
c
        write(istdpr,220)
        write(istdpr,205)
        lj=0
        do 7 k=1,irab
        izrj=izr(k)
        do 6 i=1,izrj
        lj=lj+1
    6   xi(i)=chir(lj)
        im=min0(izrj,13)
    7   write(istdpr,210) name(ir(k)),abr(k),izrj,(xi(i),i=1,im)
        write(istdpr,230) anz0, anze0
c
      end if
c
c  change to summed ionization potentials
c
      lj=0
      do 9 k=1,irab
      sum=0
      izj=izr(k)
      do 9 i=1,izj
      lj=lj+1
      sum=sum+chir(lj)
    9 chir(lj)=sum
c
c  set total number of levels for (possibly) reduced set
c
      nlvtot=lj
c
c  ***************************************
c
c  end of initialization section
c
   10 f=1.d1**fl
      t=1.d1**tl
c  test for departure of heavy elements
      if(idpco.lt.2) bdcoz=1
c  k*t, in ev
      tk=ck1*t
c
c  test whether phder has already been called
c
      if(idmu.eq.1) go to 15
c
      call phder(fl,tl,phi,hst,nosd,notd)
c
c
c  psi and derivatives
      wf= sqrt(1+f)
      psi=2*wf+log(f/((1+wf)*(1+wf)))
      psif=amm*wf
      psff=amm2*f/2/wf
c
      zf=x+2*y+anze0*z
      zf2=zf*zf
      zf3=zf*zf2
c
      ak0=ckh*zf2
c  a0**3*ne/(k*t)
      aa=caa*phi(1)/(tk*zf3)
c
c  delta mu and derivatives
c
      bmu=tk+20*ak0
      dmu=aa*bmu
      dmps=dmu-psi
      dmup(1)=dmps
      ref=phi(2)
      ret=phi(3)
      dmuf=dmu*ref*amm
      dmut=aa*(bmu*ret-20*ak0)*amm
      dmux=aa*(3*tk+20*ak0)/zf
      dmup(2)=dmuf-psif
      dmup(3)=dmut
      dmup(4)=dmux
c
      if(nosd) go to 18
      dmup(5)=(dmuf*ref+dmu*phi(4)*amm)*amm-psff
      dmup(6)=(dmut*ref+dmu*phi(5)*amm)*amm
      dmup(7)=dmux*ref*amm
      dmup(8)=aa*(20*ak0*(1-2*ret)+bmu*(phi(6)+ret*ret))*amm2
      dmup(9)=aa*(ret*(3*tk+20*ak0)-20*ak0)*amm/zf
      dmup(10)=aa*(12*tk+40*ak0)/zf2
c
      go to 18
   15 dmps=dmup(1)
c
   18 idmx=10
      if(nosd) idmx=4
      do 19 i=1,10
      anu(i)=0.e0
      ue(i)=0.e0
      anur(i)=0.e0
   19 uer(i)=0.e0
c
c
      lj=0
      tki=1.d0/tk
      bdcozl=log(bdcoz)
c
c  ***************************************
c
c  start loop over elements
c
      do 50 j=1,irab
      izoj=iz(ir(j))
      izj=izr(j)
c
c  skip elements with zero abundance
c
      if(abr(j).le.0) then
        anur(j)=0
        uer(j)=0
        xhvmn(j)=0
        lj=lj+izj
        go to 50
      end if
c
c  set limit for no ionization
c
      dptstn=dptst1
      if(j.eq.jfirst) dptstn=dptst0
c
c  set exponents phi in saha equations into xi
c
      lj0=lj
      do 20 i=1,izj
      lj=lj+1
      xhv(lj)=0
      dxp=i*dmps-chir(lj)*tki+bdcozl
   20 xi(i)=dxp
c
c  -----------------
c
c  test for complete or no ionization
c
      if(izj-1) 21,21,25
c
c  only one level
c
   21 phm=xi(1)
      if(phm) 22,22,23
c
   22 imax=0
      phm=0
c  test for no ionization
      if(xi(1).lt.-dptstn) go to 50
c
      go to 34
c
   23 imax=1
c  test for complete ionization
      if(xi(1).gt.dptst1) go to 29
      go to 34
c
c  more than one level
c
   25 ii=izj+1
      xi(ii)=0
      imax=1
      phm=xi(1)
c
c  set phm to largest xi
c
      do 26 i=2,ii
      if(xi(i).le.phm) go to 26
      imax=i
      phm=xi(i)
   26 continue
c
      if(imax.ne.izj) go to 30
c
c  test for complete ionization
c
      izj1=izj-1
      dphm=phm-xi(1)
      if(izj1.eq.1) go to 28
      do 27 i=1,izj1
   27 dphm=min(dphm,phm-xi(i))
   28 if(dphm.le.dptst1) go to 34
c
c  complete ionization
c
   29 xhv(lj)=1
      fct1=abr(j)/amr(j)
      anur(j)=fct1*izj
      uer(j)=fct1*chir(lj)
      anu(1)=anu(1)+anur(j)
      ue(1)=ue(1)+uer(j)
      xhvmn(j)=1
      go to 50
c
   30 if(imax.ne.ii) go to 34
c
c  test for no ionization
c
      do 31 i=1,izj
      if(xi(i).gt.-dptstn) go to 34
   31 continue
c
c  no ionization. skip element completely
c
      xhvmn(j)=0
c
      go to 50
c
c  ************************************
c
c  general case
c
   34 do 35 i=1,idmx
      dr(i)=0.e0
      hr(i)=0.e0
   35 gr(i)=0.e0
      if(phm.le.dptst1) dr(1)=omega(0,izoj)*exp(-phm)
      do 40 i=1,izj
      lji=lj0+i
      cchi=chir(lji)
      dxp=xi(i)-phm
      if(dxp.lt.-dptst1) go to 40
      dxp=omega(i,izoj)*exp(dxp)
      xhv(lji)=dxp
      eir(1)=1
      do 36 k=2,4
   36 eir(k)=i*dmup(k)
      eir(3)=eir(3)+cchi*amm*tki
      ii=4
      if(nosd) go to 38
      do 37 k=2,4
      do 37 l=k,4
      ii=ii+1
   37 eir(ii)=eir(k)*eir(l)+i*dmup(ii)
      eir(8)=eir(8)-amm2*cchi*tki
c
   38 do 39 k=1,idmx
      eeir=dxp*eir(k)
      dr(k)=dr(k)+eeir
      hr(k)=hr(k)+i*eeir
   39 gr(k)=gr(k)+cchi*eeir
   40 continue
c
      dr1i=1.d0/dr(1)
c  scale xhv
      do 42 i=1,izj
      lji=lj0+i
      xhv(lji)=dr1i*xhv(lji)
   42 continue
c
      fct1=abr(j)/(amr(j)*dr(1))
      hrdr=hr(1)*dr1i
      grdr=gr(1)*dr1i
      anur(j)=fct1*hr(1)
      uer(j)=fct1*gr(1)
      anu(1)=anu(1)+anur(j)
      ue(1)=ue(1)+uer(j)
c
      xhvmn(j)=hr(1)/(izj*dr(1))
c
c  derivatives
      ii=4
      do 48 k=2,4
      anu(k)=anu(k)+fct1*(hr(k)-dr(k)*hrdr)
      ue(k)=ue(k)+fct1*(gr(k)-dr(k)*grdr)
      if(nosd) go to 48
      do 45 l=k,4
      ii=ii+1
      anu(ii)=anu(ii)+fct1*(hr(ii)-(dr(k)*hr(l)+dr(l)*hr(k)+(dr(ii)
     .  -2*dr(k)*dr(l)*dr1i)*hr(1))*dr1i)
   45 ue(ii)=ue(ii)+fct1*(gr(ii)-(dr(k)*gr(l)+dr(l)*gr(k)+(dr(ii)
     .  -2*dr(k)*dr(l)*dr1i)*gr(1))*dr1i)
   48 continue
   50 continue
c
      return
  110 format(//' **** error in s/r hviona. ihvz =',i5,
     *  ' not allowed')
  200 format(///' Original heavy element data.'/)
  205 format(' Name, abundance, Z, ionization potentials:'/)
  210 format(1x,a4,f8.5,i4,13f8.2)
  220 format(///' Heavy element data after resetting:'/)
  230 format(//' now anz0 =',f10.5,'    anze0 =',f10.5//)
  250 format(/' xhv:'/(1p10e12.4))
      end
      double precision function omega(i,iz)
c  calculates statistical weight of i-th ionization stage of element
c  with number iz
      implicit double precision (a-h,o-z)
      common/hvomeg/ iom(26),iom1(20)
      common/hvomcl/ iomfll
c
      save
c
      if(i.le.1.and.iz.ge.19) go to 20
      if(i.eq.iz) go to 15
      omega=iom(iz-i)
      return
c
c  statistical weight for fully ionized atom.
c  before 5/1/84 was always set to 15, for some bizarre reason.
c  for transition introduce flag iomfll so that iomfll = 0 corresponds
c  to old situation, and iomfll .ne. 0 gives correct value omega = 1.
c
   15 omega=1
      if(iomfll.eq.0) omega=15
      return
   20 omega=iom1(2*(iz-19)+i+1)
      return
      end
      subroutine seteqs
c
c  dummy subroutine included for compatibility with dog equation  of
c  state programmes.
c
      return
      end
      block data bleqst
c  initialize data for equation of state
c
c  Modified 22/5/90 to include all levels of ionization of Ar and Fe.
c  Previously only the first 15 levels were included, the remainder
c  being forced to be unionized
c
      implicit double precision (a-h,o-z)
      character name*5
      common/eqscnt/ anh0,anhe0,ihvz,iprrad,ihmin,igndeg
      common/hvabnd/ ab(10),iab
      common/eqdpco/ frhi,bdcoh,bdcoz,idpco
      common/eqphcs/ ic,c(48)
      common/potetc/ chi(125),am(10),iz(10)
      common/hvname/ name(10)
      common /hvcntl/ icnthv,iwrthv,dptst0,dptst1
      common/hvomeg/ iom(26),iom1(20)
      common/hvomcl/ iomfll
      common/ccoulm/ epsdmu, icoulm, iclmit, iclmsd, nitdmu, epssdr
c
      data anh0,anhe0,ihvz,iprrad,ihmin /0.5,6.,1,1,0/
      data ab,iab /0.2254,0.0549,0.4987,0.0335,0.00197,0.0436,
     *  0.00403,0.0565,0.00180,0.0795,10/
      data frhi,bdcoh,bdcoz,idpco /1.,1.,1.,0/
c  coefficients for s/r phder.
      data ic,c/4,2.315472,7.128660,7.504998,2.665350,7.837752,23.507934
     .,23.311317,7.987465,9.215560,26.834068,25.082745,8.020509,3.693280
     .,10.333176,9.168960,2.668248,2.315472,6.748104,6.564912,2.132280,
     .7.837752,21.439740,19.080088,5.478100,9.215560,23.551504,19.015888
     .,4.679944,3.693280,8.859868,6.500712,1.334124,1.157736,3.770676,
     .4.015224,1.402284,8.283420,26.184486,28.211372,10.310306,14.755480
     .,45.031658,46.909420,16.633242,7.386560,22.159680,22.438048,
     .7.664928/
      data iz/6,7,8,10,11,12,13,14,18,26/
      data name/'   C','   N','   O','  Ne','  Na',
     .  '  Mg','  Al','  Si','  Ar','  Fe'/
      data am/12.00,14.01,16.00,20.17,22.99,24.31,26.98,28.08,
     .  39.94,55.84/
      data chi/11.26,24.38,47.86,64.48,391.99,489.84,
     .  14.54,29.60,47.43,77.45,97.86,551.92,666.83,
     .  13.61,35.15,54.93,77.39,113.87,138.08,739.11,871.12,
     .  21.56,41.07,63.5,97.16,126.4,157.91,207.3,239.,1196.,1360.,
     .  5.14,47.29,71.65,98.88,138.60,172.36,208.44,264.15,299.78,
     .  1465.,1646.,
     .  7.64,15.03,80.12,109.29,141.23,186.86,225.31,265.96,327.90,
     .  367.36,1761.2,2085.,
     .  5.98,18.82,28.44,119.96,153.77,190.42,241.93,285.13,330.1,
     .  398.5,441.9,2085.5,2299.,
     .  8.15,16.34,33.46,45.13,166.73,205.11,246.41,303.87,
     .  351.83,401.3,476.0,523.2,2436.,2666.,
     .  15.75,27.62,40.90,59.79,75.0,91.3,124.,143.46,421.,480.,
     .  539.5,621.1,688.5,755.5,854.4,918.,4121.,4426.,
     .  7.90,16.18,30.64,56.,79.,105.,133.,151.,235.,262.,290.,321.,
     .  355.,390.,457.,489.,1266.,1358.,1456.,1582.,1689.,1799.,1950.,
     .  2045.,8828.,9278./
      data icnthv,iwrthv /0,0/
c
c  limiting arguments in exponential in test for full or no
c  ionization in s/r hviona.
c  dptst0 is used for first level and element to ionize
c  dptst1 is used for remaining levels
c
      data dptst0,dptst1 /85.,19./
c
c  quantities for calculating statistical weights in
c  function omega.
c
      data iom/2,1,2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,10,21,28,25,6,25,
     .  28,21/
      data iom1/2,1,1,2,10,15,21,28,28,25,7,6,6,7,25,30,28,21,21,10/
c  flag for statistical weight omega for fully ionized atom.
c  before 5/1/84 omega was always set to 15, for some bizarre reason.
c  for transition introduce flag iomfll so that iomfll = 0 corresponds
c  to old situation, and iomfll .ne. 0 gives correct value omega = 1.
      data iomfll /1/
c
c  controls for Coulomb effect.
c
      data epsdmu, icoulm, iclmit, iclmsd, epssdr
     *    / 1.e-12,   0,       1,      0 , 1.e-3  /
      end
