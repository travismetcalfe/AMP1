      subroutine eqstf(fl,tl,x,y,z,nosd,notd)
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
c  modification 15/9/86: save statement added in all subroutines.
c     (needed in f77, apparently)
c
c  modification 3/1/88: dimension of dxt1 increased to 20.
c
c  Modification 6/8/88: second derivative of He ionization fractions
c     corrected. Also in commom/eqsout/ is now stored ea 
c     before rescaling, for test of derivatives.
c
c  Modified 18/5/90, adding array gmm1(4) in common/ eqstd/ containing
c  Gamma1 and derivatives.
c
c  Modification 21/5/90: include possibility of including complete
c     treatment of ionization of heavy elements, for ihvz = 4.
c     Note that this involves changing size of arrays in 
c     commons /hvredu/ and /xhv/
c
c  Modified 22/5/90 to include all levels of ionization of Ar and Fe.
c     Previously only the first 15 levels were included, the remainder
c     being forced to be unionized.
c
c  Modification 4/6/90: include average degree of ionization for
c     each element in array xhvmn(10) at the beginning of
c     common/hviond/
c  Modification 30/9/91: default value of iomfll changed from 0 to 1
c
c  Modified 5/3/03: Include common/cdgphs/ kdgeos, kdgopc, kdgeng.
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
      implicit double precision (a-h,o-z)
      logical tstl,nosd,notd,cmplio, secder
      dimension phi(30),hst(10),ex(3),ea(30),xi(30),dxt(4),dxt1(20),
     .  dxt2(4),
     .  anuh(10),ueh(10),anuhr(23),uehr(23),
     *  ehm(10),aneder(10)
      common/eqscnt/ anh0,anhe0,ihvz,iprrad,ihmin,igndeg
      common/ln10/ amm,amm2,amm3
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm
      common/dmuder/ dmu,dmuf,dmut,dmux,dmuff,dmuft,dmufx,dmutt,dmutx,
     .  dmuxx,idmu
      common/eqsout/ east(30),xii(30),dne(10),dph(20),he(20),pe(10),
     *  hi(20),pi(10),hh(20),ph(10),hr(20),pr(10)
      common/eqstd/ xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp,gmm1(4)
      common/eqstfl/ flcm
      common/eqdpco/ frhi,bdcoh,bdcoz,idpco
      common/xhminc/ xhm(10)
      common/hviond/ xhvmn(10),xhv(125)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  commons setting diagnostics and defining version
c
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
      common/cversn/ ivrbcs,ivrrhs,ivreos,ivropc,ivreng
c
      common/cdgphs/ kdgeos, kdgopc, kdgeng
c
      external bleqst, blstio
c
      equivalence(ex(1),exh)
      data kdgeos, kdgopc, kdgeng /0, 0, 0/ 
c
      save
c
      fxp(a)=exp(min(85.d0,max(a,-85.d0)))
      nuu(l,n)=((l-1)*(8-l))/2+n-l+5
c
c  *****************************************************************
c
c  set equation of state version to 0
c
      ivreos = 0
c
      secder=.not.nosd
c
c  *****************************************************************
c
c  set y
      y=1-x-z
c  store fl in common
      flcm=fl
c  number of derivatives
      ider=20
      if(notd) ider=10
      if(nosd) ider=4
      iders=min0(ider,10)
c
      f=1.d1**fl
      t=1.d1**tl
c
      call phder(fl,tl,phi,hst,nosd,notd)
c..      if(fl.le.-19) write(6,*) 'Point 1 in eqstf'
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
      if(idgeos.ge.1.and.istdpr.gt.0) then
	write(istdpr,*) 'Enter eqstf with fl, tl, x, y, z =',
     *     fl, tl, x, y, z
        write(istdpr,*) 'After call of phder, phi(1) =', phi(1)
      end if
c  k*t, in ev
      tk=ck1*t
c  a0**3*ne/(k*t)
      aa=caa*phi(1)/(tk*zf3)
c
c  delta mu and derivatives
c
      bmu=tk+20*ak0
      dmu=aa*bmu
c
      ref=phi(2)
      ret=phi(3)
      dmuf=dmu*ref*amm
      dmut=aa*(bmu*ret-20*ak0)*amm
      dmux=aa*(3*tk+20*ak0)/zf
c
      if(nosd) go to 10
      dmuff=(dmuf*ref+dmu*phi(4)*amm)*amm-psff
      dmuft=(dmut*ref+dmu*phi(5)*amm)*amm
      dmufx=dmux*ref*amm
      dmutt=aa*(20*ak0*(1-2*ret)+bmu*(phi(6)+ret*ret))*amm2
      dmutx=aa*(ret*(3*tk+20*ak0)-20*ak0)*amm/zf
      dmuxx=aa*(12*tk+40*ak0)/zf2
   10 dmu=dmu-psi
      dmuf=dmuf-psif
      idmu=1
c..      if(idgeos.eq.-2) write(6,*) '#D# eos(1)',fl, tl, dmu, dmu-54.4/tk
c  test for complete ionization of h and he
      cmplio=(dmu-54.4/tk).gt.19
c  test for use of departure coefficient
      if(idpco .ge.1) cmplio=cmplio.and.frhi.gt.1.d6
      if(cmplio) go to 31
c
c  e h, e he, e he+ and derivatives
c
      k=-10
      do 25 ia=1,3
      k=k+10
      ext=ex(ia)/tk
      eea=dmu-ext
      if(eea+30) 15,15,21
c  no ionization
   15 do 16 i=1,10
   16 ea(k+i)=0.d0
      go to 25
c
   21 eea=fxp(eea)
c  test for departure coefficient for he
      if(ia.ge.2.and.idpco.ge.2) eea=bdcoz*eea
      if(ia-2)  22,23,22
   22 eea=eea/2
c  test for departure coefficient for h
      if(idpco.le.0.or.ia.ne.1) go to 24
      if(idpco.le.2) bdcoh=frhi/eea
      eea=bdcoh*eea
      go to 24
   23 eea=2*eea
   24 ea(k+1)=eea
c  first derivatives
      ea(k+2)=dmuf*eea
      ea(k+3)=(amm*ext+dmut)*eea
      ea(k+4)=dmux*eea
      if(nosd) go to 25
c  second derivatives
      ea(k+5)=dmuff*eea+dmuf*ea(k+2)
      ea(k+6)=dmuft*eea+dmuf*ea(k+3)
      ea(k+7)=dmufx*eea+dmux*ea(k+2)
      ea(k+8)=(dmutt-amm2*ext)*eea+(amm*ext+dmut)*ea(k+3)
      ea(k+9)=dmutx*eea+dmux*ea(k+3)
      ea(k+10)=dmuxx*eea+dmux*ea(k+4)
   25 continue
c  reset e he+
      eea=ea(21)
      if(eea.eq.0) go to 35
      eea1=ea(11)
      ea(21)=eea*eea1
      ii=24
      if(nosd) go to 27
      do 26 l=1,3
      l1=l+21
      do 26 m=l,3
      ii=ii+1
      m1=m+21
   26 ea(ii)=eea*ea(ii-10)+ea(l1)*ea(m1-10)+ea(l1-10)*ea(m1)+
     .  eea1*ea(ii)
   27 do 28 i=22,24
   28 ea(i)=ea(i)*eea1+eea*ea(i-10)
   30 continue
      go to 35
c
c  x h+, x he+, x he++ and derivatives
c
c  complete ionization
   31 xi(1) =1
      xi(11)=0
      xi(21)=1
      do 33 i=2,22,10
      jj=i+8
      do 33 j=i,jj
   33 xi(j)=0
      go to 50
c
c  partial ionization
c
c  for diagnostic purposes, store current ea
c
   35 call store(ea,east,30)
c
	dnm=1+ea(1)
      xi(1)=ea(1)/dnm
c
c  to avoid over- and underflow, divide derivatives
c  of ea by dnm (modification 3/1/84)
c
      do 36 i=2,iders
   36 ea(i)=ea(i)/dnm
c
      ii=4
      do 40 l=1,3
      l1=l+1
      eal=ea(l1)
      xi(l1)=eal/dnm
      if(nosd) go to 40
      do 38 m=l,3
      m1=m+1
      ii=ii+1
   38 xi(ii)=(ea(ii)-2*eal*ea(m1))/dnm
   40 continue
c..      if(fl.le.-19) write(6,*) 'Point 2 in eqstf'
c
      dnm=1+ea(11)+ea(21)
c
c  to avoid over- and underflow, divide derivatives
c  of ea by dnm (modification 3/1/84)
c
      do 42 k=10,20,10
      do 42 i=2,iders
   42 ea(i+k)=ea(i+k)/dnm
c
      k1=11
      k2=21
   45 kd=k1-k2
      ii=k1+3
      ea1=ea(k1)
      ea2=ea(k2)
      xi(k1)=ea1/dnm
      do 48 l=1,3
      l1=k1+l
      eal1=ea(l1)
      eal2=ea(k2+l)
      anm=(1+ea2)*eal1-ea1*eal2
      xi(l1)=anm/dnm
c
c  second derivatives
c
      if(nosd) go to 48
      do 46 m=l,3
      ii=ii+1
c
c  for xi .lt. 1.e-10 zero second derivatives
c
      if(xi(k1).gt.1.e-10) go to 45500
      xi(ii)=0
      go to 46
c
45500 eam1=ea(k1+m)
      eam2=ea(k2+m)
c
c  old, incorrect form for second derivative. changes 6/8/88
c
c..      xi(ii)=(eal1*eam2+(1+ea2)*ea(ii)-eal2*eam1-ea1*ea(ii-kd)
c..     .  -2*anm*(eam1+eam2))/dnm
c
      xi(ii)=eal1*eam2-eal2*eam1+((1+ea2)*ea(ii)-ea1*ea(ii-kd)
     .  -2*anm*(eam1+eam2))/dnm
   46 continue
   48 continue
      if(k1.eq.21) go to 50
      k1=21
      k2=11
      go to 45
c
c  ***************************************
c
c  ionization of heavy elements
c
   50 if(ihvz) 51,51,52
   51 anuh(1)=anh0
      call zero(anuh(2),9)
      call zero(ueh,10)
      go to 53
   52 call hviona(fl,tl,x,y,z,nosd,notd,anuh,ueh,anuhr,uehr)
c..      if(fl.le.-19) write(6,*) 'Point 3 in eqstf'
c
c  inclusion of h-
c
   53 call zero(xhm,10)
      if(ihmin.ne.1) go to 54
      ext=exhm/tk
      eea=ext-dmu
c  test for no h-
      if(eea.le.-100) go to 54
      eea=0.5*exp(eea)
      dnm=1+ea(1)
      xhm(1)=eea/dnm
c
      ehm(1)=eea
      ehm(2)=-dmuf*eea
      ehm(3)=-(amm*ext+dmut)*eea
      ehm(4)=-dmux*eea
      if(nosd) go to 535
      ehm(5)=-dmuff*eea-dmuf*ehm(2)
      ehm(6)=-dmuft*eea-dmuf*ehm(3)
      ehm(7)=-dmufx*eea-dmux*ehm(2)
      ehm(8)=(amm2*ext-dmutt)*eea-(amm*ext+dmut)*ehm(3)
      ehm(9)=-dmutx*eea-dmux*ehm(3)
      ehm(10)=-dmuxx*eea-dmux*ehm(4)
c  derivatives of xh-
  535 dnm2=dnm*dnm
      dnm3=dnm2*dnm
      ii=4
      do 537 l=1,3
      l1=l+1
      ehml=ehm(l1)
      eal=ea(l1)
      xhm(l1)=(ehml-eea*eal)/dnm
      if(nosd) go to 537
      do 536 m=l,3
      m1=m+1
      ii=ii+1
  536 xhm(ii)=(ehm(ii)-(ehml*ea(m1)+ehm(m1)*eal+eea*ea(ii))
     *  +2*eea*eal*ea(m1))/dnm
  537 continue
c
c  ne and derivatives
c
c  combine he fractions for an in xi(20+.), and for ionization
c  energy in xi(10+.)
   54 exhc=exhe+exhep
c
      imx=10
      if(nosd) imx=4
      do 55 i=1,21,10
   55 call store(xi(i),xii(i),imx)
      xia=xi(1)
      imx=imx+20
c
      do 56 i=21,imx
      i10=i-10
      xio=exhe*xi(i10)+exhc*xi(i)
      xi(i)=2*xi(i)+xi(i10)
   56 xi(i10)=xio
c  combine h and h-
      ider2=min0(ider,10)
      do 561 l=1,ider2
  561 dxt1(l)=xi(l)-xhm(l)
c  terms needed for x-derivatives
      do 57 l=1,4
   57 dxt(l)=av*(dxt1(l)/ah-xi(20+l)/ahe)
c
      if(idgeos.ge.1.and.istdpr.gt.0) then
	write(istdpr,*) 
     *    'Setting ane(1) with av, dxt1(1), x, ah, xi(21), ',
     *    'y, ahe, z, anuh(1) ='
        write(istdpr,*) av, dxt1(1), x, ah, xi(21), y, ahe, z, anuh(1)
      end if
      ane(1)=av*(dxt1(1)*x/ah+xi(21)*y/ahe+z*anuh(1))
c
      ii=4
      do 60 l=1,3
      l1=l+1
      anel=(dxt1(l1)*x/ah+xi(l+21)*y/ahe+z*anuh(l1))*av
      tstl=l.eq.3
      if(tstl) anel=anel+dxt(1)
      ane(l1)=anel
      if(nosd) go to 60
      do 58 m=l,3
      m1=m+1
      ii=ii+1
      anelm=(dxt1(ii)*x/ah+xi(20+ii)*y/ahe+z*anuh(ii))*av
      if(tstl) anelm=anelm+dxt(m1)
      if(m.eq.3) anelm=anelm+dxt(l1)
   58 ane(ii)=anelm
   60 continue
c..      if(fl.le.-19) write(6,*) 'Point 4 in eqstf'
c
c  test that Ne is not zero
c
      if(ane(1).le.0) then
	write(istdou,300) fl, tl, x, z
	if(istdpr.gt.0.and.istdou.ne.istdpr) 
     *    write(istdpr,300) fl, tl, x, z
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
      do 60100 i=2,iders
60100 aneder(i)=ane(i)/ane(1)
c
c  the density and derivatives (note that rho(2) = dlog rho/dlog f,
c  and so on).
c
      anee=ane(1)
c
      rho(1)=phi(1)*(crho/anee)
      if(idgeos.eq.-2) write(98,*) 
     *  fl, tl, dmu, dmu-54.4/tk, phi(1), rho(1)
      if(idgeos.ge.1.and.istdpr.gt.0) 
     *  write(istdpr,*) 'Setting rho with phi(1), crho, anee =',
     *  phi(1), crho, anee
      ii=4
      jj=10
      do 63 l=1,3
      l1=l+1
      rhol=-aneder(l1)/amm
      tstl=l.le.2
      if(tstl) rhol=rhol+phi(l1)
      rho(l1)=rhol
      if(nosd) go to 63
      do 62 m=l,3
      ii=ii+1
      m1=m+1
      lm=l+m
      rholm=(aneder(l1)*aneder(m1)-aneder(ii))/amm
      if(tstl.and.m.le.2) rholm=rholm+amm*phi(2+lm)
      if(notd) go to 62
      do 61 n=m,3
      jj=jj+1
      rhd=-2*aneder(l1)*aneder(m1)*aneder(n+1)/amm
      if(l.lt.3.and.m.lt.3.and.n.lt.3) rhd=rhd+amm2*phi(4+lm+n)
   61 rho(jj)=rhd
   62 rho(ii)=rholm
   63 continue
c
c  ****************************************************
c
c  start setting total pressure and enthalpy
c
c  as initialization, zero arrays
c
c..      if(fl.le.-19) write(6,*) 'Point 5 in eqstf'
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
      xi(1)=-1/(ea(1)+1)
      xi(21)=-(2+ea(11))/(1+ea(11)+ea(21))
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
c..      if(fl.le.-19) write(6,*) 'Point 6 in eqstf'
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
      c1=amm*tk
      dxt1(2)=c1
      dph(18)=amm2*tk
      dph(19)=3*c1/zf
      dnee=dne(1)
      a03=a03*rho(1)*rho(1)
      nu=2
      nb=20
      k1=1
   75 do 77 l=2,4
   77 dxt(l-1)=dne(l)+nu*amm*dnee*rho(l)
c
      dxt1(3)=(3*tk+nb*ak0)/zf
      bmu=tk+nb*ak0
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
      if(m.eq.3.and.l.eq.3) dphlm=dphlm+dnee*(12*tk+2*nb*ak0)/zf2
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
c..      if(fl.le.-19) write(6,*) 'Point 7 in eqstf'
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
c..      if(fl.le.-19) write(6,*) 'Point 8 in eqstf'
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
c..      if(fl.le.-19) write(6,*) 'Point 9 in eqstf'
c  k*t, in ergs
      tk=ck2*t
      rhh=rho(1)
      phh=anh*rhh*tk
      ph(1)=phh
      ph(2)=amm*phh*rho(2)
      drht=1+rho(3)
      ph(3)=amm*phh*drht
      ph(4)=amm*phh*rho(4)+rhh*tk*avd1
      do 121 i=1,4
  121 pt(i)=pt(i)+ph(i)
      if(nosd) go to 125
      ph(10)=amm*(phh*rho(10)+rho(4)*(ph(4)+rhh*tk*avd1))
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
     .  +phh*rho(20)+rhh*tk*avd1*(amm*rho(4)*rho(4)+rho(10)))
c
  125 hhh=2.5*anh*tk
      hh(1)=hhh
      hh(2)=0
      hh(3)=amm*hhh
      hh(4)=2.5*tk*avd1
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
      hh4=2.5*tk*avd1
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
c..      if(fl.le.-19) write(6,*) 'Point 10 in eqstf'
      if(iprrad.eq.0) go to 150
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
      do 145 i=1,ider
  145 ht(i)=ht(i)+hr(i)
c     if(notd) go to 150
c
c  *********************************************
c
c  change to derivatives of log p
c
  150 ptt=pt(1)
c..      if(fl.le.-19) write(6,*) 'Point 11 in eqstf'
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
c..      if(fl.le.-19) write(6,*) 'Point 13 in eqstf'
      dxt(1)=rhf-dad(1)*dxtt
c
c  test for dominance of radiation pressure
c  (temporary, inelegant solution)
c
      if(dxt(1).eq.0) then
	write(istdou,400) fl,tl,x,z
	if(istdpr.gt.0.and.istdou.ne.istdpr) write(istdou,400) fl,tl,x,z
	gm1=1.33333333333333333d0
	gmm1(1)=gm1
	if(secder) call zero(gmm1(2),3)
c
      else
c
        gm1=pf/dxt(1)
        gmm1(1)=gm1
        if(secder) then
          dxt(2)=rho(5)-dad(2)*dxtt
     *      -dad(1)*(pt(6)*rho(2)+pt(3)*rho(5)
     *              -pt(5)*rho(3)-pt(2)*rho(6))
          dxt(3)=rho(6)-dad(3)*dxtt
     *      -dad(1)*(pt(8)*rho(2)+pt(3)*rho(6)
     *              -pt(6)*rho(3)-pt(2)*rho(8))
          dxt(4)=rho(7)-dad(4)*dxtt
     *      -dad(1)*(pt(9)*rho(2)+pt(3)*rho(7)
     *              -pt(7)*rho(3)-pt(2)*rho(9))
          gmm1(2)=(pt(5)-gm1*dxt(2))/dxt(1)
          gmm1(3)=(pt(6)-gm1*dxt(3))/dxt(1)
          gmm1(4)=(pt(7)-gm1*dxt(4))/dxt(1)
        else
          call zero(gmm1(2),3)
        end if
      end if
c
      tprh=rhf/dxtt
      trhp=-pf/dxtt
      rhxp=amm*x*(rho(4)-rhf*pt(4)/pf)
c  delta and derivatives
c..      if(fl.lt.-20) write(6,*) 'trph=',trph
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
      if(idgeos.le.-2) then
c..	write(6,*) '#D# phder:', fl, tl, xii1
	write(97,*) fl,tl,xii1
	if(idgeos.eq.-3) write(96,*) fl, tl, xii1, ane(1), rho(1),
     *    gmm1(1),(pt(i),i=1,10),(rho(i),i=1,10),(xhvmn(i),i=1,10)
      end if
      return
  300 format(//' ***** Error in eqstf. Ne = 0'/
     *         '       log f, log T, X, Z =',4f12.5/
     *         '       Execution terminated.')
  400 format(//' ***** Warning in eqstf.',
     *         '  Radiation pressure dominates Gamma1.'/
     *         '       log f, log T, X, Z =',4f12.5)
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
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm
      common/eqphcs/ ic,c(48)
      common/ln10/ amm,amm2,amm3
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
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
      anu=1.5d0
      mio=ic
      an32=2.5d0-ic
      kk=-10
c
      l0=1
      do 50 k=1,3
      kk=kk+10
      if(k-2) 18,16,17
c  reset fdf for k=2 and 3
   16 anu=2.5d0
      fdf=g*fdf
      go to 18
   17 mio=mio+1
      fdf=fdf*vf
   18 do 19 i=1,imax
   19 sig(i)=0.d0
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
      sfg=0.d0
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
c..      if(idgeos.eq.-2) then
c..	write(6,*) '#D# phder:', fl, tl, f, g, icase, phi(1)
c..	write(97,*) fl,tl,f,g,icase,(phi(i),i=1,10)
c..      end if
      return
      end
      subroutine hviona(fl,tl,x,y,z,nosd,notd,anu,ue,anur,uer)
c
c  Calculate ionization of heavy elements.
c
c  Modification 17/5/90: include possibility of including 
c     full ionization of heavy elements, for ihvz = 4. 
c     Note that this involves changing
c     size of arrays in commons /hvredu/ and /xhv/
c
c  Modification 4/6/90: include average degree of ionization for
c     each element in array xhvmn(10) at the beginning of
c     common/hviond/
c
c  Modified 19/10/06: Add equivalence of dmu and dmup(1). Otherwise,
c  information in common/dmuder/ is apparently not passed properly
c  into dmup for ionization calculation.
c
      implicit double precision (a-h,o-z)
      logical nosd,notd
      character name*5
      dimension anu(1),ue(1),anur(1),uer(1)
      dimension eir(10),dr(10),hr(10),gr(10),dmup(10),
     *  xi(29),phi(30),hst(30)
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm
      common/potetc/ chi(125),am(10),iz(10)
      common/hvname/ name(10)
      common/hvabnd/ ab(10),iab
      common /hvredu/ irab,jfirst,ir(10),chir(125),izr(10),abr(10),
     *  amr(10)
      common /hvcntl/ icount,iwrite,dptst0,dptst1
      common/hviond/ xhvmn(10),xhv(125)
      common/eqscnt/ anz0,anze0,ihvz,iprrad,ihmin,igndeg
      common/dmuder/ dmu,dmuf,dmut,dmux,dmuff,dmuft,dmufx,dmutt,dmutx,
     .  dmuxx,idmu
      common/eqdpco/ frhi,bdcoh,bdcoz,idpco
      common/ln10/ amm,amm2,amm3
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  Add equivalence 19/10/06
c
      equivalence(dmu,dmup(1))
c
      save
c
c  test for resetting diagnostics
c
      if(istdpr.le.0) iwrite=0
c
c  test for restricting set of elements
c
      if(icount.gt.0) go to 10
      icount=1
      if(iwrite.eq.1) then
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
          if(istdpr.gt.0.and.istdou.ne.istdpr) write(istdpr,110) ihvz
	  if(kdgeos.gt.0) then
	    kdgeos=-1
	    return
          else
            stop 'Stop in eqstf'
          end if
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
      anz0=0
      anze0=0
      do 77 j=1,irab
      anz0=anz0+abr(j)*izr(j)/amr(j)
   77 anze0=anze0+abr(j)*izr(j)
c
      if(iwrite.eq.1) then
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
        write(istdpr,230) anz0
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
      anu(i)=0.d0
      ue(i)=0.d0
      anur(i)=0.d0
   19 uer(i)=0.d0
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
c
      izoj=iz(ir(j))
      izj=izr(j)
c..      write(6,*) 'j, izj, izoj', j, izj, izoj
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
c..      write(6,*) 'complete ionization at lj =',lj
      fct1=abr(j)/amr(j)
      anur(j)=fct1*izj
      uer(j)=fct1*chir(lj)
      anu(1)=anu(1)+anur(j)
      ue(1)=ue(1)+uer(j)
      xhvmn(j)=1
      if(idgeos.eq.-2) write(6,*) '#D# hviona(1)',j, fl, tl, dphm,
     *  xhvmn(j)
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
      dr(i)=0.d0
      hr(i)=0.d0
   35 gr(i)=0.d0
      if(phm.le.dptst1) dr(1)=omega(0,izoj)*exp(-phm)
      do 40 i=1,izj
      lji=lj0+i
      cchi=chir(lji)
      dxp=xi(i)-phm
      if(dxp.lt.-dptst1) go to 40
      dxp=omega(i,izoj)*exp(dxp)
c..      print *,' j,i,dxp',j,i,dxp
      xhv(lji)=dxp
      if(abs(tl-7).le.0.1.and.idgeos.eq.-2) write(6,*) '#D# hvion(0)',
     *  j,i,cchi,xi(i)-phm,omega(i,izoj),dxp
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
c..      write(6,*) 'At i, j, lji =',i,j,lji,'  xhv =',xhv(lji)
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
      if(idgeos.eq.-2) write(6,*) '#D# hviona(2)',j, fl, tl, 
     *  (xhv(lji0+i),i=1,izj), xhvmn(j)
c
c..      print *,' j,anur(j),anu(1)',j,anur(j),anu(1)
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
  230 format(//' now anz0 =',f10.5//)
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
      data anh0,anhe0,ihvz,iprrad,ihmin,igndeg /0.5,6.,1,1,0,0/
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
c  Change default to 1, 30/9/91
      data iomfll /1/
c
c  controls for Coulomb effect.
c
      data epsdmu, icoulm, iclmit, iclmsd
     *    / 1.e-12,   0,       1,      0 /
      end
