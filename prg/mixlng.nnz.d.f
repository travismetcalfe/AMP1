      subroutine mixlng(t,hp,dr,dac,ddac,ptpg)
c
c  routine to calculate the actual temperature gradient dac
c  and its derivatives ddac(1-6) from mixing length theory.
c
c  on input t, hp and dr are temperature, pressure scale height
c  and radiative gradient.
c  note: s/r eqstf must be called before call of mixlng.
c
c  correction 4/6/1985: before this date, albd was always
c  one. in fact albd depends on etac and phc. the value one
c  is appropriate for the default values of etac and phc
c  (which were in fact used in all significant cases).
c
c  notice that albd corresponds to 27/4*lambda, where lambda
c  is as defined in thesis. Thus c1 corresponds to lambda(thesis)
c
c  addition 13/6/1985: if etac .lt. 0, use a simple exponential
c  approximation to delta - delta ad.
c
c  Modified 23/6/89: install same fudge as used in Cambridge (1986)
c  paper
c
c  Modified 6/4/01: set dyda according to approximations for small
c  and large a. Before that, the general expression was used in
c  all cases.
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
c  Modified 23/8/96, to allow call of MJM routines
c
c  Modified 18/11/98, adding ptpg to the argument list, for case
c  with turbulent pressure
c
      implicit double precision (a-h,o-z)
      dimension ddac(1),is(5),dda(5)
      common/ln10/ amm
      common/mxlcn/ cc1,cc2,etac,phc,alfa,tprfct,iconcs,imxlng,iturpr
      common/eqstd/ dum(14),rho,drho(19),ht(20),p,dp(19),cp,dcp(3),
     .  dad,ddad(3),dlt,ddlt(3),gm1
      common/rhcn/ a1,a2,a3,a4,z,nvar,ifwrt,irhtst
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
      common/opcdat/ akk,dak(3)
      common/cqhopf/ d2qhpf,tgrfct,dtgrfc(5)
      common/cnvout/ ddrad,rr,a,ddacad,amach,adrr(5),addr(5),ada(5),
     *  pturb,y,dyda
      common/noiter/ iter, ntime
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      data etacp,phcp /-1.,-1./,is/2,3,5,1,4/
c
      idgmxl=0
c
	if(idgrhs.eq.2.and.istdpr.gt.0) 
     *    write(istdpr,*) 'cc1,cc2,etac,phc,alfa', 
     *      cc1,cc2,etac,phc,alfa
c
c  test for using MJM routines
c
      if(iconcs.gt.0) then
        call mxlmjm(t,hp,dr,dac,ddac,ptpg)
	return
      end if
c
c  test for simple exponential fit to ddacad
c
      if(etac.lt.0) go to 85
c
      if(etac.eq.etacp.and.phc.eq.phcp) go to 5
c
      albd=4.5*etac*sqrt(phc)
      c1=4*albd/27d0
      c2=c1*(1.5d0-c1)
      c3=c1*(1-c1)
      c4=c3*c3*c3
      c5=9/(4*albd)
      ex=1.d0/3.d0
c
      etacp=etac
      phcp=phc
c
c  g and k
    5 g=p*(1.d0+ptpg)/(hp*rho)
      ak=cc2*t*t*t/(akk*rho)
c  the mixing length
      aml=alfa*hp
c  curly r
	if(idgrhs.eq.2.and.istdpr.gt.0) write(istdpr,*) 
     *  ' alfa,  hp, aml, rho, cp, ak, g, dlt =',
     *    alfa,  hp, aml, rho, cp, ak, g, dlt
      rr=aml*aml*rho*cp/ak
      rr=(rr/hp)*rr*g*dlt
	if(idgrhs.eq.2.and.istdpr.gt.0) write(istdpr,*) ' rr =',rr
c  delta r - delta ad
      ddrad=dr-dad
      if(ddrad.le.0) ddrad=1.e-10
c  set a
      a=2/ sqrt(rr*ddrad)/etac
	if(idgrhs.eq.2.and.istdpr.gt.0) write(istdpr,*) ' a =',a
c
c  test for asymptotic regime
c
   15 if(a.le.1.e-5) go to 20
      if(a.gt.15.) go to 30
c
c  general y
c
      gm=c1*(1.5d0/(a*a)+c2)
      x=a*(gm+ sqrt(gm*gm+c4))**ex
      y=a*(x/a-c3*a/x-c1)
      dyda=(1-y*(y+2*a))/(a*a+y*(2*a+y/c1))
      go to 40
c
c  small a
c
   20 x=(a/c5)**ex
      y=x*(1-x*x/3d0)
      dyda=ex*(x/a)*(1.d0-2.d0*x*x/3.d0)
      go to 40
c
c  large a
c
   30 x=1/(a*a)
      y=(1-x+(2d0-c5)*x*x)/a
      dyda=-x*(1.d0-3.d0*x+5.d0*(2.d0-c5)*x*x)
      dm=c5*x*x
      go to 50
c
   40 continue
	if(idgrhs.eq.2.and.istdpr.gt.0) write(istdpr,*) ' y =',y
      dm=1-y*(y+a)
   50 continue
c  grad - grad ad
      ddacad=y*(y+a)*ddrad
      dac=ddacad+dad
	if(idgrhs.eq.2.and.istdpr.gt.0) 
     *    write(istdpr,*) ' ddrad, ddacad, dac =',
     *    ddrad, ddacad, dac 
c
c  diagnostics if idgrhs = 1
c
      if(idgrhs.eq.1.and.istdpr.gt.0) write(istdpr,59901) 
     *   t,hp,dr,p,rho,rr,a,y,ddacad,dac
59901 format(' t,hp,dr,p,rho =',1p5e15.7/
     *  ' rr,a,y,ddacad,dac=',1p5e15.7)
c
c  the mach number
c
   60 amach=dm*ddrad
      amach=amach*amach
      x=aml*g*rho/p
      amach=(phc*amach/(4*rr))**ex*dlt*x*x/gm1
      amach= sqrt(amach)
c
c  derivatives
c
      ddy=(2*y+a)*dyda+y
   70 do 80 i=1,5
      j=is(i)
      if(i.ge.4) go to 72
      ddr=dak(i)+dp(i)
      dlrr=(2*dcp(i)/cp+ddlt(i)/dlt)/amm+2*dak(i)+drho(i)+3*dp(i)
      dadi=ddad(i)
      go to 75
   72 dadi=0
      ddr=0
      dlrr=0
   75 if(j.eq.4) ddr=1
      if(j.eq.1) dlrr=4
      if(j.ne.3) go to 77
      ddr=ddr-4
      dlrr=dlrr-6
c
c  set d nablaR/d tj, including contribution from atmospheric factor
c
   77 ddr=dr*(amm*ddr+dtgrfc(j)/tgrfct)
      dlrr=amm*dlrr
c
      da=-a*(dlrr*ddrad+ddr-dadi)/2
      ddac(j)=dadi+ddy*da+y*(y+a)*(ddr-dadi)
      adrr(j)=rr*dlrr
      addr(j)=ddr
      ada(j)=da/ddrad
      dda(j)=da
   80 continue
      if(idgrhs.eq.1.and.istdpr.gt.0) write(istdpr,80091) dac
80091 format(' dac =',1pe15.7)
      if(idgmxl.ge.1.and.istdpr.gt.0) 
     *  write(istdpr,'(a,1p4e14.6)') 'dac, ddac =',
     *  dac,(ddac(i),i=1,3)
      return
c
c  **************************************************
c
c  simple exponential fit
c
   85 a=alfa**(-0.8)
      xx=log10(p)-5.1
      zz=phc*xx
      ddacad=a*exp(-min(zz,50.d0))
      ddx=(1-zz)*ddacad
      ddaca1=xx*ddacad
      ddaca2=1-dad/dr
      ddacad=ddaca1*ddaca2
c
      dac=ddacad+dad
c
c  derivatives
c
      do 87 i=1,5
      j=is(i)
c
      if(i.le.3) then
        dldr=dak(i)+dp(i)
        dadi=ddad(i)
        dddaci=ddx*dp(i)
      else
        dldr=0
        dadi=0
        dddaci=0
      end if
c
      if(j.eq.4) dldr=1
      if(j.eq.3) dldr=dldr-4
c
   87 ddac(j)=dadi+dddaci*ddaca2-ddaca1*(dadi-amm*dad*dldr)/dr
c
      ddrad=dr-dad
      amach=0
      return
      end
