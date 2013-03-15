      subroutine he3abd(fl,tl,x,y,z,age,x3,anu,noder)
c  calculates he 3  abundance at t=age (in seconds), assuming
c  that reaction rates etc have been constant from time 0.
c  Initial abundance is taken from xzer3 in common/compsz/.
c  if noder=.false. first derivatives of abundance are
c  calculated and stored in x3(2-4).
c  note that s/r eqstf must be called before call of he3abd.
c
c  modified 19/12/1984 to make expansion for small ca1.
c
c  modified 29/12/1984, to include option of age .lt. 0.
c  this is used to flag for setting equilibrium he3 abundance
c  and its derivatives.
c
c  Modified 30/9/2005, to allow non-zero initial abundance.
c
      implicit double precision (a-h,o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
c  Note: engenr.n.d.incl replaced by engenr.bz.d.incl, 30/9/05
      include 'engenr.bz.d.incl'
c
      logical nosd,noder,norct
      dimension x3(*)
      dimension aa(4), an(4), beta(4), psi1(4), psi2(4)
      common/eqstd/ dumm(14),rho(10)
      common/rnrcnt/ irnrat, krnrat, irn12, irn13, irn14, 
     *  irn15c, irn15o, irn16, irn17
      common/rnratd/ al(10,krnrmx),norct
      common/consts/ av,ah,ahe
      common/rcncns/ a(knucst),b(knucst),s(knucst,isnuc),
     *  zz(knucst),zz12(2,knucst),zsmh,acno,awght(knucwg),q(knucq), 
     *  knucpm, kqvals, kawght
      common/compsz/ xzerh, yzer, xzer3, xrz12, xrz13, xrz14, xrz16
      common/ln10/ amm
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      ahe3 = awght(1)
c
c  reaction rates
c
      irnrat = 1
      nosd=.true.
      tt=10.d0**tl
      call rnrate(fl,tl,x,y,z,nosd)
c
c  test for no reactions
c
      if(norct.or.al(1,2).eq.0.or.al(1,3).eq.0) go to 50
c
c  calculate xhe3
c
      ca1=x*ahe/(y*ah*al(1,3))
      ca1=ca1*al(1,1)*al(1,2)*ca1
c
c  test for expansion
c
      if(ca1.le.1.e-4) then
        aa1=ca1*(1-0.5d0*ca1)
        aaa=1+aa1
      else
c
        aaa= sqrt(1+2*ca1)
        aa1=aaa-1
c
      end if
c
      x3eq=ahe3*y*al(1,3)*aa1/(2*ahe*al(1,2))
      anu=rho(1)*al(1,3)*y*aaa/ahe
      tau=age*anu
c
c  test for equilibrium
c
      if(tau.gt.50.or.age.lt.0) then
        coef=1
c
      else
c
        beta(1)=(aaa+1)/aa1
        ak=exp(-tau)
        ak1=1-ak
        psi1(1)=x3eq*beta(1)*ak1+(1+beta(1)*ak)*xzer3
        psi2(1)=x3eq*(beta(1)+ak)+       ak1*xzer3
        coef=psi1(1)/psi2(1)
      end if
c
      x3(1)=coef*x3eq
c
      if(idghe3.eq.1.and.istdpr.gt.0) 
     *  write(istdpr,120) 
     *  fl,tl,rho(1),al(1,2),al(1,3),anu,tau,x3eq,coef,ca1
c
c  test for derivatives
c
      if(noder) return
c
      alp=(aaa+1)/(2*ca1)
      a2=aaa/alp
      do 20 i=2,4
   20 x3(i)=x3eq*amm*(al(i,3)-al(i,2)
     *     +ca1*(al(i,1)+al(i,2)-2*al(i,3))/a2)
      x3(4)=x3(4)+x3eq*(2*ca1*(1/x+1/y)/a2-1/y)
c
c  now x3(.) contains derivatives of equilibrium abundance
c
      if(coef.eq.1) return
c
      do 22 i=2,4
   22 aa(i)=ca1*amm*(al(i,1)+al(i,2)-2*al(i,3))/aaa
      aa(4)=aa(4)+2*ca1*(1/x+1/y)/aaa
      do 24 i=2,4
   24 an(i)=anu*(amm*(rho(i)+al(i,3))+aa(i)/aaa)
      an(4)=an(4)-anu/y
c
      betfct=-2.d0/(aa1*aa1)
      do 30 i=2,4
      beta(i)=betfct*aa(i)
      psi1(i)=beta(1)*ak1*x3(i)+(ak1*x3eq+ak*xzer3)*beta(i)
     *       -beta(1)*(xzer3-x3eq)*ak*age*an(i)
   30 psi2(i)=(beta(1)+ak)*x3(i)+x3eq*beta(i)
     *       +(xzer3-x3eq)*ak*age*an(i)
c
      do 35 i=2,4
   35 x3(i)=x3(1)*(psi1(i)/psi1(1)-psi2(i)/psi2(1)+x3(i)/x3eq)
      return
c
c  no reactions. Set he3 abundance to initial value and zero derivatives
c
   50 x3(1)=xzer3
      if(noder) return
      do 55 i=2,4
   55 x3(i)=0
      return
  120 format(' in he3abd, fl, tl, rho(1),al(1,2),al(1,3), ',
     *  ' nu, tau, x3eq, coef, ca1 ='/ 2f10.5,1p8e13.5)
      end
