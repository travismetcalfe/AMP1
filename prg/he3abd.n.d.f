      subroutine he3abd(fl,tl,x,y,z,age,x3,anu,noder)
c  calculates he 3  abundance at t=age (in seconds), assuming
c  that reaction rates etc have been constant from time 0.
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
      implicit double precision (a-h,o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
c  Note: engenr.n.d.incl replaced by engenr.nnz.d.incl, 5/6/02
      include 'engenr.nnz.d.incl'
c
      logical nosd,noder,norct
      dimension x3(*),aa(4),an(4)
      common/eqstd/ dumm(14),rho(10)
      common/rnrcnt/ irnrat, krnrat, irn12, irn13, irn14, 
     *  irn15c, irn15o, irn16, irn17
      common/rnratd/ al(10,krnrmx),norct
      common/consts/ av,ah,ahe
      common/rcncns/ a(knucst),b(knucst),s(knucst,isnuc),
     *  zz(knucst),zz12(2,knucst),zsmh,acno,awght(knucwg),q(knucq), 
     *  knucpm, kqvals, kawght
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
c  calculate xhe3
      ca1=x*ahe/(y*ah*al(1,3))
      ca1=ca1*al(1,1)*al(1,2)*ca1
c
c  test for expansion
c
      if(ca1.gt.1.e-4) go to 5
      aa1=ca1*(1-0.5d0*ca1)
      aaa=1+aa1
      go to 8
c
    5 aaa= sqrt(1+2*ca1)
      aa1=aaa-1
c
    8 x3eq=ahe3*y*al(1,3)*aa1/(2*ahe*al(1,2))
      alp=(aaa+1)/(2*ca1)
      anu=rho(1)*al(1,3)*y*aaa/ahe
      tau=age*anu/2
c
c  test for equilibrium
c
      if(tau.gt.50.or.age.lt.0) then
        coef=1
c
      else
c
        ak=tanh(tau)
        coef=ak*(1+2*alp)/(1+alp*(1+ak))
      end if
c
      x3(1)=coef*x3eq
c
      if(idghe3.eq.1.and.istdpr.gt.0) then
        write(istdpr,120) fl,tl,rho(1),al(1,2),al(1,3),anu,tau,x3eq,coef
      end if
      if(noder) return
c
c  derivatives of x3
c
      a2=aaa/alp
      do 20 i=2,4
   20 x3(i)=x3eq*amm*(al(i,3)-al(i,2)+ca1*(al(i,1)+al(i,2)
     .  -2*al(i,3))/a2)
      x3(4)=x3(4)+x3eq*(2*ca1*(1/x+1/y)/a2-1/y)
      if(coef.gt.0.99999) return
c
      do 22 i=2,4
   22 aa(i)=ca1*amm*(al(i,1)+al(i,2)-2*al(i,3))/aaa
      aa(4)=aa(4)+2*ca1*(1/x+1/y)/aaa
      do 24 i=2,4
   24 an(i)=anu*(amm*(rho(i)+al(i,3))+aa(i)/aaa)
      an(4)=an(4)-anu/y
c
      dnom=(1+alp*(1+ak))**2
      dy1=age*(1+alp)*(1+2*alp)*(1+ak)*(1-ak)/dnom/2
      dy2=ak*(1-ak)*alp*alp/dnom
      do 28 i=2,4
   28 x3(i)=coef*x3(i)+x3eq*(dy1*an(i)+dy2*aa(i))
      return
c
c  no reactions. zero he3 abundance and derivatives
c
   50 x3(1)=0
      if(noder) return
      do 55 i=2,4
   55 x3(i)=0
      return
  120 format(' in he3abd, fl, tl, rho(1),al(1,2),al(1,3), ',
     *  ' nu, tau, x3eq, coef ='/ 2f10.5,1p7e13.5)
      end
