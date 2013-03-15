      subroutine he3abc(x,y,z,age,x30,alrhmn,x3,x3eq,anu,idalam,noder)
c
c  calculates he 3  abundance at t=age (in seconds), assuming
c  that reaction rates etc have been constant from time 0,
c  where abundance was x30.
c  note that s/r eqstf must be called before call of he3abd.
c  also alrhmn(1,i) must contain the relevant values of
c  rho*alam, for reaction rates. When setting derivatives,
c  alrhmn(k,i) must contain the derivatives of log10(rho*alrhmn(1,i)),
c  stored in usual way.
c
c  when age .lt. 0, x3 is returned as equilibrium abundance.
c  otherwise, equilibrium abundance is returned in x3eq.
c
c  This version of the routine he3abd has been set up to compute
c  He3 abundance in convective core. However, it may form
c  model for later revision of he3abd, for general initial
c  abundance.
c
c  Modified 21/1/96, to set also derivatives of He3 abundance,
c  when noder is false
c
c  Modified 1/10/05, correcting errors found in connection with
c  development of he3abd.bz.d.f
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
      dimension alrhmn(idalam,*),x3(*)
      dimension aa(4), an(4), beta(4), psi1(4), psi2(4)
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
c  test for no reactions
c
      if(alrhmn(1,2).eq.0.or.alrhmn(1,3).eq.0) then
	x3(1)=x30
	return
      end if
c
c  calculate xhe3
c
      ca1=x*ahe/(y*ah*alrhmn(1,3))
      ca1=ca1*alrhmn(1,1)*alrhmn(1,2)*ca1
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
      x3eq=ahe3*y*alrhmn(1,3)*aa1/(2*ahe*alrhmn(1,2))
      anu=alrhmn(1,3)*y*aaa/ahe
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
        psi1(1)=x3eq*beta(1)*ak1+(1+beta(1)*ak)*x30
        psi2(1)=x3eq*(beta(1)+ak)+       ak1*x30
        coef=psi1(1)/psi2(1)
      end if
c
      x3(1)=coef*x3eq
c
      if(idghe3.eq.1.and.istdpr.gt.0) 
     *  write(istdpr,120) anu,tau,x30,x3eq,coef
c
c  test for derivatives
c
      if(noder) return
c
      alp=(aaa+1)/(2*ca1)
      a2=aaa/alp
      do 20 i=2,4
   20 x3(i)=x3eq*amm*(alrhmn(i,3)-alrhmn(i,2)
     *     +ca1*(alrhmn(i,1)+alrhmn(i,2)-2*alrhmn(i,3))/a2)
      x3(4)=x3(4)+x3eq*(2*ca1*(1/x+1/y)/a2-1/y)
c
c  now x3(.) contains derivatives of equilibrium abundance
c
      if(coef.eq.1) return
c
      do 22 i=2,4
   22 aa(i)=ca1*amm*(alrhmn(i,1)+alrhmn(i,2)-2*alrhmn(i,3))/aaa
      aa(4)=aa(4)+2*ca1*(1/x+1/y)/aaa
      do 24 i=2,4
   24 an(i)=anu*(amm*alrhmn(i,3)+aa(i)/aaa)
      an(4)=an(4)-anu/y
c
      betfct=-2.d0/(aa1*aa1)
      do 30 i=2,4
      beta(i)=betfct*aa(i)
      psi1(i)=beta(1)*ak1*x3(i)+(ak1*x3eq+ak*x30)*beta(i)
     *       -beta(1)*(x30-x3eq)*ak*age*an(i)
   30 psi2(i)=(beta(1)+ak)*x3(i)+x3eq*beta(i)
     *       +beta(1)*(x30-x3eq)*ak*age*an(i)
c
      do 35 i=2,4
   35 x3(i)=x3(1)*(psi1(i)/psi1(1)-psi2(i)/psi2(1)+x3(i)/x3eq)
      return
  120 format(' in he3abc, nu, tau, X30, X3eq, coef ='/
     *  1p6e13.5)
      end
