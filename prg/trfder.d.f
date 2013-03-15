      subroutine trfder(fl, tl, xh, yh, zh, a,atr,trder,itr,init)
c
c  transforms thermodynamic first derivatives from initial independent
c  variables (log f, log T, X) (as used by eqstf) to
c  (ln p, ln rho, Y), as required by MJT
c
c  a contains original first derivatives, and atr returns transformed
c  derivatives. trder returns the transformation matrix.
c
c  The transformation matrix is recalculated if (fl, tl, xh) 
c  has changed since last call, or init = 1
c
c  Note: eqstf is assumed to have been called before call of trfder.
c
c  Original version: 7/4/93
c
      implicit double precision (a-h, o-z)
      dimension a(1), atr(1), trder(itr,1)
      dimension w(3,3)
      common/eqstd/ xii1(4),ane(10),rho(20),ht(20),p(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp,gmm1(4)
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev1,exhm
      common/ln10/ amm
      data flp, tlp, xhp /-1.e30,-1.e30,-1.e30/
c
      if(init.eq.1.or.flp.ne.fl.or.tlp.ne.tl.or.xhp.ne.xh) then
c
c  set matrices
c
        do 15 i=1,3
        do 12 j=1,3
   12   trder(i,j)=0
        trder(i,i)=1
        w(1,i)=amm*p(i+1)
        w(2,i)=amm*rho(i+1)
   15   w(3,i)=0
c
        yh=1-xh-zh
        w(3,3)=-1.
c
        call leq(w,trder,3,3,3,itr,err)
c
	flp=fl
	tlp=tl
	xhp=xh
c
      end if
c
c  set transformed derivatives
c
      do 30 i=1,3
      atr(i)=0
      do 30 j=1,3
   30 atr(i)=atr(i)+trder(j,i)*a(j)
c
      return
      end
