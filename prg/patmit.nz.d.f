      subroutine patmit(patmos,als,ars,ams,xh,nitmax,eps,idiag,
     *   icry)
c
c  iterates on opacity fudge factor alamop to get atmospheric
c  pressure (at bottom of atmosphere) of patmos.
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
c  Modified 19/7/95, to use general heavy-element abundance
c  as set in common/heavy/.
c
      implicit double precision (a-h,o-z)
      common/bcatms/ taumn,taumx,sbcfct,flsatm,ntau 
      common/opctcl/ alamop,zatmop,tlmopc,dtmopc,fdgopl,iwdopc,
     *  iopacm,ifdgop
      common/heavy/ zatmos, zhc, zh(1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      pl=log10(patmos)
      nit=0
      dalam=0.1
      if(alamop.le.0) alamop=1
c
      idiag1=0
      if(idiag.eq.2) idiag1=1
      if(idiag.ge.1.and.istdpr.gt.0) write(istdpr,100)
c  step in iteration
   10 nit=nit+1
      call atmos(taumn,taumx,ntau,ams,ars,als,xh,zatmos,pls,dplrs,
     *  dplls,sbcfct,idiag1,nitmax,eps,icry)
      if(icry.eq.-1) go to 20
      dpl=pl-pls
      if(idiag.ge.1.and.istdpr.gt.0) write(istdpr,110) 
     *  nit,alamop,pls,dpl
      if(nit.eq.1) go to 20
c  test for convergence of excessive number of iterations
      if(abs(dpl).le.eps) go to 30
      if(nit.eq.nitmax) go to 80
c  reset dalam
      dalam=dalam*dpl/(dplp-dpl)
c  reset alamop
   20 dplp=dpl
      alamop=alamop+dalam
      go to 10
c  iteration converged
   30 icry=1
      return
c
c  iteration failed to converge
c
   80 write(istdou,120) nit,alamop,pls,dpl
      if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,120) 
     *  nit,alamop,pls,dpl
      icry=-1
c
c  do not use matched opacity
c
   90 alamop=0
      return
  100 format(///' iterate on alamop to get given atmospheric',
     *  ' pressure'/)
  110 format(/' nit, alamop, pls, dpl =',i4,2f10.5,1pe13.5)
  120 format(///1x,10(1h*),' iteration for alamop failed to',
     *  ' converge in s/r patmit after',i4,' iterations'//
     *  11x,'final alamop, pls, dpl =',2f10.5,1pe13.5)
      end
