      subroutine setcbr(x,y,rhs,nn,iy,icscbr,iqcres)
c
c  sets up borders of convection zones (e.g. in trial model) by calling
c  s/r rhs, and afterwards sets parameters for convective overshoot
c  by calling s/r setovs.
c
c  If icscbr .ge. 1, in addition analyses gradient in core to look
c  for near instability. The action, in that case, depends on icscbr:
c    icscbr = 1: reset parameters for mixed core.
c    icscbr = 2: no action in s/r setcbr, action may be taken outside
c
c  The occurrence of this problem is flagged by returning iqcres = icscbr.
c
c  Original version: 16/1/2000
c
c  Modified: 5/6/2000, to analyze for near instability
c
c  Modified 18/4/06, setting iter1 to 0 before calling rhs
c
      implicit double precision(a-h,o-z)
      logical time0, time0p
      common/noiter/ iter1, ntime
      common/cnvout/ ddrad,rrcnv,aacnv,ddacad,amach,
     *  drr(5), ddr(5), da(5)
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor,nmxcp,qmxcp,rmxcp,
     *  qmxmin,qmxmax
      common/crxcor/ cqc, cqcp
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
      dimension x(*),y(iy,*)
      dimension zk(1),f(ivarmx),fd(ivarmx,ivarmx),zz(ivarmx),
     *  dzdy(ivarmx,ivarmx),alam(1,ivarmx),alamd(1,ivarmx,ivarmx),
     .  h(1),hd(1),ddrads(nnmax)
      common/cmtime/ age, time0
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,jnc,frcf(6),frcl(6),
     *  xrcf(6),xrcl(6)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      external rhs
c
      idiag=0
c
      if(istdpr.gt.0) write(istdpr,100) icscbr, iqcres
c
      time0p = time0
      iter = -1
      iter1p = iter1
      iter1 = 0
      iz = ivarmx
      ihd = 1
      ialam = 1
c
      do 20 n=1,nn
c
      call rhs(x(n),y(1,n),zk,zz,dzdy,ap,aq,f,fd,alam,alamd,h,hd,
     .  iz,iz,ihd,ialam,n,iter)
   20 ddrads(n)=ddrad
c
      time0 = time0p
      iter1 = iter1p
c
c  test for testing for convective-core problems
c
      if(icscbr.ge.1) then
        iqcres=0
c
c  locate local maximum in superadiabatic gradient
c
	ddrmax=-1
	nmax=-1
	if(idiag.ge.1.and.istdpr.gt.0) write(istdpr,*) 'n, q, ddrad:'
	do 30 n=nn-1,1,-1
	if(x(n).le.-2.3.and.idiag.ge.1.and.istdpr.gt.0) 
     *    write(istdpr,*) n, 10.**x(n),ddrads(n)
	if(ddrads(n).lt.0.and.x(n).le.-2.3.and.ddrads(n).gt.ddrmax.and.
     *    ddrads(n-1).le.ddrads(n).and.ddrads(n).ge.ddrads(n+1)) then
	  ddrmax=ddrads(n)
	  nmax=n
	  if(idiag.ge.1.and.istdpr.gt.0)
     *      write(istdpr,*) ' >>>>>>>>>>>>>>>>>>>>>>  Max.'
        end if
   30   continue
c
c  test for extended mixing (note that criterion is so far hard-wired;
c  that should be changed).
c
	if(nmax.gt.0.and.ddrmax.ge.-1.e-3) then
	  iqcres=icscbr
	  if(icscbr.eq.1) then
            qmxcor=10.d0**x(nmax)
	    cqc=qmxcor
	    rmxcor=10.d0**(y(1,nmax)-y(1,1))
	    nmxcor=nmax
	    frmxc=0
	    write(istdou,130) qmxcor,ddrmax
	    if(istdpr.gt.0.and.istdou.ne.istdpr)
     *        write(istdpr,130) qmxcor,ddrmax
          end if
        end if
c
      end if
c
      incovs = 0
      call setovs(x,y,nn,iy,incovs)
      return
  100 format(//' Entering s/r setcbr with icscbr =',i2,
     *  '  iqcres =',i2)
  130 format(//' ***** Warning. In s/r setcbr qmxcor has been reset to',
     *  1pe13.5/
     *         '       at local maximum = ',e13.5,' in ddrad'//)
      end
