      subroutine almmxc(x,y,iy,nn,compc,alrhmn,ialam,idalam)
c
c  calculates average rate of rho*lambda over convective core,
c  from reaction rates in alam(1,k), k = 1, ialam.
c  results are stored in alrhmn(1,k). (derivatives may later
c  be set into alrhmn(j,k)).
c  idalam is first dimension of alrhmn.
c
c  On input x and y are assumed to contain log q and evolution 
c  variables, in usual form. 
c  Also compc is the assumed uniform composition of the convective core.
c
c  Assumes that variables defining boundaries of convective core have
c  been set up in common/convpt/
c
c  Should later be modified to save computation of equation of state
c  and all that, by storing the necessary variables when available
c
c  Original version: 18/10/90.
c
c  Modified 29/5/92, to average instead over mixed core
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      implicit double precision(a-h, o-z)
c  Note: engenr.n.d.incl replaced by engenr.nnz.d.incl, 5/6/02
      include 'engenr.nnz.d.incl'
c
      parameter(naztmx = nspcmx + 3)
c
      logical nosd,notd,norct,dgnrxp
      dimension x(nn), y(iy,nn), compc(1), alrhmn(idalam,1)
      dimension yint(6)
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,jnc,frcf(6),frcl(6)
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor
      common/engcnt/ xhzlm1, xhzlm2, nxhzer
      common/rnrcnt/ irnrat, krnrat, irn12, irn13, irn14, 
     *  irn15c, irn15o, irn16, irn17
      common/rnratd/ al(10,krnrmx),norct
      common/heavy/ zh
      common/ksider/ al01,al21,aztst(naztmx),axst(naztmx)
      common/eqstd/ xi(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     .  dlt(4),gm1,tprh,trhp,rhxp
      common/anwvar/ datdum(8), yi(istrmx,1)
      common/xnwvar/ q(1)
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3,idgmxc
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      data irxprt /0/
c
c#ai#  For the moment, explicitly set irnrat = 1.
c#ai#  This needs to be fixed up later, also in s/r evolmain.
c
      irnrat = 1
c
      nosd=.true.
      notd=.true.
      dgnrxp=istdpr.gt.0.and.irxprt.gt.0
c
      if(nmxcor.ge.nn.or.nmxcor.le.0) then
	nc1=nn
      else if(frmxc.gt.1) then
        nc1=nmxcor+1
      else
        nc1=nmxcor
      end if
      xhc=compc(1)
c
c  test for interpolating to first point
c
      if(frmxc.ne.1.and.nmxcor.lt.nn.and.nmxcor.gt.0) then
	do 20 i=1,5
   20   yint(i)=frmxc*y(i,nmxcor)+(1-frmxc)*y(i,nmxcor-1)
        qli=frmxc*x(nmxcor)+(1-frmxc)*x(nmxcor-1)
	fl=yint(2)
	tl=yint(3)
	yh=1-xhc-zh
	call eqstf(fl,tl,xhc,yh,zh,nosd,notd)
c
	t=10.d0**tl
	call rnrate(fl,tl,xhc,yh,zh,nosd)
c
	q(1)=10.d0**qli
        do 22 i=1,ialam
   22   yi(i,1)=rho(1)*al(1,i)
	n1=1
      else
	n1=0
      end if
c
c  set values at intermediate points.
c
      if(dgnrxp) then
        write(istdpr,*) 'xhc, zh =',xhc, zh
        write(istdpr,*) 'n, q, fl, T, rho'
      end if
      do 30 n=nc1,nn
      n1=n1+1
      fl=y(2,n)
      tl=y(3,n)
      yh=1-xhc-zh
      call eqstf(fl,tl,xhc,yh,zh,nosd,notd)
c
      t=10.d0**tl
      call rnrate(fl,tl,xhc,yh,zh,nosd)
c
      q(n1)=10.d0**x(n)
      do 28 i=1,ialam
   28 yi(i,n1)=rho(1)*al(1,i)
c
      if(nn-n.le.5.and.dgnrxp) then
        write(istdpr,'(i5,1p4e13.5)') n, q(n1), fl, t, rho(1)
      end if
c
   30 continue
c
c  set values at central point
c
      nnc=n1+1
      tc=1.d7*aztst(2)
      tl=log10(tc)
      yh=1-xhc-zh
c
c  need to iterate to get log f at centre, starting with
c  trial from last meshpoint (already set)
c
   35 call eqstf(fl,tl,xhc,yh,zh,nosd,notd)
      dfl=log10(1.d17*aztst(1)/pt(1))/pt(2)
      if(abs(dfl).gt.1.e-8) then
	fl=fl+dfl
	if(dgnrxp) 
     *    write(istdpr,*) 'fl, dfl =',fl, dfl
	go to 35
      end if
c
      call rnrate(fl,tl,xhc,yh,zh,nosd)
c
      q(nnc)=0
      do 37 i=1,ialam
   37 yi(i,nnc)=rho(1)*al(1,i)
c
      if(dgnrxp) then
        write(istdpr,'(i5,1p3e13.5)') n, q(nnc), tc, rho(1)
      end if
c
c  integrate
c
      do 40 k=1,ialam
   40 call vinta(q,yi(k,1),yi(k+ialam,1),nnc,istrmx,istrmx)   
c
      if(dgnrxp) then
        write(istdpr,40090) (n,q(n),(yi(i,n),i=1,6),n=1,nnc)
40090   format(//' n, q, yi(1 - 6):'/(i5,1p7e12.4))
        irxprt=irxprt-1
      end if
c
c  set and print results
c
      do 50 i=1,ialam
      alrhmn(1,i)=-yi(ialam+i,nnc)/q(1)
      if(idgeng.ge.1.and.istdpr.gt.0) write(istdpr,135) i, alrhmn(1,i)
   50 continue
c
      return
  130 format(//' Output from almmxc.'/)
  135 format(' i =',i3,' alrhmn =',1pe13.5)
      end
