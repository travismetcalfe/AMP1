      subroutine rxhmxc(x,y,iy,nn,compc,qc,rxmn,drxmn,icomp,alres,icry)
c
c  calculates average rate of change rxmn of composition over convective
c  core, and drxmn = d rxmn/ dX
c
c  Also returns, in qc, mass fraction in convective core
c
c  Returns also reset luminosity variable
c  (log10(L/1.e33)) in alres, found by integrating energy generation
c
c  On input x and y are assumed to contain log q 
c  and evolution variables, in usual form. 
c  Also compc is the assumed uniform composition of the convective core.
c
c  Assumes that variables defining the boundary of mixed core have
c  been set up in common/cmxcor/ by call of s/r mixcor
c
c  Should later be modified to save computation of equation of state
c  and all that, by storing the necessary variables when available
c
c  Original version: 22/8/90.
c
c  Modified 16/10/90, to set rates of change of several elements.
c
c  Modified 11/8/91, for general treatment of CNO cycle
c
c  Modified 16/5/92, to use new setting of mixed core in 
c  common /cmxcor/
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
c  Modified 19/7/95, to use general heavy-element abundance
c  as set in common/heavy/. This has to be checked with care in
c  connection with convective mixing.
c
c  Modified 21/12/00, to return also reset luminosity variable
c  (log10(L/1.e33)) in alres, found by integrating energy generation
c  rate with mixed composition. Note that this is returned in new
c  argument to routine.
c
c  Modified 4/1/01, to use heavy-element abundance in zhc (assumed to
c  have been set to core composition in s/r cmpcvc).
c
c  Modified 3/10/02, including epsg (which must have been set previously
c  in s/r rhs) in recalculation of luminosity. Also include possible
c  luminosity shift set in alshft.
c
      implicit double precision(a-h, o-z)
      include 'engenr.bz.d.incl'
c
      parameter(idr1mx = nspcmx+3, naztmx = nspcmx+3)
      logical nosd, notd
      dimension x(nn), y(iy,nn), compc(icomp), rxmn(icomp), 
     *    drxmn(icomp), alres(1)
      dimension yint(6), ft(idr1mx,nspcmx), eps(idr1mx), dlrxmn(nspcmx),
     *  allint(1)
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev, syear
      common/totmss/ am, rs
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,jnc,frcf(6),frcl(6)
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor
      common/he3eql/ xhe3eq,ieqhe3,iequsd
      common/he3fdg/ agesh,ifdhe3
      common/engcnt/ xhzlm1, xhzlm2, nxhzer
      common/enggrv/ epsg(nnmax)
      common/clshft/ alshft
      common/heavy/ zatmos,zhc,zh(1)
      common/ksider/ al01,al21,aztst(naztmx),axst(naztmx)
      common/eqstd/ xi(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     .  dlt(4),gm1,tprh,trhp,rhxp
      common/anwvar/ datdnn(8), yi(istrmx,1)
      common/xnwvar/ q(1)
      common/cdgphs/ kdgeos, kdgopc, kdgeng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      data irxprt /1/
c
      icry=0
c
      kdgeos=1
c
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'entering rxhmxc with nn, iy, compc(1-icomp) =',
     *  nn, iy, (compc(i),i=1,icomp)
c..	write(6,*) 'n, x, y'
c..	do n=nn-10,nn
c..	  write(6,'(i5,1p6e13.5)') n,x(n),(y(i,n),i=1,5)
c..	end do
c..      write(6,*) '#D# mshstr call in rxhmxc (1)'
c..      call mshstr(x,y,in1,in,iy,-1,-1,rhs,resrc)
c
c  fudge factor for possible testing 
c
      alffct=1.0d0
c
c  switch on extra diagnostics for small convective core
c
      if(qmxcor.gt.1.e-4) then
        idiag=1
      else
        idiag=2
      end if
c
      if(alffct.ne.1.and.istdpr.gt.0) 
     * write(istdpr,'(/'' Fudge eps with factor'',f10.5)') alffct
c
c  set up variables for integration
c
      nosd=.true.
      notd=.true.
c
      if(nmxcor.ge.nn.or.nmxcor.eq.0) then
	nc1=nn
      else if(frmxc.gt.1) then
        nc1=nmxcor+1
      else
        nc1=nmxcor
      end if
      xhc=compc(1)
c
c  test for extrapolating to first point
c
      if(abs(frmxc-1.d0).ge.1.d-6.and.nmxcor.lt.nn
     *  .and.nmxcor.gt.0) then
	do 20 i=1,5
   20   yint(i)=frmxc*y(i,nmxcor)+(1-frmxc)*y(i,nmxcor-1)
        qli=frmxc*x(nmxcor)+(1-frmxc)*x(nmxcor-1)
        epsgi=frmxc*epsg(nmxcor)+(1-frmxc)*epsg(nmxcor-1)
	q(1)=10.d0**qli
	if(istdpr.gt.0) write(istdpr,'(a,i5,1p3e19.11)') 
     *    'nmxcor, frmxc, qli, q(1)', nmxcor, frmxc, qli, q(1)
	fl=yint(2)
	tl=yint(3)
	yh=1-xhc-zhc
	call eqstf(fl,tl,xhc,yh,zhc,nosd,notd)
        if(kdgeos.lt.0) go to 90
c
	t=10.d0**tl
	call engenr(fl,tl,compc,yh,zhc,eps,ft,idr1mx,nosd)
c
        do 22 i=1,icomp
	yi(2*i-1,1)=ft(1,i)
   22   yi(2*i  ,1)=ft(4,i)
	yi(2*icomp+1,1)=eps(1)+epsgi
	n1=1
      else
	n1=0
      end if
c
c  set values at intermediate points.
c
      if(irxprt.gt.0.and.istdpr.gt.0)  then
        write(istdpr,'(a,2f14.11)') 'in rxhmxc, xhc, zh =', xhc, zhc
	write(istdpr,*) 'n, q, fl, T, rho'
      end if
      do 30 n=nc1,nn
      qq=10.d0**x(n)
      if(n1.eq.0.or.qq.lt.q(1)) then
        n1=n1+1
	q(n1)=qq
        fl=y(2,n)
        tl=y(3,n)
        yh=1-xhc-zhc
        call eqstf(fl,tl,xhc,yh,zhc,nosd,notd)
        if(kdgeos.lt.0) go to 90
c
        t=10.d0**tl
        call engenr(fl,tl,compc,yh,zhc,eps,ft,idr1mx,nosd)
c
        do 28 i=1,icomp
        yi(2*i-1,n1)=ft(1,i)
   28   yi(2*i  ,n1)=ft(4,i)
	yi(2*icomp+1,n1)=eps(1)+epsg(n)
c
        if(nn-n.le.5.and.irxprt.gt.0.and.istdpr.gt.0) then
          write(istdpr,'(i5,1p4e16.9)') n, q(n1), fl, t, rho(1)
        end if
c 
      end if
c
   30 continue
c
      qc=q(1)
c
c  set values at central point
c
      nnc=n1+1
      tc=1.d7*aztst(2)
      tl=log10(tc)
      yh=1-xhc-zhc
c
c  need to iterate to get log f at centre, starting with
c  trial from last meshpoint (already set)
c
   35 call eqstf(fl,tl,xhc,yh,zhc,nosd,notd)
      if(kdgeos.lt.0) go to 90
      dfl=log10(1.d17*aztst(1)/pt(1))/pt(2)
      if(abs(dfl).gt.1.e-8) then
	fl=fl+dfl
	if(irxprt.gt.0.and.istdpr.gt.0) 
     *    write(istdpr,*) 'fl, dfl =',fl, dfl
	go to 35
      end if
c
      call engenr(fl,tl,compc,yh,zhc,eps,ft,idr1mx,nosd)
c
      q(nnc)=0
      do 37 i=1,icomp
      yi(2*i-1,nnc)=ft(1,i)
   37 yi(2*i  ,nnc)=ft(4,i)
      yi(2*icomp+1,nnc)=eps(1)
c
      if(irxprt.gt.0.and.istdpr.gt.0) then
        write(istdpr,'(i5,1p4e13.5)') n, q(nnc), fl, tc, rho(1)
      end if
c
c  integrate
c
      icomp2=2*icomp+1
      do 40 k=1,icomp2
   40 call vinta(q,yi(k,1),yi(k+icomp2,1),nnc,istrmx,istrmx)   
c
      if(irxprt.gt.1.and.istdpr.gt.0) then
          write(istdpr,40090) (n,q(n),(yi(i,n),i=1,4),
     *      abs(yi(2,n))*xhc/(abs(yi(1,n))+1.e-30),n=1,nnc)
40090   format(//' n, q, yi(1 - 4), dlrxdx:'/(i5,1p6e12.4))
        irxprt=irxprt-1
      end if
c
      if(irxprt.gt.1.and.istdpr.gt.0) 
     *  write(istdpr,'(a/(i5,1p9e11.3))')
     *  'n, q, integrands, integrals for C12 and C13:',
     *  (n, q(n), (yi(i,n),i=5,8), (yi(icomp2+i,n),i=5,8), n=1,nnc,5)
c
c  set and print average reaction rates
c
      do 50 i=1,icomp
      j=icomp2+2*i-1
      rxmn(i)=-yi(j,nnc)/q(1)
      drxmn(i)=-yi(j+1,nnc)/q(1)
      if(rxmn(i).ne.0) then
	dlrxmn(i)=compc(i)*drxmn(i)/rxmn(i) 
      else
	dlrxmn(i)=0
      end if
      if(istdpr.gt.0) write(istdpr,135) i, rxmn(i), drxmn(i), dlrxmn(i)
   50 continue
      if(istdpr.gt.0) write(istdpr,'('' rxmn(1) ='',1pe20.12)') rxmn(1)
c
c  set and possibly print integrated luminosity
c
   55 if(idiag.ge.3.and.istdpr.gt.0) then
	write(istdpr,140) (n,q(n),yi(icomp2,n),yi(2*icomp2,n),n=1,nnc)
      end if
c
      n1=nn+1
      kl=2*icomp2
      amms=amsun*am
      if(idiag.ge.2.and.istdpr.gt.0) write(istdpr,150)
      do 60 n=nnc-1,2,-1
      n1=n1-1
      q1=10.d0**x(n1)
      all=amms*(yi(kl,n)-yi(kl,nnc))
c
c  as a test, increase luminosity
c
      all=alffct*all
c
      all=all/1.d33+alshft
      if(all.gt.0) then
        alres(n1)=log10(all)
      else
        alres(n1)=-30.d0
      end if
      if(idiag.ge.2.and.istdpr.gt.0) write(istdpr,155) 
     *  n1,q1,q(n),yi(icomp2,n), y(4,n1),alres(n1)
   60 continue
c
      if(idiag.ge.2.and.istdpr.gt.0) write(istdpr,*) ' '
c
      call lir1(log10(q(1)),x,allint,y(4,1),1,iy,nn,1,inter)
      all1=alffct*amms*(yi(kl,1)-yi(kl,nnc))/1.d33
c
      delal=all1-10.d0**allint(1)+alshft
      if(idiag.ge.2.and.istdpr.gt.0) 
     *  write(istdpr,*) 'q(edge), delal =',q(1), delal
      do 65 n=n1-1,1,-1
      all=all/1.d33+alshft
      all=10.d0**y(4,n)+delal+alshft
      if(all.gt.0) then
        alres(n)=log10(all)
      else
        alres(n)=-30.d0
      end if
      alres(n)=log10(10.d0**y(4,n)+delal+alshft)
      if(idiag.ge.2.and.n1-n.le.20.and.istdpr.gt.0) 
     *  write(istdpr,157) n,10.d0**x(n),y(4,n),alres(n)
   65 continue
c
      if(idiag.ge.1.and.istdpr.gt.0) write(istdpr,160) y(4,1), alres(1)
c
c  test for unreasonably big change in surface luminosity
c
      if(abs(y(4,1)-alres(1)).ge.0.2.and.idiag.ne.3.and.istdpr.gt.0) 
     *  then
	write(istdpr,165)
	idiag=3
	go to 55
      end if
c
      return
c
c  error return
c
   90 icry = -1
      write(istdou,'(/'' ***** Error return from rxhmxc'')')
      if(istdpr.gt.0.and.istdpr.ne.istdou)
     *  write(istdpr,'(/'' ***** Error return from rxhmxc'')')
      return
c
  110 format(//' ***** error in rxhmxc.',
     *  ' Convection zone does not extend to centre'/
     *         '       Last point at n =',i6)
  130 format(//' Output from rxhmxc.'/)
  135 format(' i =',i3,' rxmn =',1pe13.5,' drxmn =',e13.5,
     *         ' dlrxmn =',e13.5)
  140 format(//' n, q, eps, int(eps):'/(i5,1p3e13.5))
  150 format(//' Reset luminosity in s/r rxhmxc.'/
     *         ' n, q(mod), q(s/r), eps, org. log L*, new log L*:'/)
  155 format(i5,1p5e13.5)
  157 format(i5,1pe13.5,26x,2e13.5)
  160 format(/' Original, reset log(L) =',1p2e15.7)
  165 format(/' ***** Unreasonable luminosity change in rxhmxc.'/
     *        '       Make additional output'/)
      end
