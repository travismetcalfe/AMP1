      subroutine atmos(taumn,taumx,ntau,ams,ars,als,x,z,pls,dplrs,
     .  dplls,sbcfct,idiag,itmax,eps,icry)
c  finds log(pressure) at taumx by integrating hydrostatic support
c  equation from taumn to taumx. mesh taken to be uniform in log(tau)
c
c  ********************************************************************
c
c  NOTE: With normal setting of input parameters, Eddington limit is
c  reach on the ZAMS at a mass of around 55 Msun.
c
c  ********************************************************************
c
c  21/9/87: implementing modifications from RECKU
c
c  modified 27/12/85 to limit size of fl in iteration to 0.5 
c
c  modified 13/3/1986 to pass trial fl on surface (fls) in
c  common /bcatms/  
c
c  modified 15/3/88 to use consistent numerical constants
c
c  modified 5/3/90, to increase size of arrays for atmospheric solution
c  from 51 to 201
c
c  Modified 19/7/95, to use general heavy-element abundance
c  as set in common/heavy/. This has to be checked with care in
c
c  Modified 19/8/96, possibly including rotational effects
c
c  Modified 12/10/98, modifying initial guess for log(f)
c  in case of convergence problems.
c
c  Modified 3/6/03, trying yet another resetting of initial guess
c  for surface log(f) in case of convergence problems in 
c  surface iteration.
c
c  Modified 20/10/04, introducing flag for testing for errors in physics
c  routines and return with flag.
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
      implicit double precision (a-h,o-z)
      logical nosd,notd,dtest,skipt,noder,qmxtst
      dimension yy(2,2),b(5,2),pratm(5,201)
      common/catmos/ natm,matm,vatm(5,201),patm(5,201)
      common/comgrt/ isprot, irotcn, a5, riner0, riner, omgrt0,
     *  omgrot(1)
      common/bcatms/ tmnbc,tmxbc,sbcfc1,fls,ntaubc,iopatm,icsrad
      common/ln10/ amm
      common/opcfdg/ alfa
      common/cder/ yd(3,201)
      common/eqstd/ xii(4),dum(10),rho(20),ht(20),p(20),dmmm(12),gm1
      common/eqstcl/ iradp,mode,dtest,skipt,noder
      common/noiter/ iter, ntime
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor,nmxcp,qmxcp,rmxcp,
     *  qmxmin,qmxmax
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev, syear
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3,idgmxc,
     *  idgsdt
c
      common/cdgrhb/ kdgrhb
      common/cdgphs/ kdgeos, kdgopc, kdgeng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      data iphsra /0/
c
c  maximum allowed change in log(f)
c
      data dflmax /0.5d0/
c
      save
c
c  extra diagnostics for small convective cores
c
      qmxtst=qmxcor.ge.1.e-8.and.qmxcor.le.1.e-4
      if(istdpr.gt.0.and.(idgbcs.ge.2.or.qmxtst)) then
        write(istdpr,*) 
     *    'Enter atmos with taumn,taumx,ntau,ams,ars,als,x,z ='
        write(istdpr,*) taumn,taumx,ntau,ams,ars,als,x,z
      end if
c
      natm=max0(ntau,1)
c
      matm=0
c
c  gravitational acceleration, possibly including rotational effect
c
      if(isprot.eq.0) then
        g=cgrav*ams/(ars*ars)
      else
        g=cgrav*ams/(ars*ars)-0.666666667d0*ars*omgrot(1)**2
      end if
c
      fcm=4*pi*ars*ars
c
      icry=1
c
      in=1
      in1=2
c  initial value for log f
      flss=fls
      flss0=flss
      ksrfit = 0
      nosd=.true.
      notd=.true.
      skipt=.false.
c
c  as an initial step, run through suitable values of fl
c  when idgbcs .ge. 2
c
      if(idgbcs.ge.2) then
        flss1=flss-2
        flss2=flss+2
        nfls=20
        dfls=(flss2-flss1)/(nfls-1)
        if(istdpr.gt.0) 
     *  write(istdpr,
     *    '(//'' n, log(f), log(rho), log(p), log(kappa), dphi:'')')
        tau=taumn
        do nfl=1,nfls
	  fl=flss1+dfls*nfl
          dphi=surdph(tau,ams,ars,als,sbcfct,g,fl,x,z,akl)
c
          if(kdgrhb.lt.0) return
c
	  if(istdpr.gt.0) write(istdpr,
     *      '(i4,4f11.6,1pe13.5)') nfl, fl, log10(rho(1)),
     *     log10(p(1)), akl, dphi
        end do
      end if
c
c  initialize counters resetting starting point
c
    5 nsrfit = 0
      fl=flss
c
      if(idiag.ge.1.and.istdpr.gt.0) write(istdpr,100)
c  tau scale
      n=1
      tau=taumn
      if(ntau.gt.1) then
        dlt=log10(taumx/taumn)/(ntau-1)
        dft=10.d0**dlt
        dlth=dlt*0.5d0
      end if
c  temperature
   10 call ttau(tau,ams,ars,als,x,z,t,qhopf,dqhopf,
     *  dtrs,dtls,dttau)
      tl=log10(t)
c  set equation of state and opacity
      y=1-x-z
      it=0
c
      deqn=1e10
      deqnmn=1e10
      flmn=0
      isecnt=0
      istep=0
c
      if(idiag.ge.1) then
	itmaxp=itmax+5
      else
	itmaxp=itmax
      end if
c
   15 continue
c
      call eqstf(fl,tl,x,y,z,nosd,notd)
      if(kdgeos.lt.0) then
        kdgrhb=-1
        return
      end if
c
      skipt=.true.
c
      rhl=log10(rho(1))
      pl=log10(p(1))
c  when iopatm = 1 call separate atmospheric opacity
      if(iopatm.eq.1) then
        call opacat(rhl,tl,x,akl,akr,akt,akx)
        if(kdgopc.lt.0) then
          kdgrhb=-1
          return
	end if
      else
c  normal opacity
        call opact(rhl,tl,x,z,akl,akr,akt,akx)
        if(kdgopc.lt.0) then
          kdgrhb=-1
          return
        end if
      end if
c
      ak=10.d0**akl
      phi=tau*g/(ak*p(1))
      akf=akr*rho(2)
c  type of iteration
      if(n.eq.1) then
c  surface boundary iteration
	dphip=dphi
        dphi=2*sbcfct*phi-1
        if(isecnt.eq.0.and.istep.eq.0) 
     *    dfl=dphi/(amm*(p(2)+2*sbcfct*phi*akf))
        if(idiag.ge.1.and.istdpr.gt.0) 
     *    write(istdpr,110) it,fl,pl,rhl,akf,dphi,dfl,xii
c
        if(abs(dfl).ge.1.and.istep.eq.0.and.it.ge.3) then
c
c  as a test, start small steps in log f
c
	  if(istdpr.gt.0.and.istdou.ne.istdpr) 
     *      write(istdpr,'(/'' Start small steps in log f''/)')
	  write(istdou,'(/'' Start small steps in log f''/)')
	  istep=1
	  fl=fl+1
	  dflst=-0.02
	  ist=1
	  go to 15
        end if
c
	if(istep.eq.1) then
	  if(dphi*dphip.lt.0.and.abs(dphi)+abs(dphip).le.0.2
     *      .and.ist.ge.2) then
	    istep=0
	    isecnt=1
          else if(abs(ist*dflst).ge.2) then
	    if(istdpr.gt.0.and.istdou.ne.istdpr) 
     *        write(istdpr,'(/'' Give up for atmos stepping'')')
            write(istdou,'(/'' Give up for atmos stepping'')')
	    go to 85
          else
	    ist=ist+1
	    dfl=dflst
	    fl=fl+dflst
            go to 15
c
	  end if
	end if
c
	if(isecnt.ge.1) then
	  if(isecnt.eq.1) 
     *      write(istdou,'(/'' Start secant iteration in log f''/)')
	  isecnt=isecnt+1
	  if(isecnt.gt.100) then
	    if(istdpr.gt.0) write(istdpr,*) 
     *        'Excessive number of secant iterations (1)'
            if(kdgrhb.gt.0) then
              kdgrhb = -1
              return
            else
	      stop 'evolatmos 1'
            end if
          end if
c..	  write(6,*) 'old dfl, dphi, dphip =',dfl,dphi,dphip
	  dfl=-dfl*dphi/(dphi-dphip)
c..	  write(6,*) 'New dfl =',dfl
	  fl=fl+dfl
	  if(abs(dfl).lt.eps) go to 40
	  go to 15
        end if
c
c  test for apparent lack of convergence
c
	if(abs(dfl).gt.dflmax) then
	  nsrfit=nsrfit+1
	  if(nsrfit.eq.4) then
            ksrfit=ksrfit+1
c
c  restart with reset initial value
c
	    if(ksrfit.eq.1) then
	      flss=flss-sign(2*dflmax,dfl)
	      isdfl=-sign(1.d0,dfl)
            else
	      isdfl=-isdfl
	      flss=flss0+2*isdfl*(ksrfit-1)*dflmax
	    end if
	    if(istdpr.gt.0) then
	      write(istdpr,*) 'ksrfit, isdfl',ksrfit, isdfl
	      write(istdpr,105) flss
	    end if
	    go to 5
	  end if
	end if
      else
c  main iteration
	deqnp=deqn
	deqn=dlth*(phi+phip)+plp-pl
c
c  temporary setting of dfl and output
c
        denom=p(2)+dlth*amm*phi*(p(2)+akf)
        dflnew=deqn/(p(2)+dlth*amm*phi*(p(2)+akf))
        if(idiag.ge.1.and.istdpr.gt.0.and.
     *    (istep.eq.1.or.abs(dflnew).ge.1))
     *    write(istdpr,110) it,fl,pl,rhl,akf,deqn,dflnew,denom
        if(abs(dflnew).ge.1.and.istep.eq.0) then
c
c  as a test, start small steps in log f
c
	  write(istdou,'(/'' Start small steps in log f''/)')
	  istep=1
	  fl=fl+1
	  dfl=-0.02
	  ist=1
	  go to 15
        end if
c
	if(istep.eq.1) then
	  if(deqn*deqnp.lt.0.and.abs(deqn)+abs(deqnp).le.0.2
     *      .and.ist.ge.2) then
	    istep=0
	    isecnt=1
          else if(abs(ist*dfl).ge.2) then
	    if(istdpr.gt.0.and.istdou.ne.istdpr) 
     *        write(istdpr,'(/'' Give up for atmos stepping'')')
            write(istdou,'(/'' Give up for atmos stepping'')')
	    go to 85
          else
	    ist=ist+1
	    fl=fl+dfl
            go to 15
c
	  end if
	end if
c
c  test for cycling of solution
c
	if((deqnp.ne.0.and.deqn/deqnp.lt.-0.5)
     *    .or.isecnt.ge.1) then
c
c  switch to secant-like iteration, possibly with undercorrection
c  (in a last desperate attempt to get on)
c
	  isecnt=isecnt+1
	  if(isecnt.gt.100) then
	    if(istdpr.gt.0) write(istdpr,*) 
     *        'Excessive number of secant iterations (2)'
            if(kdgrhb.gt.0) then
              kdgrhb = -1
              return
            else
	      stop 'evolatmos 1'
            end if
          end if
	  if(abs(deqnp).gt.abs(deqnmn)) then
	    deqnp=deqnmn
	    fl=flmn
          else
	    deqnmn=deqnp
	    flmn=fl
          end if
	  dfl=-dfl*deqn/(deqn-deqnp)
	  if(isecnt.gt.4) dfl=dfl/2
	  if(idiag.ge.1.and.istdpr.gt.0) write(istdpr,18090) 
     *      fl,deqnp, deqn
18090     format(' secant iteration. fl, deqnp, deqn =',f10.5,1p2e13.5)
        else
          dfl=deqn/(p(2)+dlth*amm*phi*(p(2)+akf))
        end if
      end if
c
c  limit size of dfl
c
      dfla=abs(dfl) 
      if(dfla.gt.dflmax) dfl=sign(dflmax,dfl) 
c
c  test for convergence
c
      fl=fl+dfl 
      if(dfla.lt.eps) then
        go to 40
      else if(it.lt.itmaxp) then
        it=it+1
	if(istdpr.gt.0.and.idiag.ge.1.and.it.gt.itmax) then
	  write(istdpr,142) it,n,tau,fl,pl,rhl,p(2),akf
        end if
        go to 15
      else
        if(istdpr.gt.0.and.istdpr.ne.istdou) 
     *    write(istdpr,140) it,n,tau,dfl,fl,pl,rhl
        write(istdou,140) it,n,tau,dfl,fl,pl,rhl
        if(dfla.gt.1.e-4) go to 85
      end if
c
c  iteration finished
c
   40 vatm(1,n)=tau
      vatm(2,n)=fl
      vatm(3,n)=tl
c  set output quantities
      pratm(1,n)=tau
      pratm(2,n)=10.**tl
      pratm(3,n)=p(1)
      pratm(4,n)=ak
      pratm(5,n)=rho(1)
c  set pressure derivatives, and mass and depth
      akp=akf/p(2)
      rhtp=rho(3)-rho(2)*p(3)/p(2)
      aktp=akt+akr*rhtp
      b(1,in)=-amm*phi*(2+aktp*dtrs)
      b(2,in)=-amm*phi*aktp*dtls
      b(3,in)=amm*phi*(1+akp)
      b(4,in)=amm*fcm*tau/ak
      b(5,in)=amm*tau/(ak*rho(1))
c
      if(n.eq.1) then
c  surface values
        do 43 i=1,2
   43   yy(i,n)=b(i,in)/b(3,in)
        if(idiag.ge.2.and.istdpr.gt.0) write(istdpr,120)
        fls=fl
        vatm(4,n)=0
        vatm(5,n)=0
      else
c  interior values
        fc1=1/(1+dlth*b(3,in))
        fc2=1-dlth*b(3,in1)
        do 47 i=1,2
   47   yy(i,in)=fc1*(yy(i,in1)*fc2+dlth*(b(i,in)+b(i,in1)))
c
        n1=n-1
        do 48 i=4,5
   48   vatm(i,n)=vatm(i,n1)+dlth*(b(i,in)+b(i,in1))
c
      end if
c
      ddt=dttau/phi
      if(idiag.ge.2.and.istdpr.gt.0) write(istdpr,130) 
     *  n,tau,pl,fl,tl,phi,yy(1,in),
     .  yy(2,in),vatm(4,n),vatm(5,n),ddt
c
      do 52 i=1,2
   52 yd(i,n)=yy(i,in)
      yd(3,n)=pl
c  set variables for pulsation
      patm(1,n)=phi/gm1
      patm(2,n)=phi*rho(2)/p(2)+rhtp*dttau
      patm(3,n)=rho(1)
      patm(4,n)=rho(1)*ak/tau
      patm(5,n)=ddt
      if(n.lt.ntau) then
c
        n=n+1
        i=in
        in=in1
        in1=i
        plp=pl
        phip=phi
        tau=tau*dft
c  extrapolate fl
	dfl=phi*dlt/p(2)
	if(abs(dfl).gt.1) dfl=sign(1.d0,dfl)
        fl=fl+dfl
        go to 10
c
      end if
c
c  finish
c
      pls=pl
      dplrs=yy(1,in)
      dplls=yy(2,in)
c
c  set fls to surface value
c
      fls=vatm(2,1)
c
      matm=natm
c  print hsra-like quantities
      if(iphsra.ne.1) return
      if(istdpr.gt.0) write(istdpr,150) (n,(pratm(i,n),i=1,5),n=1,natm)
c
      return
c  failure to converge. write diagnostics
   85 icry=-1
      if(kdgrhb.gt.0) then
        kdgrhb = -1
        write(istder,*) 'Stop in atmos'
        return
      else
        stop 'Stop in atmos'
      end if
  100 format(///' output from atmos.'//' iteration for surface p.',
     .  ' it,fl,pl,rhl,akf,dphi,dfl,xi:'/)
  105 format(//' Apparent convergence problems.'/
     *         ' Restart with surface trial log(f) =',f10.5)
  110 format(i4,3f10.5,1p3e11.3,2x,4e10.2)
  120 format(//' solution in atmosphere. n,tau,pl,fl,tl,phi,y(1-2),',
     .  ' delta mass, depth:'//)
  130 format(i4,1pe13.5,0p3f12.5,1p5e13.5,0pf10.5)
  140 format(///1x,10(1h*),' iteration failed to converge after',
     .  i10,' iterations in atmos, at n =',i4,' tau =',
     .  1pe13.5//' dfl,fl,pl,rhl =',e13.5,0p3f12.5)
  142 format(
     *  ' Convergence problems. it, n, tau, fl, pl, rhl, p(2), akf ='/
     *  2i4,1pe11.3,0p5f12.4)
  150 format(///' output as from hsra. n, tau, t, p, kappa, rho:'//
     *  (i4,1p5e13.5))
      end
      function surdph(tau,ams,ars,als,sbcfct,g,fl,x,z,akl)
c
c  determines error in surface boundary condition
c
      implicit double precision (a-h, o-z)
      logical nosd, notd
      common/eqstd/ xii(4),dum(10),rho(20),ht(20),p(20),dmmm(12),gm1
c
      common/cdgrhb/ kdgrhb
      common/cdgphs/ kdgeos, kdgopc, kdgeng
c
      call ttau(tau,ams,ars,als,x,z,t,qhopf,dqhopf,
     *  dtrs,dtls,dttau)
      tl=log10(t)
c  set equation of state and opacity
      y=1-x-z
      nosd=.true.
      notd=.true.
c
      call eqstf(fl,tl,x,y,z,nosd,notd)
      if(kdgeos.lt.0) then
        kdgrhb=-1
        return
      end if
c
      rhl=log10(rho(1))
      pl=log10(p(1))
c  opacity
      call opact(rhl,tl,x,z,akl,akr,akt,akx)
      if(kdgopc.lt.0) then
        kdgrhb=-1
        return
      end if
c
      ak=10.d0**akl
      phi=tau*g/(ak*p(1))
c
c  surface boundary condition
c
      surdph=2*sbcfct*phi-1
      return
      end
