      subroutine setdt(x,y,yp,dymx,dtmx,dtmn,age,agefns,dt,dtp,xdtden,
     *  aldtrt,dtfcrd,dtceps,iy,nn)
c
c  determines new time step dt
c
c  Modification, 17/10/02, to reduce effects of thin shell burning:
c  
c  The reduction factor is calculated as
c  
c  dtfcrd + (1 - dtfcrd)*rel_int/(rel_int + dtceps)
c  
c  where
c  
c  rel_int = (int change d x)/max(change)
c  
c  and change is the actual change in the variable considered
c  on which the time step might be decided.
c  This reduction factor is applied to the maximum change evaluated over
c  the model. When the change is concentrated in a very narrow region
c  (as is the case, e.g., for the hydrogen abundance in the hydrogen
c  shell source), rel_int is small and the reduction factor is substantially
c  below 1. Evidently, the smallest value of the reduction factor is 
c  dtfcrd.
c
c  This is (currently) applied to the changes in the hydrogen abundance
c  and luminosity.
c
c  **************************************************************************
c
c  corrected 19/7/1985 to reset ifin to 0 in last time step.
c
c  modified 2/3/88: limit increase in time step over previous
c  value to at most 20 per cent.
c
c  modified 10/8/90, to introduce factor on change in X in convective
c  core, to achieve better time resolution with the current (incorrect)
c  way of treating the convective core.
c  Factor xccfct currently hardcoded in this routine
c
c  Modified 14/8/91, to take into account the implementation of
c  the CNO cycle.
c
c  Modified 31/5/92, to check for near hydrogen exhaustion
c  in mixed core
c
c  Modified 27/7/96, adding diagnostics if idgsdt ne 0
c
c  Modified 1/7/97, testing for possible problems in reaching 
c  final age
c
c  Modified 23/9/98, changing test for growing convective core
c  (*** NOTE *** There may be problems with the definition of
c  nmxcor).
c
c  Modified 16/11/99, eliminating effect of change in X in growing
c  convective envelope. 
c  *** NOTE ***: This is currently flagged by setting ixcenv to 1
c  by hard-coding in the routine. Should perhaps be changed to
c  a control flag later.
c
c  Modified 17/10/02, introducing rescaling of changes in X and L
c  in case of narrow hydrogen-burning shell
c
c  Modified 20/6/03, applying suppression of change in growing
c  convective envelope also to log f (or log rho)
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
      implicit double precision (a-h,o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
c
      dimension x(1),y(iy,1),yp(iy,1)
      dimension dy(20),dymax(20),nymax(20),sdy(2,nnmax)
      common/rhcn/ dumrh(5),nvar
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,jnc,frcf(6),frcl(6),
     *  xrcf(6),xrcl(6)
      common/convpp/ dcsp,dccp,nfp(6),nlp(6),incp,jncp,
     *  frcfp(6),frclp(6),xrcfp(6),xrclp(6)
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor,nmxcp,qmxcp,rmxcp,
     *  qmxmin,qmxmax
      common/engcnt/ xhzlm1, xhzlm2, nxhzer, icnocs, nspec,
     *  iheccs, nspect
      common/cstyst/ iyh1, iyhe3, iyc12, iyc13, iyn14, iyo16, 
     *  iyche4, iycc12
      common/noiter/ iter, ntime
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3,idgmxc,
     *  idgsdt,idgn75,idgn76,idgn77
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      data idgsdt /0/
c
      save
c
      data init, ifin/1, 0/
      data xccfct /5.d0/
      data iopnfl /0/
c
c  hard-coded factor for test near hydrogen exhaustion
c
      data alexhs /0.8d0/
c
c hard code flag for eliminating X in growing convective envelope
c
      data ixcenv /1/
c
      if(istdpr.gt.0) write(istdpr,'(//'' Entering setdt''/)')
c
      if(istdpr.gt.0) then
        write(istdpr,*) 'In setdt limits of convective envelope are',
     *    xrclp(1), xrcl(1)
        write(istdpr,*) 'nmxcp, nmxcor =',nmxcp, nmxcor
      end if
c
      if(iopnfl.eq.0) then
        open(62,file='ttt.setdt.out',status='unknown')
        if(idgsdt.le.1) then
          write(62,140)
        else
          write(62,150)
        end if
        iopnfl=1
      end if
c
c  set test for convective core
c
      if(nmxcor.gt.0) then
       ncc=nmxcor+1
      else
       ncc=nn+1
      end if
c
c  test for approaching end
c
      if(ifin.gt.0) then
c
c  last timestep before agefns. Test that final age will be reached
c
        if(age+dt.ge.agefns) then
c
c  do not reset dt. restore ifin and init
c
          ifin=0
          init=1
          return
c
        else
c
c  in case of problems, continue as usual
c
	  ifin=0
	  init=0
c
        end if
      end if
c
      if(age+2*dt.ge.agefns) then
        dtp=dt
        dt=(agefns-age)/1.9999999d0
        ifin=1
        return
c
      end if
c
      if(dymx.ge.1) go to 80
      chmx=0.d0
      ichmx=0  
      nchmx=0  
c
      call zero(dymax,4+nspec)
c
      do 45 n=1,nn 
c
c  for testing, switch off effect of 4He etc.
c
c..      do 45 i=1,nvar
      do 45 i=1,4+nspec
      yy=y(i,n)
      yyp=yp(i,n)  
      chmxn= abs(yy-yyp)
      if(i.eq.4) then
c  
c  luminosity  
c  
        alrat=10.**(yy-y(4,1))
        chmxn=chmxn*alrat*(1+aldtrt)/(alrat+aldtrt)  
c
      else if(i.eq.5.or.(iheccs.ne.0.and.i.eq.iyche4)) then
c  
c  hydrogen or helium abundance, test for convective core 
c#ai# To avoid problems at edge of a convective core
c     (particularly if it is growing) apply factor only in
c     bulk of convectively mixed core
c     Note that change in possibly growing region
c     of convective core is suppressed below
c  
        if(n.lt.ncc+2) then
          chmxn=chmxn/(abs(yy)+xdtden) 
        else
          chmxn=xccfct*chmxn/(abs(yy)+xdtden) 
        end if
c
      else if(i.eq.6.and.ispxx3.eq.2) then
c  
c  he3 abundance
c  
        chmxn=chmxn/5.e-3
c
      end if
c
c  Neglect changes in composition or log f (log rho) in possibly 
c  growing region of convective core
c
c  Also, neglect growing convective envelope, if ixcenv = 1
c  (so far hard-coded above)
c  
      if(i.eq.2.or.i.eq.5.or.(iheccs.ne.0.and.i.eq.iyche4)) then
        if((n.le.nmxcp.and.n.ge.nmxcor-1).or.
     *    (ixcenv.eq.1
     *       .and.xrclp(1)+0.01.ge.x(n).and.x(n).ge.xrcl(1)-0.01))
     *       chmxn=0
c..	  write(6,*) 'x =',x(n),' in growing envelope'
      end if
c  
      if(chmxn.gt.chmx) then
        chmx=chmxn 
        ichmx=i
        nchmx=n
	if(i.eq.5.and.n.ge.ncc) then
	  ixccmx=1
        else
	  ixccmx=0
        end if
      end if
c
c  test for maximum in individual variable
c
      if(chmxn.gt.dymax(i)) then
        dymax(i)=chmxn
	nymax(i)=n
      end if
c
c  add to integrals of changes of X and luminosity (for possible
c  suppression in narrow shell)
c
      if(i.eq.4.or.i.eq.5) then
	i1=i-3
	if(n.eq.1) then
	  sdy(i1,n)=0.d0
        else
	  sdy(i1,n)=sdy(i1,n-1)+chmxn*(x(n-1)-x(n))
        end if
      end if
c
   45 continue 
c
c  reduce for maximum changes confined to narrow burning shell
c  Note: the reduction is not applied for the first
c
      chratl=sdy(1,nn)/dymax(4)
      chratx=sdy(2,nn)/dymax(5)
      if(dtceps.gt.0) then
        fcratl=dtfcrd+(1.d0-dtfcrd)*chratl/(chratl+dtceps)
        fcratx=dtfcrd+(1.d0-dtfcrd)*chratx/(chratx+dtceps)
      else
        fcratl=1.d0
        fcratx=1.d0
      end if
c
c#ai#  the following override is kept for testing purposes.
c#ai#  should be removed in next iteration!
c
c  Removed 17/1/03
c
c..      fcratl=0.10d0+0.75d0*chratl/(chratl+2.d-2)
c..      fcratx=0.10d0+0.75d0*chratx/(chratx+2.d-2)
      ires=0
      if(fcratl.le.0.95.and.ntime.gt.1) then
	if(istdpr.gt.0) write(istdpr,'(/a,1pe13.5)')
     *    ' Luminosity change reduced by factor',fcratl
	dymax(4)=fcratl*dymax(4)
	if(ichmx.eq.4) ires=1
      end if
c
      if(fcratx.le.0.95.and.ntime.gt.1) then
	if(istdpr.gt.0) write(istdpr,'(/a,1pe13.5)')
     *    ' Hydrogen   change reduced by factor',fcratx
	dymax(5)=fcratx*dymax(5)
	if(ichmx.eq.5) ires=1
      end if
c
c  test for resetting maximum change
c
      if(ires.eq.1) then
	chmx=0
	do i=1,4+nspec
	  if(dymax(i).gt.chmx) then
	    ichmx=i
	    chmx=dymax(i)
          end if
        end do
      end if
c
c  output full list of changes
c
      if(istdpr.gt.0) write(istdpr,'(/a/(i4,i5,1pe13.5))')
     *  ' Changes in s/r setdt. i, nmax, max. change:',
     *  (i,nymax(i),dymax(i),i=1,4+nspec)
c
      dtp=dt
      dtn=dt*dymx/chmx 
c
c  test for near, but not complete, exhaustion in mixed core
c
      if(nmxcor.gt.0.and.yp(5,nn).gt.y(5,nn).and.y(5,nn).gt.1.e-9) then 
	dtrmax=(1-alexhs)*y(5,nn)/(yp(5,nn)-y(5,nn))
	if(dtrmax.gt.0.and.dtn.gt.dtrmax*dtp) then
	  dtn=dtrmax*dtp
	  ixccmx=2
        end if
      end if
c
c  to avoid instability in the mesh determination, prevent
c  too rapid increase in step length, except for the first step
c  Note: the value of the limit may have to be adjusted.
c
      if(init.eq.0) then
        dtn=min(dtn,1.4*dt)
      else
        init=0
      end if
c
c  limit dt to be between dtmn and dtmx
c
      dt=max(dtmn,min(dtn,dtmx))
c
      if(istdpr.gt.0) then
        if(ixccmx.eq.0) then
          write(istdpr,110) chmx,ichmx,nchmx
	  if(ichmx.eq.5) write(istdpr,111) nmxcp, nmxcor
  111     format(/' nmxcp, nmxcor =',2i5)
	  npmin=max(1,nchmx-5)
	  npmax=min(nn,nchmx+5)
	  write(istdpr,112) 
     *      (n,x(n),y(ichmx,n),yp(ichmx,n),n=npmin,npmax)
  112     format(//' n, x, y, yp near maximum change:'/
     *      (i5,1p3e13.5))
        else if(ixccmx.eq.1) then
          write(istdpr,115) chmx, xccfct
        else
          write(istdpr,130) yp(5,nn),y(5,nn),alexhs
        end if
      end if
c
c  if time step is reduced drastically, output detailed
c  changes to file
c
      if(idgsdt.lt.1.and.dt/dtp.ge.0.5) return
c
      do 60 n=1,nn 
      do 55 i=1,nvar
      yy=y(i,n)
      yyp=yp(i,n)  
      dy(i)=yy-yyp
      if(i.eq.4) then
c  
c  luminosity  
c  
        alrat=10.**(yy-y(4,1))
        dy(i)=dy(i)*alrat*(1+aldtrt)/(alrat+aldtrt)  
c
      else if(i.eq.5) then
c  
c  hydrogen abundance, test for convective core 
c#ai# To avoid problems at edge of a convective core
c     (particularly if it is growing) apply factor only in
c     bulk of convectively mixed core
c  
        if(n.lt.ncc+2) then
          dy(i)=dy(i)/(abs(yy)+xdtden) 
        else
          dy(i)=xccfct*dy(i)/(abs(yy)+xdtden) 
        end if
c
      else if(i.eq.6.and.ispxx3.eq.2) then
c  
c  he3 abundance
c  
        dy(i)=dy(i)/5.e-3
c
      end if
   55 continue
      if(idgsdt.le.1) then
        write(62,'(2i5,1p15e13.5)') ntime, n, 10.d0**x(n),
     *    (dy(i),i=1,nvar)
      else
        write(62,'(2i5,1p45e13.5)') ntime, n, 10.d0**x(n),
     *    (yp(i,n),i=1,nvar),(y(i,n),i=1,nvar),
     *    (dy(i),i=1,nvar)
      end if
   60 continue
      call flush(62)
      return
c
   80 dtp=dt
      dt=dtmx  
      return
  110 format(//' in setdt maximum change is',1pe13.5,3x,
     *  ' this occurs in variable',i3,'  at point',i5) 
  115 format(//' in setdt maximum change is',1pe13.5,3x,
     *  ' in X in convective core. Enhancement factor =',0pf8.3)
  130 format(//' in setdt, dt set for near hydrogen exhaustion.',
     *  ' Xcp, Xc =',1p2e13.5,'  alpha =',e13.5)
  140 format('# Output from setdt'/'#'/
     *  '# n, q, changes in y:'/'#')
  150 format('# Output from setdt'/'#'/
     *  '# n, q, yp, y, changes in y:'/'#')
      end  
