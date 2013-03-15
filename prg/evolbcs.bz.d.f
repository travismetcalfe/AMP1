      subroutine bcs(x1,x2,y1,y2,zk,ap,aq,g,gd,iaaa,ibbb,nn)
c
c  boundary condition subroutine for evolution.
c
c  modified 20/8/1984 to include scaling of variables for
c  central conditions.
c
c  modified 19/12/1984 to correct treatment of he3 at innermost
c  mesh point.
c
c  modified 4/1/85 to use as dependent variables log(r/1.d11) and
c  log(l/1.d33).
c
c  21/9/87: implementing modifications from RECKU
c
c  ******************************************************************
c
c  31/7/91: Begin implementing changes required to generalize
c  treatment of CNO cycle.
c
c  ******************************************************************
c
c  Modified 19/7/95, to use general heavy-element abundance
c  as set in common/heavy/. 
c
c  Modified 14/7/96, to set He3 abundance from xzer3, at time0,
c  when agehe3 .lt. 0
c
c  Modified 29/5/03, reducing eps (convergence in s/r atmos) from 1.e-6
c  to 1.e-10.
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
      logical time0,nosd,notd,dtest,skipt,noder,qmxtst
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
c
      parameter(naztmx = nspcmx + 3)
c
      dimension y1(1),y2(1),zk(1),ap(1),aq(1),g(1),gd(iaaa,1),
     . az(naztmx),gh(2),dgh(2,naztmx),yg(2),isig(naztmx),xhe3fg(4)
      common/bccn/ b1,b2,b3,b4,nb
      common/bcatms/ taumn,taumx,sbcfct,flsatm,ntau
      common/rhcn/ aj(4),zdum,nvar
      common/clshft/ alshft
      common/heavy/ zatmos, zhc, zh(1)
      common/ln10/ amm
      common/cmtime/ age, time0
      common/noiter/ iter1, ntime
      common/bccomp/ xhc,xhs,xcnoc(nspcmx),xhecc(2)
      common/ksider/ aax(12)
      common/eqstd/ dummy(14),rho(20),ht(20),p(20)
      common/eqstcl/ iradp,mode,dtest,skipt,noder
      common/logf/ fl
      common/he3int/ xhe3in(4)
      common/he3fdg/ agesh,ifdhe3
      common/compsz/ xzerh, yzer, xzer3, xrz12, xrz13, xrz14, xrz16
      common/compvr/ cvr(icvrmx,1)
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec
      common/cstyst/ iyh1, iyhe3, iyc12, iyc13, iyn14, iyo16,
     *  iyche4, iycc12
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor
      common/engcnt/ xhzlm1, xhzlm2, nxhzer, icnocs, nspec,
     *  iheccs, nspect
      common/engche/ xmxrhe
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3
c
      common/cdgrhb/ kdgrhb
      common/cdgphs/ kdgeos, kdgopc, kdgeng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
c  setting of old permutation array for derivatives, valid only
c  for nb = 4, properly speaking. Note that for the moment this
c  is not used.
c
c..      data isig/2,3,5,6,1,0/
c
      data idgatm/1/
      data eps /1.d-10/
c
c  in present version, initialize flags for physics diagnostics, to
c  allow for return.
c  Flag initialized to 1: in physics routines, return rather than stop,
c  resetting flag to -1
c
      kdgeos=1
      kdgopc=1
      kdgeng=1
c
c  initialize flag for error in routine
c
      kdgrhb=1
c
c  extra diagnostics for small convective cores
c
      qmxtst=qmxcor.ge.1.e-8.and.qmxcor.le.1.e-4
c
      if(istdpr.gt.0.and.(idgbcs.ge.2.or.qmxtst)) then
	write(istdpr,105) (y1(i),i=1,nvar)
	write(istdpr,106) (y2(i),i=1,nvar)
      end if
c
c  central conditions
c
      if(time0) then
        xh=xhc
      else
        xh=y2(5)
c
c  for .not. time0 set He3. Use y2(6) if apropriate,
c  otherwise value in cvr.
c#ai# This should probably be checked more carefully.
c  Note that if ifdhe3 = 1,
c  this assignement is overridden by call to s/r he3abd below.
c
        if(ispxx3.eq.2) then
          xhe3=y2(6)
        else
          xhe3=cvr(3,nn)
        end if
	if(istdpr.gt.0) write(istdpr,*) ' X3 set to',xhe3,'  in bcs'
      end if
c
   10 if(.not.time0.and.iheccs.ne.0.and.xh.le.1.e-5
     *    .and.x2.le.xmxrhe) then
	yh=y2(iyche4)
	zhc=1-yh-xh
	if(idgbcs.ge.1.and.istdpr.gt.0) write(istdpr,*) 
     *    'Set zhc in evolbcs, X, Y, Z =',xh, yh, zhc
      else
	yh=1-xh-zhc
      end if
      nosd=.true.
      notd=.true.
      skipt=.false.
      call eqstf(y2(2),y2(3),xh,yh,zhc,nosd,notd)
      if(kdgeos.lt.0) then
        kdgrhb=-1
        return
      end if
c
c  test for setting xhe3, from fudged evolution at constant conditions
c
      if((time0.or.ifdhe3.gt.0).and.agesh.gt.0) then
c
        nosd=.true.
	idgheo=idghe3
	idghe3=0
        call he3abd(y2(2),y2(3),xh,yh,zhc,agesh,xhe3fg,anu,nosd)
	idghe3=idgheo
        xhe3=xhe3fg(1)
c
        if(idgbcs.ge.1.and.istdpr.gt.0) write(istdpr,120) xhe3
      else if(time0) then
	xhe3=xzer3
      end if
c
      az(1)=1.d-17*p(1)
      az(2)=1.d-7*1.d1**y2(3)
      az(3)=xh
      az(4)=xhe3
      az(5)=1.d1**y2(1)
c
c  set az for CNO and possibly 4He-burning variables, depending on whether 
c  or not this is time0
c
      if(icnocs.ge.1) then
	if(time0) then
	  if(icnocs.ge.1) call store(xcnoc(1), az(6), ispcno)
	  if(iheccs.ge.1) call store(xhecc(1), az(6+ispcno), 2)
        else
	  if(icnocs.ge.1) call store(y2(5+ispxx3), az(6), ispcno)
	  if(iheccs.gt.0) 
     *      call store(y2(5+ispxx3+ispcno), az(6+ispcno), 2)
        end if
      end if
      naz=5+ispcno+isphec
c
      acy=1.e-7
      nmax=10
      iter=iter1-1
c
      fl=y2(2)
c
      if(idgbcs.ge.2.and.istdpr.gt.0) write(istdpr,125) (az(i),i=1,naz)
c
      call bciter(az,gh,dgh,acy,eam,iconv,nmax,iter,nb)
c  set the boundary conditions and their derivatives
c
c  test for return
c
      if(kdgrhb.lt.0) return
c
      nb1=nb+1
      yg(1)=x2+b1
      yg(2)=y2(4)
      do 20 ial=1,2
      ig=2+ial
      g(ig)=yg(ial)-gh(ial)
c
c  restore derivatives in proper order
c  When there is a mixed core, do not set derivatives with respect
c  to composition variables
c
      do 20 ia=1,nb+1
      if(ia.le.2) then
	j = ia+1
      else if(ia.le.nb) then
	j = ia+2
      else
	j = 1
      end if
      if(nmxcor.le.0.or.j.le.4) gd(ig,j)=-dgh(ial,ia)
   20 continue
c
      gd(4,4)=1.d0
c
c     ---------------------------------------------------
c
c  surface conditions
c
      xh=y1(5)
      if(time0) xh=xhs
      yh=1-xh-zatmos
      nosd=.true.
      notd=.true.
      skipt=.false.
      call eqstf(y1(2),y1(3),xh,yh,zatmos,nosd,notd)
      if(kdgeos.lt.0) then
        kdgrhb=-1
        return
      end if
      rl=y1(1)
      pl=log10(p(1))
c
c  test for atmosphere integration
c
      if(ntau.le.1) go to 30
c  solve in atmosphere
      ams=1.d33*10.d0**b1
      ars=1.d11*10.d0**y1(1)
      als=1.d33*(10.d0**y1(4)-alshft)
      if(idgbcs.ge.1) idiag=idgatm
c
      call atmos(taumn,taumx,ntau,ams,ars,als,xh,zatmos,
     .  pls,dplrs,dplls,sbcfct,idiag,nmax,eps,icry)
c
c  test for return
c
      if(kdgrhb.lt.0) return
c
      idgatm=1
c  if not converged, use standard condition
      if(icry.eq.-1) go to 30
c  set temperature
      call ttau(taumx,ams,ars,als,xh,zatmos,ta,qhopf,dqhopf,
     *  dtrs,dtls,dttau)
      tal=log10(ta)
c  set conditions
      g(1)=y1(3)-tal
      g(2)=pl-pls
      gd(1,1)=-dtrs
      gd(1,3)=1
      gd(1,4)=-dtls
      gd(2,1)=-dplrs
      gd(2,2)=p(2)
      gd(2,3)=p(3)
      gd(2,4)=-dplls
      if(time0) then
        igd=4
      else
        igd=nvar
        gd(2,5)=p(4)
      end if
      go to 50
c
c  use simple condition
c
   30 tl=y1(3)
      rhl=log10(rho(1))
      call opact(rhl,tl,xh,zatmos,ak,rkr,rkt,rkx)
      if(kdgopc.lt.0) then
        kdgrhb=-1
        return
      end if
      akt=rkt+rho(3)*rkr
      akf=rho(2)*rkr
      akx=0
      if(xh.gt.0) akx=rkx/(amm*xh)+rho(4)*rkr
c
   35 g(1)=y1(4)-2.d0*rl-4.d0*tl-b3
      g(2)=2.d0*rl+ak+pl-b4
c
      gd(1,1)=-2.d0
      gd(1,3)=-4.d0
      gd(1,4)=1.d0
   40 gd(2,1)=2.d0
      gd(2,2)=p(2)+akf
      gd(2,3)=p(3)+akt
      if(time0) then
        igd=4
      else
        gd(2,5)=akx+p(4)
        igd=nvar
      end if
   50 continue
      if(idgbcs.lt.1.or.istdpr.le.0) return
      write(istdpr,150)
      do 60 i=1,4
   60 write(istdpr,160) i,g(i),(gd(i,j),j=1,igd)
      write(istdpr,165) yg(1), gh(1), yg(2), gh(2)
      return
  105 format(/' Entering bcs with'/
     *        ' y1: ',1p5e15.7/(5x,5e15.7))
  106 format(' y2: ',1p5e15.7/(5x,5e15.7))
  120 format(//' in s/r bcs, xhe3 has been reset to',1pe13.5)
  125 format(/' Call bciter with'/
     *        ' az:',1p5e15.7/(4x,5e15.7))
  150 format(//' output from bcs. i, g(i), gd(i,j):'/)
  160 format(i4,1pe13.5,5x,12e13.5)
  165 format(/' y(1), gh(1) =',1p2e13.5/
     *        ' y(2), gh(2) =',1p2e13.5)
      end
