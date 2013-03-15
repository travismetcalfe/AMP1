      subroutine bcsd(x1,x2,y1,y2,zk,ap,aq,g,gd,iaaa,ibbb,nn)
c
c  boundary condition subroutine for evolution.
c
c  Version for diffusion
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
c  Modified 21/7/95, to use general heavy-element abundance
c  as set in common/heavy/.
c
c  Modified 5/8/95, to allow diffusion of several elements.
c
c  Modified 21/8/95, setting zero-gradient condition on abundances
c  of diffusing non-reacting elements at centre
c
c  Corrected 5/8/97, by including engenr.nnz.d.incl, to set nspdmx etc
c
c  Modified 31/3/98, setting up central boundary condition for X
c  which includes second-derivative term (flagged by idfxbc = 2)
c  Note: derivatives of additional term are not included.
c
c  Modified 16/8/02, preparing for including 4He burning.
c
c  Modified 3/9/03, reducing eps (convergence in s/r atmos) from 1.e-6
c  to 1.e-10.
c
c  Modified 4/9/03, correcting zeroing of composition derivatives
c  for convective core.
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
     . az(naztmx),gh(2),dgh(2,naztmx),yg(2),isig(naztmx),xhe3fg(4),
     *  difx(nspdmx,6), velx(nspdmx,2,6), clamxy(4), grad(6), yr2(20)
      common/bccn/ b1,b2,b3,b4,nb
      common/bcatms/ taumn,taumx,sbcfct,flsatm,ntau
      common/rhcn/ aj(4),zdum,nvar
      common/clshft/ alshft
      common/heavy/ zatmos, zhc, zh(1)
      common/ln10/ amm
      common/totmss/ am, rscmm
      common/cmtime/ age, time0
      common/noiter/ iter1, ntime
      common/bccomp/ xhc,xhs,xcnoc(nspcmx),xhecc(2)
      common/ksider/ al01,al21,aztst(naztmx),axst(naztmx)
      common/eqstd/ dummy(14),rho(20),ht(20),p(20),cp(4),dad(4)
      common/eqstcl/ iradp,mode,dtest,skipt,noder
      common/logf/ fl
      common/he3int/ xhe3in(4)
      common/he3fdg/ agesh,ifdhe3
      common/compvr/ cvr(icvrmx,1)
      common/cstyst/ iyh1, iyhe3, iyc12, iyc13, iyn14, iyo16,
     *  iyche4, iycc12
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec, idcomp, iccomp
      common/cdiffu/ idiffus, itbdif, rhob, rmxdif, tlfdif, dtldif,
     *  ctdfmx, rctdif, rctdf1,idfxbc
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor
      common/engcnt/ xhzlm1, xhzlm2, nxhzer, icnocs, nspec,
     *  iheccs, nspect
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev, syear
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
        xh=y2(4)
c
c  for .not. time0 set He3. Use y2(7) if apropriate,
c  otherwise value in cvr.
c#ai# This should probably be checked more carefully.
c  Note that if ifdhe3 = 1,
c  this assignement is overridden by call to s/r he3abd below.
c
        if(ispxx3.eq.2) then
          xhe3=y2(5+2*idcomp)
        else
          xhe3=cvr(3,nn)
        end if
	if(idgbcs.ge.1.and.istdpr.gt.0) write(istdpr,*) 
     *    ' X3 set to',xhe3,'  in bcs'
      end if
c
c  for idiffus = 2 use value at innermost meshpoint and at the
c  same time store value in zhc for later use
c
      if(idiffus.eq.2) then
        zhc=y2(5)
      end if
c
   10 if(.not.time0.and.iheccs.ne.0.and.xh.le.1.e-5) then
	yh=y2(iyche4)
	zhc=1-yh
      else
	yh=1-xh-zhc
      end if
      nosd=.true.
      notd=.true.
      skipt=.false.
      call eqstf(y2(2),y2(3),xh,yh,zhc,nosd,notd)
c
c  test for setting xhe3, from fudged evolution at constant conditions
c
      if(time0.or.ifdhe3.gt.0) then
c
        nosd=.true.
	idgheo=idghe3
	idghe3=0
        call he3abd(y2(2),y2(3),xh,yh,zhc,agesh,xhe3fg,anu,nosd)
	idghe3=idgheo
        xhe3=xhe3fg(1)
c
        if(istdpr.gt.0.and.idgbcs.ge.1) write(istdpr,120) xhe3
      end if
c
      az(1)=1.e-17*p(1)
      az(2)=1.e-7*1.d1**y2(3)
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
	  if(icnocs.ge.1) call store(y2(4+2*idcomp+ispxx3), az(6), 
     *      ispcno)
	  if(iheccs.ne.0) 
     *      call store(y2(4+2*idcomp+ispxx3+ispcno), az(6+ispcno), 2)
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
      if(istdpr.gt.0.and.idgbcs.ge.2) write(istdpr,125) (az(i),i=1,naz)
c
      call bciter(az,gh,dgh,acy,eam,iconv,nmax,iter,nb)
c
c  set the boundary conditions and their derivatives
c
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'After call of bciter, idcomp =',idcomp
      nb1=nb+1
      yg(1)=x2+b1
      yg(2)=y2(4+idcomp)
      do 20 ial=1,2
      ig=2+idcomp+ial
      g(ig)=yg(ial)-gh(ial)
c
c  restore derivatives in proper order
c  When there is a mixed core, do not set derivatives with respect
c  to composition variables
c
      do 20 ia=1,nb+1
      if(ia.le.3) then
	j = ia+1
      else if(ia.le.nb) then
	j = ia+1+2*idcomp
      else
	j = 1
      end if
      if(nmxcor.le.0.or.j.le.3) gd(ig,j)=-dgh(ial,ia)
   20 continue
c
      gd(4+idcomp,4+idcomp)=1.d0
c
c  zero-flux or zero-gradient conditions 
c  (should be properly set up, with expansion, but certainly not now)
c
      initcp=0
      call cpydif(y2,yr2,ii1,ii2,ii3,1,idiffus,20,2,initcp)
      drad=fdrad(x2,yr2,zhc,ak,akr,akt,akx)
      if(istdpr.gt.0.and.idgbcs.ge.2) write(istdpr,*) 
     *  'In bcs, drad =',drad
      if(drad.lt.dad(1)) then
	grad(1)=drad
      else 
	grad(1)=dad(1)
      end if
      tc=10.d0**y2(3)
      rc=1.d11*10.d0**y2(1)
      amms=am*amsun
      amassc=amms*10.d0**x2
      if(istdpr.gt.0.and.idgbcs.ge.2) write(istdpr,*) 
     *  'In bcs, tc, rc, amassc, p(1), rho(1), grad(1) =',
     *  tc, rc, amassc, p(1), rho(1), grad(1)
c
c  Retain possibility of zero-velocity condition on hydrogen
c  Also, do not include derivatives of diffusion coefficients
c
      noder = .true.
      call difcff(tc,rho,p,xh,zhc,rc,amassc,grad,idiffus,
     *  itbdif,clamxy,velx,difx,nspdmx,noder)
      do 25 j=1,idcomp
      if(j.eq.1.and.idfxbc.eq.0) then
        g(4+idcomp+j)=y2(4+idcomp+j)
        gd(4+idcomp+j,4+idcomp+j)=1
      else if(j.eq.1.and.idfxbc.eq.1) then
	g(4+idcomp+j)=y2(4+idcomp+j)-velx(j,1,1)*y2(3+j)/amms
	gd(4+idcomp+j,4+idcomp+j)=1
	gd(4+idcomp+j,3+j)=-velx(j,1,1)/amms
      else if(j.eq.1) then
	rcst=10.d0**y2(1)
	qc=10.d0**x2
	g(4+idcomp+j)=y2(4+idcomp+j)-velx(j,1,1)*y2(3+j)/amms
     *     -0.666666667d0*axst(3)*rcst**2*difx(1,1)/(qc*amms*amms)
	gd(4+idcomp+j,4+idcomp+j)=1
	gd(4+idcomp+j,3+j)=-velx(j,1,1)/amms
      else
	g(4+idcomp+j)=y2(4+idcomp+j)-velx(j,1,1)*y2(3+j)/amms
     *                -velx(j,2,1)*y2(3+j)*y2(5+idcomp)
	gd(4+idcomp+j,4+idcomp+j)=1
	gd(4+idcomp+j,3+j)=-velx(j,1,1)/amms-velx(j,2,1)*y2(5+idcomp)
	gd(4+idcomp+j,5+idcomp)=-velx(j,2,1)*y2(3+j)
      end if
   25 continue
c
c     ---------------------------------------------------
c
c  surface conditions
c
      if(time0) then
	xh=xhs
      else
        xh=y1(4)
      end if
c
c  for idiffus = 2 use value at outermost meshpoint and at the
c  same time store value in zatmos for later use
c
      if(idiffus.eq.2) then
        zatmos=y1(5)
      end if
c
      yh=1-xh-zatmos
      nosd=.true.
      notd=.true.
      skipt=.false.
      call eqstf(y1(2),y1(3),xh,yh,zatmos,nosd,notd)
      rl=y1(1)
      pl=log10(p(1))
c
c  test for atmosphere integration
c
      if(ntau.le.1) go to 30
c  solve in atmosphere
      ams=1.d33*10.d0**b1
      ars=1.d11*10.d0**y1(1)
      als=1.d33*(10.d0**y1(4+idcomp)-alshft)
      if(idgbcs.ge.1) idiag=idgatm
c
      call atmos(taumn,taumx,ntau,ams,ars,als,xh,zatmos,
     .  pls,dplrs,dplls,sbcfct,idiag,nmax,eps,icry)
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
      gd(1,4+idcomp)=-dtls
      gd(2,1)=-dplrs
      gd(2,2)=p(2)
      gd(2,3)=p(3)
      gd(2,4)=p(4)
      gd(2,4+idcomp)=-dplls
      go to 50
c
c  use simple condition
c
   30 tl=y1(3)
      rhl=log10(rho(1))
      call opact(rhl,tl,xh,zatmos,ak,rkr,rkt,rkx)
      akt=rkt+rho(3)*rkr
      akf=rho(2)*rkr
      akx=0
      if(xh.gt.0) akx=rkx/(amm*xh)+rho(4)*rkr
c
   35 g(1)=y1(4+idcomp)-2.d0*rl-4.d0*tl-b3
      g(2)=2.d0*rl+ak+pl-b4
c
      gd(1,1)=-2.d0
      gd(1,3)=-4.d0
      gd(1,4+idcomp)=1.d0
      gd(2,1)=2.d0
      gd(2,2)=p(2)+akf
      gd(2,3)=p(3)+akt
      gd(2,4)=akx+p(4)
   50 continue
c
c  zero-flux conditions
c
      do 55 j=1,idcomp
      g(2+j)=y1(4+idcomp+j)
   55 gd(2+j,4+idcomp+j)=1
c
      if(idgbcs.lt.1.or.istdpr.le.0) return
      write(istdpr,150)
      do 60 i=1,4+2*idcomp
   60 write(istdpr,160) i,g(i),(gd(i,j),j=1,nvar)
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
