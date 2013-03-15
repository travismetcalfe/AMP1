      subroutine bcscen(x2,y2,nn)
c
c  driving routine for setting up central expansion coefficients.
c  based on s/r bcs.
c
c  original version: 29/12/1984.
c
c  modified 4/1/1985 to account for change in dependent variable
c
c  modified 31/7/91, generalizing to take into account CNO cycle.
c
c  Modified 19/7/95, to use general heavy-element abundance
c  as set in common/heavy/. 
c
c  Modified 14/7/96, to set He3 abundance from xzer3, at time0,
c  when agehe3 .lt. 0
c
c  Modified 12/8/02, for inclusion of 4He burning
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
      implicit double precision (a-h,o-z)
      logical time0,nosd,notd,dtest,skipt,noder
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
c  Note: engenr.n.d.incl replaced by engenr.nnz.d.incl, 5/6/02
      include 'engenr.nnz.d.incl'
c
      parameter(naztmx = nspcmx + 3)
c
      dimension y2(*),az(naztmx),gh(2),dgh(2,naztmx),xhe3fg(4)
      common/bccn/ b1,b2,b3,b4,nb
      common/heavy/ zatmos, zhc, zh(1)
      common/ln10/ amm
      common/cmtime/ age, time0
      common/bccomp/ xhc,xhs,xcnoc(nspcmx),xhecc(2)
      common/eqstd/ dummy(14),rho(20),ht(20),p(20)
      common/eqstcl/ iradp,mode,dtest,skipt,noder
      common/logf/ fl
      common/he3fdg/ agesh,ifdhe3
      common/compsz/ xzerh, yzer, xzer3, xrz12, xrz13, xrz14, xrz16
      common/compvr/ cvr(icvrmx,1)
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec
      common/engcnt/ xhzlm1, xhzlm2, nxhzer, icnocs, nspec,
     *  iheccs, npsect
      common/diagns/ idgbcs
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      if(time0) then
        xh=xhc
      else
        xh=y2(5)
c
c  for .not. time0 set He3. Use y2(6) is apropriate,
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
	if(istdpr.gt.0) 
     *  write(istdpr,*) ' X3 set to',xhe3,'  in bcscen'
      end if
c
   10 yh=1-xh-zhc
      nosd=.true.
      notd=.true.
      skipt=.false.
      call eqstf(y2(2),y2(3),xh,yh,zhc,nosd,notd)
c
c  test for setting xhe3, from fudged evolution at constant conditions
c
      if((time0.or.ifdhe3.gt.0).and.agesh.gt.0) then
c
        nosd=.true.
        call he3abd(y2(2),y2(3),xh,yh,zhc,agesh,xhe3fg,anu,nosd)
        xhe3=xhe3fg(1)
c
        if(idgbcs.ge.1.and.istdpr.gt.0) write(istdpr,120) xhe3
      else if(time0) then
	xhe3=xzer3
      end if
c
      az(1)=1.e-17*p(1)
      az(2)=1.e-7*1.d1**y2(3)
      az(3)=xh
      az(4)=xhe3
      az(5)=1.d1**y2(1)
c
c  set az for CNO variables, depending on whether or not
c  this is time0
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
c
      acy=1.e-7
      nmax=10
      iter=1
c
      fl=y2(2)
c
      call bciter(az,gh,dgh,acy,eam,iconv,nmax,iter,nb)
c
      return
  120 format(//' in s/r bcscen, xhe3 has been reset to',1pe13.5)
      end
