      program main  
c   
c  redistributes mesh in physical model 
c   
c  modified 9/1/1985 to allow use of central expansion to reset 
c  points in new model interior to first non-zero point in old  
c  model. old expansion may still be used by setting ioldex = 1.
c   
c  modified 17/5/1985. now interpolates in x**(-2)*aa(i,n)  
c  for i = 2 and 4 (these tend to a constant at x = 0), 
c  and sets aa(2,n) interior to first mesh point in original
c  model from expansion of p and other aa-s.
c   
c  modified 21/5/1985 in treatment of convection zone boundary  
c  (see s/r stcvzb) 
c   
c  modified 21/5/1985 to include contribution in stretching 
c  from superadiabatic gradient, to ensure adequate resolution very 
c  near surface.
c   
c  modified 29/3 - 1/4 1989 to include various fixes for stretching
c  polytropic model
c
c  modified 30/5/92, to include term in change in buoyancy frequency,
c  and smoothing with running mean.
c
c  modified 25/7/92, to smooth A4 near centre, particularly 
c  near hydrogen exhaustion
c
c  modified 17/2/98, to allow reading general-format amdl files
c  (including aa(6,.), say)
c
c  Modified 5/7/01, fixing criterion for discontinuity in rho, now
c  based on aa(4,.). This is, for the time being, flagged by a
c  negative ndisc. The old criterion is obtained by selecting
c  a posive ndisc.
c
c  Modified 28/1/02, adding contribution from cx at x = 0.
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90
c
      implicit double precision (a-h, o-z)
      parameter (nnmax = 10001, iaa=10)
      logical singsf
      dimension x(nnmax),aa(iaa,nnmax),xsi(nnmax),xn(nnmax),
     *  an(iaa,nnmax),
     *  data(8), nsdisc(2,100), xsdisc(2,100),
     *  ninter(2,10),nnew(2,10),nndisc(2,10),danslp(iaa)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c..      namelist/exec/ nn,icnmsh,
c..     *  icase,icvzbn,nsmth,ndisc,dlxdsc,
c..     *  cg,cx,ca,cdgr,cddgr,cdg,alphsf
c..     *  nout,cn,irsu,unew,iresa4,   
c..     *  nmodel,kmodel,itsaml,ioldex
c  defaults in exec 
c  nn: number of points in new mesh
      nn=601
      icnmsh=0  
c  icase: for icase = 1, 2, 3, 11 or 12 set standard parameters.
c     icase = 1, 2, 3 give old parameters, as used up until June 1988.
c     icase = 11 gives new, optimized p mode parameters
c     icase = 1  corresponds to p modes, old parameters
c     icase = 2  corresponds to g modes, old parameters
c     icase = 3  corresponds to intermediate modes, old parameters.
c     icase = 11 corresponds to p modes, new parameters.
c     icase = 12 corresponds to g modes, new parameters.
      icase=0   
c  icvzbn: for icvzbn .gt. 0 ensure that there is a mesh point at   
c     the lower boundary of the convective envelope.
c     for icvzbn = 1 interpolate normally in aa(i,.) for i .ne. 4,   
c     and reset aa(4,.) by linear interpolation near boundary.  
c     for icvzbn = 2 interpolate separately below and above  
c     boundary. 
      icvzbn=0  
c  nsmth: for nsmth .gt. 1, smooth integrand with Gaussian-weighted
c     running mean over nsmth points
      nsmth=0
c  ndisc: for ndisc ne 0 add abs(ndisc) points where there is a
c     discontinuity in rho
c     NOTE: ndisc lt 0 selects for new (post-05/07/01) criterion
      ndisc=0
c  dlxdsc: range in ln x over which discontinuity is distributed
      dlxdsc=1.e-3
c  dlgrmx: minimum (d log rho/d log p) for discontinuity
      dlgrmx=10
c  cacvzb: amplitude of point added at base of convection zone
c     should be used only in conjunction with smoothing
      cacvzb=0
c
      cg=1  
      cx=40.
      ca=0.01   
      cdgr=0
      cddgr=0   
c  cdg: weight for gradient in A in interior of model
c     if cdg .gt.0, limit contribution (before scaling by cdg)
c     to partial sum of other terms. If cdg .lt. 0, do not limit
c     (and use abs(cdg) as weight).
      cdg=0
c  alphsf: cuts off the singularity at the surface of a polytropic
c     model. In setting mesh, A2 is replaced by A2/(1+alphsf*A2).
      alphsf=0
      nout=100  
      cn=2  
      irsu=0
      unew=3
c  iresa4: if iresa4 = 1, reset A4 before redistribution,
c  to correct possible errors in A4 resulting from near exhaustion
c  of hydrogen
      iresa4=0
c  nmodel: for nmodel .gt. 1 read and redistribute nmodel models.   
c          otherwise only 1.
      nmodel=0  
c  kmodel: for kmodel .gt. 0, start at model no. kmodel 
      kmodel=0  
c  itsaml: if itsaml = 1 test input model for conversion errors 
      itsaml=0  
c  ioldex: if ioldex = 1, old expansion (not using second derivatives   
c     in data) is used. 
      ioldex=0  
c   
c                   ......................................
c
      nrd=0 
c
c  open files needed
c
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'Files needed: input  model on unit 2'
      if(istdpr.gt.0) 
     *  write(istdpr,*) '              output model on unit 3'
      call ofiles
      call openf(2,'o','u')
      call openf(3,'u','u')
c
c  read start of file, to check for number of variables
c
      read(2,end=80,err=90) nmod,nh,data
      close(2)
      call openf(2,'o','u')
      idata8 = idint(data(8)+0.1)
      if(idata8.ge.100) then
        ivar = 8
      else if(idata8.ge.10) then
        ivar = 6
      else
        ivar = 5
      end if
c
      ix=ivar+1
c
c  read exec
    5 continue
c..      read(5,exec,end=90)   
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'nn,icnmsh?'
      if(istdpr.gt.0) 
     *  write(istdpr,*) nn,icnmsh
      read(5,*,end=90,err=90) nn,icnmsh
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'icase,icvzbn,nsmth,ndisc,dlxdsc,dlgrmx,cacvzb?'
      if(istdpr.gt.0) 
     *  write(istdpr,*) icase,icvzbn,nsmth,ndisc,dlxdsc,dlgrmx,cacvzb
      read(5,*) icase,icvzbn,nsmth,ndisc,dlxdsc,dlgrmx,cacvzb
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'cg,cx,ca,cdgr,cddgr,cdg,alphsf?'
      if(istdpr.gt.0) 
     *  write(istdpr,*) cg,cx,ca,cdgr,cddgr,cdg,alphsf
      read(5,*) cg,cx,ca,cdgr,cddgr,cdg,alphsf
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'nout,cn,irsu,unew,iresa4   ?'
      if(istdpr.gt.0) 
     *  write(istdpr,*) nout,cn,irsu,unew,iresa4   
      read(5,*) nout,cn,irsu,unew,iresa4   
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'nmodel,kmodel,itsaml,ioldex?'
      if(istdpr.gt.0) 
     *  write(istdpr,*) nmodel,kmodel,itsaml,ioldex
      read(5,*) nmodel,kmodel,itsaml,ioldex
c   
c  test for setting standard parameters
c   
      if(icase.eq.1) then
c  old p modes  
        cg=0.001  
        cx=0.1
        ca=0.01   
        cn=2  
      else if (icase.eq.2) then
c  old g modes  
        cg=1  
        cx=40 
        ca=0.01   
        cn=2  
      else if(icase.eq.3) then
c  old intermediate modes   
        cg=0.05   
        cx=1  
        ca=0.01   
        cn=2  
      else if(icase.eq.11) then
c  new p modes  
        cg=0.001  
        cx=0.1
        ca=0.003   
        cn=2  
      else if (icase.eq.12) then
c  new g modes  
        cg=4  
        cx=40 
        ca=0.01   
        cn=2  
      end if
c
c  set cacvzb to zero, unless smoothing is applied
c
      if(nsmth.eq.0) cacvzb=0
c
c  test for old or new discontinuity formulation
c
      if(ndisc.lt.0) then
	ndisca=-ndisc
	if(istdpr.gt.0) 
     *  write(istdpr,103) ndisca
      else if(ndisc.gt.0) then
	ndisca=ndisc
	if(istdpr.gt.0) 
     *  write(istdpr,102) ndisca
      else
	ndisca=0
      end if
c
c  test for resetting ndisc
c
      if(ndisca.gt.0.and.ndisca.lt.20) ndisca=20
c   
c..      write(6,105)  
c..      write(6,exec) 
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'nn,icnmsh'
      if(istdpr.gt.0) 
     *  write(istdpr,*) nn,icnmsh
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'icase,icvzbn,nsmth,ndisc,dlxdsc,dlgrmx,cacvzb'
      if(istdpr.gt.0) 
     *  write(istdpr,*) icase,icvzbn,nsmth,ndisc,dlxdsc,dlgrmx,cacvzb
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'cg,cx,ca,cdgr,cddgr,cdg,alphsf'
      if(istdpr.gt.0) 
     *  write(istdpr,*) cg,cx,ca,cdgr,cddgr,cdg,alphsf
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'nout,cn,irsu,unew,iresa4   '
      if(istdpr.gt.0) 
     *  write(istdpr,*) nout,cn,irsu,unew,iresa4 
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'nmodel,kmodel,itsaml,ioldex'
      if(istdpr.gt.0) 
     *  write(istdpr,*) nmodel,kmodel,itsaml,ioldex
c   
c  zero model count 
c  
      if(nrd.gt.kmodel) then   
        nrd=0  
        rewind 2   
      end if   
c  
      nwrt=0   
c  
c  read model  
c  
    8 read(2,end=80,err=90) nmod,nh,data,(x(n),(aa(i,n),i=1,ivar),
     *  n=1,nh)
c  
      nrd=nrd+1
c  
c  test for correct number 
c  
      if(kmodel.gt.0.and.nrd.lt.kmodel) go to 8
c  
c  test for no model   
c  
      if(nh.le.0) then
        if(istdpr.gt.0) 
     *  write(istdpr,107) nh  
        go to 5  
      end if
c  
c  if itsaml = 1 test for conversion errors
c  
      if(itsaml.eq.1) call tstaml(x,aa,nh,iaa)   
c
c  test for resetting A4 near centre
c
      if(iresa4.eq.1) call rseta4(x,aa,nh,data,iaa)
c  
c  for icvzbn .ge. 1 add point at lower boundary of convective envelope
c  
      icvzb1=0 
c  
      if(icvzbn.ge.1) then
c  
        call stcvzb(x,aa,data,nh,iaa,ncbh) 
        if(ncbh.gt.0) then
c  
          icvzb1=icvzbn
c  
          ncbh1=ncbh-4 
          ncbh2=ncbh+4 
          xcb1=x(ncbh1)
          xcb2=x(ncbh2)
c  
c  reduce nn by 1 (to ensure that the final model has nn points)   
c  
          nn=nn-1  
        end if
      end if
c  
c  test for congruent mesh 
c  
      if(icnmsh.eq.1) then
c  
c  set xsi for congruent mesh  
c  
        cn=x(2)/(x(3)-x(2))+1
        xsi(1)=0 
        dxsi=(nn-1)/float(nh-1)  
        xsi(2)=cn+1+dxsi 
        do 12 n=1,nh
        if(n.gt.2) xsi(n)=xsi(n-1)+dxsi 
   12   aa(ix,n)=x(n) 
        go to 30
      end if
c
c  reset u?
c
      if(irsu.gt.0) aa(5,2)=unew   
c  
      if(nout.gt.0) then
        nd=max0(1,nh/nout)   
        if(istdpr.gt.0) 
     *  write(istdpr,109) 
      else
	nd=0
      end if
c  
c  set integrand   
c  
      fnr=1./(1+cg)
      cgr=fnr*cg   
      car=fnr*ca   
c  
      kdisc=0
      ndsp=-1
      do 20 n=1,nh 
      qx=aa(1,n)   
      if(x(n).eq.0) then
        php=fnr*data(5)/qx   
        pha=0
        phg=cgr*(data(6)-data(5))*qx 
        pdgr=0   
        pddgr=0   
c
c  Fixed 28/1/02, by adding cx in bracket.
c  
        dxsi=sqrt(php+abs(phg)+cx)
c
      else
c  
        x2=x(n)*x(n) 
        if(aa(2,n).ne.0) then
          aar2=aa(2,n)
          aar2=aar2/(alphsf*aar2+1)
          aar4=abs(aa(4,n))
          aar4=aar4/(alphsf*aar4+1)
c
c  to avoid problems in region of composition discontinuity,
c  limit magnitude of this term
c
	  aar4=min(aar4,10.*aar2)
        else
          aar2=1./alphsf
          aar4=1./alphsf
        end if
c
        php=fnr*aar2/(qx*x2)  
        pha=car*aar2*aar2/qx   
        phg=cgr*aar4*qx/x2
c  
        if(aa(4,n).lt.0.and.cdgr.gt.0) then
          pdgr=cdgr*aa(4,n)*aa(4,n)   
        else
          pdgr=0   
        end if
c  
        if(n.eq.1.or.cddgr.le.0) then  
          pddgr=0
        else 
          n1=n-1 
          pddgr=1.e-5*(aa(4,n)-aa(4,n1))/(x(n)-x(n1))
          pddgr=cddgr*pddgr*pddgr
        end if   
c
c  term from gradient in A4. To avoid unreasonable behaviour
c  limit basic term (before multiplication by cdg) to be
c  at most equal to remaining terms relevant to interior
c
c  Note: as a hack, do not limit, when cdg .lt. 0.
c
        if(n.eq.1.or.cdg.eq.0) then  
          pdg=0
        else 
	  partsum=php+pha+abs(phg)+cx
          n1=n-1 
          daa4=(aa(4,n)-aa(4,n1))/(x(n)-x(n1))
	  if(cdg.gt.0) then
            pdg=(1-x(n))*(1-x(n))*daa4
	    pdg=min(pdg*pdg,partsum)
          else
            pdg=(1-x(n))*(1-x(n))*daa4
	    pdg=pdg*pdg
	  end if
          pdg=abs(cdg)*pdg
        end if   
c
	dxsi=sqrt(php+pha+abs(phg)+pdgr+pddgr+pdg+cx)
c
c  set up storage indices for discontinuity
c
        if(ndisc.gt.0.and.n.gt.1) then
c
c  using old criterion
c
	  alogp=2*log(x(n)*aa(1,n))+log(aa(5,n)/(aa(2,n)*aa(3,n)))
	  if(n.gt.2) then
	    dlgrhp=log(aa(5,n)/aa(5,n-1))/(alogp-alogpp)
	    if(dlgrhp.gt.dlgrmx) then
	      if(kdisc.eq.0.or.ndsp.eq.-1.or.n.gt.ndsp+10) then
	        kdisc=kdisc+1
	        nsdisc(1,kdisc)=n-1
	        xsdisc(1,kdisc)=x(n-1)
              end if
	      ndsp=n
	      nsdisc(2,kdisc)=n
	      xsdisc(2,kdisc)=x(n)
	      dxsi=xn(n-1)
            end if
          end if
	  alogpp=alogp
c
        else if(ndisc.lt.0.and.n.gt.1) then
c
c  using new criterion
c
          dlgrhp=(aa(4,n)+aa(2,n))/(aa(2,n)*aa(3,n)+1.e-15)
          if(dlgrhp.gt.dlgrmx) then
            if(kdisc.eq.0.or.ndsp.eq.-1.or.n.gt.ndsp+10) then
              kdisc=kdisc+1
              nsdisc(1,kdisc)=n-1
              xsdisc(1,kdisc)=x(n-1)
c
c  to make sure that interval is always defined by at least
c  two points, set upper limit here
c
              nsdisc(2,kdisc)=n
              xsdisc(2,kdisc)=x(n)
            end if
            ndsp=n
	    if(n.gt.nsdisc(2,kdisc)) then
              nsdisc(2,kdisc)=n
              xsdisc(2,kdisc)=x(n)
	    end if
            dxsi=xn(n-1)
          end if
c
        end if
      end if
c  
c  test for boundary of convection zone
c
      if(cacvzb.gt.0.and.n.gt.2) then
	if(aa(4,n).gt.1.e-4.and.aa(4,n-1).lt.1.e-4.or.
     *     aa(4,n).lt.1.e-4.and.aa(4,n-1).gt.1.e-4) then
	  dxsi=dxsi+cacvzb
	  xn(n-1)=xn(ni1)+cacvzb
	  if(istdpr.gt.0) 
     *  write(istdpr,*) 'Add ',cacvzb,' at x =',x(n-1),x(n)
        end if
      end if
c
      if(nout.gt.0) then
        if(mod(n-1,nd).eq.0.and.istdpr.gt.0) 
     *    write(istdpr,110) n,x(n),php,pha,phg,pdgr,pddgr,pdg,cx,dxsi
      end if
c
   20 xn(n)=dxsi   
c
c  test for smoothing
c
      if(nsmth.gt.1) then
	call rnmean(xn,xsi,nh,1,1,nsmth)
	do 22 n=1,nh
   22   xn(n)=xsi(n)
      end if
c
      if(ndisca.gt.0.and.kdisc.gt.0) then
	if(istdpr.gt.0) 
     *  write(istdpr,120) (k,(nsdisc(i,k),i=1,2),(xsdisc(i,k),i=1,2),
     *    k=1,kdisc)
c
c  reset A4 near discontinuity 
c
	call resdsc(x,aa,nh,data,nsdisc,kdisc,iaa)
      end if
c
c  set controls for integration and interpolation
c
      if(ndisca.le.0.or.kdisc.eq.0) then
	kinter=1
	ninter(1,1)=1
	ninter(2,1)=nh
      else
	kinter=kdisc+1
	ninter(1,1)=1
	ninter(2,kinter)=nh
	do 23 k=1,kdisc
	ninter(2,k)=nsdisc(1,k)
   23   ninter(1,k+1)=nsdisc(2,k)
        if(istdpr.gt.0) 
     *  write(istdpr,'(/''ninter(1 - 2,k):''/(i2,2i5))') 
     *    (k, (ninter(i,k),i=1,2),k=1,kinter)
      end if
c
c  integrate and renormalize   
c  
      if(ndisca.le.0.or.kdisc.eq.0) then 
        call vinta(x,xn,xsi,nh,1,1)  
        xsi(1)=0 
        fnr=(nn+cn-1)/xsi(nh)
        do 24 n=1,nh 
        an(1,n)=xsi(n)
        xsi(n)=1+fnr*xsi(n)  
   24   aa(ix,n)=x(n) 
	nnew(1,1)=1
	nnew(2,1)=nn
c..	write(77,'(0pf12.7,1p3e15.7)') (x(n),xn(n),an(1,n),xsi(n),n=1,nh)
      else
c
c  separate integration on each interval
c
	xsisum=0
	do 25 k=1,kinter
	n1=ninter(1,k)
	n2=ninter(2,k)
	nhk=n2-n1+1
        call vinta(x(n1),xn(n1),xsi(n1),nhk,1,1)  
c..	if(k.eq.1) then
c..	  write(6,*) 'n, x, integrand, integral:'
c..	  write(6,'(i5,0pf10.5,1p2e13.5)') (n,x(n),xn(n),xsi(n),n=n1,n2)
c..     end if
   25   xsisum=xsisum+xsi(n2)
c
c  now set ranges, and scale xsi
c
	nnrtot=nn-kdisc*ndisca
	fnrtot=float(nnrtot)/xsisum
	nnrsum=0
	ninit=1
	do 27 k=1,kinter
	n1=ninter(1,k)
	n2=ninter(2,k)
	if(k.lt.kinter) then
	  nnk=xsi(n2)*fnrtot
	  nnrsum=nnrsum+nnk
        else
	  nnk=nnrtot-nnrsum
        end if
	if(istdpr.gt.0) 
     *  write(istdpr,*) 'k, n1, n2, xsi(n2), nnk, nnrsum',
     *    k, n1, n2, xsi(n2), nnk, nnrsum
	if(k.eq.1) then
          fnr=(nnk+cn-1)/xsi(n2)
        else
          fnr=(nnk-1)/xsi(n2)
        end if
	do 26 n=n1,n2
        an(1,n)=xsi(n)
        xsi(n)=ninit+fnr*xsi(n)  
   26   aa(ix,n)=x(n) 
	nnew(1,k)=ninit
	nnew(2,k)=ninit+nnk-1
   27   ninit=ninit+ndisca+nnk
        if(istdpr.gt.0) write(istdpr,'(/''nnew(1 - 2,k):''/(i2,2i5))') 
     *    (k, (nnew(i,k),i=1,2), k=1,kinter)
      end if
c
      if(nd.eq.1) then
	if(istdpr.gt.0) 
     *  write(istdpr,130) (xsi(n),n=1,nh)
      end if
c
c  interpolate to new mesh 
c  *********************** 
c  
   30 do 31 i=1,5  
   31 an(i,1)=aa(i,1)  
      xn(1)=x(1)  
      nst=1
      nni=nh   
      ni=1 
c
c  reset aa(2,.) and aa(4,.) to interpolate only in
c  quantities going as constants at x = 0. 
c
      if(x(1).eq.0) then
c  
        n0=2
        xs=cn+1  
        aa(2,1)=data(5)  
        aa(4,1)=data(6)-data(5)  
c
      else
        n0=1
        xs=0
      end if
c  
      do 34 n=n0,nh  
      x2=x(n)*x(n) 
      aa(2,n)=aa(2,n)/x2   
   34 aa(4,n)=aa(4,n)/x2   
c
c  test for resetting aa(2), aa(4) and aa(5) at singular surface
c
      singsf=data(7).ge.0
      if(singsf) then
        nh1=nh-1
        do 36 n=1,nh1
        x1=1.-x(n)
        aa(2,n)=aa(2,n)*x1
        aa(4,n)=aa(4,n)*x1
        if(data(7).gt.0) then
          aa(5,n)=aa(5,n)/x1**data(7)
        end if
   36   continue
        aa(2,nh)=aa(2,nh1)
        do 37 j=4,ivar
   37   aa(j,nh)=aa(j,nh1)
      end if
c  
c  prepare for separate interpolation for icvzbn = 2   
c  
      if(icvzb1.eq.2) then 
        nncbh=nh+1-ncbh
        xscbh=xsi(ncbh)
      end if   
c  
      n1=0 
c
      do 50 kint=1,kinter
c
      if(kint.eq.1) then
        nw1=n0
      else
	nw1=nnew(1,kint)
	xs=nw1-1
      end if
      nw2=nnew(2,kint)
c
      ni=1
      n1=ninter(1,kint)
      n2=ninter(2,kint)
      nhk=n2-n1+1
      nst=n1
      nni=nhk
c
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'n1, n2, nhk, nst, nni, nw1, nw2, xs',
     *  n1, n2, nhk, nst, nni, nw1, nw2, xs
c  
      do 50 n=nw1,nw2
      xs=xs+1  
c  
c  test for separate interpolation for icvzbn = 2  
c  **** Note: this has not been incorporated with discontinuities
c  
      if(icvzb1.eq.2) then 
        if(xs.le.xscbh) then   
          nst=1
          nni=ncbh 
        else   
          if(nst.eq.1) ni=1
          nst=ncbh 
          nni=nncbh
        end if 
      end if   
c
      call lir(xs,xsi(nst),an(1,n),aa(1,nst),ix,iaa,nni,ni,inter) 
      ni=2
c
      xn(n)=an(ix,n)
c
c  test that resulting xn is in fact monotonically increasing
c
      if(n.gt.nw1.and.xn(n).le.xn(n-1)) then
	if(istdpr.gt.0) 
     *  write(istdpr,135) n, xn(n), xn(n-1)
c
c  try linear interpolation in this neighbourhood instead
c
	ni=1
        call lir1(xs-1,xsi(nst),an(1,n-1),aa(1,nst),ix,iaa,nni,ni,inter) 
        call lir1(xs,xsi(nst),an(1,n),aa(1,nst),ix,iaa,nni,ni,inter) 
	xn(n-1)=an(ix,n-1)
        xn(n)=an(ix,n)
	if(istdpr.gt.0) 
     *  write(istdpr,140) n-1,xn(n-1),n,xn(n)
	if(xn(n).le.xn(n-1)) stop
      end if
c  
c  if icvzbn = 1 reset aa(4,.) from linear interpolation around
c  lower edge of convective envelope   
c  
      if(icvzb1.eq.1.and.xn(n).ge.xcb1.and.xn(n).le.xcb2) then
c  
c  locate interval in original model   
c  
   38   if(xsi(n1).le.xs.and.xs.lt.xsi(n1+1)) go to 39 
        n1=n1+1  
        go to 38   
c  
c  interpolate 
c  
   39   fct2=(xs-xsi(n1))/(xsi(n1+1)-xsi(n1))
        fct1=1-fct2  
        an41=fct1*aa(4,n1)+fct2*aa(4,n1+1)   
c  
        an(4,n)=an41 
c
      end if
c
   50 continue
c
c  test for setting variables in vicinity of discontinuities
c
      if(ndisca.gt.0.and.kdisc.gt.0) then
	do 60 k=1,kdisc
	n1=nnew(2,k)
	n2=nnew(1,k+1)
	x1=xn(n1)
	x2=xn(n2)
c
c  set ranges, steps in discontinuity region
c
        dlxds1=0.5*(x2-x1)/x2
	if(dlxdsc.gt.dlxds1) then
	  if(istdpr.gt.0) 
     *  write(istdpr,145) dlxdsc, dlxds1
        else
	  dlxds1=dlxdsc
        end if
c
	dxdsc=x2*dlxds1
	xn(n1+1)=x1+1.e-3*dxdsc
	xn(n1+3)=x1+0.5*(x2-x1-dxdsc)
	xn(n1+2)=xn(n1+3)-1.e-3*dxdsc
	xn(n1+4)=xn(n1+3)+1.e-3*dxdsc
	xn(n1+5)=xn(n1+3)+2.e-3*dxdsc
	xn(n2-1)=x2-1.e-3*dxdsc
	xn(n2-3)=x2-0.5*(x2-x1-dxdsc)
	xn(n2-2)=xn(n2-3)+1.e-3*dxdsc
	xn(n2-4)=xn(n2-3)-1.e-3*dxdsc
	xn(n2-5)=xn(n2-3)-2.e-3*dxdsc
	dxst=dxdsc/float(ndisca-9)
	do 52 n=n1+6,n2-6
   52   xn(n)=xn(n1+2)+dxst*(n-n1-5)
c
c  extrapolate to points just outside discontinuity region
c
	ns1=ninter(1,k)
	nn1=ninter(2,k)-ns1+1
	ns2=ninter(1,k+1)
	nn2=ninter(2,k+1)-ns2+1
c
	do 54 n=n1+1,n1+3
	if(istdpr.gt.0) 
     *  write(istdpr,*) 'Resetting through extrapolation at x =',xn(n)
   54   call lir1(xn(n),x(ns1),an(1,n),aa(1,ns1),ix,iaa,nn1,1,inter)
c
	do 55 n=n2-3,n2-1
	if(istdpr.gt.0) 
     *  write(istdpr,*) 'Resetting through interpolation at x =',xn(n)
   55   call lir1(xn(n),x(ns2),an(1,n),aa(1,ns2),ix,iaa,nn2,1,inter)
c
c  set on discontinuity region. Interpolate linearly, set
c  d ln rho/d ln x to a constant (remembering that we are
c  working in terms of scaled variables here)
c
	nd1=n1+3
	nd2=n2-3
	nndisc(1,k)=nd1
	nndisc(2,k)=nd2
	dlrdlx=log(an(1,nd2)*an(5,nd2)/(an(1,nd1)*an(5,nd1)))/
     *          log(xn(nd2)/xn(nd1))
	if(istdpr.gt.0) 
     *  write(istdpr,*) 'dlrdlx =',dlrdlx
	do 56 i=1,ivar
   56   danslp(i)=(an(i,nd2)-an(i,nd1))/(xn(nd2)-xn(nd1))
c
	do 58 n=nd1+1,nd2-1
	do 57 i=1,ivar
   57   an(i,n)=an(i,nd1)+danslp(i)*(xn(n)-xn(nd1))
   58   an(4,n)=-an(2,n)-dlrdlx/(xn(n)*xn(n))
c
   60   continue
c
      end if
c  
c  reset aa(2,.) and aa(4,.)   
c  
      do 62 n=n0,nn 
      x2=xn(n)*xn(n)   
      an(2,n)=x2*an(2,n)   
   62 an(4,n)=x2*an(4,n)   
c  
      do 65 n=n0,nh 
      x2=x(n)*x(n) 
      aa(2,n)=x2*aa(2,n)   
   65 aa(4,n)=x2*aa(4,n)   
c
c  reset new A_4 near discontinuity 
c
      if(ndisca.gt.0.and.kdisc.gt.0) then
	call resdsc(xn,an,nn,data,nndisc,kdisc,iaa)
      end if
c
c  test for resetting aa(2) and aa(4) at singular surface
c
      if(singsf) then
        nn1=nn-1
        do 72 n=1,nn1
        xn1=1.-xn(n)
        an(2,n)=an(2,n)/xn1
        an(4,n)=an(4,n)/xn1
        if(data(7).gt.0) then
          an(5,n)=an(5,n)*xn1**data(7)
        end if
   72   continue
        an(2,nn)=0
        an(4,nn)=0
        if(data(7).ne.0) an(5,nn)=0
c
        do 74 n=1,nh1
        x1=1.-x(n)
        aa(2,n)=aa(2,n)/x1
        aa(4,n)=aa(4,n)/x1
   74   aa(5,n)=aa(5,n)*x1**data(7)
        aa(2,nh)=0
        aa(4,nh)=0
        aa(5,nh)=0
c
      end if
c  
c  reset innermost points from expansion   
c  
      if(x(1).eq.0) call resexp(x,aa,nh,iaa,xn,an,nn,iaa,data,ioldex)
c  
c  test for setting point at convection zone boundary  
c  
      if(icvzb1.ge.1) call stcvzb(xn,an,data,nn,iaa,ncb) 
c  
      if(nout.gt.0) then
        ndn=max0(1,nn/nout)  
        if(istdpr.gt.0) 
     *  write(istdpr,160) data,(n,xn(n),(an(i,n),i=1,5),n=1,nn,ndn)   
      end if
c
c  output reset model
c
      write(3) nmod,nn,data,(xn(n),(an(i,n),i=1,ivar),n=1,nn) 
c  
c  
c  test for reading next model 
c  
      nwrt=nwrt+1  
      if(nwrt.lt.nmodel) go to 8
c  
      go to 5  
c  
c  diagnostics for end of file on model file   
c  
   80 if(istdpr.gt.0) 
     *  write(istdpr,180) 
c  
   90 continue 
      stop 
  102 format(/' **** Warning: using old discontinuity criterion,',
     *  ' ndisc =',i4)
  103 format(/' Using new discontinuity criterion, ndisc =',i4)
  107 format(//' **** error in redistrb. nh =',i10)
  109 format(///' n,x,php,pha,phg,pdgr,pddgr,pdg,cx,dxsi:'/)   
  110 format(i5,f10.7,1p8e13.5)
  120 format(//' Discontinuities in rho. k, n1, n2, x1, x2:'/
     *  (3i7,1p2e13.5))
  122 format(/' Start setting function at singularity no.',i4/
     *  ' n, x, xsi:')
  124 format(i5,1p2e13.5)
  130 format(/' xsi:'/(10f8.2))
  135 format(//' ***** error in redistrb at n =',i6/
     *         '       xn(n) = ',f12.7,'  .le. xn(n-1) =',
     *         f12.7)
  140 format(/' Try resetting with linear interpolation. n, x:'/
     *  (i5,f10.5))
  145 format(/' ***** Warning: dlxdsc reduced from',1pe13.5,
     *  ' to',e13.5)
  160 format(///' model on new mesh.'//' data:',1p8e13.5// 
     *  ' n,x,aa(1-5):'//(i5,0pf10.7,1p5e13.5))

  180 format(//' ********** end of file on d/s 2') 
      end  
      subroutine resexp(x,aa,nh,iaa,xn,an,nn,ian,data,ioldex)  
c  
c  reset innermost points in interpolated model from expansion 
c  if ioldex .ne. 1, second derivatives of p and rho in data   
c  are used.   
c  
c  for consistency with previous version of redistribution 
c  programme, also tests that an(5,n) .le. 3. this should  
c  not have any effect on recent, consistent models.   
c  
c  original version: 9/1/1985  
c  
c  modified 17/5/1985, to reset an(2,.) from expansion 
c  of p and of other aa-s  
c  
      implicit double precision (a-h, o-z)
      dimension x(nh),aa(iaa,nh),xn(nn),an(ian,nn),data(8) 
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c  
      x2=x(2)*x(2) 
      x4=x2*x2 
c  
c  test for new expansion  
c  
      if(ioldex.ne.1) go to 20 
c  
      ca1=(aa(1,2)-aa(1,1))/x2 
      cc1=0
      ca4=aa(4,2)/x2   
      cc4=0
      ca5=(aa(5,2)-aa(5,1))/x2 
      cc5=0
c  
      go to 30 
c  
c  new expansion   
c  
   20 ca1=-0.3*aa(1,1)*data(6) 
      cc1=(aa(1,2)-aa(1,1)-ca1*x2)/x4  
      ca4=data(6)-data(5)  
      cc4=(aa(4,2)-ca4*x2)/x4  
      ca5=-0.6*data(6) 
      cc5=(aa(5,2)-aa(5,1)-ca5*x2)/x4  
c  
c  expansion for aa(3,.) (same in the two cases)   
c  
   30 ca3=(aa(3,2)-aa(3,1))/x2 
      cc3=0
c  
      if(istdpr.gt.0) 
     *  write(istdpr,100) ca1,ca3,ca4,ca5,cc1,cc3,cc4,cc5 
c  
      xo2=x(2) 
c  
      amr2=data(1)/(data(2)*data(2))   
      pfct=6.6732e-8*amr2**2/(16.*atan(1.))
      cc2=0.5*aa(3,1)*data(5)  
c  
      do 40 n=1,nn 
      xx=xn(n) 
      if(xx.ge.xo2) go to 40   
      x2=xx*xx 
      x4=x2*x2 
c  
      an(1,n)=aa(1,1)+ca1*x2+cc1*x4
      an(3,n)=aa(3,1)+ca3*x2   
      an(4,n)=        ca4*x2+cc4*x4
      an(5,n)=aa(5,1)+ca5*x2+cc5*x4
c  
      an(2,n)=pfct*x2*an(1,n)**2*an(5,n)/(an(3,n)*data(3)*(1-cc2*x2))  
c  
   40 continue 
c  
c  restrict u to be .le. 3 
c  
      do 50 n=1,nn 
      if(an(5,n)-3.d0.gt.1.e-10) then
        if(istdpr.gt.0) write(istdpr,110) n,xn(n),an(5,n) 
        an(5,n)=3
      end if
   50 continue 
c  
      return   
  100 format(//' coefficients in s/r resexp:'//
     *  ' ca1, ca3, ca4, ca5 =',1p4e14.6/  
     *  ' cc1, cc3, cc4, cc5 =',4e14.6)
  110 format(' **** warning. at n =',i5,'  x =',f12.7, 
     *  '   aa(5,n) =',1pe15.7,' .gt. 3'/  
     *  '      has been reset to 3.')  
      end  
      subroutine stcvzb(x,aa,data,nn,iaa,ncb1) 
c  
c  sets mesh point at the botton edge of the convective envelope   
c  for adiabatic oscilltion model. 
c  position of edge is found by extrapolation in aa(4,.) from  
c  radiative interior. values of aa(i,.) at this point is  
c  found from linear interpolation between neighbouring points 
c  for i = 1, 2, 3 and 5, and aa(4,.) is set to zero.  
c  
c  this was introduced 21/5/1985. before that date extrapolation   
c  from radiative interior was used for all variables. 
c  
c  if extrapolated edge is outside or too close to next mesh   
c  point, this is taken to be bottom of convection zone,   
c  and no point is added   
c  
c  ncb1 returns position of bottom edge, or 0 if no edge is found  
c  
c  
      implicit double precision (a-h, o-z)
      dimension x(1),aa(iaa,1),data(1) 
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c  
c  locate bottom of convective envelope (to allow for convective overshoot,
c  define by small positive number)
c  
      ncb=0
      epscb=1.e-4
      do 10 n=2,nn 
      ncb=n
      if(aa(4,n).ge.epscb.and.aa(4,n+1).lt.epscb) go to 20 
   10 continue 
c  
c  no boundary found. write diagnostics
c  
      if(istdpr.gt.0) write(istdpr,100) 
      ncb1=0   
      return   
c  
c  extrapolate linearly to zero in aa(4,.) 
c  
   20 ncp=ncb-1
      fcp=aa(4,ncb)/(aa(4,ncb)-aa(4,ncp))  
      fcb=1.-fcp   
c  
c  set x at boundary, check for location   
c  
      ncb1=ncb+1   
      xcb=fcp*x(ncp)+fcb*x(ncb)
      if(xcb.lt.0.99*x(ncb1)+0.01*x(ncb)) go to 23 
c  
      if(istdpr.gt.0) write(istdpr,105) xcb,x(ncb1) 
      nn1=nn   
      go to 40 
c  
c  move x and aa   
c  
   23 nns=nn-ncb   
      nn1=nn+1 
      do 25 n=1,nns
      n1=nn1-n 
      n2=n1+1  
      x(n2)=x(n1)  
      do 25 i=1,5  
   25 aa(i,n2)=aa(i,n1)
c  
c  set interpolation coefficients  
c  
      fcb1=(xcb-x(ncb))/(x(ncb1)-x(ncb))   
      fcb=1-fcb1   
c  
c  set point at boundary   
c  
      x(ncb1)=xcb  
      do 30 i=1,5  
   30 aa(i,ncb1)=fcb*aa(i,ncb)+fcb1*aa(i,ncb1) 
c  
      aa(4,ncb1)=0 
c  
c  diagnostic output   
c  
   40 n1=ncb1-10   
      n2=ncb1+10   
      if(istdpr.gt.0) 
     *  write(istdpr,110) (n,x(n),(aa(i,n),i=1,5),n=n1,n2)
c  
      nn=nn1   
c  
      return   
  100 format(///' ********* no lower edge of convective envelope ',
     *  'found in s/r stcvzb') 
  105 format(///' extrapolated x at convection zone boundary =',   
     *  f10.5,'  is outside or too close to next mesh point =',
     *  f10.5/' no point added.')  
  110 format(///' points near lower boundary of convective envelope.'  
     *  //' n, x, aa(1-5):'//(i5,0pf12.7,1p5e15.7))
      end  
      subroutine tstaml(x,aa,nn,iaa)   
c  
c  test for conversion error leading to zero instead   
c  of 1
c  
      implicit double precision (a-h, o-z)
      dimension x(1),aa(iaa,1) 
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c  
      do 10 n=2,nn 
      if(x(n).ne.0) go to 10   
      if(istdpr.gt.0) write(istdpr,100) n   
      x(n)=1   
   10 continue 
      return   
  100 format(/' ********** conversion error for adiabatic',
     *  ' model. x(',i4,') = 0 has been reset to 1')   
      end  
      subroutine rseta4(x,aa,nn,data,iaa)   
c  
c  Reset A4 near centre, to correct for problems near
c  end of hydrogen burning
c  Resetting is only applied if A4/x**2 is non-monotonic
c  near centre
c  
      implicit double precision (a-h, o-z)
      dimension x(1),aa(iaa,1),data(1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  reference values
c
      axc=data(6)-data(5)
      nref=4
      axref=aa(4,nref)/x(nref)**2
c
      ireset=0
      do 10 n=2,nref-1
      ax=aa(4,n)/x(n)**2
      if((axc-ax)*(ax-axref).lt.0) ireset=1
   10 continue
c
      if(ireset.eq.1) then
	if(istdpr.gt.0) write(istdpr,110)
	axcoef=(axref-axc)/x(nref)**2
	do 20 n=2,nref-1
	x2=x(n)**2
	ax=aa(4,n)/x2
	axnew=axc+axcoef*x2
	if(istdpr.gt.0) write(istdpr,115) n, x(n), ax, axnew
   20   aa(4,n)=axnew*x2
c
      end if
c  
      return   
  110 format(//' Reset A4 near centre. ',
     *   ' n, x, old, new values of A4/x**2:'/)
  115 format(i5,1p3e13.5)
      end  
      function phidsc(x,x1,x2)
c
c  sets up function to be zero for x .le. x1 and 1 for  x .ge. x2,
c  monotomic with continuous first derivatives
c
      implicit double precision (a-h,o-z)
c
      if(x.lt.x1.or.x.gt.x2) then
	phidsc=0
      else
        xnorm=0.5*(x2*x2-x1*x1)*(x2+x1)-(x2**3-x1**3)/3.
     *    -x1*x2*(x2-x1)
        phidsc=(0.5*(x*x-x1*x1)*(x2+x1)-(x**3-x1**3)/3.
     *    -x1*x2*(x-x1))/xnorm
      end if
      return
      end
      subroutine resdsc(x,aa,nn,data,nsdisc,kdisc,iaa)
c  
c  Reset A4 near discontinuities in rho, to correct for problems
c  with numerical differentiation in evolution programme
c  
      implicit double precision (a-h, o-z)
      dimension x(1),aa(iaa,1),data(1),nsdisc(2,1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      do 30 k=1,kdisc
      n2=nsdisc(1,k)
      n1=max(n2-4,2)
      n3=nsdisc(2,k)
      n4=min(n3+4,nn-1)
      if(istdpr.gt.0) write(istdpr,110) k
c
      do 15 n=n1,n2
      if(n.eq.n1) then
	dlrhp=-(aa(2,n-1)+aa(4,n-1))
      else
	dlrhp=dlrh
      end if
      dlrh=2.*log(aa(1,n)*aa(5,n)/(aa(1,n-1)*aa(5,n-1)))/
     *        log(x(n)/x(n-1)) - dlrhp
      aa4new=-aa(2,n)-dlrh
      if(x(n)-x(n-1).lt.1.e-5*x(n)) aa4new=aa(4,n-1)
      if(istdpr.gt.0) write(istdpr,120) n, x(n), aa(4,n),aa4new
   15 aa(4,n)=aa4new

c
      do 20 n=n4,n3,-1
      if(n.eq.n4) then
	dlrhp=-(aa(2,n+1)+aa(4,n+1))
      else
	dlrhp=dlrh
      end if
      dlrh=2.*log(aa(1,n)*aa(5,n)/(aa(1,n+1)*aa(5,n+1)))/
     *        log(x(n)/x(n+1)) - dlrhp
      aa4new=-aa(2,n)-dlrh
      if(x(n+1)-x(n).lt.1.e-5*x(n)) aa4new=aa(4,n+1)
      if(istdpr.gt.0) write(istdpr,120) n, x(n), aa(4,n),aa4new
   20 aa(4,n)=aa4new
c
   30 continue
      return
  110 format(/' Reset A_4 at discontinuity no.',i3/
     *  ' n, x, old A_4, new A_4:')
  120 format(i5,1p3e13.5)
      end
