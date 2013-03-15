      subroutine srdist(i_paramset, ierr_param, i_inout, 
     *  x, aa, data, xn, an, nh, nn, ivar, iaa, ian)
c
c  Redistributes mesh in model for adiabatic pulsations
c
c  Modified 22/9/08, introducing s/r test_mono to test monotonicity of 
c  xsi.
c
      implicit double precision (a-h, o-z)
      include '../adipr.incl'
      parameter (iax = 10)
      logical singsf
      dimension x(*),aa(iaa,*),xn(nnmax),an(ian,*),data(8)
      dimension xsi(nnmax),  nsdisc(2,100), xsdisc(2,100),
     *  ninter(2,10),nnew(2,10),nndisc(2,10),danslp(iax)
      common/comgrp/ isprtp, irotcp, omgrtp(1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      common/cstdio_in/ istdin_in, istdpr_in
      common/cstdio_def/ istdin_def, istdou_def, istdpr_def, istder_def
c
      data idsinp, idsoup /0, 0/
      data init_paramset /0/
c
      save
c
c..      namelist/exec/ nn,icnmsh,
c..     *  icase,icvzbn,nsmth,ndisc,dlxdsc,
c..     *  cg,cx,ca,cdgr,cddgr,cdg,alphsf
c..     *  nout,cn,irsu,unew,iresa4,   
c..     *  nmodel,kmodel,itsaml,ioldex
c
      write(istdou,'(//
     *  '' --------------------------------------------------''/
     *  '' Entering srdist''/)')
      if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,'(//
     *  '' --------------------------------------------------''/
     *  '' Entering srdist''/)')
      ierr_param=0
c
c  test for resetting printed output
c  A fairly preliminary reset, probably to be improved
c
      if(i_paramset.eq.1.and.istdpr.ne.istdpr_def.and.istdpr.gt.0) then
        write(istdpr,'(//'' Change istdpr to'',i3)') istdpr
	write(istdou,'(/'' Opening standard print on unit'',i3,
     *    '' file '',a40)') istdpr, 'ttt.redistrb.out'
        open(istdpr,file='ttt.redistrb.out',status='unknown')
      end if
c
c  test for skipping parameter input
c
      if(i_paramset.eq.1.and.init_paramset.ne.0) go to 6
c
c  defaults in exec 
c  nn: number of points in new mesh
      nn=601
      icnmsh=0  
c  icase: for icase = 1, 2, 3, 11, 12, 21, 22, 23, 24 set standard parameters.
c     icase = 1, 2, 3 give old parameters, as used up until June 1988.
c     icase = 11 gives new, optimized p mode parameters
c     icase = 1  corresponds to p modes, old parameters
c     icase = 2  corresponds to g modes, old parameters
c     icase = 3  corresponds to intermediate modes, old parameters.
c     icase = 11 corresponds to p modes, new parameters.
c     icase = 12 corresponds to g modes, new parameters.
c  Revision 3/5/10, including cg .lt. 0 (see below).
c  For icase .ge. 20, set new standard nsmth, ndisc, dlxdsc
c     icase = 21  corresponds to p modes, revised parameters
c     icase = 22  corresponds to g modes, revised parameters
c     icase = 23  corresponds to intermediate modes, revised parameters.
c     icase = 24  corresponds to modes in red giants
c
c  Note that cg .lt. 0 is used to flag for ensuring that a reasonable
c  fraction of points correspond to p modes (fraction hardcoded in
c  call of s/r testcg).
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
c  cg: abs(cg) weight for g-mode part of mesh
c     if cg .lt. 0 ensure a minimum numbers of meshpoint in p-mode region
c     The fraction is hard-coded to 0.1 in ppfrac below.
      cg=1  
c  cx: constant part of mesh function
      cx=40.
c  ca: weight to near-surface term
      ca=0.01   
c  cdgr, cddgr: additional weight to term depending on A4 (of somewhat
c     uncertain significance for now)
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
c
c  adda4: scaling of term to add to A(4,.)/x**2 in g-mode term
c     Actual term added is adda4*A(4,nr)/x(nr)**2, where nr-2 is
c     the innermost convectively stable point.
      adda4=0.d0
c  accrmn: Minimum value of A(4,.) used in cg term (set .gt. 0 to
c     increase number of points in convective core
      accrmn=0
c  nout: number of meshpoints in printed output
      nout=100  
c  cn: control of innermost non-zero meshpoint
      cn=2
c  irsu, unew: if irsu .gt. 0, reset U at innermost non-zero meshpoint to
c  unew
      irsu=0
      unew=3
c  iresa4: if iresa4 = 1, reset A4 before redistribution,
c  to correct possible errors in A4 resulting from near exhaustion
c  of hydrogen
c  if iresa4 = -1, reset data(6) on input, to correct for apparent
c  errors in setting amdl file (added 7/3/08)
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
c idsin, idsout: input and output units for file I/O with i_inout = 1
      idsin=2
      idsout=3
c   
c                   ......................................
c
c  set up files
c
      if(i_inout.eq.1) then
        if(istdpr.gt.0) then
          write(istdpr,*) 
     *      'Files needed: input  model on unit idsin (default 2)'
          write(istdpr,*) 
     *      '              output model on unit idsout (default 3)'
        end if
        call ofiles
      end if
c
c  read exec
    5 continue
c..      read(5,exec,end=90)   
      write(istdou,*) 'nn,icnmsh?'
      write(istdou,*) nn,icnmsh
      read(istdin,*,end=90,err=90) nn,icnmsh
      write(istdou,*) 'icase,icvzbn,nsmth,ndisc,dlxdsc,dlgrmx,cacvzb?'
      write(istdou,*) icase,icvzbn,nsmth,ndisc,dlxdsc,dlgrmx,cacvzb
      read(istdin,*) icase,icvzbn,nsmth,ndisc,dlxdsc,dlgrmx,cacvzb
      write(istdou,*) 'cg,cx,ca,cdgr,cddgr,cdg,alphsf,adda4,accrmn?'
      write(istdou,*) cg,cx,ca,cdgr,cddgr,cdg,alphsf,adda4,accrmn
      read(istdin,*) cg,cx,ca,cdgr,cddgr,cdg,alphsf,adda4,accrmn
      write(istdou,*) 'nout,cn,irsu,unew,iresa4   ?'
      write(istdou,*) nout,cn,irsu,unew,iresa4   
      read(istdin,*) nout,cn,irsu,unew,iresa4   
      write(istdou,*) 'nmodel,kmodel,itsaml,ioldex,idsin,idsout?'
      write(istdou,*) nmodel,kmodel,itsaml,ioldex,idsin,idsout
      read(istdin,*) nmodel,kmodel,itsaml,ioldex,idsin,idsout
c
c  store nn for later calls
c
      nn_res=nn
c
c  set flag for test of cg
c
      if(cg.lt.0) then
	itstcg = 1
	cg=abs(cg)
      else
	itstcg=0
      end if
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
      else if (icase.ge.20) then
c
c  Revised standard parameters, from 3/5/10
c
	nsmth=9
	ndisc=-60
	dlxdsc=0.0002
c
	if(icase.eq.21) then
c  Revised p modes
          cg=0.001  
          cx=1
          ca=0.005   
	  cdg=0.01
          cn=2  
        else if(icase.eq.22) then
c  Revised g modes (so far not changed)
          cg=4  
          cx=40 
          ca=0.01   
          cn=2  
        else if(icase.eq.23) then
c  Revised intermediate modes 
          cg=0.1  
          cx=1
          ca=0.005   
	  cdg=0.01
          cn=2  
c  Modes in red giants
        else if(icase.eq.24) then
          cg=-0.0001  
          cx=5
          ca=0.005   
	  cdg=0.01
          adda4=0.5
          cn=2  
        end if
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
	if(istdpr.gt.0) write(istdpr,103) ndisca
      else if(ndisc.gt.0) then
	ndisca=ndisc
	if(istdpr.gt.0) write(istdpr,102) ndisca
      else
	ndisca=0
      end if
c
c  test for resetting ndisc
c
      if(ndisca.gt.0.and.ndisca.lt.20) ndisca=20
c
c  ***************************************************************
c
c  In initial call with parameter setting, return after input of
c  parameters. Actual calculation is carried out with the subsequent
c  calls.
c
      if(i_paramset.eq.1.and.init_paramset.eq.0) then
c
        init_paramset = 1
        write(istdou,'(//'' End parameter input'')')
        if(istdpr.gt.0) write(istdpr,'(//'' End parameter input'')')
        ierr_param=1
        if(istdpr.gt.0.and.istdou.ne.istdpr) write(istdpr,'(//
     *    '' Exiting srdist''/
     *    '' --------------------------------------------------''/)')
        write(istdou,'(//
     *    '' Exiting srdist''/
     *    '' --------------------------------------------------'')')
        return
      end if
c
c  ***************************************************************
c
c  entry point after skipping parameter input
c
    6 continue
c
c  reset nn
c
      nn=nn_res
c
      if(istdpr.gt.0) then
        write(istdpr,*) 'nn,icnmsh'
        write(istdpr,*) nn,icnmsh
        write(istdpr,*) 'icase,icvzbn,nsmth,ndisc,dlxdsc,dlgrmx,cacvzb'
        write(istdpr,*) icase,icvzbn,nsmth,ndisc,dlxdsc,dlgrmx,cacvzb
        write(istdpr,*) 'cg,cx,ca,cdgr,cddgr,cdg,alphsf,adda4,accrmn'
        write(istdpr,*) cg,cx,ca,cdgr,cddgr,cdg,alphsf,adda4,accrmn
        write(istdpr,*) 'nout,cn,irsu,unew,iresa4   '
        write(istdpr,*) nout,cn,irsu,unew,iresa4 
        write(istdpr,*) 'nmodel,kmodel,itsaml,ioldex,idsin,idsout'
        write(istdpr,*) nmodel,kmodel,itsaml,ioldex,idsin,idsout
      end if
c
    7 if(i_inout.eq.1) then
c
c  possibly reset input and output files
c
	if(idsin.ne.idsinp) then
          nrd=0
          call openfc(idsin,idsinp,'o','u')
c
c  read start of file, to check for number of variables
c
          read(idsin,end=80,err=90) nmod,nh,data
          close(idsin)
          call openf(idsin,'o','u')
          idata8 = idint(data(8)+0.1)
          if(idata8.ge.100) then
            ivar = 8
          else if(idata8.ge.10) then
            ivar = 6
          else
            ivar = 5
          end if
c
	end if
c
        call openfc(idsout,idsoup,'u','u')
c   
c  zero model count 
c  
        if(nrd.gt.kmodel) then   
          nrd=0  
          rewind idsin   
        end if   
c  
        nwrt=0   
c  
c  read model  
c  
    8   read(idsin,end=80,err=90) nmod,nh,data,(x(n),(aa(i,n),i=1,ivar),
     *    n=1,nh)
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
          if(istdpr.gt.0) write(istdpr,107) nh  
          go to 5  
        end if
c  
c  if itsaml = 1 test for conversion errors
c  
        if(itsaml.eq.1) call tstaml(x,aa,nh,iaa)   
c
      end if
c
      ix=ivar+1
      if(isprtp.ne.0) then
	ixrot=ix+1
      else
	ixrot=ix
      end if
c
c  test for resetting A4 near centre
c
      if(iresa4.eq.1) call rseta4(x,aa,nh,data,iaa)
c
c  test for resetting second derivative of density
c
      if(iresa4.eq.-1) then
        dat6=(aa(4,2)+aa(2,2))/x(2)**2
        if(istdpr.gt.0) write(istdpr,'(//
     *    '' ***** Warning. data(6) reset from'',1pe13.5,
     *    '' to'',e13.5)') data(6),dat6
        data(6)=dat6
      end if
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
	aa(ix,n)=x(n) 
	if(isprtp.ne.0) aa(ixrot,n)=omgrtp(n)
   12   continue
        go to 30
      end if
c
c  reset u?
c
      if(irsu.gt.0) aa(5,2)=unew   
c
c  test for error in A(1,.) at centre
c
      if(x(1).eq.0.and.aa(1,1).le.aa(1,2)) then
	xs2=x(2)*x(2)
	xs3=x(3)*x(3)
	a1c=(xs3*aa(1,2)-xs2*aa(1,3))/(xs3-xs2)
	if(istdpr.gt.0) 
     *    write(istdpr,'(//'' ***** Warning. Central A1 changed from'',
     *    1pe13.5,'' to '',e13.5)') aa(1,1),a1c
	aa(1,1)=a1c
      end if
c
c  test for resetting cg to secure reasonable p-mode resolution
c
      if(itstcg.gt.0) then 
	ppfrac=0.1
	call testcg(x,aa,iaa,nh,ppfrac,cg)
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
c
c  test for setting term to add to A(4,.)/x**2
c
      if(adda4.gt.0) then
	nrcz=0
	do n=2,nh
	  if(aa(4,n).gt.0.and.nrcz.eq.0) nrcz=n
        end do
	cadda4=adda4*aa(4,nrcz+2)/x(nrcz+2)**2
	if(istdpr.gt.0) 
     *    write(istdpr,'(//'' Addition term to A(4,.)/x**2 set to'',
     *    1pe13.5)') cadda4
c
      else
	cadda4=0.d0
      end if
c
      if(nout.gt.0) then
        nd=max0(1,nh/nout)   
        if(istdpr.gt.0) write(istdpr,109) 
      else
	nd=0
      end if
c
c  start loop setting integrand
c
      do 20 n=1,nh 
      qx=aa(1,n)   
      if(x(n).eq.0) then
        php=fnr*data(5)/qx   
        pha=0
        phg=cgr*(data(6)-data(5)+cadda4)*qx 
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
          if(accrmn.gt.0) aar4=max(aar4,accrmn)
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
        phg=cgr*(aar4/x2+cadda4)*qx
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
c  Changed 16/2/05, to set always upper limit one point higher
c  than present point in discontinuity region.
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
              nsdisc(2,kdisc)=n+1
              xsdisc(2,kdisc)=x(n+1)
            end if
            ndsp=n
	    if(n+1.gt.nsdisc(2,kdisc)) then
              nsdisc(2,kdisc)=n+1
              xsdisc(2,kdisc)=x(n+1)
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
     *      write(istdpr,*) 'Add ',cacvzb,' at x =',x(n-1),x(n)
        end if
      end if
c
      if(nout.gt.0.and.istdpr.gt.0) then
        if(mod(n-1,nd).eq.0) 
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
     *    write(istdpr,120) (k,(nsdisc(i,k),i=1,2),(xsdisc(i,k),i=1,2),
     *    k=1,kdisc)
c
c  reset A4 near discontinuity 
c
	call resdsc(x,aa,nh,data,nsdisc,kdisc,iaa,ndisc)
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
     *    write(istdpr,'(/''ninter(1 - 2,k):''/(i2,2i5))') 
     *    (k, (ninter(i,k),i=1,2),k=1,kinter)
      end if
c
c  integrate and renormalize   
c  
      if(ndisca.le.0.or.kdisc.eq.0) then 
        call vinta(x,xn,xsi,nh,1,1)  
        call test_mono(x,xn,xsi,nh)  
        xsi(1)=0 
        fnr=(nn+cn-1)/xsi(nh)
        do 24 n=1,nh 
        an(1,n)=xsi(n)
        xsi(n)=1+fnr*xsi(n)  
	aa(ix,n)=x(n) 
	if(isprtp.ne.0) aa(ixrot,n)=omgrtp(n)
   24   continue
	nnew(1,1)=1
	nnew(2,1)=nn
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
        call test_mono(x(n1),xn(n1),xsi(n1),nhk)  
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
          if(nnk.le.1) then
            write(istdou,'(/'' ***** Warning: in srdist nnk ='',i3,
     *        '' for n1, n2 ='',2i5/
     *                      ''       nnk reset to 2'')') nnk, n1, n2
            if(istdpr.gt.0.and.istdpr.ne.istdou)
     *        write(istdpr,'(/'' ***** Warning: in srdist nnk ='',i3,
     *        '' for n1, n2 ='',2i5/
     *                      ''       nnk reset to 2'')') nnk, n1, n2
            nnk=2
          end if
	  nnrsum=nnrsum+nnk
        else
	  nnk=nnrtot-nnrsum
        end if
	if(istdpr.gt.0) 
     *    write(istdpr,*) 'k, n1, n2, xsi(n2), nnk, nnrsum',
     *    k, n1, n2, xsi(n2), nnk, nnrsum
	if(k.eq.1) then
          fnr=(nnk+cn-1)/xsi(n2)
        else
          fnr=(nnk-1)/xsi(n2)
        end if
	do 26 n=n1,n2
        an(1,n)=xsi(n)
        xsi(n)=ninit+fnr*xsi(n)  
	aa(ix,n)=x(n) 
	if(isprtp.ne.0) aa(ixrot,n)=omgrtp(n)
   26   continue
	nnew(1,k)=ninit
	nnew(2,k)=ninit+nnk-1
   27   ninit=ninit+ndisca+nnk
        if(istdpr.gt.0) write(istdpr,'(/''nnew(1 - 2,k):''/(i2,2i5))') 
     *    (k, (nnew(i,k),i=1,2), k=1,kinter)
      end if
c
      if(nd.eq.1.and.istdpr.gt.0) write(istdpr,130) (xsi(n),n=1,nh)
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
      if(istdpr.gt.0) write(istdpr,*) 
     *  'n1, n2, nhk, nst, nni, nw1, nw2, xs',
     *  n1, n2, nhk, nst, nni, nw1, nw2, xs
c
      if(istdpr.gt.0) then
        write(istdpr,*) 'xsi(nst) =',xsi(nst)
        write(istdpr,*) 'n, x, xsi, a4'
        write(istdpr,'(i5,0p2f10.5,1pe13.5)') 
     *    (n,x(n),xsi(n),aa(4,n)*x(n)*x(n),n=nst,nst+10)
      end if
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
      call lir(xs,xsi(nst),an(1,n),aa(1,nst),ixrot,iaa,nni,ni,inter) 
      ni=2
c
      xn(n)=an(ix,n)
      if(isprtp.ne.0) then
	omgrtp(n)=an(ixrot,n)
      else
	omgrtp(n)=0.d0
      end if
      if(n-nw1.le.5.and.istdpr.gt.0) 
     *  write(istdpr,*) n, xs, xn(n), an(4,n)*xn(n)*xn(n)
c
c  test that resulting xn is in fact monotonically increasing
c
      if(n.gt.nw1.and.xn(n).le.xn(n-1)) then
	if(istdpr.gt.0) write(istdpr,135) n, xn(n), xn(n-1)
c
c  try linear interpolation in this neighbourhood instead
c
	ni=1
        call lir1(xs-1,xsi(nst),an(1,n-1),aa(1,nst),ixrot,iaa,nni,
     *    ni,inter) 
        call lir1(xs,xsi(nst),an(1,n),aa(1,nst),ixrot,iaa,nni,ni,inter)
	xn(n-1)=an(ix,n-1)
        xn(n)=an(ix,n)
	if(istdpr.gt.0) write(istdpr,140) n-1,xn(n-1),n,xn(n)
	if(xn(n).le.xn(n-1)) then
	  write(istdou,142) 
	  if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdou,142) 
	  ierr_param=-1
	  go to 95
        end if
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
      if(istdpr.gt.0) write(istdpr,'(//'' End normal pass''/)')
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
	  if(istdpr.gt.0) write(istdpr,145) dlxdsc, dlxds1
        else
	  dlxds1=dlxdsc
        end if
c
	dxdsc=x2*dlxds1
	do n=1,5
	  xn(n1+n)=x1+n*1.e-3*dxdsc
	  xn(n2-n)=x2-n*1.e-3*dxdsc
	end do
	dxst=(xn(n2-5)-xn(n1+5))/float(n2-n1-10)
	do 52 n=n1+6,n2-6
   52   xn(n)=xn(n1+5)+dxst*(n-n1-5)
c
c  extrapolate to points just outside discontinuity region
c
	ns1=ninter(1,k)
	nn1=ninter(2,k)-ns1
	ns2=ninter(1,k+1)+1
	nn2=ninter(2,k+1)-ns2
c
	if(istdpr.gt.0) 
     *    write(istdpr,'(/'' Extrapolate from x range'',1p2e17.9)')
     *    x(ns1),x(ns1+nn1-1)
	do 54 n=nnew(1,k),n1+3
	if(n.gt.n1.or.xn(n).gt.x(ns1+nn1-1)) then
          call lir1(xn(n),x(ns1),an(1,n),aa(1,ns1),ixrot,iaa,
     *      nn1,1,inter)
          if(istdpr.gt.0) write(istdpr,*) 
     *      'Resetting at n, x , A4 =',n, xn(n), an(4,n)
        end if
   54   continue
c
	if(istdpr.gt.0) then
          write(istdpr,'(/'' Extrapolate from x range'',1p2e17.9)')
     *      x(ns2),x(ns2+nn2-1)
          write(istdpr,'('' x, A4 at first two points:''/
     *      (1pe17.9,e13.5))')
     *      x(ns2),aa(4,ns2),x(ns2+1),aa(4,ns2+1)
	end if
	do 55 n=n2-3,nnew(2,k+1)
	if(n.lt.n2-1.or.xn(n).lt.x(ns2+1)) then
          call lir1(xn(n),x(ns2),an(1,n),aa(1,ns2),ixrot,
     *      iaa,nn2,1,inter)
          if(istdpr.gt.0) 
     *      write(istdpr,*) 'Resetting at n, x, A4/x^2, A4 =',
     *      n, xn(n), an(4,n), an(4,n)*xn(n)**2
        end if
   55   continue
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
	if(istdpr.gt.0) write(istdpr,*) 'dlrdlx =',dlrdlx
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
c  Switched off 10/6/03
c
c..      if(ndisca.gt.0.and.kdisc.gt.0) then
c..	   call resdsc(xn,an,nn,data,nndisc,kdisc,iaa,ndisc)
c..      end if
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
      if(x(1).eq.0) call resexp(x,aa,nh,iaa,xn,an,nn,ian,data,ioldex)
c  
c  test for setting point at convection zone boundary  
c  
      if(icvzb1.ge.1) call stcvzb(xn,an,data,nn,ian,ncb) 
c  
      if(nout.gt.0.and.istdpr.gt.0) then
        ndn=max0(1,nn/nout)  
        write(istdpr,160) data,(n,xn(n),(an(i,n),i=1,5),n=1,nn,ndn)   
      end if
c
c  output reset model
c
      if(i_inout.eq.1) 
     *  write(idsout) nmod,nn,data,(xn(n),(an(i,n),i=1,ivar),n=1,nn) 
c
c  return if not reading and writing models
c
      if(i_inout.ne.1) then
	go to 95
c
      else
c  
c  test for reading next model 
c  
        nwrt=nwrt+1  
        if(nwrt.lt.nmodel) go to 7
c  
        go to 5  
      end if
c  
c  diagnostics for end of file on model file   
c  
   80 if(istdpr.gt.0) write(istdpr,180) 
c  
   90 continue 
      ierr_param=-1
c
c  common return point, to enable suitable closing of files etc.
c
   95 if(istdpr.gt.0.and.istdou.ne.istdpr) write(istdpr,'(//
     *  '' Exiting srdist''/
     *  '' --------------------------------------------------''/)')
      write(istdou,'(//
     *  '' Exiting srdist''/
     *  '' --------------------------------------------------'')')
      if(istdpr.ne.istdpr_def.and.istdpr.gt.0) then
        close(istdpr)
      end if
c
      return 
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
  142 format(/' ***** Error. Problems with non-monotonicity in sredist')
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
c  for adiabatic oscillation model. 
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
      if(x(n).eq.0) then
        if(istdpr.gt.0) write(istdpr,100) n   
        x(n)=1   
      end if
   10 continue 
      return   
  100 format(/' ********** conversion error for adiabatic',
     *  ' model. x(',i4,') = 0 has been reset to 1')   
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
      subroutine resdsc(x,aa,nn,data,nsdisc,kdisc,iaa,ndisc)
c  
c  Reset A4 near discontinuities in rho, to correct for problems
c  with numerical differentiation in evolution programme
c
c  Modified 16/2/05, to set simply d ln rho/d ln x constant over the
c  last 4 meshpoints around the discontinuity
c  
      implicit double precision (a-h, o-z)
      dimension x(1),aa(iaa,1),data(1),nsdisc(2,1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      do 50 k=1,kdisc
      n2=nsdisc(1,k)
      n1=max(n2-4,2)
      n3=nsdisc(2,k)
      n4=min(n3+4,nn-1)
      if(istdpr.gt.0) write(istdpr,110) k
c
c  test for using old or new formulation
c
      if(ndisc.gt.0) then
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
   15   aa(4,n)=aa4new

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
   20   aa(4,n)=aa4new
c
      else
c
c  use new formulation with constant d ln rho/d ln x
c
        dlrh=log(aa(1,n1)*aa(5,n1)/(aa(1,n2)*aa(5,n2)))/
     *        log(x(n1)/x(n2))
        do 25 n=n1,n2
        aa4new=-aa(2,n)-dlrh
        if(x(n)-x(n-1).lt.1.e-5*x(n)) aa4new=aa(4,n-1)
        if(istdpr.gt.0) write(istdpr,120) n, x(n), aa(4,n),aa4new
   25   aa(4,n)=aa4new

c
        dlrh=log(aa(1,n4)*aa(5,n4)/(aa(1,n3)*aa(5,n3)))/
     *        log(x(n4)/x(n3))
        do 30 n=n4,n3,-1
        aa4new=-aa(2,n)-dlrh
        if(x(n+1)-x(n).lt.1.e-5*x(n)) aa4new=aa(4,n+1)
        if(istdpr.gt.0) write(istdpr,120) n, x(n), aa(4,n),aa4new
   30   aa(4,n)=aa4new
c
      end if
c
   50 continue
      return
  110 format(/' Reset A_4 at discontinuity no.',i3/
     *  ' n, x, old A_4, new A_4:')
  120 format(i5,1pe17.9,2e13.5)
      end
      subroutine testcg(x,aa,iaa,nh,ppfrac,cg)
c
c  resets cg to ensure that roughly the fraction ppfrac of the points
c  are in the p-mode region
c
c  Original version: 9/3/05
c
      implicit double precision (a-h, o-z)
      include '../adipr.incl'
      parameter(iw = 4)
      dimension x(1), aa(iaa,1)
      dimension w(iw,nnmax)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      do n=2,nh
	w(1,n)=sqrt(aa(2,n)/(aa(1,n)*x(n)*x(n)))
	w(2,n)=sqrt(abs(aa(4,n))*aa(1,n)/(x(n)*x(n)))
      end do
      w(1,1)=w(1,2)
      w(2,1)=w(2,2)
c
      call vinta(x,w(1,1),w(3,1),nh,iw,iw)
      call vinta(x,w(2,1),w(4,1),nh,iw,iw)
c
c  test for limit on cg
c
      cgmax=(1-ppfrac)*w(3,nh)/(ppfrac*w(4,nh))
      cgmax=cgmax*cgmax
      if(istdpr.gt.0) then
        write(istdpr,*) 'w(3,nh),w(4,nh),cgmax,cg'
        write(istdpr,*)  w(3,nh),w(4,nh),cgmax,cg 
      end if
      if(cg.gt.cgmax) then
	if(istdpr.gt.0) write(istdpr,110) cg, cgmax, cgmax
	cg=cgmax
      end if
      return
      if(cg.gt.cgmax) then
	if(istdpr.gt.0) write(istdpr,110) cg, cgmax, cgmax
	cg=cgmax
      end if
      return
  110 format(//' ***** Warning: cg = ',1pe13.5,' gt cgmax =',e13.5/
     *         '                cg reset to', e13.5)
      end
      subroutine test_mono(x,xn,xsi,nh)  
c
c  test that xsi resulting from integration is monotonic. If not,
c  reset with trapezoidal integration
c
      implicit double precision (a-h, o-z)
      include '../adipr.incl'
      dimension x(*), xn(*), xsi(*)
      dimension xsi1(nnmax)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      n=2
      dxsi=0.d0
      ireset=0
c
   10 n=n+1
      if(xsi(n).le.xsi(n-1)) then
c
c  find region of monotonicity
c
	ns=n-1
	do n1=n,nh
	  if(xsi(n1).gt.xsi(ns)) then
	    nf=n1
	    go to 20
          end if
        end do
        write(istdou,'(/
     *    '' **** Warning. In test_mono no monotonic point found''/
     *    ''      Reset to surface''/)')
        if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,'(/
     *    '' **** Warning. In test_mono no monotonic point found''/
     *    ''      Reset to surface''/)')
        nf=nh
c
c  reset on this interval
c
   20   continue
        do n1=ns+1,nf
	  xsi1(n1)=xsi1(n1-1)+0.5d0*(x(n1)-x(n1-1))*(xn(n1-1)+xn(n1))
        end do
        ireset=1
	n=nf
	dxsi=xsi1(nf)-xsi(nf)
      else
	xsi1(n)=xsi(n)+dxsi
      end if
      if(n.lt.nh) go to 10
c
c  transfer result, diagnostic output
c
      if(ireset.eq.1.and.istdpr.gt.0) 
     *  write(istdpr,'(//'' n, x, old xsi, new xsi:''/)')
      do n=3,nh
	if(xsi1(n).ne.xsi(n)) 
     *    write(istdpr,'(i5,f12.7,1p2e13.5)') n, x(n), xsi(n), xsi1(n)
	xsi(n)=xsi1(n)
      end do
      return
      end
