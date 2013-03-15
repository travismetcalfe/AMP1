      subroutine conmxc(x,y,iy,jc,nn,nv,ireset)
c
c  mixes material in convective core, making composition homogeneous.
c  abundances must be in y(4+i,n), i=1,nv. iy is first dimension of y.
c  Boundaries, interpolation coefficients are assumed to be given
c  in nf(jc), nl(jc), frcf(jc), frcl(jc).
c  However, unlike s/r conmix, this routine shifts outer boundary
c  of core to point of marginal stability with mixed composition.
c
c  If inc = 0 on input (so that convective regions have not
c  been set, or are non-existent), routine searches for convective
c  core and stores relevant parameters for inc = 1.
c
c  Original version: 10/8/90
c
c  Modified 24/8/90, to use trapezoidal rule integration to set
c  mixed composition (to avoid coupling to composition outside
c  the core, and setting ddacad with unmixed composition outside
c  core.
c
c  Modified for option of not resetting boundaries, unless they
c  have not been reset previously. This is determined by the
c  variable ireset, which is passed through the argument
c  (modification 16/10/90).
c
c  Note that when evolution of convective core abundances is
c  treated more or less correctly in programme, it is probably
c  not such a good idea to reset also in conmix, since
c  there are problems in the definition of the convection
c  zone boundary.
c
c  Modified 19/7/95, to use general heavy-element abundance
c  as set in common/heavy/. This has to be checked with care in
c  connection with convective mixing.
c
c  Modified 18/4/00, to drop resetting if central point is
c  found to be stable. Before some peculiar (and likely erroneous)
c  resetting based on convective-core boundary set in s/r rhs 
c  was used.
c
c  Modified 14/8/02, correcting storage from/to aztst and in cvr.
c  Prepare for inclusion of 4He burning.
c
c  Modified 30/5/03, correcting call of fdrad, replacing zh(n1) by zh(n)
c
c  Modified 22/8/03, changing diagnostics when unmixed composition
c  leads to instability (as expected at onset of hydrogen-burning
c  core convection)
c
c  Modified 11/11/03, including effect of composition gradient in
c  test for convection, for case with semiconvection (flagged by
c  idiffc1 .gt. 0)
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      implicit double precision (a-h,o-z)
c  Note: engenr.n.d.incl replaced by engenr.bz.d.incl, 14/8/02
      include 'engenr.bz.d.incl'
c
      dimension x(1),y(iy,1)
      dimension yr(20)
      common/anwvar/ data(8), yi(istrmx,1)
      common/xnwvar/ xi(1)
      common/heavy/ zatmos, zhc, zh(1)
      common/compvr/ cvr(icvrmx,1)
      common/cnvout/ ddrad,rrcnv,aacnv,ddacad,amach,
     *  drr(5), ddr(5), da(5),pturb,ycnv,dyda, dtxh
      common/cstcvr/ icvh1, icvhe4, icvhe3, icvc12, icvc13, icvn14, 
     *  icvo16
      common/cstyst/ iyh1, iyhe3, iyc12, iyc13, iyn14, iyo16, 
     *  iyche4, iycc12
      common/cstazt/ iah1, iahe3, iac12, iac13, ian14, iao16, 
     *  iache4, iacc12
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec
      common/engcnt/ xhzlm1, xhzlm2, nxhzer, icnocs, nspec,
     *  iheccs, nspect
      common/eqstd/ xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp,gmm1(4)
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,jnc,frcf(6),frcl(6)
      common/cxhder/ xhder(2,nnmax)
      common/cdiffu/ idiffus, itbdif, rhob, rmxdif, tlfdif, dtldif,
     *  ctdfmx, rctdif, rctdf1, idfxbc, ismdif, nsmdif, idiffc1
      common/ksider/ al0,al2,aztst(5)
      common/ln10/ amm
      common/cmtime/ age
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      data iopnfl /0/
c
      if(istdpr.gt.0) then
	write(istdpr,*) 
     *    'Entering conmxc with nv =',nv,' ireset =',ireset
        write(istdpr,*) 'jc, inc =', jc, inc
      end if
c
c  test for convective regions not set
c
      if(inc.le.0) then
        nmn=1
        nlim=1
      else
        nmn=nf(jc)
        nlim=max0(1,nmn - 40)
      end if
      if(istdpr.gt.0) write(istdpr,*) 'nmn, nlim =',nmn, nlim
c
c  set limits for extra output
c
      if(inc.gt.0.and.nl(inc).eq.nn) then
        ncp1=nf(inc)-3
        ncp2=nf(inc)+3
	if(istdpr.gt.0) write(istdpr,*) 'n, log q, X'
      else
        ncp1=0
        ncp2=0
      end if
c
c  set integrands
c
      xi(1)=0
      if(ispxx3.eq.1) then
        ishft = 3
      else
        ishft = 2
      end if
      yi(1,1)=aztst(3)
      do 10 i=2,nv
   10 yi(i,1)=aztst(ishft+i)
c
      n1=1
      do 15 n=nn,nlim,-1
      n1=n1+1
      xi(n1)=10.d0**x(n)
c
      if(n.ge.ncp1.and.n.le.ncp2.and.istdpr.gt.0) then
	write(istdpr,14090) n, x(n), y(5,n)
14090   format(i5,1pe13.5,0pf11.6)
      end if
c
      do 15 i=1,nv
   15 yi(i,n1)=y(4+i,n)
c
      ntot=n1
c
c  step through to locate convection zone boundary
c
      n1=1
      do 22 i=1,nv
   22 yi(i+nv,1)=0
      nc=0
c
c  extra output
c
      if(ncp1.lt.ncp2) then
        if(istdpr.gt.0) write(istdpr,*) 
     *    ' n, log q, log T, log rho, X, log kappa, ddad'
      end if
      do 40 n=nn,nlim,-1
      n1=n1+1
      call store(y(1,n),yr,6)
c
c  set mixed composition up to this point, using trapezoidal 
c  integration
c
      do 32 i=1,nv
      yi(i+nv,n1)=yi(i+nv,n1-1)+
     *            0.5d0*(xi(n1)-xi(n1-1))*(yi(i,n1)+yi(i,n1-1))
   32 yr(4+i)=yi(i+nv,n1)/xi(n1)
c
      if(n.ge.ncp1.and.n.le.ncp2.and.istdpr.gt.0)
     *  write(istdpr,'(a, 2i5,1p4e13.5)') 
     *  ' n, n1, etc ', n, n1, xi(n1),xi(n1-1), yi(i,n1),yi(i,n1-1)
c
      ddadp=ddad
c
c  before 30/5/03 fdrad was called with zh(n1); that must be
c  obviously wrong.
c
      drad=fdrad(x(n),yr,zh(n),ak,akr,akt,akx) 
      if(idiffc1.eq.0) then
        ddad=drad-dad(1)
      else
	dtxh=fdtxh(x(nlim),y(1,nlim),y(5,nlim),xhder(2,nlim))
        ddad=drad-dad(1)-dtxh
      end if
c
c  test for output
c
      if(n.ge.ncp1.and.n.le.ncp2) then
	tl=yr(3)
	rhl=log10(rho(1))
	akl=log10(ak)
	xh=yr(5)
	if(istdpr.gt.0) then
	  write(istdpr,20090) n, x(n), tl, rhl, xh, akl, ddad
20090     format(i5,6f11.6)
          write(istdpr,*) 'yi:',(yi(i,n1),i=1,2*nv)
        end if
      end if
c
c  test for convection at innermost point
c
      if(n1.eq.2) then
        if(ddad.le.0) then
c
c  test whether convection zones have been set previously
c  otherwise quit
c
          if(inc.le.0) then
            return
          else
            write(istdou,105) ddad
            if(istdpr.gt.0.and.istdou.ne.istdpr) write(istdpr,105) ddad
	    return
          end if
        end if
c
c  test for change of sign
c
      else if(ddadp.ge.0.and.ddad.lt.0) then
c
c  boundary found. set interpolation coefficients etc
c  Repeat determination of ddad with unmixed composition
c
	if(istdpr.gt.0) write(istdpr,*) 'boundary found at n =',n
        call store(y(1,n),yr,6)
        drad1=fdrad(x(n),yr,zh(n),ak,akr,akt,akx) 
	ddad1=drad1 - dad(1)
c
c  test for output
c
        if(n.ge.ncp1.and.n.le.ncp2) then
	  tl=yr(3)
	  rhl=log10(rho(1))
	  akl=log10(ak)
	  xh=yr(5)
	  if(istdpr.gt.0) write(istdpr,20090) 
     *      n, x(n), tl, rhl, xh, akl, ddad1
        end if
c
c  test that with original composition layer is still stable.
c  otherwise use mixed value
c  Note that this is generally as expected (22/8/03)
c
c  BUT THE LOGICS OF THIS, ORIGINALLY INVERTED, SEEMED RATHER ODD
c
	if(istdpr.gt.0) write(istdpr,*) 'ddad1 =',ddad1
	if(ddad1.ge.0) then
	  if(istdpr.gt.0) write(istdpr,112) x(n), tl, rhl, xh
        else
	  ddad=ddad1
        end if
c
        nc=n+1
	if(istdpr.gt.0) write(istdpr,*) 'ddad, ddadp', ddad, ddadp
        fc=ddad/(ddad-ddadp)
        go to 45
      end if
c
   40 continue
c
c  write diagnostics for no boundary found. 
c  For simplicity, use boundary from s/r rhs.
c
      write(istdou,110) 
      if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,110) 
c
   42 write(istdou,115) 
      if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,115) 
      nc=nmn
      fc=frcf(jc)
      go to 50
c
c  reset boundary quantities in convpt
c
   45 xcprev=frcf(jc)*x(nmn)+(1-frcf(jc))*x(nmn-1)
      xcnew =fc*x(nc)+(1-fc)*x(nc-1)
      if(istdpr.gt.0) then
	write(istdpr,*) 'nc, fc, x(nc), 1-fc, x(nc-1)',
     *    nc, fc, x(nc), 1-fc, x(nc-1)
        write(istdpr,120) nmn, frcf(jc), xcprev, nc, fc, xcnew
      end if
c
c  test for new setting of convective core
c
      if(inc.le.0) then
        inc=1
        jc=1
        nl(jc)=nn
        frcl(jc)=1
        nf(jc)=nc
        frcf(jc)=fc
      else if(ireset.eq.1) then
        nf(jc)=nc
        frcf(jc)=fc
      end if
c
      nmn=nc
c
c  reset composition
c
   50 qc=frcf(jc)*10.d0**x(nmn)+(1-frcf(jc))*10.d0**x(nmn-1)
      if(ireset.eq.1) then
        nc1=nn-nmn+2
	if(istdpr.gt.0) write(istdpr,*) 
     *    'Resetting composition: nc1, qc =',nc1, qc
        do 60 i=1,nv
        yc=(fc*yi(i+nv,nc1)+(1-fc)*yi(i+nv,nc1+1))/qc
	if(istdpr.gt.0) write(istdpr,*) 
     *    'i, fc, yi(i+nv,nc1),yi(i+nv,nc1+1)), yc =',
     *    i, fc, yi(i+nv,nc1),yi(i+nv,nc1+1), yc
	if(i.eq.1.and.yc.gt.y(5,1)) then
	  if(istdpr.gt.0) write(istdpr,152) yc,yc-y(5,1)
	  yc=y(5,1)
        end if
        do 60 n=nmn,nn
        y(4+i,n)=yc
c
c#ai#  Storage in cvr is still not very elegant and might need
c      fixing.
c
        if(i.eq.1) then
	  cvr(1,n)=yc
	  cvr(2,n)=1-yc-zh(n)
        else if(i.eq.iyche4-4) then
          cvr(icvhe4,n)=yc
        else if(i.eq.iycc12-4) then
          cvr(icvc12,n)=yc
        else if(i.eq.nv) then
          call setcno(y(1,n),n)
	end if
   60   continue
c
        aztst(3)=y(5,nn)
        do 70 i=2,nv
   70   aztst(ishft+i)=y(4+i,nn)
c
      else
	if(istdpr.gt.0) write(istdpr,160)
      end if
c
      return
  102 format('# entering conmxc. age =',1pe15.7/'#'/
     *  '# n, q, original X, mixed X'/'#')
  105 format(//' ***** Error in s/r conmxc. At innermost point ddad =',
     *  1pe13.5,' is negative'/
     *         '       Skip resetting of composition')
  110 format(//' ***** Error in s/r conmxc. No boundary found.')
  112 format(//' *** Problem in conmxc. Point stable with mixed ',
     *  'composition, unstable with unmixed composition'/
     *  ' log q =',f11.6,' log T =',f11.6,' log rho =',f11.6,
     *  ' unmixed X =',f11.6)
  115 format(/'       Use boundary set in s/r rhs')
  120 format(//' Outer boundary of convective core reset in conmxc.'/
     *  ' Old n, frc, xc =',i5,0pf10.5,1pe13.5/
     *  ' New n, frc, xc =',i5,0pf10.5,1pe13.5)
  150 format(/' resetting X(',i1,'). n, log q, Xi, Xbari'/
     *  (i5,1p3e13.5))
  152 format(/' ***** Warning in s/r conmxc:'/
     *        '       Xc = ',f15.8,' exceeds surface value by',1pe13.5/
     *        '       Xc reset to surface value')
  160 format(' *** Composition not reset in s/r conmxc')
      end
