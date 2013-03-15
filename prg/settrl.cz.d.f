      subroutine settrl(idstrl,itrial,xmdtrl,ipartr,nrdtr1,nvartr,
     *  idatmd,nspcm1,nbcpr1,x,y,datmod,ndtmod,bccoef,iy,iformr,nnt,
     *  nrdtmr,nidtmr,nbccfr,nspcmr,icnotr,ihectr,isprtt,icry)
c
c  Subroutine for setting trial solution.
c  *************************************
c
c  Also restores bccoef in form appropriate for current
c  value of nspcmx
c
c  Note: when reading single-precision data (for itrial .lt. 0)
c  assume old format, with nvartr .le. 6, 31 variables in 
c  datmod and 54 variables in bccoef.
c  datmdr and ndtmdr are only used for reading old format
c  and hence can be dimensioned to fit that.
c
c  Double-precision data may be read either on old or new
c  form, automatically decided by s/r rdemdl.
c
c  Original version: 13/8/91
c
c  Modified 13/7/05, entering with (possibly) non-integer model
c  number xmdtrl. If xmdtrl is not an integer, interpolate
c  linearly between model no. nmdtrl = int(xmdtrl) and nmdtrl+1
c
c  Modified 31/3/06, to include ipartr and isprtt in parameter list.
c  ipartr may force setting of rotational flag and angular velocity
c  if available and isprtt returns in rotation flag from the input
c  model.
c
c  Modified 20/7/06 to use internal storage of evolution sequence 
c  if iastr = -1 and isetos = 1
c
      implicit double precision(a-h, o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.cz.d.incl'
c
      logical mod_int
      real xrds, yrds, datrds,bccrds
      dimension x(*), y(iy,*), datmod(*), ndtmod(*), bccoef(*)
      dimension datrds(31), ndtrds(3), bccrds(54), datmdr(31),
     *  ndtmdr(3),bccfst(1), x2(nnmax), y2(ivarmx,nnmax), yint(ivarmx),
     *  datmod2(nrdtmd), ndtmod2(nidtmd), bccoef2(nbcprv)
      common/comgrt/ isprot, irotcn, a5, riner0, riner, omgrt0,
     *  omgrot(1)
      common/sooner/ xrds(1)
      common/work/ yrds(6,1)
      common/cevlio/ isetos, iastr
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      equivalence(datmdr(29), ndtmdr(1)),(datrds(29), ndtrds(1)),
     *  (bccfst(1),xrds(1),x2(1)),(yrds(1,1),y2(1,1))
c
      data nrdtrl /-1/
c
c  initialize icasrd so that first read tests for type of model
c
      data icasrd / 0/
c
      save 
c
      if(istdpr.gt.0) 
     *  write(istdpr,'(/'' Enter settrl with itrial, xmdtrl, nrdtr1 ='',
     *  i4,f10.5,i4/)') itrial, xmdtrl, nrdtr1
c
c  if nrdtr1 on input is negative, reset nrdtrl
c
      if(nrdtr1.lt.0) nrdtrl = -1
c
c  test for interpolation
c
      nmdtrl=xmdtrl+1.d-10
      if(abs(xmdtrl-nmdtrl).gt.1.d-7) then
	mod_int=.true.
	frc_int2=xmdtrl-nmdtrl
	frc_int1=1.d0-frc_int2
      else
	mod_int=.false.
      end if
c
c  test for getting model from internal storage
c
      if(isetos.eq.1.and.iastr.eq.-1) then
	call mread_emdl(nmdtrl,x,y,datmod,ndtmod,bccoef,iy,
     *    nnt,nrdtmr,nidtmr,nbccfr,nvartr,icry)
	if(icry.lt.0) return
	iformr=1
        go to 45
      end if
c
c  read model in standard (evolution) form from d/s idstrl
c
      if(nrdtrl.lt.0) then
        rewind idstrl
        nrdtrl=0
      end if
c
   20 if(nmdtrl.gt.0) then
c
        if(nmdtrl.lt.nrdtrl) then
          rewind idstrl
          nrdtrl=0
	else if(nmdtrl.eq.nrdtrl) then
	  go to  45
        end if
c
      end if
c
c  test for single-precision or double precision read
c
      if(itrial.lt.0) then
c
        if(idatmd.eq.0) then
          read(idstrl,end=40) nnt,(xrds(n),(yrds(j,n),j=1,nvartr),
     *      n=1,nnt),(bccrds(i),i=1,54)
        else
          read(idstrl,end=40) nnt,(datrds(i),i=1,31),(xrds(n),
     .      (yrds(j,n),j=1,nvartr),n=1,nnt),(bccrds(i),i=1,54)
	  do 26 i=1,29
   26     datmod(i)=datrds(i)
	  do 27 i=1,3
   27     ndtmod(i)=ndtrds(i)
        end if
c
        nrdtrl=nrdtrl+1
c
	do 28 i=1,54
   28   bccoef(i)=bccrds(i)
c
	do 29 n=1,nnt
	x(n)=xrds(n)
	do 29 i=1,nvartr
   29   y(i,n)=yrds(i,n)
c
c  set iformr = 0 to flag for old format
c
	iformr = 0
c
      else
c
        call rdemdl(idstrl,x,y,datmod,ndtmod,bccoef,iy,iformr,
     *    nnt,nrdtmr,nidtmr,nbccfr,nvartr,icasrd,icryrd)
c
	icasrd = -1
c
	if(icryrd.eq.-2) then
	  stop 'Error in s/r settrl'
        else if(icryrd.eq.-1) then
	  go to 40
        else if(icryrd.eq.1) then
c
c  file was rewound. Reset nrdtrl
c
	  nrdtrl = 1
        else
          nrdtrl=nrdtrl+1
        end if
c
      end if
c
      if(nmdtrl.gt.0) then
	go to 20
      else
	go to 45
      end if
c
   40 if(nmdtrl.lt.10000) then
c
c  error for end of file on trial data set
c
        write(istdou,103) nmdtrl,nrdtrl
        if(istdpr.gt.0.and.istdpr.ne.istdou) 
     *    write(istdpr,103) nmdtrl,nrdtrl
        rewind idstrl
        nrdtrl=0
	icry = -1
	nrdtr1=0
	return
c
      end if
c
c  test for setting storage parameters, for read of old format
c
   45 if(iformr.eq.0) then
c
        nrdtmr = 28
        nidtmr = 3
        nbccfr = 54
c
        icnotr = 0
	ihectr = 0
        nspcmr = 0
c
      else
c
c  set flag for CNO and 4He cases and parameter for composition storage
c
        icnotr = ndtmod(5)
        nspcmr = ndtmod(7)
	if(nidtmr.ge.80) then
	  ihectr = ndtmod(61)
        else
	  ihectr = 0
        end if
      end if
c
c  test for interpolating between two models in sequence.
c  If so, read next model or extract from internal storage
c
      if(mod_int) then
	if(isetos.eq.1.and.iastr.eq.-1) then
c
c  internal storage
c
	  call mread_emdl(nmdtrl+1,x2,y2,datmod2,ndtmod2,bccoef2,ivarmx,
     *      nnt,nrdtmr,nidtmr,nbccfr,nvartr,icryrd)
        else
c
c  read next
c
          call rdemdl(idstrl,x2,y2,datmod2,ndtmod2,bccoef2,ivarmx,
     *      iformr,nnt,nrdtmr,nidtmr,nbccfr,nvartr,icasrd,icryrd)
c
	end if
c
        if(icryrd.ne.0) then
	  write(istdou,'(//'' ***** Error in s/r settrl.''/
     *                     ''       No subsequent model for'',
     *      '' interpolation with xmdtrl ='',f10.5/)') xmdtrl
	  icry=-1
	  return
        end if
	nrdtrl=nrdtrl+1
        call trlint(x,y,datmod,ndtmod,bccoef,x2,y2,datmod2,bccoef2,
     *    frc_int1, frc_int2, nnt, nvartr, nrdtmr, nbccfr, iy, ivarmx)
c
      end if
c
c  Test for writing datmod
c
      if(idatmd.eq.1.and.istdpr.gt.0) then
	if(mod(nrdtmr,5).eq.0) then
	  krec = nrdtmr/5
        else
	  krec = nrdtmr/5+1
	end if
	write(istdpr,135) 
	do 60 k=1,krec
	i1=1+5*(k-1)
	i2=min0(i1+4,nrdtmr)
        write(istdpr,140) i1,(datmod(i),i=i1,i2)
   60   continue
c
	if(mod(nidtmr,5).eq.0) then
	  krec = nidtmr/5
        else
	  krec = nidtmr/5+1
	end if
	write(istdpr,145) 
	do 65 k=1,krec
	i1=1+5*(k-1)
	i2=min0(i1+4,nidtmr)
        write(istdpr,150) i1,(ndtmod(i),i=i1,i2)
   65   continue
c
      end if
c
c  test for rescaling central boundary condition coefficients and
c  resetting solution.
c
      call sclsol(y,nnt,iy,bccoef,datmod,ndtmod)
c
c  possibly reset storage of variables in common/ksider/
c
      call rsksdr(bccoef,bccfst,nspcmr,nspcmx)
      call store(bccfst,bccoef,nbcprv)
c
      icry = 0
c
c  test for storing angular velocity in omgrot if ipartr = 1
c
      if(ndtmod(35).gt.0.and.ipartr.eq.1) then
	isprot=mod(ndtmod(35),1000)
	isprtt=ndtmod(35)
	if(ndtmod(35).gt.1000) then
	  do n=1,nnt
	    omgrot(n)=y(nvartr,n)
          end do
	  if(istdpr.gt.0) then
	    write(istdpr,'(//
     *      '' Angular velocity set from input y, with isprot ='',i3/)')
     *      isprot
            write(istdpr,'(//'' n, omega:''/(i5,1pe13.5))')
     *        (n,omgrot(n),n=1,nnt,20)
          end if
        end if
      end if
c
      nrdtr1=nrdtrl
      if(istdpr.gt.0) write(istdpr,'(/'' Exit settrl''/)')
      return
  103 format(////10(1h*),' model no',i10,' is not on trial d/s.',
     *  ' Last model read had number',i10)
  135 format(/' Data for trial model.'/' datmod:')
  140 format(i3,':',1p5e13.5)
  145 format(/' ndtmod:')
  150 format(i3,':',5i13)
      end
      subroutine trlint(x,y,datmod,ndtmod,bccoef,x2,y2,datmod2,bccoef2,
     *  frc_int1, frc_int2, nn, nvartr, nrdtmr, nbccfr, iy, iy2)
      implicit double precision(a-h, o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.cz.d.incl'
      parameter(kdisc_max=10)
      logical set_disc
c
      dimension x(*), y(iy,*), datmod(*), ndtmod(*),bccoef(*), 
     *  x2(*), y2(iy2,*), datmod2(*), bccoef2(*)
      dimension ndisc1(kdisc_max),ndisc2(kdisc_max),
     *  xnew(nnmax),ynew(ivarmx,nnmax),
     *  pl(nnmax),pl2(nnmax),plnew(nnmax),yint(ivarmx)
c
c
c  Interpolate between models at two adjacent timesteps.
c  
c  Test for presence of discontinuity in composition (from
c  growing convective core) and take appropriate action.
c
c  NOTE: Includes valid test with treatment of variables involved in
c  diffusion and settling. However, it is striking that a negative
c  dX/dq is encountered.
c
c  Original version: 15/7/05
c
c  Modified 9/4/06, for use with diffusion and settling.
c
      common/eqstd/ xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp,gmm1(4)
      common/heavy/ zatmos, zhc, zh(1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      data idiag /1/
c
      if(istdpr.le.0) idiag=0
      if(idiag.ge.1) write(istdpr,'(/'' Entering trlint''/)')
c
      if(idiag.ge.2) then
        write(istdpr,'(a,2f15.11)') 
     *    'Enter trlint with frc_int1, frc_int2 =', frc_int1, frc_int2
        write(istdpr,
     *    '(/'' Model 1. n, x, y2, y3, y5:''/(i5,1p4e18.10))')
     *      (n,x(n),y(2,n),y(3,n),y(5,n),n=1,nn,20)
        write(istdpr,
     *    '(/'' Model 2. n, x, y2, y3, y5:''/(i5,1p4e18.10))')
     *    (n,x2(n),y2(2,n),y2(3,n),y2(5,n),n=1,nn,20)
      end if
c
c  test for composition discontinuities
c
      idisc1=0
      do n=nn-5,5,-1
        dxh=abs(y(5,n+1)-y(5,n))
        if(dxh.ge.1.e-6.and.dxh/(x(n)-x(n+1)).gt.
     *    5*abs(y(5,n+4)-y(5,n-4))/(x(n-4)-x(n+4)).and.
     *    idisc1.lt.kdisc_max) then
          if(idiag.gt.0) write(istdpr,'(//
     *      '' Discontinuity found on Model 1 at x ='',1pe13.5/
     *      '' n, x, X:''//(i5,1pe13.5,0pf13.8))') 
     *      x(n),(i,x(i),y(5,i),i=n-4,n+4)
          idisc1=idisc1+1
          ndisc1(idisc1)=n
        end if
      end do
c
      idisc2=0
      do n=nn-5,5,-1
        dxh=abs(y2(5,n+1)-y2(5,n))
        if(dxh.ge.1.e-6.and.dxh/(x2(n)-x2(n+1)).gt.
     *    5*abs(y2(5,n+4)-y2(5,n-4))/(x2(n-4)-x2(n+4)).and.
     *    idisc2.lt.kdisc_max) then
          if(idiag.gt.0) write(istdpr,'(//
     *      '' Discontinuity found on Model 2 at x ='',1pe13.5/
     *      '' n, x, X:''//(i5,1pe13.5,0pf10.5))') 
     *      x2(n),(i,x2(i),y2(5,i),i=n-4,n+4)
          idisc2=idisc2+1
          ndisc2(idisc2)=n
        end if
      end do
c
c  test for presence for and management of discontinuities
c  So far only deal with innermost discontinuity.
c  Certainly needs to be fixed.
c
      set_disc=.false.
      if(idisc1.ne.idisc2) then
        write(istdou,'(/''  ***** problem in trlint. idisc1 ='',i2,
     *    '' not eq idisc2 ='',i2/
     *                 ''        Ignore discontinuities''/)')
     *    idisc1, idisc2
        if(istdpr.gt.0.and.istdpr.ne.istdou)
     *    write(istdpr,'(/''  ***** problem in trlint. idisc1 ='',i2,
     *    '' not eq idisc2 ='',i2/
     *                 ''        Ignore discontinuities''/)')
     *    idisc1, idisc2
      else if(idisc1.gt.0) then
        nd1=ndisc1(1)
        nd2=ndisc2(1)
        xd1=x(nd1)
        xd2=x2(nd2)
        if(abs(xd1-xd2).gt.0.1) then 
          write(istdou,
     *      '(/''  ***** problem in trlint. xdisc1 ='',1pe13.5,
     *      '' too far from xdisc2 ='',e13.5/
     *                 ''        Ignore discontinuities''/)')
     *      xd1, xd2
        if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,
     *      '(/''  ***** problem in trlint. xdisc1 ='',1pe13.5,
     *      '' too far from xdisc2 ='',e13.5/
     *                 ''        Ignore discontinuities''/)')
     *      xd1, xd2
        else
          set_disc=.true.
        end if
      end if
c
c  Initially set up interpolated model on one mesh
c  use mesh in model closest to target
c  Interpolation near a possible discontinuity could be refined
c
      if(frc_int1.ge.0.5) then
        do n=1,nn
          xnew(n)=x(n)
          call lir(xnew(n),x2,yint,y2,nvartr,ivarmx,nn,n,inter)
          do i=1,nvartr
            ynew(i,n)=frc_int1*y(i,n)+frc_int2*yint(i)
          end do
        end do
c
      else
        do n=1,nn
          xnew(n)=x2(n)
          call lir(xnew(n),x,yint,y,nvartr,iy,nn,n,inter)
          do i=1,nvartr
            ynew(i,n)=frc_int1*yint(i)+frc_int2*y2(i,n)
          end do
        end do
      end if
c
c  test for resetting near discontinuities
c
      if(set_disc) then
c
c  set pressure in models for later resetting of log(f) (or whatever)
c  Could be more efficient by restricting to region near discontinuity.
c  ALSO badly needs fixing of usage of zh here!
c  For the time being use fixed zh, from initial value
c
        zhh=datmod(1)
c
        do n=1,nn
          fl=y(2,n)
          tl=y(3,n)
          xh=y(5,n)
          yh=1-xh-zhh
          call eqstf(fl,tl,xh,yh,zhh,nosd,notd)
          pl(n)=log10(pt(1))
        end do
c
        do n=1,nn
          fl=y2(2,n)
          tl=y2(3,n)
          xh=y2(5,n)
          yh=1-xh-zhh
          call eqstf(fl,tl,xh,yh,zhh,nosd,notd)
          pl2(n)=log10(pt(1))
        end do
c
c  set interpolated pressure on xnew
c
        do n=1,nn
          call lir(xnew(n),x,yint,pl,1,1,nn,n,inter)
          plnew(n)=frc_int1*yint(1)
        end do
        do n=1,nn
          call lir(xnew(n),x2,yint(1),pl2,1,1,nn,n,inter)
          plnew(n)=plnew(n)+frc_int2*yint(1)
        end do
c
c  set location of and range around discontinuity on chosen mesh
c
        xdi=frc_int1*xd1+frc_int2*xd2
        if(idiag.gt.0) write(istdpr,'(//
     *    '' Interpolated discontinuity at x ='',1pe13.5)')  xdi
        xdm=min(xd1,xd2)
        xdp=max(xd1,xd2)
        do n=1,nn
          if(xnew(n).gt.xdp) ndp=n
          if(xnew(n).gt.xdi) ndi=n
          if(xnew(n).gt.xdm) ndm=n
        end do
        ndp=ndp-3
        ndm=ndm+3
        ndi=ndi+1
c
c  reset composition around discontinuity (for now just hydrogen;
c  obviously has to be fixed)
c
c  Also, assumes that discontinuity is further out in model 2 than in
c  model 1, probably
c
        if(idiag.gt.0) 
     *    write(istdpr,'(/'' Reset X. n, x, old, new X:'')')
        dxhm=ynew(5,ndm)-y2(5,ndm)
        ni=1
        do n=ndp,ndi-1
          xhold=ynew(5,n)
          call lir(xnew(n),x,yint(1),y(5,1),1,iy,nn,ni,inter)
          ni=ni+1
          if(n.eq.ndp) dxhp=ynew(5,ndp)-yint(1)
          ynew(5,n)=yint(1)+dxhp
          if(idiag.gt.0) write(istdpr,'(i5,1pe13.5,0p2f10.5)')
     *      n,xnew(n),xhold,ynew(5,n)
        end do
        if(idiag.gt.0) write(istdpr,'('' '')')
        ni=0
        do n=ndm,ndi,-1
          xhold=ynew(5,n)
          call lir(xnew(n),x2,yint(1),y2(5,1),1,iy2,nn,ni,inter)
          ni=ni+1
          if(n.eq.ndm) dxhm=ynew(5,ndm)-yint(1)
          ynew(5,n)=yint(1)+dxhm
          if(idiag.gt.0) write(istdpr,'(i5,1pe13.5,0p2f10.5)')
     *      n,xnew(n),xhold,ynew(5,n)
        end do
c
c  Finally reset log(f) or log(rho) such that to recover log(p) with
c  new composition
c
        if(idiag.gt.0) write(istdpr,
     *    '(//'' Iterate for log f. n, x, old, new log f:'')')
        do n=ndp,ndm
c
c  iterate to determine new log f
c
          tl=ynew(3,n)
          xh=ynew(5,n)
          yh=1-xh-zhh
          fl=ynew(2,n)
          nit=0
c
   10     call eqstf(fl,tl,xh,yh,zhh,nosd,notd)
          dfl=(plnew(n)-log10(pt(1)))/pt(2)
          fl=fl+dfl
          nit=nit+1
          if(nit.ge.10) then
            write(istdou,*) 'Iteration in trlint unconverged. dfl =',
     *        dfl
            if(idiag.gt.0.and.istdpr.ne.istdou)
     *        write(istdou,*) 'Iteration in trlint unconverged. dfl =',
     *          dfl
            fl=y(2,n)
          else if(abs(dfl).ge.1.d-8) then
            go to 10
          end if
c
          if(idiag.gt.0) write(istdpr,'(i5,1pe13.5,0p2f10.5)')
     *      n, x(n), ynew(2,n),fl
          ynew(2,n)=fl
        end do
      end if
c
c  restore variables
c
      do n=1,nn
        x(n)=xnew(n)
        do i=1,nvartr
          y(i,n)=ynew(i,n)
        end do
      end do
c
      if(idiag.gt.0)
     *  write(istdpr,'(//'' Interpolated model: n, x, y(1-6)''/
     *  (i3,1p7e11.3))') (n, x(n), (y(i,n),i=1,6),n=1,nn)
c
      do i=1,nrdtmr
        datmod(i)=frc_int1*datmod(i)+frc_int2*datmod2(i)
      end do
      if(idiag.gt.0) 
     *  write(istdpr,'(//'' Interpolated datmod:''/
     *  (1p5e13.5))') (datmod(i),i=1,nrdtmr)
c
      do i=1,nbccfr
        bccoef(i)=frc_int1*bccoef(i)+frc_int2*bccoef2(i)
      end do
      if(idiag.gt.0) 
     *  write(istdpr,'(//'' Interpolated bccoef:''/
     *  (1p5e13.5))') (bccoef(i),i=1,nbccfr)
      return
      end
