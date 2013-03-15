      subroutine settrl(idstrl,itrial,nmdtrl,nrdtr1,nvartr,idatmd,
     *  nspcmx,nbcprv,x,y,datmod,ndtmod,bccoef,iy,iformr,nnt,
     *  nrdtmr,nidtmr,nbccfr,nspcmr,icnotr,ihectr,icry)
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
      implicit double precision(a-h, o-z)
      real xrds, yrds, datrds,bccrds
      dimension x(*), y(iy,*), datmod(*), ndtmod(*), bccoef(*)
      dimension datrds(31), ndtrds(3), bccrds(54), datmdr(31),
     *  ndtmdr(3),bccfst(1)
      common/sooner/ xrds(1)
      common/work/ yrds(6,1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      equivalence(datmdr(29), ndtmdr(1)),(datrds(29), ndtrds(1)),
     *  (bccfst(1),xrds(1))
c
      data nrdtrl /-1/
c
c  initialize icasrd so that first read tests for type of model
c
      data icasrd / 0/
c
      save icasrd
c
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'Enter settrl with itrial, nmdtrl, nrdtr1 =',
     *  itrial, nmdtrl, nrdtr1
c
c  if nrdtr1 on input is negative, reset nrdtrl
c
      if(nrdtr1.lt.0) nrdtrl = -1
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
   40 if(nmdtrl.lt.1000) then
c
c  error for end of file on trial data set
c
        write(istdou,103) nmdtrl,nrdtrl
        if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,103) 
     *    nmdtrl,nrdtrl
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
c  Test for writing datmod
c
      if(idatmd.eq.1) then
	if(mod(nrdtmr,5).eq.0) then
	  krec = nrdtmr/5
        else
	  krec = nrdtmr/5+1
	end if
	if(istdpr.gt.0) write(istdpr,135) 
	do 60 k=1,krec
	i1=1+5*(k-1)
	i2=min0(i1+4,nrdtmr)
        if(istdpr.gt.0) write(istdpr,140) i1,(datmod(i),i=i1,i2)
   60   continue
c
	if(mod(nidtmr,5).eq.0) then
	  krec = nidtmr/5
        else
	  krec = nidtmr/5+1
	end if
	if(istdpr.gt.0) write(istdpr,145) 
	do 65 k=1,krec
	i1=1+5*(k-1)
	i2=min0(i1+4,nidtmr)
        if(istdpr.gt.0) write(istdpr,150) i1,(ndtmod(i),i=i1,i2)
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
      nrdtr1=nrdtrl
      return
  103 format(////10(1h*),' model no',i10,' is not on trial d/s.',
     *  ' Last model read had number',i10)
  135 format(/' Data for trial model.'/' datmod:')
  140 format(i3,':',1p5e13.5)
  145 format(/' ndtmod:')
  150 format(i3,':',5i13)
      end
