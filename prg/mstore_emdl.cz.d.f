      subroutine mstore_emdl(ntime, nn, nrdtmd, nidtmd, nvar, nbccf,
     *  datmod, ndtmod, x, y, bccoef, iy)
      implicit double precision (a-h, o-z)
      dimension datmod(*), ndtmod(*), x(*), y(iy,*), bccoef(*)
c
c  Stores evolution models internally, for later use in setting
c  up selected oscillation results or gong files.
c
c  Note that although the storage arrays are treated internally
c  as 1-dimensional, they have the following equivalent meaning:
c
c  store_emdl(i,n,k)
c  datmod(j,k)
c  ndtmod(j,k)
c  bccoef(j,k)
c
c  Here k numbers the timestep (with k = 1 for the initial timestep,
c  typically the ZAMS model; hence k = ntime+1, where ntime is
c  the timestep number passed as input to the routine)
c
c  n numbers the meshpoints in the model
c  i = 1 corresponds to x (store_datmod(1,n,.) = x(n))
c  i > 1 corresponds to y(i-1,.) (store_datmod(i,n,.) = y(i-1,n))
c  j numbers the variables in datmod, ndtmod or bccoef
c
c  The (implicit) dimensioning of the storage arrays is defined as
c
c  store_emdl(nvar+1,nn,*)
c  datmod(nrdtmd,*)
c  ndtmod(nidtmd,*)
c  bccoef(nbccf,*)
c
c  where nvar, nn, nrdtmd, nidtmd and nbccf are passed as arguments.
c  Thus these variables must not change during an evolution run.
c
c  If entered with ntime .lt. 0, store model in istore_emdl+1, etc.,
c  i.e., as the subsequent model.
c
c  istore_emdl, istore_datmod, istore_ndtmod and istore_bc
c  are set to the number of the latest model stored in a given
c  sequence. (Thus they should all be identical; they are kept
c  as separate variables for possible future flexibility.)
c
c  Original version: 20/7/06
c
c  commons for internal storage of model variables
c
      common/cstore_emdl/ nstore_emdl, istore_emdl, store_emdl(1)
      common/cstore_datmod/ nstore_datmod, istore_datmod, 
     *  store_datmod(1)
      common/cstore_ndtmod/ nstore_ndtmod, istore_ndtmod, 
     *  lstore_ndtmod(1)
      common/cstore_bc/ nstore_bc, istore_bc, store_bc(1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  set location for next model
c
      if(ntime.lt.0) then
	kstore_emdl  =istore_emdl+1
	kstore_datmod=istore_datmod+1
	kstore_ndtmod=istore_ndtmod+1
	kstore_bc    =istore_bc+1
      else
	kstore_emdl  =ntime+1
	kstore_datmod=ntime+1
	kstore_ndtmod=ntime+1
	kstore_bc    =ntime+1
      end if
c
c  test for adequate storage space
c
c  Note: stop to be replaced by output of model to file,
c  change to file output
c
      itot_emdl=  (nvar+1)*nn*kstore_emdl
      itot_datmod=nrdtmd*kstore_datmod
      itot_ndtmod=nidtmd*kstore_ndtmod
      itot_bc=    nbccf*kstore_bc
c
      if(itot_emdl.gt.nstore_emdl) then
	write(istdou,110) kstore_emdl, itot_emdl, nstore_emdl, 
     *    nvar+1, nn
	if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,110) 
     *    kstore_emdl, itot_emdl, nstore_emdl, nvar+1, nn
	stop 'mstore_emdl'
      else if(itot_datmod.gt.nstore_datmod) then
	write(istdou,120) kstore_datmod, itot_datmod, nstore_datmod, 
     *    nrdtmd
	if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,120) 
     *    kstore_datmod, itot_datmod, nstore_datmod, nrdtmd
	stop 'mstore_emdl'
      else if(itot_ndtmod.gt.nstore_ndtmod) then
	write(istdou,130) kstore_ndtmod, itot_ndtmod, nstore_ndtmod, 
     *    nidtmd
	if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,130) 
     *    kstore_ndtmod, itot_ndtmod, nstore_ndtmod, nidtmd
	stop 'mstore_emdl'
      else if(itot_bc.gt.nstore_bc) then
	write(istdou,140) kstore_bc, itot_bc, nstore_bc, nnbccf
	if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,140) 
     *    kstore_bc, itot_bc, nstore_bc, nidtmd
	stop 'mstore_emdl'
      end if
c
c  store emdl
c
      j=(kstore_emdl-1)*nn*(nvar+1)
      do n=1,nn
	j=j+1
	store_emdl(j)=x(n)
	do i=1,nvar
	 j=j+1
	 store_emdl(j)=y(i,n)
        end do
      end do
      if(istdpr.gt.0) write(istdpr,
     *  '(/'' Store emdl with nn, nvar ='',i5,i3/
     *     '' Number of variables used in store_emdl:'',i10)')
     *  nn, nvar, j
c
c  store datmod
c
      j0=(kstore_datmod-1)*nrdtmd
      do j=1,nrdtmd
	store_datmod(j0+j)=datmod(j)
      end do
      if(istdpr.gt.0) write(istdpr,
     *  '(/'' Store datmod with nrdtmd ='',i5/
     *     '' Number of variables used in store_datmod:'',i10)')
     *  nrdtmd, j0+nrdtmd
c
c  store ndtmod
c
      j0=(kstore_ndtmod-1)*nidtmd
      do j=1,nidtmd
	lstore_ndtmod(j0+j)=ndtmod(j)
      end do
      if(istdpr.gt.0) write(istdpr,
     *  '(/'' Store ndtmod with nidtmd ='',i5/
     *     '' Number of variables used in store_ndtmod:'',i10)')
     *  nidtmd, j0+nidtmd
c
c  store bccoef
c
      j0=(kstore_bc-1)*nbccf
      do j=1,nbccf
	store_bc(j0+j)=bccoef(j)
      end do
      if(istdpr.gt.0) write(istdpr,
     *  '(/'' Store bccoef with nbccf ='',i5/
     *     '' Number of variables used in store_bc:'',i10)')
     *  nbccf, j0+nbccf
c
      istore_emdl  =kstore_emdl
      istore_datmod=kstore_datmod
      istore_ndtmod=kstore_ndtmod
      istore_bc    =kstore_bc
      return
  110 format(/' ***** Error in s/r mstore_emdl for model no',i5/
     *        '       For emdl, required space =',i10,
     *        ' exceeds available space =', i10/
     *        '       nvar + 1 =',i3,' nn =',i5)
  120 format(/' ***** Error in s/r mstore_emdl for model no',i5/
     *        '       For datmod, required space =',i10,
     *        ' exceeds available space =', i10/
     *        '       ndtmod =',i4)
  130 format(/' ***** Error in s/r mstore_emdl for model no',i5/
     *        '       For ndtmod, required space =',i10,
     *        ' exceeds available space =', i10/
     *        '       nitmod =',i4)
  140 format(/' ***** Error in s/r mstore_emdl for model no',i5/
     *        '       For bccoef, required space =',i10,
     *        ' exceeds available space =', i10/
     *        '       nbccf =',i4)
      end
      subroutine mread_emdl(nmodrd,x,y,datmod,ndtmod,bccoef,iy,
     *  nn,nrdtmd,nidtmd,nbccf,nvar,icry)
c
c  Extracts evolution model no nmodrd from in internal storage
c  Note: model numbering starts at 1.
c
c  Storage is defined in comments to s/r istore_emdl
c
c  Original version: 20/7/06
c
c  commons for internal storage of model variables
c
      implicit double precision(a-h, o-z)
      dimension x(*),y(iy,*),datmod(*),ndtmod(*),bccoef(*)
      common/cstore_emdl/ nstore_emdl, istore_emdl, store_emdl(1)
      common/cstore_datmod/ nstore_datmod, istore_datmod, 
     *  store_datmod(1)
      common/cstore_ndtmod/ nstore_ndtmod, istore_ndtmod, 
     *  lstore_ndtmod(1)
      common/cstore_bc/ nstore_bc, istore_bc, store_bc(1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      icry=0
c
c  test for model in range of stored models
c
      if(nmodrd.gt.istore_emdl) then
	write(istdou,110) nmodrd, istore_emdl
	if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,110) 
     *    nmodrd, istore_emdl
	icry=-1
	return
      else if(nmodrd.gt.istore_datmod) then
	write(istdou,120) nmodrd, istore_datmod
	if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,120) 
     *    nmodrd, istore_datmod
	icry=-1
	return
      else if(nmodrd.gt.istore_ndtmod) then
	write(istdou,130) nmodrd, istore_ndtmod
	if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,130) 
     *    nmodrd, istore_ndtmod
	icry=-1
	return
      else if(nmodrd.gt.istore_bc) then
	write(istdou,140) nmodrd, istore_bc
	if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,140) 
     *    nmodrd, istore_bc
	icry=-1
	return
      end if
c
c  extract emdl
c
      j=(nmodrd-1)*nn*(nvar+1)
      do n=1,nn
	j=j+1
	x(n)=store_emdl(j)
	do i=1,nvar
	 j=j+1
	 y(i,n)=store_emdl(j)
        end do
      end do
c
c  extract datmod
c
      j0=(nmodrd-1)*nrdtmd
      do j=1,nrdtmd
	datmod(j)=store_datmod(j0+j)
      end do
c
c  extract ndtmod
c
      j0=(nmodrd-1)*nidtmd
      do j=1,nidtmd
	ndtmod(j)=lstore_ndtmod(j0+j)
      end do
c
c  extract bccoef
c
      j0=(nmodrd-1)*nbccf
      do j=1,nbccf
	bccoef(j)=store_bc(j0+j)
      end do
c
      return
  110 format(/' ***** Error in s/r mread_emdl.'/
     *        '       Requested model no =',i5,
     *        ' exceeds number of emdl records =',i5)
  120 format(/' ***** Error in s/r mread_emdl.'/
     *        '       Requested model no =',i5,
     *        ' exceeds number of datmod records =',i5)
  130 format(/' ***** Error in s/r mread_emdl.'/
     *        '       Requested model no =',i5,
     *        ' exceeds number of ndtmod records =',i5)
  140 format(/' ***** Error in s/r mread_emdl.'/
     *        '       Requested model no =',i5,
     *        ' exceeds number of bccoef records =',i5)
      end
