      subroutine rdemdl(ids,x,y,datmod,ndtmod,bccoef,iy,iform,
     *  nn,nrdtmd,nidtmd,nbccf,nvar,icase,icry)
c
c  Reads evolution model file on unit ids, possibly testing for old
c  or new format. 
c
c  If icase = 0, test on case, after rewinding file.
c  If icase = 1, assume old format
c  If icase = 2, assume new format
c  If icase = -1, assume same format as in previous call
c    (if no format has been set, rewind and test)
c
c  icry is returned as 
c       1 for completion with rewind but without error
c       0 for completion without rewind and error
c      -1 for end of file 
c      -2 for error in read.
c
c  Original version: 13/8/91
c
      implicit double precision(a-h, o-z)
      dimension x(1), y(iy,1), datmod(1), ndtmod(1), bccoef(1)
      dimension datmrd(31),ndtmrd(3)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      equivalence(datmrd(29), ndtmrd(1))
      data icasel /-1/
c
      save
c
c  initialize icry to 0 (optimistically)
c
      icry = 0
c
c  test for testing for case
c
      if(icase.eq.0.or.(icase.eq.-1.and.icasel.eq.-1)) then
	icry = 1
	rewind ids
	read(ids,end=80,err=90) iread
	rewind ids
	if(iread.le.0) then
	  icasel = 2
        else
	  icasel = 1
        end if
      else if(icase.eq.-1) then
c
      else if(icase.eq.1.or.icase.eq.2) then
	icasel = icase
      else
	if(istdpr.gt.0) write(istdpr,110) icase
	icry = -2
	return
      end if
c
c  now read next model, depending on icasel
c
      if(icasel.eq.1) then
c
c  old format
c
	read(ids,end=80,err=90) nn,(datmrd(i),i=1,31),
     *    (x(n),(y(i,n),i=1,6),n=1,nn),(bccoef(i),i=1,54) 
        do 20 i=1,28
   20   datmod(i)=datmrd(i)
        do 25 i=1,3
   25   ndtmod(i)=ndtmrd(i)
c
c  set storage parameters for this case
c
	iform = 0
	nrdtmd = 28
	nidtmd = 3
	nbccf  = 54
	nvar   = 6
c
      else
c
c  new case
c
        read(ids,end=80,err=90) iform,nn,nrdtmd,nidtmd,nvar,nbccf,
     *    (datmod(i),i=1,nrdtmd),(ndtmod(i),i=1,nidtmd),
     *    (x(n),(y(i,n),i=1,nvar),n=1,nn),
     *    (bccoef(i),i=1,nbccf)
c
      end if
      return
c
c  diagnostics and return for end of file or error
c
   80 icry = -1
      if(istdpr.gt.0) write(istdpr,120) ids
      return
c
   90 icry = -2
      if(istdpr.gt.0) write(istdpr,130) ids
      return
  110 format(//' ***** Error in s/r rdemdl. icase = ',i4,
     *  ' not allowed')
  120 format(//' ***** EOF encountered in s/r rdemdl on unit',i4)
  130 format(//' ***** Error encountered in s/r rdemdl on unit',i4)
      end
