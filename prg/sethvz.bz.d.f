      subroutine sethvz(x,y,nn,iy,datmod,ndtmod,z,isetzh,idshvz)
c
c  Sets heavy-element abundance into common/heavy/, depending on
c  isetzh:
c  isetzh = 0: set uniform composition at z
c  isetzh .gt. 0: set composition based on values read from 
c  d/s idshvz and interpolated to mesh x = log10(m/M)
c  for isetzh = 2, scale z(h), such that surface value is z.
c
c  isetzh = -1: Set Z from trial model, provided in y, datmod and
c  ndtmod
c  isetzh = -2: With helium burning, set Z from X and Y in trial model,
c  otherwise set Z from datmod(1) in trial model
c
c  The data are (at least for isetzh = 1 and 2) assumed to be in 
c  ASCII form, given as log10(m/M), Z.
c
c  Original version: 19/7/95
c
c  Modified 11/6/99, to allow setting Z from trial model in cases
c  with heavy-element settling
c
c  Modified 8/10/02, allowing setting Z from X and Y in trial model,
c  in case of 4He burning
c
      implicit double precision(a-h, o-z)
      parameter(nnrmax=5000)
      dimension x(*),y(iy,*),datmod(*),ndtmod(*)
      dimension xr(nnrmax), zr(nnrmax)
      common/heavy/ zatmos, zhc, zh(1)
      common/cstyst/ iyh1, iyhe3, iyc12, iyc13, iyn14, iyo16, 
     *  iyche4, iycc12
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      data idshvp /-1/
c
      save
c
      if(isetzh.eq.0) then
	do 15 n=1,nn
   15   zh(n)=z
	zatmos = z
	zhc = z
      else if(isetzh.eq.1.or.isetzh.eq.2) then
	if(idshvz.ne.idshvp) then
	  call skpcom(idshvz)
	  n=1
   20     if(n.gt.nnrmax) then
	    write(istdou,110) n,nnrmax
	    if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,110) 
     *        n,nnrmax
	    stop 'Stop 1 in sethvz'
          end if
	  read(idshvz,*,end=25) xr(n),zr(n)
c..	  write(6,*) n, xr(n), zr(n)
	  n=n+1
	  go to 20
c
   25     nnr=n-1
	  idshvp=idshvz
        end if
c
c  carry out interpolation to model mesh
c
	do 30 n=1,nn
        call lir(x(n),xr,zh(n),zr,1,1,nnr,n,inter)
c..	write(6,*) 'new:',n,x(n),zh(n)
   30   continue
c
c  for isetzh = 2, scale to obtain surface value of z
c
	if(isetzh.eq.2) then
	  zscale=z/zh(1)
	  do 35 n=1,nn
   35     zh(n)=zscale*zh(n)
	end if
c
	zatmos=zh(1)
	zhc=zh(nn)
      else if(isetzh.eq.-1) then
c
c  set Z from trial model
c
	idiffus=ndtmod(21)

c  
	if(idiffus.le.1) then
	  zhh=datmod(1)
	  if(istdpr.gt.0) write(istdpr,120) idiffus, zhh
	  do 40 n=1,nn
   40     zh(n)=zhh
        else
c
	  idcomp=ndtmod(28)
	  iccomp=ndtmod(29)
	  iz=6+iccomp
	  do 45 n=1,nn
   45     zh(n)=y(iz,n)
c
	end if
c
	zatmos=zh(1)
	zhc=zh(nn)
      else if(isetzh.eq.-2) then
c
c  set Z from X and Y in trial model (if 4He burning is included)
c
	if(iyche4.le.0) then
	  zhh=datmod(1)
	  if(istdpr.gt.0) write(istdpr,130) iyche4, zhh
	  do 50 n=1,nn
   50     zh(n)=zhh
        else
c
	  do 55 n=1,nn
   55     zh(n)=1.d0-y(5,n)-y(iyche4,n)
c
	end if
c
	zatmos=zh(1)
	zhc=zh(nn)
c
      else
c  
c  diagnostics
c
	write(istdou,180) isetzh
	if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,110) isetzh
	stop 'Stop 2 in sethvz'
      end if
c
      return
  110 format(' ***** Error in reading Z(q) in s/r sethvz.'/
     *       '       n =',i5,' exceeds nmax =',i5)
  120 format(' ***** Warning: idiffus = ',i2,' for isetzh = -1'/
     *       '       Z set to',1pe13.5,' from datmod')
  130 format(' ***** Warning: iyche4 = ',i2,' for isetzh = -1'/
     *       '       Z set to',1pe13.5,' from datmod')
  180 format(' ***** Error in s/r sethvz. isetzh =',i4,
     *  ' not implemented')
      end
