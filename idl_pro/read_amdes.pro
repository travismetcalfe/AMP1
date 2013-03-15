pro read_amdes,lun,x,y,nst,cs,file=kfile,double=kdouble,diag=kdiag, $
  mcase=mcase,nskip=nskip,nrdmax=nrdmax,modes=modes,resorder=resorder, $
  swap=swap

;  reads several binary eigenfunctions from file.
;  if file is set, opens file and reads first eigenfunction
;  if the unit number is set in lun, this is used,
;  otherwise lun is set to 1.
;  if mcase is set and is equal to 2, read reduced set
;  Otherwise full set is assumed
;  If nskip is set, return kernels at every nskip-th point
;  If keyword nrdmax is set, read at most nrdmax eigenfunctions.
;  If keyword modes is set, select only those modes that match modes
;  in modes. modes is an array whose first two columns must contain
;  degree and order. It is assumed to be ordered in the same
;  way as the eigenfunction file.

;  if file is not set, read is from unit lun, which is assumed
;  to be open

;  if double is set, the data are in double precision,
;  otherwise in single precision

;  If resorder is set, reset order for l = 1, by counting down from
;  top (so far very naively)

;  Note: it is assumed that x is the same for all eigenfunctions.
;  x(0:nn-1) returns the radial mesh (r/R)
;  nst+1 returns the number of modes
;  The eigenfunctions are returned in y(0:ivar-1,0:nn-1,0:nst)
;  cs(0:49,0:nst) return the grand summaries for the modes.

;  x
;  
common c_read_amde,init,nn,xx

numamdemax=1000

if keyword_set(mcase) then icase = mcase else icase = 1

if keyword_set(nrdmax) then mrdmax = nrdmax else mrdmax = numamdemax

if keyword_set(kfile) then begin
	if keyword_set(lun) eq 0 then lun=1
	if keyword_set(swap) then swap = 1 else swap = 0
	openr,lun,kfile,/f77_unformatted,swap_endian=swap
        init=1
endif else begin
	init=0
endelse

set_modes = keyword_set(modes)

if set_modes then imodes = n_elements(modes(0,*)) else imodes = 1

iker=0
im=0


if EOF(lun) then begin
	print,' ****** EOF on unit',lun,'^G'
	n=-1
	return
endif

n=1l
if keyword_set(kdouble) then begin
  csr=dblarr(50) 
  cs =dblarr(50,mrdmax+1)
endif else begin
  csr=fltarr(50)
  cs =fltarr(50,mrdmax+1)
endelse

point_lun,-lun,pos

if icase eq 2 then begin
	if init eq 1 then begin
	  readu,lun,n
	  point_lun,lun,pos
          if keyword_set(kdouble) then x=dblarr(n) else x=fltarr(n)
	  readu,lun,n,x
	  nn=n
	  xx=x
	endif
        if keyword_set(kdouble) then yr=dblarr(2,nn) else yr=fltarr(2,nn)
	iy=2
	readu,lun,csr,yr
	x=xx
endif else begin
	readu,lun,csr,n
	if keyword_set(kdouble) then a=dblarr(7,n) else a=fltarr(7,n)
	iy=6
	point_lun,lun,pos
	readu,lun,csr,n,a
	x=a(0,*)
	yr=a(1:6,*)
endelse

if keyword_set(nskip) then nsk=nskip else nsk = 1

if nsk gt 1 then begin
	nst=1+(n-1)/nsk
	sel=nsk*indgen(nst)
	x=x(sel)
endif else begin
	nst=n
	sel=indgen(n)
endelse

;  set storage y array and set first mode

if keyword_set(kdouble) then y=dblarr(iy,nst,mrdmax+1) else $
			     y=fltarr(iy,nst,mrdmax+1) 
im  = 0

if set_modes then begin
    l=csr(17)
    iord=csr(18)
    iorder=order_array([l,iord],modes(0:1,im))
    if iorder eq 1 then begin
      while iorder eq 1 and im lt imodes do begin
        print,modes(0:1,im), format = $
          '(" **** Error: mode l, order =",2f6.1," not found")'
        im=im+1
        if im lt imodes then iorder=order_array([l,iord],modes(0:1,im))
      endwhile
    endif
    if iorder eq 0 then begin
      cs(*,0)=csr
      y(*,*,0)=yr(*,sel)
      if keyword_set(kdiag) then $
        print,' In read_amde l, n, freq  =', csr(17),csr(18),csr(26)
      im=im+1
    endif
endif else begin
  cs(*,0)=csr
  y(*,*,0)=yr(*,sel)
  if keyword_set(kdiag) then $
    print,' In read_amde l, n, freq  =', csr(17),csr(18),csr(26)
endelse

imd = 0

while (NOT EOF(lun) and imd lt mrdmax and im lt imodes) do begin
  if icase eq 2 then begin
	readu,lun,csr,yr
  endif else begin
	readu,lun,csr,n,a
	yr=a(1:6,*)
  endelse
  if set_modes then begin
      l=csr(17)
      iord=csr(18)
      iorder=order_array([l,iord],modes(0:1,im))
      if iorder eq 1 then begin
        while iorder eq 1 and im lt imodes do begin
          print,modes(0:1,im), format = $
            '(" **** Error: mode l, order =",2f6.1," not found")'
          im=im+1
          if im lt imodes then iorder=order_array([l,iord],modes(0:1,im))
        endwhile
      endif
      if iorder eq 0 then begin
        imd=imd+1
        cs(*,imd)=csr
        y(*,*,imd)=yr(*,sel)
        if keyword_set(kdiag) then $
           print,' In read_amde l, n, freq  =', csr(17),csr(18),csr(26)
        im=im+1
      endif
  endif else begin
      imd=imd+1
      cs(*,imd)=csr
      y(*,*,imd)=yr(*,sel)
      if keyword_set(kdiag) then $
         print,' In read_amde l, n, freq  =', csr(17),csr(18),csr(26)
  endelse
endwhile

cs=cs(*,0:imd)
y=y(*,*,0:imd)

if keyword_set(resorder) then begin
  s1=where(cs(17,*) eq 1,ns1)
  cs1=cs(*,s1)
  for i=ns1-2,0,-1 do begin
    if cs1(18,i+1) ne 1 then $
      new_order=cs1(18,i+1)-1 else $
      new_order=-1
    if new_order ne cs1(18,i) then begin
      print,' **** change order from '+string(cs1(18,i),format='(i4)')+ $
                               '  to '+string(new_order,format='(i4)')
      cs1(18,i)=new_order
    endif
  endfor
  cs(*,s1)=cs1
endif

nst=imd

return
end
