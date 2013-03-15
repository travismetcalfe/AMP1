pro write_amdes,lun,x,y,nst,cs,file=kfile,mcase=mcase

;  writes adiabatic eigenfunctions to file kfile
;  if mcase set set and equal to 2, assume reduced set,
;  otherwise full set is assumed
;  nst+1 is the total number of eigenfunctions to be written

if keyword_set(mcase) then icase = mcase else icase = 1

openw,lun,kfile,/f77_unformatted

nn=long(n_elements(x))
print,'nn =',nn

if icase eq 1 then begin
  iprec=size(x)
  iprec=iprec(2)
  if iprec eq 4 then yw=fltarr(7,nn) else yw=dblarr(7,nn)
endif else begin
  writeu,lun,nn,x
endelse

for i=0,nst do begin
  if icase eq 1 then begin
    yw(0,*)=x 
    yw(1:6,*)=y(*,*,i)
    writeu,lun,cs(*,i),nn,yw
  endif else begin
    yw=y(*,*,i)
    writeu,lun,cs(*,i),yw
  endelse
endfor
close,lun

return
end

