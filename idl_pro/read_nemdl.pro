pro read_nemdl,lun,nn,ndata,datmod,idatmd,bc,x,y,file=kfile,number=knum, $
    swap=swap
;  reads unformatted evolution model from unit lun.
;  This is based on new format for evolution model output.
;  if file is set, opens file and reads first model
;  if the unit number is set in lun, this is used,
;  otherwise lun is set to 1.
;  If swap is set, use swap_endian, for reading, e.g., linux files
;  from bigcat

;  Original version: 8/5/92

common numemdl, nread

if keyword_set(kfile) then begin
	if keyword_set(lun) eq 0 then lun=1
	if keyword_set(swap) then swap = 1 else swap = 0
	openr,lun,kfile,/f77_unformatted,swap_endian=swap
	nread=0
endif

if EOF(lun) then begin
	print,' ****** EOF on unit',lun,''
	nn=-1
	return
endif

point_lun,-lun,pos

iform = 1l
nn=1l
nvar=6l
nrdtmd=31l
nidtmd=0l
nbccf=54l

;  test for format

readu,lun,iform
point_lun,lun,pos
nread=nread+1

if iform gt 0 then begin
  nn=iform
  datmod=dblarr(nrdtmd)
  bc=dblarr(nbccf)
  x=dblarr(nn)
  y=dblarr(nvar,nn)
  xy=dblarr(nvar+1,nn)
  readu,lun,nn,datmod,xy,bc
  x=xy(0,*)
  y=xy(1:nvar,*)
endif else begin
  readu,lun,iform,nn,nrdtmd,nidtmd,nvar,nbccf
  datmod=dblarr(nrdtmd)
  idatmd=lonarr(nidtmd)
  bc=dblarr(nbccf)
  x=dblarr(nn)
  y=dblarr(nvar,nn)
  xy=dblarr(nvar+1,nn)
  point_lun,lun,pos
  readu,lun,iform,nn,nrdtmd,nidtmd,nvar,nbccf,datmod,idatmd,xy,bc
  x=transpose(xy(0,*))
  y=xy(1:nvar,*)
endelse

if keyword_set(knum) then begin
  while nread lt knum do begin
	if iform gt 0 then begin
  	  readu,lun,nn,datmod,xy,bc
	endif else begin
  	  readu,lun,iform,nn,nrdtmd,nidtmd,nvar,nbccf,datmod,idatmd,xy,bc
	endelse
	nread=nread+1
  endwhile
  x=transpose(xy(0,*))
  y=xy(1:nvar,*)
endif

if iform gt 0 then $
  print,'nread, nn, nvar =',nread, nn $
else $
  print,'nread, nn, nrdtmd, nidtmd, nvar =',nread, nn, nrdtmd, nidtmd, nvar

ndata=[nrdtmd,nidtmd,nvar,nbccf,iform]
return
end
