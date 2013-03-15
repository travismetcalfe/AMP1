pro read_nemdls,lun,nn,ndata,datmod,idatmd,bc,x,y,file=kfile,number=knum, $
    nrdmax=nrdmax,swap=swap
;  reads several unformatted evolution models from unit lun.
;  This is based on new format for evolution model output.
;  if file is set, opens file and reads first model
;  if the unit number is set in lun, this is used,
;  otherwise lun is set to 1.
;  number defines the number of the first model read
;  nrdmax may be set to the maximum number of models. Otherwise
;  read is to EOF or to an internally set maximum (numemdlmax)
;  If swap is set, use swap_endian, for reading, e.g., linux files
;  from bigcat

;  Original version: 19/12/99

common numemdl, nread

numemdlmax=200

if keyword_set(nrdmax) then mrdmax = nrdmax else mrdmax = numemdlmax

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
  datmodr=dblarr(nrdtmd)
  bcr=dblarr(nbccf)
  xr=dblarr(nn)
  yr=dblarr(nvar,nn)
  xy=dblarr(nvar+1,nn)
  readu,lun,nn,datmodr,xy,bcr
  xr=xy(0,*)
  yr=xy(1:nvar,*)
endif else begin
  readu,lun,iform,nn,nrdtmd,nidtmd,nvar,nbccf
  datmodr=dblarr(nrdtmd)
  idatmdr=lonarr(nidtmd)
  bcr=dblarr(nbccf)
  xr=dblarr(nn)
  yr=dblarr(nvar,nn)
  xy=dblarr(nvar+1,nn)
  point_lun,lun,pos
  readu,lun,iform,nn,nrdtmd,nidtmd,nvar,nbccf,datmodr,idatmdr,xy,bcr
  xr=transpose(xy(0,*))
  yr=xy(1:nvar,*)
endelse

if keyword_set(knum) then begin
  while nread lt knum do begin
	if iform gt 0 then begin
  	  readu,lun,nn,datmod,xy,bcr
	endif else begin
  	  readu,lun,iform,nn,nrdtmd,nidtmd,nvar,nbccf,datmodr,idatmdr,xy,bcr
	endelse
	nread=nread+1
  endwhile
  xr=transpose(xy(0,*))
  yr=xy(1:nvar,*)
endif

if iform gt 0 then $
  print,'nread, nn, nvar =',nread, nn $
else $
  print,'nread, nn, nrdtmd, nidtmd, nvar =',nread, nn, nrdtmd, nidtmd, nvar

ndata=[nrdtmd,nidtmd,nvar,nbccf,iform]

x=dblarr(nn,mrdmax+1)
y=dblarr(nvar,nn,mrdmax+1)

x(*,0)=xr
y(*,*,0)=yr
datmod=dblarr(nrdtmd,mrdmax+1)
bc=dblarr(nbccf,mrdmax+1)
datmod(*,0)=datmodr
bc(*,0)=bcr

if iform le 0 then begin
  idatmd=lonarr(nidtmd,mrdmax+1)
  idatmd(*,0)=idatmdr
endif

ird = 0

while (NOT EOF(lun) and ird lt mrdmax) do begin
  ird=ird+1
  if iform gt 0 then begin
    readu,lun,nn,datmod,xy,bcr
  endif else begin
    readu,lun,iform,nn,nrdtmd,nidtmd,nvar,nbccf,datmodr,idatmdr,xy,bcr
    idatmd(*,ird)=idatmdr
  endelse
  x(*,ird)=transpose(xy(0,*))
  y(*,*,ird)=xy(1:nvar,*)
  datmod(*,ird)=datmodr
  bc(*,ird)=bcr
endwhile

x=x(*,0:ird)
y=y(*,*,0:ird)
datmod=datmod(*,0:ird)
bc=bc(*,0:ird)
if iform le 0 then idatmd=idatmd(*,0:ird)

return
end
