pro read_ngong,lun,cdata,nmod,nn,ndata,datmod,idatmd,datgng,bc,yvar, $
   file=kfile, number=knum,diag=diag,new=new,quiet=quiet,swap=swap
;  reads unformatted GONG model from unit lun.
;  if file is set, opens file and reads first model
;  if the unit number is set in lun, this is used,
;  otherwise lun is set to 1.
;  If new is set, force read of new format, otherwise test 
;  on data
;  If swap is set, use swap_endian, for reading, e.g., linux files
;  from bigcat

;  Modified 5/5/92, to allow reading also new format of models

common numgong, nread

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

if keyword_set(diag) then idiag = diag else idiag = 0

point_lun,-lun,pos

cdata=make_array(4,/string,value=string(replicate(32b,80)))

iform = 1l
nmod = 1l
nn=1l
nvar=1l
nrdtmd=31l
nidtmd=0l
ndtgng=30l
nbccf=54l

;  test for format

if keyword_set(new) then inew = 1 else inew = 0

if inew eq 0 then begin
  readu,lun,cdata,nmod,iform
  if idiag ge 1 then print,cdata,nmod,iform
  point_lun,lun,pos
  if iform lt 0 then inew = 1
endif

if inew eq 0 then begin
  readu,lun,cdata,nmod,nn,nvar
  if idiag ge 1 then print,cdata,nmod,nn,nvar
  nread=nread+1
  datmod=dblarr(nrdtmd)
  datgng=dblarr(ndtgng)
  bc=dblarr(nbccf)
  yvar=dblarr(nvar,nn)
  point_lun,lun,pos
  readu,lun,cdata,nmod,nn,nvar,datmod,datgng,bc,yvar
endif else begin
  readu,lun,cdata,nmod,iform,nn,nrdtmd,nidtmd,ndtgng,nvar,nbccf
  if idiag ge 1 then print,cdata,nmod,iform,nn,nrdtmd,nidtmd,ndtgng,nvar,nbccf
  nread=nread+1
  datmod=dblarr(nrdtmd)
  idatmd=lonarr(nidtmd)
  datgng=dblarr(ndtgng)
  bc=dblarr(nbccf)
  yvar=dblarr(nvar,nn)
  point_lun,lun,pos
  readu,lun,cdata,nmod,iform,nn,nrdtmd,nidtmd,ndtgng,nvar,nbccf, $
    datmod, idatmd, datgng, bc, yvar
endelse

if keyword_set(knum) then begin
  while nread lt knum do begin
	if inew eq 0 then begin
  	  readu,lun,cdata,nmod,nn,nvar,datmod,datgng,bc,yvar
	endif else begin
  	  readu,lun,cdata,nmod,iform,nn,nrdtmd,nidtmd,ndtgng,nvar,nbccf, $
    	    datmod, idatmd, datgng, bc, yvar
	endelse
	nread=nread+1
  endwhile
endif

if not keyword_set(quiet) then begin
  if inew eq 0 then $
    print,'nread, nn, nvar =',nread, nn, nvar $
  else $
    print,'nread, nn, nrdtmd, nidtmd, ndtgng, nvar =', $
       nread, nn, nrdtmd, nidtmd, ndtgng, nvar
endif
  
ndata=[nrdtmd,nidtmd,ndtgng,nvar,nbccf]
return
end
