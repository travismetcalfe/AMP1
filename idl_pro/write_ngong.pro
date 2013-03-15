pro write_ngong,lun,cdata,nmod,nn,ndata,datmod,idatmd,datgng,bc,yvar, $
   file=kfile, diag=diag
;  writes unformatted GONG model to unit lun.
;  if file is set, opens file 
;  if the unit number is set in lun, this is used,
;  otherwise lun is set to 1.

;  Original version: 1/7/03

common numgong, nwr

if keyword_set(kfile) then begin
	if keyword_set(lun) eq 0 then lun=1
	openw,lun,kfile,/f77_unformatted
	nwr=0
endif

if keyword_set(diag) then idiag = diag else idiag = 0

iform=-2l

writeu,lun,cdata,nmod,iform,nn,ndata, $
    datmod, idatmd, datgng, bc, yvar
nwr=nwr+1
  
; Note: ndata=[nrdtmd,nidtmd,ndtgng,nvar,nbccf]

return
end
