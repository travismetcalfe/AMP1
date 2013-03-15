pro write_nemdl,lun,nn,ndata,datmod,idatmd,bc,x,y,file=kfile
;  writes unformatted evolution model to unit lun.
;  This is based on new format for evolution model output.
;  if file is set, opens file and writes first model
;  if the unit number is set in lun, this is used,
;  otherwise lun is set to 1.
;  It is assumed that ndata contains [nrdtmd,nidtmd,nvar,nbccf,iform]

;  Original version: 25/3/96

if keyword_set(kfile) then begin
	if keyword_set(lun) eq 0 then lun=1
	openw,lun,kfile,/f77_unformatted
endif

nrdtmd=ndata(0)
nidtmd=ndata(1)
nvar=ndata(2)
nbccf=ndata(3)
iform=ndata(4)

xy=dblarr(nvar+1,nn)
xy(0,*)=x
xy(1:nvar,*)=y

writeu,lun,iform,nn,nrdtmd,nidtmd,nvar,nbccf,datmod,idatmd,xy,bc

return
end
