pro write_amdl,lun,nmodel,nn,data,aa,file=kfile
;  writes unformatted adiabatic oscillation model to unit lun,
;  if file is set, opens file and writes model
;  if the unit number is set in lun, this is used,
;  otherwise lun is set to 1.

common numamdl, nread

if keyword_set(kfile) then begin
	if keyword_set(lun) eq 0 then lun=1
	openw,lun,kfile,/f77_unformatted
	nread=0
endif

writeu,lun,nmodel,nn,data,aa


return
end
