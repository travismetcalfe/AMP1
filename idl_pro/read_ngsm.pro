pro read_ngsm,lun,cs,ics,file=kfile,swap=swap

;  reads binary grand summary from file, including both the real part
;  (cs(0:37)) and the integer part (ics(0:7))
;  if file is set, opens file and reads first grand summary
;  if the unit number is et in lun, this is used,
;  otherwise lun is set to 1.

;  if file is not set, read is from unit lun, which is assumed
;  to be open

; On EOF, cs(0) and cs(1) are returned as -1.

;  data assumed to be in double precision


if keyword_set(kfile) then begin
	if keyword_set(lun) eq 0 then lun=1
	if keyword_set(swap) then swap = 1 else swap = 0
	openr,lun,kfile,/f77_unformatted,swap_endian=swap
endif

cs=dblarr(38) 
icsr=lonarr(24)

if EOF(lun) then begin
;	print,' ****** EOF on unit',lun,''
	cs(0:1)=-1
	return
endif

point_lun,-lun,pos

readu,lun,cs,icsr
ics=icsr(0:7)

return
end
