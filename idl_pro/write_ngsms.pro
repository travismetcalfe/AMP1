pro write_ngsms,gsms,igsms,file=file
;
;  writes grand summaries defined in gsms and igsms to file file

openw,99,file,/f77_unformatted

n_gsms=n_elements(gsms(0,*))
n_igsms=n_elements(igsms(0,*))

igsms_full=lonarr(24)
if n_gsms ne n_igsms then begin
  print,' *** Warning: n_gsms =',n_gsms,'  n_igsms =',n_igsms
  n_gsms=min([n_gsms,n_igsms])
endif

for i=0,n_gsms-1 do begin
  igsms_full(0:7)=igsms(0:7,i)
  writeu,99,gsms(*,i),igsms_full
endfor

close,99

return
end
