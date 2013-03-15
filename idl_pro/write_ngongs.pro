pro write_ngongs,lun,cdatas,nmod,nn,ndata,datmods,idatmds,datgngs,bcs,yvars, $
   file=kfile
;  writes unformatted GONG models to unit lun.
;  if file is set, opens file.
;  if the unit number is set in lun, this is used,
;  otherwise lun is set to 1.

;  Original version: 1/7/03.

common numgong, nwr

if keyword_set(kfile) then cfile=kfile else cfile=0

nmods=n_elements(datmods(0,*))-1

write_ngong,lun,cdatas(*,0),nmod,nn,ndata,datmods(*,0),idatmds(*,0), $
  datgngs(*,0),bcs(*,0),yvars(*,*,0), file=cfile

n=1

while(n le nmods) do begin 
  write_ngong,lun,cdatas(*,n),nmod,nn,ndata,datmods(*,n),idatmds(*,n), $
    datgngs(*,n),bcs(*,n),yvars(*,*,n)
  n=n+1 
endwhile

return
end
