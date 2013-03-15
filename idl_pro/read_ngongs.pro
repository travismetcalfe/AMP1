pro read_ngongs,lun,cdatas,nmod,nn,ndata,datmods,idatmds,datgngs,bcs,yvars, $
   file=kfile, number=knum,nmax=nmax,swap=swap
;  reads unformatted GONG models from unit lun.
;  if file is set, opens file.
;  if the unit number is set in lun, this is used,
;  otherwise lun is set to 1.
;  If new is set, force read of new format, otherwise test 
;  on data
;  If swap is set, use swap_endian, for reading, e.g., linux files
;  from bigcat

;  Original version: 6/7/00.

common numgong, nread

if keyword_set(kfile) then cfile=kfile else cfile=0

if keyword_set(swap) then swap = 1 else swap = 0

read_ngong,lun,cdata,nmod,nnr,ndata,datmod,idatmd,datgng,bc,yvar, $
   file=cfile,swap=swap

if keyword_set(nmax) then n_max=nmax else n_max=500

nn=nnr

n=0
cdatas=strarr(4,n_max)
yvars=dblarr(ndata(3),nn,n_max)
datmods=dblarr(ndata(0),n_max)
idatmds=lonarr(ndata(1),n_max)
datgngs=dblarr(ndata(2),n_max)
bcs=dblarr(ndata(4),n_max)

while(nnr gt 0 and n lt n_max) do begin 
  cdatas(*,n)=cdata
  yvars(*,*,n)=yvar 
  datmods(*,n)=datmod 
  idatmds(*,n)=idatmd 
  datgngs(*,n)=datgng 
  bcs(*,n)=bc
  read_ngong,lun,cdata,nmod,nnr,ndata,datmod,idatmd,datgng,bc,yvar,/quiet, $
    swap=swap
  n=n+1 
endwhile


yvars=yvars(*,*,0:n-1)
datmods=datmods(*,0:n-1)
idatmds=idatmds(*,0:n-1)
datgngs=datgngs(*,0:n-1)

return
end
