pro read_ngsms,a,file=file,modpar=modpar, gsms=gsms, igsms=igsms, $
   nrdmax=nrdmax,swap=swap
;
;  If nrdmax is set, at most ndrmax modes are read, otherwise
;  the entire summary is read
;
;  reads grand summaries on file file
;  returns:
;  a(0,*): radius
;  a(1,*): l
;  a(2,*): order
;  a(3,*): sigma^2
;  a(4,*): nu (microHz)
;  a(5,*): E

;  If keyword modpar is set, in addition returns model parameters in
;  modpar:

;  modpar(0,*): mass
;  modpar(1,*): radius
;  modpar(2,*): central pressure
;  modpar(3,*): central density
;  modpar(4,*): scaled central second derivative of pressure
;  modpar(5,*): scaled central second derivative of density
;  modpar(6,*): polytropic index

;  If keyword gsms is set, in addition stores entire
;  real grand summary in gsms

;  If keyword igsms is set, in addition stores entire
;  integer grand summary in igsms

a=dblarr(6,10000)
if keyword_set(modpar) then modpar=dblarr(7,10000)
if keyword_set(gsms) then gsms=dblarr(38,10000)
if keyword_set(igsms) then igsms=lonarr(8,10000)
if keyword_set(swap) then swap = 1 else swap = 0

if keyword_set(nrdmax) then imax = nrdmax else imax = 10000

read_ngsm,1,cs,ics,file=file,swap=swap
i=-1

while cs(1) gt 0 and i lt imax do begin
  i=i+1
  a(0,i)=cs(2)
  a(1,i)=cs(17)
  a(2,i)=cs(18)
  a(3,i)=cs(20)
  a(4,i)=1000.*cs(36)
  a(5,i)=cs(23)
  if keyword_set(modpar) then modpar(0:6,i)=cs(1:7)
  if keyword_set(gsms) then gsms(*,i)=cs
  if keyword_set(igsms) then igsms(*,i)=ics
  read_ngsm,1,cs,ics,swap=swap
endwhile

if i ge 0 then begin
  a=a(*,0:i)
  if keyword_set(modpar) then modpar=modpar(*,0:i)
  if keyword_set(gsms) then gsms=gsms(*,0:i)
  if keyword_set(igsms) then igsms=igsms(*,0:i)
endif else begin
  a=a(*,0:0)
  a(0,0)=-1.d0
  if keyword_set(modpar) then modpar=modpar(*,0:0)
  if keyword_set(gsms) then gsms=gsms(*,0:0)
  if keyword_set(igsms) then igsms=igsms(*,0:0)
endelse
close,1
return
end
