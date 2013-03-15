      subroutine wdtmod(ids,datmod,nrdtmd,ndtmod,nidtmd,ihead)
c
c  outputs datmod and ndtmod on d/s ids.
c  If ihead = 1, set "#" at start of line, for header
c
c  Original version: 15/8/92
c
      implicit double precision(a-h,o-z)
      character*80 form1, form2, form3, form4
      dimension datmod(1),ndtmod(1)
c
      if(ihead.ne.1) then
        form1='(/''  datmod:'')'
        form2='(i3,'':'',1p5e13.5)'
        form3='(/'' ndtmod:'')'
        form4='(i3,'':'',5i13)'
      else
        form1='(''#''/''#  datmod:'')'
        form2='(''#'',i3,'':'',1p5e13.5)'
        form3='(''#'',/''# ndtmod:'')'
        form4='(''#'',i3,'':'',5i13)'
      endif
c
      if(mod(nrdtmd,5).eq.0) then
	krec = nrdtmd/5
      else
	krec = nrdtmd/5+1
      end if
      write(ids,form1) 
      do 30 k=1,krec
      i1=1+5*(k-1)
      i2=min0(i1+4,nrdtmd)
   30 write(ids,form2) i1,(datmod(i),i=i1,i2)
c
      if(mod(nidtmd,5).eq.0) then
	krec = nidtmd/5
      else
	krec = nidtmd/5+1
      end if
      write(ids,form3) 
      do 35 k=1,krec
      i1=1+5*(k-1)
      i2=min0(i1+4,nidtmd)
   35 write(ids,form4) i1,(ndtmod(i),i=i1,i2)
c
      return
      end
