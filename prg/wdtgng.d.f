      subroutine wdtgng(ids,datgng,ndtgng,ihead)
c
c  outputs datgng on d/s ids.
c  If ihead = 1, set "#" at start of line, for header
c
c  Original version: 15/8/92
c
      implicit double precision(a-h,o-z)
      character*80 form1, form2
      dimension datgng(1)
c
      if(ihead.ne.1) then
        form1='(/''  datgng:'')'
        form2='(i3,'':'',1p5e13.5)'
      else
        form1='(''#''/''#  datgng:'')'
        form2='(''#'',i3,'':'',1p5e13.5)'
      endif
c
      if(mod(ndtgng,5).eq.0) then
	krec = ndtgng/5
      else
	krec = ndtgng/5+1
      end if
      write(ids,form1) 
      do 40 k=1,krec
      i1=1+5*(k-1)
      i2=min0(i1+4,ndtgng)
   40 write(ids,form2) i1,(datgng(i),i=i1,i2)
      return
      end
