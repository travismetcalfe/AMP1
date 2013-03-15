      character*(*) function str_replace(s, s1, s2, ierr)
c
c  replaces first occurrence of string s1 in string s by string s2
c  ignoring trailing blanks
c
c  If s1 is not found in s, ierr is returned as -1, and the function
c  returns s.
c
c  Original version: 11/2/06
c
      character*(*) s, s1, s2
      character*120 srest, snew
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  find possible location of s1 in s
c
      ls =length(s)
      ls1=length(s1)
      ls2=length(s2)
c
      snew=s
      lsnew=ls
      ierr=0
      do i=1,ls
	if(s(i:i+ls1-1).eq.s1(1:ls1)) then
	  srest=s(i+ls1:ls)
	  snew(i:i+ls2-1)=s2(1:ls2)
	  lsnew=ls+ls2-ls1
	  snew(i+ls2:lsnew)=srest
	  go to 20
        end if
      end do
c
c  diagnostic output
c
      ierr=-1
      write(istdou,100) s1(1:ls1)
      write(istdou,110) s(1:ls)
      if(istdpr.gt.0.and.istdpr.ne.istdou) then
        write(istdpr,100) s1(1:ls1)
        write(istdpr,110) s(1:ls)
      end if
c
   20 str_replace=snew(1:lsnew)
      return
  100 format(/' ***** Warning. String ',a)
  110 format( '   not found in string ',a)
      end
