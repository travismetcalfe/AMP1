      subroutine leqpr(a,b,nn,mm,ia,ib,err)   
c   
c  driver routine for leq, which prints equations and
c  solution if idgtnr in common/cdgtnr/ is .ge. 2
c   
c     nn - dimension of segment of a to be used  
c     mm - number of right hand columns of b to be used  
c     ia - the total number of rows in large array a
c     ib - the total number of rows in large array b
c
      implicit double precision (a-h, o-z)
      dimension a(ia,1), b(ib,1)
      common/cdgtnr/ idgtnr
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  test for printing equations
c
      if(idgtnr.ge.2.and.istdpr.gt.0) then
        write(istdpr,100)
	do 15 i=1,nn
   15   write(istdpr,110) (a(i,j),j=1,nn),(b(i,j),j=1,mm)
      end if
c 
      call leq(a,b,nn,mm,ia,ib,err)   
c
c  test for printing solution
c
      if(idgtnr.ge.2.and.istdpr.gt.0) then
        write(istdpr,120)
	do 25 i=1,nn
   25   write(istdpr,110) (b(i,j),j=1,mm)
      end if
      return
  100 format(' leq equations:')
  110 format(1p10e13.5)
  120 format(' leq solution:')
      end
