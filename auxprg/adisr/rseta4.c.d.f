      subroutine rseta4(x,aa,nn,data,iaa)   
c  
c  Reset A4 near centre, to correct for problems near
c  end of hydrogen burning
c  Resetting is only applied if A4/x**2 is non-monotonic
c  near centre
c
c  Orignial version: 25/7/92
c  
c  Modified 8/6/03, suppressing resetting in convective core. In this
c  case, also reset data(6) to ensure consistency.
c
      implicit double precision (a-h, o-z)
      dimension x(*),aa(iaa,*),data(*)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  reference values
c
      axc=data(6)-data(5)
      nref=5
      axref=aa(4,nref)/x(nref)**2
c
      ireset=0
      do 10 n=2,nref-1
      ax=aa(4,n)/x(n)**2
      if((axc-ax)*(ax-axref).lt.0) ireset=1
   10 continue
c
      if(ireset.eq.1) then
c
c  suppress resetting in case of convective core
c
	if(axref.lt.0) then
	  write(istdou,105) 
	  if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,105)
	  datan6=data(5)+axref
	  write(istdou,107) data(6), datan6
	  if(istdou.ne.istdpr.and.istdpr.gt.0) 
     *      write(istdpr,107) data(6), datan6
	  data(6)=datan6
	  return
        end if
c
	write(istdpr,110)
	axcoef=(axref-axc)/x(nref)**2
	do 20 n=2,nref-1
	x2=x(n)**2
	ax=aa(4,n)/x2
	axnew=axc+axcoef*x2
	write(istdpr,115) n, x(n), ax, axnew
   20   aa(4,n)=axnew*x2
c
      end if
c  
      return   
  105 format(//
     *  ' Resetting of A4 near centre suppressed in convective core')
  107 format(/ 'Reset data(6). Old, new values =',1p2e13.5)
  110 format(//' Reset A4 near centre. ',
     *   ' n, x, old, new values of A4/x**2:'/)
  115 format(i5,1p3e13.5)
      end  
