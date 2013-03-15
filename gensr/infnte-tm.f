      integer function infnte(x)
c
c  test whether the argument is a valid floating-point number.
c  In that case, returns 0. Otherwise returns 1,  and (for now) 
c  writes a message
c
c  Original version: 7/7/00
c
c  Modified 12/7/02, using Travis Metcalfe trick
c
      implicit double precision (a-h, o-z)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      if(x.eq.x) then
	infnte=0
      else
	if(istdpr.gt.0) write(istdpr,100) x
	infnte=1
      end if
      return
  100 format(/' ****** Argument to infnte is ',e13.5)
      end
