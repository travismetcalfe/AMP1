      subroutine gid(x,pp,f)
c  same as GH v12 version
      implicit double precision(a-h,o-z)
      dimension f(4)
      h=x-1.d0
      g=-pp+pp*x-1.d0
      e=pp*x+1.d0
      f(1)=-1.d0
      f(2)= 1.d0
      f(3)=-h*h*( 3.d0+pp+2.d0*pp*x)/(e*e)
      f(4)=-x*x*(-3.d0*pp+2.d0*pp*x-3.d0)/(g*g)
      return
      end
