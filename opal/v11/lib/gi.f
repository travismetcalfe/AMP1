      subroutine gi(x,pp,f)
c same as GH v12 version
      implicit double precision(a-h,o-z)
      dimension f(4)
      h=1.d0-x
      f(1)=h
      f(2)=x
      f(3)=h*h*h/(pp*x+1.d0) 
      f(4)=x*x*x/(pp*h+1.d0) 
      return
      end
