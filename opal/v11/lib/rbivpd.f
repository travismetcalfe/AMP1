      subroutine rbivpd(u,v,n,ndim,m,mdim,x,y,a,pp,z,zx,zy,iflag)
c  based on GH v12 version
      implicit double precision (a-h,o-z)
      real a
      dimension x(n),y(m),a(ndim,mdim,4,4)
      dimension g(4),h(4),hx(4),gy(4)
      data i,j /2*1/
c     if(n.lt.2.or.m.lt.2)then
c       iflag=1
c       return
c     endif
      call inttwo(x,n,y,m,u,v,i,j,iflag)
      if(iflag.ne.0)return
      a1=1.d0/(x(i+1)-x(i))
      a2=1.d0/(y(j+1)-y(j))
      ux=(u-x(i))*a1
      vy=(v-y(j))*a2
      z =0.d0
      zx=0.d0
      zy=0.d0
      call gi (ux,pp,h)
      call gid(ux,pp,hx)
      call gi (vy,pp,g)
      call gid(vy,pp,gy)
      do 20 k=1,4
         do 10 l=1,4
            z =z +dble(a(i,j,k,l))*h(k)*g(l)
            zx=zx+dble(a(i,j,k,l))*hx(k)*g(l)*a1
            zy=zy+dble(a(i,j,k,l))*h(k)*gy(l)*a2
c#dble#            z =z +a(i,j,k,l)*h(k)*g(l)
c#dble#            zx=zx+a(i,j,k,l)*hx(k)*g(l)*a1
c#dble#            zy=zy+a(i,j,k,l)*h(k)*gy(l)*a2
  10     continue
  20  continue
      return
      end
