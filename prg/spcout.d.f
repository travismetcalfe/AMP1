      subroutine spcout(x,y,yp,nn,iy,intr,agey,datout,ispcpr)
c
c  template for user-supplied routine for special output
c  x:  log q
c  y:  solution at actual time step
c  yp: solution at previous time step
c  nn: number of mesh points
c  intr: printing step used in other output
c  agey: age, in years
c  datout(1-31): global output parameters
c  ispcpr: user-supplied parameter
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
      implicit double precision (a-h,o-z)
      dimension x(nn),y(iy,nn),yp(iy,nn),datout(31)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  this is just a dummy routine
c
      if(istdpr.gt.0) write(istdpr,100) ispcpr
      return
  100 format(//' **** s/r spcout called with ispcpr =',i4,
     *  ' no action taken')
      end
