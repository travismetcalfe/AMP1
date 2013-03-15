      subroutine wremdl(ids,x,y,datmod,ndtmod,bccoef,iy,iform,
     *  nn,nrdtmd,nidtmd,nbccf,nvar)
c
c  Writes evolution model file on unit ids, assumign new format
c
c  Originial version: 20/7/92
c
      implicit double precision(a-h, o-z)
      dimension x(1), y(iy,1), datmod(1), ndtmod(1), bccoef(1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      write(ids) iform,nn,nrdtmd,nidtmd,nvar,nbccf,
     *  (datmod(i),i=1,nrdtmd),(ndtmod(i),i=1,nidtmd),
     *  (x(n),(y(i,n),i=1,nvar),n=1,nn),
     *  (bccoef(i),i=1,nbccf)
c
      return
      end
