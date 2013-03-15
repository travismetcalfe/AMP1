      program runadi
c
c  main programme for adiabatic pulsations
c
c  quantities set in parameter statement:
c  nnmax: maximum number of mesh points in integration
c
c  Note: as presently set up, nnmax is set also in s/r nrkint and
c  eigin4
c
c  Modified 20/1/95 to increase first dimension of aa to 10.
c  Note: with present version of code this must also be implemented
c  in s/r adipls and rhs
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      parameter (nnmax = 10000, iaa = 10, iaa1 = 10, iy = 8)
      parameter (nnmax1 = nnmax+1, nnmax2 = nnmax+10,
     *  nnwwin=20*nnmax2)
      common/rhsdat/ dt1(29),aa(iaa,nnmax) /xarra/ x(nnmax)
     *  /xarr1/ x1(nnmax1) /xarr2/ x2(nnmax1)
     *  /worksp/ aa1(iaa1,nnmax)   /yyyyyy/ y(iy,nnmax)  
     *  /yyyyri/ yri(4,nnmax)
     *  /sysord/ sysyso(1)
     *  /work/ wwnrk(20,nnmax2)
      common/nrmchk/ nnww, ncfac, irdist
      common/wrkleq/ wwwlll(1500)
      common/cderst/ derc(6,nnmax)  /cintst/ aintc(6,nnmax)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      write(istdou,100) nnmax
c
c  for consistency with previous usage, set nnww to 0 for the moment.
c  A more reasonable value would be
c     nnww = nnwwin
c
      nnww = 0
c
      call adipls
      stop
  100 format(//61('*')//
     *  ' In this version, the maximum number of mesh points is',i7//
     *  61('*'))
      end
      subroutine leqdet(a,b,nn,mm,ia,ib,err)
c
c  Double precision version.
c  +++++++++++++++++++++++++
c
c  Dated: 10/3/90
c
      implicit double precision (a-h,o-z)
      dimension a(100),b(100)
c
      call leq(a,b,nn,mm,ia,ib,err)
      return
      end
      subroutine spcout
      return
      end
      subroutine modmod
      return
      end
      subroutine resdif
      return
      end
