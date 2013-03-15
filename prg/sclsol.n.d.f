      subroutine sclsol(y,nn,iy,bc,datmod,ndtmod)
c
c  rescales coefficients in central boundary conditions,
c  following redifinition of variables on 20/8 1984.
c  rescaling is done if i7 = int(icase/1.d7) = 0 and bc(4) (the
c  central temperature variable) .gt. 1.d4.
c
c  also resets dependent variables y(1,n) and y(4,n) to be
c  log(r/1.d11) and log(l/1.d33), respectively, if i7 .le. 1.
c  i7 is reset to 2 to flag for change of variables.
c
c  original version from 4/1/1985.
c
c  Modified 13/8/91, for change in input/output formats
c
c  Dated 11/3/90
c
      implicit double precision (a-h,o-z)
      dimension y(iy,nn),bc(*),datmod(*),ndtmod(*),gm(5)
c
      save
c
      data gm /1.e-17,1.e-7,1.,1.,1./
c
c  test for rescaling
c
      icase=ndtmod(1)
      i7=icase/10 000 000
      if(i7.ne.0.or.bc(4).le.1.d4) go to 50
c
c  rescale boundary condition coefficients
c
      bc(2)=1.d22*bc(2)
c
      do 10 i=3,7
      i1=i-2
      i2=i+5
      bc(i)=gm(i1)*bc(i)
   10 bc(i2)=1.d22*gm(i1)*bc(i2)
c
      bc(13)=1.d5*bc(13)
      bc(14)=1.d15*bc(14)
c
      do 15 i=15,17
      gmi=1./gm(i-14)
      i2=i+3
      bc(i)=1.d15*(gmi*bc(i))
   15 bc(i2)=gmi*bc(i2)
c
      do 20 i=1,5
      gmi=1./gm(i)
      i1=15+i
      do 20 ip=1,3
      i1=i1+5
   20 bc(i1)=gmi*bc(i1)
c
      bc(37)=1.d22*bc(37)
c
      do 25 i=38,40
      gmi=1./gm(i-37)
      i2=i+3
      i3=i+6
      bc(i)=gmi*bc(i)
      bc(i2)=1.d22*(gmi*bc(i2))
   25 bc(i3)=gmi*bc(i3)
c
      i1=46
      do 30 ia=1,6
      do 30 i=1,5
      i1=i1+1
      if(i1.eq.55) go to 50
   30 bc(i1)=gm(i)*bc(i1)
c
c  test for resetting y(1,n) and y(4,n)
c
   50 if(i7.ge.2) return
c
      do 55 n=1,nn
      y(1,n)=y(1,n)-11
   55 y(4,n)=y(4,n)-33
c
c  reset icase
c
   60 icase=icase+(2-i7)*10 000 000
      ndtmod(1)=icase
c
      return
      end
