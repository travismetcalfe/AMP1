      subroutine sclcbc(bc,datmod)
c
c  rescales coefficients in central boundary conditions,
c  following redifinition of variables on 20/8 1984.
c  rescaling is done if icase .lt. 1.d7 and bc(4) (the
c  central temperature variable) .gt. 1.d4.
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
      implicit double precision (a-h,o-z)
      dimension bc(*),datmod(*),gm(5)
      equivalence (dt,idt)
c
      save
c
      data gm /1.e-17,1.e-7,1.,1.,1./
c
c  test for rescaling
c
      dt=datmod(29)
      i7=idt/10 000 000
      if(i7.ne.0.or.bc(4).le.1.d4) return
c
c  rescale
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
c  reset icase
c
   50 idt=idt+10 000 000
      datmod(29)=dt
c
      return
      end
