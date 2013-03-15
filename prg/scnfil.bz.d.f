      subroutine scnfil(ids,iw,ic)
c  list evolution models
c 
c  corrected 16/12/85 for rescaling of central boundary
c  condition quantities
c 
c  Modified 16/3/89, storing x and y in commons /sooner/ and /work/
c     Note: these now have to be dimensioned in calling programme.
c
c  Modified 13/8/91, to take into account new output format.
c
c
      implicit double precision (a-h,o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
      dimension datmod(nrdtmd),bc(nbcprv),ndtmod(nidtmd)
      common/sooner/ x(1)
      common/work/ y(ivarmx,1)
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev, syear
c
      save
c
      rewind ids
      go to (1,2), ic
    1 write(iw,100) ids
      go to 5
    2 write(iw,102) ids
c
    5 nmod=0
c
   10 read(ids,end=90) iform,nn,nrdtmr,nidtmr,nvarrd,nbccf,
     *  (datmod(i),i=1,nrdtmr),(ndtmod(i),i=1,nidtmr),
     *  (x(n),(y(i,n),i=1,nvarrd),n=1,nn),
     *  (bc(i),i=1,nbccf)
      nmod=nmod+1
      tef=10.**y(3,1)
      ams=datmod(23)/amsun
c
c  test for case
c
      go to (20,30), ic
c
   20 pc=bc(3)*1.d17  
      tc=bc(4)*1.d7
      xc=bc(5)
      write(iw,110) nmod,ams,datmod(2),datmod(1),datmod(22),datmod(24),
     *  tef,datmod(25),pc,tc,xc,ndtmod(1)  
      go to 10
c  test for setting of patmos
c
c  note that as a result of resetting  of datmod (23/11/82), datmod(21)
c  is non-zero and patmos is in datmod(20) instead of in datmod(19)
c
   30 iptmos=19
      if(datmod(21).ne.0) iptmos=20
      patmos=datmod(iptmos)
      write(iw,120) nmod,ams,datmod(22),datmod(2),datmod(1),datmod(3),
     *  patmos,datmod(24),datmod(25),ndtmod(1)
      go to 10
   90 continue
      if(nmod.eq.0) write(iw,190) ids
      return
  100 format('1evolution sequence on d/s ',i3//
     *  ' nmod,mass,surface x,z,age,radius,',
     *  'eff. temperature,luminosity,pc,tc,xc,icase:'/)
  102 format('1evolution sequence on d/s ',i3//
     *  ' nmod,mass,age,surface x,z,alfa,patmos,radius,',
     *  'luminosity,icase:'/)
  110 format(i4,f11.6,2f10.6,1p6e12.5,0pf10.6,i10)
  120 format(i4,f12.6,1pe13.5,0p3f10.6,1p3e15.7,i10)
  190 format(//1x,10('*'),' no models on d/s',i3)
      end
