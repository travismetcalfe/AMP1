      program opares
c
c     compares opacities by reading a model in GONG format
c     and calling the opacity interpolation routines again
c     by using the models temperature and density
c     (see X22::D:/hg/opacity/jcd/opares.f)
c
c     witten first in May 1997 when comparing JCD's model with
c     models by the Munich people (garsom5.bin) used by
c     Armin Weiss and Schattel (?)
c
c     hg: 07/10/02: adopted for using 'my' (HG's) GONG format
c                   (as used in s/r gongot in programme em_v3.7v9b.f)
c                   and comparing OPAL95 with Frank Piper's low-Z
c                   tables
c
      implicit double precision(a-h,o-z)
c
      character*12 tabnam
      character*80 cdata(4)
      dimension datmod(100),datgng(20),idatmd(20),bccoef(60)
      dimension var(30,4000),opalg(4000),residu(4000)
c
      open(10,file='gong.dat',form='unformatted')
      read(10)(cdata(j),j=1,4),nmod,iform,nn,nrdtmd,nidtmd,ndtgng,
     +      nvar,nbccf,(datmod(j),j=1,nrdtmd),(idatmd(j),j=1,nidtmd),
     +      (datgng(j),j=1,ndtgng),(bccoef(j),j=1,nbccf),
     +      ((var(j,k),j=1,nvar),k=nn,1,-1)
      close(10)
c     print *,'nn=',nn
c     print *,'nvar=',nvar
c
      print *,cdata(1)
      print *,cdata(2)
      print *,cdata(3)
      print *,cdata(4)
      print *,'Z=     ',datmod(1)
      print *,'X=     ',datmod(2)
      print *,'alpha= ',datmod(3)
      print *,'Mass=  ',datmod(23)
      print *,'Radius=',datmod(24)
      print *,'Lumino=',datmod(25)
      print *,''
c
c     just for Sun-mode, read Z
c     open(13,file='fgong.Z',form='unformatted')
c     read(13)(var(17,j),j=1,nn)
c     close(13)
c
      data iexp/0/
      data tabnam /'OPINTPATH_02'/
      data iorder /4/
      data imode /-2/          ! splines with electron conduction disabled
c
c     initialize optables (read them)
      call maceps(drelpr)
      print *,'initialize tables'
      call opinit(drelpr,iorder,tabnam,imode)
c
      print *,'compute residuals...'
c     print *,'zval,xval,rlg,tlg,val,opa'
      do i=1,nn
        tlg=var(4,i)
        rlg=dlog10(var(10,i))-3.d0*tlg+18.d0
        xval=var(6,i)
        zval=datmod(1)
        call opints(xval,zval,tlg,rlg,opalg(i),opr,opt,opx,opz,iexp,ier)
        residu(i)=1.0d1**opalg(i)-var(12,i)
        print'(i5,4f9.6,1p3e12.5)',i,zval,xval,rlg,tlg,
     #        var(12,i),1.0d1**opalg(i),residu(i)/var(12,i)
      enddo  
c
c     write results
      open(20,file='residuals.bin',form='unformatted')
      write(20)(var(4,i),i=1,nn),
     #         (var(12,i),i=1,nn),
     #         (opalg(i),i=1,nn),
     #         (residu(i),i=1,nn)
      close(20)
c
      stop
      end
