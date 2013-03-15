      program main
c
c  finds differences between models on GONG format.
c
c  original version: 29/12/87
c
c
c  Important note: Before 3/8/88, calculated delta GAMMA1. Also
c  **************  delta ln c was set incorrectly.
c
c  Modified 13/6/89, to bring it in line with diff-fgong. Specifically
c   - Now allows comparisons at fixed r/R, q or r/Rs
c   - Uses same special treatment of innermost meshpoints.
c   - computes delta (ln ad. grad) and delta (ln delta) instead of
c     delta (ad. grad.) and delta (delta)
c
c  Modified 5/2/91 (for the third time!) to include difference of
c  Gamma = d ln p / d ln rho. Also added log q to output
c  Changed delta A output from delta (ln A) to (delta A)/Vg.
c
c  Modified 12/2/91, to allow differences at fixed pressure and
c  at fixed m/Ms (note that the normal q is normalized at the
c  photosphere).
c
c  Modified 13/8/91, to take into account changed format of
c  GONG models.
c
c  Modified to allow differences between models in same sequence
c
c  Modified 3/8/92, to disallow extrapolation for icase = 2
c  (quick fix that should be reconsidered carefully)
c
c  Modified 31/5/94 to allow differences at fixed m/M(fixed r/R).
c
c  Modified 12/5/97, to disallow extrapolation for icase = 5
c  (quick fix that should be reconsidered carefully)
c
c
c  Double precision version
c  ++++++++++++++++++++++++
c
      implicit double precision (a-h, o-z)
      real yvarrs, datgns, bcs, datmrs
      include 'evolpr.incl'
      parameter(ndiff=20)
      character*80 cdata,fin1,fin2, fout, fgjd
      character dyname*20, ccase*4, cprec*6
      dimension idss(2),nmodel(2),iprec(2),datmod(nrdtmx,2),
     *  yvar(igvrmx,nnmax,2),cprec(2),cdata(4,2),
     *  yvarrd(igvrmx,nnmax),datgng(ndtgmx,2),bc(nbccmx,2),
     *  yvarrs(20,nnmax),datgns(30),bcs(54),
     *  nvarrd(2),nnrd(2),
     *  xab1(nnmax),xab2(nnmax),xr1(nnmax),xr2(nnmax),
     *  xq1(nnmax),xq2(nnmax),vg1(nnmax),vg2(nnmax),
     *  xabc1(nnmax),xabc2(nnmax),yint(igvrmx),yvarc1(igvrmx,20),
     *  dy(ndiff),xmax(ndiff),dymin(ndiff),dymax(ndiff),daymax(ndiff),
     *  dyname(ndiff),ccase(7),
     *  datmrs(31),idatrs(3),idatmd(nidtmx,2),
     *  nrdtmd(2),nidtmd(2),ndtgng(2),nbccf(2)
c
      equivalence(datmrs(29),idatrs(1))
c
      data dyname /
     *  'delta(ln q)','delta(ln f)',
     *  'delta(ln T)','delta(ln L)',
     *  'delta(X)','delta(ln p)',
     *  'delta(ln rho)','delta(ln kappa)',
     *  'delta(ln epsilon)','delta(ln GAMMA1)',
     *  'delta(ln ad. grad)','delta(ln delta)',
     *  'delta(ln cp)','(delta A)/Vg',
     *  'delta(ln c)','delta(ln Ne)',
     *  'delta(ln rX)','delta(ln r)',
     *  'delta(ln GAMMA)','delta(ln X3)' /
      data ccase /'r/R', 'q', 'r/Rs', 'r', 'p', 'm/Ms', 'm/Mr'/
      data cprec /'single', 'double'/
c
      powten(x)=10.d0**(min(max(x,-36.d0),36.d0))
c
      cgrav=6.67232e-8
c
      write(6,*) 'Enter first file name'
      read(5,'(a)') fin1
      open(2,file=fin1,status='old',form='unformatted')
      idss(1) = 2
c
      write(6,*) 'Enter second file name'
      read(5,'(a)') fin2
c
c  test for differences in same model sequence
c
      if(fin2.ne.fin1) then
        open(3,file=fin2,status='old',form='unformatted')
	idss(2)=3
      else
	idss(2)=2
      end if
c
      nmodel(1)=1
      nmodel(2)=1
      write(6,*) 'Enter model numbers on first and second file'
      read(5,*) nmodel
c
      iprec(1)=2
      iprec(2)=2
      write(6,*) 'Enter precision of first and second model'
      write(6,*) '1 for single precision, 2 for double precision'
      read(5,*) iprec
c
      write(6,100) nmodel(1),fin1,nmodel(2),fin2
c
      write(6,*) 'Enter output file'
      read(5,'(a)') fout
      open(10,file=fout,status='unknown')
      write(10,102) nmodel(1),fin1,nmodel(2),fin2
c
c  file for summary
c
      open(12,file='dgong-sum',status='unknown')
      write(6,103)
      write(12,100) nmodel(1),fin1,nmodel(2),fin2
c
      write(6,*) 'Enter 1 for differences at fixed r/R, 2 for fixed q,'
      write(6,*) '3 for differences at fixed r/r(last point),'
      write(6,*) '4 for differences at fixed r'
      write(6,*) '5 for differences at fixed p'
      write(6,*) '6 for differences at fixed m/M(last point)'
      write(6,*) '7 for differences at fixed m/M(fixed r/R)'
      icase=1
      read(5,*) icase
c
      if(icase.ne.7) then
        write(6,104) ccase(icase)
        write(10,105) ccase(icase)
        write(12,104) ccase(icase)
      else
	write(6,*) 'Enter r/R for point of fixing mass'
	read(5,*) rmsfix
        write(6,106) rmsfix
        write(10,107) rmsfix
        write(12,106) rmsfix
      end if
c
      write(6,*) 'Enter 2 for diagnostics from interpolation'
      itest=0
      read(5,*) itest
c
      write(6,*) 'Enter 1 for output on Guenther et al form'
      read(5,*) igjdpr
      if(igjdpr.eq.1) then
        write(6,*) 'Enter file name for Guenther et al output'
        read(5,'(a)') fgjd
        open(14,file=fgjd,status='unknown')
        write(14,102) nmodel(1),fin1,nmodel(2),fin2
	if(icase.ne.7) then
          write(14,105) ccase(icase)
        else
          write(14,107) rmsfix
	end if
        write(14,108)
      end if
c
      do 25 ids1=1,2
      ids=idss(ids1)
      if(ids1.eq.1.or.ids.eq.3) then
        nrd=0
        icasrd=0
      end if
c
   10 if(iprec(ids1).eq.1) then
c
c  single precision read. Assume old format
c
	read(ids,end=25) (cdata(i,ids1),i=1,4),
     *    nmod,nn,nvar,(datmrs(i),i=1,31),
     *    (datgns(i),i=1,30),(bcs(i),i=1,54),
     *    ((yvarrs(i,n),i=1,nvar),n=1,nn)
c
        do 11 i=1,28
   11   datmod(i,ids1)=datmrs(i)
        do 12 i=1,3
   12   idatmd(i,ids1)=idatrs(i)
c
	do 13 i=1,30
   13   datgng(i,ids1)=datgns(i)
	do 14 i=1,54
   14   bc(i,ids1)=bcs(i)
c
      else
c
        call rdgong(ids,cdata(1,ids1),nmod,datmod(1,ids1),
     *    idatmd(1,ids1),datgng(1,ids1),bc(1,ids1),yvarrd,
     *    igvrmx,iform,nn,nrdtmd(ids1),nidtmd(ids1),ndtgng(ids1),
     *    nbccf(ids1),nvar,icasrd,icry)
        if(icry.eq.-1) then
	  go to 25
        else if(icry.eq.-2) then
	  stop 'Read error'
        end if
c
	icasrd = -1
c
c
      end if
c
c  store such that r increases with n
c  storage depends on precision parameter
c
      if(iprec(ids1).eq.1) then
        if(yvarrs(1,1).lt.yvarrs(1,nn)) then
          do 18 n=1,nn
          do 18 i=1,nvar
   18     yvar(i,n,ids1)=yvarrs(i,n)
	else
          do 19 n=1,nn
          n1=nn+1-n
          do 19 i=1,nvar
   19     yvar(i,n1,ids1)=yvarrs(i,n)
	end if
      else
        if(yvarrd(1,1).lt.yvarrd(1,nn)) then
          do 20 n=1,nn
          do 20 i=1,nvar
   20     yvar(i,n,ids1)=yvarrd(i,n)
	else
          do 21 n=1,nn
          n1=nn+1-n
          do 21 i=1,nvar
   21     yvar(i,n1,ids1)=yvarrd(i,n)
	end if
      end if
c
      nvarrd(ids1)=nvar
c
c  for icase = 7, chop model to stop at fixed radius
c
      if(icase.eq.7) then
	nrfix=-1
	do 22 n=1,nn
        if(yvar(2,n,ids1).le.-35) then
          xr1(n)=0
        else 
          xr1(n)=10.d0**yvar(2,n,ids1)/(1.e-11*datgng(2,ids1))
	end if
	if(nrfix.eq.-1.and.xr1(n).ge.rmsfix) nrfix=n
   22   continue
c
        if(nrfix.eq.-1) then
	  write(6,110) rmsfix, xr1(nn)
	  stop
        end if
c
c  test for interpolating to set last point
c
        if(xr1(nrfix).gt.rmsfix) then
	  call lir(rmsfix,xr1,yint,yvar(1,1,ids1),nvar,igvrmx,nn,
     *             1,inter)
	  do 23 i=1,nvar
   23     yvar(i,nrfix,ids1)=yint(i)
        end if
c
	nn=nrfix
      end if
c
      nnrd(ids1)=nn
      nrd=nrd+1
      if(nrd.lt.nmodel(ids1)) go to 10
c
   25 continue
c
      nvar=min0(nvarrd(1),nvarrd(2))
      if(nvar.le.18) then
        ndiffa=15
      else
        ndiffa=ndiff
      end if
c
c  reset log q and log L, and set r/R for the two models
c  also set Vg and Gamma
c
      fctvg1=cgrav*datmod(23,1)/datmod(24,1)
      fctvg2=cgrav*datmod(23,2)/datmod(24,2)
c
      cq=2
      dq=1./log(10.d0)
c
      do 28 n=1,nnrd(1)
      if(yvar(2,n,1).le.-35) then
        xr1(n)=0
      else if(icase.ne.3) then
        xr1(n)=10.d0**yvar(2,n,1)/(1.e-11*datgng(2,1))
      else
        xr1(n)=10.d0**(yvar(2,n,1)-yvar(2,nnrd(1),1))
      end if
      if(icase.eq.1.or.icase.eq.3) then
        xab1(n)=xr1(n)
        xabc1(n)=xr1(n)*xr1(n)
      else if(icase.eq.4) then
        xab1(n)=10.d0**yvar(2,n,1)
        xabc1(n)=xab1(n)*xab1(n)
      else if(icase.eq.5) then
        xab1(n)=-log10(yvar(9,n,1))
      else
c
c  test for differences at fixed q or fixed m/Ms
c
        if(icase.eq.2.or.xr1(n).eq.0) then
          xlogq=yvar(1,n,1)
        else
          xlogq=yvar(1,n,1) - yvar(1,nnrd(1),1)
        end if
c..        write(6,*) xlogq
        if(xr1(n).gt.0) then
          xab1(n)=
     *      (powten(xlogq/3)+cq)*xlogq/(dq-xlogq)
          xabc1(n)=powten(0.6666667*xlogq)
        else
          xab1(n)=-cq
          xabc1(n)=0
        end if
      end if
c
      xq1(n)=yvar(1,n,1)
      yvar(1,n,1)=powten(yvar(1,n,1))
      yvar(2,n,1)=1.e11*powten(yvar(2,n,1))
      if(xr1(n).gt.0) then
        vg1(n)=fctvg1*yvar(1,n,1)*yvar(10,n,1)/
     *         (yvar(14,n,1)*yvar(9,n,1)*xr1(n))
        yvar(21,n,1)=yvar(14,n,1)/(1+yvar(18,n,1)/vg1(n))
      else
        vg1(n)=0
        yvar(21,n,1)=0
      end if
c..      write(79,*) xr1(n), vg1(n)
   28 continue
c
      do 30 n=1,nnrd(2)
      if(yvar(2,n,2).le.-35) then
        xr2(n)=0
      else if(icase.ne.3) then
        xr2(n)=10.d0**yvar(2,n,2)/(1.e-11*datgng(2,2))
      else
        xr2(n)=10.d0**(yvar(2,n,2)-yvar(2,nnrd(2),2))
      end if
      if(icase.eq.1.or.icase.eq.3) then
        xab2(n)=xr2(n)
        xabc2(n)=xr2(n)*xr2(n)
      else if(icase.eq.4) then
        xab2(n)=10.d0**yvar(2,n,2)
        xabc2(n)=xab2(n)*xab2(n)
      else if(icase.eq.5) then
        xab2(n)=-log10(yvar(9,n,2))
      else
c
c  test for differences at fixed q or fixed m/Ms
c
        if(icase.eq.2.or.xr2(n).eq.0) then
          xlogq=yvar(1,n,2)
        else
          xlogq=yvar(1,n,2) - yvar(1,nnrd(2),2)
        end if
c..        write(6,*) xlogq
        if(xr2(n).gt.0) then
          xab2(n)=
     *      (powten(xlogq/3)+cq)*xlogq/(dq-xlogq)
          xabc2(n)=powten(0.6666667*xlogq)
        else
          xab2(n)=-cq
          xabc2(n)=0
        end if
      end if
      xq2(n)=yvar(1,n,2)
      yvar(1,n,2)=powten(yvar(1,n,2))
      yvar(2,n,2)=1.e11*powten(yvar(2,n,2))
      if(xr2(n).gt.0) then
        vg2(n)=fctvg2*yvar(1,n,2)*yvar(10,n,2)/
     *         (yvar(14,n,2)*yvar(9,n,2)*xr2(n))
        yvar(21,n,2)=yvar(14,n,2)/(1+yvar(18,n,2)/vg2(n))
      else
        vg2(n)=0
        yvar(21,n,2)=0
      end if
   30 continue
c
c  set internal value of nvar, allowing for variables added
c
      nvarl=21
c
c  set special variables in innermost part of model 1, for
c  special interpolation
c
      do 32 n=1,20
      do 32 i=1,nvarl
   32 yvarc1(i,n)=yvar(i,n,1)
c
      if(icase.eq.1.or.icase.eq.3.or.icase.eq.4) then
        yvarc1(1,1)=(4.1887902e33/datmod(23,1))*yvar(10,1,1)
        yvarc1(11,1)=4.1887902*yvar(10,1,1)*yvar(13,1,1)
        yvarc1(18,1)=1.e22*(datgng(13,1)/(yvar(9,1,1)*yvar(14,1,1))
     *    -datgng(14,1)/yvar(10,1,1))
        do 34 n=2,20
        yvarc1( 1,n)=1.e33*yvar( 1,n,1)/yvar(2,n,1)**3
        yvarc1( 2,n)=yvar(2,n,1)**2
        yvarc1(11,n)=      yvar(11,n,1)/yvar(2,n,1)**3
   34   yvarc1(18,n)=1.e22*yvar(18,n,1)/yvar(2,n,1)**2
c
c..        write(6,122) (n,xabc1(n),yvarc1(1,n),yvarc1(11,n),yvarc1(18,n),
c..     *    n=1,20)
c
      else if(icase.eq.2.or.icase.eq.6) then
        yvarc1(2,1)=
     *    1.e-11*(datmod(23,1)/(4.1887902*yvar(10,1,1)))**0.3333333
        yvarc1(11,1)=1.e-33*datmod(23,1)*yvar(13,1,1)
        yvarc1(18,1)=(datgng(13,1)/(yvar(9,1,1)*yvar(14,1,1))
     *    -datgng(14,1)/yvar(10,1,1))*
     *    (datmod(23,1)/(4.1887902*yvar(10,1,1)))**0.6666667
        do 36 n=2,20
        yvarc1(1,n)=yvar(1,n,1)**0.6666667
        yvarc1(2,n)=1.e-11*yvar(2,n,1)/yvar(1,n,1)**0.3333333
        yvarc1(11,n)=1.e-33*yvar(11,n,1)/yvar(1,n,1)
   36   yvarc1(18,n)=      yvar(18,n,1)/xabc1(n)
c
c..        write(6,123) (n,xabc1(n),yvarc1(2,n),yvarc1(11,n),yvarc1(18,n),
c..     *    n=1,20)
      end if
c
      do 45 ids1=1,2
      write(10,113) ids1, cprec(iprec(ids1))
      write(10,115) (cdata(i,ids1),i=1,4)
      write(10,116) (datmod(i,ids1),i=1,nrdtmd(ids1))
      write(10,117) (idatmd(i,ids1),i=1,nidtmd(ids1))
      write(10,118) (datgng(i,ids1),i=1,ndtgng(ids1))
   45 write(10,119) (bc(i,ids1),i=1,nbccf(ids1))
c
c  step through model 2, interpolating in model 1
c
      amm=log(10.d0)
      write(10,150) ndiffa,(i,dyname(i),i=1,ndiffa)
      write(10,155)
c
      nint=0
c
      do 60 n=1,nnrd(2)
c
c  interpolation
c
      if(icase.ne.5.or.xab2(n).gt.xab1(4)) then
        nint=nint+1
        call lir(xab2(n),xab1,yint,yvar(1,1,1),nvarl,igvrmx,nnrd(1),
     *    nint,inter)
c
      else
c
c  for innermost points, do interpolation in (r/R)**2 or q**(2/3).
c  to get more accurate values for q and L
c  inside innermost non-zero meshpoint use linear interpolation
c
	write(6,*) 'special interpolation at n =',n
        nint=0
        if(xab2(n).gt.xab1(2)) then
          call lir(xabc2(n),xabc1,yint,yvarc1(1,1),nvarl,igvrmx,20,
     *      nint,inter)
        else
          call lir1(xabc2(n),xabc1,yint,yvarc1(1,1),nvarl,igvrmx,20,
     *      nint,inter)
        end if
c
c  reset variables, depending on icase
c
        if(n.eq.1) then
          yint( 1) = 0
          yint( 2) = 0
          yint(11) = 0
          yint(18) = 0
        else if(icase.eq.1.or.icase.eq.3.or.icase.eq.4) then
          yint( 2)=sqrt(yint(2))
          yint( 1)=1.e-33*yint( 1)*yint(2)**3
          yint(11)=       yint(11)*yint(2)**3
          yint(18)=1.e-22*yint(18)*yint(2)**2
        else
          yint(1)=yint(1)**1.5
          yint(2)=1.e11*yint(2)*yint(1)**0.3333333
          yint(11)=1.e33*yint(11)*yint(1)
          yint(18)=     yint(18)*xabc2(n)
        end if
      end if
c
c  test for interpolation errors
c
      if(inter.lt.0) go to 60
c
c  disallow extrapolation for icase = 2 or 5
c
      if(inter.eq.0.and.(icase.eq.2.or.icase.eq.5)) go to 60
c
      if(itest.eq.2) then
        write(6,130) 'yvar(i,n,2):',(yvar(i,n,2),i=1,nvar)
        write(6,130) 'yint:',yint
      end if
c
      if(yint(1).gt.0.and.yvar(1,n,2).gt.0) then
        dy(1)=-log(yint(1)/yvar(1,n,2))
      else
        dy(1)=0
      end if
      dy(2)=-amm*(yint(3)-yvar(3,n,2))
      dy(3)=-amm*(yint(4)-yvar(4,n,2))
      if(yint(11).gt.0.and.yvar(11,n,2).gt.0) then
        dy(4)=-log(yint(11)/yvar(11,n,2))
      else
        dy(4)=0
      end if
      dy(5)=     -yint(6)+yvar(6,n,2)
      dy(6)=-log(yint(9)/yvar(9,n,2))
      dy(7)=-log(yint(10)/yvar(10,n,2))
      if(yint(12).gt.0.and.yvar(12,n,2).gt.0) then
        dy(8)=-log(yint(12)/yvar(12,n,2))
      else
        dy(8)=0
      end if
      if(abs(yint(13)).gt.0.and.abs(yvar(13,n,2)).gt.1.e-5) then
        dy(9)=-log(abs(yint(13)/yvar(13,n,2))+1.e-20)
      else
        dy(9)=0
      end if
      if(yint(14).gt.0.and.yvar(14,n,2).gt.0) then
        dy(10)=-log(yint(14)/yvar(14,n,2))
      else
        dy(10)=0
      end if
      dy(11)=-log(yint(15)/yvar(15,n,2))
      dy(12)=-log(yint(16)/yvar(16,n,2))
      dy(13)=-log(yint(17)/yvar(17,n,2))
      if(vg2(n).ne.0) then
        dy(14)=(yvar(18,n,2)-yint(18))/vg2(n)
      else
        dy(14)=0
      end if
      dy(15)=0.5d0*(dy(6)+dy(10)-dy(7))
c
c  test for later set, including 19 and 20
c
      if(nvar.ge.20) then
        if(abs(yint(19)).gt.0.and.abs(yvar(19,n,2)).gt.1.e-30) then
          dy(16)=-log(yint(19)/yvar(19,n,2))
        else
          dy(16)=0
        end if
        if(abs(yint(20)).gt.0.and.abs(yvar(20,n,2)).gt.1.e-30) then
          dy(17)=-log(abs(yint(20)/yvar(20,n,2))+1.e-20)
        else
          dy(17)=0
        end if
      end if
      if(yint(2).gt.0.and.yvar(2,n,2).gt.0) then
        dy(18)=-log(yint(2)/yvar(2,n,2))
      else
        dy(18)=0
      end if
      if(abs(yint(21)).gt.0.and.abs(yvar(21,n,2)).gt.0) then
        dy(19)=-log(abs(yint(21)/yvar(21,n,2))+1.e-20)
      else
        dy(19)=0
      end if
      if(abs(yint(7)).gt.0.and.abs(yvar(7,n,2)).gt.0) then
        dy(20)=-log(abs(yint(7)/yvar(7,n,2))+1.e-20)
      else
        dy(20)=0
      end if
c
c  set extreme values
c
      if(n.eq.1) then
        do 35 i=1,ndiffa
        xmax(i)=xr2(1)
        dymax(i)=dy(i)
        dymin(i)=dy(i)
   35   daymax(i)=abs(dy(i))
c
      else
c
        do 37 i=1,ndiffa
        if(abs(dy(i)).gt.daymax(i)) then
          xmax(i)=xr2(n)
          daymax(i)=abs(dy(i))
        end if
        if(dy(i).gt.dymax(i)) then
          dymax(i)=dy(i)
        else if(dy(i).lt.dymin(i)) then
          dymin(i)=dy(i)
        end if
   37   continue
c
      end if
c
      write(10,160) xr2(n),xq2(n),(dy(i),i=1,ndiffa)
c
c  test for Guenther et al output
c
      if(igjdpr.eq.1) then
        rr=yint(2)
        drho=yint(10)*dy(7)
        dc=sqrt(yint(9)*yint(14)/yint(10))*dy(15)
        dgamma1=yint(14)*dy(10)
        dkappa=yint(12)*dy(8)
        write(14,165) n, xr2(n), rr, drho, dc, dgamma1, dkappa
      end if
   60 continue
c
c  output maximum values
c
      write(6,170) (i,dyname(i),xmax(i),dymin(i),dymax(i),i=1,ndiffa)
      write(12,170) (i,dyname(i),xmax(i),dymin(i),dymax(i),i=1,ndiffa)
c
   90 continue
c
      stop
  100 format(//' differences between GONG models.'/
     *  ' First  model: no.',i3,' on file ',a60/
     *  ' Second model: no.',i3,' on file ',a60)
  102 format('#'/'#'/'# differences between GONG models.'/
     *  '# First  model: no.',i3,' on file ',a60/
     *  '# Second model: no.',i3,' on file ',a60)
  103 format(//' summary written to file dgong-sum')
  104 format(//' differences at fixed ',a)
  105 format('#'/'#'/'# differences at fixed ',a)
  106 format(//' differences at fixed m/M(r/R =',f10.6,')')
  107 format('#'/'#'/'# differences at fixed m/M(r/R =',f10.6,')')
  108 format('#'/'# n, r/R, r, delta rho, delta c, delta Gamma_1,',
     *  ' delta kappa'/'#')
  110 format(//' ***** Error in diff-gong for icase = 7:'/
     *         '       rmsfix =',f10.6,' .gt. Rs/R =',f10.6)
  113 format('#'/'#  model no. ',i2/'#'/
     *  '# ',a6,' precision data')
  115 format('#'/'#'/'# cdata:'/('#',a79))
  116 format('#'/'#'/'# datmod:'/('#',1p5e15.6))
  117 format('#'/'#'/'# idatmd:'/('#',5i15))
  118 format('#'/'#'/'# datgng:'/('#',1p5e15.6))
  119 format('#'/'#'/'# bc:'/('#',1p5e15.6))
  122 format(//' n, xabc1(n), varc1(1,n), varc1(11,n), varc1(18,n):'/
     *  (i4,1p4e14.6))
  123 format(//' n, xabc1(n), varc1(2,n), varc1(11,n), varc1(18,n):'/
     *  (i4,1p4e14.6))
  130 format(/a10/(1p10e12.4))
  140 format(//' nn =',i6,'   nvar =',i4//)
  150 format('#'/'#'/'# output format:'/
     *  '# r/R(2), log q(2), dy(1-',i2,'),       where'/'#'/
     *  ('# dy(',i2,'): ',a20))
  155 format('#')
  160 format(0pf12.6,1pe13.5,20e11.3)
  165 format(i5,f12.7,1pe13.5,4e12.4)
  170 format(//' extreme differences.'/
     *  '   variable                  xextreme          extremes:'/
     *  '                                            min        max'//
     *  (i4,1x,a20,0pf12.6,1p2e13.5))
      end
