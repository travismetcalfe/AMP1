      program main
c
c  set formatted file of GONG model variables, as defined in
c  GONG model manual, on the basis of GONG model set.
c  note that variables is stored in same order as evolution
c  variables, i.e. with surface at point 1.
c
c  format:
c
c  write(idsgng) (cdata(i),i=1,4),nmod,nn,nvar,(datmod(i),i=1,31),
c    (datgng(i),i=1,30),(bc(i),i=1,54),
c    ((yvar(i,n),i=1,nvar),n=1,nn)
c
c  original version: 20/03/88
c
c  Modified 17/8/92, to take into account extended format of
c  binary output
c
c  Modified 29/11/93, to include epsilon_g as variable 19
c
c  Modified 5/5/95, including He3 and CNO abundances, depending
c  on parameter. If so, reset ivers to 200
c
c  Modified 12/12/95, including also Z. Reset version number to 210
c
c  Modified 17/12/95, giving possibility of including derivatives of
c      Gamma_1. If so, set version number to 250.
c  
c  Double precision version
c  ++++++++++++++++++++++++
c
      include 'evolpr.incl'
      parameter(ivrout=30)
      implicit double precision (a-h, o-z)
      character*80 cdata,fin, fout, fhead, head, fgam1
      dimension datmod(nrdtmx),ndtmod(nidtmx),yvar(igvrmx,nnmax),
     *  cdata(4),datgng(ndtgmx),bc(nbccmx),
     *	 dgam1(13,nnmax),datout(15), varout(ivrout), head(4)
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .	ct,crho,cpe,che,ca03,caa,ckh,car,ergev1,exhm
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm,
     *	cgrav, amsun, echar, ergev, syear, iver
      common/ln10/ amm
c
      data head /4*' '/
c
      call setcns
c
      write(6,*) 'Enter input file name'
      read(5,'(a)') fin
      open(2,file=fin,status='old',form='unformatted')
c
      write(6,*) 'Enter output file name'
      read(5,'(a)') fout
      open(3,file=fout,status='unknown')
c
      write(6,*) 'Enter number of model on file'
      read(5,*) nrdout
      write(6,100) fin
c
      write(6,*) 'Enter file name for header'
      read(5,'(a)') fhead
      open(4,file=fhead,status='old')
c
      read(4,'((a))',end=5) (head(i),i=1,4)
    5 continue
c
      write(6,*) 'Enter case number (determining output)'
      write(6,*) '   case 0: old (pre-1995) format'
      write(6,*) '   case 1: including He3 and CNO abundances'
      write(6,*) '   case 2: including He3 and CNO abundances and Z'
      write(6,*) '   case 3: including He3 and CNO abundances and Z,'
      write(6,*) '           as well as Gamma_1 derivatives'
      read(5,*) iout
c
      nrd=0
      icase=0
c
   10 call rdgong(2,cdata,nmod,datmod,ndtmod,datgng,bc,yvar,
     *  igvrmx,iform,nn,nrdtmd,nidtmd,ndtgng,nbccf,nvar,icase,icry)
      icase=-1
      if(icry.eq.-1) then
	go to 20
      else if(icry.eq.-2) then
	go to 90
      end if
c
      nrd=nrd+1
c
      if(nrd.lt.nrdout) go to 10
c
   20 write(6,102) nrd
c
      write(6,105) cdata
      call wdtmod(6,datmod,nrdtmd,ndtmod,nidtmd,0)
      call wdtgng(6,datgng,ndtgng,0)
      call wrtbc(6,bc,nbccf,0)
      write(6,140) nn,nvar
c
c  set version number etc.
c
      iconst=15
      if(iout.eq.0) then
        ivar=20
        ivers=100
      else if(iout.eq.1) then
        ivar=25
	icnocs=ndtmod(5)
        ivers=200
      else if(iout.eq.2) then
        ivar=25
	icnocs=ndtmod(5)
        ivers=210
      else if(iout.eq.3) then
        ivar=30
	icnocs=ndtmod(5)
        ivers=250
      else
	write(6,145) iout
	stop
      end if
c
c  test for reading gamma_1 coefficients
c
      if(iout.ge.3) then
	write(6,*) 'Enter file of Gamma_1 coefficients'
	read(5,'(a)') fgam1
	open(10,file=fgam1,status='old')
	call skpcom(10)
	n=1
   22   read(10,*,end=23) (dgam1(i,n),i=1,13)
	n=n+1
	go to 22
   23   continue
	nngam=n-1
	if(nngam.ne.nn) then
	  write(6,150) nngam, nn
	  stop
        end if
      end if
c
c  start output to file
c
      do 25 i=1,iconst
   25 datout(i)=0
      do 27 i=1,ivar
   27 varout(i)=0
c
c  set constants
c
      datout(1) = datmod(23)
      datout(2) = datmod(24)
      datout(3) = datmod(25)
      datout(4) = datmod(1)
      datout(5) = datmod(2)
      datout(6) = datmod(3)
c
      datout(7) = 1.0d0/(2*datmod(4)*sqrt(datmod(5)))
      datout(8) = datmod(4)**2/4
c
      datout(9) = datmod(14)
      if(ndtmod(3).le.1) then
	datout(10) = datmod(15)
      else
	datout(10) = 1
      end if
c
      rsq = datmod(24)*datmod(24)
c
      datout(11) = rsq*datgng(13)/datgng( 6)
      datout(12) = rsq*datgng(14)/datgng(10)
c
      write(3,160) (head(i),i=1,4)
      write(3,165) nn, iconst, ivar, ivers
      write(3,170) (datout(i),i=1, iconst)
c
c   step through model for variable output
c
      do 40 n=1,nn
      varout(1) = 1.d11*10.d0**yvar(2,n)
      varout(2) = amm*yvar(1,n)
      varout(3) = yvar(8,n)
      varout(4) = yvar(9,n)
      varout(5) = yvar(10,n)
      varout(6) = yvar(6,n)
      varout(7) = yvar(11,n)
      varout(8) = yvar(12,n)
      varout(9) = yvar(13,n)
      varout(10) = yvar(14,n)
      varout(11) = yvar(15,n)
      varout(12) = yvar(16,n)
      varout(13) = yvar(17,n)
      varout(14) = amu*yvar(19,n)
      varout(15) = yvar(18,n)
      varout(16) = yvar(20,n)
      if(iout.lt.2) then
        varout(17) = datout(2)-varout(1)
      else
        varout(17) = 1-yvar(21,n)-varout(6)
        varout(18) = datout(2)-varout(1)
      end if
      varout(19) = yvar(29,n)
c
      if(iout.gt.0) then
	varout(21) = yvar(7,n)
	if(icnocs.eq.0) then
	  varout(24) = datmod(52)*datmod(1)
	else if(icnocs.eq.1) then
	  varout(24) = yvar(22,n)
	  varout(25) = yvar(23,n)
	else if(icnocs.eq.2.or.icnocs.eq.4) then
	  varout(22) = yvar(22,n)
	  varout(23) = yvar(23,n)
	  varout(24) = yvar(24,n)
	  varout(25) = yvar(25,n)
	else
	  varout(22) = yvar(22,n)
	  varout(23) = yvar(23,n)
	  varout(24) = yvar(24,n)
        end if
      end if
c
      if(iout.ge.3) then
	varout(26)=dgam1(9,n)
	varout(27)=dgam1(8,n)
	varout(28)=dgam1(10,n)
      end if
c
   40 write(3,170) (varout(i),i=1,ivar)
c
   90 continue
c
      stop
  100 format(//' Set formatted GONG models from file ',a60)
  102 format(///' Formatted output of record no',i5)
  105 format(//' cdata:'/(a))
  140 format(//' nn =',i6,'   nvar =',i4//)
  145 format(//' ***** iout = ',i4,' not implemented')
  150 format(//' ***** Error. nngam = ',i5,' not equal to nn =',i5)
  160 format(a)
  165 format(4i10)
  170 format(1p5e16.9)
      end
