      subroutine enghec(fl,t,x,z,epshec,fthec,ifthec, secder)
c
c  Sets contributions to energy generation rate from triple-alpha
c  and 12C burning, as well as rates of change of 4He and 12C abundances.
c
c  On input, fl is log(f) (or whatever),
c  t is temperature, and x(.) is an array containing abundances;
c  for now, 4He abundance is assumed to be in x(nspec+1) and 12C
c  abundance in x(nspec+2), where nspec is the number of elements
c  involved in H burning, set by s/r engcse.
c  [NEEDS later change for consistency with rest of programme??]
c  z is heavy element abundances.
c
c  Subroutine controlled by iheccs in common /engcnt/. 4He burning
c  is included in the energy generation if iheccs .ne. 0.
c  If iheccs .lt. 0 12C burning is switched off.
c
c  Returns epshec(1) as energy generation from 4He and 12C burning
c  epshec(2, ...) as first derivatives of log(epshec).
c  Finally, fthec(1,k), k = 1,2 are rates of change of 4He and 12C,
c  respectively,
c  and fthec(2 - , k) are first derivatives of fthec(1,k),
c  in the usual manner, including all nspec + 2 abundances.
c
c  Total and partial 4He and 12C energy generation rates are set in
c  common /cephec/, in the following variables (in ergs/sec):
c
c  epsche:    total rate from 4He and 12C
c  eppche(1): Rate associated with reaction 3 4He -> 12C
c  eppche(2): Rate associated with reaction 12C + 4He -> 16O
c
c  NOTE: for now only first derivatives are set. However, argument
c  secder (to flag for setting second derivatives, possibly by numerical
c  differentiation) is kept for consistency.
c
c  Original date: 26/7/02
c
      implicit double precision (a-h,o-z)
      logical secder, norche, nosd
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
      parameter (idermx = ((nspcmx+3)*(nspcmx+4))/2)
c
      dimension x(*),epshec(*),fthec(ifthec,*)
      dimension qhec(2), echtl1(idermx), echtl2(idermx)
c
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
      common/consts/ av,ah,ahe
      common/eqstd/ xii(4),ane(10),rho(10)
      common/rnrhed/ alhe(10,10),norche
      common/rcncns/ a(knucst),b(knucst),s(knucst,isnuc),
     *  zz(knucst),zz12(2,knucst),zsmh,acno,awght(knucwg),q(knucq), 
     *  knucpm, kqvals, kawght
      common/engcnt/ xhzlm1, xhzlm2, nxhzer, icnocs, nspec, 
     *  iheccs, nspect
      common/he3eql/ xhe3eq,ieqhe3,iequsd
      common/cengcs/ ixc12, ixc13, ixn14, ixn15, ixo16, ixo17
      common/rnrcnt/ irnrat, krnrat, irn12, irn13, irn14, 
     *  irn15c, irn15o, irn16, irn17
      common/compos/ xseth, yset, xset3, xset12, xset13, xset14, xset16
      common/ln10/ amm,amm2
      common/cnofrc/ fcno, xtlcno
      common/cephec/ epsche, eppche(2)
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  Q values (in MeV) from Angulo et al.
c
      data qhec /7.274d0, 7.162d0/
c
c  statement function defining storage of second derivatives
c
      jjsder(i, j) = (i*(2*nspect+5-i))/2 +j +1
c
c..      write(6,'(a,0pf12.7,1p2e15.7,0p2f10.7)') 
c..     *  'Enter enghec with fl, t, rho, x =',fl, t, rho(1), x(1),x(2)
c
c  if iheccs = 0, return zero energy generation
c
      if(iheccs.eq.0) then
	norche=.true.
	eppche(1)=0
	eppche(2)=0
	epsche=0
	do 10 i=1,nspec+5
	epshec(i)=0
	do 10 k=1,2
   10   fthec(i,k)=0
c
	return
c
      end if
c
      xh = x(1)
c
      yh = x(nspec+1)
      xset12 = x(nspec+2)
c
c
      call zero(eppche,2)
c
      if(istdpr.gt.0.and.(idgeng.lt.0.or.idgeng.ge.1)) then
        write(istdpr,*) 'Entering enghec. rho, T, Xi', rho(1), t,
     *    (x(i),i=1,4)
      end if
c
c  set reaction rates
c
      tl=log10(t)
      nosd=.true.
      call rnrhec(fl,tl,xh,yh,z,nosd)
c
c  if no reactions, return zero energy generation
c
      if(norche) then
	eppche(1)=0
	eppche(2)=0
	epsche=0
	do 15 i=1,nspec+5
	epshec(i)=0
	do 15 k=1,2
   15   fthec(i,k)=0
c
	return
c
      end if
c
c  if iheccs .lt. 0 set 12C + 4He rate to zero
c
      if(iheccs.lt.0) alhe(2,1)=0
c
      yfact=x(nspec+1)/ahe
      yfact3=yfact*yfact*yfact
      cfact=x(nspec+2)/awght(2)
c..      write(6,'(a,1p5e13.5)') 'yfact, etc:',yfact,cfact,av,alhe(1,1),
c..     *  alhe(2,1)
c
c  set energy generation rate, resetting from MeV to erg
c
      echtl1(1)=1.d6*ergev*qhec(1)*rho(1)*rho(1)*yfact3*av*alhe(1,1)
      echtl2(1)=1.d6*ergev*qhec(2)*rho(1)*cfact*yfact*av*alhe(2,1)
      eppche(1)=echtl1(1)
      eppche(2)=echtl2(1)
      epsche=eppche(1)+eppche(2)
      epshec(1)=epsche
c..      write(6,*) 'eppche(1-2):',eppche
c
c  set logarithmic first derivatives of echtl1(1) and echtl2(1)
c  wrt (log f, log T, X_i)
c
      call zero(echtl1(2),idermx-1)
      call zero(echtl2(2),idermx-1)
      do 20 i=2,4
      echtl1(i)=2*rho(i)+alhe(1,i)
   20 echtl2(i)=rho(i)+alhe(2,i)
c
      echtl1(nspec+4)=3.d0/(amm*x(nspec+1))
      echtl2(nspec+4)=1.d0/(amm*x(nspec+1))
      echtl2(nspec+5)=1.d0/(amm*x(nspec+2))
c
c  set logarithmic first derivatives of epshec
c  wrt (log f, log T, X_i)
c
      do 25 i=2,nspec+5
   25 epshec(i)=(echtl1(1)*echtl1(i)+echtl2(1)*echtl2(i))/epshec(1)
c
c  set rates of change of 4He and 12C abundances
c
      ftct11=-3.d0*rho(1)*rho(1)*yfact3*alhe(1,1)
      ftct12=-rho(1)*cfact*yfact*alhe(2,1)
      fthec(1,1)=ahe*(ftct11+ftct12)
      ftct21=rho(1)*rho(1)*yfact3*alhe(1,1)
      ftct22=-rho(1)*cfact*yfact*alhe(2,1)
      fthec(1,2)=awght(2)*(ftct21+ftct22)
c
c  set first derivatives of ft(1,1) and ft(1,2) 
c  wrt (log f, log T, X_i)
c
      call zero(fthec(2,1),ifthec-1)
      call zero(fthec(2,2),ifthec-1)
      do 30 i=2,4
      fthec(i,1)=amm*ahe*(ftct11*(2.d0*rho(i)+alhe(1,i))+
     *                    ftct12*(rho(i)+alhe(2,i)))
   30 fthec(i,2)=amm*awght(2)*(ftct21*(2.d0*rho(i)+alhe(1,i))+
     *                         ftct22*(rho(i)+alhe(2,i)))
c
      if(x(nspec+1).gt.0) then
        fthec(nspec+4,1)=ahe*(3.d0*ftct11+ftct12)/x(nspec+1)
        fthec(nspec+4,2)=awght(2)*(3.d0*ftct21+ftct22)/x(nspec+1)
      end if
      if(x(nspec+2).gt.0) then
        fthec(nspec+5,1)=ahe*ftct12/x(nspec+2)
        fthec(nspec+5,2)=awght(2)*ftct22/x(nspec+2)
      end if
      if(idgeng.eq.-99.and.istdpr.gt.0) write(istdpr,*) 'in enghec:',
     *  nspec,awght(2),ftct22,x(nspec+2),fthec(nspec+5,2)
c
      return
      end
