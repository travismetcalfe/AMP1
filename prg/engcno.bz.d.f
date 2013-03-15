      subroutine engcno(t,x,z,alt,epscno,fthcno,ft,idalt,ift,
     *  ider1, ider2, ider, iderh, nspech, nspcmh, secder)
c
c  Sets contributions to energy generation rate and rate of change 
c  of X from CNO cycle, as well as rates of change of abundances
c  in CNO cycle.
c
c  On input, t is temperature, x(.) is an array containing abundances
c  z is heavy element abundances, and alt must have been set
c  to reaction rates and derivatives of lambda.
c  In particular, nspech (= 1 or 2) is the number of elements amongst
c  H, 3He that are included.
c  (Note that the total number of elements, nspect, is passed in 
c  common/engcnt/ and is assumed to have been set by call of s/r engcse.
c  Also, storage indices ider1, ider2, ider, iderh, nspech and nspcmh 
c  must have been set (this is done automatically when engcno is
c  called from s/r engenr)).
c
c  Returns epscno(1) as energy generation from CNO cycle and
c  epscno(2, ...) as derivatives of log(epscno).
c  Also fthcno(1) is contribution to rate of change of X,
c  and fthcno(2, ...) are derivatives of fthcno(1).
c  Finally, ft(1,k + nspech) is rate of change of k-th CNO element,
c  and ft(2 - , k + nspech) are derivatives of ft(1,k + nspech).
c
c  Total and partial CNO energy generation rates are set in
c  common /cepcno/, in the following variables (in ergs/sec):
c
c  epscnc:    total CNO rate
c  eppcno(1): Rate associated with reaction 12C + 1H
c  eppcno(2): Rate associated with reaction 13C + 1H
c  eppcno(3): Rate associated with reaction 14N + 1H
c  eppcno(4): Rate associated with reaction 16O + 1H
c
c  Original date: 26/7/91.
c
c  Modified 15/7/96, to include the full case (icnocs = 4), albeit 
c  still without second derivatives
c
c  Modified 24/7/96, to include output of total and partial CNO
c  energy generation rates in commin /cepcno/
c
c  Corrected 5/8/97, by including engenr.nnz.d.incl, to set nspdmx etc
c
c  Modified 1/8/02, changing storage to use nspect instead of nspec,
c  to account for inclusion of helium burning
c
c  Modified 16/6/04, allowing for case xtlcno = 0. In this case, 
c  simply return with zero energy generation and rate of change of
c  CNO abundances
c
      implicit double precision (a-h,o-z)
      logical secder
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.nnz.d.incl'
      parameter (idermx = ((nspcmx+3)*(nspcmx+4))/2)
c
      dimension x(*),alt(idalt,*),epscno(*),fthcno(*),ft(ift,*)
      dimension dfth(4), gam15(10), al15(10), qn14(10), epstil(idermx),
     *  fthtil(idermx), fcnotl(idermx, nspcmx)
c
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
      common/noiter/ iter, ntime
      common/consts/ av,ah,ahe
      common/eqstd/ xii(4),ane(10),rho(10)
      common/rnratd/ al(10,krnrmx),norct
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
      common/cepcno/ epscnc, eppcno(4)
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
c  statement function defining storage of second derivatives
c
      jjsder(i, j) = (i*(2*nspect+5-i))/2 +j +1
c
      xh = x(1)
c
      call zero(eppcno,4)
      epscnc=0
c
c  for xtlcno (assumed to flag no CNO elements) simply return
c
      if(xtlcno.eq.0) return
c
      if(icnocs.eq.0) then
c
c  equilibrium CNO energy generation rate
c
        epscno(1)=av*q(6)*rho(1)*xh*z*fcno*al(1,5)/(ah*acno)
c..	if(log10(t).gt.7.3.and.ntime.le.3) 
c..     *    write(6,*) 'engcno', rho(1), t, q(6), xh, z, fcno, epscno(1)
        do 15 i=2,4
   15   epscno(i)=rho(i)+al(i,5)
        epscno(4)=epscno(4)+1/(amm*xh)
        if(iderh.eq.5) epscno(5)=0
c
c  reset from MeV to erg
c
        epscno(1)=1.d6*ergev*epscno(1)
c
c  contribution to rate of change of X
c
        ct4=4*fcno*z/acno
        cc4=ct4*xh
        fthcno(1)=-cc4*rho(1)*alt(1,5)
        do 20 i=2,4
        dfth(i)=-cc4*alt(i,5)
        if(i.eq.4) dfth(i) = dfth(i) - ct4*alt(1,5)
   20   fthcno(i)=amm*rho(i)*fthcno(1)+rho(1)*dfth(i)
c
        if(secder) then
c
c  setting second derivatives
c
          nspcm1=nspect-1
          jjxx =jjsder(3, 3)
          ii=4
          jj=ider1
          do 30 i=2,4
          do 25 j=i,4
          ii=ii+1
          jj=jj+1
          epscno(jj) = rho(ii) + al(ii,5)
          ddfth=0
          if(i.eq.4) ddfth=ddfth-alt(j,5)
          if(j.eq.4) ddfth=ddfth-alt(i,5)
   25     fthcno(jj) = amm*(rho(ii)*fthcno(1)+rho(i)*fthcno(j)+
     *      rho(j)*rho(1)*dfth(i))+rho(1)*(-cc4*alt(ii,5)+ct4*ddfth)
   30     jj=jj+nspcm1
c
          epscno(jjxx)=epscno(jjxx)-1/(amm*xh*xh)
        end if
c
c  in this case, set just the assumed N14 abundance in common/compos/
c
        xset14 = fcno*z/acno
c
	epscnc=epscno(1)
c
        return
c
      else if(icnocs.ne.1.and.icnocs.ne.4) then
c
c  diagnostic message 
c
        write(istdou,190) icnocs
        if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,190) icnocs
        stop 'engcno'
      end if
c
c  start section handling general CNO cases
c
      if((idgeng.eq.-1.or.idgeng.ge.1).and.istdpr.gt.0) then
        write(istdpr,*) 'Entering engcno. rho, T, Xi', rho(1), t,
     *    (x(i),i=1,4)
        write(istdpr,*) 'ieqhe3,icnocs, irn14,irn15o,irn15c, ixn14',
     *    ieqhe3,icnocs, irn14,irn15o,irn15c, ixn14
       end if
c
c  test for setting alpha15 and gamma15.
c
      if(icnocs.ne.3) then
        al15(1)=al(1,irn15o)+al(1,irn15c)
        if(al15(1).gt.0) then
          gam15(1)=al(1,irn15o)/al15(1)
          alph15 = al(1,irn15c)/al15(1)
          do 35 i=2,4
          al15(i) = alph15*al(i,irn15c)+gam15(1)*al(i,irn15o)
   35     gam15(i) = amm*gam15(1)*(al(i,irn15o) - al15(i))
c
        else
c
          alph15=1
          do 36 i=1,10
   36     gam15(i)=0
c
        end if
c
        if(secder.and.al15(1).gt.0) then
          ii = 4
          do 37 i=2,4
          do 37 j=i,4
          ii=ii+1
          al15(ii)=(al(i,irn15o) - al(i,irn15c))*gam15(j)
     *             +alph15*al(ii,irn15c)+gam15(1)*al(ii,irn15o)
   37     gam15(ii)=amm*(gam15(i)*(al(j,irn15o) - al15(j))
     *                  +gam15(1)*(al(ii,irn15o) - al15(ii)))
        end if
c
      end if
c
      if(idgeng.ne.0.and.istdpr.gt.0) 
     *  write(istdpr,*) 'ider, nspect, idermx',
     *  ider,nspect,idermx
c
c  zero all intermediate derivatives
c
      call zero(fthtil,ider)
      call zero(epstil,ider)
      call zero(fcnotl,nspect*idermx)
c
c  set brackets and derivatives, depending on icnocs
c
      if(icnocs.eq.1) then
        xtln14 =  x(ixn14)/awght(4)
	if(xtlcno.ge.1.e-12) then
          xtlo16 = xtlcno - xtln14
        else
          xtlo16 = 0.
        end if
c
        if(idgeng.ge.1.and.istdpr.gt.0) 
     *    write(istdpr,*) 'xtln14, xtlo16', xtln14, xtlo16
        cft14 = awght(4)/ah
        fcnotl(1, ixn14) = cft14*(-gam15(1)*al(1,irn14)*xtln14
     *    + al(1,irn16)*xtlo16)
	if(idgeng.eq.-1.and.istdpr.gt.0) then
	  write(istdpr,*)
     *     'cft14, gam15(1), al(1,irn14), xtln14, al(1,irn16), xtlo16,',
     *       'fcnotl(1, ixn14)'
          write(istdpr,*)
     *      cft14, gam15(1), al(1,irn14), xtln14, al(1,irn16), xtlo16,
     *      fcnotl(1, ixn14) 
	end if
        fthtil(1) = -2*((2 + gam15(1))*al(1,irn14)*xtln14
     *    + al(1,irn16)*xtlo16)
        cepstl = av/ah
        qn14(1) = alph15*q(6) + gam15(1)*(q(9) + q(11))
        epstil(1) = cepstl*(qn14(1)*al(1,irn14)*xtln14
     *    + (q(12) + q(13))*al(1,irn16)*xtlo16)
        eppcno(3) = cepstl*qn14(1)*al(1,irn14)*xtln14
        eppcno(4) = cepstl*(q(12) + q(13))*al(1,irn16)*xtlo16
        if(idgeng.eq.-10.and.istdpr.gt.0) write(istdpr,*) 
     *    'epstil(1),cepstl,qn14(1),al(1,irn14),xtln14',
     *    'q(12),q(13),al(1,irn16),xtlo16',
     *    epstil(1),cepstl,qn14(1),al(1,irn14),xtln14,
     *    q(12),q(13),al(1,irn16),xtlo16
c
c  first derivatives
c
        do 40 i = 2,4
        fcnotl(i, ixn14) = cft14*(
     *    -(gam15(i)*al(1,irn14)+gam15(1)*alt(i,irn14))*xtln14
     *    + alt(i,irn16)*xtlo16)
        fthtil(i) = -2*(((2 + gam15(1))*alt(i,irn14) + 
     *    gam15(i)*al(1,irn14))*xtln14 + alt(i,irn16)*xtlo16)
        qn14(i) = gam15(i)*(q(9) + q(11) - q(6))
   40   epstil(i) = cepstl*((qn14(i)*al(1,irn14) + 
     *    qn14(1)*alt(i,irn14))*xtln14
     *    + (q(12) + q(13))*alt(i,irn16)*xtlo16)
c
c  derivatives wrt N14 abundance
c
        idn14 = ixn14 + 3
        fcnotl(idn14, ixn14) = 
     *    -cft14*(gam15(1)*al(1,irn14) + al(1,irn16))/awght(4)
        fthtil(idn14) = 
     *    -2*((2 + gam15(1))*al(1,irn14) - al(1,irn16))/awght(4)
        epstil(idn14) = cepstl*(qn14(1)*al(1,irn14) - 
     *    (q(12) + q(13))*al(1,irn16))/awght(4)
c
c  test for setting second derivatives
c
        if(secder) then
          nspcm1 = nspect - 1
          ii=4
          jj=ider1
          do 44 i = 2,4
          do 42 j = i,4
          ii = ii+1
          jj = jj+1
          fcnotl(jj,ixn14) = cft14*(
     *      -(gam15(ii)*al(1,irn14) + gam15(i)*alt(j,irn14) + 
     *      gam15(j)*alt(i,irn14)+gam15(1)*alt(ii,irn14))*xtln14
     *      + alt(ii,irn16)*xtlo16)
          fthtil(jj) = -2*(((2 + gam15(1))*alt(ii,irn14) + 
     *      gam15(i)*alt(j,irn14) + gam15(j)*alt(i,irn14) + 
     *      gam15(ii)*al(1,irn14))*xtln14 + alt(ii,irn16)*xtlo16)
          qn14(ii) = gam15(ii)*(q(9) + q(11) - q(6))
   42     epstil(jj) = cepstl*((qn14(ii)*al(1,irn14) + 
     *      qn14(i)*alt(j,irn14) + qn14(j)*alt(i,irn14) +
     *      qn14(1)*alt(ii,irn14))*xtln14
     *      + (q(12) + q(13))*alt(ii,irn16)*xtlo16)
c
c  derivatives involving N14 abundance
c
          jjdn14 = jjsder(i-1, ixn14 + 2)
          fcnotl(jjdn14, ixn14) = 
     *      -cft14*(gam15(1)*alt(i,irn14) + gam15(i)*al(1,irn14) + 
     *      alt(i,irn16))/awght(4)
          fthtil(jjdn14) = 
     *      -2*((2 + gam15(1))*alt(i,irn14) + gam15(i)*al(1,irn14) - 
     *      alt(i,irn16))/awght(4)
          epstil(jjdn14) = cepstl*(qn14(i)*al(1,irn14) + 
     *      qn14(1)*alt(i,irn14) - 
     *      (q(12) + q(13))*alt(i,irn16))/awght(4)
   44     jj = jj + nspcm1
        end if
c
c..     write(6,910) 'fcnotl(.,ixn14)',(fcnotl(k,ixn14),k=1,ider)
c..     write(6,910) 'epstil',(epstil(k),k=1,ider)
c..     write(6,910) 'fthtil',(fthtil(k),k=1,ider)
c
c  set N14 and O16 abundances in common/compos/
c
        xset14 = x(ixn14)
        xset16 = xtlo16*awght(6)
c
c  end for icnocs = 1
c
      else if(icnocs.eq.4) then
c..	write(6,*) 'ixc12, x(ixc12), awght(2)',ixc12, x(ixc12), awght(2)
        xtlc12 =  x(ixc12)/awght(2)
        xtlc13 =  x(ixc13)/awght(3)
        xtln14 =  x(ixn14)/awght(4)
	if(xtlcno.ge.1.e-12) then
          xtlo16 = xtlcno - xtlc12 - xtlc13- xtln14
        else
          xtlo16 = 0.
        end if
c
        if(idgeng.ge.1.and.istdpr.gt.0) 
     *    write(istdpr,*) 'xtln14, xtlo16', xtln14, xtlo16
        cft12 = awght(2)/ah
        cft13 = awght(3)/ah
        cft14 = awght(4)/ah
        fcnotl(1, ixc12) = cft12*(-al(1,irn12)*xtlc12
     *    + alph15*al(1,irn14)*xtln14)
        fcnotl(1, ixc13) = cft13*(al(1,irn12)*xtlc12
     *    - al(1,irn13)*xtlc13)
        fcnotl(1, ixn14) = cft14*(al(1,irn13)*xtlc13-al(1,irn14)*xtln14
     *    + al(1,irn16)*xtlo16)
	if(idgeng.eq.-1.and.istdpr.gt.0) then
          write(istdpr,*)
     *     'cft14, gam15(1), al(1,irn14), xtln14, al(1,irn16), xtlo16,',
     *       'fcnotl(1, ixn14)'
          write(istdpr,*)
     *     cft14, gam15(1), al(1,irn14), xtln14, al(1,irn16), xtlo16,
     *       fcnotl(1, ixn14) 
	end if
        fthtil(1) = -al(1,irn12)*xtlc12 -al(1,irn13)*xtlc13 
     *     -2*al(1,irn14)*xtln14 - 2*al(1,irn16)*xtlo16
        cepstl = av/ah
        qn14(1) = q(9) + alph15*q(10) + gam15(1)*q(11)
        epstil(1) = cepstl*(q(7)*al(1,irn12)*xtlc12
     *    + q(8)*al(1,irn13)*xtlc13 + qn14(1)*al(1,irn14)*xtln14
     *    + (q(12) + q(13))*al(1,irn16)*xtlo16)
        eppcno(1) = cepstl*q(7)*al(1,irn12)*xtlc12
        eppcno(2) = cepstl*q(8)*al(1,irn13)*xtlc13 
        eppcno(3) = cepstl*qn14(1)*al(1,irn14)*xtln14
        eppcno(4) = cepstl*(q(12) + q(13))*al(1,irn16)*xtlo16
c
c  first derivatives
c
c..	if(t.ge.1.e7) then
c..	  write(6,*) 'In engcno, T, X, Z =',t, x, z
c..	  write(6,*) 'cft12, xtlc12, ixc12, fcnotl(1, ixc12)',
c..     *                cft12, xtlc12, ixc12, fcnotl(1, ixc12)
c..	  write(6,*) 'al(1,irn14), xtln14', al(1,irn14), xtln14
c..	  write(6,*) 
c..     *      'i, alt(i,irn12), gam15(i), alt(i,irn14), fcnotl(i,ixc12)'
c..	end if
        do 65 i = 2,4
        fcnotl(i, ixc12) = cft12*(-alt(i,irn12)*xtlc12
     *    +(-gam15(i)*al(1,irn14)+alph15*alt(i,irn14))*xtln14)
c..	if(t.gt.1.e7) 
c..     *    write(6,*) 
c..     *      i, alt(i,irn12), gam15(i), alt(i,irn14), fcnotl(i,ixc12)
        fcnotl(i, ixc13) = cft13*(alt(i,irn12)*xtlc12
     *    -alt(i,irn13)*xtlc13)
        fcnotl(i, ixn14) = cft14*(alt(i,irn13)*xtlc13
     *    -alt(i,irn14)*xtln14 + alt(i,irn16)*xtlo16)
        fthtil(i) = -alt(i,irn12)*xtlc12 - alt(i,irn13)*xtlc13
     *     -2*alt(i,irn14)*xtln14 -2*alt(i,irn16)*xtlo16
        qn14(i) = gam15(i)*(q(11) - q(10))
   65   epstil(i) = cepstl*(q(7)*alt(i,irn12)*xtlc12
     *    + q(8)*alt(i,irn13)*xtlc13 
     *    + (qn14(i)*al(1,irn14) + qn14(1)*alt(i,irn14))*xtln14
     *    + (q(12) + q(13))*alt(i,irn16)*xtlo16)
c
c..      if(t.gt.1.e7) stop
c
c  derivatives wrt CNO abundances
c
        idc12 = ixc12 + 3
        idc13 = ixc13 + 3
        idn14 = ixn14 + 3
	fcnotl(idc12, ixc12)= -cft12*al(1,irn12)/awght(2)
	fcnotl(idn14, ixc12)=  cft12*alph15*al(1,irn14)/awght(4)
	fcnotl(idc12, ixc13)=  cft13*al(1,irn12)/awght(2)
	fcnotl(idc13, ixc13)= -cft13*al(1,irn13)/awght(3)
	fcnotl(idc12, ixn14)= -cft14*al(1,irn16)/awght(2)
	fcnotl(idc13, ixn14)=  
     *    cft14*(al(1,irn13) - al(1,irn16))/awght(3)
	fcnotl(idn14, ixn14)= 
     *   -cft14*(al(1,irn14) + al(1,irn16))/awght(4)
        fthtil(idc12) = -(al(1,irn12)-2*al(1,irn16))/awght(2)
        fthtil(idc13) = -(al(1,irn13)-2*al(1,irn16))/awght(3)
        fthtil(idn14) = -2*(al(1,irn14) - al(1,irn16))/awght(4)
        epstil(idc12) = cepstl*
     *    (q(7)*al(1,irn12) - (q(12) + q(13))*al(1,irn16))/awght(2)
        epstil(idc13) = cepstl*
     *    (q(8)*al(1,irn13) - (q(12) + q(13))*al(1,irn16))/awght(3)
        epstil(idn14) = cepstl*
     *    (qn14(1)*al(1,irn14) - (q(12) + q(13))*al(1,irn16))/awght(3)
c
c  test for setting second derivatives
c
        if(secder) then
c
c  Second derivatives not implemented for icnocs = 4
c
        end if
c
c..     write(6,910) 'fcnotl(.,ixn14)',(fcnotl(k,ixn14),k=1,ider)
c..     write(6,910) 'epstil',(epstil(k),k=1,ider)
c..     write(6,910) 'fthtil',(fthtil(k),k=1,ider)
c
c  set CNO abundances in common/compos/
c
        xset12 = x(ixc12)
        xset13 = x(ixc13)
        xset14 = x(ixn14)
        xset16 = xtlo16*awght(6)
      end if
c
c  set derivatives of final rates of change of CNO abundances
c  contribution to rate of change of X and energy generation
c
      do 70 k=nspech+1, nspec
   70 call mltder(fcnotl(1,k), rho, xh, ft(1,k), ider1, ider, secder)
c..      write(6,910) 'ft(.,nspec)',(ft(k,nspec),k=1,ider)
c
      call mltder(fthtil, rho, xh, fthcno, ider1, ider, secder)
      call mltder(epstil, rho, xh, epscno, ider1, ider, secder)
c..      write(6,910) 'fthcno',(fthcno(k),k=1,ider)
c..      write(6,910) 'Before transformation, epscno',
c..     *  (epscno(k),k=1,ider)
c
c
c  change to ergs/(g*sec), transform to logarithmic derivatives
c
      if(epscno(1).gt.0) then
        do 72 i=2, ider1
   72   epscno(i) = epscno(i)/(amm*epscno(1))
c
        if(secder) then
          jj = ider1
          do 74 i=2, ider1
          do 74 j=i, ider1
          jj = jj+1
   74     epscno(jj) = epscno(jj)/(amm*epscno(1)) 
     *      - amm*epscno(i)*epscno(j)
        end if
      else
	do 75 i=2,ider
   75   epscno(i)=0
      end if
c
      epscno(1)=1.d6*ergev*epscno(1)
c
c  set final quantities in common /cepcno/
c
      epscnc=epscno(1)
      do 80 i=1,4
   80 eppcno(i)=1.d6*ergev*rho(1)*xh*eppcno(i)
c
      if(idgeng.eq.-10.and.istdpr.gt.0) 
     *  write(istdpr,*) 'epscno =',epscno(1)
c
      
      if(idgeng.eq.-42.and.istdpr.gt.0) then 
        write(istdpr,*) 'After transformation, fthcno:'
        write(istdpr,*) (fthcno(i),i=1,ider1)
      end if
c
      return
  190 format(//' ***** Error in s/r engcno. icnocs =',i5,
     *  ' not inplemented')
  910 format(/a60/(1p7e11.3))
      end
      subroutine mltder(atil, rho, xh, a, ider1, ider, secder)
c
c  Sets product rho*X*atil and its derivatives, stored in usual
c  manner, from atil and rho.
c  ider1 - 1 is number of first derivatives.
c
c  Original version: 29/7/91
c
      implicit double precision (a-h, o-z)
      logical secder
      dimension atil(*), rho(*), a(*)
      common/ln10/ amm,amm2
c
      istore = 0
      if(istore.eq.1) then
        do 10 i=1,ider
   10   a(i)=atil(i)
        return
      end if
c
      a(1) = rho(1)*xh*atil(1)
      do 15 i = 2,ider1
      if(i.le.4) then
        a(i) = rho(1)*xh*atil(i) + amm*rho(i)*a(1)
      else
        a(i) = rho(1)*xh*atil(i) 
      end if
   15 continue
      a(4) = a(4) + rho(1)*atil(1)
c
c  test for second derivatives
c
      if(secder) then
        jj = ider1
        ii = 4
        do 30 i = 2, ider1
        do 30 j = i, ider1
        jj = jj+1
        a(jj) = rho(1)*xh*atil(jj)
        if(i.le.4) then
          if(i.eq.4) a(jj) = a(jj) + rho(1)*atil(j)
          a(jj) = a(jj) + amm*rho(i)*a(j)
          if(j.le.4) then
            ii = ii+1
            a(jj) = a(jj) + 
     *         amm*(rho(ii)*a(1) + rho(1)*rho(j)*xh*atil(i))
            if(j.eq.4) a(jj) = a(jj) + rho(1)*atil(i)
            if(i.eq.4) a(jj) = a(jj) + amm*rho(1)*rho(j)*atil(1)
          end if
        end if
   30   continue
      end if
      return
      end
