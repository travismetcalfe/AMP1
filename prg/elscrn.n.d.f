      subroutine elscrn(fl, tl, x, y, z, iscren, 
     *  zt, uw, zzscr, thte, nosd)
c
c  Calculates electron screening for nuclear reactions,
c  using weak or intermediate screening as determined by iscren.
c
c  iscren = 1: use weak screening in original formulation, with
c     thte = thtec (as given in input)
c  iscren = 2: use weak screening, with consistent thte
c  iscren = 3: use intermediate screening, as done by Bahcall.
c  iscren = 4: use intermediate screening in Mitler formulation,
c     as used by Turck-Chieze.
c
c  Note: eqstf must have been called before call of elscrn
c
c  uw(1) returns screening function, and zzscr(k) returns charge
c  factor for reaction no. k, such that the screening factor
c  is fscr(k) = exp(zzscr(k)*uw(1))
c  uw(2 - 10) return d log fscr / d t(i) = (1/amm) d uw / d t(i)
c  and d2 log fscr / d t(i) d t(j) = (1/amm) d2 uw / d t(i) d t(j)
c  where amm = ln(10)
c
c  zt(1) returns zeta, and zt(2 - 10) return d log zeta/d t(i)
c  and d2 log zeta/d t(i)d t(j), where, as usual, 
c  t(i) = (log f, log T, X). Also, log is to base 10.
c
c  thte(1) returns theta_e (as defined by Grabowske et al), and
c  thte(2 - 10) return d log thte / d t(i) and 
c  d2 log thte/d t(i) d t(j). Note that as thte is a function of
c  psi, and hence f, alone, most of the derivatives are zero.
c  However, for consistency the same storage of derivatives 
c  as for other thermodynamic quantities is maintained
c
c  Orignial version: 19/2/93
c
c  Modified 17/4/98, to include Mitler form of weak screening
c  (see Dzitko et al. 1995; ApJ 447, 428)
c  This is done by setting up zzscr, appropriately.
c  Note that consistent derivatives have not yet been set.
c  These would have to be defined in an extension to zzscr.
c
c  Modified 18/4/98, to set fl if not provided, from log(rho) and
c  N_e.
c
      implicit double precision (a-h,o-z)
      logical noder, nosd, secder
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
c  Note: engenr.n.d.incl replaced by engenr.nnz.d.incl, 5/6/02
      include 'engenr.nnz.d.incl'
      parameter (idermx = ((nspcmx+3)*(nspcmx+4))/2, idalt = 10,
     *  krnrm1 = krnrmx + 1, iptdat = krnrm1 - 12)
      dimension zt(*), uw(*), zzscr(*), thte(*)
      dimension iptrnr(krnrm1, 4),phi(30), hst(30)
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm
      common/rnrcnt/ irnrat, krnrat, irn12, irn13, irn14, 
     *  irn15c, irn15o, irn16, irn17
      common/rcncns/ a(knucst),b(knucst),s(knucst,isnuc),
     *  zz(knucst),zz12(2,knucst),zsmh,acno,awght(knucwg),q(knucq), 
     *  knucpm, kqvals, kawght
      common/eqstd/ xii(4),ane(10),rho(10)
      common/degfct/ thtec
      common/ln10/ amm,amm2
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      external blengr
c
c  set pointing indices to data in common /rcncns/
c
c                  1  2  3  4  5  6  7  8   9  10  11  12
      data iptrnr /1, 2, 3, 4, 5, 0, 0, 0,  0,  0,  0,  0, iptdat*0,
     *             1, 2, 3, 4, 5, 0, 8, 9, 10, 11, 10, 11, iptdat*0,
     *             1, 2, 3, 4, 5, 0, 6, 7,  0,  0,  0,  0, iptdat*0,
     *             1, 2, 3, 4, 5, 0, 6, 7,  8,  9, 10, 11, iptdat*0/
c
c  for the moment, hard-code bint, akint and z3bm1
c
      data bint, akint, z3bm1 /0.86, 0.38, 1.55/
c
      iscnuc=mod(iscren,10)
      secder=.not.nosd
c
c  test for no screening
c
      if(iscnuc.eq.0) then
        do 10 i=1,10
	thte(i)=0
        zt(i)=0
   10   uw(i)=0
        do 12 k=1,krnrat
   12   zzscr(k)=0
	return
      end if
c
c  test for setting theta
c
      if(iscnuc.eq.1) then
	thte(1)=thtec
	do 14 i=2,10
   14   thte(i)=0
c
      else
c
c  call phder to set degeneracy functions. 
c  Note: for some reason the original implementation (Dec. 93)
c  set a very low temperature to use nonrelativistic formulation.
c  This has been replaced (18/4/98) by the actual temperature
c
c  Possibly iterate for log(f), if not provided
c
	ttph=1.d3
	tlph=tl
	ttph=10.d0**tlph
	flph=fl
	rhol=log10(rho(1))
	noder=.false.
c
c  test for case
c
	if(abs(flph-rhol).ge.1.e-3) then
c
c  assume that fl contains log(f) (test may need refinement)
c
	  call phder(flph, tlph, phi, hst, noder, nosd)
c
        else
c
c  iterate for appropriate log(f)
c
	  rho00=2.3152d0
	  ts=10.d0**tl/ct
          ff=rho(1)*ane(1)/(rho00*crho*ts**1.5)
	  flph=log10(ff)
	  nitfl=0
c
   15     call phder(flph, tlph, phi, hst, noder, nosd)
	  rhoph=crho*phi(1)/ane(1)
          dlogrh=log10(rho(1)/rhoph)
	  if(nitfl.gt.10) then
	    write(istdou,120) tlph,log10(rho(1)),flph,dlogrh
	    if(istdpr.gt.0.and.istdpr.ne.istdou) 
     *          write(istdpr,120) tlph,log10(rho(1)),flph,dlogrh
	  else if(abs(dlogrh).gt.1.e-8) then
	    flph=flph+dlogrh/phi(2)
	    go to 15
          end if
        end if
c
	ff=10.d0**flph
	ff1=1.d0+10.d0**flph
	sqff1=sqrt(ff1)
        thte(1)=phi(2)/sqff1
c
c  set derivatives
c
        do 16 i=2,10
   16   thte(i)=0
	thte(2)=(phi(4) - 0.5*amm*phi(2)*ff/ff1)/sqff1
	thte(5)=(phi(7) - amm*phi(4)*ff/ff1
     *         + amm*amm*phi(2)*ff*(0.75*ff/ff1 - 0.5)/ff1)/sqff1
      end if
c
c  set zeta 
c
      zzt=x/ah+4*y/ahe+zsmh*z+thte(1)*ane(1)/av
      zt(1)= sqrt(zzt)
      zzt=1/(2*zzt*amm)
      do 20 k=2,4
   20 zt(k)=(thte(k)*ane(1) + thte(1)*ane(k))*zzt/av
      zt(4)=zt(4)+(1/ah-4/ahe)*zzt
      if(secder) then
        ii=4
        do 24 i=2,4
        do 24 j=i,4
        ii=ii+1
   24   zt(ii)=(thte(1)*ane(ii) + thte(i)*ane(j) + thte(j)*ane(i)
     *        +thte(ii)*ane(1))*zzt/av - 2*zt(i)*zt(j)/amm
      end if
c
      t = 10.d0**tl
      t62= sqrt(t)/1000
      t623=t62*t62*t62
c
c  set uw and zzscr, depending on iscnuc
c  For test, we currently set *either* weak or intermediate
c  screening. Later a suitable transition must be implemented.
c
      if(iscnuc.ne.3) then
c
c  weak screening
c
        uw(1)=0.188d0* sqrt(rho(1))*zt(1)/t623
        do 30 k=2,4
   30   uw(k)=uw(1)*(rho(k)/2+zt(k))
        uw(3)=uw(3)-1.5d0*uw(1)
        if(secder) then
          ii=4
          do 35 i=2,4
          do 35 j=i,4
          ii=ii+1
   35     uw(ii)=uw(1)*(rho(ii)/2 + zt(ii)) + amm*uw(i)*uw(j)/uw(1)
        end if
c
c  charge factor
c
	do 37 k=1,krnrat
	if(k.ne.6) then
	  kpar=iptrnr(k, irnrat)
	  zzscr(k)=zz(kpar)
        end if
   37   continue
c
c  test for applying Mitler correction
c
	if(iscnuc.eq.4) then
c
c  apply Mitler correction (note: derivatives have to be added)
c
          zeta0=0.56372d0*sqrt(rho(1))*zt(1)**3/(t623*ane(1)/av)
	  do 38 k=1,krnrat
	  if(k.ne.6) then
	    kpar=iptrnr(k, irnrat)
	    zeta1=zz12(1,kpar)*zeta0
	    zeta2=zz12(2,kpar)*zeta0
	    zetas=zeta1+zeta2
c
c  test for using expansion
c
	    if(zetas.le.1.e-3) then
	      fctmit=1-zetas/6.d0
            else
	      pow=5.d0/3.d0
	      fctmit=0.9d0*((1+zetas)**pow - (1+zeta1)**pow
     *               - (1+zeta2)**pow+1)/(zeta1*zeta2)
	    end if
	    zzscr(k)=fctmit*zzscr(k)
          end if
   38     continue
c
	end if
c
      else
c
c  Bahcall intermediate screening
c
	b31=3*bint-1
	b32=3*bint-2
	bm22=2-2*bint
        ztb31= x/ah+2.d0**b31*y/ahe+z3bm1*z
	etaint=ztb31/(zt(1)**b32*(ane(1)/av)**bm22)
	uw(1)= akint*etaint*(0.188d0* sqrt(rho(1))/t623)**bint
c
c  derivatives 
c
        do 40 k=2,4
   40   uw(k)=uw(1)*(bint*rho(k)/2 - b32*zt(k) 
     *       - bm22*ane(k)/(amm*ane(1)))
        uw(3)=uw(3)-1.5d0*bint*uw(1)
	dztb31=1/ah - 2.d0**b31/ahe
	uw(4)=uw(4) + dztb31/(amm*ztb31)
        if(secder) then
          ii=4
          do 45 i=2,4
          do 45 j=i,4
          ii=ii+1
   45     uw(ii)=uw(1)*(bint*rho(ii)/2 - b32*zt(ii)
     *	        - bm22*(ane(ii) - ane(i)*ane(j)/ane(1))/(amm*ane(1)))
     *          + amm*uw(i)*uw(j)/uw(1)
	  uw(10)=uw(10) - uw(1)*(dztb31*dztb31)/(amm*ztb31*ztb31)
        end if
c
c  charge factor
c
	do 55 k=1,krnrat
	if(k.ne.6) then
	  kpar=iptrnr(k, irnrat)
	  zzscr(k)=((zz12(1,kpar)+zz12(2,kpar))**(1+bint)
     *      - zz12(1,kpar)**(1+bint) - zz12(2,kpar)**(1+bint))
        end if
   55   continue
c
      end if
c
      if(idgeng.ge.2.and.istdpr.gt.0) then
        write(istdpr,*) 'x,ah,y,ahe,zsmh,thte(1),ane(1),av'
        write(istdpr,*) x,ah,y,ahe,zsmh,thte(1),ane(1),av
        write(istdpr,220) zt(1), uw(1)
      end if
c
      return
c
  120 format(//' **** Error in s/r elscrn. Iternation for log(f)',
     *  ' failed to converge.'/
     *  '             log(T), log(rho), last log(f), error:',
     *  1p4e13.5)
  220 format(' electron screening with zeta =',1pe13.5,
     *  '  Uw = ',e13.5)
      end
