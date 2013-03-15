      subroutine setcns
c
c  sets the physical and mathematical constants used in the programme
c
c  numerical constants from CODATA report, with additional values
c  for atomic masses and ionization potentials. this is
c  the GONG set of values.
c
c  corresponds to iver = 2
c
c  modified on 8/12/87 to approach a consistent setting based on
c  minimum set of basic constants. so far keep original values,
c  but print consistent values
c
c  modified 15/3/88 to use consistent numerical constants
c  note that a version number iver has been added in common/fconst/
c  to indicate which set of values have been used.
c
c  also include constants for simple equation of state, with
c  no partial degeneracy, for use for GONG models.
c
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 9/3/90
c
c  Modified 24/8/96, to read in file name for Livermore tables,
c  setting EOS table version number
c  Also sets up table version number based on list coded in
c  routine
c
c  Modified 24/3/04, to read in two Livermore Z tables and
c  associated values of Z, for linear Z interpolation using
c  esac-liv01z.nnz.d.f
c
      implicit double precision(a-h,o-z)
      parameter (nztab=2)
      character*80 finliv
      character fin*60, cindx*9
      dimension indx(6), cindx(6)
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev1,exhm
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev, syear, iver
      common/eqphsm/ rho00, f0
      common/ln10/ amm,amm2,amm3
      common/eqstd/ ccc1(90)
      common/eqsout/ ccc2(210)
      common/dmuder/ ccc3(11)
      common/ctablv/ iztab, zinliv(nztab), finliv(nztab)
      common/cversn/ ivrbcs,ivrrhs,ivreos,ivropc,ivreng
      common/ctabvr/ ivteos
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      data init/0/
      data indx, nindx  /101,  102, 103, 111, 112, 113, 6/
      data cindx        /'data.Z2','data.Z18', 'data.Z1', 
     *  '2001.Z2','2001.Z18','2001.Z1'/
c
      save
c
      idiag=1
c
      if(istdpr.gt.0) write(istdpr,100)
c
c  version number for constants
c
      iver = 2
c
c  *****************************************************************
c
c  set equation of state version number to 8 for Livermore interpolation
c
      ivreos = 8
c
c  *****************************************************************
c
c  test for reading name of Livermore tables
c
      if(init.eq.0) then
c
	init=1
c
	if(istdpr.gt.0) write(istdpr,*) 
     *    'Enter file names for Livermore tables'
	read(istdin,'(a)') finliv
	read(istdin,*) zinliv
c
c  set index number
c
        do 20 i=1,nindx
        if(index(finliv(1),cindx(i)).gt.0) then
	  ivteos=indx(i)
	  go to 25
        end if
   20   continue
c
        ivteos=100
	if(istdpr.gt.0) write(istdpr,105) finliv,ivteos,zinliv
   25   continue
c
        if(idiag.gt.0.and.istdpr.gt.0) write(istdpr,107) 
     *    finliv, ivteos, zinliv
c
      end if
c
c  ********************************************************************
c
      pi=4.d0*atan(1.d0)
      ebase=exp(1.d0)
c
c  atomic mass unit
      amu=1.6605402d-24
c  electron mass
      ame=9.1093897d-28
c  speed of light
      clight=2.99792458d10
c  planck's constant
      planck=6.6260755d-27
c  boltzman's constant, cgs
      boltzm=1.380658d-16
c  gravitational constant (chosen from planetary orbit data, assuming
c  msun = 1.989d33. consistent, but not identical, with CODATA value).
      cgrav=6.67232d-8
c  solar mass 
      amsun=1.989d33
c  electron charge, in ESU
      echar=4.8032068d-10
c  number of ergs in 1 ev
      ergev=1.60217733d-12
c  number of seconds in a year
c  (should be replaced with second value)
      syear =3.155692597d7
c
c  start on derived values, largely
c
      if(istdpr.gt.0) write(istdpr,110)
c  avogadro's number (should be 1/amu)
      avp=6.02217d23
      av=1.0/amu
      if(istdpr.gt.0) write(istdpr,120) 'av',avp,av
c  atomic weights of h and he
      ah=1.007825d0
      ahe=4.002603d0
c  average atomic weight for heavy elements.
c  before 4/1/84 was given value 17.8, but was in fact
c  most often reset to the current value in calling programme.
      az=16.389d0
c  auxiliary quantities for eggleton fudge calculation
      avda=av*(1/ah-2/ahe)
      avd1=av*(1/ah-1/ahe)
c  boltzmann's constant in ev/deg
      ck1p=8.6170837d-5
      ck1=boltzm/ergev
      if(istdpr.gt.0) write(istdpr,120) 'ck1',ck1p,ck1
c  the same, in ergs/deg
      ck2=boltzm
c  ionization potentials
      exh=13.595d0
      exhe=24.580d0
      exhep=54.403d0
c  constants for transition to starred variables
      ctp=1.686304d-10
      crhop=1.759547d30
      cpep =1.440588d24
      chep =8.187265d-7
c
      alamc=planck/(ame*clight)
      alamc3=alamc**3
      ct=boltzm/(ame*clight**2)
      crho=8*pi/alamc3
      che=ame*clight**2
      cpe=crho*che
      if(istdpr.gt.0) then
        write(istdpr,120) 'ct',ctp,ct
        write(istdpr,120) 'crho',crhop,crho
        write(istdpr,120) 'cpe',cpep,cpe
        write(istdpr,120) 'che',chep,che
      end if
c  constants for pressure ionization
      ca03p=2.147d-24
      caap=1.759547d30*ca03p
      ckhp=13.5d0
c
      a0=planck/(2*pi*echar)
      a0=(a0/ame)*a0
c
c  for the moment, set additional fudge factor by hand at this
c  point. need to think about a more suitable way to do it later
c
      efffac=15
      ca03=efffac*a0**3
      caa=crho*ca03
      ckh=exh
c
c  change from ev to ergs and include factor av**2 to compensate for
c  redefinition of dne (see notes of 3/1/84)
c
      ca03p=av*1.602192d-12*ca03p*av
      ca03=av*ergev*ca03*av
c
      if(istdpr.gt.0) then
        write(istdpr,120) 'ca03',ca03p,ca03
        write(istdpr,120) 'caa',caap,caa
        write(istdpr,120) 'ckh',ckhp,ckh
      end if
c
c  the radiation constant
c
      carp=7.5647d-15
      car=boltzm/(clight*planck)
      car=8*pi**5*boltzm*car**3/15
      if(istdpr.gt.0) write(istdpr,120) 'car',carp,car
c  number of ergs in 1 ev
      ergev1=ergev
c  ionization potential for h-
      exhm=0.754d0
c  ln 10
      amm=log(1.d1)
      amm2=amm*amm
      amm3=amm2*amm
c
c constant for phder, with no partial degeneracy
c
      rho00=sqrt(2*pi)*ebase**2/8
      f0=av/(rho00*crho*ct**1.5d0)
c
c  set commons from s/r eqstf to zero
      call zero(ccc1,90)
      call zero(ccc2,210)
      call zero(ccc3,10)
      return
  100 format(//2x,75('*')//
     *  ' Set constants from CODATA Report.'/' See Cohen & Taylor,',
     *  ' Rev. Mod. Phys., vol. 59, 1121 (1987)'/
     *  ' ++++ Double precision version'/)
  105 format(//
     *  ' ***** Warning in s/r setcns-liv. Table file not in index'/
     *  '       File names: '/a80/a80/
     *  '       Version number set to',i5/
     *  '       Z =',2f10.5)
  107 format(//' EOS data to be read from files:'/1x,a/1x,a/
     *  ' Version number set to',i5/
     *  ' Z =',2f10.5)
  110 format(/' variable   old value    current, consistent value'/)
  120 format(a10,1p2e18.7)
      end
