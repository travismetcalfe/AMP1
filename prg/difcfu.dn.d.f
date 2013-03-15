      subroutine difcfu(t,rho,p,xh,zh,r,amass,grad,dad,
     *  totmss, age, itbdif, rdift, dift, noder)
c
c  Sets approximations to turbulent diffusion coefficient.
c  This version is designed for user modification. 
c  It is called by the programme when itbdif .gt. 100.
c  On the first call, additional parameters may be read from 
c  input file; this happens after all other input parameters to
c  the programme have been read.
c
c  The routine is called at each meshpoint in the model, during the
c  iteration to solve the evolution equations.
c
c  On input:
c
c    t: temperature, in K.
c    rho(1-20) is density array; rho(1) is density, in g/cm**3
c    p(1-10) is pressure array; p(1) is pressure, in dyn/cm**2
c    xh: hydrogen abundance
c    zh: heavy-element abundance s
c    r: distance from centre, in cm
c    amass: mass inside point,in g
c    grad: temperature gradient array; grad(1) = d log T/d log p
c    dad: adiabatic temperature gradient array; 
c         dad(1) = (d log T/d log p)_ad
c    totmss: total mass of star, in g
c    age: age of star, in years
c    rhob:  density at base of convective envelope, in g/cm**3
c    itbdif: parameter (.gt.100) which may be used to control
c            calculation
c    noder: logical variable, not currently used
c
c  Returns results in arrays rdift and dift;
c  only rdift(1) and dift(1) are currently used:
c
c  rdift(1): turbulent diffusion coefficient D, in cm**2/sec,
c            relevant for diffusion equation written in terms of r
c  dift(1):  curly D = ((4*pi*r**2*rho)**2)*D, i.e., coeffient
c            relevant for diffusion equation written in terms of m
c          
c  Original version: 3/11/2000
c
c  ***************************************************************
c
c  As an example, this version of the routine assumes a diffusion
c  coefficient of the form
c  
c  D = D_0 * exp( -((T - T_0)/delta T)**2)
c
c  where D_0, T_0 and delta T are read in
c
c  ***************************************************************
c
      implicit double precision(a-h, o-z)
      logical noder
      dimension rho(1), p(1), grad(1), dad(1), rdift(1), dift(1)
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev1,exhm
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,inccor,frcf(6),frcl(6)
      common/convvr/ cvvarf(7,6), cvvarl(7,6)
      common/cdiffu/ idiffus, itbdf1, rhob, rmxdif, tlfdif, dtldif,
     *  ctdfmx, rctdif, rctdf1
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev, syear, iver
      common/ln10/ amm,amm2,amm3
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      data idiag /0/
c
c  initialize flag for parameter input
c
      data intpar /1/
c
      save
c
c  internal function to calculate exp without over- or underflow
c
      fxp(a)=exp(min(85.d0,max(a,-85.d0)))
c
      if(intpar.eq.1) then
c
c  Read parameters from input
c
        if(istdpr.gt.0) write(istdpr,110) 
        read(istdin,*) d0, tmp0, deltmp
c
        if(istdpr.gt.0) write(istdpr,120) itbdif, d0, tmp0, deltmp
        open(99,file='difcfu.log',status='unknown')
        write(99,120) itbdif, d0, tmp0, deltmp
	close(99)
c
	intpar=0
c
      end if
c
      rdift(1) = d0*fxp(-((t-tmp0)/deltmp)**2)
c
c  set diffusion coefficient wrt mass
c
      dmdr=4*pi*r*r*rho(1)
      dift(1)=dmdr*dmdr*rdift(1)
      return
  110 format(//' Input parameter for s/r difcfu'/)
  120 format(//
     *  ' In s/r difcfu, itbdif =',i5/
     *  ' The following parameters have been set:'//
     *  ' D_0:     ', 1pe13.5, ' cm**2/sec'/
     *  ' T_0:     ', e13.5,' K'/
     *  ' delta T: ', e13.5,' K'/)
      end
