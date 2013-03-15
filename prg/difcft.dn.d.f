      subroutine difcft(t,rho,p,xh,zh,r,amass,grad,dad,
     *  totmss, age, itbdif, rdift, dift, noder)
c
c  sets approximations to turbulent diffusion coefficient,
c
c  On input t is temperature, rho is density array, p is pressure array,
c  (both stored in standard form for thermodynamic variables)
c  xh and zh are hydrogen and heavy-element abundances, r is
c  distance from centre, amass is mass inside point (in g),
c  grad is temperature gradient array, stored as for velx and difx
c  (see below), dad is adiabatic gradient array, totmss is total
c  mass of star, age is age (in years), rhob is density at
c  base of convection zone.
c
c  Returns turbulent diffusion (in cm**2/sec) in rdift(1),
c  and coefficient appropriate for mass (curly d) in
c  dift(1). Derivatives may be added later.
c
c  The type of diffusion depends on itbdif.
c 
c  itbdif = 1: Proffitt match to Pinsonneault et al
c  (see Proffitt and Michaud, ApJ 380, 238).
c
c  itbdif = 2: simple power law in density
c
c  itbdif = 3: gaussian decreasing with distance to 
c  convectively mixed core, gaussian width parameter rctdif
c
c  itbdif = 4: constant and then gaussian decreasing with distance to 
c  convectively mixed core, constant over rctdf1, 
c  gaussian width parameter rctdif 
c
c  The width parameters may be either in cm or in units of the
c  pressure scale height at the edge of the convective core.
c  If rctdif and rctdf1 are above 1000, a value in cm is assumed.
c  Otherwise the distances are set to rctdif*hpcor and
c  rctdf1*hpcor, where hpcor is the smaller of the pressure scale
c  height and the radius of the core (to avoid problems when the 
c  pressure scale height becomes unphysically small).
c
c  itbdif = 5: no action in this routine. diffusion coefficient
c  scaled by ctdfmx in s/r difcff, to mimick Schatzman et al.
c
c If noder is true, no derivatives are set
c
c  Original version: 19/8/92 
c
c  Modified 30/12/99, to allow specifying width parameters in units
c  of pressure scale height.
c
c  Modified 7/8/04, resetting rhob to value at base of overshoot region
c  in case with overshoot.
c
      implicit double precision(a-h, o-z)
      logical noder, setder
      dimension rho(1), p(1), grad(1), dad(1), rdift(1), dift(1)
      dimension tfact(4), dmdr(6), difxr(6)
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev1,exhm
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,inccor,frcf(6),frcl(6)
      common/convvr/ cvvarf(7,6), cvvarl(7,6)
      common/cnvovs/ icnvos, jcnvos, clcovs, cldovs, alphos,
     *  rczl, rczlfx, rcnvos, qlcnos, rhobos
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
      setder=.false.
c
c..      write(6,*) 'Enter difcft with rmxdif =',rmxdif
c
c  test for resetting rhob with value at base of overshoot region
c
      if(icnvos.eq.1.or.icnvos.eq.2) rhob=rhobos
c
c  test for setting distances for turbulent mixing
c
      if(itbdif.ge.3) then
	hpcor=min(cvvarf(7,inccor),rmxdif)
	if(rctdif.ge.1.e3) then
	  rctd=rctdif
        else
	  rctd=rctdif*hpcor
        end if
	if(rctdf1.ge.1.e3) then
	  rctd1=rctdf1
        else
	  rctd1=rctdf1*hpcor
        end if
	if(idiag.gt.0.and.mod(idiag,10).eq.1.and.istdpr.gt.0) then
          write(istdpr,*) 'In difcft, hpcor, rctd =',hpcor,rctd
        end if
	idiag=idiag-1
      end if
c
      if(itbdif.eq.1) then
	age1=max(age,5.d7)
	qq=amass/totmss
	if(age.lt.3.d7.or.qq.lt.0.5) then
	  rdift(1)=0
        else
	  ddadf=min(0.2d0,max(dad(1)-grad(1),0.015d0))
	  agef1=10.d0**(11.94-0.95*log10(age))
	  if(qq.gt.0.95) then
	    rdift(1)=5.*ddadf*agef1
          else
	    agef2=10.d0**(-3.16+0.48*log10(age))
	    rdift(1)=5.*ddadf*agef1*10.d0**((qq-0.95)*agef2)
          end if
        end if
c
      else if(itbdif.eq.2) then
	if(rhob.le.0) then
	  rdift(1)=ctdfmx
        else
	  rdift(1)=ctdfmx*(max(1.d0,rho(1)/rhob))**(-3.)
        end if
      else if(itbdif.eq.3) then
	if(rmxdif.le.0.or.r.le.rmxdif) then
	  rdift(1)=0
        else
	  xx=(r-rmxdif)/rctd
	  rdift(1)=ctdfmx*exp(-min(xx*xx,100.d0))
        end if
      else if(itbdif.eq.4) then
	if(rmxdif.le.0) then
	  rdift(1)=0
        else if(r-rmxdif.le.rctd1) then
	  rdift(1)=ctdfmx
        else
	  xx=(r-rmxdif-rctd1)/rctd
	  rdift(1)=ctdfmx*exp(-min(xx*xx,100.d0))
        end if
c..	write(6,'(a,1p4e13.5)') 'r,rmxdif,xx,rdift(1)',r,rmxdif,xx,rdift(1)
      else
	rdift(1)=0
      end if
c
c  set diffusion coefficient wrt mass
c
      dmdr(1)=4*pi*r*r*rho(1)
      dift(1)=dmdr(1)*dmdr(1)*rdift(1)
      return
      end
