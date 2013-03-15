      subroutine difcff(t,rho,p,xh,zh,r,amass,grad,
     *  idiffus,itbdif,clamxy,velx,difx,idmdif,noder)
c
c  sets settling velocity velx, diffusion coefficient difx and
c  Coulomb logarithm clamxy for diffusion.
c
c  On input t is temperature, rho is density array, p is pressure array,
c  (both stored in standard form for thermodynamic variables)
c  xh and zh are hydrogen and heavy-element abundances, r is
c  distance from centre, amass is mass inside point (in g),
c  grad is temperature gradient array, stored as for velx and difx
c  (see below). 
c
c  if itbdif gt 0 sets turbulent diffusivity, from s/r difcft.
c  Currently try not to include derivatives of turbulent diffusivity.
c
c  Returns Coulomb logarithm in clamxy(1) and derivatives wrt
c  log10 f, log10 T and X in clamxy(2-4) (note that log Lambda is
c  a thermodynamic quantity and hence derivatives are
c  set accordingly).
c
c  Returns settling-velocity coefficients in velx(i,1,1) and velx(1,2,1),
c  diffusion coefficient c  in difx(i,1), and derivatives of velx and 
c  difx wrt log10 r, log10 f, log10 T, log10 L, X and Z in velx(i,k,2 - 7) 
c  and difx(i,2 - 7), for i running over diffusing species
c  (in conformance with storage used, e.g., in s/r mixlng;
c  note that further transformation is needed to conform with
c  order of variables in set solved in diffusion case).
c
c  More precisely, if V_i and X_i are diffusion velocity and 
c  abundance of species i, we have
c
c                                     d X_i
c  4 pi r^2 rho V_i X_i = - difx(i,1) ----- - velx(i,1,1) X_i -
c                                      d m
c
c                         - velx(i,2,1) M Y_H X_i
c
c  where Y_i = 4 pi r^2 rho V_i X_i M^(-1) and Y_H is Y_i for i = 1,
c  i.e., hydrogen. For i = 1, the term in velx(1,2,1) is evidently
c  absent.
c
c  In addition, for output purposes we set in common/cdiffv/ the
c  variables cdifr(.), cvelr(.,.) and cdiftr, such that
c
c                        d X_i
c  V_i X_i = - cdifr(i) ----- - cvelr(i,1) X_i - cvelr(i,2) M X V_H X_i
c                         d r
c
c  and cdiftr is the turbulent contribution to cdifr(i).
c
c  If tlfdif (in common/cdiffu/) is gt 0, suppress settling at
c  log T lt tlfdif
c
c  If noder is true, no derivaties are set
c
c  Based on expressions in Michaud & Proffitt 
c  (Proc. IAU Colloq. 137: Inside the stars;
c  ASP Conf. Ser. vol. 40, p. 246; 1993), obtained in part from
c  a bit of code obtained from Proffitt
c  (see /hosts/big_scr/usr/jcd/evolprg/dp/diffusion/f.proffitt.2008).
c  Note that the expression for velx(2,2,1) contains an additional 
c  term -0.23 not in that paper. The origin of that is unclear.
c
c  Original version: 15/8/92 
c
c  Modified 20/5/93, to enable suppression of diffusion 
c  at low temperature
c
c  Modified 5/8/95, to allow diffusion of several species.
c
c  Corrected 5/8/97, by including engenr.nnz.d.incl, to set nspdmx etc
c
c  Corrected 27/12/99, including suppression in settling velocity
c  for elements other than H when tlfdif gt 0
c
c  Modified 3/11/00, allowing for call of user-supplied routine
c  difcfu when itbdif .gt. 100
c
c  Modified 24/10/03, adding s/r difscn to set diffusion coefficient
c  for semiconvection.
c
c  Modified 26/10/03, preparing for the inclusion of Z derivatives
c  (not consistently implemented so far)
c
      implicit double precision(a-h, o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.nnz.d.incl'
      logical noder, setder, time0, lastmd
      dimension rho(1), p(1), grad(1), clamxy(1), 
     *  velx(idmdif,2,1), difx(idmdif,1)
      dimension tfact(7), dmdr(7), difxr(nspdmx,7), diftb(7), rdiftb(7),
     *  velt1(7), veltr2(7), velt2(7)
      common/eqstd/ xi(4),ane(10),rho1(20),ht(20),pt(20),cp(4),dad(4),
     .  dlt(4),gm1,tprh,trhp,rhxp
      common/totmss/ am, rs
      common/cmtime/ age, time0, lastmd
      common/cdiffu/ idiffs, itbdf1, rhob, rmxdif, tlfdif, dtldif,
     *  ctdfmx, rctdif, rctdf1
      common/cdiffv/ cdifr(nspdmx),cvelr(nspdmx,2),cdiftr
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev1,exhm
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
      data bconst /-1.d0/
      data diftxr /0.d0/
c
      setder=.not.noder
c
c  test for suppressing settling at low temperature
c
      tl=log10(t)
      if(tlfdif.gt.tl) then
	xx=(tlfdif-tl)/dtldif
	supfct=exp(-min(100.d0,xx*xx))
      else
	supfct=1
      end if
c
c  test for setting turbulent diffusivity
c
      if(itbdif.gt.0) then
	ageyr=age/syear
	totlms=am*amsun
	qq=amass/totlms
c..        write(6,'(a,1p10e13.5)') 
c..     *    'in difcff, t, rho, etc =',t, rho(1), totlms, ageyr, qq
c..	write(6,*) 'Calling difcft with itbdif =',itbdif
        if(itbdif.le.100) then
          call difcft(t,rho,p,xh,zh,r,amass,grad,dad,
     *      totlms, ageyr, itbdif, rdiftb, diftb, noder)
        else
          call difcfu(t,rho,p,xh,zh,r,amass,grad,dad,
     *      totlms, ageyr, itbdif, rdiftb, diftb, noder)
	end if
        diftxr=rdiftb(1)
c..        write(6,*) 'diftxr set to',diftxr
      else
	diftxr=0
      end if
c
c  as a temporary fudge, suppress turbulent diffusion if
c  hydrogen abundance is small
c
      if(xh.le.1.e-3) diftxr=0
c
      if(bconst.lt.0) then
        bconst=0.9375d0*sqrt(0.4d0*amu/pi)*boltzm**2.5d0/echar**4
	if(istdpr.gt.0) then
          write(istdpr,*) 'amu, pi, boltzm, echar',
     *      amu, pi, boltzm, echar
          write(istdpr,*) 'bconst =',bconst
        end if
      end if
c
      clamxy(1)=-19.95d0-0.5d0*log(rho(1))-0.5d0*log((xh+3.d0)/2.d0)
     *  +1.5d0*log(t)
      if(clamxy(1).lt.1) then
	clamxy(1)=0.5
	if(setder) call zero(clamxy(2),3)
      else if(setder) then
	clamxy(2)=-0.5d0*amm*rho(2)
	clamxy(3)=amm*(-0.5d0*rho(3)+1.5d0)
	clamxy(4)=-0.5d0*amm*rho(4)-0.5d0/(xh+3d0)
      end if
c
c  first set coefficients for hydrogen
c
      tfact(1)=bconst*t**2.5d0/(clamxy(1)*(0.7d0+0.3d0*xh))
c..      write(6,*) 'tfact(1) =',tfact(1)
      if(setder) then
	tfact(2)=-tfact(1)*clamxy(2)/clamxy(1)
	tfact(3)=tfact(1)*(-clamxy(3)/clamxy(1) + 2.5d0*amm)
	tfact(4)=-tfact(1)*(clamxy(4)/clamxy(1)+0.3d0/(0.7d0+0.3*xh))
      end if
c
      velx(1,1,1)=-4d0*pi*tfact(1)*
     *  (1.25d0+1.125d0*grad(1))*(rho(1)*cgrav*amass/p(1))*(1-xh-zh)
      dmdr(1)=4d0*pi*rho(1)*r*r
      difxr(1,1)=(tfact(1)/rho(1))*(3d0+xh)/((1+xh)*(3d0+5d0*xh))
c
c  for itbdif = 5, scale diffusion coefficient (but not
c  settling velocity) by ctdfmx, to mimick Schatzman et al.
c
      if(itbdif.eq.5) difxr(1,1)=ctdfmx*difxr(1,1)
c
      difx(1,1)=dmdr(1)*dmdr(1)*(difxr(1,1)+diftxr)
c
c  include possible factor for suppressing settling at low
c  temperature
c
      velx(1,1,1)=supfct*velx(1,1,1)
c
c  no second term for hydrogen
c
      velx(1,2,1)=0
c
      if(setder) then
        velx(1,1,2)=velx(1,1,1)*1.125d0*grad(2)/(1.25d0+1.125d0*grad(1))
        velx(1,1,3)=velx(1,1,1)*(tfact(2)/tfact(1)+
     *          1.125d0*grad(3)/(1.25d0+1.125d0*grad(1))+
     *          amm*(rho(2)-p(2)))
        velx(1,1,4)=velx(1,1,1)*(tfact(3)/tfact(1)+
     *          1.125d0*grad(4)/(1.25d0+1.125d0*grad(1))+
     *          amm*(rho(3)-p(3)))
        velx(1,1,5)=velx(1,1,1)*1.125d0*grad(5)/(1.25d0+1.125d0*grad(1))
        velx(1,1,6)=velx(1,1,1)*(tfact(4)/tfact(1)+
     *          1.125d0*grad(6)/(1.25d0+1.125d0*grad(1))+
     *          amm*(rho(4)-p(4))-1.d0/(1-xh-zh))
        velx(1,1,7)=0
c
	do 25 i=2,7
   25   velx(1,2,i)=0
c
	dmdr(2)=2.d0*amm*dmdr(1)
	dmdr(3)=amm*dmdr(1)*rho(2)
	dmdr(4)=amm*dmdr(1)*rho(3)
	dmdr(5)=0
	dmdr(6)=amm*dmdr(1)*rho(4)
	dmdr(7)=0
	difxr(1,2)=0
	difxr(1,3)=difxr(1,1)*(tfact(2)/tfact(1) - amm*rho(2))
	difxr(1,4)=difxr(1,1)*(tfact(3)/tfact(1) - amm*rho(3))
	difxr(1,5)=0
	difxr(1,6)=difxr(1,1)*(tfact(4)/tfact(1) - amm*rho(4) +
     *   1.d0/(3d0+xh) - 1.d0/(1d0+xh) - 5.d0/(3d0+5d0*xh))
	difxr(1,7)=0
	do 30 i=2,7
   30   difx(1,i)=dmdr(1)*(2.d0*dmdr(i)*(difxr(1,1)+diftxr)+
     *            dmdr(1)*difxr(1,i))
      end if
c
c  store values in common/cdiffv/
c
      cdifr(1)=difxr(1,1)+diftxr
      cvelr(1,1)=velx(1,1,1)/dmdr(1)
      cvelr(1,2)=0
      cdiftr=diftxr
c
      if(idiffus.eq.1) return
c
c  coefficients for other elements, depending on idiffus
c
c  For the moment, represent heavy-element diffusion by oxygen
c
      zhv = 8.d0
      amhv = 16.d0
c
      axzs=dsqrt((ah*amhv)/(ah+amhv))
      ayzs=dsqrt((ahe*amhv)/(ahe+amhv))
      cxz=log(dexp(1.2d0*(clamxy(1)+log(2.d0/zhv)))+1.0d0)/1.2d0
      cyz=log(dexp(1.2d0*(clamxy(1)-log(zhv)))+1.0d0)/1.2d0
      cc1=xh*(axzs*cxz-ayzs*cyz)
      cc=cc1+ayzs*cyz
c
      tfact(1)=bconst*t**2.5d0
c
c leading coefficient for many terms 
c
      fccvz=0.8944272d0/(cc*zhv*zhv)
c
c  diffusion coefficient
c
      difxr(2,1)=fccvz*tfact(1)/rho(1)
c
c  add turbulent part (as for hydrogen)
c
      difx(2,1)=dmdr(1)*dmdr(1)*(difxr(2,1)+diftxr)
c
      vlt1fc=4d0*pi*tfact(1)*rho(1)*cgrav*amass/p(1)
      grt1fc=0.54d0*(4.75d0*xh+2.25d0)/(clamxy(1)+5.0d0)
      velt1(1)=-vlt1fc*(fccvz*(1.0d0+zhv-amhv*0.25d0*(5.d0*xh+3.d0)) -
     *   grt1fc*grad(1))
c
      veltr2(1)=-fccvz*(tfact(1)/rho(1))*
     *          (2.d0*zhv+5.d0*(1.d0+xh))/((1.d0+xh)*(5.d0*xh+3.d0))
      velt2(1)=dmdr(1)*dmdr(1)*veltr2(1)
c
c  first settling coefficient
c
      velx(2,1,1)=velt1(1)
     *   - (velx(1,1,1)/difx(1,1))*xh*velt2(1)
c
c  second settling coefficient
c
      velx(2,2,1) = velt2(1)/difx(1,1) + cc1/(xh*cc) - 0.23d0
c
c  When setting derivatives, ignore derivatives of cxz etc,
c  and be sloppy with composition derivatives
c
      if(setder) then
c
	do 35 i=2,7
   35   tfact(i)=0
c
	tfact(3)=2.5d0*amm*tfact(1)
c
	difxr(2,2)=0
	difxr(2,3)=-amm*difxr(2,1)*rho(2)
	difxr(2,4)=difxr(2,1)*(tfact(3)/tfact(1)-amm*rho(3))
	difxr(2,5)=0
	difxr(2,6)=-amm*difxr(2,1)*rho(4)
	difxr(2,7)=0
c
	veltr2(2)=0
	veltr2(3)=-amm*veltr2(1)*rho(2)
	veltr2(4)=veltr2(1)*(tfact(3)/tfact(1)-amm*rho(3))
	veltr2(5)=0
	veltr2(6)=-amm*veltr2(1)*rho(4)
	veltr2(7)=0
c
	do 40 i=2,7
        difx(2,i)=dmdr(1)*(2.d0*dmdr(i)*(difxr(2,1)+diftxr)+
     *            dmdr(1)*difxr(2,i))
   40   velt2(i)=dmdr(1)*(2.d0*dmdr(i)*veltr2(1)+
     *            dmdr(1)*veltr2(i))
c
	difx(2,7)=0.d0
c
        velt1(2)=vlt1fc*grt1fc*grad(2)
        velt1(3)=amm*velt1(1)*(rho(2)-p(2))+vlt1fc*grt1fc*grad(3)
        velt1(4)=velt1(1)*(tfact(3)/tfact(1)+amm*(rho(3)-p(3)))
     *          +vlt1fc*grt1fc*grad(4)
        velt1(5)=vlt1fc*grt1fc*grad(5)
        velt1(6)=amm*velt1(1)*(rho(4)-p(4))+vlt1fc*grt1fc*grad(6)
        velt1(7)=0
	do 45 i=2,7
        velx(2,1,i)=velt1(i)
     *     - (velx(1,1,i)/difx(1,1))*xh*velt2(1)
     *     + (velx(1,1,1)/difx(1,1))*(difx(1,i)/difx(1,1))*xh*velt2(1)
     *     - (velx(1,1,1)/difx(1,1))*xh*velt2(i)
   45   velx(2,2,i) = velt2(i)/difx(1,1) 
     *     - (velt2(1)/difx(1,1))*(difx(1,i)/difx(1,1))
c
c  add X-derivative
c 
	velx(2,1,6)=velx(2,1,6) - (velx(1,1,1)/difx(1,1))*velt2(1)
      end if
c
c  include suppression factor, if tlfdif gt 0
c
      if(tlfdif.gt.tl) then
	do 50 j=1,2
	do 50 i=1,6
   50   velx(2,j,i)=supfct*velx(2,j,i)
      end if
c
c  store values in common/cdiffv/
c
      cdifr(2)=difxr(2,1)+diftxr
      cvelr(2,1)=velx(2,1,1)/dmdr(1)
      cvelr(2,2)=-velx(2,2,1)
c
      return
      end
      subroutine difscn(t,rho,p,xh,zh,r,amass,grad,ak,dtxh,difscx,
     *  noder)
c
c  sets diffusion coefficient difscx from semiconvection.
c  So far derivatives are set only for dependence on nabla.
c  Other derivatives need to be added (sigh!)
c
c  On input t is temperature, rho is density array, p is pressure array,
c  (both stored in standard form for thermodynamic variables)
c  xh and zh are hydrogen and heavy-element abundances, r is
c  distance from centre, amass is mass inside point (in g),
c  grad is temperature gradient array, ak is log10(opacity) and
c  dtxh is composition contribution to superadiabatic gradient.
c
c  Based on expressions quoted by Richard et al. (2001; ApJ 558, 377).
c
c  Original version: 23/10/03
c
      implicit double precision(a-h, o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
      logical noder, setder, time0, lastmd
      dimension rho(1), p(1), grad(1), difscx(1)
      dimension difsc(7)
      common/eqstd/ xi(4),ane(10),rho1(20),ht(20),pt(20),cp(4),dad(4),
     .  dlt(4),gm1,tprh,trhp,rhxp
      common/totmss/ am, rs
      common/cmtime/ age, time0, lastmd
      common/cdiffu/ idiffs, itbdf1, rhob, rmxdif, tlfdif, dtldif,
     *  ctdfmx, rctdif, rctdf1
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev1,exhm
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
c  Factor in Richard et al. expression
c
      data aldfsc /1.d-2/
c
      setder=.not.noder
c
      cfact=4.d0*pi*aldfsc*car*clight/9.d0
      akk=10.d0**ak
      fact=cfact*t**3/(cp(1)*akk*rho(1)*rho(1))
      difsc(1)=fact*(grad(1)-dad(1))/(dad(1)+dtxh-grad(1))
c
      dmdr=4d0*pi*rho(1)*r*r
      difscx(1)=dmdr*dmdr*difsc(1)
c..      write(6,*) 'SEMIC. R',t,cp(1),akk,rho(1),
c..     *  fact,grad(1),dad(1),dtxh,dmdr,difscx(1)
c
c  test for setting incomplete derivatives (considering only variations
c  in grad)
c
      if(setder) then
	fscgrd=1.d0/(grad(1)-dad(1))+1.d0/(dad(1)+dtxh-grad(1))
	do i=2,7
	  difscx(i)=difscx(1)*fscgrd*grad(i)
	  difsc(i)=difsc(1)*fscgrd*grad(i)
        end do
      end if
      return
      end
