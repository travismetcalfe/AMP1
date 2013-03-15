c*****************************************************************************
c
	subroutine mxlmjm (t,hp,dr,dac,ddac,ptpg)
c
c  Routine to calculate the actual temperature gradient "dac"
c  and its derivatives "ddac(1-3)" from CM formulations:
c		a-(1991)	b-(1992)	c-(1996)
c  See subroutine "init_conv" for how to choose one of these
c  formulations in particular.
c  The use of this subroutine assumes that "init_conv" was called
c  ONCE at the begining (necessarily before this one is used).
c
c  This version is valid both for the envelope and evolution codes
c  since it incorporates a parameter which establishes which of
c  those is being used and writes the derivatives accordingly.
c
c  Modified 29/7/96 to correct setting of turbulent velocity and
c  turbulent pressure; see note from MJM in f.mjm.960729
c
c  Modified 18/11/98, adding ptpg to the argument list, for case
c  with turbulent pressure
c  (Note: this still need to be made consistent with the computation
c  of turbulent pressure in MJM routine).
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 7/96
c
c  Modified 23/8/96, to interface directly with JC-D code.
c
c  Corrected 6/5/03, to fix several errors (related to log(10.0d0))
c  in the calculation of the derivatives. These have now been confirmed
c  using numerical derivatives. The array tsmjm was introduced for this
c  test.
c
	implicit double precision (a-h,o-z)
	character aconv*1,aconvtype*3
	parameter (error=1.0d-14,imax=200)
	dimension ddac(1),is(5),dt(5),dlxl(5),ddr(5),dlb(5)
	common/ln10/ amm
	common/mxlcn/ cc1,cc2,etac,phc,alfa,tprfct,iconcs,imxlng,iturpr
	common/eqstd/ dum(14),rho,drho(19),ht(20),p,dp(19),cp,dcp(3),
     *		dad,ddad(3),dlt,ddlt(3),gm1
	common/rhcn/ a1,a2,a3,a4,zheavy,nvar,ifwrt,irhtst
	common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
	common/opcdat/ akk,dak(3)
	common/cnvout/ ddrad,rr,a,ddacad,bamach,adrr(5),addr(5),ada(5)
	common/ctsmjm/ tsmjm(6,6)
c
	common/crhsvr/qx,x
	common/totmss/xm,xr
c
c - Convection specific "common" input:
c
c  The variable "ncodeversion" indicates which version of the
c  code is being used:	0 - envelope code
c			1 - evolution code
c  While the vector "is" gives the position of the derivatives for output
	common /cmconva/ncodeversion,is
c
c  The options passed on by "cmconv0" are,
c		aconv='h' - mixing length is proportional to Hp
c		aconv='z' - mixing length is given by depth (CM 1991)
c       	aconv='f' - mixing length varies according to f (CM 1992)
c  and
c		aconvtype='new' - new CM theory (1996)
c		aconvtype='old' - Old CM theory (1991)
c		aconvtype='par' - parametrization by MCDT (1996)
c
	common /cmconv0/aconv,aconvtype
c
c  Block "cmconv1" is only needed if aconv='z'
	common /cmconv1/xlgs,fstab
c
c  This block provides the parameters of the formulation as given by the
c  authors:
	common /cmconv2/xa1,xa2,em,en,ep
	common /cmconv3/zc,zd,ze,zf,zp,zq,zr,zt
	common /cmconv4/xk0,zg,zh,zi,zj,zk,zl,zo,zu,zv
c
c - Common defining standard input and output
c
	common/cstdio/ istdin, istdou, istdpr, istder
c
	save
c
	data etacp,phcp /-1.,-1./
	data dt/0.0d0,1.0d0,0.0d0,0.0d0,0.0d0/
c
c  Testing for the initialization of the parameters in "init_conv":
	data initconv/0/
	idgmjm=0
	if(idgmjm.eq.2.and.istdpr.gt.0) then
          write(istdpr,*) 'Enter mixlng with t,hp,dr =',t,hp,dr
          write(istdpr,*) 'qx, x, xm, xr =', qx,x,xm,xr
          write(istdpr,*) 'rho, akk, p, cp =', rho, akk, p, cp
          write(istdpr,*) 'dak =',dak
          write(istdpr,*) 'dp =',(dp(i),i=1,3)
          write(istdpr,*) 'drho =',(dp(i),i=1,3)
          write(istdpr,*) 'dcp =',(dcp(i),i=1,3)
          write(istdpr,*) 'ddad =',(ddad(i),i=1,3)
          write(istdpr,*) 'ddlt =',(ddlt(i),i=1,3)
	end if
	if (initconv.eq.0) then
		initconv=1
		call init_conv
		if(istdpr.gt.0) write(istdpr,*) 
     *            ' After initialization of MJM initconv'
	endif
c
	do i=1,6
	  do j=1,6
	    tsmjm(i,j)=0
          end do
        end do
c
	if(idgrhs.eq.2.and.istdpr.gt.0) 
     *    write(istdpr,*) 'cc1,cc2,etac,phc,alfa', 
     *    cc1,cc2,etac,phc,alfa
	if(idgrhs.eq.2.and.istdpr.gt.0) write(istdpr,*) 
     *		' alfa,  hp, aml, rho, cp, ak, g, dlt =',
     *		alfa,  hp, aml, rho, cp, ak, g, dlt
c-------
c
c  delta r - delta ad
	ddrad=dr-dad
	if(ddrad.lt.0) ddrad=1.e-10
c
	temp1=cc1/(9.0d0*cc2)
	temp2=dlt*rho/(2.0d0*p)
	temp3=cp*akk*rho*rho*qx*hp*hp/(x*x*t*t*t)
	bcoef=temp1*temp1*temp2*temp3*temp3
	y0=4.0d0*bcoef
c
c  Here "xu" should be the top of the convection zone and not "1*R"!
	xu=1.0d0*xr
	z=xu-x
c
	b1=xa1*y0**em
	b2=y0*xa2
	if(idgmjm.eq.2.and.istdpr.gt.0) 
     *    write(istdpr,*) 'y0, xa2, b2', y0, xa2, b2
	ys0=log(40.5d0*y0)
	if (aconv.eq.'f') then
		xl=z/hp
	else
		xl=alfa
	endif
c
	iter=0
	s1=1.0d-30
	s2=xl**4*ddrad
 10	iter=iter+1
	g1=gs_conv(alfa,s1,ddrad,xu,x,hp,b1,b2,ys0)
	g2=gs_conv(alfa,s2,ddrad,xu,x,hp,b1,b2,ys0)
	if (g1*g2.gt.0.) then
		if (iter.gt.30) then
			write (*,*) 'WARNING: No root yet in CONV !'
			write (*,*) x,s1,s2
			pause
			iter=0
		endif
		s1=s1/10.0d0
		s2=s2*10.0d0
		goto 10
	endif
c
	iter=0
 20	iter=iter+1
	sm=0.5d0*(s1+s2)
	gm=gs_conv(alfa,sm,ddrad,xu,x,hp,b1,b2,ys0)
	if (gm*g1.lt.0.0d0) then
		s2=sm
		g2=gm
	else
		s1=sm
		g1=gm
	endif
	err=abs((s2-s1)/(s1+s2))
	if (iter.gt.imax) then
		write (*,*) 'ERROR: Too many iterations in CONV !'
		write (*,*) s1,sm,s2
		pause
		iter=0
	endif
	if (err.gt.error) goto 20
	s=0.5d0*(s1+s2)
c-------
c  Now the results we need:
c
	if (aconv.eq.'f') then
		ys=log(40.5d0*y0*s)
		fss=fs_conv(ys)
		xl=(z/hp)*exp(-alfa*log(fss))
	else if (aconv.eq.'z') then
		xl=((1.0d0+alfa*1.0d-4)*xu-x)/hp
	else
		xl=alfa
	endif
c  grad - grad ad
	ddacad=s/xl**4
	dac=dad+ddacad
   	if(idgmjm.eq.2.and.istdpr.gt.0) 
     *    write(istdpr,*) ' ddrad, ddacad, dac =',
     *		ddrad, ddacad, dac
c
c  mach number (amach=vt/c) and turbulent pressure ratio (ptp=Pt/P):
	if (aconvtype.eq.'new') then
	        bs=40.5d0*y0*s
		f3=zg*bs*bs/(1.0d0*sqrt(1.0d0+zh*bs*bs))
		zkbszl=zk*bs**zl
		f4=zi+zj*(-1.0d0+zkbszl)/(1.0d0+zkbszl)
		ptfac=(xk0/1.5d0)**3/(40.5d0*y0*xl**2)
		amach=ptfac*f3*f4/gm1
		zvbszl=zv*bs**zl
		f5=zo+zu*(-1.0d0+zvbszl)/(1.0d0+zvbszl)
		ptp=ptfac*f3*f5
	else if (aconvtype.eq.'par') then
c  These relations are a very rough approximation;
	        amach=sqrt(2.0d0*s/gm1)/xl
		ptp=s/xl
	else
		amach=sqrt(2.0d0*s/gm1)/xl
		ptp=2.0d0*s/(5.0d0*xl)
	endif
c  necessary in an output common block, but not used here,
c are "curly r" and the quantity "a":
	g=p*(1.d0+ptpg)/(hp*rho)
	ak=cc2*t*t*t/(akk*rho)
c  the mixing length ("xl" has been defined above!)
	aml=xl*hp
c  curly r
	rr=aml*aml*rho*cp/ak
	rr=(rr/hp)*rr*g*dlt
	if(idgrhs.eq.2.and.istdpr.gt.0) write(istdpr,*) ' rr =',rr
c  set a
	a=2.0d0/(sqrt(rr*ddrad)*etac)
	if(idgrhs.eq.2.and.istdpr.gt.0) write(istdpr,*) ' a =',a
c
c  diagnostics if idgrhs = 1
c
	if(idgrhs.eq.1.and.istdpr.gt.0) write(istdpr,59901) 
     *		t,hp,dr,p,rho,b,y0,s,ddacad,dac
59901	format(' t,hp,dr,p,rho =',1p5e15.7/
     *		' b,y0,s,ddacad,dac=',1p5e15.7)
c
c-------
c  Derivatives:
c
	temp0=1.0d0+b2*s
	temp1=en*ep*b2*s*temp0**(en-1.0d0)/(temp0**en-1.0d0)
	temp2=xl**4*ddrad-s
        temp3=s/temp2+em+1.0d0+temp1
	if(idgmjm.eq.2.and.istdpr.gt.0) then
	  write(istdpr,*) 'en, ep, b2, s, temp0 =',en, ep, b2, s, temp0
          write(istdpr,*) 'temp1, temp2 =',temp1, temp2
        end if
c
	do 102 jk=1,5
c  The derivatives of "xl" are zero only for aconv='h', otherwise they
c  are not. So the expression for the derivative of $log(s)$ should
c  include this extra term. But we shall leave them equal to zero for
c  the time being! (To messy otherwise!)
		dlxl(jk)=0.0d0
 102		continue
c
c  Derivatives for: (log f, log T, X);
	do 101 jk=1,3
		ik=is(jk)
		ddr(jk)=dr*log(10.0d0)*(dak(jk)+dp(jk)-4.0d0*dt(jk))
		dlcp=(dcp(jk)/cp)/log(10.0d0)
		dldlt=(ddlt(jk)/dlt)/log(10.0d0)
c..		dlb(jk)=dldlt+2.0d0*dlcp+2.0d0*dak(jk)+drho(jk)+
		dlb(jk)=dldlt+2.0d0*dlcp+2.0d0*dak(jk)+drho(jk)
     *			-6.0d0*dt(jk)+3.0d0*dp(jk)
		ddac(ik)=((xl**4/temp2)*(4.0d0*ddrad*dlxl(jk)+ddr(jk)
     *			-ddad(jk))-log(10.0d0)*(em+temp1)*dlb(jk))/temp3
		ddac(ik)=ddad(jk)+s*(ddac(ik)-4.d0*dlxl(jk))/xl**4
c
c  store in tsmjm
c
	       tsmjm(1,1)=dr
	       tsmjm(1,ik+1)=ddr(jk)
	       tsmjm(2,1)=log10(cp)
	       tsmjm(2,ik+1)=dlcp
	       tsmjm(3,1)=log10(dlt)
	       tsmjm(3,ik+1)=dldlt
	       tsmjm(4,1)=log10(bcoef)
	       tsmjm(4,ik+1)=dlb(jk)
	       tsmjm(5,1)=dad
	       tsmjm(5,ik+1)=ddad(jk)
	       tsmjm(6,1)=xl
	       tsmjm(6,ik+1)=dlxl(jk)
 101		continue
c
	jk=4
	ik=is(jk)
c  Derivative with respect to radius;
	ddr(jk)=0.0d0
	dlb(jk)=4.0d0
	ddac(ik)=((xl**4/temp2)*(4.0d0*ddrad*dlxl(jk)+ddr(jk))
     *		-log(10.0d0)*(em+temp1)*dlb(jk))/temp3
	ddac(ik)=s*(ddac(ik)-4.d0*dlxl(jk))/xl**4
	tsmjm(4,ik+1)=dlb(jk)
c
	jk=5
	ik=is(jk)
	if (ncodeversion.ne.0) then
c  Derivative with respect to luminosity;
		ddr(jk)=dr*log(10.0d0)
		dlb(jk)=0.0d0
		ddac(ik)=((xl**4/temp2)*(4.0d0*ddrad*dlxl(jk)
     *			+ddr(jk))-log(10.0d0)*(em+temp1)*dlb(jk))/temp3
		ddac(ik)=s*(ddac(ik)-4.d0*dlxl(jk))/xl**4
	        tsmjm(1,ik+1)=ddr(jk)
	else
c  Derivative with respect to mass;
		ddr(jk)=dr*(-1.0d0)*log(10.0d0)
		dlb(jk)=-2.0d0
		ddac(ik)=((xl**4/temp2)*(4.0d0*ddrad*dlxl(jk)
     *			+ddr(jk))-log(10.0d0)*(em+temp1)*dlb(jk))/temp3
		ddac(ik)=s*(ddac(ik)-4.d0*dlxl(jk))/xl**4
	endif
c
      if(idgrhs.eq.1.and.istdpr.gt.0) write(istdpr,80091) dac
      if(idgmjm.eq.1.and.istdpr.gt.0) write(istdpr,'(a,1p4e14.6)') 
     *  'dac, ddac =',dac, (ddac(i),i=1,3)
80091 format(' dac =',1pe15.7)
c
      return
      end
c
c**************************************************************************
c
	double precision function gs_conv (xa,s,dra,xu,x,hp,b1,b2,ys0)
c This subroutine defines the fucntion that allow us to find the reduced
c convective efficiency "s", corresponding to the zero of "gs_conv(s)=0".
c
	implicit double precision (b-h,o-z)
	character aconv*1,aconvtype*3
	common /cmconv0/aconv,aconvtype
	common /cmconv2/xa1,xa2,em,en,ep
	common /cmconv3/zc,zd,ze,zf,zp,zq,zr,zt
	common /cmconv4/xk0,zg,zh,zi,zj,zk,zl,zo,zu,zv
c
c Mixing length over pressure scale height:
c
	if (aconv.eq.'f') then
		ys=ys0+log(s)
		fss=fs_conv(ys)
		xl=((xu-x)/hp)*exp(-xa*log(fss))
	else if (aconv.eq.'z') then
		xl=((1.0d0+xa*1.0d-4)*xu-x)/hp
	else
		xl=xa
	endif
c
c Correction term ofr the new theory (1996):
c
	if (aconvtype.eq.'new') then
	        bs=40.5d0*(b2/xa2)*s
		f2=1.0d0+zc*bs**zp/(1.0d0+zd*bs**zq)+
     *			ze*bs**zr/(1.0d0+zf*bs**zt)
		b11=b1*f2
	else
		b11=b1
	endif
c
c Actual function (see Monteiro's thesis):
c
	gs_conv=1.0d0-s*(1.0d0+b11*s**em*((1.0d0+b2*s)**en-1.0d0)**ep)/
     *		(xl**4*dra)
c
	return
	end
c
c****************************************************************************
c
	double precision function fs_conv(x)
c This subroutine provides the values of "f(s)" from the table given
c in CM (1992) which enters the definition of the mixing length.
c It is only used if "aconv='f'".
c
	implicit double precision (b-h,o-z)
	parameter (np=1,na=2*np,nf=8,nint=10)
	dimension xlgs(nf),fstab(nf),xa(nint),ya(nint)
	common /cmconv1/xlgs,fstab
c
c Location of the point:
c
	if (x.lt.xlgs(1)) then
		i1=0
	else if (x.gt.xlgs(nf)) then
		i1=8
	else
		i1=1
		i2=8
 10		if (i2-i1.gt.1) then
			im=int((i1+i2)/2)
			if (x.gt.xlgs(im)) then
				i1=im
			else
				i2=im
			endif
			goto 10
		endif
	endif
c
c Interpolation (linear if np=1):
c
	j=min(max(i1+1-np,1),nf+1-2*np)
	do 50 k=1,na
		xa(k)=xlgs(k+j-1)
		ya(k)=fstab(k+j-1)
 50		continue
	call interp_conv (xa,ya,na,x,y,dy)
	fs_conv=y
c
	return
	end
c
c***************************************************************************
c
	subroutine interp_conv (xi,yi,n,x0,y,dy)
c Last changed: Jan 95
c Subroutine for polynomial interpolation (degree "n-1") of the
c "n" points "(xi,yi)", at the value "x0". The result is "y".
c The set of points is reduced to the interval [0,1], both in
c "x" and "y" (only important for n>4).
c
c See subroutine POLINT in Numerical Recipies, pag. 82.
c
        implicit double precision (b-h,o-z)
        parameter (np=10)
        dimension xi(np),yi(np),xa(np),ya(np),c(np),d(np)
c
        xmax=0.0d0
        ymax=0.0d0
        do 1 i=1,n
                xmax=max(abs(xi(i)),xmax)
                ymax=max(abs(yi(i)),ymax)
 1              continue
        if (xmax.eq.0.0d0) xmax=1.0d0
        if (ymax.eq.0.0d0) ymax=1.0d0
        do 2 i=1,n
                xa(i)=xi(i)/xmax
                ya(i)=yi(i)/ymax
 2              continue
        x=x0/xmax
c
        ns=1
        dif=abs(x-xa(1))
        do 11 i=1,n
                dift=abs(x-xa(i))
                if (dift.lt.dif) then
                        ns=i
                        dif=dift
                endif
                c(i)=ya(i)
                d(i)=ya(i)
 11             continue
        y=ya(ns)
        ns=ns-1
        do 13 m=1,n-1
                do 12 i=1,n-m
                        h0=xa(i)-x
                        hp=xa(i+m)-x
                        w=c(i+1)-d(i)
                        den=h0-hp
                        if (den.eq.0.0d0) then
			  write (*,*) 'ERROR: Points with same x ',
     *					'in INTERP !'
				stop
			endif
                        den=w/den
                        d(i)=hp*den
                        c(i)=h0*den
 12                     continue
                if (2*ns.lt.n-m) then
                        dy=c(ns+1)
                else
                        dy=d(ns)
                        ns=ns-1
                endif
                y=y+dy
 13             continue
c
        y=y*ymax
        dy=dy*ymax
c
        return
        end
c
c****************************************************************************
c
	subroutine init_conv
c This is the setup routine for the convection subroutine implementing 
c the CM formulations (1991,1992 and 1996) and also the MCDT (1996)
c PTG parametrization.
c
c This subroutine has to be called ONCE in the begining before the
c convection subroutine is ever used.
c
c  In this version, the choice of methods is based on the 
c  variables in common/mxlcn/, which are used to initialize
c  MJM character flags.
c
c In order to define what will be used it is necessary to establish
c the following variables:
c	ncodeversion - which code is being used
c	aconv - which mixing length to use
c	aconvtype - which formulation
c and	(betac,rem) - parameters of the PTG formulation
c		     (necessary for aconvtype='par')
c The default values are: ncodeversion=1, aconv='h', aconvtype='old',
c			  betac=1, em=-1
c
	implicit double precision (b-h,o-z)
	double precision alfa
	character aconv*1,aconvtype*3
	dimension xlgs(8),fstab(8),is(5)
c
c  The options passed on by "cmconv0" define the mixing length to
c  be used and the formulation (see below).
        common/mxlcn/ cc1,cc2,etac,phc,alfa,tprfct,iconcs,imxlng,iturpr
	common /cmconva/ncodeversion,is
	common /cmconv0/aconv,aconvtype
	common /cmconv1/xlgs,fstab
	common /cmconv2/xa1,xa2,em,en,ep
	common /cmconv3/zc,zd,ze,zf,zp,zq,zr,zt
	common /cmconv4/xk0,zg,zh,zi,zj,zk,zl,zo,zu,zv
	common /cmconv5/rem,betac
c
c  common defining standard input and output
c
        common/cstdio/ istdin, istdou, istdpr, istder
c
c--The possible types of mixing length available are:
c	aconv='h' - mixing length is proportional to Hp
c	aconv='z' - mixing length is given by depth
c       aconv='f' - mixing length varies according to f (CM 1992)
c
c  Default value;
	data aconv/'h'/
c
c--The formulations available are:
c	aconvtype='new' - New CM theory (1996)
c	aconvtype='old' - Old CM theory (1991)
c	aconvtype='par' - Parametrized temperature gradient
c			  as given by MCDT, A&A (1996)
c
c  Default value;
	data aconvtype/'old'/
c
c--Default values for the case of having "aconvtype='par'";
	data betac/1.0d0/,rem/-1.0d0/
c
c--The code version being used is indicated in the variable
c		ncodeversion =	0  for the envelope code 
c				1  for the evolution code (default value)
	data ncodeversion/1/
c
c  The vector "is" defines the ordering of the derivatives in the
c  output of subroutine "mixlng", depending on what code is beign used.
	if (ncodeversion.ne.0) then
		is(1)=2
		is(2)=3
		is(3)=5
		is(4)=1
		is(5)=4
	else
		is(1)=1
		is(2)=4
		is(3)=5
		is(4)=2
		is(5)=3
	endif
c
c-------
c..        write(6,*) 'Enter mxlmjm'
c..	write(6,*) 'cc1,cc2,etac,phc,alfa,tprfct,iconcs,imxlng,iturpr:'
c..	write(6,*) cc1,cc2,etac,phc,alfa,tprfct,iconcs,imxlng,iturpr
c
c  set character flags based on values in common/mxlcn/
c
	if(iconcs.eq.1) then
	  aconvtype='par'
	else if(iconcs.eq.2) then
	  aconvtype='old'
	else if(iconcs.eq.3) then
	  aconvtype='new'
	else
	  if(istdpr.gt.0) write(istdpr,120) iconcs
	  stop 'init_conv 1'
	end if
c
	if(imxlng.eq.0) then
	  aconv='h'
	else if(imxlng.eq.1) then
	  aconv='z'
	else if(imxlng.eq.2) then
	  aconv='f'
	else
	  if(istdpr.gt.0) write(istdpr,130) imxlng
	  stop 'init_conv 2'
	end if
c
	if(iconcs.eq.1) then
	  betac=phc
	  rem=etac
	end if
c
	if(istdpr.gt.0) 
     *    write(istdpr,*) 'Enter init_conv with aconv, aconvtype =',
     *		aconv, '  ', aconvtype
c
	if (aconvtype.eq.'new') then
c - Parameters as given in Canuto (ApJ, 1996):
		if(istdpr.gt.0) write(istdpr,*) 
     *            'Initialize for Canuto 1996'
		xk0=1.8d0
c   Basic (F1):
		xaa=10.8654d0
		xbb=0.00489073d0
		em=0.149888d0
c		en=0.189238d0
		en=0.5d0*(1.0d0-2.0d0*em)/(2.0d0-em)
c		ep=1.85011d0
		ep=2.0d0-em
c
		xa1=xaa*40.5d0**em*(xk0/1.5d0)**3
		xa2=xbb*40.5d0
c   Corrections:
c          F2;
		zc=0.0108071d0
		zd=0.00301208d0
		ze=0.000334441d0
		zf=0.000125d0
		zp=0.72d0
		zq=0.92d0
		zr=1.2d0
		zt=1.5d0
c   turbulent pressure and velocity:
c          F3;
                zg=0.00101392d0
                zg=zg*(xk0/1.5d0)**3
                zh=0.000017848d0
c          F4;
                zi=6.39899d0
                zj=2.246815d0
                zk=0.000777055d0
                zl=0.868589d0
c          F5;
		zo=1.49168d0
		zu=0.45185d0
		zv=0.00111378d0
c
	else if (aconvtype.eq.'par') then
c - This is the parametrization in terms of $\beta_c$ and $m$ as given
c   in MCDT - A&A (1996):
		if(istdpr.gt.0) write(istdpr,*) 
     *            'Initialize for MCDT 1996'
		em=rem
        	en=0.5d0*(1.0d0-2.0d0*em)/(2.0d0-em)
        	ep=2.0d0-em
c
		xa1=9.0d0/8.0d0
		xa1=xa1*(betac**(5.0d0-4.0d0*em)*((1.0d0-2.0d0*em)/
     *			(2.0d0-em))**(2.0d0*em*em-5.0d0*em+2.0d0)/
     *			2.0d0**(2.0d0*em*em+em-1.0d0))**(1.0d0/3.0d0)
		xa2=(((2.0d0-em)/(1.0d0-2.0d0*em))**(2.0d0-em)/
     *			(betac**2*2.0d0**(em+1.0d0)))**(2.0d0/3.0d0)
c   If "betac=1" and "em=-1" (the default values) then
c   "xa1=9/8" and "xa2=1" corresponding to the MLT formulation!
c
	else
c - Parameters as given in Canuto & Mazzitelli (A&A, 1991):
		if(istdpr.gt.0) write(istdpr,*) 
     *            'Initialize for CM 1991'
        	xa1=24.868d0
        	xa2=0.097666d0
        	em=0.14972d0
c       	en=0.18931d0
        	en=0.5d0*(1.0d0-2.0d0*em)/(2.0d0-em)
c       	ep=1.8503d0
        	ep=2.0d0-em
	endif
c
	if (aconv.eq.'f') then
c - Function defined in Canuto & Mazzitelli (A&A, 1992) as a table
c   which enters the definition of the mixing length:
		do 10 i=1,8
			xlgs(i)=2.0d0*dfloat(i)
 10			continue
		fstab(1)=1.1375d0
		fstab(2)=0.9260d0
		fstab(3)=0.7325d0
		fstab(4)=0.6296d0
		fstab(5)=0.5818d0
		fstab(6)=0.5581d0
		fstab(7)=0.5428d0
		fstab(8)=0.5492d0
	endif
c
	return
  120   format(//' ***** Error in s/r init_conv. iconcs =',i4,
     *  ' not allowed')
  130   format(//' ***** Error in s/r init_conv. imxlng =',i4,
     *  ' not allowed')
	end
c
