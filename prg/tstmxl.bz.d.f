      subroutine tstmxl(x,y,n)
c
c  Routine to test derivatives from s/r mixlng (for now,
c  in the version calling MJM routines) at the point x=x(n),
c  y=y(.,n)
c
c  Original version: 6/5/03
c
      implicit double precision (a-h,o-z)
      logical nosd,notd,conv
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
c
      parameter(idr1mx = nspcmx+3, idermx = ((nspcmx+3)*(nspcmx+4))/2,
     *  ivrmx1 = ivarmx+1, nspcm2=2*nspcmx, ivrmx4=4*ivarmx,
     *  iyfdmx = 1+ivarmx*(3+2*ivarmx))
c
      dimension y(*),ddac(5),f(5),y0(5),ddac0(5),dacst(20),
     *  tsmjmst(6,20),tsmjm0(6,6)
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev, syear
      common/ln10/ amm
      common/rhcn/ a1,a2,a3,a4,zdummy,nvar,ifwrt,irhtst,inentr,
     *  ii1,ii2,ii3,icomp
      common/crhsvr/ qx, radius
      common/heavy/ zatmos, zhc, zh(1)
      common/prtvar/ rhl,ak,akt,akp,akx,epsnuc(idr1mx),dr,dac,conv,
     *  convos
      common/caddvr/ addvar(5,nnmax)
      common/cnvout/ ddrad,rrcnv,aacnv,ddacad,amach, 
     *  drr(5), ddr(5), da(5),pturb,ycnv,dyda
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,jnc,frcf(6),frcl(6)
      common/eqstcl/ iradp,mode,dtest,skipt,noder
      common/eqprcl/ dmeqpr(12),idiag,iscskp
      common/eqstd/ xi(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     .  dlt(4),gm1,tprh,trhp,rhxp
      common/opcdat/ akk(4)
      common/opcxdr/ akxa
       common/ctsmjm/ tsmjm(6,6)
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3,idgmxc,
     *  idgsdt,idgn75,idgn76,idgn77
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      data init /0/
c
c
      if(istdpr.gt.0) write(istdpr,*) 'Enter tstmxl'
c
c  initialize turbulent pressure to zero
c
      pturb=0
      ptpg=0
c
      qx=10.d0**x
c
      eps=1.d-5
      ichange=0
      jchange=0
      kchange=5
      do i=1,kchange
	y0(i)=y(i)
      end do
      amms=amsun*am
      amass=amms*qx
c
c  start loop over changes to y
c
   10 jchange=jchange+1
      do 12 i=1,kchange
   12 y(i)=y0(i)
      if(jchange.gt.1) then
	ichange=jchange/2
	if(mod(jchange,2).eq.0) then
          y(ichange)=y(ichange)-eps
        else
          y(ichange)=y(ichange)+eps
        end if
      end if
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'jchange, ichange, y:',jchange,ichange,
     *  (y(i),i=1,kchange)
c
c  end setting new value
c
c  calculate log(rho) and log(kappa)
c
      fl=y(2)
      tl=y(3)
      rlred=y(1)
      radius=1.d11*10.d0**y(1)
      allred=y(4)
      xh=y(5)
      yh=1-xh-zh(n)
      zhh=zh(n)
c
      nosd=.false.
      notd=.true.
c
      call eqstf(fl,tl,xh,yh,zhh,nosd,notd)
      pl=log10(pt(1))
c
c  store log p in addvar
c
      addvar(1,n)=pl
c
      p=pt(1)
      t=10.d0**tl
      r=1.d11*1.d1**rlred
      rhl=log10(rho(1))
      call opact(rhl,tl,xh,zhh,ak,rkr,rkt,rkx)
c
      akf=rkr*rho(2)
      akt=rkt+rkr*rho(3)
      akx=akxa+rkr*rho(4)
c
      tgrfct=1.d0
      alshft=0.d0
c
c  set the right hand sides
c
      f(1)=a1*10.d0**(x-3.d0*rlred-rhl)
      f(2)=-a2*10.d0**(2.d0*x-4.d0*rlred-pl)
      f(3)=-a3*tgrfct*10.d0**(ak+x-4.d0*(tl+rlred))*
     *        ((10.d0**allred)-alshft)
c
c  test for convection
c
      dr=f(3)/f(2)
      dac=dr
      ddrad=dr-dad(1)
c
c  pressure scale height
c
      hhp=-r*f(1)/f(2)
c
      conv=dr.gt.dad(1)
      ddacad=dr-dad(1)
      ddarad=dr-dad(1)
c
      if(.not.conv) then
	stop 'tstmxl called in radiative region'
      end if
c
c  store opacity in array for mixlng
c
      akk(1)=10.d0**ak
      akk(2)=akf
      akk(3)=akt
      akk(4)=akx
      call mixlng(t,hhp,dr,dac,ddac,ptpg)
c
      dacst(jchange)=dac
      do k=1,6
	tsmjmst(k,jchange)=tsmjm(k,1)
      end do
      if(jchange.eq.1) then
	do k=1,kchange
          ddac0(k)=ddac(k)
	  do i=1,6
	    tsmjm0(i,k)=tsmjm(i,k+1)
          end do
        end do
      end if
c
      if(jchange.lt.2*kchange+1) go to 10
c
c  output analytical and numerical derivatives
c
      if(istdpr.gt.0) then
        write(istdpr,110) n, x, y0
        write(istdpr,120) dacst(1),eps
        do 30 k=1,kchange
        ddack=(dacst(2*k+1)-dacst(2*k))/(2.d0*eps)
   30   write(istdpr,130) k, ddac0(k),ddack,ddac0(k)-ddack
c
c  derivatives of auxiliary variables from tsmjm
c
        do 40 i=1,6
        write(istdpr,140) i,tsmjmst(i,1),eps
        do 40 k=1,kchange
        dtsmjmk=(tsmjmst(i,2*k+1)-tsmjmst(i,2*k))/(2.d0*eps)
   40   write(istdpr,130) k, tsmjm0(i,k),dtsmjmk,tsmjm0(i,k)-dtsmjmk
      end if
c
      return
  110 format(//' Testing s/r mixlng at n =',i5,'  x =',1pe13.5/
     *  ' y =',1p5e13.5)
  120 format(/' Derivatives of dac, with central value =',1pe13.5/
     *  ' k, analytical deriv., numerical deriv. with eps =',
     *  1pe13.5)
  130 format(i3,1p3e15.7)
  140 format(/' Derivatives of tsmjm(',i1,'), with central value =',
     *  1pe13.5/
     *  ' k, analytical deriv., numerical deriv. with eps =',1pe13.5)
      end
