      subroutine eqstf(rl,tl,x,y,z,nosd,notd)
c
c
c
c========= emulator for Livermore tables ============
c
c
c  Double precision version interface, Livermore routines single precision
c  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  Dated 28/11/94
c
c  Modified 9/9/95: When temperature is below 5.e3, call esac at
c  5.e3 and scale p and H by temperature ratio
c  Add derivatives wrt X.
c
c  Modified 11/9/94, to set centralized first derivatives, second
c  derivatives, when nosd is false.
c
c  Modified 24/8/96, to call double-precision version of esac
c
c  Modified 15/3/05, replacing ht(2) by htdrh in calculation of htdt
c
      implicit double precision (a-h, o-z)
      logical nosd,notd
c
      dimension ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),dlt(4)
      dimension der(10,7)
c
      common/ln10/ amm,amm2,amm3
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev1,exhm
      common/eqscnt/ anh0, anhe0, ihvz, iprrad, ihmin
      common/eqstd/ xii1(4),ane,rho,ht,pt,cp,dad,dlt,gm1,tprh,trhp,
     .  rhxp,gmm1(4)
c
      common/e/esact,eos(10)
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
c..      write(61,'(5f12.7)') rl, tl, x, y, z
c
c  minimum temperature in tables
c
      t6mins=0.005001d0
c
c  steps in derivatives
c
c
      dtl=0.002
      drhl=0.01
      dx=0.005
c
c
c===================================================================
c              call livermore routine
c===================================================================
c
      ztabs=z
      iorder=10
      if(iprrad.eq.0) then
        irad = 0
      else
        irad = 1
      end if
c
c  step over central point and one- or two-sided derivatives
c
      if(nosd) then
	nodir = 4
      else
	nodir = 7
      end if
c
      do 30 ider=1,nodir
      if(ider.eq.1) then
        rh=10.d0**rl
        tt=10.d0**tl
        xs=x
      else if(ider.eq.2) then
        rh=10.d0**(rl+drhl)
        tt=10.d0**tl
        xs=x
      else if(ider.eq.3) then
        rh=10.d0**rl
        tt=10.d0**(tl+dtl)
        xs=x
      else if(ider.eq.4) then
        rh=10.d0**rl
        tt=10.d0**tl
        xs=x+dx
      else if(ider.eq.5) then
        rh=10.d0**(rl-drhl)
        tt=10.d0**tl
        xs=x
      else if(ider.eq.6) then
        rh=10.d0**rl
        tt=10.d0**(tl-dtl)
        xs=x
      else if(ider.eq.7) then
        rh=10.d0**rl
        tt=10.d0**tl
        xs=x-dx
      end if
      t6s=1.e-6*tt
      rs=rh
c
c  test for temperature too low
c
      if(t6s.gt.t6mins) then
        call esac(xs,ztabs,t6s,rs,iorder,irad)
      else
        call esac(xs,ztabs,t6mins,rs,iorder,irad)
	trat=t6s/t6mins
	eos(1)=trat*eos(1)
	eos(2)=trat*eos(2)
	eos(4)=trat*eos(4)
      end if
c
c  test for error
c
      if(eos(1).lt.0) then
	pt(1)=eos(1)
	return
      end if
c
      ptt   = eos(1)*1.d12
      ptdrh = eos(6)
      ptdt  = eos(7)
      htt   = eos(2)*1.d12+ptt/rh
      dadd  = 1./eos(9)
      dltt  = eos(7)/eos(6)
      alpha = 1./eos(6)
c  Note: here assume that energy is still in 1.e12 erg, T is T6
      cpp   = 1.e6*eos(5)+ptt*dltt*dltt/(rh*tt*alpha)
      htdrh = 1.e12*amm*rh*eos(4)+amm*(ptt/rh)*(eos(6)-1)
      htdt  = amm*cpp*tt+htdrh*dltt
c
      if(ider.eq.1) then
        call zero(xii1,90)
c
        rho(1)  = rh
        pt(1)   = ptt
	pt(2)   = ptdrh
	pt(3)   = ptdt
        ht(1)   = htt
	ht(2)   = htdrh
	ht(3)   = htdt
	cp(1)   = cpp
        dad(1)  = dadd
        dlt(1)  = dltt
        gm1     = eos(8)
        gmm1(1) = eos(8)
        tprh    = 1./eos(7)
        trhp    = -1./dlt(1)
c
      end if
      der(1,ider) = htt
      der(2,ider) = log10(ptt)
      der(3,ider) = cpp
      der(4,ider) = dadd
      der(5,ider) = dltt
      der(6,ider) = eos(8)
      der(7,ider) = ptdrh
      der(8,ider) = ptdt
      der(9,ider) = htdrh
      der(10,ider) = htdt
c
   30 continue
c
c  set derivatives
c
      if(nosd) then
c
        do 40 k=1,6
        der(k,2)=(der(k,2)-der(k,1))/drhl
        der(k,3)=(der(k,3)-der(k,1))/dtl
   40   der(k,4)=(der(k,4)-der(k,1))/dx
c
      else
c
c  first set second derivatives wrt X, before der(.,4) is overwritten
c
	htdxx  =(der(1,4)+der(1,7)-2*der(1,1))/(dx*dx)
	ptdxx  =(der(2,4)+der(2,7)-2*der(2,1))/(dx*dx)
c
        do 45 k=1,10
c
        der(k,2)=(der(k,2)-der(k,5))/(2*drhl)
        der(k,3)=(der(k,3)-der(k,6))/(2*dtl)
        der(k,4)=(der(k,4)-der(k,7))/(2*dx)
   45   continue
c
      end if
c
      rho(2)  = 1.d0
c..      write(71,'(2f10.5,1p8e13.5)') 
c..     *  rl, tl, pt(2), pt(3), ht(2), ht(3), der(2,2), 
c..     *  der(2,3), der(1,2), der(1,3)
      do 50 i=2,4
      cp(i)   = der(3,i)
      dad(i)  = der(4,i)
      dlt(i)  = der(5,i)
   50 gmm1(i) = der(6,i)
c
      ht(4)   = der(1,4)
      pt(4)   = der(2,4)
c
      if(.not.nosd) then
	pt(5) = der(7,2)
	pt(6) = der(8,2)
	pt(7) = der(7,4)
	pt(8) = der(8,3)
	pt(9) = der(8,4)
	pt(10) = ptdxx
c
	ht(5) = der(9,2)
	ht(6) = der(10,2)
	ht(7) = der(9,4)
	ht(8) = der(10,3)
	ht(9) = der(10,4)
	ht(10) = htdxx
c
      end if
c
      if(.not.notd) then
c
c  approximate treatment of third derivatives
c
        amu1 = 2*x/ah + 3*y/ahe + (anh0 + 1./az)*z
        b = (2/ah - 3/ahe)
c
        pt(20) = 2*((b/amu1)**3)/amm
c
        ht(17) = amm3*ht(1)
        ht(18) = amm2*ht(1)*b/amu1
c
      end if
c
c  set Ne as for fully ionized gas (only relevant, probably,
c  in electron screening in the core)
c
      ane(1)=av*(1.d0*x/ah+2.d0*y/ahe+z/anh0)
      ane(4)=av*(1.d0/ah-2.d0/ahe)
c
      rhxp=-amm*x*pt(4)/pt(2)
c
c..      write(6,*) 'rl,tl,rho(1),pt(1),cp(1),dad(1),dlt(1),gmm1(1)'
c..      write(6,*) rl,tl,rho(1),pt(1),cp(1),dad(1),dlt(1),gmm1(1)
c
      return
      end
      subroutine seteqs
c
c  dummy subroutine included for compatibility with dog equation  of
c  state programmes.
c
      return
      end
