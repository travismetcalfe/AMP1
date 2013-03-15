      subroutine rnrhec(fl, tl, x, y, z, nosd)
c
c  Set relevant reaction rates for He burning.
c
c  Data from Angulo et al. (1999; Nucl. Phys. A656, 3 -- 183).
c
c  Results are stored in alhe(k,i), k = 1, ..., 10, where alhe(1,i)
c  gives the reaction rate and k = 2, ..., 10 give the derivatives
c  of log10(alhe(1,i)), stored in the usual way. alhe is transmitted
c  in common/rnrhed/.
c
c  The reaction rates are stored in the following order
c
c    alhe(1,.): 4He(4He 4He,gamma)12C (this sets N_A^2 <sigma v>)
c    alhe(2,.): 12C(4He, gamma)16O    (this sets N_A <sigma v>)
c
c  For diagnostic purposes, the individual contributions to alhe(2,1),
c  (but not the derivatives) are stored in
c
c  common/rnrccd/ alcc(4),
c
c  in the order
c
c    alcc(1): E1 contribution
c    alcc(2): E2 contribution
c    alcc(3): res contribution
c    alcc(4): total reaction rate
c
c  Note: normalization by Angulo et al. is consistent with 
c  Fowler et al. (1975) and hence with usage in s/r rnrate.
c
c  Note: so far second derivatives are not set
c
c  For internal calculations, set derivatives of the various
c  contributions to the reaction rates, rather than their logarithms
c
c  Original version: 18/7/02
c
      implicit double precision (a-h, o-z)
      logical nosd, norche
      include 'engenr.bz.d.incl'
c..      dimension alhe1(10),alhe2(10),alcc1(10),alcc2(10),alcc3(10)
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
      common/eqstd/ xii(4),ane(10),rho(10)
      common/rnrhed/ alhe(10,10), norche
      common/rnrccd/ alcc(4)
      common/ln10/ amm,amm2
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
C
      fxp(x)=exp(min(200.d0,max(-200.d0,x)))
c
c  lower limit for temperature (note that a more sophisticated 
c  condition may be required)
c
      tlowlm = 60.d6
c
      do 10 j=1,2
      do 10 i=1,10
   10 alhe(j,i)=0
c
      t=10.d0**tl
c
c  test for no reactions
c
      if(t.lt.tlowlm) then
	norche = .true.
	return
      else
	norche = .false.
      end if

c..   call zero(alhe1,10)
c..   call zero(alhe2,10)
c..   call zero(alcc1,10)
c..   call zero(alcc2,10)
c..   call zero(alcc3,10)
c
      t9=t/1.d9
      t13=t9**(1.d0/3.d0)
      ti13=1.d0/t13
      ti23=ti13*ti13
      ti32=t9**(-1.5d0)
c
c  Pseudorate for 4He(4He,gamma)...
c
c..   alhe1(1)=2.43d9*ti23*fxp(-13.490*ti13-t9*t9/0.0225)*(1+74.5d0*t9)
c..  *         +6.09d5*ti32*fxp(-1.054d0/t9)
      flhe11=2.43d9*ti23*fxp(-13.490*ti13-t9*t9/0.0225)
      flhe12=6.09d5*ti32*fxp(-1.054d0/t9)
      alhe1=flhe11*(1+74.5d0*t9)+flhe12
      dalhe1=flhe11*((-2.d0/3.d0+(13.490/3.d0)*ti13-2*t9*t9/0.0225)*
     *       (1+74.5d0*t9)
     *        + 74.5d0*t9)
     *      +flhe12*(-1.5d0+1.054d0/t9)
c..      alhe1=1
c..      dalhe1=0
c
c  Pseudorate for 8Be(4He,gamma)...
c
c..   alhe2(1)=2.76d7*ti23*fxp(-23.570d0*ti13 -t9*t9/0.16d0)*
c..  *         (1+5.47d0*t9+326d0*t9*t9)
c..  *         +130.7d0*ti32*fxp(-3.38/t9)+2.51d4*ti32*fxp(-20.307/t9)
      flhe21=2.76d7*ti23*fxp(-23.570d0*ti13 -t9*t9/0.16d0)
      flhe22=130.7d0*ti32*fxp(-3.38/t9)
      flhe23=2.51d4*ti32*fxp(-20.307/t9)
      alhe2=flhe21*(1+5.47d0*t9+326d0*t9*t9)+flhe22+flhe23
      dalhe2=flhe21*((-(2.d0/3.d0)+(23.570d0/3.d0)*ti13-t9*t9/0.08d0)*
     *    	 (1+5.47d0*t9+326d0*t9*t9)+5.47d0*t9+652d0*t9*t9)
     *	    +flhe22*(-1.5d0+3.38d0/t9)
     *      +flhe23*(-1.5d0+20.307d0/t9)
c..      alhe2=1
c..      dalhe2=0
c
c Rate for 4He(4He 4He,gamma)12C
c
      if(t9.le.0.03d0) then
	fct=3.07d-16*(1-29.1d0*t9+1.308d3*t9*t9)
	dfct=-3.07d-16*t9*(29.1d0-2.616d3*t9)
      else
c
c  Note: coefficient in following expression changed from 3.44d-16 to
c  3.4685d-16, to ensure continuity.
c
	fct=3.4685d-16*(1+0.0158d0*t9**(-0.65d0))
	dfct=-3.4685d-16*0.0158d0*0.65d0*t9**(-0.65d0)
      end if
c
      alhe(1,1)=alhe1*alhe2*fct
      alhe(1,3)=dalhe1/alhe1+dalhe2/alhe2+dfct/fct
c
c  12C(4He, gamma)16O
c
c..   alcc1=6.66d7/(t9*t9)*fxp(-32.123*ti13-(t9/4.6)**2)*
c..  *          (1+t9*(2.54+t9*(1.04-0.226*t9)))+
c..  *          1.39d3*ti32*fxp(-28.93/t9)
      flcc11=6.66d7/(t9*t9)*fxp(-32.123*ti13-(t9/4.6)**2)
      flcc12=1.39d3*ti32*fxp(-28.93/t9)
      alcc1=flcc11*(1+t9*(2.54+t9*(1.04-0.226*t9)))+flcc12
      dalcc1=flcc11*((-2.d0+(32.123d0/3.d0)*ti13-2*(t9/4.6)**2)*
     *             (1+t9*(2.54+t9*(1.04-0.226*t9)))
     *               +t9*(2.54+t9*(2.08-0.678*t9)))
     *      +flcc12*(-1.5d0+28.93d0/t9)
c..   alcc2=6.56d7/(t9*t9)*fxp(-32.123*ti13-(t9/1.3)**2)*
c..  *          (1+t9*(9.23-t9*(13.7-7.4*t9)))
      flcc21=6.56d7/(t9*t9)*fxp(-32.123*ti13-(t9/1.3)**2)
      alcc2=flcc21*(1+t9*(9.23-t9*(13.7-7.4*t9)))
      dalcc2=flcc21*((-2.d0+(32.123d0/3.d0)*ti13-2*(t9/1.3)**2)*
     *                (1+t9*(9.23-t9*(13.7- 7.4*t9)))+
     *                   t9*(9.23-t9*(27.4-22.2*t9)))
      alcc3=19.2d0*t9*t9*fxp(-26.9d0/t9)
      dalcc3=alcc3*(2.d0+26.9d0/t9)
c
      alhe(2,1)=alcc1+alcc2+alcc3
c..      write(6,'(a,1p4e13.5)') 'alhe(2,1):',alcc1,alcc2,alcc3,alhe(2,1)
      alhe(2,3)=(dalcc1+dalcc2+dalcc3)/alhe(2,1)
c
      return
      end
