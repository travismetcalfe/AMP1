c Set of modules for the computation of nuclear reaction rates 
c using NACRE database
c
c started 22/07/05 by Michael Bazot
c
c
c
c Subroutine snrnacre sets the constants for reaction rate 
c calculation. It is called by s/r srncns when ivreng=8. 
c The constant are those given in Angulo et al. 1999 
c (see http://pntpm.ulb.ac.be/nacre.htm).
c
c
c Subroutine cnrnacre computes the reaction rates according
c to the analytical approximations given in Angulo et al. 1999.
c 
c problem with reaction 10 : O16(p,g)F17 -> shape is not
c 1/T923 but T9**-0.82... To be fixed in crnacre?


c------------------------------------------------------------------

      subroutine snrnacre

c  Note that the coefficients are for computing Na*lambda(i,j),
c  where Na is Avogadro's number and lambda(i,j) is the rate
c  of reactions per pair (i,j)

      implicit double precision (a-h,o-z)
c  Note: engenr.n.d.incl replaced by engenr.nnz.d.incl, 5/6/02
      include 'engenr.nnz.d.incl'
c

      common/rcncns/ a(knucst),b(knucst),s(knucst,isnuc),
     *  zz(knucst),zz12(2,knucst),zsmh,acno,awght(knucwg),q(knucq), 
     *  knucpm, kqvals, kawght
      

c
c  Set NACRE (Angulo et al. 1999) parameter set
c
c        write(istdpr,105)
        a(1)=4.08e-15
        b(1)=3.381
        s(1,1)=3.82
        s(1,2)=1.51
        s(1,3)=0.144
        s(1,4)=-1.14e-2
c
        a(2)=5.59e10
        b(2)=12.277
        s(2,1)=-0.135
        s(2,2)=2.54e-2
        s(2,3)=-1.29e-3
        s(2,4)=0.0000        
c
        a(3)=5.46e6
        b(3)=12.827
        s(3,1)=-0.307
        s(3,2)=8.81e-2
        s(3,3)=-1.06e-2
        s(3,4)=4.46e-4
c
        a(4)=2.61e5
        b(4)=10.264
        s(4,1)=-5.11e-2
        s(4,2)=4.68e-2
        s(4,3)=-6.60e-3
        s(4,4)=3.12e-4
c
        a(5)=4.83e7
        b(5)=15.231
        s(5,1)=-2.00
        s(5,2)=3.41
        s(5,3)=-2.43
        s(5,4)=0.0000

c
c  Z1 =   1.0  Z2 =   6.0 A1 =   1.00782  A2 =  12.00000
c
        a(6)=  2.00E+07
        b(6)=   13.692
        s(6,1)=    9.89
        s(6,2)=    -59.8
        s(6,3)=    266
        s(6,4)=0.0000
 
c
c  Z1 =   1.0  Z2 =   6.0 A1 =   1.00782  A2 =  13.00335
c

        a(7)=  9.57e7
        b(7)=   13.720
        s(7,1)=    3.56
        s(7,2)=  0.0000
        s(7,3)=  0.0000
        s(7,4)=  0.0000
c 
c
c  Z1 =   1.0  Z2 =   7.0 A1 =   1.00782  A2 =  15.00011
c
        a(8)=  1.12E+1
        b(8)=   15.253
        s(8,1)=    4.95
        s(8,2)=    143
        s(8,3)=    0.0000
        s(8,4)=    0.0000
c
c  Z1 =   1.0  Z2 =   7.0 A1 =   1.00782  A2 =  15.00011
c

        a(9)=  1.08E+09
        b(9)=   15.524
        s(9,1)=    6.15
        s(9,2)=    16.4
        s(9,3)=    0.0000
        s(9,4)=    0.0000

         
           
c        endif
c
c  Z1 =   1.0  Z2 =   8.0 A1 =   1.00782  A2 =  15.99492
c
        a(10)=  7.37E+07
        b(10)=   16.696
        s(10,1)=    0.0000
        s(10,2)=    0.0000
        s(10,3)=    0.0000
        s(10,4)=    0.0000
c

c problem with reaction 10 : O16(p,g)F17 -> shape is not
c 1/T923 but T9**-0.82... To be fixed in crnacre?

c
c  O17(H1,He4)N14
c
        a(11)=  9.20d8
        b(11)=   16.715
        s(11,1)=    -80.31
        s(11,2)=    2211
        s(11,3)=    0.0000
        s(11,4)=    0.0000
c
c Test 
c

        write(*,66) s(1,3),a(4),a(10),b(10)
 66     format('Some coefficientts, s(1,3)',e13.5,' a(4) :',e13.5,
     &  ' a(10) :',e13.5,' b(10) :',e13.5,' ...Seems to work')

c
c End test
c

        return
        end

c-------------------------------------------------------------------

      subroutine cnrnacre(tl,ider,norct1,secder,zzscr) 

c Started 30/07/05 by Michael Bazot

c Calculate nuclear reaction rates and their derivatives for 
c with coefficients from Angulo et al. (1999). The major change
c compared to the classical procedure in s/r engenr is that the 
c nuclear reaction rates is given with a polynomial part, function
c of T. In the previous works (Fowler et al. 1975,...), this polynom
c was given as a function of T**(1/3).

c Summary of the used parameters :
c
c 
c
c
c
c
c
c
c
c
c
c
c
c Note : 
c The reaction rate for the reaction O16(p,g)F17 does not have 
c the same form as the other one. It is thus treated separately.
c


      implicit double precision (a-h,o-z)

      include 'engenr.nnz.d.incl'
      parameter (krnrm1 = krnrmx + 1, iptdat = krnrm1 - 12)

c
c  cut-off parameter for exponential
c
      parameter(coexp = 120.d0)
c
      logical norct1,secder
      dimension tpw(4)
      dimension iptrnr(krnrm1,4),zzscr(krnrmx)


      common/rnratd/ al(10,krnrmx),norct
      common/rnrcnt/ irnrat, krnrat, irn12, irn13, irn14, 
     *  irn15c, irn15o, irn16, irn17
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
      common/rcncns/ a(knucst),b(knucst),s(knucst,isnuc),
     *  zz(knucst),zz12(2,knucst),zsmh,acno,awght(knucwg),q(knucq), 
     *  knucpm, kqvals, kawght
      common/rnrout/ zt(10),uw(10),rs(10),ee(20),alr(10),
     *  dd(10),sr(10),f(10),alc(10)

c
c
c                  1  2  3  4  5  6  7  8   9  10  11  12
      data iptrnr /1, 2, 3, 4, 5, 0, 0, 0,  0,  0,  0,  0, iptdat*0,
     *             1, 2, 3, 4, 5, 0, 8, 9, 10, 11, 10, 11, iptdat*0,
     *             1, 2, 3, 4, 5, 0, 6, 7,  0,  0,  0,  0, iptdat*0,
     *             1, 2, 3, 4, 5, 0, 6, 7,  8,  9, 10, 11, iptdat*0/


      save
c
      fxp(aa)= exp(min(coexp,max(aa,-coexp)))
c

      t=10.d0**tl
      t9=t/1.d9
      tpw(1)=t9
      t13=t9**0.3333333333333333d0
      do 5 l=2,4               
 5       tpw(l)=t9*tpw(l-1)

      do 50 k=1,krnrat       
c
c  skip reaction Be7(e-,nu)Li7
c
         if(k.eq.6) go to 50 

c
c  set index for reaction parameters 
c
         kpar = iptrnr(k,irnrat) 

         if(idgeng.ge.2) then  
            write(istdpr,221) k, kpar, a(kpar), b(kpar), zzscr(k), t13
         end if
c
c Test 
c
            write(*,221) k, kpar, a(kpar), b(kpar), zzscr(k), t9
c
c End test
c
         exx=b(kpar)/t13        ! exponential part has the same form in every case

c
c  test for no reactions
c
         if(a(kpar).le.0) stop 'rnrate' 
         aexx=log(a(kpar))-exx  ! exponential part has the same form in every case
c
c
c
         write(*,*) 'aexx = ',aexx,' coexp = ',coexp
c
c End test
c
         

         if(aexx.le.-coexp) then 
c
c Test 
c
            write(*,*) 'Null reaction rate'
c
c End test
c


c
c  zero reaction rates and logarithmic derivatives
c
            do 36 i=1,ider
 36            al(i,k)=0
               go to 50

            end if

c
c Test 
c
               write(*,*) 'k before skipping O16 reaction : ',k
c
c End test
c

c
c reaction O16(p,g)F17 treated separately
c
            if (kpar.ne.10) then


c
c  set reaction rates
c
            norct1=.false.
c     
            sum1=1
            sum2=0
            sum3=0
      
            do 38 l=1,4               
               t1=s(kpar,l)*tpw(l)      
               sum1=sum1+t1            
               sum2=sum2+l*t1   ! for d(ln alt)/d(ln T)
 38            sum3=sum3+l*l*t1 ! for d2(ln alt)/d(ln T)2

      
               if(idgeng.ge.2) then 
                  write(istdpr,222) k, aexx, sum1, sum2, sum3, tpw(2)
               end if


               alt=fxp(aexx)*sum1/t13**2 
               dalt=exx/3-0.6666666666666667d0+sum2/sum1 !d(ln alt)/d(ln T) 
               if(secder) then
                  ddalt=-amm*(exx/9+(sum2*sum2/sum1-sum3)/sum1)
               end if


               al(1,k)= exp(zzscr(k)*uw(1))*alt !setting of reaction rates including the screening factor
               do 45 i=2,ider                   
 45               al(i,k)=zzscr(k)*uw(i)
                  al(3,k)=al(3,k)+dalt
                  if(secder) then
                     al(8,k)=al(8,k)+ddalt
                  end if
c     
c     reaction O16(p,g)F17 
c     
               else
         
                  alt=fxp(aexx)*t**(-0.82)
                  dalt=exx/3-0.82 !d(ln alt)/d(ln T) 
                  if(secder) then
                     ddalt=-amm*(exx/9)
                  end if
         
                  al(1,k)= exp(zzscr(k)*uw(1))*alt
                  do 46 i=2,ider                   
 46                  al(i,k)=zzscr(k)*uw(i)
                     al(3,k)=al(3,k)+dalt
                     if(secder) then
                        al(8,k)=al(8,k)+ddalt
                     end if
c
c Test           
c            
                     write(*,*) 'reaction O16(p,g)F17 reaction rate =', 
     &                al(1,k)

c
c End test
c

                  end if
 50            continue

c
c Test
c
                     
                     write(*,*) 'Reaction rates in cnrnacre'
                     do j=1,11
                        write(*,*) al(1,j)
                     enddo

c
c End test
c

 221           format(' start setting lambda for k =',i3,'  kpar =',i3/
     *  ' a(k), b(k), zzscr(k), t9 =',1p4e13.5)
 222           format(' k =',i3,'  aexx, sum1, sum2, sum3=',
     *  1p4e13.5)
               
      return
      end
