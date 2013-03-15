      subroutine setf4(nspe)
c
c  Note: commons free, map and spes have been renamed freecl,
c  mapcl and spescl to avoid conflict when calling routine
c  from MHD table emulator.
c                                           1/5/90       jcd
c
c  Modified 11/5/90 to allow arbitrary order of expansion for tau
c  in s/r f4mhd (modification necessary because of unavailability
c  of real*16).                                          jcd
c
c  Modified 5/6/90 to reset constants with JC-D values for consistency.
c  Note: this requires that s/r setcns be called before calling setf4.
c  In JC-D usage, setf4 would normally be called from setcns.
c
c  Modified 20/6/03: Include common/cdgphs/ kdgeos, kdgopc, kdgeng.
c  If, on input, kdgeos .gt. 0, do not stop on fatal error, but
c  return with kdgeos = -1.
c
c======================================================================
c
      implicit real*8 (a-h,o-z)
c
      parameter ( mchem =  3, mspes = 18, mz = 9, mion = mz + 1 )
      parameter ( mlam  = mspes - mchem -1 )
      parameter ( mfe   =   4, mfd    =  5)
c
      common /cstmhd/
     .                                camu   ,               cc     ,
     .                                ce     ,               ch     ,
     .                                ck     ,               cme    ,
     .                                cpi    ,               cevw   ,
     .                                cx     ,               cf4    ,
     .                                carad
c
      common /degmhd/
     .                                eta    ,               detdn  ,
     .                                detdt  ,               detdv  ,
     .                                d2etdn2,               d2etdnt,
     .                                d2etdnv,               d2etdt2,
     .                                d2etdtv,               d2etdv2,
     .                                exeta  ,               thet   ,
     .                                dthet  ,               d2thet ,
     .                                fd    (mfd   )
      common /freecl/
     .                         e     (mfe   ),        p     (mfe   ),
     .                         f     (mfe   ),               d2fdt2 ,
     .                                d2fdtv ,               d2fdv2 ,
     .                         fscr  (mspes ),        dfdn  (mspes ),
     .                         d2fdnt(mspes ),        d2fdnv(mspes ),
     .                         d2fdn2( mspes ,mspes )
      common /mapcl /
     .                                nchem  ,               nspes  ,
     .                                ise
      common /spescl/
     .                         zs    (mspes ),        zsq   (mspes ),
     .                         zmask (mspes ),        sn    (mspes )
c
c  commons of fundamental constants set by JC-D routine setcns.
c  for consistency, replace WD values by these.
c
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev1,exhm
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev, syear, iver
c
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
      common/cdgphs/ kdgeos, kdgopc, kdgeng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      if(idgeos.ge.1.and.istdpr.gt.0) write(istdpr,*) 'Calling setf4'
c
      camu  =  1.6605655d-24
      cc    =  2.9979246d+10
      ce    =  4.8032426d-10
      ch    =  6.6261764d-27
c
      ck    =  1.3806624d-16
      cme   =  9.1095345d-28
      cpi   =  3.141592654d0
      cevw  =  8.0654653d+03
c
c  resetting with JC-D values
c
      camu  =  amu
      cc    =  clight
      ce    =  echar
      ch    =  planck
      ck    =  boltzm
      cme   =  ame
      cpi   =  pi
c
c  Note: the meaning of cevw is currently unclear. On the other
c  hand, it appears not to be used. Should be fixed up.
      
c
      cf4   = 2.0d0 * dsqrt( cpi ) * ce**3 / ( 3.0d0 * dsqrt( ck ) )
      cx    = 3.0d0 * cf4 / ck
c
      carad =  7.56567  d-15
c
c  JC-D value
c  
      carad =  car
c
      nspes =  nspe
      ise   =  nspe
c
      zmask(1) = 0.d0
      zmask(2) = 1.d0
      zmask(3) = 0.d0
      zmask(4) = 1.d0
      zmask(5) = 1.d0
c
      zs(1) = 0.d0
      zs(2) = 1.d0
      zs(3) = 0.d0
      zs(4) = 1.d0
      zs(5) = 2.d0
c
      zsq(1) = 0.d0
      zsq(2) = 1.d0
      zsq(3) = 0.d0
      zsq(4) = 1.d0
      zsq(5) = 4.d0
c
c................... following values should never be used; if none the less,
c................... they ought to provoke a crash of the program!
      zmask(6) = -1.d200
      zs(6)    = -1.d200
      zsq(6)   = -1.d200
c
      return
      end
c
c======================================================================
      subroutine f4der(rhol,tl,sn,nspe,f4,e4,p4,h4,d2f4r2,d2f4t2,
     .           d2f4rt,df4,d2f4i,d2f4f,d2f4x,d2f4r,d2f4t,
     .           dp4i,dp4dr,dp4dt,dp4dx,dh4i,dh4dr,dh4dt,dh4dx,npar)
c======================================================================
c
c     nspe  : number of particles actually used
c     npar  : first dimension of d2f4i,d2f4f (=number of ionization
c             degrees = number of number fractions)
c
c     f4    : Coulomb configurational free energy (Debye-Huckel)
c     e4    : internal energy
c     p4    : pressure
c     h4    : enthalpy 
c
c     df4   : derivatives with respect to reaction parameters ("Saha equations")
c     d2f4i : derivatives of df4 with respect to ionization degrees
c     d2f4f : derivatives of df4 with respect to number fractions
c     d2f4x : derivatives of df4 with respect to X (Y varying accordingly,
c             and at fixed ionization degrees)
c
c     d2f4r : derivatives of df4 with respect to log10 rho
c     d2f4t : derivatives of df4 with respect to log10 t
c     d2f4r2: derivative  of  f4 with respect to log10 rho ** 2
c     d2f4t2: derivative  of  f4 with respect to log10 t   ** 2
c     d2f4rt: derivative  of  f4 with respect to log10 rho and log10 t
c
c     in both d2f4i and d2f4f the element (i,j) denotes the derivative
c     of df4(i) with respect to parameter j (ionization degree or number
c     fraction)
c
c     dp4i  : derivatives of p4 with respect to ionization degrees
c     dp4dr : derivative  of p4 with respect to log10 rho
c     dp4dt : derivative  of p4 with respect to log10 t
c     dp4dx : derivative  of p4 with respect to X (Y varying accordingly
c             and at fixed ionization degrees)
c
c     dh4i  : derivatives of h4 with respect to ionization degrees
c     dh4dr : derivative  of h4 with respect to log10 rho
c     dh4dt : derivative  of h4 with respect to log10 t
c     dh4dx : derivative  of h4 with respect to X (Y varying accordingly
c             and at fixed ionization degrees)
c
c=========================================================================
c
c  Modified 5/11/95 to allow switching off tau correction. This
c  is controlled by setting the integer variable notau in 
c  common/cnttau/ to 1, in
c
c=========================================================================
c
      implicit real*8 (a-h,o-z)
c
      parameter ( mchem =  3, mspes = 18, mz = 9, mion = mz + 1 )
      parameter ( mlam  = mspes - mchem -1 )
      parameter ( mfe   =   4, mfd    =  5)
c
      common /freecl/
     .                         e     (mfe   ),        p     (mfe   ),
     .                         f     (mfe   ),               d2fdt2 ,
     .                                d2fdtv ,               d2fdv2 ,
     .                         fscr  (mspes ),        dfdn  (mspes ),
     .                         d2fdnt(mspes ),        d2fdnv(mspes ),
     .                         d2fdn2( mspes ,mspes )
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev1,exhm
c
      dimension sn(nspe),df4(npar),d2f4i(npar,npar),d2f4f(npar,npar),
     .          d2f4x(npar),d2f4r(npar),d2f4t(npar),
     .          dp4i(npar),dh4i(npar)
c
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
      common/cdgphs/ kdgeos, kdgopc, kdgeng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      data umod /2.302585092994046 d0/
c
c...assume ah = 1.008, ahe = 4.0026
      data hovhe,heovh /0.251836,3.97083/
c
c  reset with JC-D values
c
      hovhe = ah/ahe
      heovh = ahe/ah
c
      t   = 10.d0**tl
      vol = 10.d0**(-rhol)
c
      call f4n(t,vol,sn,nspe)
c
c  test for return with error
c
      if(kdgeos.lt.0) return
c
c     f4    : Coulomb configurational free energy (Debye-Huckel)
c     e4    : internal energy
c     p4    : pressure
c     p4    : enthalpy
c
      f4     =  f(4)
      e4     =  e(4)
      p4     =  p(4)
      h4     =  e(4) + p(4)*vol
c
c     d2f4r2: derivative of  f4 with respect to log10 rho ** 2
c     d2f4t2: derivative of  f4 with respect to log10 t   ** 2
c     d2f4rt: derivative of  f4 with respect to log10 rho and log10 t
c
      d2f4r2 =  d2fdv2*vol*vol*umod*umod
      d2f4rt =  d2fdtv*  t*vol*umod*umod
      d2f4t2 =  d2fdt2*  t*  t*umod*umod
c
c     df4   : derivatives with respect to reaction parameters ("Saha equations")
c
      df4(1) =  dfdn(2) + dfdn(6)
      df4(2) =  dfdn(4) + dfdn(6)
      df4(3) = -dfdn(4) + dfdn(5) + dfdn(6)
      if(idgeos.eq.-101.and.tl.le.4) write(istdpr,*) '#D# in f4der',
     *  tl, dfdn(2),dfdn(6), df4(1)
c
c     d2f4r : derivatives of df4 with respect to log10 rho
c
      d2f4r(1) = -(             d2fdnv(2) + d2fdnv(6))*vol*umod
      d2f4r(2) = -(             d2fdnv(4) + d2fdnv(6))*vol*umod
      d2f4r(3) = -(-d2fdnv(4) + d2fdnv(5) + d2fdnv(6))*vol*umod
c
c     d2f4t : derivatives of df4 with respect to log10 t
c
      d2f4t(1) =  (             d2fdnt(2) + d2fdnt(6))*t  *umod
      d2f4t(2) =  (             d2fdnt(4) + d2fdnt(6))*t  *umod
      d2f4t(3) =  (-d2fdnt(4) + d2fdnt(5) + d2fdnt(6))*t  *umod
c
      toth   =  sn(1) + sn(2)
      tothe  =  sn(3) + sn(4) + sn(5)
c
c     d2f4i : derivatives of df4 with respect to ionization degrees
c
      d2f4i(1,1) = toth *(-d2fdn2(2,1)+d2fdn2(2,2)+d2fdn2(2,6)
     .                    -d2fdn2(6,1)+d2fdn2(6,2)+d2fdn2(6,6))
      d2f4i(1,2) = tothe*(-d2fdn2(2,3)+d2fdn2(2,4)+d2fdn2(2,6)
     .                    -d2fdn2(6,3)+d2fdn2(6,4)+d2fdn2(6,6))
      d2f4i(1,3) = tothe*(-d2fdn2(2,3)+d2fdn2(2,5)+d2fdn2(2,6)*2.d0
     .                    -d2fdn2(6,3)+d2fdn2(6,5)+d2fdn2(6,6)*2.d0)
c
      d2f4i(2,1) = toth *(-d2fdn2(4,1)+d2fdn2(4,2)+d2fdn2(4,6)
     .                    -d2fdn2(6,1)+d2fdn2(6,2)+d2fdn2(6,6))
      d2f4i(2,2) = tothe*(-d2fdn2(4,3)+d2fdn2(4,4)+d2fdn2(4,6)
     .                    -d2fdn2(6,3)+d2fdn2(6,4)+d2fdn2(6,6))
      d2f4i(2,3) = tothe*(-d2fdn2(4,3)+d2fdn2(4,5)+d2fdn2(4,6)*2.d0
     .                    -d2fdn2(6,3)+d2fdn2(6,5)+d2fdn2(6,6)*2.d0)
c
      d2f4i(3,1) = toth *( d2fdn2(4,1)-d2fdn2(4,2)-d2fdn2(4,6)
     .                    -d2fdn2(5,1)+d2fdn2(5,2)+d2fdn2(5,6)
     .                    -d2fdn2(6,1)+d2fdn2(6,2)+d2fdn2(6,6))
      d2f4i(3,2) = tothe*( d2fdn2(4,3)-d2fdn2(4,4)-d2fdn2(4,6)
     .                    -d2fdn2(5,3)+d2fdn2(5,4)+d2fdn2(5,6)
     .                    -d2fdn2(6,3)+d2fdn2(6,4)+d2fdn2(6,6))
      d2f4i(3,3) = tothe*( d2fdn2(4,3)-d2fdn2(4,5)-d2fdn2(4,6)*2.d0
     .                    -d2fdn2(5,3)+d2fdn2(5,5)+d2fdn2(5,6)*2.d0
     .                    -d2fdn2(6,3)+d2fdn2(6,5)+d2fdn2(6,6)*2.d0)
c
      fach       =  sn(1)*sn(1)/toth
c
c  Note: temporary fudge to avoid problems when sn(3) or sn(4)
c  are zero. Need correction from WD           10/5/90
c
      if(sn(3).ne.0) then
        ff2      =  sn(4)/sn(3)
      else
        ff2      =  1
      end if
      if(sn(4).ne.0) then
        ff3      =  sn(5)/sn(4)
      else
        ff3      =  1
      end if
      ff22       =  ff2*ff2
      fac2       =  1.d0 + ff2
      fac3       =  1.d0 + ff3
      fache      =  sn(3)*sn(3)/tothe
      dhed2      = -fac3*fache
      dhed3      = -ff2*fache
      dhepd2     =  fache
      dhepd3     = -ff22*fache
      dhe2d2     =  ff3*fache
      dhe2d3     =  ff2*fac2*fache
      ded2       =  (1.d0 + ff3*2.d0)*fache
      ded3        = ff2*(2.d0 + ff2)*fache
c
c     d2f4f : derivatives of df4 with respect to number fractions
c
      d2f4f(1,1) = fach*(-d2fdn2(2,1)+d2fdn2(2,2)+d2fdn2(2,6)
     .                   -d2fdn2(6,1)+d2fdn2(6,2)+d2fdn2(6,6))
      d2f4f(2,1) = fach*(-d2fdn2(4,1)+d2fdn2(4,2)+d2fdn2(4,6)
     .                   -d2fdn2(6,1)+d2fdn2(6,2)+d2fdn2(6,6))
      d2f4f(3,1) = fach*( d2fdn2(4,1)-d2fdn2(4,2)-d2fdn2(4,6)
     .                   -d2fdn2(5,1)+d2fdn2(5,2)+d2fdn2(5,6)
     .                   -d2fdn2(6,1)+d2fdn2(6,2)+d2fdn2(6,6))
c
      d2f4f(1,2) =    ( d2fdn2(2,3)+d2fdn2(6,3))*dhed2
     .              + ( d2fdn2(2,4)+d2fdn2(6,4))*dhepd2
     .              + ( d2fdn2(2,5)+d2fdn2(6,5))*dhe2d2
     .              + ( d2fdn2(2,6)+d2fdn2(6,6))*ded2
      d2f4f(2,2) =    ( d2fdn2(4,3)+d2fdn2(6,3))*dhed2
     .              + ( d2fdn2(4,4)+d2fdn2(6,4))*dhepd2
     .              + ( d2fdn2(4,5)+d2fdn2(6,5))*dhe2d2
     .              + ( d2fdn2(4,6)+d2fdn2(6,6))*ded2
      d2f4f(3,2) =    (-d2fdn2(4,3)+d2fdn2(5,3)+d2fdn2(6,3))*dhed2
     .              + (-d2fdn2(4,4)+d2fdn2(5,4)+d2fdn2(6,4))*dhepd2
     .              + (-d2fdn2(4,5)+d2fdn2(5,5)+d2fdn2(6,5))*dhe2d2
     .              + (-d2fdn2(4,6)+d2fdn2(5,6)+d2fdn2(6,6))*ded2
c
      d2f4f(1,3) =    ( d2fdn2(2,3)+d2fdn2(6,3))*dhed3
     .              + ( d2fdn2(2,4)+d2fdn2(6,4))*dhepd3
     .              + ( d2fdn2(2,5)+d2fdn2(6,5))*dhe2d3
     .              + ( d2fdn2(2,6)+d2fdn2(6,6))*ded3
      d2f4f(2,3) =    ( d2fdn2(4,3)+d2fdn2(6,3))*dhed3
     .              + ( d2fdn2(4,4)+d2fdn2(6,4))*dhepd3
     .              + ( d2fdn2(4,5)+d2fdn2(6,5))*dhe2d3
     .              + ( d2fdn2(4,6)+d2fdn2(6,6))*ded3
      d2f4f(3,3) =    (-d2fdn2(4,3)+d2fdn2(5,3)+d2fdn2(6,3))*dhed3
     .              + (-d2fdn2(4,4)+d2fdn2(5,4)+d2fdn2(6,4))*dhepd3
     .              + (-d2fdn2(4,5)+d2fdn2(5,5)+d2fdn2(6,5))*dhe2d3
     .              + (-d2fdn2(4,6)+d2fdn2(5,6)+d2fdn2(6,6))*ded3
c
c     d2f4x : derivatives of df4 with respect to X (Y varying accordingly)
c
      dnhdx    =  (toth + tothe*heovh)/toth
      dnhedx   = -(toth*hovhe + tothe)/tothe
c
      d2f4x(1) = dnhdx *((d2fdn2(2,1)+d2fdn2(6,1))* sn(1) +
     .                   (d2fdn2(2,2)+d2fdn2(6,2))* sn(2) +
     .                   (d2fdn2(2,6)+d2fdn2(6,6))* sn(2))+
     .           dnhedx*((d2fdn2(2,3)+d2fdn2(6,3))* sn(3) +
     .                   (d2fdn2(2,4)+d2fdn2(6,4))* sn(4) +
     .                   (d2fdn2(2,5)+d2fdn2(6,5))* sn(5) +
     .                   (d2fdn2(2,6)+d2fdn2(6,6))*(sn(4)+2.d0*sn(5)))
c
      d2f4x(2) = dnhdx *((d2fdn2(4,1)+d2fdn2(6,1))* sn(1) +
     .                   (d2fdn2(4,2)+d2fdn2(6,2))* sn(2) +
     .                   (d2fdn2(4,6)+d2fdn2(6,6))* sn(2))+
     .           dnhedx*((d2fdn2(4,3)+d2fdn2(6,3))* sn(3) +
     .                   (d2fdn2(4,4)+d2fdn2(6,4))* sn(4) +
     .                   (d2fdn2(4,5)+d2fdn2(6,5))* sn(5) +
     .                   (d2fdn2(4,6)+d2fdn2(6,6))*(sn(4)+2.d0*sn(5)))
c
      d2f4x(3) = dnhdx *((-d2fdn2(4,1)+d2fdn2(5,1)+d2fdn2(6,1))*
     .                     sn(1) +
     .                   (-d2fdn2(4,2)+d2fdn2(5,2)+d2fdn2(6,2))*
     .                     sn(2) +
     .                   (-d2fdn2(4,6)+d2fdn2(5,6)+d2fdn2(6,6))* 
     .                     sn(2))+
     .           dnhedx*((-d2fdn2(4,3)+d2fdn2(5,3)+d2fdn2(6,3))*
     .                     sn(3) +
     .                   (-d2fdn2(4,4)+d2fdn2(5,4)+d2fdn2(6,4))*
     .                     sn(4) +
     .                   (-d2fdn2(4,5)+d2fdn2(5,5)+d2fdn2(6,5))*
     .                     sn(5) +
     .                   (-d2fdn2(4,6)+d2fdn2(5,6)+d2fdn2(6,6))*
     .                    (sn(4)+2.d0*sn(5)))
c
c     dp4i  : derivatives of p4 with respect to ionization degrees
c     dp4dr : derivative  of p4 with respect to log10 rho
c     dp4dt : derivative  of p4 with respect to log10 t
c     dp4dx : derivative  of p4 with respect to X (Y varying accordingly
c             and at fixed ionization degrees)
c
      dp4i(1) = - toth *( -d2fdnv(1)+d2fdnv(2)+d2fdnv(6))
      dp4i(2) = - tothe*( -d2fdnv(3)+d2fdnv(4)+d2fdnv(6))
      dp4i(3) = - tothe*( -d2fdnv(3)+d2fdnv(5)+d2fdnv(6)*2.d0)
c
      dp4dr   =    umod  * vol * d2fdv2
      dp4dt   =  - umod  *   t * d2fdtv
c
      dp4dx   = -dnhdx * (d2fdnv(1)* sn(1) +
     .                    d2fdnv(2)* sn(2) +
     .                    d2fdnv(6)* sn(2))
     .          -dnhedx* (d2fdnv(3)* sn(3) +
     .                    d2fdnv(4)* sn(4) +
     .                    d2fdnv(5)* sn(5) +
     .                    d2fdnv(6)*(sn(4)+2.d0*sn(5)))
c
c     ... analogous for h4 ...
c
      dh4i(1) =   toth *((-dfdn  (1)+dfdn  (2)+dfdn  (6))            
     .                  -(-d2fdnt(1)+d2fdnt(2)+d2fdnt(6))*t 
     .                  -(-d2fdnv(1)+d2fdnv(2)+d2fdnv(6))*vol  )
c
      dh4i(2) =   tothe*((-dfdn  (3)+dfdn  (4)+dfdn  (6))            
     .                  -(-d2fdnt(3)+d2fdnt(4)+d2fdnt(6))*t 
     .                  -(-d2fdnv(3)+d2fdnv(4)+d2fdnv(6))*vol  )
c
      dh4i(3) =   tothe*((-dfdn  (3)+dfdn  (5)+dfdn  (6)*2.d0)            
     .                  -(-d2fdnt(3)+d2fdnt(5)+d2fdnt(6)*2.d0)*t 
     .                  -(-d2fdnv(3)+d2fdnv(5)+d2fdnv(6)*2.d0)*vol  )
c
      dh4dr   =    umod * vol * t * d2fdtv  +  vol * dp4dr
      dh4dt   =  - umod *  t  * t * d2fdt2  +  vol * dp4dt
c
      dh4dx   = dnhdx *((dfdn(1)-t*d2fdnt(1))* sn(1) +
     .                  (dfdn(2)-t*d2fdnt(2))* sn(2) +
     .                  (dfdn(6)-t*d2fdnt(6))* sn(2))+
     .          dnhedx*((dfdn(3)-t*d2fdnt(3))* sn(3) +
     .                  (dfdn(4)-t*d2fdnt(4))* sn(4) +
     .                  (dfdn(5)-t*d2fdnt(5))* sn(5) +
     .                  (dfdn(6)-t*d2fdnt(6))*(sn(4) + 2.d0*sn(5)))+
     .          vol * dp4dx
c
      return
      end
c
c======================================================================
      subroutine f4n(t,vol,snn,nspe)
c======================================================================
c
c     derivatives with respect to number abundances.
c     calls f4 of MHD package and prepares quantities not provided by f4
c     (degeneracy-related stuff).
c
c === all units c.g.s. and degrees Kelvin
c
c
c  Modified 20/6/03: Include common/cdgphs/ kdgeos, kdgopc, kdgeng.
c  If, on input, kdgeos .gt. 0, do not stop on fatal error, but
c  return with kdgeos = -1.

      implicit real*8 (a-h,o-z)
c
      parameter ( mchem =  3, mspes = 18, mz = 9, mion = mz + 1 )
      parameter ( mlam  = mspes - mchem -1 )
      parameter ( mfe   =   4, mfd    =  5)
c
      dimension snn(nspe)
c
      common /cstmhd/
     .                                camu   ,               cc     ,
     .                                ce     ,               ch     ,
     .                                ck     ,               cme    ,
     .                                cpi    ,               cevw   ,
     .                                cx     ,               cf4    ,
     .                                carad
c
      common /degmhd/
     .                                eta    ,               detdn  ,
     .                                detdt  ,               detdv  ,
     .                                d2etdn2,               d2etdnt,
     .                                d2etdnv,               d2etdt2,
     .                                d2etdtv,               d2etdv2,
     .                                exeta  ,               thet   ,
     .                                dthet  ,               d2thet ,
     .                                fd    (mfd   )
      common /freecl/
     .                         e     (mfe   ),        p     (mfe   ),
     .                         f     (mfe   ),               d2fdt2 ,
     .                                d2fdtv ,               d2fdv2 ,
     .                         fscr  (mspes ),        dfdn  (mspes ),
     .                         d2fdnt(mspes ),        d2fdnv(mspes ),
     .                         d2fdn2( mspes ,mspes )
      common /mapcl /
     .                                nchem  ,               nspes  ,
     .                                ise
      common /spescl/
     .                         zs    (mspes ),        zsq   (mspes ),
     .                         zmask (mspes ),        sn    (mspes )
c
      common/cdgphs/ kdgeos, kdgopc, kdgeng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      do 5 is=1,nspes
 5    sn(is) = snn(is)
c
      ckt   = ck * t
c
      call neweta(t,vol,sn(ise),ier)
      if(ier.ne.0) then
         write(istdou,*) 
     *     'error in s/r f4n. failure in neweta. f4 not called'
         if(istdpr.gt.0.and.istdou.ne.istdpr) write(istdpr,*) 
     *     'error in s/r f4n. failure in neweta. f4 not called'
	 if(kdgeos.gt.0) then
	   kdgeos=-1
	   return
         else
           stop 'f4n'
         end if
      end if
c
      call f4mhd(t,vol)
c
      return
      end
c
      subroutine neweta (t,vol,sne,ier)
c
c***********************************************************************
c     calculate degeneracy parameter and its                           c
c     derivatives. evaluate theta and its derivatives                  c
c     ................................
c     (modified from MHD package)
c     ................................
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
c
      parameter (mfd  =  5)
c
      common /degmhd/
     .                                eta    ,               detdn  ,
     .                                detdt  ,               detdv  ,
     .                                d2etdn2,               d2etdnt,
     .                                d2etdnv,               d2etdt2,
     .                                d2etdtv,               d2etdv2,
     .                                exeta  ,               thet   ,
     .                                dthet  ,               d2thet ,
     .                                fd    (mfd   )
      common /cstmhd/
     .                                camu   ,               cc     ,
     .                                ce     ,               ch     ,
     .                                ck     ,               cme    ,
     .                                cpi    ,               cevw   ,
     .                                cx     ,               cf4    ,
     .                                carad
c
      data   niter,   conv
     .    /    20, 1.d-10 /
c
c........... the following factorizing in gse is made to avoid
c........... overflows on simple machines with ranges of about 1.e37
c
      gse  = 2.*cpi*(cme/ch)*(ck/ch)
      gse  = gse**1.5
c
      ceta = sqrt( cpi ) / ( 4.0 * gse )
c
c..      data   niter,   conv
c..     .    /     20, 1.e-10 /
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialization                                                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ier = 0
      rhs = ceta * sne   / ( t**1.5 * vol )
ccc      write(6,*) ' rhs,ceta = ',rhs,ceta
c
c........... initialization to a few per cent accuracy
c........... (see ref. in dappen, astron. astrphys., 1980)
c
      if(rhs.lt.1.d-6) then
            eta = dlog(rhs)
      else if(rhs.lt.100.d0) then
            g   =rhs
            g2  =rhs*rhs
            et0 =g+(g2/(4.45+2.257*g+g2))*dexp((1.32934*g)**0.66666667)
            eta = dlog(et0)
      else
            eta = (1.32934*rhs)**0.66666667
      end if
      eta = eta + 0.120782
c
ccc      write(6,*) 'eta,rhs ini: ',eta,rhs
      call ferdir( eta, fd )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     newton-raphson iteration
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      iter  = 0
    1 deta  = 2.0d0 * ( rhs - fd(2) ) / fd(1)
      eta = eta + deta
      call ferdir( eta, fd )
      if( dabs(deta) .le. conv ) go to 3
      iter = iter + 1
      if( iter .le. niter ) go to 1
c
c     failure to converge
      write  ( 6, 2 ) t, vol, sne, eta
    2 format ( ' nonconvergence of degeneracy parameter ' /
     .         ' t =',1pg10.2 ',vol =',g10.2 ,'ne =',g10.3,
     .         ' eta   =',g12.4 )
      ier = 1
      return
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     convergence                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
    3 exeta = dexp( - eta )
c
      thet   =  0.5d0 *   fd(1) / fd(2)
      dthet  = thet * ( fd(4) / fd(1) - thet )
      d2thet = thet * ( fd(5) / fd(1) - 3.0d0 * dthet - thet**2 )
c
      detdv   = - 1.0d0 / ( thet * vol )
      detdt   = - 1.5d0 / ( thet * t )
      detdn   =   1.0d0 / ( thet * sne )
      d2etdn2 = - detdn**2 * fd(4) / fd(1)
      d2etdnt = - dthet * detdn * detdt / thet
      d2etdnv = - dthet * detdn * detdv / thet
      d2etdtv = - dthet * detdt * detdv / thet
      d2etdt2 = - ( 1.0d0 / t   + dthet * detdt / thet ) * detdt
      d2etdv2 = - ( 1.0d0 / vol + dthet * detdv / thet ) * detdv
c
      if(iter.gt.5) write (6,*)' slow convergence in neweta: iter,eta',
     . ',t,vol = ',iter,eta,t,vol
c
ccc      write(6,*) 'eta,fd ',eta,fd
      return
      end
c
      subroutine ferdir( x, fd )
c     ................................
c     (unmodified from MHD package)
c     ................................
      implicit real*8 (a-h,o-z)
c
c***********************************************************************
c     calculate fermi-dirac integrals                                  c
c                                                                      c
c     fd(1) = f    (x)                                                 c
c              -1/2                                                    c
c                                                                      c
c     fd(2) = f    (x)                                                 c
c               1/2                                                    c
c                                                                      c
c     fd(3) = f    (x)                                                 c
c               3/2                                                    c
c                                                                      c
c     fd(4) = f'   (x)                                                 c
c              -1/2                                                    c
c                                                                      c
c     fd(5) = f"   (x)                                                 c
c              -1/2                                                    c
c***********************************************************************
c
      dimension fd(5)
      dimension p1(5),p2(5),p3(5),p4(5),p5(5),p6(5),p7(5),p8(5),p9(5)
      dimension q1(5),q2(5),q3(5),q4(5),q5(5),q6(5),q7(5),q8(5),q9(5)
c
      data p1
     ./-1.25331 41288 20d+0, -1.72366 35577 01d+0, -6.55904 57292 58d-1,
     . -6.34228 31976 82d-2, -1.48838 31061 16d-5/
      data q1
     ./+1.00000 00000 00d+0, +2.19178 09259 80d+0, +1.60581 29554 06d+0,
     . +4.44366 95274 81d-1, +3.62423 22881 12d-2/
      data p2
     ./-3.13328 53055 70d-1, -4.16187 38522 93d-1, -1.50220 84005 88d-1,
     . -1.33957 93751 73d-2, -1.51335 07001 38d-5/
      data q2
     ./+1.00000 00000 00d+0, +1.87260 86759 02d+0, +1.14520 44465 78d+0,
     . +2.57022 55875 73d-1, +1.63990 25435 68d-2/
      data p3
     ./-2.34996 39854 06d-1, -2.92737 36375 47d-1, -9.88309 75887 38d-2,
     . -8.25138 63795 51d-3, -1.87438 41532 23d-5/
      data q3
     ./+1.00000 00000 00d+0, +1.60859 71091 46d+0, +8.27528 95308 80d-1,
     . +1.52232 23828 50d-1, +7.69512 04750 64d-3/
      data p4
     ./+1.07381 27694 00d+0, +5.60033 03660 00d+0, +3.68822 11270 00d+0,
     . +1.17433 92816 00d+0, +2.36419 35527 00d-1/
      data q4
     ./+1.00000 00000 00d+0, +4.60318 40667 00d+0, +4.30759 10674 00d-1,
     . +4.21511 32145 00d-1, +1.18326 01601 00d-2/
      data p5
     ./+6.78176 62666 00d-1, +6.33124 01791 00d-1, +2.94479 65177 20d-1,
     . +8.01320 71141 90d-2, +1.33918 21294 00d-2/
      data q5
     ./+1.00000 00000 00d+0, +1.43740 40039 70d-1, +7.08662 14845 00d-2,
     . +2.34579 49473 50d-3, -1.29449 92883 50d-5/
      data p6
     ./+1.15302 13402 00d+0, +1.05915 58972 00d+0, +4.68988 03095 00d-1,
     . +1.18829 08784 00d-1, +1.94387 55787 00d-2/
      data q6
     ./+1.00000 00000 00d+0, +3.73489 53841 00d-2, +2.32484 58137 00d-2,
     . -1.37667 70874 00d-3, +4.64663 92781 00d-5/
      data p7
     ./-8.22255 93300 00d-1, -3.62036 93450 00d+1, -3.01538 54100 00d+3,
     . -7.04987 15790 00d+4, -5.69814 59240 00d+4/
      data q7
     ./+1.00000 00000 00d+0, +3.93568 98410 00d+1, +3.56875 62660 00d+3,
     . +4.18189 36250 00d+4, +3.38513 89070 00d+5/
      data p8
     ./+8.22449 97626 00d-1, +2.00463 03393 00d+1, +1.82680 93446 00d+3,
     . +1.22265 30374 00d+4, +1.40407 50092 00d+5/
      data q8
     ./+1.00000 00000 00d+0, +2.34862 07659 00d+1, +2.20134 83743 00d+3,
     . +1.14426 73596 00d+4, +1.65847 15900 00d+5/
      data p9
     ./+2.46740 02368 40d+0, +2.19167 58236 80d+2, +1.23829 37907 50d+4,
     . +2.20667 72496 80d+5, +8.49442 92003 40d+5/
      data q9
     ./+1.00000 00000 00d+0, +8.91125 14061 90d+1, +5.04575 66966 70d+3,
     . +9.09075 94630 40d+4, +3.89960 91564 10d+5/
c
c
      if( x .gt. 4.0d0) go to 2
      if( x .gt. 1.0d0) go to 1
c
c
      y   = dexp(x)
c
      p   = y**2*(   p1(1) + y*(   p1(2) + y*(    p1(3) + y*(   p1(4)
     .                     + y*    p1(5)))))
      dp  = y**2*(2.*p1(1) + y*(3.*p1(2) + y*( 4.*p1(3) + y*( 5.*p1(4)
     .                     + y* 6.*p1(5)))))
      d2p = y**2*(4.*p1(1) + y*(9.*p1(2) + y*(16.*p1(3) + y*(25.*p1(4)
     .                     + y*36.*p1(5)))))
c
      q   = q1(1) +y*(q1(2) + y*(   q1(3) + y*(   q1(4) + y*    q1(5))))
      dq  =        y*(q1(2) + y*(2.*q1(3) + y*(3.*q1(4) + y* 4.*q1(5))))
      d2q =        y*(q1(2) + y*(4.*q1(3) + y*(9.*q1(4) + y*16.*q1(5))))
c
      fd(1) = 1.7724 53850 90552 d0*y + p/q
c
      fd(2) = y*(0.8862 26925 45276 d0 + y*
     .          (p2(1) + y*(p2(2) + y*(p2(3) + y*(p2(4) + y*p2(5)))))/
     .          (q2(1) + y*(q2(2) + y*(q2(3) + y*(q2(4) + y*q2(5))))))
c
      fd(3) = y*(1.3293 40388 17914 d0 + y*
     .          (p3(1) + y*(p3(2) + y*(p3(3) + y*(p3(4) + y*p3(5)))))/
     .          (q3(1) + y*(q3(2) + y*(q3(3) + y*(q3(4) + y*q3(5))))))
c
      fd(4) = 1.7724 53850 90552 d0*y + (dp*q - p*dq)/q**2
c
      fd(5) = 1.7724 53850 90552 d0*y +
     .               ((d2p*q - p*d2q)*q - 2.*(dp*q -p*dq)*dq)/q**3
c
      return
c
c
    1 p   =  p4(1)+ x*(p4(2) + x*(   p4(3) + x*(   p4(4) + x*   p4(5))))
      dp  =            p4(2) + x*(2.*p4(3) + x*(3.*p4(4) + x*4.*p4(5)))
      d2p =                       2.*p4(3) + 6.*x*(p4(4) + x*2.*p4(5))
c
      q   =  q4(1)+ x*(q4(2) + x*(   q4(3) + x*(   q4(4) + x*   q4(5))))
      dq  =            q4(2) + x*(2.*q4(3) + x*(3.*q4(4) + x*4.*q4(5)))
      d2q =                       2.*q4(3) + 6.*x*(q4(4) + x*2.*q4(5))
c
      fd(1) = p/q
c
      fd(2) = (p5(1)  + x*(p5(2) + x*(p5(3) + x*(p5(4) + x*p5(5)))))/
     .        (q5(1)  + x*(q5(2) + x*(q5(3) + x*(q5(4) + x*q5(5)))))
c
      fd(3) = (p6(1)  + x*(p6(2) + x*(p6(3) + x*(p6(4) + x*p6(5)))))/
     .        (q6(1)  + x*(q6(2) + x*(q6(3) + x*(q6(4) + x*q6(5)))))
c
      fd(4) = (dp *q - p*dq )/q**2
c
      fd(5) = (d2p*q - p*d2q)/q**2 - 2.*fd(4)*dq/q
c
      return
c
c
    2 root  = dsqrt(x)
      xsq   = x**2
      y     = 1.0d0/xsq
c
      p   =        y * (p7(1) + y*(   p7(2) + y*(   p7(3) + y*(    p7(4)
     .                        + y*    p7(5)))))
      dp  = (-2.*y/x)* (p7(1) + y*(2.*p7(2) + y*(3.*p7(3) + y*( 4.*p7(4)
     .                        + y* 5.*p7(5)))))
      d2p =  4.*y**2 * (p7(1) + y*(4.*p7(2) + y*(9.*p7(3) + y*(16.*p7(4)
     .                        + y*25.*p7(5))))) - dp/x
c
      q   = q7(1) + y * (q7(2) +y*(   q7(3) +y*(   q7(4) +y*    q7(5))))
      dq  = (-2.*y/x) * (q7(2) +y*(2.*q7(3) +y*(3.*q7(4) +y* 4.*q7(5))))
      d2q =  4.*y**2  * (q7(2) +y*(4.*q7(3) +y*(9.*q7(4) +y*16.*q7(5))))
     .                -  dq/x
c
      fd(1) = root * (2.0d0 + p/q)
c
      fd(2) = x*root*(0.66666 66666 66667 d0 + y*
     .        (p8(1) + y*(p8(2) + y*(p8(3)   + y*(p8(4) + y*p8(5)))))/
     .        (q8(1) + y*(q8(2) + y*(q8(3)   + y*(q8(4) + y*q8(5))))))
c
      fd(3) = xsq*root*(0.4d0   + y*
     .        (p9(1) + y*(p9(2) + y*(p9(3) + y*(p9(4) + y*p9(5)))))/
     .        (q9(1) + y*(q9(2) + y*(q9(3) + y*(q9(4) + y*q9(5))))))
c
      fd(4) =    fd(1)/(2.d0*x) +   root*(dp*q - p*dq)/q**2
c
      fd(5) = ( -fd(1)/x + fd(4))/(2.d0*x)
     .           + root* (dp *q - p*dq )/(2.*x*q**2)
     .           + root*((d2p*q - p*d2q) - 2.*dq*(dp*q - p*dq)/q)/q**2
c
      return
      end
c
      subroutine f4mhd(t,vol)
c     ................................
c     (modified from MHD package)
c     ................................
      implicit real*8 (a-h,o-z)
c
c***********************************************************************
c     free energy of coulomb interactions and derivatives              c
c***********************************************************************
c
c  Modified 5/11/95 to allow switching off tau correction. This
c  is controlled by setting the integer variable notau in 
c  common/cnttau/ to 1, in
c
      parameter ( mchem =  3, mspes = 18, mz = 9, mion = mz + 1 )
      parameter ( mlam  = mspes - mchem -1 )
      parameter ( mfe   =   4, mfd    =  5)
c
c  parameters for tau expansion. xtautr gives transition between
c  expansion and direct expression.
c  With nmax = 15, xtautr = 0.1, the error at the transition point
c  is around 1.e-13 in tau, 1.e-10 in dtau and 1.e-9 in d2tau.
c
      parameter ( nmax = 15, xtautr = 0.1)
c
c-------------------------------------------------------------------
c
      common /cstmhd/
     .                                camu   ,               cc     ,
     .                                ce     ,               ch     ,
     .                                ck     ,               cme    ,
     .                                cpi    ,               cevw   ,
     .                                cx     ,               cf4    ,
     .                                carad
c
      common /degmhd/
     .                                eta    ,               detdn  ,
     .                                detdt  ,               detdv  ,
     .                                d2etdn2,               d2etdnt,
     .                                d2etdnv,               d2etdt2,
     .                                d2etdtv,               d2etdv2,
     .                                exeta  ,               thet   ,
     .                                dthet  ,               d2thet ,
     .                                fd    (mfd   )
      common /freecl/
     .                         e     (mfe   ),        p     (mfe   ),
     .                         f     (mfe   ),               d2fdt2 ,
     .                                d2fdtv ,               d2fdv2 ,
     .                         fscr  (mspes ),        dfdn  (mspes ),
     .                         d2fdnt(mspes ),        d2fdnv(mspes ),
     .                         d2fdn2( mspes ,mspes )
      common /mapcl /
     .                                nchem  ,               nspes  ,
     .                                ise
      common /spescl/
     .                         zs    (mspes ),        zsq   (mspes ),
     .                         zmask (mspes ),        sn    (mspes )
c
      common/cnttau/ notau
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      data notau /0/
c
c
c  for running on Suns, take out quadruple precision.
c
c  *** Note: we may need to fix this up in expansion later, by
c      including more terms
c
c..      real*16 dpx, dp1, dp2, dp3
c
      dimension   dxdn(mspes), d2xdnt(mspes), d2xdnv(mspes),
     .            d2xdn2(mspes, mspes)
      dimension   ctau(nmax), cdtau(nmax), cd2tau(nmax)
c
      equivalence (dxdn  , dzdn  ), (d2xdnt, d2zdnt), (d2xdnv, d2zdnv),
     .            (d2xdn2, d2zdn2)
c
      data initcf /0/
c
      save initcf, ctau, cdtau, cd2tau
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     sums over charges                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      sumn = dmax1( sdot(nspes - 1, zmask, 1, sn, 1), 1.d-70 )
      zn   = dmax1( sdot(nspes - 1, zs   , 1, sn, 1), 1.d-70 )
      znt  = dmax1( sdot(nspes - 1, zsq  , 1, sn, 1) + sn(ise)*thet,
     .                                              1.d-70 )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     parameter x and its derivatives                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      frat = fd(2) / fd(3)
      x    = cx * frat * (zn/sumn) * dsqrt( znt ) / dsqrt( vol*t**3 )
c
c-----------------------------------------------------------------------
c     zero everything
c-----------------------------------------------------------------------
c
      do 2 is = 1, nspes
      fscr  (is) = 0.0d0
      dxdn  (is) = 0.0d0
      d2xdnt(is) = 0.0d0
      d2xdnv(is) = 0.0d0
      do 1 js = 1, nspes
      d2xdn2(is, js) = 0.0d0
    1 continue
    2 continue
c
c-----------------------------------------------------------------------
c     first derivatives                                                c
c-----------------------------------------------------------------------
c
      do 3 is = 1, nspes - 1
      dxdn(is) = x * ( zs(is)/zn - zmask(is)/sumn+0.5d0*zsq(is)/znt )
    3 continue
c
      dxdn(ise) = x * ( 0.5d0* (thet + sn(ise)*dthet*detdn)/znt
     .                - 1.5d0* frat * detdn
     .                + 1.0d0 / sn(ise) )
c
      dxdt = x * ( detdt*(thet - 1.5d0*frat + 0.5d0*sn(ise)*dthet/znt)
     .           - 1.5d0/t   )
      dxdv = x * ( detdv*(thet - 1.5d0*frat + 0.5d0*sn(ise)*dthet/znt)
     .           - 0.5d0/vol )
c
c-----------------------------------------------------------------------
c     second derivatives
c-----------------------------------------------------------------------
c
      do 5 js = 1, nspes - 1
c
      if ( zs (js) .ne. 0.0d0 ) then
            do 4 is = 1, js
            d2xdn2(is, js) = ( dxdn(is) * dxdn(js) / x
     .                   + x * (       zmask(is) * zmask(js) / sumn**2
     .                         -       zs   (is) * zs   (js) / zn**2
     .                       - 0.5d0 * zsq  (is) * zsq  (js) / znt**2) )
    4       continue
c
            if( js .gt. 1 ) call scopy( js - 1, d2xdn2( 1,js), 1,
     .                                          d2xdn2(js, 1), mspes )
      d2xdnt(js) = dxdt * dxdn(js) / x
     .               - 0.5d0 * x * zsq(js)*sn(ise) * dthet*detdt/znt**2
      d2xdnv(js) = dxdv * dxdn(js) / x
     .               - 0.5d0 * x * zsq(js)*sn(ise) * dthet*detdv/znt**2
      end if
    5 continue
c
      do 6 is = 1, nspes - 1
      d2xdn2(is, ise) = ( dxdn(is) * dxdn(ise) / x
     .   - 0.5d0 * x * zsq(is) * (thet + sn(ise)*dthet*detdn)/znt**2 )
    6 continue
      call scopy( nspes - 1, d2xdn2(1,ise), 1, d2xdn2(ise,1), mspes )
c
      d2xdn2(ise, ise) =
     .    dxdn(ise)**2/x
     .    - x * ( 1.0d0/sn(ise)**2
     .       + 1.5d0 * frat * (d2etdn2 + (thet - 1.5d0*frat)*detdn**2)
     .         - 0.5d0 * ( d2thet * sn(ise) * detdn**2
     .                  + dthet * (2.d0*detdn + sn(ise)*d2etdn2)
     .                - (thet + sn(ise)*dthet*detdn)**2/znt ) / znt  )
c
      d2xdnt(ise) = dxdt*dxdn(ise)/x + x * (
     .  -1.5d0*frat* (d2etdnt + (thet - 1.5d0*frat)*detdn*detdt)
     .  +0.5d0*( d2thet * sn(ise)*detdn*detdt
     .       + dthet * (detdt + sn(ise)*d2etdnt)
     .     - (thet+sn(ise)*dthet*detdn)*sn(ise)*dthet*detdt/znt)/znt )
c
      d2xdnv(ise) = dxdv*dxdn(ise)/x + x * (
     .  -1.5d0*frat* (d2etdnv + (thet - 1.5d0*frat)*detdn*detdv)
     .  +0.5d0*( d2thet * sn(ise)*detdn*detdv
     .       + dthet * (detdv + sn(ise)*d2etdnv)
     .     - (thet+sn(ise)*dthet*detdn)*sn(ise)*dthet*detdv/znt)/znt )
c
      d2xdt2 = dxdt**2/x + x * ( 1.5d0/t**2
     .     + d2etdt2 * ( thet - 1.5d0*frat + 0.5d0*sn(ise)*dthet/znt )
     .       + detdt**2 * ( dthet - 1.5d0*frat*(thet - 1.5d0 * frat)
     .         + 0.5d0*sn(ise)*(d2thet - sn(ise)*dthet**2/znt)/znt ) )
c
      d2xdtv = dxdt*dxdv/x + x * (
     .       d2etdtv * ( thet - 1.5d0*frat + 0.5d0*sn(ise)*dthet/znt )
     .    + detdt * detdv * ( dthet - 1.5d0*frat*(thet - 1.5d0 * frat)
     .         + 0.5d0*sn(ise)*(d2thet - sn(ise)*dthet**2/znt)/znt ) )
c
      d2xdv2 = dxdv**2/x + x * ( 0.5d0/vol**2
     .     + d2etdv2 * ( thet - 1.5d0*frat + 0.5d0*sn(ise)*dthet/znt )
     .     + detdv**2 * ( dthet - 1.5d0*frat*(thet - 1.5d0 * frat)
     .         + 0.5d0*sn(ise)*(d2thet - sn(ise)*dthet**2/znt)/znt ) )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tau and its derivatives                                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c..      if ( x .gt. 1.0d-1 ) then
c
c  test for not including tau correction
c
      if(notau.eq.1) then
	tau   = 1.d0
	dtau  = 0.d0
	d2tau = 0.d0
c
      else
c
c
        if ( x .gt. xtautr ) then
            dpx =   x
            dp1 =   3.d0 * ( log(1.d0 + dpx)
     .                     - dpx*(1.d0 - 0.5d0*dpx) ) / dpx**3
            dp2 = - 3.d0 * ( dp1 - 1.d0/(1.d0 + dpx)    ) / dpx
            dp3 =   - ( 4.d0*dp2 + 3.d0/(1.d0 + dpx)**2 ) / dpx
            tau   = dp1
            dtau  = dp2
            d2tau = dp3
        else
c..          tau  = ((((((((x/11. - 1.d0/10.)*x +  1.d0/9.)*x - 1.d0/8.)*x
c..     .                         + 1.d0/7. )*x -  1.d0/6.)*x + 1.d0/5.)*x
c..     .                         - 1.d0/4. )*x +  1.d0/3.)*3.d0
c..          dtau = (((((((8.d0*x/11. - 7.d0/10.)*x +2.d0/3.)*x-5.d0/8.)*x
c..     .                              + 4.d0/7.)*x -1.d0/2.)*x+2.d0/5.)*x
c..     .                              - 1.d0/4. )*3.d0
c..          d2tau= ((((((56.d0*x/11. - 21.d0/5.)*x+10.d0/3.)*x-5.d0/2.)*x
c..     .                             + 12.d0/7.)*x-1.d0)*x+2.d0/5.)*3.d0
c
c  test for setting coefficients
c
          if(initcf.eq.0) then
c
c  set coefficients for expansion
c
            isg=6*mod(nmax,2)-3
            do n=nmax,3,-1
              ctau(n)  =dfloat(isg)/dfloat(n)
              cdtau(n) =dfloat(isg*(n-3))/dfloat(n)
              cd2tau(n)=dfloat(isg*(n-3)*(n-4))/dfloat(n)
              isg=-isg
            end do
            initcf = 1
          end if
c
c  do expansion as do loop, to allow arbitrary high order
c
          tau=ctau(nmax)
          do n=nmax-1,3,-1
            tau=tau*x+ctau(n)
          end do
c
          dtau=cdtau(nmax)
          do n=nmax-1,4,-1
            dtau=dtau*x+cdtau(n)
          end do
c
          d2tau=cd2tau(nmax)
          do n=nmax-1,5,-1
            d2tau=d2tau*x+cd2tau(n)
          end do
        end if
       end if
c
c..      write(6,*) 'x, tau, etc =',x,tau, dtau, d2tau
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     free energy, pressure, internal energy                           c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c..      write(6,*) 'x, tau, dtau',x, tau,dtau
c     free energy
      f(4) = - cf4 * tau * dsqrt( znt**3 / (t * vol) )
c
c     pressure
      p(4) = f(4) * ( 0.5d0/vol-1.5d0* sn(ise) * dthet * detdv / znt
     .                          - dxdv * dtau / tau )
c
c     internal energy
      e(4) = f(4) * ( 1.5d0/t -1.5d0* sn(ise) * dthet * detdt / znt
     .                          - dxdt * dtau / tau ) * t
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     derivatives of f4                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c-----------------------------------------------------------------------
c     first derivatives                                                c
c-----------------------------------------------------------------------
c
      df4dt =   f(4) * ( -0.5d0/ t + 1.5d0*sn(ise) * dthet * detdt / znt
     .                            + dxdt * dtau / tau )
      df4dv = - p(4)
c
      do 7 is = 1, nspes - 1
      fscr(is) = f(4)*(1.5d0*zsq(is)/znt + dxdn(is)*dtau/tau)
      dfdn(is) = fscr(is)
    7 continue
      fscr(ise) = f(4)*( 1.5d0* ( thet + sn(ise)*dthet*detdn )/znt
     .                   + dxdn(ise) * dtau/tau )
      dfdn(ise) = fscr(ise)
c
c-----------------------------------------------------------------------
c     second derivatives                                               c
c-----------------------------------------------------------------------
c
      do 9 js = 1, nspes - 1
c
      if ( zs(js)  .ne. 0.d0 ) then
            do 8 is = 1, js
            d2fdn2(is, js) = ( fscr(is) * fscr(js) / f(4)
     .            + ( (d2tau/tau - (dtau/tau)**2) * dxdn(is) * dxdn(js)
     .            + d2xdn2(is, js) * dtau / tau
     .            - 1.5d0* zsq(is) * zsq(js) / znt**2 ) * f(4) )
    8       continue
            if( js .gt. 1 ) call scopy( js - 1, d2fdn2(1 ,js), 1,
     .                                          d2fdn2(js, 1), mspes )
c
            d2fdnt(js) = fscr(js) * df4dt / f(4)
     .              + ( (d2tau/tau - (dtau/tau)**2) * dxdn(js) * dxdt
     .              + d2xdnt(js) * dtau/tau
     .              - 1.5d0*zsq(js)*sn(ise)*dthet*detdt/znt**2)*f(4)
c
            d2fdnv(js) = fscr(js) * df4dv / f(4)
     .              + ( (d2tau/tau - (dtau/tau)**2) * dxdn(js) * dxdv
     .              + d2xdnv(js) * dtau/tau
     .              - 1.5d0*zsq(js)*sn(ise)*dthet*detdv/znt**2)*f(4)
      end if
    9 continue
c
      do 10 is = 1, nspes - 1
      d2fdn2(is, ise) = ( fscr(is) * fscr(ise) / f(4)
     .    + f(4)*( (d2tau/tau - (dtau/tau)**2) * dxdn(is) * dxdn(ise)
     .    + d2xdn2(is, ise) * dtau / tau
     .    -1.5d0*zsq(is)*(thet+sn(ise)*dthet*detdn)/znt**2))
   10 continue
      call scopy( nspes - 1, d2fdn2(1,ise), 1, d2fdn2(ise,1), mspes )
c
      d2fdn2(ise, ise) = fscr(ise)**2 / f(4)
     .    + f(4) * ( (d2tau/tau - (dtau/tau)**2) * dxdn(ise)**2
     .    + d2xdn2(ise, ise) * dtau / tau
     .    + 1.5d0*( d2thet * sn(ise) * detdn**2
     .    + dthet * ( 2.d0*detdn + sn(ise) * d2etdn2 )
     .    - (thet + sn(ise)*dthet*detdn)**2/znt) / znt )
c
      d2fdnt(ise) = fscr(ise)*df4dt/f(4) + f(4) *
     . ( (d2tau/tau - (dtau/tau)**2)*dxdt*dxdn(ise)
     . + d2xdnt(ise) * dtau/tau
     . + 1.5d0 *( d2thet*sn(ise)*detdt*detdn
     . +dthet*(detdt + sn(ise)*d2etdnt)
     . -(thet+sn(ise)*dthet*detdn)*sn(ise)*dthet*detdt/znt)/znt )
c
      d2fdnv(ise) = fscr(ise)*df4dv/f(4) + f(4) *
     . ( (d2tau/tau - (dtau/tau)**2)*dxdv*dxdn(ise)
     . + d2xdnv(ise) * dtau/tau
     . + 1.5d0 *( d2thet*sn(ise)*detdv*detdn
     . +dthet*(detdv + sn(ise)*d2etdnv)
     . -(thet+sn(ise)*dthet*detdn)*sn(ise)*dthet*detdv/znt)/znt )
c
      d2fdtv = df4dt*df4dv/f(4) + f(4) *
     . ( (d2tau/tau - (dtau/tau)**2)*dxdt*dxdv
     . + d2xdtv * dtau/tau
     . + 1.5d0 * sn(ise) * ( d2thet*detdt*detdv
     . + dthet*(d2etdtv-sn(ise)*dthet*detdt*detdv/znt))/znt )
c
      d2fdt2 = df4dt**2/f(4) + f(4) *
     . ( 0.5d0 / t**2
     . + (d2tau/tau - (dtau/tau)**2)*dxdt**2
     . + d2xdt2 * dtau/tau
     . + 1.5d0 * sn(ise) * ( d2thet*detdt**2
     . + dthet*(d2etdt2-sn(ise)*dthet*detdt**2/znt))/znt )
c
      d2fdv2 = df4dv**2/f(4) + f(4) *
     . ( 0.5d0 / vol**2
     . + (d2tau/tau - (dtau/tau)**2)*dxdv**2
     . + d2xdv2 * dtau/tau
     . + 1.5d0 * sn(ise) * ( d2thet*detdv**2
     . + dthet*(d2etdv2-sn(ise)*dthet*detdv**2/znt))/znt )
c
      return
      end
