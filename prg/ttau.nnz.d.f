      subroutine ttau(tau,ams,ars,als,x,z,t,qhopf,dqhopf,
     *  dtrs,dtls,dttau)
c  for given tau, and mass (ams), radius (ars), luminosity (als)
c  and composition (x,z) sets temperature and dtrs = d log t/d log rs,
c  dtls = d log t/d log ls, dttau = d log t/d log tau, using dog analyti
c  fit to hsra
c
c  modified 15/3/88 to use consistent numerical constants
c
c  Modified 10/6/98, redefining qatmos for consistency with RT 
c  version of routine.
c
c  Modified 7/4/01, adding second derivative of qhopf in common/cqhopf/
c
      implicit double precision (a-h,o-z)
      common/qatmos/ iqfit,iqqbc,qq(7)
      common/cqhopf/ d2qhpf
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev, syear
c
c  common giving derived constants, and constants used for
c  equation of state
c
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev1,exhm
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3,idgmxc,
     *  idgsdt
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      fxp(a)=exp(min(85.d0,max(a,-85.d0)))
c
      idiag=0
c
      if(tau.lt.0) then
	if(idiag.gt.0) write(istdou,110) tau
	taul=0.d0
      else
	taul=tau
      end if
c
      e1=qq(2)* fxp(-qq(3)*taul)
      e2=qq(4)* fxp(-qq(5)*taul)
      qhopf=qq(1)+e1+e2
      dqhopf=-qq(3)*e1-qq(5)*e2
c..      if(taul.le.-10) write(6,*) 'dqhopf,qq(3),e1,qq(5),e2',
c..     *  dqhopf,qq(3),e1,qq(5),e2
      d2qhpf=qq(3)*qq(3)*e1+qq(5)*qq(5)*e2
      if(ars.le.0) return
c
      tt=taul+qhopf
      cteff=0.75/(pi*clight*car)
      t=(cteff*als*tt/(ars*ars))**0.25d0
      dtrs=-0.5d0
      dtls=0.25d0
      dttau=taul*(1+dqhopf)/(4*tt)
      return
  110 format(/' ***** Warning. tau =',1pe13.5,' negative in ttau'/
     *        '       Replace by tau = 0')
      end
