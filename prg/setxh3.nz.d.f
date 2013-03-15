      subroutine setxh3(x,y,xhe3,nn,iy,ix3,iche3)
c
c  set he3 abundance throughout solar model. if iche3
c  (in argument list) is 2, equilibrium value is set.
c  otherwise, if ifdhe3 (in common/he3fdg/) is 1, the
c  fudged abundance for age agesh is used.
c
c  original version: 29/12/1984
c
c  Modified 19/7/95, to use general heavy-element abundance
c  as set in common/heavy/. 
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
      implicit double precision (a-h,o-z)
      logical noder
      dimension x(1),y(iy,1),xhe3(ix3,1),xhe3fg(4)
c
      common/he3fdg/ agesh,ifdhe3
      common/heavy/ zatmos, zhc, zh(1)
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
c  test for case
c
   10 if(iche3.eq.2) go to 20
c
      if(ifdhe3.ne.1) go to 90
c
c  fudge. use age in agesh
c
      ageact=agesh
      if(istdpr.gt.0) write(istdpr,100)
      go to 30
c
c  he3 in equilibrium. set age to -1 as flag.
c
   20 ageact=-1.
      if(istdpr.gt.0) write(istdpr,110)
c
c  step through model
c
   30 noder=.true.
c
      do 50 n=1,nn
      fl=y(2,n)
      tl=y(3,n)
      xh=y(5,n)
      yh=1-xh-zh(n)
      idghe3=0
c
      call eqstf(fl,tl,xh,yh,zh(n),noder,noder)
      call he3abd(fl,tl,xh,yh,zh(n),ageact,xhe3fg,anu,noder)
c
   50 xhe3(1,n)=xhe3fg(1)
c
      idghe3=0
c
      return
c
c  diagnostics for incorrect values of iche3 and ifdhe3
c
   90 write(istdou,120) iche3,ifdhe3
      if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,120) 
     *  iche3,ifdhe3
      return
  100 format(//' set fudged he3 abundance in model')
  110 format(//' set equilibrium he3 abundance in model')
  120 format(//' ***** warning. s/r setxh3 called with iche3 =',
     *  i5,'  ifdhe3 =',i5/14x,' no action taken')
      end
