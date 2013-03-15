      subroutine setovs(x,y,nn,iy,incovs)
c
c  Set parameters required for treatment of convective overshoot,
c  from convection-zone parameters defined in s/r rhs.
c  If incovs = 1, initialize parameters suitably, for no overshoot
c  (this may be changed later, by setting convection-zone parameters
c  for model)
c
c  qmxove, rmxove and nmxove are set to the innermost point in the
c  envelope overshoot region
c  qmxovc, rmxovc and nmxovc are set to the outermost point in the
c  core overshoot region
c
c  Original version: 5/1/00.
c
      implicit double precision (a-h, o-z)
      logical concor
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
c
      dimension x(1), y(iy,1)
c
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,jnc,frcf(6),frcl(6),
     *  xrcf(6),xrcl(6)
      common/convvr/ cvvarf(icvvar,6), cvvarl(icvvar,6)
      common/covsec/ icsove, icsovc, alpove, alpovc, 
     *  qmxove, qmxovc, rmxove, rmxovc, imxove, imxovc, 
     *  nmxove, nmxovc
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      data idiag /1/
c
      if(istdpr.le.0) idiag=0
c
c  Test for initialization
c
      if(incovs.eq.1.or.inc.eq.0) then
	if(istdpr.gt.0) write(istdpr,110)
	qmxove = -1
	qmxovc = -1
	rmxove = -1
	rmxovc = -1
	nmxove = -1
	nmxovc = -1
	imxove = -1
	imxovc = -1
	return
      end if
c
c  set flag for convective core
c
      concor = nl(inc).eq.nn
c
      rs=1.d11*10.d0**y(1,1)
c
c  test for initializing parameters for envelope overshoot
c
      if(icsove.lt.1.or.alpove.le.0) then
	qmxove = -1
	rmxove = -1
	nmxove = -1
	imxove = -1
c
      else
c
c  start initialization
c
c  set relevant convection zone (the last that is not the core)
c
	if(concor) then
	  imxove = inc-1
        else
	  imxove = inc
	end if
c
c  radius (in cm) at the edge of the overshoot region
c
	rmxove = rs*cvvarl(2,imxove)-alpove*cvvarl(7,imxove)
	rmxovl = log10(rmxove/1.d11)
c
c  locate point
c
	do 20 n=nl(imxove)-1,nn
	if(y(1,n).lt.rmxovl) go to 25
	nmxove = n
   20   continue
c
   25   qmxove = 10.d0**x(nmxove)
c
	if(idiag.eq.1) write(istdpr,115) nmxove, rmxove, qmxove
c
      end if
c
c  test for initializing parameters for core overshoot
c
      if(icsovc.lt.1.or.alpovc.le.0.or..not.concor) then
	qmxovc = -1
	rmxovc = -1
	nmxovc = -1
	imxovc = -1
c
      else
c
c  start initialization
c
c  set relevant convection zone 
c
	imxovc = inc
c
c  set scale of overshoot, with warning if core is small
c
	rccore=rs*cvvarf(2,imxovc)
	hpcore=cvvarf(7,imxovc)
	if(rccore.lt.hpcore) then
	  if(istdpr.gt.0) write(istdpr,120) rccore, hpcore
	  sccore = rccore
        else
	  sccore = hpcore
        end if
c
c  radius (in cm) at the edge of the overshoot region
c
	rmxovc = rccore + alpovc*sccore
	rmxovl = log10(rmxovc/1.d11)
c
c  locate point
c
	do 30 n=nf(imxovc)+1,1,-1
	if(y(1,n).gt.rmxovl) go to 35
	nmxovc = n
   30   continue
c
   35   qmxovc = 10.d0**x(nmxovc)
c
	if(idiag.eq.1) write(istdpr,125) nmxovc, rmxovc, qmxovc
c
      end if
c
      return
c
  110 format(//' Initialize overshoot parameters in s/r setovs',
     *         ' for initial model, to no overshoot'/)
  115 format(/' In setovs, nmxove, rmxove, qmxove =',
     *        i5, 1pe13.5,0pf12.7)
  120 format(/' ***** Warning in s/r setovs. rcore =',1pe13.5,' le ',
     *        ' Hp =',e13.5/
     *        '       Overshoot scale set to rcore.'/)
  125 format(/' In setovs, nmxovc, rmxovc, qmxovc =',
     *        i5, 1pe13.5,0pf12.7)
      end
