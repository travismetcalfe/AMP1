      subroutine eqstrh(rhl, tl, xh, yh, zh, nosd, notd, fl, nit)
c
c  equation of state routine, with independent variables
c    rhl = log(rho)
c    tl  = log(T)
c
c  Iterates, using eqstf, to determine log(f), and hence set equation
c  of state variables. 
c  log(f) is returned in fl.
c  nit returns number of iterations.
c
c  Accuracy is set in eps below, and should be adjusted for
c  computer.
c
c  If fl .le. -1.e10 on input, a trial value is set in the routine.
c  Otherwise the value of fl provided is used as trial.
c
c  As for s/r eqstf, s/r setcns must be called before calling
c  eqstrh.
c
c  Note: to obtain the trial value of fl, a fit to mue**-1, as a
c  function of log T, is used, in the region of partial ionization.
c  this is currently based on solar composition and
c  log T - log rho relation, and hence may give problems under
c  different circumstances.
c
c  Original version: 07/08/87
c
c  This version is set up for use with Coulomb term version
c  of eqstf. When starting from trial value of fl, make initial
c  iteration for fl without Coulomb terms.
c
c  Note that present this is the case regardless of input fl
c
c  Version from 1/5/90.
c
c  Modified 21/6/03: Include common/cdgphs/ kdgeos, kdgopc, kdgeng.
c  If, on input, kdgeos .gt. 0, do not stop on fatal error, but
c  return with kdgeos = -1.
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 9/3/90
c
c
      implicit double precision(a-h,o-z)
      logical nosd, notd
      dimension flini(2), wt(2)
      common/eqstd/ xii(4), ane(10), rho(20), ht(20), p(20), cp(4),
     *  dad(4), dlt(4), gm1, tprh, trhp, rhxp
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
      common/ccoulm/ epsdmu, icoulm, iclmit, iclmsd, nitdmu, epssdr
c
      common/cdgphs/ kdgeos, kdgopc, kdgeng
c
      save
c
      data eps /1.e-10/
c..      data kdgeos, kdgopc, kdgeng /0, 0, 0/
c
c  store original value of Coulomb case flag
c
      icoulp=icoulm
c
c  test for setting trial
c
c
c  try to use wd setting of trial fl in all cases
c
      icase = 2
      icoulm=-1
      call inteff(tl,pgl,rhl,flini,wt,icase,iextr)
      fl = wt(1)*flini(1) + wt(2)*flini(2)
c..      write(6,*) 'tl, rhl =', tl, rhl
c..      write(6,*) 'trial fl set to', fl
      if(fl.le.-1.e10) then
c..        xt=max(0.d0,tl-3.76119)
c..        xt=xt*xt
c..        xmufl=3.8*xt/(xt+0.035)-3.83015
c..        fl=7.829+xmufl+rhl-1.5*tl
c
      end if
c
c  start loop for iteration
c
      nit=0
   10 call eqstf(fl, tl, xh, yh, zh, nosd, notd)
c
      if(rho(1).le.0.or.rho(2).eq.0) then
        write(istdou, 105) fl, tl, rhl, rho(1)
        if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr, 105) 
     *    fl, tl, rhl, rho(1)
        if(istdpr.gt.0) call dmpeqs
	if(kdgeos.gt.0) then
	  kdgeos=-1
	  return
        else
	  stop 'eqstrh'
        end if
      end if
      rhli=log10(rho(1))
      dfl=(rhl - rhli)/rho(2)
      nit=nit+1
c
c  limit size of dfl
c
      if(abs(dfl).gt.0.4) dfl=sign(0.4d0,dfl)
c
c  stupid statement added to make the code run on tuc47, 6/12/04
c
      dflstr=dfl
      if(idgeos.ge.1.and.istdpr.gt.0) write(istdpr,*) 
     *  nit,rhl, tl, fl, dfl
c
c  test for continuing iteration
c
      if(abs(dfl).lt.eps) then
        go to 20
      else if(nit.le.60) then
        fl = fl+dfl
        go to 10
      else
c
c  diagnostics for failed convergence
c
        write(istdou,110) rhl, tl, fl, dfl
        if(istdpr.ne.istdou.and.istdpr.gt.0)
     *    write(istdpr,110) rhl, tl, fl, dfl
        fl = -1.e11
        go to 20
c
      end if
c
c  this is the end
c
   20 continue
c
c  test for repeat with Couilomb term
c
      if(icoulm.ne.icoulp) then
        icoulm=icoulp
        nit=0
        go to 10
      end if
c..      write(6,*) 'converged fl =',fl
      return
  105 format(//' **** error in eqstrh. rho .le. 0 or rho(2) = 0 for'/
     *  ' log f =',f10.5,'  log T =',f10.5,' log rho =', f10.5,
     *  ' rhoi =',1pe13.5)
  110 format(//' ***** Iteration failed in s/r eqstrh.'/
     *  7x,'log(rho) =',f10.5,'  log(T) =',f10.5,
     *  '  Last log(f) =',f10.5,'  dfl =',1pe13.5)
      end
      subroutine inteff(tl,pgl,rl,flini,wt,icase,iextr)
c>>>>>
c>>>>> finds starting value of flini for eff routines if the argument
c>>>>> is either log10(gas pressure) or log10(density). iextr is
c>>>>> set to 1 if extrapolation outside the table had to be made.
c>>>>>
c>>>>> based on eff-pre-computation with x=.73,y=.25,z=.02
c>>>>>
c>>>>> contains no common statements. parameters are local.
c>>>>>
c===== icase = 1: input=(tl,pgl),  rl ignored
c===== icase = 2: input=(tl,rl ), pgl ignored
c
      implicit double precision(a-h,o-z)
      parameter (nm = 11,nml = 5,nmh = 6,ntl = nml*nm,nth = nmh*nm)
c
      dimension tlog(nm),rlog(nm,nm),pglog(nm,nm),flog(nm,nm)
      dimension rll(ntl),rlh(nth),pgll(ntl),pglh(nth),fll(ntl),flh(nth)
      dimension flini(2),wt(2)
c
      equivalence(rlog (1,1),rll (1)) ,(rlog (1,nmh),rlh (1))
      equivalence(pglog(1,1),pgll(1)) ,(pglog(1,nmh),pglh(1))
      equivalence(flog (1,1),fll (1)) ,(flog (1,nmh),flh (1))
c
      data tlog/
     .    3.40000,   3.80000,   4.20000,   4.60000,   5.00000,
     .    5.40000,   5.80000,   6.20000,   6.60000,   7.00000,
     .    7.40000/
c
      data rll/
     .  -15.93806, -14.59077, -13.11393, -11.71751, -10.38835,
     .   -8.97767,  -7.59375,  -6.16023,  -4.74661,  -3.34048,
     .   -1.93629,
     .  -14.39806, -13.08988, -11.78384, -10.46730,  -9.15835,
     .   -7.85767,  -6.54375,  -5.24023,  -3.92661,  -2.62048,
     .   -1.31627,
     .  -12.79806, -11.57078, -10.36165,  -9.14343,  -7.92835,
     .   -6.71767,  -5.50375,  -4.28023,  -3.06660,  -1.85045,
     .   -0.63617,
     .  -11.19806, -10.06214,  -8.95086,  -7.82325,  -6.69833,
     .   -5.57766,  -4.45374,  -3.33021,  -2.20654,  -1.08028,
     .    0.04407,
     .   -9.59806,  -8.56014,  -7.52659,  -6.49503,  -5.46791,
     .   -4.43760,  -3.40366,  -2.37003,  -1.33614,  -0.30961,
     .    0.72372/
      data rlh/
     .   -7.99806,  -7.05242,  -6.11177,  -5.17967,  -4.23200,
     .   -3.29677,  -2.35276,  -1.41861,  -0.47419,   0.46015,
     .    1.40286,
     .   -6.39806,  -5.55147,  -4.69828,  -3.84877,  -3.00142,
     .   -2.15675,  -1.30499,  -0.46287,   0.38404,   1.23722,
     .    2.08971,
     .   -4.79807,  -4.04144,  -3.28452,  -2.53132,  -1.76802,
     .   -1.01603,  -0.25685,   0.49360,   1.25430,   2.00696,
     .    2.76640,
     .   -3.19831,  -2.53173,  -1.86318,  -1.20584,  -0.53947,
     .    0.12553,   0.78489,   1.45261,   2.11985,   2.77988,
     .    3.44720,
     .   -1.59748,  -1.02039,  -0.44723,   0.12189,   0.69032,
     .    1.26763,   1.83971,   2.41156,   2.98379,   3.55092,
     .    4.12641,
     .    0.00202,   0.48728,   0.96517,   1.44458,   1.92401,
     .    2.40470,   2.88056,   3.36421,   3.84159,   4.32331,
     .    4.80062/
c
      data pgll/
     .   -4.71634,  -2.70110,  -0.79422,   1.01898,   2.74815,
     .    4.55883,   6.34275,   8.17627,   9.98989,  11.79601,
     .   13.60020,
     .   -3.17634,  -1.30952,   0.53582,   2.26909,   3.97815,
     .    5.67883,   7.39275,   9.09627,  10.80989,  12.51602,
     .   14.22022,
     .   -1.57634,   0.09093,   1.95691,   3.59094,   5.20815,
     .    6.81883,   8.43275,  10.05627,  11.66990,  13.28603,
     .   14.90028,
     .    0.02366,   1.56549,   3.35740,   4.90076,   6.43816,
     .    7.95883,   9.48275,  11.00628,  12.52992,  14.05612,
     .   15.58052,
     .    1.62366,   3.06194,   4.76523,   6.22483,   7.66836,
     .    9.09886,  10.53279,  11.96637,  13.40014,  14.82664,
     .   16.26059/
      data pglh/
     .    3.22366,   4.56932,   6.08520,   7.53751,   8.90120,
     .   10.23927,  11.58325,  12.91715,  14.26164,  15.59761,
     .   16.94023,
     .    4.82375,   6.07052,   7.37716,   8.83489,  10.11675,
     .   11.37463,  12.62823,  13.87314,  15.12352,  16.37610,
     .   17.63041,
     .    6.42746,   7.58888,   8.76598,  10.05827,  11.32914,
     .   12.52233,  13.69271,  14.83754,  16.00047,  17.15703,
     .   18.32207,
     .    8.15380,   9.30678,  10.50206,  11.73099,  12.68725,
     .   13.70433,  14.76934,  15.84437,  16.91982,  17.98796,
     .   19.06506,
     .   10.79939,  11.94400,  13.08449,  13.22093,  14.17459,
     .   15.14145,  16.10082,  17.05997,  18.01959,  18.97151,
     .   19.93379,
     .   13.96914,  14.93979,  14.51362,  15.31529,  16.11750,
     .   16.92229,  17.71977,  18.53028,  19.33067,  20.13707,
     .   20.93387/
c
      data fll/
     .  -15.20875, -12.62213, -11.68335, -10.85416, -10.12500,
     .   -9.31434,  -8.53047,  -7.69708,  -6.88379,  -6.07850,
     .   -5.27639,
     .  -13.66875, -11.40213, -10.35335,  -9.60416,  -8.89500,
     .   -8.19434,  -7.48047,  -6.77708,  -6.06379,  -5.35850,
     .   -4.65639,
     .  -12.06875, -10.50213,  -8.93335,  -8.28416,  -7.66500,
     .   -7.05434,  -6.44047,  -5.81708,  -5.20379,  -4.58850,
     .   -3.97639,
     .  -10.46875,  -9.61213,  -7.54335,  -6.98416,  -6.43500,
     .   -5.91434,  -5.39047,  -4.86708,  -4.34379,  -3.81850,
     .   -3.29639,
     .   -8.86875,  -8.40213,  -6.15335,  -5.66416,  -5.20500,
     .   -4.77434,  -4.34047,  -3.90708,  -3.47379,  -3.04850,
     .   -2.61639/
      data flh/
     .   -7.26875,  -6.92213,  -4.97335,  -4.35416,  -3.97500,
     .   -3.63434,  -3.29047,  -2.95708,  -2.61379,  -2.27850,
     .   -1.93639,
     .   -5.66875,  -5.42213,  -4.11335,  -3.09416,  -2.77500,
     .   -2.50434,  -2.25047,  -2.00708,  -1.75379,  -1.49850,
     .   -1.24639,
     .   -4.06875,  -3.91213,  -3.30335,  -2.10416,  -1.64500,
     .   -1.40434,  -1.20047,  -1.03708,  -0.87379,  -0.71850,
     .   -0.55639,
     .   -2.46875,  -2.40213,  -2.26335,  -1.25416,  -0.26500,
     .   -0.17434,  -0.11047,  -0.03708,   0.03621,   0.10150,
     .    0.17361,
     .   -0.85875,  -0.88213,  -0.90335,   1.24584,   1.20500,
     .    1.17566,   1.13953,   1.10292,   1.06621,   1.02150,
     .    0.98361,
     .    0.92125,   0.77787,   3.15665,   2.99584,   2.83500,
     .    2.67566,   2.50953,   2.35292,   2.18621,   2.02150,
     .    1.84361/
c
      it=-999
      ip=-999
      ir=-999
      iextr=0
c========== select isotherms for linear inter(extra)polation
      do 11 i=1,nm
      if(tl.lt.tlog(i)) goto 20
      it=i
 11   continue
         iextr=1
         it=nm-1
c
  20  if(it.eq.-999) then
         iextr=1
         it=1
      end if
c========== tl-part of arguments
      it1=it+1
      x0 =tlog(it)
      x1 =tlog(it1)
      x2 =x0
      x3 =x1
      x  =  tl
c
c................................. icase = 1 ............................
      if(icase.eq.2) goto 200
c
c========== select pressure points for linear inter(extra)polation
      do 21 i=1,nm
      if(pgl.lt.pglog(it,i)) goto 30
      ip=i
  21  continue
         iextr=1
         ip=nm-1
c
  30  if(ip.eq.-999) then
         iextr=1
         ip=1
      end if
c
c========== define the three function values for inter(extra)polation
      ip1=ip+1
      y0 =pglog(it,ip)
      y1 =pglog(it1,ip)
      y2 =pglog(it,ip1)
      y3 =pglog(it1,ip1)
      z0 =flog(it,ip)
      z1 =flog(it1,ip)
      z2 =flog(it,ip1)
      z3 =flog(it1,ip1)
c
c========== define arguments
      y  = pgl
      goto 1000
c................................. icase = 2 ............................
c
c========== select density points for linear inter(extra)polation
 200  do 31 i=1,nm
      if(rl.lt.rlog(it,i)) goto 40
      ir=i
  31  continue
         iextr=1
         ir=nm-1
c
  40  if(ir.eq.-999) then
         iextr=1
         ir=1
      end if
c
c========== define the three function values for inter(extra)polation
      ir1=ir+1
      y0 =rlog(it,ir)
      y1 =rlog(it1,ir)
      y2 =rlog(it,ir1)
      y3 =rlog(it1,ir1)
      z0 =flog(it,ir)
      z1 =flog(it1,ir)
      z2 =flog(it,ir1)
      z3 =flog(it1,ir1)
c
c========== define argument
      y  =  rl
c
c========== call bilinear interpolation
c
1000  continue
c
c....... lower triangle
      call  bilin (x0,y0,z0,x1,y1,z1,x2,y2,z2,x,y,z)
      flini(1) = z
c
c....... upper triangle
      call  bilin (x3,y3,z3,x2,y2,z2,x1,y1,z1,x,y,z)
      flini(2) = z
c
c....... weights (quite arbitrary)
      wlow = 1./( (x-x0)**2 + (y-y0)**2 + 1.e-5 )
      whig = 1./( (x-x3)**2 + (y-y3)**2 + 1.e-5 )
      wtot = wlow + whig
      wlow = wlow/wtot
      whig = whig/wtot
      wt(1) = wlow
      wt(2) = whig
c
      return
      end
c
c
      subroutine bilin(x0,y0,z0,x1,y1,z1,x2,y2,z2,x,y,z)
c
c------- performs bilinear interpolation of the function f,
c------- given on three arbitray points (x0,y0),(x1,y1),(x2,y2),
c------- where the respective function values are z0,z1,z2.
c
      implicit double precision(a-h,o-z)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      x10=x1-x0
      x20=x2-x0
      y10=y1-y0
      y20=y2-y0
      z10=z1-z0
      z20=z2-z0
c
      det=x10*y20-y10*x20
c
      if (det.ne.0) then
c
        dzdx=(z10*y20-z20*y10)/det
        dzdy=(z20*x10-z10*x20)/det
c
        z = z0 + (x-x0)*dzdx + (y-y0)*dzdy
c
      else
        write(istdou,8000) x10,x20,y10,y20
        if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,8000) 
     *    x10,x20,y10,y20
        stop 'bilin'
      end if
      return
8000  format(/' collinear points in s/r bilin. error stop.',
     . ' x1-x0,x2-x0,y1-y0,y2-y0 = ',/1x,1p4g15.6/)
      end
      subroutine dmpeqs
c  dumps commons from s/r eqstf
      implicit double precision (a-h, o-z)
      common/eqstd/ c1(90)
      common/eqsout/ c2(210)
      common/dmuder/ c3(10)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      write(istdpr,105)
      write(istdpr,110) c1
      write(istdpr,115) c2
      write(istdpr,120) c3
      return
  105 format(///' output from dmpeqs:'/1x,20(1h*))
  110 format(/' common/eqstd/:'/1p4e13.5/(10e13.5))
  115 format(/' common/eqsout/:'/(1p10e13.5))
  120 format(/' common/dmuder/:'/1p10e13.5)
      end
