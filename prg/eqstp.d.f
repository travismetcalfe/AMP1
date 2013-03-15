      subroutine eqstp(pl, tl, xh, yh, zh, nosd, notd, fl, nit)
c
c  equation of state routine, with independent variables
c    pl = log(p)
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
c  If fl .le. -1.d10 on input, a trial value is set in the routine.
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
c
      save
c
      data eps /1.e-10/
c
c  test for setting trial
c
      if(fl.le.-1.d10) then
        icase = 1
        call inteff(tl,pl,rhl,flini,wt,icase,iextr)
        fl = wt(1)*flini(1) + wt(2)*flini(2)
c
      end if
c
c  start loop for iteration
c
      nit=0
   10 continue
c..      write(6,*) 'In eqstp, call eqstf with fl, tl, xh =',fl, tl, xh
      call eqstf(fl, tl, xh, yh, zh, nosd, notd)
c
      pli=log10(p(1))
      dfl=(pl - pli)/p(2)
      nit=nit+1
c
c  limit size of dfl
c
      if(abs(dfl).gt.0.4) dfl=sign(0.4d0,dfl)
      if(idgeos.ge.1.and.istdpr.gt.0) write(istdpr,*) 
     *  nit,pl, tl, fl, dfl
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
        write(istdou,110) pl, tl, fl, dfl
        if(istdpr.gt.0.and.istdpr.ne.istdou)
     *    write(istdpr,110) pl, tl, fl, dfl
        fl = -1.d11
        go to 20
c
      end if
c
c  this is the end
c
   20 continue
c..      write(6,*) 'converged fl =',fl
      return
  110 format(//' ***** Iteration failed in s/r eqstp.'/
     *  7x,'log(p) =',f10.5,'  log(T) =',f10.5,
     *  '  Last log(f) =',f10.5,'  dfl =',1pe13.5)
      end
