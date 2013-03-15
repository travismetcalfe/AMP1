      program lmmin 

c     x = initial guess
c     yy = observations
c     sig = errors
c     ia = 1 = fit, 0 = no fit

      implicit double precision (a-h,o-z)

      double precision dat(36),data(36),x_initial, yi(10000), sig(10000)
      double precision one, y_initial, z_initial, chisq, retx(7)
c      double precision rng(3)
      integer npar,IARGC,i, nsel, sub(10000), ia(36), nobs, nip, nse
      integer ell(10000)
      character*80 par_string

c     common things
      integer nfreq, nonseis,subscripts(10000),iwflag,isflag,idif
      double precision rng(3), f0
      common/olcdint/ nfreq,nonseis,subscripts,iwflag,isflag,idif
      common/olcddbl/ rng


c     get data from the file.... these are initial parameters + data
      call mrqdata (ia, npar, dat, nse, yi, sig, ell)
c      print *,'out of mrqdata'
c      stop
c     get the subscripts that we need
c      call rd_subs (nsel, sub)   

c     if we give initial conditions at beginning, use these
c      nip = IARGC()
c      if (nip.GE.6) then
c         do i=1,nip
c            CALL GETARG(i,par_string)
c            READ(par_string,*) dat(i)
c         enddo
c      endif
 
c
c     Change to read in the initial He abundance
c
c      one = 1.0
c      y_initial = dat(2)
c      z_initial = dat(3)
c      x_initial = one - y_initial - z_initial
c      dat(2) = x_initial

c     nobs will be dimensiton of YMOD array, we only use ndata of these
      nobs = 1200

c     MINIMIZE HERE:c
c     INPUT: DAT = parameters, YI = observed data, sig = errors, 
c     ia= parameters to fit (0,1), SUB = subscripts of observables to use (#)
c     npar = n_parameters, nobs = total number of observables

      call reducemrq(dat, yi,sig,ia,subscripts, npar, nobs, chisq, retx)


      end

