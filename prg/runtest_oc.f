      program runtest



      implicit double precision (a-h,o-z)

      double precision data(36),yy(10000), chi, k
c      double precision rng(3)
      integer npar,IARGC,i, nsel, subs(10000), nobs, sub(10000)
      character*80 par_string

      integer nfreq, nonseis,subscripts(10000),iwflag,isflag,idif
      double precision rng(3), f0, sys(1200), error, sig2i, dyi

      integer ia(36), nse, ell(10000)
      double precision dat(36), yi(10000), sig(10000), syserror(1200)
      common/olcdint/ nfreq,nonseis,subscripts,iwflag,isflag,idif
      common/olcddbl/ rng

      npar = IARGC()
      if (npar.LT.1) then
         write(*,*)
         write(*,*) "Usage: runevol <mass> <Y> <Z> <alpha> <age> <rot>"
         write(*,*)
	 stop
      endif
      do i=1,npar
         CALL GETARG(i,par_string)
         READ(par_string,*) data(i)
      enddo


      if (npar.eq.7) then 
         npar = npar - 1
      endif

      nobs = 1200
c
c     this is the subroutine that can be called from anywhere
c     i need a wrapper to call this!! 
c
c     get the data
      call mrqdata (ia, npar, dat, nse, yi, sig, ell)
c     run the model
      call evalmodel(npar, data, nobs, nse, yy)
      
      
      do i = 1, nobs
         sub(i) = 0
      enddo
      do i = 1, nse
         sub(subscripts(i)) = 1
      enddo

c     save the model observables to a file
      open (13, file = "model2", status = "OLD")
 2    format (i5, i5, d16.9)
      do i = 1, nobs
         write(13, 2) i, sub(i), yy(i)
      enddo
      close (13)

      call surfacesigi(nobs, syserror)
      k = 1
      do i = 1, nobs
         if (sub(i).ne.0) then 
            sys(i) = syserror(k)
            k = k + 1
         endif
      enddo

c
c     Calculate chi^2
c
      k = 1
      chi = 0.d0
      do i = 1, nse
            error = sqrt(sig(i)*sig(i)+0.25*syserror(i)*syserror(i))
            sig2i = 1.d0/error**2.
            dyi = yi(i) - yy(subscripts(i))
            chi = chi + dyi*dyi*sig2i
      enddo

      print *,'Reduced chi sqd',chi/nse
      end
 
