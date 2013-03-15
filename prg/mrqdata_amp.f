      subroutine mrqdata (ia, npar, xx, nobs, yy, sig, ell)


      implicit none
      integer ma, ia(36), i, sub(10000), nobs, x1, npar, ell(10000), x4
      integer col1int, nip
      double precision  yy(10000),sig(10000), x3, xx(36), x2
      character a*55, col1str*1, par_string*80
      integer nsel
      double precision  rng(3)

      real ntemp, n(50)
      integer ltemp, l(50)
      integer nonseis, nfreq, subscripts(10000), iwflag, isflag, idif
      character a1*10, a2*10
      double precision f0, timestep


      common/olcdint/ nfreq,nonseis,subscripts,iwflag,isflag,idif
      common/olcddbl/ rng
c
c  read obs.dat file
c
         nonseis=0
         open(55,file='obs.dat',status='old')
         read(55,*) nobs, iwflag, isflag, idif
         do i=1,nobs+5
            read(55,*,end=50) col1str,yy(i),sig(i)
            if (col1str .eq. "T") then
               ell(i)=5
c               T_obs=yy(i)
c               Terr=sig(i)
               nonseis=nonseis+1
            elseif (col1str .eq. "G") then
               ell(i)=6
c               G_obs=yy(i)
c               Gerr=sig(i)
               nonseis=nonseis+1
            elseif (col1str .eq. "M") then
               ell(i)=7
c               M_obs=yy(i)
c               Merr=sig(i)
               nonseis=nonseis+1
            elseif (col1str .eq. "L") then
               ell(i)=8
c               L_obs=yy(i)
c               Lerr=sig(i)
               nonseis=nonseis+1
            elseif (col1str .eq. "R") then
               ell(i)=9
c               R_obs=yy(i)
c               Rerr=sig(i)
               nonseis=nonseis+1
            else
               read(col1str,*) col1int
               ell(i)=col1int
            endif
         enddo
 50      close(55)
         nfreq = nobs
         nobs = nfreq + nonseis
         
      npar=5
      do i = 1, npar 
         ia(i) = 1
      enddo
    
c
c     I need to read the paramlog file to get the information 
c     and assign the indices values
c
c
      close(11)
 
      open ( 11, FILE = 'paramlog.001', STATUS = 'OLD' )
 3    format (a55)
c
c     read in the first set of lines until we get close to < 1 time step
c
      do i = 1, 7
         read (11, 3) a
      enddo
c
c     now read and check the size of the step. if it is less than abs(1) then
c     we can move on to reading params
c
      do i = 1, 10
         read ( 11, *) timestep
         if (abs(timestep).lt.1) goto 23
      enddo
      
 23   read (11, 3) a
c 
c     read in parameter
c
      read (11, *) a1, xx(1), xx(2), xx(3), xx(4), xx(5)
      do i = 1, 3 
         read (11, 3) a
      enddo
c      print *,xx(1), xx(2), xx(3), xx(4), xx(5)
c 
c     read in f0
c
      read (11, *) a1, a2

c
c     read in last 3 lines before frequency information
c
      do i = 1, 3
         read (11, 3) a
      enddo
c
c     read in the frequency information, mainly l and n
c
      do i = 1, nfreq
         read (11, *) ltemp, ntemp
         l(i) = ltemp
         n(i) = ntemp
      enddo
      close(11)
c
c     assign the indices  (subs)
c
      do i = 1, nfreq
         subscripts(i) = 100 + 100*l(i)+n(i)
      enddo

     
      do i = nfreq + 1, nfreq + nonseis
         if (ell(i) .eq. 5) then 
            subscripts(i) = 4
         elseif (ell(i) .eq. 6) then 
            subscripts(i) = 24
         elseif (ell(i) .eq. 7) then 
            subscripts(i) = 21
         elseif (ell(i) .eq. 8) then
            subscripts(i) = 27
         elseif (ell(i) .eq. 9) then
            subscripts(i) = 3
         endif
      enddo

c
c     hardwire a few things for code that are used in mrqdata_oc.f version only
c     we do not need them here but need to save them to /common
c
               
      rng(1) = 1400.d0
      rng(2) = 2400.d0
      rng(3) = 3800.d0

      end




      subroutine wrpar(npar, par)
      integer npar, ia(36), i ,onpar, x1
      double precision par(36), xx(100), x2
      character*55 a 
 
 1    format (i5)
      open (12, FILE = 'PARAMETERS', STATUS = 'OLD')
 3    format (d16.10)
 5    format (a55)

      read (12, 1) onpar
      do i = 1, onpar
         read (12, 3) x2
         xx(i) = x2
      enddo
      read (12, 5) a
      do i = 1, onpar 
         read (12, 1) x1
         ia(i) = x1
c         print *,x1
      enddo


      close(12)

      a = 'c ...... '
      open (14, FILE = 'PARAMETERS', STATUS = 'OLD')
      write (14, 1) npar
      do i = 1, npar
         write (14, 3) par(i)
      enddo
      write (14, 5) a
      do i = 1, npar 
         write (14,1) ia(i)
      enddo

      close(14)
      end


      subroutine wria(npar, iai)
      integer npar, ia(36), i ,onpar, x1, iai(36)
      double precision par(36), xx(100), x2
      character*55 a 
 
 1    format (i5)
      open (12, FILE = 'PARAMETERS', STATUS = 'OLD')
 3    format (d16.10)
 5    format (a55)

      read (12, 1) onpar
      do i = 1, onpar
         read (12, 3) x2
         xx(i) = x2
      enddo
      read (12, 5) a
      do i = 1, onpar 
         read (12, 1) x1
         ia(i) = x1
c         print *,x1
      enddo


      close(12)

      a = 'c ...... '
      open (14, FILE = 'PARAMETERS', STATUS = 'OLD')

      write (14, 1) npar
      do i = 1, npar
         write (14, 3) xx(i)
      enddo
      write (14, 5) a
      do i = 1, npar 
         write (14,1) iai(i)
      enddo

      close(14)
      end
 
 

      

