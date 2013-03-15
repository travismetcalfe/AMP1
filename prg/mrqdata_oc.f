      subroutine mrqdata (ia, npar, xx, nobs, yy, sig, ell)


      implicit none
      integer ma, ia(36), i, sub(10000), nobs, x1, npar, ell(10000), x4
      double precision  yy(10000),sig(10000), x3, xx(36), x2
      integer col1int, nip
      character a*55, col1str*1, par_string*80
      integer nsel,  subs(10000)
      double precision  rng(3)

      real ntemp, n(50)
      integer ltemp, l(50)
      integer nonseis, nfreq, subscripts(10000), iwflag, isflag, idif
      character a1*10, a2*10
      double precision f0

      common/olcdint/ nfreq,nonseis,subscripts,iwflag,isflag,idif
      common/olcddbl/ rng

c
c  read obs.dat file to get data, errors and ell
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

c
c     GET INITIAL PARAMETERS & IA
c
      open (12, FILE = 'PARAMETERS', STATUS = 'OLD')
 13    format (d16.10)
 5    format (a55)

      read (12, 1) npar
      do i = 1, npar
         read (12, 13) x2
         xx(i) = x2
      enddo
      read (12, 5) a
      do i = 1, npar 
         read (12, 1) x1
         ia(i) = x1
      enddo

      close(12)

c
c     OR READ THEM FROM THE COMMAND LINE
c
      
      nip = IARGC()
      do i=1,nip
         CALL GETARG(i,par_string)
         READ(par_string,*) xx(i)
      enddo

c
c     read the param.log file to get "details"
c
c

      close(11)
 
      open ( 11, FILE = 'param.log', STATUS = 'OLD' )
 3    format (a55)
   
c
c     read in the first set of lines until we get to the l and n
c
      do i = 1, 7
         read (11, 3) a
      enddo

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

 
c      get DETAILS about whether to dosurface correction and 
c     range of frequencies to use
c

      open (11, FILE = 'DETAILS', STATUS = 'OLD')
 1    format (i5)
c
c     read number of observables and then the subscripts
c
      read (11, 3) a
      read (11, 1) nsel
      read (11, 3) a
      do i = 1, nsel
         read (11, 1) x1
         subs(i) = x1
      enddo
      read (11, 3) a
c
c     read the frequency ranges .... to calculate average sep
c
 2    format (f9.5)
      do i = 1, 3
         read (11, 2) x2
         rng(i) = x2
      enddo
      read (11, 3) a

c
c     read info on surface effects  0 = none, 1 = yes
c
c 4    format (i1)
c      read (11, 4) dosurf
c      read (11, 3) a
c
c     read info on f0 term for surface effects
c
c 5    format (f16.9)
c      read (11, 5) f0
      close(11)


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
 
 


c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c     READ THE FILE "DETAILS" TO GET VARIOUS INPUTS TO RUN THE CODE...
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c
c
c This routine should be obselete now, because I pass this information 
c in a common block
c

      subroutine rdsubs(nsel, subs, rng, f0)

c
c     this routine will read in a file with information about which 
c     observables to use and which range of frequencies to use to 
c     calculate the average separation as well as a flag to tell it 
c     to turn on surface effects.
c
      implicit none
      integer i, subs(10000), nsel, x1
      double precision  x2, rng(3), f0
      character*70 a

       close(11)
     
      open (11, FILE = 'DETAILS', STATUS = 'OLD')
 3    format (a70)
 1    format (i5)
c
c     read number of observables and then the subscripts
c
      read (11, 3) a
      read (11, 1) nsel
      read (11, 3) a
      do i = 1, nsel
         read (11, 1) x1
         subs(i) = x1
      enddo
      read (11, 3) a
c
c     read the frequency ranges .... to calculate average sep
c
 2    format (f9.5)
      do i = 1, 3
         read (11, 2) x2
         rng(i) = x2
      enddo
      read (11, 3) a

c
c     read info on surface effects  0 = none, 1 = yes
c
c 4    format (i1)
c      read (11, 4) dosurf
c      read (11, 3) a
c
c     read info on f0 term for surface effects
c
c 5    format (f16.9)
c      read (11, 5) f0
c      close(11)
 
      end

