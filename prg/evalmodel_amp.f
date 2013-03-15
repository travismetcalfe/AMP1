      subroutine grmm_slave
c ---------------------------------------------
c evalmodel slave program
c ---------------------------------------------
      implicit none

      include 'mpif.h'

      integer myid, ierr, status(MPI_STATUS_SIZE)
      integer master, msgtype, trial, npar, nobs, nsel
      double precision data(32), yy(10000)

c ---------------------------------------------
c identify this slave task
c ---------------------------------------------
      call mpi_comm_rank( MPI_COMM_WORLD, myid, ierr )

c ---------------------------------------------
c listen for a new job
c ---------------------------------------------
      master = 0
 25   msgtype = 1

c ---------------------------------------------
c receive data from master host
c ---------------------------------------------
      call mpi_recv( trial, 1, MPI_INTEGER, master,
     +               msgtype, MPI_COMM_WORLD, status, ierr )
      if (trial .EQ. -1) goto 99
      call mpi_recv( npar, 1, MPI_INTEGER, master,
     +               msgtype, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( data, npar, MPI_DOUBLE_PRECISION, master,
     +               msgtype, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( nobs, 1, MPI_INTEGER, master,
     +               msgtype, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( nsel, 1, MPI_INTEGER, master,
     +               msgtype, MPI_COMM_WORLD, status, ierr )

c ---------------------------------------------
c perform calculations with data
c ---------------------------------------------
      call evalmodel (npar, data, nobs, nsel, yy)

c ---------------------------------------------
c send result to master host
c ---------------------------------------------      
      msgtype = 2 
      call mpi_send( myid, 1, MPI_INTEGER, master,
     +               msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( trial, 1, MPI_INTEGER, master,
     +               msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( nobs, 1, MPI_INTEGER, master,
     +               msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( yy, nobs, MPI_DOUBLE_PRECISION, 
     +       master, msgtype, MPI_COMM_WORLD, ierr )

c ---------------------------------------------
c go back for more work
c ---------------------------------------------
      goto 25

 99   return
      end

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c     TAKE INPUT PARAMETERS AND RUN THE MODEL, GET OUTPUT OBSERVABLES
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine evalmodel (npar, data, nobs, nsel, yy)



      implicit double precision (a-h,o-z)

      include 'engenr.cz.d.incl'
      double precision yy(nobs),fitness,userff,data(npar), age, R_Ro
      double precision Teff, L_Lo, xc, qc, gg, msol, logg, pi
      double precision num, den, meandensity, mbol, ff(1000), rng(3)
      double precision avgs(12), f0, xinitial
      integer nsel, npar, dpar, subs(nobs), mfit, nparam
      integer subscripts(10000), nonseis, nfreq,iwflag,isflag,idif
c     other declarations to make ff_amp work
      logical quiet

c
c  common with evolution sequence RESULTS!
c
      common /csum_param/ icsum, nstep, csum_st(icsum_max, nstep_max)
      common /csum_indiv/ icsum_ind, nstep_ind,
     *  csum_ind(icsum_max, nstep_max) 
c      common /xmodage/ age
c
c  common with pulsation RESULTS! (degree, order, frequency, inertia)
c
      common/cobs_param/ icobs_st, nobs_st, obs_st(10,1)
c
c      common/olcdint/ rng, subs, f0, nfreq
      common/olcdint/ nfreq,nonseis,subscripts,iwflag,isflag,idif
      common/olcddbl/ f0, rng

c     other declarations to make ff_amp work
      common /verbosity/ quiet
      quiet = .FALSE.
c     If we don't want to generate output:
      quiet = .TRUE.   

      nparam = 5
      fitness = userff(nparam,data,0)
      gg = 6.673e-5
      msol = 1.98892e30
      logg = gg*csum_ind(1,nstep_ind)*msol/(csum_ind(3,nstep_ind)**2.)
      pi = 2.*acos(0.0)
      num = csum_ind(1,nstep_ind)*msol
      den = (4.*pi*(csum_ind(3,nstep_ind)/100)**3.)/3.
      meandensity = num/den

      call seismo(ff)
c
c     get average separations
c
      call avgsep(ff, rng, avgs) 
c
c
c     assign global parameters (notes p 25, same positions)
c
c     csum(1) = mass (solar masses)
c     csum(2) = age (years)
c     csum(3) = surface radius (cm) ... but i divide by 6.96d+10 for output
c     csum(4) = effective temperature (K)
c     csum(5) = surface luminosity (erg/sec) .. i divide by 3.846d+33
c     csum(6) = depth of outermost convection zone (units of surface radius)
c     csum(7) = p_c central pressure
c     csum(8) = T_c central temperature
c     csum(9) = X_c central hydrogen
c     csum(10) = X_3,c
c     csum(11) = rho_c central density
c     csum(12) = epsilon_c 
c     csum(13) = kappa_c
c     csum(14) = nabla_rad,c
c     csum(15) = nabla_ad,c
c     csum(16) = m_c/M (mass in convectively mixed core
c     csum(17) = r_c/R (fraction radius of convectively mixed core)
c     csum(18) = X(14 N)  abundance
c     csum(19) = X(16 O) abundance
c     
c     yy(20) = log_10 (z/x)
c     yy(21) = log_10 (z/x) - log_10(z/x)_solar
c     yy(22) = distance (parsec)
c     yy(23) = parallax (mas)
c     yy(24) = log g  ;;;; g = GM/r**2: units alog10(cm**3 kg**-1 s**-2)
c     yy(25) = density (mean) ;;; M/(4/3 pi r^3) (my calc is kg/m^3)
c     yy(26) = bolometric magnitude
c     yy(27) = luminosity (1 J/s = 1e-7 erg/s = 1 W) ;solar = 3.24e26J/s
c     yy(28) = xinitial
c

      do i = 1, 19
         yy(i) = csum_ind (i, nstep_ind)
      enddo
      yy(3) = yy(3)/6.96d+10
      yy(5) = yy(5)/3.846d+33
      xinitial = 1.d0 - data(3) - data(2)
      yy(20) = log10(data(2)/xinitial)
      yy(21) = log10(data(2)/xinitial) - alog10(0.02/0.71)
      yy(22) = data(npar+1)
      yy(23) = 1.0/yy(22) *1000.0
      yy(24) = log10(logg)
      yy(25) = meandensity 
      mbol=42.36-5*log10(yy(3))-10*log10(yy(4)) 
      yy(26) = mbol
      yy(27) = csum_ind(5, nstep_ind)*1e-7

c
c     these are the initial parameters
c
c     yy(30) = data(1) 
c     yy(31) = data(2)
c     yy(32) = data(3) 
c     yy(33) = data(4)
c     yy(34) = data(5)
c     yy(35) = xinitial
c

      yy(30) = data(1) 
      yy(31) = data(2)
      yy(32) = data(3) 
      yy(33) = data(4)
      yy(34) = data(5)
      yy(35) = xinitial

c
c     assign average separations 
c
c     yy(50) = Delta_0 (lower range of frequencies)
c     yy(51) = Delta_0 (higher range of frequencies)
c     yy(52) = Delta_1 (lower range of frequencies)
c     yy(53) = Delta_1 (higher range of frequencies)
c     yy(54) = Delta_2 (lower range of frequencies)
c     yy(55) = Delta_2 (higher range of frequencies)
c     yy(56) = Delta_3 (lower range of frequencies)
c     yy(57) = Delta_4 (higher range of frequencies)
c     yy(58) = delta_02 (lower range of frequencies)
c     yy(59) = delta_02 (higher range of frequencies)
c     yy(60) = delta_13 (lower range of frequencies)
c     yy(61) = delta_13 (higher range of frequencies)
c
c
c
      do i = 1, 12
         yy(i+49) =  avgs(i) 
      enddo

c
c     assign frquencies 
c
c     yy(100:199) = l = 0
c     yy(200:299) = l = 1
c     yy(300:399) = l = 2
c     yy(400:499) = l = 3
c     yy(500:599) = Delta l = 0
c     yy(600:699) = Delta l = 1
c     yy(700:799) = Delta l = 2
c     yy(800:899) = Delta l = 3
c     yy(900:999) = delta l = 02
c     yy(1000:1099) = delta l = 13

      do i = 1, 999
         yy(i+100) = ff(i)
      enddo
      
c     call 'surface' to modify (or not) the frequencies due to surface effects
      call surface (nsel, yy) 

      end
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c     SAVE THE FREQUENCIES INTO THE ARRAY FF
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine seismo(ff)

      implicit double precision (a-h,o-z)

      include 'engenr.cz.d.incl'
      double precision ff(1000)
      real k
c
c  common with pulsation RESULTS! (degree, order, frequency, inertia)
c
      common/cobs_param/ icobs_st, nobs_st, obs_st(10,1)
c
c      save done_setups, /xmodage/

c     assign individual frequencies
      do i = 1, nobs_st
         if (obs_st(1,i).eq.0) then
            ff(obs_st(2,i)) = obs_st(3,i)
         endif
         if (obs_st(1,i).eq.1) then
            ff(obs_st(2,i)+100) = obs_st(3,i)
         endif
         if (obs_st(1,i).eq.2) then
            ff(obs_st(2,i)+200) = obs_st(3,i)
         endif
         if (obs_st(1,i).eq.3) then
            ff(obs_st(2,i)+300) = obs_st(3,i)
         endif
      enddo

c
c     then frequqnecy separations ....
c     NOTE: ff(403) = v_0,n - v_0,n-1 (ff(3)-f(2))
c
        do i = 1, 48
            ff(i+400) = ff(i)-ff(i-1)
            ff(i+500) = ff(i+100)-ff(i-1+100)
            ff(i+600) = ff(i+200)-ff(i-1+200)
            ff(i+700) = ff(i+300)-ff(i-1+300)
            ff(i+800) = ff(i) - ff(i+200-1)
            ff(i+900) = ff(i+100) - ff(i+300-1)

      enddo


      end
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c     GET THE AVERAGE SEPARATION VALUES BASED ON INPUT RANGES OF FREQ.
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c     This will never be used in the amp version
c
c
      subroutine avgsep(ff, rng, avgs)

      implicit none
      double precision  ff(1000), rng(3), avgs(12), df, sf
      integer k, i 

c
c 1.    l = 0 average separations in 2 ranges
c
      df = 0.0
      k = 0
      do i = 2, 49
         if ((ff(i).le.rng(2)).AND.(ff(i).ge.rng(1))) then
            if ((ff(i-1).le.rng(2)).AND.(ff(i-1).ge.rng(1))) then        
               df = df + ff(i+400)
               k = k + 1
            endif
         endif
      enddo

      avgs(1) = df/float(k)

      df = 0.0
      k = 0
      do i = 2, 49
         if ((ff(i).le.rng(3)).AND.(ff(i).ge.rng(2))) then
            if ((ff(i-1).le.rng(3)).AND.(ff(i-1).ge.rng(2))) then        
               df = df + ff(i+400)
               k = k + 1
            endif
         endif
      enddo
      avgs(2) = df/float(k)
      
c
c 2.    small freq separations for l = 0
c

      
      df = 0.0
      k = 0
      do i = 2, 49
         if ((ff(i).le.rng(2)).AND.(ff(i).ge.rng(1))) then
         if ((ff(i-1+200).le.rng(2)).AND.(ff(i-1+200).ge.rng(1))) then        
               df = df + ff(i+800)
               k = k + 1
            endif
         endif
      enddo
      avgs(9) = df/float(k)

      df = 0.0
      k = 0
      do i = 2, 49
         if ((ff(i).le.rng(3)).AND.(ff(i).ge.rng(2))) then
         if ((ff(i-1+200).le.rng(3)).AND.(ff(i-1+200).ge.rng(2))) then        
               df = df + ff(i+800)
               k = k + 1
            endif
         endif
      enddo
      avgs(10) = df/float(k)


c
c 3.    l = 1 average separations in 2 ranges
c
      df = 0.0
      k = 0
      do i = 2, 49
         if ((ff(i+100).le.rng(2)).AND.(ff(i+100).ge.rng(1))) then
         if ((ff(i-1+100).le.rng(2)).AND.(ff(i-1+100).ge.rng(1))) then        
               df = df + ff(i+500)
               k = k + 1
            endif
         endif
      enddo
 
      avgs(3) = df/float(k)

      df = 0.0
      k = 0
      do i = 2, 49
         if ((ff(i+100).le.rng(3)).AND.(ff(i+100).ge.rng(2))) then
         if ((ff(i-1+100).le.rng(3)).AND.(ff(i-1+100).ge.rng(2))) then        
               df = df + ff(i+500)
               k = k + 1
            endif
         endif
      enddo
      avgs(4) = df/float(k)
      
c
c 4.    small freq separations for l = 1
c

      
      df = 0.0
      k = 0
      do i = 2, 49
         if ((ff(i+100).le.rng(2)).AND.(ff(i+100).ge.rng(1))) then
         if ((ff(i-1+300).le.rng(2)).AND.(ff(i-1+300).ge.rng(1))) then        
               df = df + ff(i+900)
               k = k + 1
            endif
         endif
      enddo
      avgs(11)= df/float(k)
       df = 0.0
      k = 0
      do i = 2, 49
         if ((ff(i+100).le.rng(3)).AND.(ff(i+100).ge.rng(2))) then
         if ((ff(i-1+300).le.rng(3)).AND.(ff(i-1+300).ge.rng(2))) then        
               df = df + ff(i+900)
               k = k + 1
            endif
         endif
      enddo
      avgs(12) = df/float(k)


c
c 5.    l = 2 average separations in 2 ranges
c
      df = 0.0
      k = 0
      do i = 2, 49
         if ((ff(i+200).le.rng(2)).AND.(ff(i+200).ge.rng(1))) then
         if ((ff(i-1+200).le.rng(2)).AND.(ff(i-1+200).ge.rng(1))) then        
               df = df + ff(i+600)
               k = k + 1
            endif
         endif
      enddo

      avgs(5) = df/float(k)

      df = 0.0
      k = 0
      do i = 2, 49
         if ((ff(i+200).le.rng(3)).AND.(ff(i+200).ge.rng(2))) then
         if ((ff(i-1+200).le.rng(3)).AND.(ff(i-1+200).ge.rng(2))) then        
               df = df + ff(i+600)
               k = k + 1
            endif
         endif
      enddo
      avgs(6) = df/float(k)

c
c 6.    l = 3 average separations in 2 ranges
c
      df = 0.0
      k = 0
      do i = 2, 49
         if ((ff(i+300).le.rng(2)).AND.(ff(i+300).ge.rng(1))) then
         if ((ff(i-1+300).le.rng(2)).AND.(ff(i-1+300).ge.rng(1))) then        
               df = df + ff(i+700)
               k = k + 1
            endif
         endif
      enddo

      avgs(7) = df/float(k)

      df = 0.0
      k = 0
      do i = 2, 49
         if ((ff(i+300).le.rng(3)).AND.(ff(i+300).ge.rng(2))) then
         if ((ff(i-1+300).le.rng(3)).AND.(ff(i-1+300).ge.rng(2))) then        
               df = df + ff(i+700)
               k = k + 1
             endif
         endif
      enddo
      avgs(8) = df/float(k)

      end


c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c     IF DOSURF NE 0 THEN WE NEED TO CALCUALTE SOME SURFACE EFFECTS...
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine surface(nobs, ymod)

c
c     Calculate the surface effects and check the flags iwflag and isflag
c     to write the correct systematic error to the common block, which 
c     will be used in the chi^2 calculation
c

      implicit none
      integer nobs
      double precision ymod(nobs)
c
c     data & sub parameters
c
      integer npar, ia(36), nsel, ell(10000), subs(10000)
      integer ndata
      double precision dat(36), yi(10000), sig(10000)
c
c     parameters within routine
c
      integer i, num_obs, num_fm, num_rnu 
      double precision dnu_obs, dnu_calc, sum_obs, sum_fm, fm, sum_n
      double precision sum_top, sum_bot, n0, sum_freq
      double precision bexp, c1, c2, r, a0, sys, syserror(nobs)

c
c     parameters from common
c
      integer nfreq, nonseis,subscripts(10000),iwflag,isflag,idif
      double precision f0, rng(3), surfsigi(100), surfcorr(100)

      common/olcdint/ nfreq,nonseis,subscripts,iwflag,isflag,idif
      common/olcddbl/ f0, rng

c     surfsigi are the systematic errors to use for chi^2, surfcorr is just 
c     to save the surface correction if iwflag eq 0
      common/olcsurf/ surfsigi, surfcorr



c     GET DATA, (I need: NOBS, YI, SIG, ELL)
      call mrqdata (ia, npar, dat, ndata, yi, sig, ell)

      c1 = 1.2616d0
      c2 = 0.2616d0
      bexp = 4.823d0

c
c     get n0, to calculate dnu_obs
c
      sum_n = 0
      num_rnu = 0
      sum_freq = 0.
      do i = 1, ndata
         if (ell(i).eq.0) then 
           sum_n = sum_n + (i-100)
           num_rnu = num_rnu + 1
           sum_freq = sum_freq + yi(i)
        endif
      enddo
      n0 = sum_n / float(num_rnu)
      f0 = sum_freq / float(num_rnu)

      num_obs = 0
      sum_obs=0.
      sum_top=0.   
      sum_bot = 0.
      do i = 1, ndata
         if (ell(i).eq.0) then 
            sum_obs = sum_obs + (yi(i)/f0)**bexp
            num_obs = num_obs  + 1
            sum_top = sum_top + (yi(i)-f0)*((i-100)-n0)
            sum_bot = sum_bot + ((i-100)-n0)*((i-100)-n0)
         endif
      enddo
      dnu_obs = sum_top / sum_bot
c
c     calculate the n0, fm for dnu_calc
c
      sum_n = 0
      sum_fm = 0.0d0
      num_rnu = 0
      do i = 1, ndata
         if ((subscripts(i).gt.100).AND.(subscripts(i).le.199)) then
            sum_n = sum_n + (i-100)
            sum_fm = sum_fm + ymod(subscripts(i))
            num_rnu = num_rnu + 1
         endif
      enddo     
      n0 = sum_n / float(num_rnu)
      fm = sum_fm / float(num_rnu)

      num_rnu = 0
      sum_top = 0.
      sum_bot = 0.
      do i = 1, ndata
         if ((subscripts(i).gt.100).AND.(subscripts(i).le.199)) then
            sum_top = sum_top + (ymod(subscripts(i)) - fm)*((i-100)-n0)
            sum_bot = sum_bot + ((i-100)-n0)*((i-100)-n0)
            num_rnu = num_rnu + 1
         endif
      enddo
      dnu_calc = sum_top/sum_bot
c
c     calculate the r and a0 term
c
       r = 1./(c1 * (fm/f0) - c2 *(Dnu_calc/Dnu_obs))
       a0 = ((f0 - fm) * float(num_obs))/ sum_obs


c
c     calculate systematic offset to model frequencies
c     and the new weights for the errors
c
        do i = 1, ndata
         if ((subscripts(i).gt.100).AND.(subscripts(i).le.499)) then
            sys = a0*(ymod(subscripts(i))/f0)**bexp
            ymod(subscripts(i)) = ymod(subscripts(i)) + sys
            surfsigi(i)  =  sys
            surfcorr(i) = sys
            if (iwflag.eq.0) surfcorr(i) = 0.d0
         endif
      enddo
      

c
c     Now check for if surface effects are turned on or not (isflag), and 
c     whether to include them in the chi^ 2 calculation (iwflag)
c     And save to common
c
      if (isflag.eq.0) then 
         do i=1,ndata
            surfsigi(i) = 0.0d0
            surfcorr(i) = 0.0d0
         enddo
      endif


      end
   
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c     READ IN THE SURFACE ERRORS...
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine surfacesigi(nobs, sigi)

c
c     this routine has changed.  it used to read in the surface
c     correction from 
c     a file, but now it just copies it from the common and
c     puts it in the vector sigi
c
 
      implicit none
      integer i, nobs
      double precision  surfsigi(100), surfcorr(100), sigi(nobs)

      common/olcsurf/ surfsigi, surfcorr


c     assign sigi the surface correction
      do i = 1, 100
         sigi(i) = surfsigi(i)
      enddo


      end
