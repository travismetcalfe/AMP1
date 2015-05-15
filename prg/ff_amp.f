      function userff ( npar, data, myid )
c
c        npar = number of parameters
c  data(npar) = scaled parameter values [0.0 -> 1.0]
c        myid = parallel slave number [1 -> Nproc]
c
      implicit double precision (a-h,o-z)

      integer npar,myid,ell(100),obs,calc,high,low,match(100),tstep
      integer n(100),peak,gap,nn,nmin,nmax,num_rat,num_r02
      double precision data(36),freq(100),spacing(400),err(100),target
      double precision step,log_z,avg_sys(5),Tresid,Lresid,ominusc
      double precision Teff,L_Lo,T_obs,L_obs,Terr,Lerr,offset(400)
      double precision R_Ro,Rerr,Rresid,n0,sum_n,sum_top,sum_bot
      double precision mixed(100),nu,rnu,lowf,highf,lowI,highI
      double precision rinertia,inertia,G_obs,Gerr,logG,Gresid
      double precision M_obs,Merr,M_H,Mresid,Dnu_in,penalty,chisq(4)
      double precision nu_calc(100,5),nu_obs(100,5),nu_err(100,5)
      double precision Dnu0_calc(100),Dnu0_obs(100),Dnu0_err(100)
      double precision Dnu1_calc(100),Dnu1_obs(100),Dnu1_err(100)
      double precision d02_calc(100),d02_obs(100),d02_err(100)
      double precision d13_calc(100),d13_obs(100),d13_err(100)
      double precision d01_calc,d01_obs,d01_err,d10_calc,d10_obs,d10_err
      double precision r01_calc,r01_obs,r01_err,r10_calc,r10_obs,r10_err
      double precision r02_calc,r02_obs,r02_err,chisq_rat,chisq_r02
      double precision r13_calc,r13_obs,r13_err

      logical isnan, quiet
      real userff
      character fname*11, header_par*80, tfile*80, col1str*1
      logical done_setups
      data spacing / 400*0.0 /
      data header_par, icase, tfile, amtg1, amtg2, damtg,
     +     ztg1, ztg2, nztg, xhtg1, xhtg2, dxhtg
     +    / "param", 1, "emdl.rho.0100.Z2.01.s.0", 0.75, 1.75, 0.05,
     +      0.002, 0.05, 8, 0.70, 0.70, 0.0 /
      data bexp,c1,c2 / 4.823d0, 1.2616d0, 0.2616d0 /

      real h1(195),h2(195),hfrq(195),w(100)
      data h1   /0.21538462,0.21923077,0.22307692,0.22692308,0.23076923,
     +0.23461538,0.23846154,0.24230769,0.24615385,0.25000000,0.25384615,
     +0.25769231,0.26153846,0.26538462,0.26923077,0.27307692,0.27692308,
     +0.28076923,0.28461538,0.28846154,0.29230769,0.29615385,0.30000000,
     +0.30384615,0.30769231,0.31153846,0.31538462,0.31923077,0.32307692,
     +0.32692308,0.33076923,0.33461538,0.33846154,0.34230769,0.34615385,
     +0.35000000,0.35384615,0.35769231,0.36153846,0.36538462,0.36923077,
     +0.37307692,0.37692308,0.38076923,0.38461538,0.38846154,0.39230769,
     +0.39615385,0.40000000,0.40384615,0.40769231,0.41153846,0.41538462,
     +0.41923077,0.42307692,0.42692308,0.43076923,0.43461538,0.43846154,
     +0.44230769,0.44615385,0.45000000,0.45384615,0.45769231,0.46153846,
     +0.46538462,0.46923077,0.47307692,0.47692308,0.48076923,0.48461538,
     +0.48846154,0.49230769,0.49615385,0.50000000,0.50384615,0.50769231,
     +0.51153846,0.51538462,0.51923077,0.52307692,0.52692308,0.53076923,
     +0.53461538,0.53846154,0.54230769,0.54615385,0.55000000,0.55384615,
     +0.55769231,0.56153846,0.56538462,0.56923077,0.57307692,0.57692308,
     +0.58076923,0.58461538,0.58846154,0.59230769,0.59615385,0.60000000,
     +0.60384615,0.60769231,0.61153846,0.61538462,0.61923077,0.62307692,
     +0.62692308,0.63076923,0.63461538,0.63846154,0.64230769,0.64615385,
     +0.65000000,0.65384615,0.65769231,0.66153846,0.66538462,0.66923077,
     +0.67307692,0.67692308,0.68076923,0.68461538,0.68846154,0.69230769,
     +0.69615385,0.70000000,0.70384615,0.70769231,0.71153846,0.71538462,
     +0.71923077,0.72307692,0.72692308,0.73076923,0.73461538,0.73846154,
     +0.74230769,0.74615385,0.75000000,0.75384615,0.75769231,0.76153846,
     +0.76538462,0.76923077,0.77307692,0.77692308,0.78076923,0.78461538,
     +0.78846154,0.79230769,0.79615385,0.80000000,0.80384615,0.80769231,
     +0.81153846,0.81538462,0.81923077,0.82307692,0.82692308,0.83076923,
     +0.83461538,0.83846154,0.84230769,0.84615385,0.85000000,0.85384615,
     +0.85769231,0.86153846,0.86538462,0.86923077,0.87307692,0.87692308,
     +0.88076923,0.88461538,0.88846154,0.89230769,0.89615385,0.90000000,
     +0.90384615,0.90769231,0.91153846,0.91538462,0.91923077,0.92307692,
     +0.92692308,0.93076923,0.93461538,0.93846154,0.94230769,0.94615385,
     +0.95000000,0.95384615,0.95769231,0.96153846/
      data h2    / 0.00000000,-0.00455772,-0.00903524,-0.01329270,
     +-0.01718160,-0.02054358,-0.02321156,-0.02501226,-0.02576128,
     +-0.02526680,-0.02334288,-0.02002764,-0.01555160,-0.01016232,
     +-0.00412720, 0.00226774, 0.00872352, 0.01492120, 0.02052464,
     + 0.02518050, 0.02852280, 0.03040114, 0.03107208, 0.03085108,
     + 0.03007360, 0.02908872, 0.02826868, 0.02800254, 0.02869776,
     + 0.03077850, 0.03468896, 0.04077690, 0.04890160, 0.05878094,
     + 0.07011540, 0.08259342, 0.09588608, 0.10964700, 0.12351788,
     + 0.13711920, 0.15006528, 0.16194732, 0.17240160, 0.18108090,
     + 0.18762200, 0.19165558, 0.19279632, 0.19065094, 0.18481008,
     + 0.17485440, 0.16035044, 0.14086122, 0.11614536, 0.08625170,
     + 0.05126880, 0.01128648,-0.03359552,-0.08326970,-0.13762308,
     +-0.19653040,-0.25986088,-0.32748066,-0.39919636,-0.47465768,
     +-0.55346880,-0.63522096,-0.71948280,-0.80581728,-0.89376224,
     +-0.98284750, -1.0725826, -1.1624640, -1.2520333, -1.3413417,
     + -1.4307163, -1.5204981, -1.6110391, -1.7027030, -1.7958640,
     + -1.8909072, -1.9882330, -2.0882444, -2.1913627, -2.2978479,
     + -2.4077872, -2.5212605, -2.6383430, -2.7591106, -2.8836432,
     + -3.0120212, -3.1443261, -3.2806372, -3.4210378, -3.5655402,
     + -3.7139550, -3.8660590, -4.0216160, -4.1803884, -4.3421316,
     + -4.5065847, -4.6734886, -4.8425771, -5.0135644, -5.1862016,
     + -5.3602912, -5.5356662, -5.7121556, -5.8896823, -6.0680033,
     + -6.2469693, -6.4264277, -6.6061893, -6.7860946, -6.9659130,
     + -7.1456134, -7.3250962, -7.5042946, -7.6831757, -7.8617063,
     + -8.0398185, -8.2174435, -8.3945471, -8.5710596, -8.7469460,
     + -8.9221356, -9.0965930, -9.2702464, -9.4430599, -9.6149237,
     + -9.7858007, -9.9556165, -10.124371, -10.291951, -10.458281,
     + -10.623284, -10.786805, -10.948573, -11.108350, -11.265855,
     + -11.420842, -11.573059, -11.722253, -11.868124, -12.010450,
     + -12.148924, -12.283596, -12.414682, -12.542359, -12.666976,
     + -12.788765, -12.907923, -13.024775, -13.139531, -13.252569,
     + -13.364068, -13.474338, -13.583438, -13.691516, -13.798639,
     + -13.905000, -14.010628, -14.115681, -14.220319, -14.324575,
     + -14.428660, -14.532655, -14.636420, -14.739769, -14.842558,
     + -14.944684, -15.045864, -15.145989, -15.244860, -15.342363,
     + -15.438203, -15.532260, -15.624230, -15.713711, -15.800293,
     + -15.883702, -15.963422, -16.039122, -16.110415, -16.176815,
     + -16.237925, -16.293388, -16.343088, -16.387536, -16.427158,
     + -16.462584, -16.494256, -16.522724, -16.548599, -16.572399,
     + -16.594705/

      include 'engenr.cz.d.incl'
      parameter(iaa_adi=11)
      character trailer_par*80, id*1, settrailer*80,
     *  strcompr*80, file*80, filess*80, par_trial*80, set_trial*80,
     *  in_evol*80, in_rdist*80, in_adi*80, paramfile*80
      dimension xr_adi(nnmax), aar_adi(iaa_adi,nnmax)
      common/cvr_param/ par_am, par_z, par_agefin, par_rsfin, 
     *  par_alsfin, par_zxsfin, par_xmdtrl,
     *  par_xxh, par_fdgopl, par_tlopfg, par_dtlopf, par_tlopf1,
     *  par_fcno, par_agehe3, par_xrz12, par_xrz13, par_xrz14,
     *  par_xrz16, par_tprfct, par_alfa, par_etac, par_phc, par_alpove,
     *  par_alpovc, par_velrot
      common/cvi_param/ ipar_nt, ipar_icsove, ipar_icsovc, ipar_isprot,
     *  ipar_icmout, ipar_istosc, ipar_isetos
      common/cadr_param/
     *  para_el, para_els1, para_dels, para_dfsig1, para_dfsig2,
     *  para_sig1, para_sig2, para_dfsig, para_eltrw1, para_eltrw2,
     *  para_sgtrw1, para_sgtrw2
      common/cadi_param/
     *  ipara_nsel, ipara_nsig1, ipara_nsig2, ipara_itrsig, ipara_nsig,
     *  ipara_istsig, ipara_inomd1, ipara_iscan, ipara_irotkr
      common/trial_param/ par_trial
      common/trl_param/ trailer_par
      common/cofile/ nfiles, idsfil(99), file(99), iopen(99), filess(99)
      common/nmbmsh/ nn_evol, nn_adi, ivar_adi
      common/xnwvar/ x_adi(1)
      common/anwvar/ data_adi(8), aa_adi(istrmx,1)
c
c  dataset designators
c
      common/cdsdsg/ idstrl,idshvz,idszab,idsevl,idssum,idscen,idsnov,
     *  idsoov,idsgng,idsgsm,idstm1,idstm2,iducen,iddgm1
c
      common/cversn/ ivrbcs,ivrrhs,ivreos,ivropc,ivreng
c
c  commons defining input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      common/cstdio_in/ istdin_in, istdpr_in
      common/cevlio/ isetos, iastr
c
c  common with evolution sequence RESULTS!
c
      common /csum_param/ icsum, nstep, csum_st(icsum_max, nstep_max)
      common /csum_indiv/ icsum_ind, nstep_ind,
     *  csum_ind(icsum_max, nstep_max) 
      common /verbosity/ quiet
      common /xmodage/ age,R_Ro,L_Lo,Teff,M_H,chisq
c
c  common with pulsation RESULTS! (degree, order, frequency, inertia)
c
      common/cobs_param/ icobs_st, nobs_st, obs_st(10,1)
c
      save done_setups, /xmodage/
c
c  set defaults for debugging flags
c
      iwflag=0
      isflag=1
      idif=1
      Dnu_in=0.0
      iovs=0
      a_ovs=0.0
      if (myid .gt. 0) then
c
c  read obs.dat file
c
         nonseis=0
         open(55,file='obs.dat',status='old')
         read(55,*) num_obs,iwflag,isflag,idif,Dnu_in,iovs,a_ovs
         do obs=1,num_obs+5
            read(55,*,end=50) col1str,freq(obs),err(obs)
            if (col1str .eq. "T") then
               ell(obs)=5
               T_obs=freq(obs)
               Terr=err(obs)
               nonseis=nonseis+1
            elseif (col1str .eq. "G") then
               ell(obs)=6
               G_obs=freq(obs)
               Gerr=err(obs)
               nonseis=nonseis+1
            elseif (col1str .eq. "M") then
               ell(obs)=7
               M_obs=freq(obs)
               Merr=err(obs)
               nonseis=nonseis+1
            elseif (col1str .eq. "L") then
               ell(obs)=8
               L_obs=freq(obs)
               Lerr=err(obs)
               nonseis=nonseis+1
            elseif (col1str .eq. "R") then
               ell(obs)=9
               R_obs=freq(obs)
               Rerr=err(obs)
               nonseis=nonseis+1
            else
               read(col1str,*) col1int
               ell(obs)=col1int
            endif
         enddo
 50      close(55)
c
c  calcalate mean observed radial frequency
c
         num_rnu=0
         sum_freq=0.
         do obs=1,num_obs
            if (ell(obs) .eq. 0) then
               sum_freq = sum_freq + freq(obs)
               num_rnu = num_rnu + 1
            endif
         enddo
         f0 = sum_freq / float(num_rnu)
c
c  assign relative n values
c
         n(1) = 11
         sum_n = n(1)
         peak = INT(num_rnu/2+0.5)
         Dnu_approx = freq(peak+1) - freq(peak)
         do obs=2,num_rnu
            n(obs) = n(obs-1) + 1
            Dnu = freq(obs) - freq(obs-1)
            if (Dnu .GT. 1.5*Dnu_approx) then
               gap = INT(Dnu/Dnu_approx+0.5)
               n(obs) = n(obs-1) + gap
            endif
            sum_n = sum_n + n(obs)
         enddo
         n0 = sum_n / float(num_rnu)
c
c  least squares fit for observed large frequency spacing
c
         sum_obs=0.
         sum_top=0.
         sum_bot=0.         
         do obs=1,num_rnu
            sum_obs = sum_obs + (freq(obs)/f0)**bexp
            sum_top = sum_top + (freq(obs)-f0)*(n(obs)-n0)
            sum_bot = sum_bot + (n(obs)-n0)*(n(obs)-n0)
         enddo
         Dnu_obs = sum_top / sum_bot
         if (Dnu_in .gt. 0) Dnu_obs=Dnu_in
         Dnu_fit = Dnu_obs + 1.
         target = 0.0
      endif
c
c  scale the input parameters
c
      if (npar .GE. 1) then
	 par_am = data(1)
         if (myid .gt. 0) par_am = 1.00*data(1)+0.75         
         am_new=par_am
      else
         par_am = 1.00
         am_new=par_am
      endif
      if (npar .GE. 2) then
         par_z = data(2)
         log_z = 1.4*data(2)-2.7
         if (myid .gt. 0) par_z = 10.**(log_z)
      else
         par_z = 0.02
      endif
      if (npar .GE. 3) then
         par_y = data(3)
         if (myid .gt. 0) par_y = 0.10*data(3)+0.22
         par_xxh = 1.-par_y-par_z
      else
         par_xxh = 0.71
      endif
      if (npar .GE. 4) then
         par_alfa = data(4)
         if (myid .gt. 0) par_alfa = 2.0*data(4)+1.0
      else
         par_alfa = 1.9
      endif
      if (npar .GE. 5) then
         age_new = data(5)
         ipar_nt = 300
      else
         age_new = 0.d0
         ipar_nt = 300
      endif
c
c  check for parameters out of allowed range
c
      if (myid .eq. 0) then
         if ((par_am.lt.0.75).or.(par_am.gt.1.75)) 
     +      stop "M out of allowed range"
         if ((par_z.lt.0.0019).or.(par_z.gt.0.051)) 
     +      stop "Z out of allowed range"
         if ((par_y.lt.0.22).or.(par_y.gt.0.32)) 
     +      stop "Y out of allowed range"
         if ((par_alfa.lt.1.0).or.(par_alfa.gt.3.0)) 
     +      stop "alpha out of allowed range"
      endif
c
c  set number of time steps based on mass
c
      if (myid .gt. 0) then
         if (par_am.lt.1.09) 
     +        ipar_nt=160
         if (par_am.ge.1.09.and.par_am.le.1.2) 
     +        ipar_nt=320
         if (par_am.gt.1.2.and.par_am.lt.1.6) 
     +        ipar_nt=int(500.-150.*par_am)
         if (par_am.ge.1.6) 
     +        ipar_nt=260
      endif
c
c  direct verbose output to run.log
c
      istdou = 26
      istdpr = 26
      istder = 26
c      open(istdou,file='run.log',status='unknown')
      open(istdou,file='/dev/null',status='unknown')
c
c  hardwire input filenames
c
 1    format(".evolrin.",i3.3)
      write(in_evol,1) myid
 2    format(".redistrbin.",i3.3)
      write(in_rdist,2) myid
 3    format(".adiplsin.",i3.3)
      write(in_adi,3) myid
c
c  run setups on first call only
c
      if (done_setups) goto 10
c
c  flags for internal passing of emdl, non-verbose mode
c
      if (quiet) then
         iastr = -1
         istdpr_in = -1
      endif
c
c  create input files
c
      call write_evol(myid,idif,iovs,a_ovs)
      call write_rdist(myid)
      call write_adi(myid)
c
c  general setups of storage and initial parameters
c
      call setups_package(in_evol, in_rdist, in_adi, ierr_param)
      if(ierr_param.lt.0) then
        userff=0.0
        goto 99
      endif
c
c  test for existence of trial grid and read in catalogue
c
      call test_gridz(icase, itrial_err)
c
c  If not available, set grid of trial models
c
      if(itrial_err.lt.0) 
     *  call trial_gridz(amtg1, amtg2, damtg, ztg1, ztg2, nztg,
     *  xhtg1, xhtg2, dxhtg, tfile, icase, header_par)
c
c  initialize with zero age model                                              
c
      xmod_new=-1
      nt_max=-1
      age_max=0.1d0
      interpol=1
      newseq=1

      call set_modage(am_new, par_z, icase, xmod_new, age_new,
     *    interpol, nt_max, age_max, newseq, ierr_param)
c
      done_setups = .TRUE.
c
c  evolve the model track (entry point after first call)
c
   10 continue

   11 format(".paramlog.",i3.3)
      write(paramfile,11) myid
      open(55,file=paramfile,status='unknown')

   12 format("data: ",4(1x,F4.2),1x,E10.3)
      write(55,12) data(1),data(2),data(3),data(4),data(5)
      call flush(55)

      xmod_new=-1
      nt_max=-1
      age_max=age_new
      interpol=1
      newseq=1
      
      if (myid .gt. 0) then
         nt_max=ipar_nt
         age_max=0
         age_new=-1
         interpol=0
      endif

      call set_modage(am_new, par_z, icase, xmod_new, age_new, 
     *    interpol, nt_max, age_max, newseq, ierr_param)
      if(ierr_param.lt.0) then
        userff=0.0
        goto 99
      endif
c
c  SINGLE MODEL: calculate the mean model large frequency separation
c
      num_nu = 0
      sum_Dnu = 0.
      do mode=15,25
         l=obs_st(1,mode)
         l1=obs_st(1,mode+1)
         if ((l.eq.l1).and.(l.eq.0)) then
            Dnu = obs_st(3,mode+1) - obs_st(3,mode)
            if(Dnu.le.0) goto 13
            sum_Dnu = sum_Dnu + Dnu
            num_nu = num_nu + 1
         endif
 13      continue
      enddo
      Dnu_calc = sum_Dnu / num_nu
      age = csum_ind(2,nstep_ind)
      R_Ro = csum_ind(3,nstep_ind)/6.96d+10
      Teff = csum_ind(4,nstep_ind)
      L_Lo = csum_ind(5,nstep_ind)/3.846d+33
      logG = log10((6.67232d-8*csum_ind(1,nstep_ind)*1.989d+33)/
     +     (csum_ind(3,nstep_ind)*csum_ind(3,nstep_ind)))
      M_H = log10(csum_ind(22,nstep_ind)/
     +            csum_ind(21,nstep_ind))+1.61d0
      chisq_r = 1./Dnu_calc
c
c  OPTIMIZATION RUN
c
      if (myid .gt. 0) then
c
c  store pulsation parameters and set for faster Dnu_calc
c
         rtmp_para_el = para_el
         itmp_ipara_nsel = ipara_nsel
         itmp_ipara_iscan = ipara_iscan
         para_el = 0
         ipara_nsel = 0
         ipara_iscan = 200
c
c  select model from track using binary decision tree
c
         do i=1,400
            offset(i)=0.d0
            spacing(i)=0.d0
         enddo
         xmod_old = 0
         Dnu_calc = 10000.
         ntree = int(log10(float(ipar_nt))/log10(2.))+2

         do num_step=1,ntree

            step = int(float(nt_max)/(2.**num_step)) + 1
c
c calculate full frequency list for final step
c
            if (num_step .eq. ntree) then
               para_el = rtmp_para_el
               ipara_nsel = itmp_ipara_nsel
               ipara_iscan = itmp_ipara_iscan
            endif

            if (Dnu_calc .gt. Dnu_fit) then
               if (num_step .eq. ntree) then
                  high = int(xmod_old)
                  low = high+1
                  step=(target-offset(high))/(offset(low)-offset(high))
                  if (abs(step) .ge. 1) step=(Dnu_fit-spacing(high))/
     +                                  (spacing(low)-spacing(high))
               endif
               xmod_new = xmod_old + step
               if (int(xmod_new) .gt. nt_max) xmod_new = nt_max
            elseif (Dnu_calc .lt. Dnu_fit) then
               if (num_step .eq. ntree) then
                  low = int(xmod_old)
                  high = low-1
                  step=(offset(low)-target)/(offset(low)-offset(high))
                  if (abs(step) .ge. 1) step=(spacing(low)-Dnu_fit)/
     +                                  (spacing(low)-spacing(high))
               endif
               xmod_new = xmod_old - step
               if (int(xmod_new) .lt. 1) xmod_new = 1
            endif

 15         newseq=0

            call set_modage(am_new, par_z, icase, xmod_new, age_new, 
     *           interpol, nt_max, age_max, newseq, ierr_param)
            if(ierr_param.lt.0) then
               userff=0.0
               goto 99
            endif
c
c  calculate the mean model large frequency separation (far from obs)
c
            num_nu = 0
            sum_Dnu = 0.
            do mode=15,25
               l=obs_st(1,mode)
               l1=obs_st(1,mode+1)
               if ((l.eq.l1).and.(l.eq.0)) then
                  Dnu = obs_st(3,mode+1) - obs_st(3,mode)
                  if(Dnu.le.0) goto 23
                  sum_Dnu = sum_Dnu + Dnu
                  num_nu = num_nu + 1
               endif
 23            continue
            enddo
            Dnu_calc = sum_Dnu / num_nu
c
c  generate index of matching frequencies (model near obs)
c
            if (num_step .ge. ntree-5) then

               do obs=1,num_obs
                  best = 10000.
                  do calc=1,nobs_st
                     l=obs_st(1,calc)
                     n(1)=obs_st(2,calc)
                     if((l.eq.ell(obs)).and.(n(1).gt.0)) then
                        resid = freq(obs) - obs_st(3,calc)
                        if (abs(resid) .lt. best) then
                           best=abs(resid)
                           ominusc=resid
                           mode=calc
                        endif
                     endif
                  enddo
                  match(obs)=mode
                  if (obs .eq. 1) then
                     tstep = int(xmod_new)
                     offset(tstep) = ominusc
                  endif
               enddo
c
c  calculate mean model radial frequency and n-value
c
               sum_n=0.
               sum_fm=0.
               do obs=1,num_rnu
                  mode = match(obs)
                  sum_n = sum_n + obs_st(2,mode)
                  sum_fm = sum_fm + obs_st(3,mode)
               enddo
               n0 = sum_n / float(num_rnu)
               fm = sum_fm / float(num_rnu)
c  
c  least squares fit for model large frequency spacing
c
               sum_top=0.
               sum_bot=0.
               do obs=1,num_rnu
                  mode = match(obs)
                  sum_top=sum_top + 
     +            (obs_st(3,mode)-fm)*(obs_st(2,mode)-n0)
                  sum_bot=sum_bot + 
     +            (obs_st(2,mode)-n0)*(obs_st(2,mode)-n0)
               enddo
               Dnu_calc = sum_top / sum_bot

            endif

            spacing(int(xmod_new)) = Dnu_calc
            xmod_old = xmod_new
            age = csum_ind(2,nstep_ind)
            R_Ro = csum_ind(3,nstep_ind)/6.96d+10
            Teff = csum_ind(4,nstep_ind)
            L_Lo = csum_ind(5,nstep_ind)/3.846d+33
            logG = log10((6.67232d-8*csum_ind(1,nstep_ind)*1.989d+33)/
     +             (csum_ind(3,nstep_ind)*csum_ind(3,nstep_ind)))
            M_H = log10(csum_ind(22,nstep_ind)/
     +                  csum_ind(21,nstep_ind))+1.61d0
c            write(55,*) step,xmod_new,Dnu_calc,age,Teff,L_Lo,R_Ro
c            call flush(55)

         enddo

      endif

      write(55,'("params: ",F6.4,1X,F7.5,1X,F6.4,1X,F4.2,1X,E21.16)')
     + par_am,par_z,par_y,par_alfa,age
      write(55,'("suffix: ",A30)') trailer_par
      write(55,'("Teff, L_Lo, R_Ro: ",F7.1,1X,F6.3,1X,F6.3)') 
     + Teff,L_Lo,R_Ro
      write(55,'("logG, M_H: ",F6.3,1X,F7.3)') logG,M_H
      call flush(55)
c
c  FINAL INTERPOLATED MODEL
c
      if (myid .gt. 0) then
c
c  calculate correction for empirical surface effect
c
         r = 1./(c1*(fm/f0) - c2*(Dnu_calc/Dnu_obs))
         a0 = ((f0 - fm)*float(num_rnu)) / sum_obs
c
c  populate arrays for solar surface correction
c
         if (isflag .eq. 2) then
            frqac = par_am / (R_Ro*R_Ro*sqrt(Teff/5778.d0))*5200.d0
            do i=1,195
               hfrq(i)=frqac*h1(i)
            enddo
            do mode=match(1),match(num_obs)
               nu = obs_st(3,mode)
c  interpolate correction to matching model frequencies
               do i=1,195
                  if (hfrq(i).lt.nu) then
                     lowf=hfrq(i)
                     lowI=h2(i)
                  endif
                  if (hfrq(i).gt.nu) then
                     highf=hfrq(i)
                     highI=h2(i)
                     goto 39
                  endif
               enddo
 39            continue
               w(mode) = lowI + (highI-lowI)*(nu-lowf)/(highf-lowf)
            enddo
c  set parameters of correction
            aa=0.d0
            bb=0.d0
            cc=0.d0
            dd=0.d0
            ee=0.d0
            f0=0.d0
            do i=1,num_obs
               mode=match(i)
               f0=f0+freq(i)
               aa=aa+(obs_st(3,mode)/err(i))**2
               bb=bb+(w(mode)/err(i))**2
               cc=cc+obs_st(3,mode)*w(mode)/err(i)**2
               dd=dd+obs_st(3,mode)*freq(i)/err(i)**2
               ee=ee+w(mode)*freq(i)/err(i)**2
            enddo
            r=(dd*bb-ee*cc)/(aa*bb-cc**2)
            amp=(aa*ee-cc*dd)/(aa*bb-cc**2)
            f0=f0/num_obs
c  interpolate correction to f0
            do i=1,195
               if (hfrq(i).lt.f0) then
                  lowf=hfrq(i)
                  lowI=h2(i)
               endif
               if (hfrq(i).gt.f0) then
                  highf=hfrq(i)
                  highI=h2(i)
                  goto 49
               endif
            enddo
 49         continue
            h_f0 = lowI + (highI-lowI)*(f0-lowf)/(highf-lowf)
            a0 = amp*h_f0
         endif

         write(55,'("f0, Dnu_obs: ",F8.3,1X,F7.3)') f0,Dnu_obs
         write(55,'("fm, Dnu_calc: ",F8.3,1X,F7.3)') fm,Dnu_calc
         write(55,'(" r: ",F8.6)') r
         write(55,'("a0: ",F8.4)') a0
         call flush(55)
c
c  calculate mode inertia ratio for mixed l>0 modes
c
         do mode=match(1),match(num_obs)
            l=obs_st(1,mode)
            if ((l.ge.1).and.(obs_st(2,mode).gt.0)) then
               nu = obs_st(3,mode)
               inertia = obs_st(4,mode)
c
c  linearly interpolate inertia of radial modes to each l>0 frequency
c
               do calc=match(1),match(num_obs)
                  l=obs_st(1,calc)
                  if ((l.eq.0).and.(obs_st(2,mode).gt.0)) then
                     rnu = obs_st(3,calc)
                     rinertia = obs_st(4,calc)
                     if (rnu .LT. nu) then
                        lowf = rnu
                        lowI = rinertia
                     endif
                     if (rnu .GT. nu) then
                        highf = rnu
                        highI = rinertia
                        goto 59
                     endif
                  endif
               enddo
 59            continue
               rinertia = lowI + (highI-lowI)*(nu-lowf)/(highf-lowf)
               mixed(mode) = rinertia / inertia
            endif

         enddo
c
c  match to the individual frequencies 
c
         sum_rsq = 0.
         do obs=1,num_obs
            mode = match(obs)
            nn=obs_st(2,mode)
            l=obs_st(1,mode)
            ll=l+1
c
c  calculate surface correction and scale for mixed l>0 modes
c
            sys = a0*(obs_st(3,mode)/f0)**bexp
            if (isflag .eq. 2) then
               sys = amp*w(mode)
            endif
            if (l.ge.1) then
               sys=sys*mixed(mode)
            endif
c
c  apply surface correction (or not)
c
            fcalc_sys = obs_st(3,mode)
            if (isflag .ge. 1) then
               fcalc_sys=fcalc_sys + sys
            endif
c
c  weight the fit (or not)
c
            if (iwflag .eq. 0) then
               weight = err(obs)
            elseif (iwflag .eq. 4) then
               weight = sqrt(err(obs)*err(obs) + 0.11662225*sys*sys)
            else
               weight = sqrt(err(obs)*err(obs) + 0.25*sys*sys)
            endif

            resid = (freq(obs)-fcalc_sys)/weight
            write(55,'(2(I2,1X),F8.2,1X,F5.2,1X,F8.3,2(1X,F7.3))') 
     + l,nn,freq(obs),err(obs),fcalc_sys,sys,resid
            call flush(55)
            sum_rsq = sum_rsq + (resid*resid)
c
c  populate arrays of matching observed and calculated modes
c
            nu_calc(nn,ll)=fcalc_sys
            nu_obs(nn,ll)=freq(obs)
            nu_err(nn,ll)=err(obs)
         enddo

         chisq_seis = sum_rsq/float(num_obs)
c
c  observed and calculated large and small separations
c
         nmax=0
         nmin=100
         do obs=1,num_obs
            mode = match(obs)
            nn=obs_st(2,mode)
            l=obs_st(1,mode)
            ll=l+1
c
c  calculate Dnu0(n)
c
            if (l .eq. 0) then
               Dnu0_calc(nn) = nu_calc(nn,ll) - nu_calc(nn-1,ll)
               Dnu0_obs(nn) = nu_obs(nn,ll) - nu_obs(nn-1,ll)
               Dnu0_err(nn) = SQRT(nu_err(nn,ll)*nu_err(nn,ll) + 
     +                   nu_err(nn-1,ll)*nu_err(nn-1,ll)) / Dnu0_obs(nn)
               if (nn .lt. nmin) nmin=nn
               if (nn .gt. nmax) nmax=nn
            endif
c
c  calculate Dnu1(n)
c
            if (l .eq. 1) then
               Dnu1_calc(nn) = nu_calc(nn,ll) - nu_calc(nn-1,ll)
               Dnu1_obs(nn) = nu_obs(nn,ll) - nu_obs(nn-1,ll)
               Dnu1_err(nn) = SQRT(nu_err(nn,ll)*nu_err(nn,ll) + 
     +                   nu_err(nn-1,ll)*nu_err(nn-1,ll)) / Dnu1_obs(nn)
            endif
         enddo
c
c  match to frequency ratios r01, r10, r02, r13
c
         l0=1
         l1=2
         l2=3
         l3=4
         num_rat=0
         num_r02=0
         sum_rsq = 0.
         sum_r02 = 0.
         do nn=nmin+1,nmax-1
c
c  calculate d01 and r01
c
            d01_calc = 0.125*(nu_calc(nn-1,l0)-4*nu_calc(nn-1,l1)+
     +             6*nu_calc(nn,l0)-4*nu_calc(nn,l1)+nu_calc(nn+1,l0))
            d01_obs = 0.125*(nu_obs(nn-1,l0)-4*nu_obs(nn-1,l1)+
     +            6*nu_obs(nn,l0)-4*nu_obs(nn,l1)+nu_obs(nn+1,l0))
            d01_err = 0.125*SQRT(nu_err(nn-1,l0)*nu_err(nn-1,l0)+
     +           16*nu_err(nn-1,l1)*nu_err(nn-1,l1)+
     +           36*nu_err(nn,l0)*nu_err(nn,l0)+
     +           16*nu_err(nn,l1)*nu_err(nn,l1)+
     +           nu_err(nn+1,l0)*nu_err(nn+1,l0)) / d01_obs

            r01_calc = d01_calc/Dnu1_calc(nn)
            r01_obs = d01_obs/Dnu1_obs(nn)
            r01_err = r01_obs*SQRT(d01_err*d01_err +
     +                Dnu1_err(nn)*Dnu1_err(nn))
            if (r01_obs.gt.0 .and. r01_obs.lt.1) then
               resid = (r01_obs-r01_calc)/r01_err
               write(55,'("r01",1X,F7.2,3(1X,F10.8))')
     + nu_obs(nn,l0),r01_obs,r01_err,r01_calc
               call flush(55)
               if (nn .lt. nmax) then
                  num_rat = num_rat + 1
                  sum_rsq = sum_rsq + (resid*resid)
               endif
            endif
c
c  calculate d10 and r10
c
            d10_calc = -0.125*(nu_calc(nn-1,l1)-4*nu_calc(nn,l0)+
     +           6*nu_calc(nn,l1)-4*nu_calc(nn+1,l0)+nu_calc(nn+1,l1))
            d10_obs = -0.125*(nu_obs(nn-1,l1)-4*nu_obs(nn,l0)+
     +          6*nu_obs(nn,l1)-4*nu_obs(nn+1,l0)+nu_obs(nn+1,l1))
            d10_err = 0.125*SQRT(nu_err(nn-1,l1)*nu_err(nn-1,l1)+
     +           16*nu_err(nn,l0)*nu_err(nn,l0)+
     +           36*nu_err(nn,l1)*nu_err(nn,l1)+
     +           16*nu_err(nn+1,l0)*nu_err(nn+1,l0)+
     +           nu_err(nn+1,l1)*nu_err(nn+1,l1)) / d10_obs

            r10_calc = d10_calc/Dnu0_calc(nn+1)
            r10_obs = d10_obs/Dnu0_obs(nn+1)
            r10_err = r10_obs*SQRT(d10_err*d10_err +
     +                 Dnu0_err(nn+1)*Dnu0_err(nn+1))
            if (r10_obs.gt.0 .and. r10_obs.lt.1) then
               resid = (r10_obs-r10_calc)/r10_err
             write(55,'("r10",1X,F7.2,3(1X,F10.8))')
     + nu_obs(nn,l1),r10_obs,r10_err,r10_calc
               call flush(55)
               if (nn .lt. nmax) then
                  num_rat = num_rat + 1
                  sum_rsq = sum_rsq + (resid*resid)
               endif
            endif
c
c  calculate d02 and r02
c
            d02_calc(nn) = nu_calc(nn,l0) - nu_calc(nn-1,l2)
            d02_obs(nn) = nu_obs(nn,l0) - nu_obs(nn-1,l2)
            d02_err(nn) = SQRT(nu_err(nn,l0)*nu_err(nn,l0) + 
     +           nu_err(nn-1,l2)*nu_err(nn-1,l2)) / d02_obs(nn)
            r02_calc = d02_calc(nn)/Dnu1_calc(nn)
            r02_obs = d02_obs(nn)/Dnu1_obs(nn)
            r02_err = r02_obs*SQRT(d02_err(nn)*d02_err(nn) +
     +                Dnu1_err(nn)*Dnu1_err(nn))
            if (r02_obs.gt.0 .and. r02_obs.lt.1) then
               resid = (r02_obs-r02_calc)/r02_err
               write(55,'("r02",1X,F7.2,3(1X,F10.8))')
     + nu_obs(nn,l0),r02_obs,r02_err,r02_calc
               call flush(55)
               if (nn .lt. nmax) then
                  num_r02 = num_r02 + 1
                  sum_r02 = sum_r02 + (resid*resid)
               endif
            endif
c
c  calculate d13 and r13
c
            d13_calc(nn) = nu_calc(nn,l1) - nu_calc(nn-1,l3)
            d13_obs(nn) = nu_obs(nn,l1) - nu_obs(nn-1,l3)
            d13_err(nn) = SQRT(nu_err(nn,l1)*nu_err(nn,l1) + 
     +           nu_err(nn-1,l3)*nu_err(nn-1,l3)) / d13_obs(nn)
            r13_calc = d13_calc(nn)/Dnu0_calc(nn+1)
            r13_obs = d13_obs(nn)/Dnu0_obs(nn+1)
            r13_err = r13_obs*SQRT(d13_err(nn)*d13_err(nn) +
     +                Dnu0_err(nn+1)*Dnu0_err(nn+1))
            if (r13_obs.gt.0 .and. r13_obs.lt.1) then
               resid = (r13_obs-r13_calc)/r13_err
               write(55,'("r13",1X,F7.2,3(1X,F10.8))')
     + nu_obs(nn,l1),r13_obs,r13_err,r13_calc
               call flush(55)
               if (nn .lt. nmax) then
                  num_r02 = num_r02 + 1
                  sum_r02 = sum_r02 + (resid*resid)
               endif
            endif

         enddo

         chisq_rat = sum_rsq/float(num_rat)
         chisq_r02 = sum_r02/float(num_r02)
c
c  add residuals from Teff, L_Lo, and R_Ro
c
         sum_rsq=0.
         do obs=num_obs+1,num_obs+5
            if (ell(obs) .eq. 5) then
               Tresid = (T_obs - Teff)/Terr
               sum_rsq = sum_rsq + (Tresid*Tresid)
            endif
            if (ell(obs) .eq. 6) then
               Gresid = (G_obs - logG)/Gerr
               sum_rsq = sum_rsq + (Gresid*Gresid)
            endif
            if (ell(obs) .eq. 7) then
               Mresid = (M_obs - M_H)/Merr
               sum_rsq = sum_rsq + (Mresid*Mresid)
            endif
            if (ell(obs) .eq. 8) then
               Lresid = (L_obs - L_Lo)/Lerr
               sum_rsq = sum_rsq + (Lresid*Lresid)
            endif
            if (ell(obs) .eq. 9) then
               Rresid = (R_obs - R_Ro)/Rerr
               sum_rsq = sum_rsq + (Rresid*Rresid)
            endif
         enddo
c
c  obs: nu's + L + T | par: M + Z + Y + a + t
c
         chisq_spec = sum_rsq/float(nonseis)
         if (par_y .lt. 0.248) then
            penalty = 100.*(0.248-par_y)
c         penalty = 100.*(par_xxh-0.752 + 2.4*par_z)
            chisq_spec = chisq_spec + penalty*penalty
         endif
         write(55,'("chisq(seis,r010,r02,spec): ",4(2X,F7.3))')
     +            chisq_seis,chisq_rat,chisq_r02,chisq_spec
         call flush(55)
         if (iwflag .eq. 3) then
            chisq_r = (chisq_seis+chisq_rat+chisq_spec)/3.
         elseif (iwflag .eq. 2) then
            chisq_r = 0.25*(chisq_seis+chisq_rat+chisq_r02+chisq_spec)
         else
            chisq_r = 0.66666667*chisq_seis + 0.33333333*chisq_spec
         endif
         if (isflag .eq. 0) then
            anchor = (freq(1) - obs_st(3,match(1)))**2/err(1)**2
            chisq_r = chisq_rat*float(num_rat) + 
     +                chisq_spec*float(nonseis) + 
     +                chisq_r02*float(num_r02) + anchor
            chisq_r = chisq_r/(num_rat+num_r02+nonseis+1 - 5) 
         endif

      endif
      chisq(1) = chisq_seis
      chisq(2) = chisq_rat
      chisq(3) = chisq_r02
      chisq(4) = chisq_spec

      userff = 1./chisq_r
      if (isnan(userff)) userff=0.0

      if (myid .gt. 0) then
         write(55,'("fitness: ",F12.8)') userff
         call flush(55)
      else
         write(55,'("Delta_nu: ",F12.8)') userff
         call flush(55)
      endif

 99   close(istdou)
      close(55)

      return
      end

***************************************************************************
      subroutine write_evol(myid,idif,iovs,a_ovs)

      implicit double precision (a-h, o-z)
      integer length
      character*80 eprgdir
      character*12 fname
      common/cvr_param/ par_am, par_z, par_agefin, par_rsfin,
     *  par_alsfin, par_zxsfin, par_xmdtrl,
     *  par_xxh, par_fdgopl, par_tlopfg, par_dtlopf, par_tlopf1,
     *  par_fcno, par_agehe3, par_xrz12, par_xrz13, par_xrz14,
     *  par_xrz16, par_tprfct, par_alfa, par_etac, par_phc, par_alpove,
     *  par_alpovc, par_velrot
      common/cvi_param/ ipar_nt, ipar_icsove, ipar_icsovc, ipar_isprot,
     *  ipar_icmout, ipar_istosc, ipar_isetos
      common/cevlio/ isetos, iastr

 1    format(".evolrin.",i3.3)
      write(fname,1) myid
      open(55,file=fname,status='unknown')

      itmass = INT(10000.*(par_am)+0.5)

      write(55,'(A)') "8 'ttt.l3b.prt'"
      write(55,'(A)') "2 'emdl.0100.Z2.01.s.0'"
 2    format("3 'emdl/emdl.d",i3.3,"'")
      write(55,2) myid
      write(55,'(A)') "12 'Z.proffitt'"
c should work on Stampede
      call getenv("EPRGDIR",eprgdir)
c fallback for Stampede
      if (length(eprgdir).eq.0) 
     + eprgdir="/work/01038/gridamp/evolpack"
c 3    format("13 '",a,"/opac/ghwd-v11.gn93_ax94'")
c     A&F94 (above), Ferg05 (below)
 3    format("13 '",a,"/opac/ghwd-v11.gn93_ax05'")
      write(55,3) eprgdir(1:length(eprgdir))
      write(55,'(A)') "14 'hvabund.g91'"

      if (iastr .eq. -1) then
         write(55,'(A)') "15 '0'"
         write(55,'(A)') "16 '0'"
         write(55,'(A)') "17 '0'"
         write(55,'(A)') "18 '0'"
         write(55,'(A)') "10 '0'"
         write(55,'(A)') "11 '0'"
         write(55,'(A)') "19 '0'"
         write(55,'(A)') "22 '0'"
         write(55,'(A)') "67 '0'"
         write(55,'(A)') "99 '0'"
      else
         write(55,'(A)') "15 'amdl/amdl.d'"
         write(55,'(A)') "16 '0'"
         write(55,'(A)') "17 'gong/gong.d'"
         write(55,'(A)') "18 'gong/gsum.d'"
         write(55,'(A)') "10 '0'"
         write(55,'(A)') "11 'csum.d'"
         write(55,'(A)') "19 'cusum.d'"
         write(55,'(A)') "22 'amdl/dgamma1.d'"
         write(55,'(A)') "67 'ttt/prt.d'"
         write(55,'(A)') "99 'evol-file.log'"
      endif

      write(55,'(A)') "-1 ''"

      write(55,'(A)') "1"
 4    format(a,"/liv-eos.05/EOS2005.Z020")
      write(55,4) eprgdir(1:length(eprgdir))
      write(55,'(A)') "0.02"

      if (iovs.eq.0) then 
         write(55,'(A)') 
     + " dsn.mod.tri.eos.opa.eng.rot.dif.bcs.con.int.msh.tst.out.dgn"
c include overshoot (below)
      elseif (iovs.eq.1) then 
         write(55,'(A)')
     +" dsn.mod.tri.eos.opa.eng.ovs.rot.dif.bcs.con.int.msh.tst.out.dgn"
      endif

c                  istcon istart
      write(55,'(A)') " 0 1"
      write(55,'(A)') ",,,,,,,,,,,,,,,,,,,,,,,,,,,,"
      write(55,'(A)') ",,,,,,,,,,,,,,,"
      write(55,'(A)') ",,19,22,,,,,,,,,,,,,,"
      if (iastr .eq. -1) then
         write(55,'(A)') ",,-1,,,,,,,,,,,,,,,,,,"
      else
         write(55,'(A)') ",,67,,,,,,,,,,,,,,,,,,"
      endif
c                    #M          #Z
 5    format(" 0.0",i5.5,"0e2 ",f7.5," 601   Model S")
      write(55,5) itmass,par_z
c                    #nt
 6    format(" .0 ",i3.3,"  1.e15 0")
      write(55,6) ipar_nt

      write(55,'(A)') " 0 7 6.959900E+10 3.846000E+33 ,, 1.00E-06"
c                       #X
 7    format(" 3 1 2 ",f5.4,"  1 2 0 2  Model S  Liv05")
      write(55,7) par_xxh

      write(55,'(A)') " 6 1 0 0 0"
      write(55,'(A)') " 1 2 1 1 1"
      write(55,'(A)') " 18 .73307800 3.5000 5.5 -12.0 .0"
      write(55,'(A)') " 8.0019 -8.0 6.0 -5.0 13 ,4,,"
      write(55,'(A)') " 0 .0 6.3019 .3001 5.0,,"
      write(55,'(A)') " 0 0 .0 1.0 2.00E-02"
c      write(55,'(A)') " .2490 1 1.0 2 5.000E+07 1 6 0    Model S"
c BP95 rates (above), NACRE rates (below)
      write(55,'(A)') " .2490 1 1.0 2 5.000E+07 1 8 0    Model S"
      write(55,'(A)') " ,,,,0.51537,,,,,,,,,,"
      write(55,'(A)') " ,,,,,,,,,,,,,,,,,,"
      write(55,'(A)') " ,,,,,,,,,,,,,,,,,,"
c              #alpha
 8    format(2x,f4.2," .15713485 2.0     Model S, Liv05, v. 11")
      write(55,8) par_alfa
c                       imixcr
      write(55,'(A)') " 1000112,,,,,,,,,,,,,,"
c include overshoot
      if (iovs.eq.1) then
      write(55,'(A)') " ,,,,,,,,,,,,,,,,,,,,,,,,"
c                      #alpha_ov
 9    format(1x,",1,,",f4.2)
      write(55,9) a_ovs
      endif
c                          #v_rot
      write(55,'(A)') " 0, 2.e5"
c             #idif (0=off, 1=He, 2=He+Z)
 10   format(1x,i1,1x,"0,,,2000.,,,,,,,,,,,,,,,,,,,,,,,,")
      write(55,10) idif
c                     #ismdif
      write(55,'(A)') " 0,,,,,,,,,,"
      write(55,'(A)') ".5000 .5000 .5000 0.5 1.000 0.5,,,,,,,,,,,,,"
c     + "  , , , ,.5,1,.5, ,1,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,"
      write(55,'(A)') 
c     #ntaubc       #tmnbc
     + " 61 1.0 1.0 1.00E-03 .41169819 -9.5 5   1    Model S"
c     + " 81 1.0 1.0 1.00E-04 .41169819 -9.5 5   1    Model S"
      write(55,'(A)') 
     + " 1.0361 -.3134 2.44799995 -.29589999 30.0,,,,,,,,,,"
      write(55,'(A)') " 1  2"
c                                nit
      write(55,'(A)') " 1.00E-08 30 .3000 6 0 0 0"
      write(55,'(A)') " .5000 .5000 .5000 1.0 .5000 1.0,,,,,,,,,,,,,"
      write(55,'(A)') " 1 0"
      write(55,'(A)') " 21 .021 1.0 5.0 3.0 5. .3001   0.1   3"
      write(55,'(A)') " ,,500.,,,,,,,,,,,,,,,,,,,,,,"
      write(55,'(A)') " 1 0     std"
      write(55,'(A)') 
     + " 2.00E+12 2.00E+16 1.00E+09 .0500 .1000 100.0 1.0,,,,"
      write(55,'(A)') "  0 10 0 0 0 20 0"
      write(55,'(A)') " 1 0 0 0 0"
      write(55,'(A)') " 1 0 0 0 0 0 1 0"
      write(55,'(A)') " 0 0  0"
      write(55,'(A)') " 0 0"

      close(55)

      return
      end

***************************************************************************
      subroutine write_rdist(myid)

      implicit double precision (a-h, o-z)

      character*15 fname

 1    format(".redistrbin.",i3.3)
      write(fname,1) myid
      open(55,file=fname,status='unknown')

      write(55,'(A)') "2402,,,"
      write(55,'(A)') ",,9,60,0.0002,,,"
      write(55,'(A)') "0.001,1,0.005,,,0.01,,,,,,,,,,,"
      write(55,'(A)') "30,,,,,,,,,"
      write(55,'(A)') ",,,,,,,,,,,,,,,"

      close(55)

      return
      end

***************************************************************************
      subroutine write_adi(myid)

      implicit double precision (a-h, o-z)

      character*13 fname
      common/cevlio/ isetos, iastr

 1    format(".adiplsin.",i3.3)
      write(fname,1) myid
      open(55,file=fname,status='unknown')

      if (iastr .eq. -1) then
         write(55,'(A)') "39 '0'"
         write(55,'(A)') "32 '0'"
         write(55,'(A)') "31 '0'"
         write(55,'(A)') "35 '0'"
      else
         write(55,'(A)') "4 'amdl/amde'"
         write(55,'(A)') "39 'ttt/ttt.adipls.prt'"
         write(55,'(A)') "32 'osc/rotker'"
         write(55,'(A)') "31 'osc/agsm'"
         write(55,'(A)') "35 'ttt/ttt.adipls.ssm'"
      endif

      write(55,'(A)') "-1 ''"
      write(55,'(A)') "dsn.mod.osc.cst.int.out"
      write(55,'(A)') "31,35,4,0,,,,,,,,,,,,,,,,,,"
      write(55,'(A)') ",,,,,,,"
      write(55,'(A)') ",,,,,,,,,,,,,,,,,,,,,,"
      write(55,'(A)') "   ,4,0,1,,,,,,,"
      write(55,'(A)') "    1,  4,   ,    1,10,,,,,,,,"
      write(55,'(A)') ",2,500,2500,,,,,,,,,,,,"
      write(55,'(A)') ",,,,,,,,,,,,,,,,"
      write(55,'(A)') "6.67232e-8"
      write(55,'(A)') ",,0,,"
      write(55,'(A)') "1,,,,,,,,,,,,,,,"
      write(55,'(A)') "1,1,0.99,,,,15,,,,,,,,,,"
      write(55,'(A)') ",,,,,-1,,,,,,,,,,,,,,"
      write(55,'(A)') "39,10,,12,,,,,,,,,"
      write(55,'(A)') "1,,,,3,,,,,,,,,,,"
      write(55,'(A)') "11,,,,,,,,,,"
      write(55,'(A)') "10010,,,,,-5000,100,,,,,,,,,,,,,,,,,,,,,,,,"

      close(55)

      return
      end

***************************************************************************
      function isnan(f)

      implicit none

      real f
      logical isnan

       if (f.eq.f) then
         isnan=.FALSE.
       else
         isnan=.TRUE.
       endif

       return
       end

***************************************************************************

