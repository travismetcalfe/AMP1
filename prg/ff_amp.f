      function userff ( npar, data, myid )
c
c        npar = number of parameters
c  data(npar) = scaled parameter values [0.0 -> 1.0]
c        myid = parallel slave number [1 -> Nproc]
c
      implicit double precision (a-h,o-z)

      integer npar,myid,ell(100),obs,calc,high,low,match(100),tstep
      integer n(100),peak,gap
      double precision data(36),freq(100),spacing(400),err(100),target
      double precision step,log_z,avg_sys(5),Tresid,Lresid,ominusc
      double precision Teff,L_Lo,T_obs,L_obs,Terr,Lerr,offset(400)
      double precision R_Ro,Rerr,Rresid,n0,sum_n,sum_top,sum_bot
      double precision mixed(100),nu,rnu,lowf,highf,lowI,highI
      double precision rinertia,inertia,G_obs,Gerr,logG,Gresid
      double precision M_obs,Merr,M_H,Mresid,Dnu_in

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
      common /xmodage/ age
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
      if (myid .gt. 0) then
c
c  read obs.dat file
c
         nonseis=0
         open(55,file='obs.dat',status='old')
         read(55,*) num_obs,iwflag,isflag,idif,Dnu_in
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
      call write_evol(myid,idif)
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
            write(55,*) step,xmod_new,Dnu_calc,age,Teff,L_Lo,R_Ro
            call flush(55)

         enddo

      endif

      write(55,'("params: ",F6.4,1X,F7.5,1X,F6.4,1X,F4.2,1X,E21.16)')
     + par_am,par_z,par_y,par_alfa,age
      write(55,'("suffix: ",A30)') trailer_par
      write(55,*) "Teff, L_Lo, R_Ro: ",Teff,L_Lo,R_Ro
      write(55,*) "logG, M_H: ",logG,M_H
      call flush(55)
c
c  FINAL INTERPOLATED MODEL
c
      if (myid .gt. 0) then
c
c  calculate correction for surface effects
c
         r = 1./(c1*(fm/f0) - c2*(Dnu_calc/Dnu_obs))
         a0 = ((f0 - fm)*float(num_rnu)) / sum_obs

         write(55,*) "f0, Dnu_obs: ",f0,Dnu_obs
         write(55,*) "fm, Dnu_calc: ",fm,Dnu_calc
         write(55,*) " r: ",r
         write(55,*) "a0: ",a0
         call flush(55)
c
c  calculate mode inertia ratio for mixed l=1 modes
c
         do mode=match(1),match(num_obs)
            l=obs_st(1,mode)
            if ((l.eq.1).and.(obs_st(2,mode).gt.0)) then
               nu = obs_st(3,mode)
               inertia = obs_st(4,mode)
c
c  linearly interpolate inertia of radial modes to each l=1 frequency
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
                        goto 49
                     endif
                  endif
               enddo
 49            continue
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
            l=obs_st(1,mode)
c
c  scale surface correction for mixed l=1 modes
c
            sys = a0*(obs_st(3,mode)/f0)**bexp
            if (l.eq.1) then
               sys=sys*mixed(mode)
            endif
c
c  apply surface correction (or not)
c
            fcalc_sys = obs_st(3,mode)
            if (isflag .eq. 1) then
               fcalc_sys=fcalc_sys + sys
            endif
c
c  weight the fit (or not)
c
            if (iwflag .eq. 0) then
               weight = err(obs)
            else
               weight = sqrt(err(obs)*err(obs) + 0.25*sys*sys)
            endif

            resid = (freq(obs)-fcalc_sys)/weight
            write(55,*) ell(obs),obs_st(2,mode),freq(obs),
     +           err(obs),fcalc_sys,sys,resid
            call flush(55)
            sum_rsq = sum_rsq + (resid*resid)
         enddo

         chisq_seis = sum_rsq/float(num_obs)
         sum_rsq=0.
c
c  add residuals from Teff, L_Lo, and R_Ro
c
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
         write(55,'("chisq(seis,spec): ",F6.2,2X,F6.2)')
     +            chisq_seis,chisq_spec
         call flush(55)
         chisq_r = 0.66666667*chisq_seis + 0.33333333*chisq_spec
c         chisq_r = sum_rsq/float(num_obs+nonseis-5)

      endif

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
      subroutine write_evol(myid,idif)

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
c should work on Kraken
      call getenv("EPRGDIR",eprgdir)
c fallback for Frost
      if (length(eprgdir).eq.0) eprgdir="/home/gridamp/evolpack"
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

      write(55,'(A)') 
     + " dsn.mod.tri.eos.opa.eng.rot.dif.bcs.con.int.msh.tst.out.dgn"
c include overshoot (below)
c     +" dsn.mod.tri.eos.opa.eng.ovs.rot.dif.bcs.con.int.msh.tst.out.dgn"
c
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
c      write(55,'(A)') " ,,,,,,,,,,,,,,,,,,,,,,,,"
c                           #alpha_ov
c      write(55,'(A)') " ,1,,0.25"
c                          #v_rot
      write(55,'(A)') " 0, 2.e5"
c             #idif (0=off, 1=He, 2=He+Z)
 9    format(1x,i1,1x,"0,,,2000.,,,,,,,,,,,,,,,,,,,,,,,,")
      write(55,9) idif
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
         write(55,'(A)') "39 'ttt/ttt.adipls.prt'"
         write(55,'(A)') "32 'osc/rotker'"
         write(55,'(A)') "31 'osc/agsm'"
         write(55,'(A)') "35 'ttt/ttt.adipls.ssm'"
      endif

      write(55,'(A)') "-1 ''"
      write(55,'(A)') "dsn.mod.osc.cst.int.out"
      write(55,'(A)') "31,35,,0,,,,,,,,,,,,,,,,,,"
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
      write(55,'(A)') "1,,,,,,,,,,,,,,,"
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

