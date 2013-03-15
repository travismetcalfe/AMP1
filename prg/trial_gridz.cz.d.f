      subroutine trial_gridz(amtg1, amtg2, damtg, ztg1, ztg2, nztg,
     *  xhtg1, xhtg2, dxhtg, tfile, icase, header_par)
c
c  set grid of trial models between masses amtg1 and amtg2, with step
c  damtg; heavy-element abundance between ztg1 and ztg2 in nztg steps
c  (step in hydrogen abundance may be implemented later)
c  tfile should give file name for the original trial model,
c  which will be used as point of departure.
c
c  icase is the the case number for the run; for the trial models 
c  this is augmented by 1000 to flag for trial-model output
c
c  Original version: 5/10/06
c
c  Modified 19/12/07, increasing size of arrays in common/cofile/
c
      implicit double precision(a-h, o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.cz.d.incl'
      parameter(naztmx = nspcmx+3, iymax = 2*ivarmx+1, 
     *  kaxst = 3+naztmx, nbmax = nspcmx+2,
     *  nam_max=201, nz_max=10)
      character header_par*(*), tfile*(*)
      character trailer_par*80, id*1, settrailer*80,
     *  strcompr*80, file*80, par_trial*80, tg_files*80, cata_file*80,
     *  int2str*4, header_parc*80, trial_par_seq*80
      dimension xtr(nnmax), ytr(iymax,nnmax), datmod_tr(nrdtmd),
     *  ntdmod_tr(nidtmd), bccoef_tr(nbcprv)
      common/cvr_param/ par_am, par_z, par_agefin, par_rsfin, 
     *  par_alsfin, par_zxsfin, par_xmdtrl,
     *  par_xxh, par_fdgopl, par_tlopfg, par_dtlopf, par_tlopf1,
     *  par_fcno, par_agehe3, par_xrz12, par_xrz13, par_xrz14,
     *  par_xrz16, par_tprfct, par_alfa, par_etac, par_phc, par_alpove,
     *  par_alpovc, par_velrot
      common/cvi_param/ ipar_nt, ipar_icsove, ipar_icsovc, ipar_isprot,
     *  ipar_icmout, ipar_istosc, ipar_isetos
      common/trial_param/ par_trial
      common/trl_param/ trailer_par
      common/par_control/ itcase, icase_trailer, xseq_trailer, 
     *  header_parc, trial_par_seq
      common/cofile/ nfiles, idsfil(99), file(99), iopen(99)
      common/ctgam_param/ n_tg, nztgl, nxtgl, icase_tg, 
     *  amtg(nam_max,nz_max),
     *  ztg(nam_max,nz_max), tg_files(nam_max,nz_max)
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev, syear
c
c  dataset designators
c
      common/cdsdsg/ idstrl,idshvz,idszab,idsevl,idssum,idscen,idsnov,
     *  idsoov,idsgng,idsgsm,idstm1,idstm2,iducen
c
      common/cversn/ ivrbcs,ivrrhs,ivreos,ivropc,ivreng
      common/cevlio/ isetos, iastr
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      common/cstdio_in/ istdin_in, istdpr_in
c
c  hard-code idiag (mainly for output of parameters)
c
      idiag=0
      if(istdpr_in.lt.0) idiag=0
c
      i_paramset = 1
c
      iastr_orig=iastr
      iastr=1
      nt_orig=ipar_nt
      ipar_nt=0
c
c  determine parameters for trial model
c
      open(99,file=tfile,status='old',form='unformatted')
      icase_tr=0
      call rdemdl(99,xtr,ytr,datmod_tr,ntdmod_tr,bccoef_tr,iymax,iform,
     *  nn_tr, nrdtmd_tr, nidtmd_tr, nbccf_tr, nvar_tr, icase_tr, icry)
      close(99)
      tmass=datmod_tr(23)/amsun
      txh=datmod_tr(2)
      tzh=datmod_tr(1)
c
      par_trial=tfile
c
      if(nztg.le.0) then
        nztgl=1
	ztgg=par_z
      else if(nztg.eq.1) then
	nztgl=1
	dztgl=0
	ztgg=ztg2
      else
	nztgl=nztg
	dztgl=log10(ztg2/ztg1)/(nztgl-1)
      end if
      if(nztgl.gt.nz_max) then
	write(istdou,'(/ '' ***** Error. no of Z points in trial ='',
     *    i3,'' exceeds maximum ='',i3/)') nztgl, nz_max
	if(istdpr.ne.istdou.and.istdpr.gt.0)
     *    write(istdou,'(/ '' ***** Error. no of Z points in trial ='',
     *      i3,'' exceeds maximum ='',i3/)') nztgl, nz_max
	stop 'in trial_gridz'
      end if
c
c  set mesh in masses for trial grid
c  and locate starting point
c
      n_tg=(amtg2 - amtg1)/damtg + 1.5
c
      if(n_tg.gt.nam_max) then
	write(istdou,'(/ '' ***** Error. no of M points in trial ='',
     *    i3,'' exceeds maximum ='',i3/)') n_tg, nam_max
	if(istdpr.ne.istdou.and.istdpr.gt.0)
     *    write(istdou,'(/ '' ***** Error. no of M points in trial ='',
     *      i3,'' exceeds maximum ='',i3/)') n_tg, nam_max
	stop 'in trial_gridz'
      end if
      dam_start=amtg2
      do n=1,n_tg
        amtg_n=amtg1+damtg*(n-1)
        if(abs(amtg_n-tmass).le.dam_start) then
          dam_start=abs(amtg_n-tmass)
          n_start=n
        end if
	do i=1,nztgl
          amtg(n,i)=amtg_n
        end do
      end do
c
c  set mesh in heavy-element abundance for trial grid
c
      if(ntgl.eq.1) then
	do n=1,n_tg
	  ztg(n,1)=ztgg
	end do
c
      else
c
c  also locate starting point in Z
c
        dz_start=ztg2
        do i=1,nztgl
	  ztgg=ztg1*10.d0**(dztgl*(i-1))
	  do n=1,n_tg
	    ztg(n,i)=ztgg
	  end do
          if(abs(ztgg-tzh).le.dz_start) then
            dz_start=abs(ztgg-tzh)
            i_start=i
          end if
        end do
      end if
c
      icase_tg=icase+1000
      icase_trailer=icase_tg
c
c  initial pass through masses, at Z(i_start)
c
      par_z=ztg(1,i_start)
      do n=n_start,1,-1
        par_am=amtg(n,1)
        trailer_par=
     *    strcompr(settrailer(par_am, par_z, icase_tg, -1.d0, itcase))
        call mnevol(i_paramset, ierr_param)
        if(idiag.gt.0) call print_par(header_par)
        if(ierr_param.lt.0) stop 'Error in call of mnevol'
c 
c  get prefix for emdl file for use in subsequent trial model file
c
        call stfile(idsevl,nfsevl)
        par_trial=strcompr(file(nfsevl)//trailer_par)
        tg_files(n,i_start)=par_trial
      end do
c
      if(n_start.lt.n_tg) then
        par_trial=tg_files(n_start,i_start)
        do n=n_start+1,n_tg
          par_am=amtg(n,1)
          trailer_par=
     *      strcompr(settrailer(par_am, par_z, icase_tg, -1.d0, itcase))
          call mnevol(i_paramset, ierr_param)
          if(idiag.gt.0) call print_par(header_par)
          if(ierr_param.lt.0) stop 'Error in call of mnevol'
c 
c  get prefix for emdl file for use in subsequent trial model file
c
          call stfile(idsevl,nfsevl)
          par_trial=strcompr(file(nfsevl)//trailer_par)
          tg_files(n,i_start)=par_trial
        end do
      end if
c
c  test for other Z values
c
      if(i_start.lt.nztgl) then
	do i=i_start+1,nztgl
          par_z=ztg(1,i)
          par_trial=tg_files(1,i-1)
          do n=1,n_tg
            par_am=amtg(n,i)
            trailer_par=
     *        strcompr(settrailer(par_am, par_z, icase_tg, -1.d0, 
     *        itcase))
            call mnevol(i_paramset, ierr_param)
            if(idiag.gt.0) call print_par(header_par)
            if(ierr_param.lt.0) stop 'Error in call of mnevol'
c 
c  get prefix for emdl file for use in subsequent trial model file
c
            call stfile(idsevl,nfsevl)
            par_trial=strcompr(file(nfsevl)//trailer_par)
            tg_files(n,i)=par_trial
          end do
        end do
      end if
c
      if(i_start.gt.1) then
	do i=i_start-1,1,-1
          par_z=ztg(1,i)
          par_trial=tg_files(1,i+1)
          do n=1,n_tg
            par_am=amtg(n,i)
            trailer_par=
     *        strcompr(settrailer(par_am, par_z, icase_tg, -1.d0, 
     *        itcase))
            call mnevol(i_paramset, ierr_param)
            if(idiag.gt.0) call print_par(header_par)
            if(ierr_param.lt.0) stop 'Error in call of mnevol'
c 
c  get prefix for emdl file for use in subsequent trial model file
c
            call stfile(idsevl,nfsevl)
            par_trial=strcompr(file(nfsevl)//trailer_par)
            tg_files(n,i)=par_trial
          end do
        end do
      end if
c
c  output catalogue of trial models 
c
      iunit=99
      cata_file=strcompr('trial.cata.'//int2str(icase_tg))
      open(iunit,file=cata_file,status='unknown')
c
      write(iunit,'(a, i6/a/a/a/a)') 
     *  '#  Trial models, case =',icase_tg,'#',
     *  '#  ',
     *  '#  no of masses, no of Z-values',
     *  '#  case, M/Msun, Z, X, file:','#'
      write(iunit,'(2i4)') n_tg, nztgl
      do i=1,nztgl
        do n=1,n_tg
          ll=length(tg_files(n,i))
          write(iunit,'(i6,1p3e15.7,1x,a)') icase_tg, amtg(n,i), 
     *      ztg(n,i), par_xxh, tg_files(n,i)(1:ll)
        end do
      end do
      close(iunit)
c
c  restore number of timesteps and flag for file output
c
      ipar_nt = nt_orig
      iastr=iastr_orig
      return
      end
      subroutine test_gridz(icase, itrial_err)
c
c  test for presence of trial grid with appropriate icase and,
c  if so, read in catalogue
c
c  Original version: 20/12/04
c
      implicit double precision(a-h, o-z)
      parameter(nam_max=201, nz_max=10)
      character strcompr*80, cata_file*80, int2str*4, tg_files*80
      logical file_exists
      common/ctgam_param/ n_tg, nztgl, nxtgl, icase_tg, 
     *  amtg(nam_max,nz_max),
     *  ztg(nam_max,nz_max), tg_files(nam_max,nz_max)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      icase_tg=icase+1000
c
c  read catalogue of trial models 
c
      cata_file=strcompr('trial.cata.'//int2str(icase_tg))
c..      i=system('ls -l '//cata_file)
c..      if (i.ne.0) then
      inquire(file=cata_file,exist=file_exists)
      if(file_exists) then
        iunit=99
        open(iunit,file=cata_file,status='old')
        call skpcom(iunit)
        read(iunit,'(2i4)') n_tg, nztgl
c
c  test for sufficient size of arrays
c
        if(nztgl.gt.nz_max) then
	  write(istdou,'(/ '' ***** Error. no of Z points in trial ='',
     *      i3,'' exceeds maximum ='',i3/)') nztgl, nz_max
	  if(istdpr.ne.istdou.and.istdpr.gt.0)
     *      write(istdou,'(/ '' ***** Error. no of Z points in trial ='',
     *        i3,'' exceeds maximum ='',i3/)') nztgl, nz_max
	  stop 'in trial_gridz'
        end if
        if(n_tg.gt.nam_max) then
	  write(istdou,'(/ '' ***** Error. no of M points in trial ='',
     *      i3,'' exceeds maximum ='',i3/)') n_tg, nam_max
	  if(istdpr.ne.istdou.and.istdpr.gt.0)
     *      write(istdou,'(/ '' ***** Error. no of M points in trial ='',
     *        i3,'' exceeds maximum ='',i3/)') n_tg, nam_max
	  stop 'in trial_gridz'
        end if
c
        read(iunit,'(i6,1p3e15.7,1x,a)') 
     *    ((icase_rd, amtg(n,i), ztg(n,i), xxh_trial,
     *    tg_files(n,i),n=1,n_tg), i=1,nztgl)
c
        if(istdpr.gt.0) then
          write(istdpr,'('' Trial data read in'')')
	  do i=1,nztgl
	    do n=1,n_tg
	      ltg=length(tg_files(n,i))
              write(istdpr,'(1p2e13.5,1x,a)')
     *          amtg(n,i), ztg(n,i), tg_files(n,i)(1:ltg)
	    end do
	  end do
	end if
        itrial_err=0
        close(iunit)
      else
        ll=length(cata_file)
        write(istder,'('' Did not find file '',a)') cata_file(1:ll)
        itrial_err=-1
        return
      end if
      return
      end
      character*(*) function set_trialz(am, z)
c
c  sets file name for closest trial model for given mass am and
c  heavy-element abundance z 
c
c  Original version: 13/12/04
c
      implicit double precision(a-h, o-z)
      parameter(nam_max=201, nz_max=10)
      character tg_files*80
      common/ctgam_param/ n_tg, nztgl, nxtgl, icase_tg, 
     *  amtg(nam_max,nz_max),
     *  ztg(nam_max,nz_max), tg_files(nam_max,nz_max)
c
c  dataset designators
c
      common/cdsdsg/ idstrl,idshvz,idszab,idsevl,idssum,idscen,idsnov,
     *  idsoov,idsgng,idsgsm,idstm1,idstm2,iducen
c
      common/cversn/ ivrbcs,ivrrhs,ivreos,ivropc,ivreng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      dam_fit=1.d10
      do n=1,n_tg
        if(abs(amtg(n,1)-am).le.dam_fit) then
          dam_fit=abs(amtg(n,1)-am)
          n_mfit=n
        end if
      end do
c
      if(nztgl.eq.1) then
	i_zfit=1
      else
        dz_fit=1.d10
        do i=1,nztgl
          if(abs(ztg(1,i)-z).le.dz_fit) then
            dz_fit=abs(ztg(1,i)-z)
            i_zfit=i
          end if
        end do
      end if
c
      set_trialz = tg_files(n_mfit, i_zfit)
c
      return
      end
