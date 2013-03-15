      program main
c
c  NOTE: all setups have been moved to s/r setups_main.
c
c  Modified 19/12/07, increasing size of arrays in common/cofile/
c
      implicit double precision (a-h,o-z)
      character trailer_par*80, id*1, header_par*80, settrailer*80,
     *  strcompr*80, file*80, filess*80, par_trial*80, set_trial*80
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
      common/cofile/ nfiles, idsfil(99), file(99), iopen(99), filess(99)
c
c  dataset designators
c
      common/cdsdsg/ idstrl,idshvz,idszab,idsevl,idssum,idscen,idsnov,
     *  idsoov,idsgng,idsgsm,idstm1,idstm2,iducen,iddgm1
c
      common/cversn/ ivrbcs,ivrrhs,ivreos,ivropc,ivreng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      call setups_main
c
c  input header and case for files for parameter outputs
c
      read(istdin,*) header_par
      read(istdin,*) icase
c
c  input trial mass
      read(istdin,*) tmass
c
c  input masses for setting grid of trial model
c
      read(istdin,*) amtg1, amtg2, damtg
c
      i_paramset = 1
c
      call mnevol(i_paramset, ierr_param)
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'Exit first call of mnevol'
c
c  test for existence of appropriate trial grid and read in
c  catalogue
c
      call test_grid(icase, itrial_err)
c
c If not available, set grid of trial models
c
      if(itrial_err.lt.0) 
     *  call trial_grid(amtg1, amtg2, damtg, tmass, icase, header_par)
c
c  loop over `random' masses to test setting of trial
c
   20 read(istdin,*,end=90) am_new
	if(istdpr.gt.0) write(istdpr,110)
        par_am=am_new
        par_trial=set_trial(par_am, par_z)
        trailer_par=strcompr(settrailer(par_am, par_z, icase, -1.d0, 0))
        call mnevol(i_paramset, ierr_param)
	call print_par(header_par)
        call dump_csum('csum1'//trailer_par, 'bcsum1'//trailer_par)
	if(ierr_param.lt.0) stop 'Error in call of mnevol'
c 
        go to 20
c
   90 continue
      stop
  110 format(//1x,70('*')//
     *  ' New call of mnevol'//
     *  1x,70('*')//)
      end
      subroutine trial_grid(amtg1, amtg2, damtg, tmass, icase, 
     *  header_par)
c
c  set grid of trial models between masses amtg1 and amtg2, with step
c  damtg; tmass should give the mass of the original trial model,
c  which will be used as point of departure.
c
c  icase is the the case number for the run; for the trial models 
c  this is augmented by 1000 to flag for trial-model output
c
c  Original version: 13/12/04
c
      implicit double precision(a-h, o-z)
      character header_par*(*)
      character trailer_par*80, id*1, settrailer*80,
     *  strcompr*80, file*80, par_trial*80, tg_files*80, cata_file*80,
     *  int2str*4, filess*80
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
      common/cofile/ nfiles, idsfil(99), file(99), iopen(99), filess(99)
      common/ctgam_param/ n_tg, icase_tg, amtg(1000),tg_files(1000)
c
c  dataset designators
c
      common/cdsdsg/ idstrl,idshvz,idszab,idsevl,idssum,idscen,idsnov,
     *  idsoov,idsgng,idsgsm,idstm1,idstm2,iducen,iddgm1
c
      common/cversn/ ivrbcs,ivrrhs,ivreos,ivropc,ivreng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      i_paramset = 1
c
      nt_orig=ipar_nt
      ipar_nt=0
c
c  set mesh in masses for trial grid
c  and locate starting point
c
      n_tg=(amtg2 - amtg1)/damtg + 1.5
      dam_start=amtg2
      do n=1,n_tg
        amtg(n)=amtg1+damtg*(n-1)
        if(abs(amtg(n)-tmass).le.dam_start) then
          dam_start=abs(amtg(n)-tmass)
          n_start=n
        end if
      end do
c
      icase_tg=icase+1000
c
      do n=n_start,1,-1
        par_am=amtg(n)
        trailer_par=
     *    strcompr(settrailer(par_am, par_z, icase_tg, -1.d0, 0))
        call mnevol(i_paramset, ierr_param)
	call print_par(header_par)
	if(ierr_param.lt.0) stop 'Error in call of mnevol'
c 
c  get prefix for emdl file for use in subsequent trial model file
c
        call stfile(idsevl,nfsevl)
        par_trial=strcompr(file(nfsevl)//trailer_par)
        tg_files(n)=par_trial
        write(istder,*) 'End setting',par_trial
      end do
c
      if(n_start.lt.n_tg) then
        par_trial=tg_files(n_start)
        do n=n_start+1,n_tg
          par_am=amtg(n)
          trailer_par=
     *      strcompr(settrailer(par_am, par_z, icase_tg, -1.d0, 0))
          call mnevol(i_paramset, ierr_param)
	  call print_par(header_par)
	  if(ierr_param.lt.0) stop 'Error in call of mnevol'
c 
c  get prefix for emdl file for use in subsequent trial model file
c
          call stfile(idsevl,nfsevl)
          par_trial=strcompr(file(nfsevl)//trailer_par)
          tg_files(n)=par_trial
          write(istder,*) 'End setting',par_trial
        end do
      end if
c
c  output catalogue of trial models 
c
      iunit=99
      cata_file=strcompr('trial.cata.'//int2str(icase_tg))
      open(iunit,file=cata_file,status='unknown')
c
      write(iunit,'(a, i6/a/a/a)') 
     *  '#  Trial models, case =',icase_tg,'#',
     *  '#  case, M/Msun, Z, file:','#'
      do i=1,n_tg
	ll=length(tg_files(i))
	write(iunit,'(i6,1p2e15.7,1x,a)') icase_tg, amtg(i), par_z,
     *    tg_files(i)(1:ll)
      end do
      close(iunit)
c
c  restore number of timesteps
c
      ipar_nt = nt_orig
      return
      end
      subroutine test_grid(icase, itrial_err)
c
c  test for presence of trial grid with appropriate icase and,
c  if so, read in catalogue
c  This surely will have to be substantially extended later
c
c  Original version: 20/12/04
c
      implicit double precision(a-h, o-z)
      character strcompr*80, cata_file*80, int2str*4, tg_files*80
      logical file_exists
      common/ctgam_param/ n_tg, icase_tg, amtg(1000),tg_files(1000)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      icase_tg=icase+1000
c
c  output catalogue of trial models 
c
      cata_file=strcompr('trial.cata.'//int2str(icase_tg))
c..      i=system('ls -l '//cata_file)
c..      if (i.ne.0) then
      inquire(file=cata_file,exist=file_exists)
      if(file_exists) then
        iunit=99
        open(iunit,file=cata_file,status='old')
	call skpcom(iunit)
	n=1
   20   read(iunit,'(i6,1p2e15.7,1x,a)',end=25) 
     *    icase_rd, amtg(n), z_rd, tg_files(n)
        n=n+1
	go to 20
c
   25   n_tg=n-1
	if(istdpr.gt.0) 
     *  write(istdpr,'('' Trial data read in''/(1pe13.5,1x,a))')
     *    (amtg(n), tg_files(n), n=1,n_tg)
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
      character*(*) function set_trial(am, z)
c
c  sets file name for closest trial model for given mass am and
c  heavy-element abundance z (note: z is so far not used).
c
c  Original version: 13/12/04
c
      implicit double precision(a-h, o-z)
      character tg_files*80
      common/ctgam_param/ n_tg, icase_tg, amtg(1000),tg_files(1000)
c
c  dataset designators
c
      common/cdsdsg/ idstrl,idshvz,idszab,idsevl,idssum,idscen,idsnov,
     *  idsoov,idsgng,idsgsm,idstm1,idstm2,iducen,iddgm1
c
      common/cversn/ ivrbcs,ivrrhs,ivreos,ivropc,ivreng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      dam_fit=1.d10
      do n=1,n_tg
        if(abs(amtg(n)-am).le.dam_fit) then
          dam_fit=abs(amtg(n)-am)
          n_fit=n
        end if
      end do
c
      set_trial = tg_files(n_fit)
c
      return
      end
