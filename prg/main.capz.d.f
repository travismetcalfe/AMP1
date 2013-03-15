      program main
c
c  NOTE: all setups have been moved to s/r setups_main.
c
c  Modified 19/12/07, increasing size of arrays in common/cofile/
c
      implicit double precision (a-h,o-z)
      include 'engenr.bz.d.incl'
      parameter(iaa_adi=11)
      character trailer_par*80, id*1, header_par*80, settrailer*80,
     *  strcompr*80, file*80, filess*80, par_trial*80, set_trial*80,
     *  in_evol*80, in_rdist*80, in_adi*80
      dimension xr_adi(nnmax), aar_adi(iaa_adi,nnmax)
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
      common/nmbmsh/ nn_evol, nn_adi, ivar_adi
      common/xnwvar/ x_adi(1)
      common/anwvar/ data_adi(8), aa_adi(istrmx,1)
      common/cobs_param/ icobs_st, nobs_st, obs_st(10,1)
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
      common/cstdio_in/ istdin_in, istdpr_in
c
c  input files for parameter input to evolution, redistribution
c  and pulsations
c
      read(istdin,*) in_evol
      read(istdin,*) in_rdist
      read(istdin,*) in_adi
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
c  general setups of storage and initial parameters
c
      call setups_package(in_evol, in_rdist, in_adi, ierr_param)
      if(ierr_param.lt.0) stop 'Error in setups_package'
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
   20 read(istdin,*,end=90) am_new, nmod
	if(am_new.le.0) go to 30
	if(istdpr.gt.0) write(istdpr,110)
        par_am=am_new
        par_trial=set_trial(par_am, par_z)
	ipar_icmout=1
	ipar_istosc=1
        trailer_par=strcompr(settrailer(par_am, par_z, icase, -1.d0, 0))
        call mnevol(i_paramset, ierr_param)
	call print_par(header_par)
        call dump_csum('csum1'//trailer_par, 'bcsum1'//trailer_par)
	if(ierr_param.lt.0) stop 'Error in call of mnevol'
c
c  redistribute mesh in last model
c
	i_paramset=1
	i_inout=0
	call srdist(i_paramset, ierr_param, i_inout,
     *    x_adi, aa_adi, data_adi, xr_adi, aar_adi, nn_adi,
     *    nnr_adi, ivar_adi, istrmx, iaa_adi)
	if(ierr_param.lt.0) stop 'Error in call of redistrb'
c
c  calculate oscillations in last model
c
	i_paramset=1
	i_inout=0
	call adipls(i_paramset, ierr_param, i_inout,
     *    xr_adi, aar_adi, data_adi, nnr_adi, ivar_adi, iaa_adi)

	if(ierr_param.lt.0) stop 'Error in call of adipls'
	call dump_obs(strcompr('osc/obs'//trailer_par))
c
c  set model no. nmod and compute frequencies for this model
c
c  set trial model from previous output
c
	call stfile(idsevl,nfsevl)
	par_trial=strcompr(file(nfsevl)//trailer_par)
	xmdtrl_orig=par_xmdtrl
	par_xmdtrl=nmod
	nt_orig=ipar_nt
	ipar_nt=1
	isetos_orig=ipar_isetos
	ipar_isetos=1
        trailer_par=strcompr(settrailer(par_am, par_z, icase,
     *    par_xmdtrl, 0))
        call mnevol(i_paramset, ierr_param)
	if(istdpr_in.ge.0) then
	  call print_par(header_par)
          call dump_csum('csum1'//trailer_par, 'bcsum1'//trailer_par)
	end if
	if(ierr_param.lt.0) stop 'Error in call of mnevol'
c
	ipar_nt=nt_orig
	ipar_isetos=isetos_orig
	par_xmdtrl=xmdtrl_orig
c
c  redistribute mesh in last model
c
	i_paramset=1
	i_inout=0
	call srdist(i_paramset, ierr_param, i_inout,
     *    x_adi, aa_adi, data_adi, xr_adi, aar_adi, nn_adi,
     *    nnr_adi, ivar_adi, istrmx, iaa_adi)
	if(ierr_param.lt.0) stop 'Error in call of redistrb'
c
c  calculate oscillations in last model
c
	i_paramset=1
	i_inout=0
	call adipls(i_paramset, ierr_param, i_inout,
     *    xr_adi, aar_adi, data_adi, nnr_adi, ivar_adi, iaa_adi)

	if(ierr_param.lt.0) stop 'Error in call of adipls'
	call dump_obs('osc/obs'//trailer_par)
c 
        go to 20
c
c  start tests of set_modage
c
   30 read(istdin,*,end=90) am_new, age_new
c
c  test for input giving age or xmod
c
	if(age_new.lt.0) then
	  xmod_new=-1.
	  nt_max=-age_new
	  write(istdou,'(//
     *      '' Call set_modage just to set evolution sequence''/)')
	else if(age_new.lt.1.d4) then
	  xmod_new=age_new
	  age_new=-1.
	  nt_max=0
	  write(istdou,'(//'' Call set_modage with am ='',f10.5,
     *      '' xmod ='', f10.5/)') am_new, xmod_new
        else
	  xmod_new=-1.
	  nt_max=0
	  write(istdou,'(//'' Call set_modage with am ='',f10.5,
     *      '' age ='', 1pe13.5/)') am_new, age_new
	end if
	ipar_nt=1000
	interpol=1
	age_max=0
	newseq=0
        call set_modage(am_new, par_z, icase, xmod_new, age_new, 
     *    interpol, nt_max, age_max, newseq, ierr_param)
        if(ierr_param.lt.0) stop 'Error in call of set_modage'
        go to 30
c
   90 continue
      stop
  110 format(//1x,70('*')//
     *  ' New call of mnevol'//
     *  1x,70('*')//)
      end
