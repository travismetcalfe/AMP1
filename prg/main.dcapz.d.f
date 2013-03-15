      program main
c
c  NOTE: all setups have been moved to s/r setups_main.
c
c  Modified 19/12/07, increasing size of arrays in common/cofile/
c
      implicit double precision (a-h,o-z)
      include 'engenr.cz.d.incl'
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
	ipar_nt=300
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
