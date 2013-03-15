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
     *  in_evol*80, in_rdist*80, in_adi*80, tfile*80
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
      common/par_control/ itcase, icase_trailer, xseq_trailer, 
     *  header_par
      common/cofile/ nfiles, idsfil(99), file(99), iopen(99), filess(99)
      common/nmbmsh/ nn_evol, nn_adi, ivar_adi
      common/xnwvar/ x_adi(1)
      common/anwvar/ data_adi(8), aa_adi(istrmx,1)
      common/cobs_param/ icobs_st, nobs_st, obs_st(10,1)
      common/cevlio/ isetos, iastr
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
c  input file name with starting evolution model for trial grid
c
      read(istdin,*) tfile
c
c  input case for trial grid, and masses and heavy-element abundances 
c  for possibly setting grid of trial model
c  If icaset .le. 0 use case for present run.
c
      read(istdin,*) icaset, amtg1, amtg2, damtg, ztg1, ztg2, nztg
c
c  general setups of storage and initial parameters
c
      call setups_package(in_evol, in_rdist, in_adi, ierr_param)
      if(ierr_param.lt.0) stop 'Error in setups_package'
c
      itcase=1
c
c  test for existence of appropriate trial grid and read in
c  catalogue
c
      if(icaset.le.0) icaset=icase
c
      call test_gridz(icaset, itrial_err)
c
c If not available, set grid of trial models
c
      xhtg1=0.7
      xhtg2=0.7
      dxhtg=0
      if(itrial_err.lt.0) then
        if(istdpr.gt.0)  then
	  write(istdpr,*) 'Call trial_gridz with'
          write(istdpr,*) 'amtg1, amtg2, damtg, ztg1, ztg2, nztg,',
     *      'xhtg1, xhtg2, dxhtg, tfile ='
          write(istdpr,*) amtg1, amtg2, damtg, ztg1, ztg2, nztg,
     *      xhtg1, xhtg2, dxhtg, tfile 
        end if
        call trial_gridz(amtg1, amtg2, damtg, ztg1, ztg2, nztg,
     *    xhtg1, xhtg2, dxhtg, tfile, icaset, header_par)
      end if
c
c  start setting models of specified parameters
c
   30 read(istdin,*,end=90) am_new, age_new, z_new, xxh_new
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
	par_z=z_new
	par_xxh=xxh_new
	ipar_nt=1000
	ipar_nt=59
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
