      subroutine set_modage(am, az, icase, xmod, age, 
     *  interpol, nt_max, age_max, newseq, ierr_param)
c
c  Set up model quantities for model characterized by mass am 
c  (in solar units), heavy-element abundance az and real model number
c  xmod or age age (in years).
c
c  Note: initial calls of mnevol, srdist and adipls must be made 
c  before calling set_modage, to read in initial input parameters.
c
c  Tests if appropriate evolution sequence has already just been
c  calculated and, if not, calculate sequence.
c  If newseq .ne. 0, force calculation of new sequence (e.g., if
c  other parameters have been reset).
c
c  If xmod .ge. 0 set model number xmod (using linear interpolation 
c  if xmod is not integer), otherwise set model of given age. In that case,
c  if interpol = 1, interpolate in sequence to model with specified age;
c  else select model with age closest to the specified age.
c
c  If nt_max .gt. 0, run evolution sequence for nt_max steps.
c  Else, if age_max .gt. 0 run sequence until first model older than
c  age_max.
c  If neither of these conditions are satisfied, set number of models
c  based on xmod, if set, or maximum age on the basis of the desired age.
c
c  Original version: 12/7/05
c
c  Modified 19/12/07, increasing size of arrays in common/cofile/
c
      implicit double precision (a-h,o-z)
      include 'engenr.cz.d.incl'
      parameter(iaa_adi=11)
      character trailer_par*80, trailer_par_seq*80, id*1, 
     *  header_par*80, settrailer*80,
     *  strcompr*80, file*80, filess*80, par_trial*80, set_trialz*80
      logical just_seq
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
     *  header_par, trailer_par_seq
      common/cofile/ nfiles, idsfil(99), file(99), iopen(99), filess(99)
      common/nmbmsh/ nn_evol, nn_adi, ivar_adi
      common/cevlio/ isetos, iastr
      common /csum_param/ icsum, nstep, csum_st(icsum_max, nstep_max)
      common /csum_indiv/ icsum_ind, nstep_ind,
     *  csum_ind(icsum_max, nstep_max)
      common/xnwvar/ x_adi(1)
      common/anwvar/ data_adi(8), aa_adi(istrmx,1)
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
      common/cstdio_in/ istdin_in, istdpr_in
c
      data am_prev, az_prev /0.d0, 0.d0/
      data init /0/
      data iseq /0/
c
      save
c
c  scale factor for testing age range
c
      agescl=1.01d0
c
      if(istdpr.gt.0) write(istder,
     *  '(''Enter set_modage, am ='',f10.5,
     *  ''  xmod ='',f10.5,''  age ='',1pe13.5)') am, xmod, age
c
c  test for just setting evolution sequence
c
      just_seq=xmod.lt.0.and.age.lt.0
c
      i_paramset=1
c
      if(init.eq.0) then
	init=1
	age_last=0.d0
      else
	age_last=csum_st(2,nstep)
      end if
c
      par_am=am
      par_z=az
c
      icmout_orig=ipar_icmout
      istosc_orig=ipar_istosc
c
c  test for making new sequence
c
      if(just_seq.or.newseq.ne.0.or.abs(am_prev-am).gt.1.d-8.or.
     *  abs(az_prev-az).gt.1.d-8.or.age.ge.age_last/agescl) then
c
c  compute evolution sequence
c
        write(istdou,'(//'' Start new evolution sequence with M ='',
     *    f10.5,''  Z ='',f10.6/)') am, az
	if(nt_max.gt.0) then
	  ipar_nt=nt_max
	  par_agefin=1.d20
        else if(age_max.gt.0) then
	  par_agefin=-age_max
        else if(xmod.ge.0) then
	  ipar_nt=xmod+3
	  par_agefin=1.d20
        else
	  par_agefin=-agescl*age
	end if
c
c  use setting the input istosc .ge. 10 to flag for forcing 
c  full oscillation output for evolution-sequence part of calculation
c
	ipar_icmout=1
	if(istosc_orig.ge.10) then
	  ipar_istosc=1
        else
	  ipar_istosc=0
        end if
c
c  clear output file for evolution models to force reopen
c
        if(iastr.gt.0) then
	  call stfile(idsevl,nfsevl)
          filess(nfsevl)=''
        end if
c
        par_trial=set_trialz(par_am, par_z)
	icase_trailer=iseq+icase
	xseq_trailer=-1.d0
        trailer_par=strcompr(settrailer(am, az, icase_trailer, 
     *    xseq_trailer, itcase))
        trailer_par_seq=trailer_par
        call mnevol(i_paramset, ierr_param)
        if(istdpr_in.ge.0) call dump_csum('csum1'//trailer_par,
     *    'bcsum1'//trailer_par)
	if(ierr_param.lt.0) then
	  write(istder,*) 'Error in call of mnevol'
          ipar_icmout=icmout_orig
          ipar_istosc=istosc_orig
	  return
        end if
	am_prev=am
	az_prev=az
c
c  test for just setting evolution sequence
c
	if(just_seq) then 
          ipar_icmout=icmout_orig
          ipar_istosc=istosc_orig
	  return
       end if
c
      end if
c
c  test for type of model setting
c
      if(xmod.ge.0) then 
	xmdtrl=xmod
c
c  test for interpolation in age
c
      else if(interpol.eq.1) then
c
c  locate models spanning age in sequence and set
c  corresponding xmdtrl
c
        nmod=-1
        do n=1,nstep
	  if(age.lt.csum_st(2,n).and.nmod.eq.-1) nmod=n-1
        end do
        xmdtrl=nmod
     *    +(age-csum_st(2,nmod))/(csum_st(2,nmod+1)-csum_st(2,nmod))	
c
c  test that bracketing ages have been found
c
	if(nmod.eq.-1) then
	  write(istdou,'(//
     *      ''  ***** Error in set_modage.''/
     *      ''        Age ='',1pe13.5,
     *      '' exceeds maximum age in sequence ='',e13.5)')
     *      age, csum_st(2,nstep)
	  ierr_param=-1
	  return
        end if
c
        write(istdou,'(//'' Target age ='',1pe13.5,
     *    '' Interpolate with xmdtrl ='',0pf10.4/
     *    '' between model no.'',i4,'' age ='',1pe13.5,
     *    '' and model no.'',i4,'' age ='',1pe13.5/)')
     *    age,xmdtrl,nmod,csum_st(2,nmod),nmod+1,csum_st(2,nmod+1)
c
      else
c
c  locate model closest in age in sequence
c
        dmin_age=age
        do n=1,nstep
	  if(abs(age-csum_st(2,n)).lt.dmin_age) then
	    dmin_age=abs(age-csum_st(2,n))
	    nmod=n
          end if
        end do
	xmdtrl=nmod

c
        write(istdou,'(//'' Target age ='',1pe13.5,
     *    '' Closest age ='',e13.5,'' for model no.'',i4)')
     *    age,csum_st(2,nmod),nmod
c
      end if
c
c  set model no. nmod and compute frequencies for this model
c
c  set trial model from previous output
c
      if(iastr.gt.0) then
	call stfile(idsevl,nfsevl)
        par_trial=strcompr(file(nfsevl)//trailer_par_seq)
      end if
      xmdtrl_orig=par_xmdtrl
      par_xmdtrl=xmdtrl
      nt_orig=ipar_nt
      ipar_nt=1
      ipar_icmout=1
      ipar_istosc=1
      isetos_orig=ipar_isetos
      ipar_isetos=1
      icase_trailer=icase
      xseq_trailer=par_xmdtrl
      trailer_par=strcompr(settrailer(am, az, icase_trailer,
     *    par_xmdtrl, itcase))
      call mnevol(i_paramset, ierr_param)
      if(istdpr_in.ge.0) call dump_csum('csum1'//trailer_par,
     *  'bcsum1'//trailer_par)
      if(ierr_param.lt.0) then
	write(istder,*) 'Error in call of mnevol'
	return
      end if
c
      ipar_nt=nt_orig
      ipar_icmout=icmout_orig
      ipar_istosc=istosc_orig
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
      if(ierr_param.lt.0) then
	write(istder,*) 'Error in call of redistrb'
	return
      end if
c
c  calculate oscillations
c
      i_paramset=1
      i_inout=0
      call adipls(i_paramset, ierr_param, i_inout,
     *    xr_adi, aar_adi, data_adi, nnr_adi, ivar_adi, iaa_adi)
c
      if(ierr_param.lt.0) then 
	write(istder,*) 'Error in call of adipls'
	return
      end if
      if(istdpr_in.ge.0) call dump_obs('osc/obs'//trailer_par)
      return
      end
