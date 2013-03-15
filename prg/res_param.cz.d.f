      subroutine res_param(i_param,
     *  am, z, agefin, rsfin, alsfin, zxsfin, xmdtrl, xxh, fdgopl,
     *  tlopfg, dtlopf, tlopf1, fcno, agehe3, xrz12, xrz13, xrz14,
     *  xrz16, tprfct, alfa, etac, phc, alpove, alpovc, velrot,
     *  nt, icsove, icsovc, isprot, icmout, istosc, isetos)
c
c  set or reset parameters to be iterated in calls of mnevol.
c
c  For i_param = 1, store parameters in par_am, par_z, etc. from 
c  parameters in arguments. Also, stores trial file name in par_trial
c  in common/trial_param/
c
c  For i_param = 2, reset parameters in arguments from par_am, ...
c  Also, reopen file for trial-model input with filename par_trial
c  set in common/trial_param/
c
c  Original version: 11/10/04
c
c  Modified 19/12/07, increasing size of arrays in common/cofile/
c
      implicit double precision(a-h, o-z)
      character par_trial*80, file*80
      common/cvr_param/ par_am, par_z, par_agefin, par_rsfin, 
     *  par_alsfin, par_zxsfin, par_xmdtrl,
     *  par_xxh, par_fdgopl, par_tlopfg, par_dtlopf, par_tlopf1,
     *  par_fcno, par_agehe3, par_xrz12, par_xrz13, par_xrz14,
     *  par_xrz16, par_tprfct, par_alfa, par_etac, par_phc, par_alpove,
     *  par_alpovc, par_velrot
      common/cvi_param/ ipar_nt, ipar_icsove, ipar_icsovc, ipar_isprot,
     *  ipar_icmout, ipar_istosc, ipar_isetos
      common/trial_param/ par_trial
      common/cofile/ nfiles, idsfil(99), file(99), iopen(99)
c
c  dataset designators
c
      common/cdsdsg/ idstrl,idshvz,idszab,idsevl,idssum,idscen,idsnov,
     *  idsoov,idsgng,idsgsm,idstm1,idstm2,iducen
      common/cevlio/ iseto1, iastr
c
c  common defining standard input and output, standard error
c
      common/cstdio/ istdin, istdou, istdpr, istder
c 
      save
c
      if(istdpr.gt.0) then
        write(istdpr,*) 'Enter res_param with i_param =',i_param
        write(istdpr,*) 'xxh, par_xxh =',xxh, par_xxh
      end if
c
c  test for storing original parameters in common /cvr_param/
c
      if(i_param.eq.1) then
        par_am = am
        par_z = z
        par_agefin = agefin
        par_rsfin = rsfin
        par_alsfin = alsfin
        par_zxsfin = zxsfin
        par_xmdtrl = xmdtrl
        par_xxh = xxh
        par_fdgopl = fdgopl
        par_tlopfg = tlopfg
        par_dtlopf = dtlopf
        par_tlopf1 = tlopf1
        par_fcno = fcno
        par_agehe3 = agehe3
        par_xrz12 = xrz12
        par_xrz13 = xrz13
        par_xrz14 = xrz14
        par_xrz16 = xrz16
        par_tprfct = tprfct
        par_alfa = alfa
        par_etac = etac
        par_phc = phc
        par_alpove = alpove
        par_alpovc = alpovc
        par_velrot = velrot
c
        ipar_nt = nt
        ipar_icsove = icsove
        ipar_icsovc = icsovc
        ipar_isprot = isprot
        ipar_icmout = icmout
        ipar_istosc = istosc
        ipar_isetos = isetos
c
c  store initial trial model file in par_trial
c
        if(istdpr.gt.0) write(istdpr,*) 'Before call of stfile'
        call stfile(idstrl,nfstrl)
        if(nfstrl.gt.0) par_trial = file(nfstrl)
c
        if(istdpr.gt.0) write(istdpr,*) 'par_trial set to',par_trial
        write(istdou,'(/'' par_trial set to '',a50/)') par_trial
c
c  test for setting parameters from /cvr_param/ into internal 
c  parameters
c
      else if(i_param.eq.2) then
        am = par_am
        z = par_z
        agefin = par_agefin
        rsfin = par_rsfin
        alsfin = par_alsfin
        zxsfin = par_zxsfin
        xmdtrl = par_xmdtrl
        xxh = par_xxh
        fdgopl = par_fdgopl
        tlopfg = par_tlopfg
        dtlopf = par_dtlopf
        tlopf1 = par_tlopf1
        fcno = par_fcno
        agehe3 = par_agehe3
        xrz12 = par_xrz12
        xrz13 = par_xrz13
        xrz14 = par_xrz14
        xrz16 = par_xrz16
        tprfct = par_tprfct
        alfa = par_alfa
        etac = par_etac
        phc = par_phc
        alpove = par_alpove
        alpovc = par_alpovc
        velrot = par_velrot
c
        nt = ipar_nt
        icsove = ipar_icsove
        icsovc = ipar_icsovc
        isprot = ipar_isprot
        icmout = ipar_icmout
        istosc = ipar_istosc
        isetos = ipar_isetos
c
c  set new trial input model
c  This may well need elaboration later
c
	if(iastr.gt.0.or.isetos.ne.1) then
          close(idstrl)
          write(istdou,'(/'' Open '',a50,'' as trial''/)') par_trial
          if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,*) 
     *      'Open ',par_trial,' as trial'
          open(idstrl,file=par_trial,form='unformatted',status='old')
        end if
      end if
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'On exit xxh, par_xxh =',xxh, par_xxh
      return 
      end
