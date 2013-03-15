      subroutine setups_package(in_evol, in_rdist, in_adi, ierr_param)
c
c  setups and initial parameter input for full combined package
c  in_evol, in_rdist and in_adi must be file names for the input files
c  for the evolution, mesh distribution and oscillation packages,
c  assumed already to be stripped of comment lines.
c  If one of these file names is given as '' (an empty string)
c  the corresponding input (and possibly setup) is skipped.
c
      implicit double precision (a-h,o-z)
      character *(*) in_evol, in_rdist, in_adi
      character*80 header_par, trailer_par_seq
      dimension x_adi(1), aa_adi(5,1), data_adi(1), 
     *  xr_adi(1), aar_adi(5,1)
      common/par_control/ itcase, icase_trailer, xseq_trailer, 
     *  header_par, trailer_par_seq
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      write(istdou,'(/'' Entering s/r setups_package''/)')
      if(istdpr.gt.0.and.istdpr.ne.istdou)
     *  write(istdpr,'(/'' Entering s/r setups_package''/)')
c
c  initialize header_par (often later overwritten)
c
      header_par='ttt.par'
c
      itcase=0
c
      il_evol= length(in_evol)
      il_rdist=length(in_rdist)
      il_adi=  length(in_adi)
c
      if(il_evol.gt.0) call setups_main
      if(il_adi.gt.0)  call setups_adi
c
      i_paramset = 1
      istdin_orig=istdin
c
c  Initial calls to read parameter input files
c
      if(il_evol.gt.0) then
c
c  Evolution parameters
c
        istdin=21
        open(istdin,file=in_evol,status='old')
        call mnevol(i_paramset, ierr_param)
        close(istdin)
        if(ierr_param.lt.0) then
          write(istder,'(/'' Error in call of mnevol''/)')
	  istdin=istdin_orig
          return
        end if
c
      end if
c
      if(il_rdist.gt.0) then
c
c  Redistribution parameters
c
        istdin=21
        open(istdin,file=in_rdist,status='old')
        ivarmd=5
        i_inout=0
        call srdist(i_paramset, ierr_param, i_inout,
     *    x_adi, aa_adi, data_adi, xr_adi, aar_adi, nn_adi,
     *    nnr_adi, ivarmd, istrmx, iaa_adi)
        close(istdin)
        if(ierr_param.lt.0) then
          write(istder,'(/'' Error in call of srdist''/)')
	  istdin=istdin_orig
          return
        end if
c
      end if
c
      if(il_adi.gt.0) then
c
c  Pulsation parameters
c
        istdin=21
        open(istdin,file=in_adi,status='old')
        call adipls(i_paramset, ierr_param, i_inout,
     *    xr_adi, aar_adi, data_adi, nnr_adi, ivarmd, iaa_adi)
        close(istdin)
        if(ierr_param.lt.0) then
          write(istder,'(/'' Error in call of adipls''/)')
	  istdin=istdin_orig
          return
        end if
c
      end if
c
      write(istdou,'(/'' Exiting s/r setups_package''/)')
      if(istdpr.gt.0.and.istdpr.ne.istdou)
     *  write(istdpr,'(/'' Exiting s/r setups_package''/)')
c
      istdin=istdin_orig
      return
      end
