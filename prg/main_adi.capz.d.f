      program main
c
c  Main programme for calculation of oscillation frequencies (and possibly
c  rotational splittings for a given model on amdl form).
c  This includes the mesh redistribution.
c
c  The file name for the amdl file is read in directly and the
c  file is assumed to contain a single model, to be used at each
c  mesh point. Reading the model is carried out by readml_sim,
c  included in this file.
c
c  Results are produced in obs file, replacing `amdl' by `obs'
c  in the name of the model file or, if amdl_file does not contain
c  `amdl', in the file obs.no_name.
c
c  Original version: 10/2/06
c
c  Modified 19/12/07, increasing size of arrays in common/cofile/
c
      implicit double precision (a-h,o-z)
      include 'engenr.bz.d.incl'
      parameter(iaa_adi=11)
      character trailer_par*80, strcompr*80, str_replace*80, 
     *  file*80, filess*80, 
     *  in_evol*80, in_rdist*80, in_adi*80, obs_file*80, amdl_file*80
      dimension xr_adi(nnmax), aar_adi(iaa_adi,nnmax)
c
      common/trl_param/ trailer_par
      common/cofile/ nfiles, idsfil(99), file(99), iopen(99), filess(99)
      common/nmbmsh/ nn_evol, nn_adi, ivar_adi
      common/xnwvar/ x_adi(1)
      common/anwvar/ data_adi(8), aa_adi(istrmx,1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  hardcode trailer in output file names
c
      trailer_par='.no_name'
c
c  input files for parameter input to redistribution
c  and pulsations
c
      read(istdin,*) in_rdist
      read(istdin,*) in_adi
c
c  read input file for amdl model, assumed to be of form 
c
      read(istdin,*) amdl_file
c
c  general setups of storage and initial parameters
c
      call setups_package('', in_rdist, in_adi, ierr_param)
      if(ierr_param.lt.0) stop 'Error in setups_package'
c
c  read model
c
      call readml_sim(amdl_file, x_adi, aa_adi, data_adi, nn_adi,
     *  ivarmd, istrmx, ierr_readml)
c
c  redistribute mesh in model
c
      i_paramset=1
      i_inout=0
      call srdist(i_paramset, ierr_param, i_inout,
     *  x_adi, aa_adi, data_adi, xr_adi, aar_adi, nn_adi,
     *  nnr_adi, ivarmd, istrmx, iaa_adi)
      if(ierr_param.lt.0) stop 'Error in call of redistrb'
c
c  set rotation rate
c
      icontr=0
      call set_rotation(xr_adi, nnr_adi, icontr)
c
c  calculate oscillations in model
c
      i_paramset=1
      i_inout=0
      call adipls(i_paramset, ierr_param, i_inout,
     *  xr_adi, aar_adi, data_adi, nnr_adi, ivarmd, iaa_adi)

      if(ierr_param.lt.0) stop 'Error in call of adipls'
      obs_file=str_replace(amdl_file,'amdl/','osc/',ierr_str)
      obs_file=str_replace(obs_file,'amdl','obsr',ierr_str)
      if(ierr_str.lt.0) obs_file=strcompr('osc/obsr'//trailer_par)
      call dump_obs(strcompr(obs_file))
c
   90 continue
      stop
      end
      subroutine readml_sim(amdl_file, x_adi, aa_adi, data_adi, nn_adi,
     *  ivarmd, iaa_adi, ierr_readml)
c
c  simplified read of model for adiabatic pulsation.
c
      implicit double precision (a-h, o-z)
      character*(*) amdl_file
      dimension x_adi(1), aa_adi(iaa_adi,1), data_adi(1)
      dimension mlname(4)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      ifind = 1
      xmod = 0
      mdintg = 0
      ids = 99
      idsp = 0
      in = 1
      nprmod = 0
      call addfil(ids, amdl_file)
c
      call readml(ifind, xmod, mdintg, ids, idsp, in, 
     *  data_adi, x_adi, aa_adi, iaa_adi, nn_adi,
     *  nprmod, ivarmd, mlname, nwmod, icry)
c
      close(ids)
c
      if(icry.eq.1) then
	ierr_readml = 0
      else
	ierr_readml = -1
      end if
c
      return
      end
      subroutine mnevol(i_paramset, ierr_param)
c
c  dummy subroutine for this pulsation-only package
c
      return
      end
