      subroutine print_par(header_par)
c
c  Prints parameters to file defined by header_par and trailer_par
c
      implicit double precision (a-h, o-z)
      character*(*) header_par
      character*80 trailer_par, file
      common/cvr_param/ par_am, par_z, par_agefin, par_rsfin, 
     *  par_alsfin, par_zxsfin, par_xmdtrl,
     *  par_xxh, par_fdgopl, par_tlopfg, par_dtlopf, par_tlopf1,
     *  par_fcno, par_agehe3, par_xrz12, par_xrz13, par_xrz14,
     *  par_xrz16, par_tprfct, par_alfa, par_etac, par_phc, par_alpove,
     *  par_alpovc, par_velrot
      common/cvi_param/ ipar_nt, ipar_icsove, ipar_icsovc, ipar_isprot,
     *  ipar_icmout, ipar_istosc, ipar_isetos
      common/trl_param/ trailer_par
c
c  common defining standard input and output, standard error
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  set complete file name.
c
      file=header_par(1:length(header_par))//trailer_par
      iunit=99
      open(iunit,file=file,form='formatted',status='unknown')
      write(iunit,110) par_am, par_z, par_agefin, par_rsfin, par_alsfin,
     *  par_zxsfin, par_xmdtrl, par_xxh, par_fdgopl, 
     *  par_tlopfg, par_dtlopf, par_tlopf1,
     *  par_fcno, par_agehe3, par_xrz12, par_xrz13, par_xrz14,
     *  par_xrz16, par_tprfct, par_alfa, par_etac, par_phc, par_alpove,
     *  par_alpovc, par_velrot, ipar_nt, ipar_icsove, ipar_icsovc, 
     *  ipar_isprot, ipar_icmout, ipar_istosc, ipar_isetos
c
      close(iunit)
      return
  110 format(//
     *   ' Parameters set in common/cvr_param/ and common/cvi_param/:'//
     *   '  par_am     = ',1pe13.5/
     *   '  par_z      = ',1pe13.5/
     *   '  par_agefin = ',1pe13.5/
     *   '  par_rsfin  = ',1pe13.5/
     *   '  par_alsfin = ',1pe13.5/
     *   '  par_zxsfin = ',1pe13.5/
     *   '  par_xmdtrl = ',1pe13.5/
     *   '  par_xxh    = ',1pe13.5/
     *   '  par_fdgopl = ',1pe13.5/
     *   '  par_tlopfg = ',1pe13.5/
     *   '  par_dtlopf = ',1pe13.5/
     *   '  par_tlopf1 = ',1pe13.5/
     *   '  par_fcno   = ',1pe13.5/
     *   '  par_agehe3 = ',1pe13.5/
     *   '  par_xrz12  = ',1pe13.5/
     *   '  par_xrz13  = ',1pe13.5/
     *   '  par_xrz14  = ',1pe13.5/
     *   '  par_xrz16  = ',1pe13.5/
     *   '  par_tprfct = ',1pe13.5/
     *   '  par_alfa   = ',1pe13.5/
     *   '  par_etac   = ',1pe13.5/
     *   '  par_phc    = ',1pe13.5/
     *   '  par_alpove = ',1pe13.5/
     *   '  par_alpovc = ',1pe13.5/
     *   '  par_velrot = ',1pe13.5//
     *   ' ipar_nt     = ',i13/
     *   ' ipar_icsove = ',i13/
     *   ' ipar_icsovc = ',i13/
     *   ' ipar_isprot = ',i13/
     *   ' ipar_icmout = ',i13/
     *   ' ipar_istosc = ',i13/
     *   ' ipar_isetos = ',i13/)
      end
