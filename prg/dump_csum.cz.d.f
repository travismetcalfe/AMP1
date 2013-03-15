      subroutine dump_csum(csum_file, bcsum_file)
c
c  outputs complete csum and bcsum files based on
c  evolution sequence stored in common /csum_param/ or common/csum_indiv/,
c  depending on isetos, to files csum_file and bcsum_file
c
      implicit double precision(a-h, o-z)
      character*(*) csum_file, bcsum_file
      character cxcnoc*30
c
      include 'engenr.cz.d.incl'
      dimension cxcnoc(4)
c
      common/engcnt/ xhzlm1, xhzlm2, nxhzer, icnocs, nspec,
     *  iheccs, nspect
      common /csum_param/ icsum, nstep, csum_st(icsum_max, nstep_max)
      common /csum_indiv/ icsum_ind, nstep_ind,
     *  csum_ind(icsum_max, nstep_max)
      common/cevlio/ isetos, iastr
c
c  common defining standard input and output, standard error
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      data cxcnoc /'X14,c, X16,c',
     *             'X12,c, X13,c, X14,c, X16,c',
     *             'X12,c, X13,c, X14,c',
     *             'X12,c, X13,c, X14,c, X16,c'/
c
c  test for no models
c
      if(nstep.le.0) then
	write(istdou,105) 
	if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,105)
	return
      end if
c
      open(99,status='unknown',file=csum_file,form='formatted')
c
      alsun = 3.846d33
c
      if(icnocs.lt.1) then
        write(99,110)
      else
        lxcnoc=length(cxcnoc(icnocs))
        write(99,115) cxcnoc(icnocs)(1:lxcnoc)
      end if
c
      if(isetos.eq.0) then
        do n = 1, nstep
          write(99,120) (csum_st(i,n),i=1,icsum)
        end do
      else
        do n = 1, nstep_ind
          write(99,120) (csum_ind(i,n),i=1,icsum_ind)
        end do
      end if
c
      close(99)
c
      open(99,status='unknown',file=bcsum_file,form='formatted')
c
      write(99,130)
c
      if(isetos.eq.0) then
        do n = 1, nstep
          write(99,140) n, csum_st(1,n), csum_st(2,n)/1.d9, 
     *      csum_st(3,n), csum_st(4,n), csum_st(5,n)/alsun, 
     *      csum_st(9,n), csum_st(16,n)
        end do
      else
        do n = 1, nstep_ind
          write(99,140) n, csum_ind(1,n), csum_ind(2,n)/1.d9, 
     *      csum_ind(3,n), csum_ind(4,n), csum_ind(5,n)/alsun, 
     *      csum_ind(9,n), csum_ind(16,n)
        end do
      end if
      close(99)
c
  105 format(//' ***** Error in s/r dump_csum. No models stored'/)
  110 format('# Summary of evolution calculation.'/'#'/
     *  '# M/Msun, age, R, Teff, L, Dcz, pc, Tc, Xc, X3c,',
     *  ' rhoc, epsilonc, kappac, dradc, dadc, qcc, xcc, Yc',
     *  ' Xs, Zs:'/'#')
  115 format('# Summary of evolution calculation.'/'#'/
     *  '# M/Msun, age, R, Teff, L, Dcz, pc, Tc, Xc, X3c,',
     *  ' rhoc, epsilonc, kappac, dradc, dadc, qcc, xcc, Yc, ',a,
     *  ', Xs, Zs:'/'#')
  120 format(f8.3,1p2e13.5,0pf10.2,1pe13.5,0pf10.5,1p7e13.5,
     *  0p2f10.5,1p8e13.5)
  130 format('# n, M/Msun, age (Gyr), R, Teff, L/Lsun, Xc, qc'/
     *  '#')
  140 format(i5, f9.4, f9.5, 1pe14.5, 0pf8.1, f10.4, 1pe13.5, e12.4)
      end
