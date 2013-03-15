      subroutine extcno(xcno, n)
c
c  Extracts abundances of active CNO elements from cvr(.,n) 
c  and stores them in xcno(i), i = 1, ... icvcno
c
c  Uses storage indices set in common/cstcvr/ by s/r incomp.
c
c  Original version: 5/8/02
c
      implicit double precision (a-h,o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.nnz.d.incl'
      parameter (idermx = ((nspcmx+3)*(nspcmx+4))/2)
c
      dimension xcno(*)
c
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
      common/engcnt/ xhzlm1, xhzlm2, nxhzer, icnocs, nspec,
     *  iheccs, nspect
      common/cengcs/ ixc12, ixc13, ixn14, ixn15, ixo16, ixo17
      common/compos/ xseth, yset, xset3, xset12, xset13, xset14, xset16
      common/compvr/ cvr(icvrmx, 1)
      common/cstcvr/ icvh1, icvhe4, icvhe3, icvc12, icvc13, icvn14, 
     *  icvo16
      common/cstyst/ iyh1, iyhe3, iyc12, iyc13, iyn14, iyo16, 
     *  iyche4, iycc12
      common/cstazt/ iah1, iahe3, iac12, iac13, ian14, iao16, 
     *  iache4, iacc12
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c
c  start section handling general CNO cases
c
      if(icnocs.eq.1.or.icnocs.eq.11) then
c
        xcno(1) = cvr(icvn14,n)
        xcno(2) = cvr(icvo16,n)
c
      else if(icnocs.eq.2.or.icnocs.eq.4) then
c
        xcno(1) = cvr(icvc12,n)
        xcno(2) = cvr(icvc13,n)
        xcno(3) = cvr(icvn14,n)
        xcno(4) = cvr(icvo16,n)
c
      else if(icnocs.eq.3) then
c
        xcno(1) = cvr(icvc12,n)
        xcno(2) = cvr(icvc13,n)
        xcno(3) = cvr(icvn14,n)
c
      else
c
c  diagnostic message 
c
        write(istdou,190) icnocs
        if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdou,190) icnocs
        stop 'extcno'
      end if
c
      return
  190 format(//' ***** Error in s/r extcno. icnocs =',i5,
     *  ' not implemented')
      end
