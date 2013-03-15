      subroutine setcno(y, n)
c
c  Sets CNO abundances in xset12, ... and, if n gt 0, in cvr(4 - ..., n)
c  based on composition data from y, the standard set of
c  dependent variables (for case without diffusion, as used for
c  model output).
c
c  Uses storage indices set in common/csycst/ by s/r incomp.
c
c  Original version: 15/5/2000
c
      implicit double precision (a-h,o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.nnz.d.incl'
      parameter (idermx = ((nspcmx+3)*(nspcmx+4))/2)
c
      dimension y(1)
c
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
      common/consts/ av,ah,ahe
      common/rcncns/ a(knucst),b(knucst),s(knucst,isnuc),
     *  zz(knucst),zz12(2,knucst),zsmh,acno,awght(knucwg),q(knucq), 
     *  knucpm, kqvals, kawght
      common/engcnt/ xhzlm1, xhzlm2, nxhzer, icnocs, nspec,
     *  iheccs, nspect
      common/he3fdg/ agesh,ifdhe3,iche30,iche31
      common/cengcs/ ixc12, ixc13, ixn14, ixn15, ixo16, ixo17
      common/compos/ xseth, yset, xset3, xset12, xset13, xset14, xset16
      common/cstcvr/ icvh1, icvhe4, icvhe3, icvc12, icvc13, icvn14,
     *  icvo16
      common/cstyst/ iyh1, iyhe3, iyc12, iyc13, iyn14, iyo16, iyche4,
     *  iycc12
      common/compvr/ cvr(icvrmx, 1)
      common/cnofrc/ fcno, xtlcno
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
      xh = y(5)
c
c  storage contribution from 3He
c
      if(iche30.eq.1) then
	ishxx3=1
      else
	ishxx3=0
      end if
c
      if(icnocs.eq.0) then
c
c  in this case, set just the assumed N14 abundance in common/compos/
c
        xset14 = fcno*z/acno
c
        return
c
      else if(icnocs.ne.1.and.icnocs.ne.4) then
c
c  diagnostic message 
c
        write(istdou,190) icnocs
        if(istdpr.ne.istdou.and.istdpr.gt.0) write(istdpr,190) icnocs
        stop 'setcno'
      end if
c
c  start section handling general CNO cases
c
      if(icnocs.eq.1) then
        xtln14 =  y(iyn14)/awght(4)
        xtlo16 = xtlcno - xtln14
c
c  set N14 and O16 abundances in common/compos/
c
        xset14 = y(iyn14)
        xset16 = xtlo16*awght(6)
c 
	if(n.gt.0) then
	  cvr(icvn14,n) = xset14
	  cvr(icvo16,n) = xset16
        end if
c
c  end for icnocs = 1
c
      else if(icnocs.eq.4) then
c
        xtlc12 =  y(iyc12)/awght(2)
        xtlc13 =  y(iyc13)/awght(3)
        xtln14 =  y(iyn14)/awght(4)
        xtlo16 = xtlcno - xtlc12 - xtlc13- xtln14
c
c  set CNO abundances in common/compos/
c
        xset12 = y(iyc12)
        xset13 = y(iyc13)
        xset14 = y(iyn14)
        xset16 = xtlo16*awght(6)
c
	if(n.gt.0) then
	  cvr(icvc12,n) = xset12
	  cvr(icvc13,n) = xset13
	  cvr(icvn14,n) = xset14
	  cvr(icvo16,n) = xset16
        end if
c
      end if
c
      return
  190 format(//' ***** Error in s/r setcno. icnocs =',i5,
     *  ' not inplemented')
  910 format(/a60/(1p7e11.3))
      end
