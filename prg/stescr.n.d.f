      subroutine stescr
c
c  initialize scratch directory for evolution output of various
c  sort (formerly placed in scrdir as set in include file)
c
c  The scratch root, including the system and the user dependent
c  parts, are set by call of s/r setscr
c
      character*80 scrusr, scrdir
      common/cscrtc/ scrdir, lscr
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      call setscr(scrusr)
      lscr=length(scrusr)
      scrdir=scrusr(1:lscr)//'/evolprg/'
      lscr=length(scrdir)
      write(istdou,100) scrdir(1:lscr)
      return
  100 format(//' Now scratch directory has been set to'/1x,a//)
      end
