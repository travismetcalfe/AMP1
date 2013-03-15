      subroutine rsksdr(bcp,bc,nspcp,nspc)
c
c  Copies contents of common/ksider/ stored in bcp with
c  nspc = nspcp into bc, with nspc = nspc.
c  If nspcp .le. 0, assume instead that bcp is stored as in
c  old (pre-July 1991) version of programme.
c
c  Original version: 10/8/91.
c
      implicit double precision(a-h, o-z)
      logical diag
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
c  Note: engenr.n.d.incl replaced by engenr.nnz.d.incl, 5/6/02
      include 'engenr.nnz.d.incl'
c
      dimension bcp(1), bc(1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  for the moment, hardcode idiag
c
      idiag=0
c
      diag = idiag.ge.1.and.istdpr.gt.0
c
c  It is not possible to store larger array into smaller
c 
      if(nspcp.gt.nspc) then
	write(istdou,105) nspcp, nspc
	if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,105) 
     *    nspcp, nspc
	stop 'rsksdr'
      end if
c
      if(diag) write(istdpr,107) nspcp, nspc
c
c  set total number of elements to be set in bc
c
      nstore = 12 + 5*idthm + nspc*(4+nspc)
c
c  to initialize, zero bc
c
      call zero(bc,nstore)
c
      nazt = nspc+3
      nb = nspc+2
c
      nxxh = nspc*nb
      if(nspcp.gt.0) then
        naztp = nspcp+3
        nbp = nspcp+2
	nxxhp = nspcp*nbp
      else
	naztp = 5
	nbp = 5
	nxxhp = 15
      end if
c
      ist = 1
      bc(ist) = bcp(ist)
      if(diag) write(istdpr,110) 'al0 ', ist, ist
      ist = 2
      bc(ist) = bcp(ist)
      if(diag) write(istdpr,110) 'al2 ', ist, ist
      ist = ist + 1
      istp = ist 
      nvar = nazt
      nvarp = naztp
      call store(bcp(istp),bc(ist),nvarp)
      if(diag) write(istdpr,120) 'azt ',istp,istp+nvarp-1,
     *  ist,ist+nvar-1
      ist = ist + nvar
      istp = istp + nvarp
      call store(bcp(istp),bc(ist),nvarp)
      if(diag) write(istdpr,120) 'ax  ',istp,istp+nvarp-1,
     *  ist,ist+nvar-1
      ist = ist + nvar
      istp = istp + nvarp
      bc(ist)=bcp(istp)
      if(diag) write(istdpr,110) 'pt  ', istp, ist
      ist = ist + 1
      istp = istp + 1
      bc(ist)=bcp(istp)
      if(diag) write(istdpr,110) 'tt  ', istp, ist
      ist = ist + 1
      istp = istp + 1
      nvar = idthm
      nvarp = idthm
      call store(bcp(istp),bc(ist),nvarp)
      if(diag) write(istdpr,120) 'th  ',istp,istp+nvarp-1,
     *  ist,ist+nvar-1
      ist = ist + nvar
      istp = istp + nvarp
      call store(bcp(istp),bc(ist),nvarp)
      if(diag) write(istdpr,120) 'xt  ',istp,istp+nvarp-1,
     *  ist,ist+nvar-1
      ist = ist + nvar
      istp = istp + nvarp
      nvar = nxxh
      nvarp = nxxhp
      call store(bcp(istp),bc(ist),nvarp)
      if(diag) write(istdpr,120) 'xxh ',istp,istp+nvarp-1,
     *  ist,ist+nvar-1
      ist = ist + nvar
      istp = istp + nvarp
      bc(ist)=bcp(istp)
      if(diag) write(istdpr,110) 'alt ', istp, ist
      ist = ist + 1
      istp = istp + 1
      bc(ist)=bcp(istp)
      if(diag) write(istdpr,110) 'altt', istp, ist
      ist = ist + 1
      istp = istp + 1
      nvar = idthm
      nvarp = idthm
      call store(bcp(istp),bc(ist),nvarp)
      if(diag) write(istdpr,120) 'alh ',istp,istp+nvarp-1,
     *  ist,ist+nvar-1
      ist = ist + nvar
      istp = istp + nvarp
      call store(bcp(istp),bc(ist),nvarp)
      if(diag) write(istdpr,120) 'alhh',istp,istp+nvarp-1,
     *  ist,ist+nvar-1
      ist = ist + nvar
      istp = istp + nvarp
      call store(bcp(istp),bc(ist),nvarp)
      if(diag) write(istdpr,120) 'albb',istp,istp+nvarp-1,
     *  ist,ist+nvar-1
      return
  105 format(//' ***** Error in s/r rsksdr. nspcp =',i4,
     *  ' .gt. nspc =',i4)
  107 format(//' Restore variables in ksider. Old nspcmx =',i4,
     *  '  New nspcmx =',i4//
     *  ' Variable   Old storage    New storage'/)
  110 format(3x,a4,5x,i4,12x,i4)
  120 format(3x,a4,3x,i4,' - ',i4,5x,i4,' - ',i4)
      end
