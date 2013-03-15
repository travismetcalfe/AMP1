      subroutine cpydif(y1,y2,ii1,ii2,ii3,nn,idiffus,iy,icase,init)
c
c  copies dependent variables between output form and diffusion
c  solution form, from y1 to y2
c  If icase = 1, y1 is in output form and y2 in diffusion form
c  If icase = 2, y1 is in diffusion form and y2 in output form
c
c  Also note that Y-variables (diffusion velocity) are not
c  set in this routine. This must be done separately
c
c  Original version 16/8/92
c
c  Modified 5/8/95, to include general number of diffusing species.
c
c  Modified 31/10/00, to use routine also for transformation of 
c  luminosity, including a shift in the luminosity variable used
c  for computing when the central hydrogen abundance is very small
c  (hard-coded to .le. 1.e-4). This is reset only when init = 1,
c  which is assumed to require icase = 1.
c
      implicit double precision(a-h, o-z)
      dimension y1(iy,1), y2(iy,1)
c
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec, idcomp, iccomp
      common/clshft/ alshft
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      data xhclim /1.d-4/
      data alsfct /1.d-4/
c
      icomps=idcomp+iccomp
      idcom2=idcomp+idcomp
c
c  try switching off alshft setting here.
c..c
c..c  test for setting luminosity shift
c..c
c..      if(init.eq.1.and.icase.eq.1) then
c..	write(6,*) 'In cpydif,', y1(5,nn),xhclim,alsfct,y1(4,1)
c..	if(y1(5,nn).lt.xhclim) then
c..	  alshft=alsfct*10.d0**y1(4,1)
c..	else
c..	  alshft=0.d0
c..        end if
c..	if(istdpr.gt.0) write(istdpr,110) alshft
c..      end if
c
c  test for case
c
      if(icase.eq.1) then
	do 20 n=1,nn
	do 15 i=1,3
   15   y2(i,n)=y1(i,n)
	y2(4,n)=y1(5,n)
	if(idcomp.gt.1) then
	  do 16 j=2,idcomp
   16     y2(3+j,n)=y1(4+iccomp+j,n)
	end if
	y2(4+idcomp,n)=y1(4,n)
c..	if(alshft.eq.0) then
c..	  y2(4+idcomp,n)=y1(4,n)
c..        else
c..	  y2(4+idcomp,n)=log10(10.d0**y1(4,n)+alshft)
c..	end if
	do 17 j=1,idcomp
   17   y2(4+idcomp+j,n)=y1(4+icomps+j,n)
	if(iccomp.gt.0) then
	  do 18 j=1,iccomp
   18     y2(4+idcom2+j,n)=y1(5+j,n)
	end if
   20   continue
c
      else
	do 30 n=1,nn
	do 25 i=1,3
   25   y2(i,n)=y1(i,n)
	y2(4,n)=y1(4+idcomp,n)
c..	if(alshft.eq.0) then
c..	  y2(4,n)=y1(4+idcomp,n)
c..        else
c..	  y2(4,n)=log10(10.d0**y1(4+idcomp,n)-alshft)
c..	end if
	y2(5,n)=y1(4,n)
	if(iccomp.gt.0) then
	  do 27 j=1,iccomp
   27     y2(5+j,n)=y1(4+idcom2+j,n)
	end if
	if(idcomp.gt.1) then
	  do 28 j=2,idcomp
   28     y2(4+iccomp+j,n)=y1(3+j,n)
	end if
	do 29 j=1,idcomp
   29   y2(4+icomps+j,n)=y1(4+idcomp+j,n)
   30   continue
c
	xhc=y2(5,nn)
c
      end if
      return
  110 format(//' In cpydif, alshft reset to ',1pe13.5/)
      end
