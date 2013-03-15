      subroutine pushy(x,y,zk,yp,zkp,dt,dtp,eta,ii,kk,nn,id)
c
c  estimates trial solution at time dt after y and zk by linear
c  extrapolation, with correction factor eta, and sets into yp, zkp.
c  eta = 1 gives full extrapolation
c  eta = 0 gives yp = y
c  routine must be supplied with yp, zkp at time dtp prior to y, zk.
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
c  New version with simplified (non-Gough) logics. From 7/11/90
c
c  Modified 9/10/02 to test for negative abundances. This needs
c  to be changed to allow also for diffusion, later.
c
      implicit double precision (a-h,o-z)
      dimension x(1),y(id,1),zk(1),yp(id,1),zkp(1)
      common/clshft/ alshft
      common/heavy/ zatmos, zhc, zh(1)
      common/cstyst/ iyh1, iyhe3, iyc12, iyc13, iyn14, iyo16, 
     *  iyche4, iycc12
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec, idcomp, iccomp
      common/engcnt/ xhzlm1, xhzlm2, nxhzer, icnocs, nspec,
     *  iheccs, nspect
      common/cdiffu/ idiffus, itbdif, rhob, rmxdif, tlfdif, dtldif,
     *  ctdfmx, rctdif, rctdf1,idfxbc, ismdif, nsmdif
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
      if(istdpr.gt.0) write(istdpr,'(/'' Entering pushy'')')
      idiag=2
      beta=eta*dt/dtp
      alpha=1.d0+beta
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'In pushy, alpha, beta =',alpha,beta
      ipr=4
      if(idiag.ge.2.and.istdpr.gt.0) 
     *  write(istdpr,*) ' Results for i =', ipr
c
      do 20 n=1,nn
      do 20 i=1,ii
      ynew=alpha*y(i,n)-beta*yp(i,n)
      if(idiag.ge.2.and.(mod(n,5).eq.0.or.n.ge.nn-5).and.i.eq.ipr
     *  .and.istdpr.gt.0) 
     *    write(istdpr,'(i5,1p4e15.7)') n, x(n), yp(i,n),y(i,n),ynew
   20 yp(i,n)= alpha*y(i,n)-beta*yp(i,n)
c
c  test for resetting lumininosity with alshft (where extrapolation
c  is not formally correct)
c
      if(alshft.gt.0) then
	alshfl=log10(alshft)
	do n=1,nn
	  if(y(4,n).le.alshfl+0.1) then
	    yp(4,n)=y(4,n)
	  end if
	end do
      end if
c
c  reset composition to avoid negative values (needs fixing with 
c  diffusion)
c
      if(idiffus.eq.0) then
        do 25 n=1,nn
        if(yp(5,n).lt.xhzlm1) then
          if(iheccs.ne.0) then
            yp(iyche4,n)=yp(iyche4,n)+yp(5,n)-1.d-10
          end if
          yp(5,n)=1.d-10
        end if
   25   continue
      end if
c
      if(kk.eq.0) return
      do 30 k=1,kk
   30 zkp(k)=alpha*zk(k)-beta*zkp(k)
      return
      end
