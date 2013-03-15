      subroutine tstcon(x,y,nn,nvar,iy,eacon,eamcon,ieacon,init,iter)
c
c  test convergence of evolution calculation, after possible resetting
c  of convective core etc. (which might affect tnrkt convergence measure).c
c
c  For init = 1, just store solution in internal array yp for later
c  comparison.
c
c  For now just consider case without diffusion. 
c
c  Original version: 18/11/05
c
c  Modified 18/7/07 to use routine also to restrict magnitude of change
c  Maximum value hardcoded as dyymax, for now. 
c  Also, test is restricted to
c  y(1) - y(3) (which may be most likely to cause problems).
c
c  Modified 9/2/08, to store yp in common for resetting when
c  mesh is reset with mrefine.
c
      implicit double precision (a-h, o-z)
      logical time0
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.cz.d.incl'
      dimension x(1), y(iy,1), eacon(ieacon,3)
      dimension dy(ivarmx),ymean(ivarmx),dymax(ivarmx),
     *  ndmax(ivarmx), nrsmax(ivarmx),dyrsmx(ivarmx)
c
      common/cmtime/ age, time0, lastmd
      common/c_tstcon/ yp(ivarmx,nnmax)
c
c  common defining standard input and output, standard error
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  maximum allowed change
c
      data dyymax /0.5d0/
c
      save
c
      if(init.eq.1) then
	do i=1,nvar
	  do n=1,nn
	    yp(i,n)=y(i,n)
	  end do
	end do
	return
      end if
c
c  set average and maximum changes
c
      ireset=0
      do i=1,nvar
	dy(i)=0.d0
	dymax(i)=0.d0
	ymean(i)=0.d0
        nrsmax(i)=0
        dyrsmx(i)=0.d0
	do n=1,nn
	  dyy=y(i,n)-yp(i,n)
c 
c  Test for restricting change 
c
	  dyabs=abs(dyy)
          if(i.le.3.and.dyabs.ge.dyymax) then
            ireset=1
            dyy1=sign(dyymax,dyy)
            if(dyabs.gt.abs(dyrsmx(i))) then
              nrsmax(i)=n
              dyrsmx(i)=dyy
            end if
            dyabs=dyymax
            dyy=dyy1
            y(i,n)=yp(i,n)+dyy
          end if
	  dy(i)=dy(i)+dyabs
	  ymean(i)=ymean(i)+abs(y(i,n))
	  if(dyabs.gt.abs(dymax(i))) then
	    ndmax(i)=n
	    dymax(i)=dyy
          end if
        end do
	dy(i)=dy(i)/nn
	ymean(i)=ymean(i)/nn
      end do
c
c  diagnostic output for reset
c
      if(ireset.eq.1) then
        write(istdou,'(/'' ***** Excessive changes reset to'',f10.5/
     *         ''       i, max. original change, n at maximum:''/)')
     *    dyymax
        if(istdpr.gt.0.and.istdpr.ne.istdou) 
     *    write(istdpr,'(/'' ***** Excessive changes reset to'',f10.5/
     *         ''       i, max. original change, n at maximum:''/)')
     *    dyymax
        do i=1,ivarmx
          if(nrsmax(i).gt.0) then
            write(istdou,'(2i5,1pe13.5)') i, nrsmax(i), dyrsmx(i)
            if(istdpr.gt.0.and.istdpr.ne.istdou) 
     *        write(istdpr,'(2i5,1pe13.5)') i, nrsmax(i), dyrsmx(i)
          end if
        end do
      end if
c
c  store in eacon. Note that the first 4 quantities are logarithmic
c  and hence relative differences to not make sense
c
      do i=1,4
	eacon(i,1)=dy(i)
	eacon(i,2)=dymax(i)
	eacon(i,3)=ndmax(i)
      end do
      do i=5,nvar
	eacon(i,1)=dy(i)/ymean(i)
	eacon(i,2)=dymax(i)
	eacon(i,3)=ndmax(i)
      end do
c
c  set average change over log r, log f, log T and log L
c
      eamcon=amean(eacon(1,1),4)
c
c  output
c
      write(istdou,'(/'' Convergence test. Iteration no'', i3,
     *  '' average change eamcon ='',1pe13.5/)') iter, eamcon
      if(istdpr.gt.0) then
        if(istdpr.ne.istdou) 
     *    write(istdpr,'(/'' Convergence test. Iteration no'', i3,
     *    '' average change eamcon ='',1pe13.5/)') iter, eamcon
	write(istdpr,'('' eacon:'')')
        if(time0) then
          npvar = 4
        else
          npvar = nvar
        end if
	do k=1,2
	  write(istdpr,'(1p10e11.3)') (eacon(i,k),i=1,npvar)
        end do
        write(istdpr,'(10i11)') (ndmax(i),i=1,npvar)
        write(istdpr,'(/)') 
      end if
c
c  finally store solution for comparing with next iteration
c
      do i=1,nvar
        do n=1,nn
          yp(i,n)=y(i,n)
        end do
      end do
      return
      end
