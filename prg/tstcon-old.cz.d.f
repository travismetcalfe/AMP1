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
      implicit double precision (a-h, o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.cz.d.incl'
      dimension x(1), y(iy,1), eacon(ieacon,3)
      dimension yp(ivarmx,nnmax),dy(ivarmx),ymean(ivarmx),dymax(ivarmx),
     *  ndmax(ivarmx)
c
c  common defining standard input and output, standard error
c
      common/cstdio/ istdin, istdou, istdpr, istder
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
      do i=1,nvar
	dy(i)=0.d0
	dymax(i)=0.d0
	ymean(i)=0.d0
	do n=1,nn
	  dyy=y(i,n)-yp(i,n)
	  dyabs=abs(dyy)
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
     *  '' average change eamcon ='',1pe13.5/)') iter, eamcom
      if(istdpr.gt.0) then
        if(istdpr.ne.istdou) 
     *    write(istdpr,'(/'' Convergence test. Iteration no'', i3,
     *    '' average change eamcon ='',1pe13.5/)') iter, eamcom
	write(istdpr,'('' eacon:'')')
	do k=1,2
	  write(istdpr,'(1p10e11.3)') (eacon(i,k),i=1,nvar)
        end do
        write(istdpr,'(10i11)') (ndmax(i),i=1,nvar)
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
