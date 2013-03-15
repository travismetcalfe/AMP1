

      subroutine reducemrq (ixx, iyy, isig, iia, isub, 
     +     ma, nobs, chisq, retx)

      implicit none
      double precision iyy(10000), isig(10000)
      double precision ixx(36), retx(7), f0
      integer iia(36), isub(10000)

c     now define for the rest of this program... reduce sizes
      integer ma, nobs
      double precision yi(nobs), si(nobs), xx(ma)
      integer ia(ma), subs(nobs)


c     internal parameters
      integer i, mfit, ndata, niter

c     other internal arrays etc
      double precision chisq, covar(ma, ma), alpha (ma ,ma), alambda

c
c     rescale arrays
c

      do i = 1, ma 
         xx(i) = ixx(i)
         ia(i) = iia(i)
      enddo

c
c    expand data and errors and subs into large arrays of dimension nobs
c      
      do i = 1, nobs
         subs(i) = 0
         yi(i) = 0.d0
         si(i) = 0.d0
      enddo
      do i = 1, nobs
          subs(isub(i)) = 1
          yi(isub(i)) = iyy(i)
          si(isub(i)) = isig(i)
      enddo
 
c
c
c     get dimensions that are used
c

      mfit = 0 
      do i = 1, ma 
         if (ia(i).ne.0) mfit = mfit+1
      end do


      ndata = 0 
      do i = 1, nobs 
         if (subs(i).ne.0) ndata = ndata+1
      end do

 1    format (a25)   
 2    format (i1,x,f9.7,x,f9.7,x,f9.7,x,f7.5,x,d16.10,x)  
      open (16, FILE='svdlm.log', STATUS='UNKNOWN', ACCESS='append')
      write (16, 2) 0,xx(1),xx(2),xx(3),xx(4),xx(5)
      close (16)

c
c     initialize for iteration
c     
      alambda = -1.
      chisq = 1000.
      do niter = 1, 4
         print *,"   "
         print *," iteration ", niter
         print *," "
         print *,"****************************************"
         print *,"parameters (Mass, Z, Y, alpha, Age) & Chi^2" 
         print *,xx(1),xx(2),xx(3),xx(4),xx(5),' ',chisq
         print *,"****************************************"
         print *,"   "
         if (chisq.ge.0.001) then 
c     xx = dim ma, yi = dim nobs, si = dim nobs, subs = dim nobs
c     ndata = actual number of data we have should be equal to where
c     sub eq 1.
c     ia = dim ma 

            call mrqmin (xx, yi, si, ia, subs, ma, nobs, mfit, ndata,
     +           covar, alpha, chisq, alambda)
         endif
      enddo
c
c     write some things to terminal
c
      print *,"   "
      print *,"****************************************"
      print *,'Fit parameters are ', xx
      print *,'with chisq', chisq
      print *,"****************************************"
      print *,"   "


      open (16, FILE='svdlm.log', STATUS='UNKNOWN', ACCESS='append')
 3    format (i1,x,f9.7,x,f9.7,x,f9.7,x,f7.5,x,d16.10,x,f9.7)  
      write (16, 3) 1,xx(1),xx(2),xx(3),xx(4),xx(5),chisq

      close (16)
      do i = 1, 7 
        retx(i) = xx(i)
      enddo   


         
      return
      end

c
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c
c
      subroutine mrqmin (x, y, sig, ia, sub, ma, nobs, mfit, ndata,  
     +       covar, alpha, chisq, alambda)

c
c     This subroutine calculates a few things to decide whether to keep
c     the latest set of parameters.  It updates them if chisq is lower
c     with the newer set of parameters.  
c     the important outputs are alpha, chisq, alambda
c     It calls the routine mrqcof which actually evaluates the func and der.


c     parameters variables

      implicit none
      integer ma, ndata, ia(ma),  nobs, mfit
      double precision alambda, x(ma)
      double precision chisq, alpha(ma,ma)
      double precision y(nobs), sig(nobs),  covar(ma,ma)
c     ones that will be saved
      double precision da(ma), beta(ma)
      double precision atry(ma)
      integer j,k,l, sub(nobs),i
      double precision ochisq
c     these are for the SVD fit for covar matrix
      double precision cvr(mfit, mfit), cv(mfit, mfit), dar(mfit)
      double precision odo(mfit), inv(mfit,mfit), cc(ma)
c     these are to save the 'save' items
      double precision sda(20), sbeta(20), satry(20), f0
      save ochisq, sda, sbeta, satry


c
c     i want to keep atry, da, beta mfit.  im going to keep
c
      do i = 1, ma
         da(i) = sda(i)
         beta(i) = sbeta(i)
         atry(i) = satry(i)
      enddo
      
c     initialization 
      if (alambda.lt.0.0) then 
         alambda = 0.001
         call mrqcof(x , y, sig, ia, sub, ma, nobs,  alpha, beta,
     +         chisq)
        chisq =chisq/ndata
        ochisq = chisq
         do j = 1, ma
            atry (j) = x(j)
         end do
      endif
 

      do j = 1, mfit
         do k = 1, mfit
            covar(j,k) = alpha(j,k) 
         end do
         covar(j,j) = alpha(j,j)* (1.+alambda)
         da(j) = beta(j)
      end do
c
c    1.   here i need to reduce all of the matrices for the inversion
c     the inversion uses SVD
c
c     2.  SVD:  cvr gets replaced by the u so dont need it anymore
c     odo contains the solution: in are mfit, cvr, dar
c     return: cv, odo, inv
c     inv = the inverse; 
c     
      
      do i = 1, mfit
         do j = 1, mfit
            cvr(i,j) = covar(i,j)
         enddo
         dar(i) = da(i)
      enddo

      call singvd (mfit, mfit, cvr, dar, cv, odo, inv, 1)

      do i = 1, mfit
         do j = 1, mfit
            covar(i,j) = inv(i,j)
         enddo
         da(i) = odo(i)
      enddo


      if (alambda.eq.0.0) then 
         call covsrt(covar, nobs, ma, ia, mfit)
         call covsrt(alpha, nobs, ma, ia, mfit)
         return
      endif

c
c     test new set (ATRY); here I have a constant multiplying for the 
c     age derivative
c     this is beause I multiplied the dyda * cc to avoid singular matrices
c     so now to need take this constant back into account
c
      cc(1) = 1.0d0
      cc(2) = 1.0d0
      cc(3) = 1.0d0
      cc(4) = 1.0d0
      cc(5) = 1.0d9
      cc(6) = 1.0d0
      cc(7)= 1.0d0
      
      j = 0
      do l = 1, ma
         if (ia(l).ne.0) then 
            j = j + 1
            atry(l) = x(l) + odo(j)*cc(l)
         endif
      enddo

      call mrqcof(atry, y, sig, ia, sub, ma, nobs, 
     +     covar, da, chisq)
      chisq =chisq/ndata

c
c     save and rerrange some arrays for next iteration
c 



      if (chisq.lt.ochisq) then
         alambda = 0.1*alambda
         ochisq = chisq
         do j=1, mfit
            do k = 1, mfit
               alpha(j,k) = covar(j,k)
            enddo
            beta(j) = da(j)
         enddo
         do l = 1, ma
            x(l) = atry(l)
         enddo
      else
         alambda = 10.*alambda
         chisq = ochisq
      endif

      do i = 1, ma
         sda(i) = da(i)
         sbeta(i) = beta(i)
         satry(i) = atry(i)
      enddo

      return
  
      end

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mrqcof(x , y, sig, ia, sub, ma,  nobs, 
     +     alpha, beta, chisq)

c
c     this subroutine calculates the matrices alpha and beta using information
c     from the models (namely derivatives).  Important RETURN items are
c     alpha, beta, chisq
c
c     INPUT: x(ma) = parameters, y = data, sig = errors, ia= [0,1,0..] telling
c     which parameters to use, sub = [0,1,...]telling which observables to use.
c

      implicit none
      integer ma,  nobs, mfit, nonseis, nfreq
      double precision chisq,  alpha(ma, ma), beta(ma), ymodb(nobs)
      double precision y(nobs), dyda(nobs, ma), sig(nobs), ymod(nobs)
      integer i, j, k, l, m, sub(nobs), ia(ma), ndata, isb(nobs)
      double precision dyi, sig2i, wti, x(ma), rx(ma), dp(ma), dx
      double precision  ymoda(nobs), dyr(38, 5), xxx(38, 5)
      double precision v(5,5), w(5), cc(ma), ct, f0, chisq1, chisq2
      double precision sysi, sys(nobs), syserror(nobs), error

cccc  initialize some arrays

      mfit = 0 
      do i = 1, ma 
         if (ia(i).ne.0) mfit = mfit+1
      end do
      ndata = 0 
      do i = 1, nobs
         if (sub(i).ne.0) ndata = ndata+1
      end do

c     save x for derivatives
      do i = 1, ma
         rx(i) = x(i)
      enddo

      do i = 1, mfit
         do j = 1, i 
            alpha (i, j ) = 0.
         end do
         beta (i) = 0.
      end do

c
c     EFUNC  calculates the model observables at X and the derivatives
c     INPUT: x(ma) = parameters with ia(ma) = [0,1...] which to use
c            sub(nobs) = [0,1,0 ...]which observables to use.
c     RETURN: ymod (nobs): the returned model observables at X
c             dyda(nobs, ma): the derivatives at x

      call efunc(ma, ia, x, nobs, sub, ymod, dyda)
      

      
c
c     NEED TO READ IN THE SURFACE CORRECTION and expand it to SYS
c     if iwflag is 0 then SYSERROR will be = 0.0, 
c      ......      1 then SYSERROR will be = values from "SURFACE"
c
      call surfacesigi(nobs, syserror)

      k = 1
      do i = 1, nobs
         if (sub(i).ne.0) then 
            sys(i) = syserror(k)
            k = k + 1
         endif
      enddo

c
c  calculate the alpha (j, k) & beta (j) & chisq component
c
      chisq1 = 0.d0
      chisq2 = 0.d0
      nonseis = 0
      nfreq = 0
      do i = 1, nobs
         if (sub(i).ne.0) then 
c     ONLY PLACE i NEED To use the systematic error
            error = sqrt(sig(i)*sig(i)+0.25*sys(i)*sys(i))
            sig2i = 1.d0/error**2.
            dyi = y(i) - ymod(i)
            j = 0
            do l = 1, ma
               if (ia(l).ne.0) then
                  j= j+1
                  wti = dyda(i,l)*sig2i
                  k = 0
                  do m = 1 , l
                     if (ia(m).ne.0) then
                        k = k+1
                        alpha(j,k) = alpha(j,k) + wti*dyda(i,m)
                     endif
                  end do
                  beta(j) = beta(j) + dyi*wti
               endif
            end do
            if (i.ge.50) then 
               chisq1 = chisq1 + dyi*dyi*sig2i            
               nfreq = nfreq + 1
            endif
            if (i.lt.50) then 
               chisq2 = chisq2 + dyi*dyi*sig2i
               nonseis = nonseis + 1
            endif
         endif
      end do
      
      chisq2 = chisq2/float(nonseis)
      chisq1 = chisq1/float(nfreq)
      chisq = 0.5*(chisq1 + chisq2)

c
c fill in symmetric side
c
      do j = 2, mfit
         do k = 1, j-1
            alpha(k,j) = alpha(j,k)
         end do
      end do


      return
      end



cccc  this is the programme to do COVSRT, it just expands the matrix
cccc  back into the original dimensions of the problem.
cccc  

      subroutine covsrt (covar, npc, ma, ia, mfit)

      
      implicit none
      integer i, j, k, ma, mfit, npc
      double precision covar(npc, npc), swap
      integer ia(ma)

       do i = mfit + 1, ma
         do j = 1, i 
            covar(i,j) = 0
            covar(j,i) = 0
         enddo
      enddo
      k = mfit
      do j = ma, 1, -1 
         if (ia(j).ne.0) then
            do i = 1, ma
               swap = covar(i,k)
               covar(i,k) = covar(i,j)
               covar(i,j) = swap
            enddo
            do i = 1, ma
               swap = covar(k,i)
               covar(k,i) = covar(j,i)
               covar(j,i) = swap
            enddo
            k = k-1
         endif
      enddo

      return 
      end

