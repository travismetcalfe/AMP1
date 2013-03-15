
      program gridmin

c     x = initial guess
c     yy = observations
c     sig = errors
c     ia = 1 = fit, 0 = no fit


c
c     This runs several minimizations starting at different Y and Z values
c     It loops through Z first by testing different values, saves the 
c     best and then goes through Y values and does the same.
c     After each minimization, the orlagh.log file is updated with
c     the result
c
      implicit double precision (a-h,o-z)

      double precision dat(36),data(36),x_initial, yi(10000), sig(10000)
      double precision one, y_initial, z_initial
      integer npar,IARGC,i, nsel, sub(10000), ia(36), nobs, nip, nse
      integer ell(10000), j
      character*80 par_string

c     internal arrays for gridmin
      double precision dpar(7),spar(7),nchi(9),zdiff1(5),zdiff2(4)
      double precision chisq, retx(7), ochisq, fchisq
      data (zdiff1(i), i=1,5)
     +/ -0.002,-0.001,0.,0.001,0.002 /
      data (zdiff2(i), i=1,4)
     +/ -0.0004,-0.0002,+0.0002,+0.0004/

c     common things
      integer nfreq, nonseis,subscripts(10000),iwflag,isflag,idif
      double precision rng(3), f0
      common/olcdint/ nfreq,nonseis,subscripts,iwflag,isflag,idif
      common/olcddbl/ rng

c     get data from the file.... these are initial parameters + data
      call mrqdata (ia, npar, dat, nse, yi, sig, ell)

c 
c     save the original parameters
c
      do i = 1, 7
        spar(i) = dat(i)
      enddo
      ochisq = 100000.
      nobs = 1200
      one = 1.0
 
c     set up the initial parameters and minimize for 5 Z values
      do i = 1, 5
         do j = 1,7     
            dat(j) = spar(j)
         enddo       
         dat(2) = dat(2) + zdiff1(i)    
         call reducemrq(dat,yi,sig,ia,subscripts, npar,nobs,chisq,retx)
         if (chisq.lt.ochisq) then
            szz = zdiff1(i)
            ochisq = chisq
            do j=1,7
               dpar(j) = retx(j)
            enddo
         endif  
      enddo
         
      do i = 1, 4
         do j=1,7 
            dat(j) = spar(j)
         enddo
         dat(2) = spar(2) + szz + zdiff2(i)
         call reducemrq(dat,yi,sig,ia,subscripts, npar,nobs,chisq,retx)
         if (chisq.lt.ochisq) then
            do j=1,7
               dpar(j) = retx(j) 
            enddo    
            ochisq = chisq
         endif   
      enddo
      
      print *,'new parameterss are'
      print *,dpar(1),dpar(2),dpar(3),dpar(4),dpar(5)

c
c    now i have the best  z answer, now lets do the same but for Y
c
      do i = 1, 7
         spar(i) = dpar(i)
      enddo
c
      do i = 1, 5
         do j = 1,7     
            dat(j) = spar(j)
         enddo       
         dat(3) = dat(3) + zdiff1(i)*5.    
         call reducemrq(dat,yi,sig,ia,subscripts, npar,nobs,chisq,retx)
         if (chisq.lt.ochisq) then
            szz = zdiff1(i)*5.
            ochisq = chisq
            do j=1,7
               dpar(j) = retx(j)
            enddo    
         endif       
      enddo
      do i = 1, 4
         do j=1,7 
            dat(j) = spar(j)
         enddo
         dat(3) = dat(3) + szz + zdiff2(i)*10.
         call reducemrq(dat,yi,sig,ia,subscripts, npar,nobs,chisq,retx)
         if (chisq.lt.ochisq) then
            do j=1,7
               dpar(j) = retx(j)
            enddo    
            ochisq = chisq
         endif   
      enddo
c
cc     test to see if we have a lower chisq, and if not just go to end
cc     otherwise go through z once more
c
c

      do i = 1, 7 
          spar(i) = dpar(i)
      enddo   

c     go to z again if we are still getting better answers
      do i = 1, 5
         do j=1,7 
            dat(j) = spar(j)
         enddo
         dat(2) = spar(2) + zdiff1(i)/5.
         call reducemrq(dat,yi,sig,ia,subscripts, npar,nobs,chisq,retx)
         if (chisq.lt.ochisq) then
            do j=1,7
               dpar(j) = retx(j) 
            enddo    
            ochisq = chisq
         endif   
      print *,dpar(1), dpar(2), dpar(3), dpar(4),dpar(5), ochisq
      print *,dat(1),dat(2), dat(3),dat(4), dat(5)
      enddo




   10 continue 
      
      print *,' =================================================='
      print *,'final answers are'
      print *,dpar(1), dpar(2), dpar(3), dpar(4), dpar(5),'='
      print *,'with chisq', ochisq
      print *,' =================================================='
      
 1    format (a25)   
 2    format (i10,x,f9.7,x,f9.7,x,f9.7,x,f7.5,x,d16.10,x)  
      open (16, FILE='orlagh.log', STATUS = 'OLD', ACCESS='append')
      write (16, 2) 'Final: ',dpar(1),dpar(2),dpar(3),dpar(4),dpar(5)
      close (16)
      
      end

