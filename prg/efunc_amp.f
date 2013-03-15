c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c     CALCULATE MODEL OBSERVABLES AND DERIVATIVES @ INPUT PARAMETERS
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine efunc (ma, ia, x, nobs, sub, ymod, dyda)

      implicit none
      integer ma, nobs, sub(nobs), ia(ma)
      double precision x(ma), dyda(nobs, ma), ymod(nobs)

c     internal
      integer n, i, j, ndata
      double precision dp(ma), ymod2(nobs,15), rx(ma,15)
      double precision cc(ma), age_tmp
c
c     define the parameter steps to make for numerical derivatives
c      
      ndata = 0
      do i=1,nobs
         if (sub(i).ne.0) ndata=ndata+1
      enddo

      dp(1) = 0.001d0
      dp(2) = 0.001d0
      dp(3) = 0.001d0
      dp(4) = 0.01d0
      dp(5) = 0.005d9
      dp(6) = 1.0d0
      dp(7) = 1.0d0
c
c populate array with all parameter sets
c
      do n=1,2*ma+1
         do i=1,ma
            if (ia(i).ne.0) then
               rx(i,n) = x(i)
               if (n .eq. i) rx(i,n) = x(i) + dp(i)
               if (n .eq. i+ma) rx(i,n) = x(i) - dp(i)
            endif
         enddo
      enddo
c      
c evaluate models in parallel
c
      call mpi_evalmodel (ma, rx, nobs, ndata, ymod2)

      cc(1) = 1.0d0
      cc(2) = 1.0d0
      cc(3) = 1.0d0
      cc(4) = 1.0d0
      cc(5) = 1.0d9
      cc(6) = 1.0d0
      cc(7) = 1.0d0
      do i=1,ma
         do j=1,nobs
            if (sub(j).ne.0) then 
               dyda(j,i) = (ymod2(j,i) - ymod2(j,i+ma))/(2.d0*dp(i))
c scale the derivatives to avoid inversion errors.
               dyda(j,i) = dyda(j,i)*cc(i)
            endif
c put central model into ymod array
            if (i.eq.ma) ymod(j) = ymod2(j,2*ma+1)
         enddo
      enddo

      return
      end

c**********************************************************************

      subroutine mpi_evalmodel (ma, rx, nobs, ndata, ymod2)
c ---------------------------------------------
c     parallel evaluation using MPI
c ---------------------------------------------
      implicit none

      include 'mpif.h'
      
      integer i,j,ma,nobs,ndata
      integer ndone, nproc, ierr, nspawn
      integer msgtype, job, num_jobs, trial, slave
      integer npar, status(MPI_STATUS_SIZE)
      double precision data(32),rx(ma,15),ymod2(nobs,15),yy(nobs)
      logical receiving

c ---------------------------------------------
c     initialize counter
c ---------------------------------------------
      ndone = 0

c ---------------------------------------------
c     determine number of processors
c ---------------------------------------------
      call mpi_comm_size( MPI_COMM_WORLD, nproc, ierr )

c ---------------------------------------------
c     run jobs on slave nodes only
c ---------------------------------------------
      nspawn=nproc-1

c ---------------------------------------------
c     send an initial job to each node
c ---------------------------------------------
      do job=1,nspawn
         trial = job
         slave = job
c fill data with parameters from rx
         do i=1,ma
            data(i) = rx(i,trial)
         enddo
         call sendjob(trial,slave,ma,data,nobs,ndata)
      enddo

      num_jobs=2*ma+1
      do job=1,num_jobs
c ---------------------------------------------
c     listen for responses
c ---------------------------------------------
 25      msgtype = 2 
         call mpi_iprobe( MPI_ANY_SOURCE, msgtype, MPI_COMM_WORLD,
     +                    receiving, status, ierr )

         if (receiving) then
c ---------------------------------------------
c     get data from responding node
c ---------------------------------------------
            call mpi_recv( slave, 1, MPI_INTEGER, MPI_ANY_SOURCE,
     +                     msgtype, MPI_COMM_WORLD, status, ierr )
            call mpi_recv( trial, 1, MPI_INTEGER, slave,
     +                     msgtype, MPI_COMM_WORLD, status, ierr )
            call mpi_recv( nobs, 1, MPI_INTEGER, slave,
     +                     msgtype, MPI_COMM_WORLD, status, ierr )
            call mpi_recv( yy, nobs, MPI_DOUBLE_PRECISION, slave,
     +                     msgtype, MPI_COMM_WORLD, status, ierr )

            write(*,*) "job",trial,rx(1,trial),rx(2,trial),
     +                 rx(3,trial),rx(4,trial),rx(5,trial)

c fill ymod2 with results from yy
            do j=1,nobs
               ymod2(j,trial) = yy(j)
            enddo
            ndone = ndone + 1

c ---------------------------------------------
c     send new job to responding node
c ---------------------------------------------
 140        if (ndone .LE. (num_jobs-nspawn)) then
               trial = job + nspawn
c fill data with parameters from rx
               do i=1,ma
                  data(i) = rx(i,trial)
               enddo
               call sendjob(trial,slave,ma,data,nobs,ndata)
            endif
            goto 100
         endif

c ---------------------------------------------
c     return to listen again or move on
c ---------------------------------------------
         if (.NOT. receiving) goto 25

         goto 199
 100     continue
      enddo

c ---------------------------------------------
c     ready for next set of jobs
c ---------------------------------------------

 199  return
      end

c**********************************************************************
      subroutine sendjob(trial,slave,npar,data,nobs,nsel)

      implicit none

      include 'mpif.h'

      integer trial, slave, npar, nobs, nsel, ierr, msgtype
      double precision data(32)

      msgtype = 1
      call mpi_send( trial, 1, MPI_INTEGER, slave,
     +               msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( npar, 1, MPI_INTEGER, slave,
     +               msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( data, npar, MPI_DOUBLE_PRECISION, slave,
     +               msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( nobs, 1, MPI_INTEGER, slave,
     +               msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( nsel, 1, MPI_INTEGER, slave,
     +               msgtype, MPI_COMM_WORLD, ierr )

      return
      end

c**********************************************************************
