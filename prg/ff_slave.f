      subroutine ff_slave 
c ---------------------------------------------
c fitness function slave program
c ---------------------------------------------
      implicit none

      include 'mpif.h'

      integer myid, ierr, status(MPI_STATUS_SIZE)
      integer master, msgtype, trial, n, i
      double precision data(32),result,age,R_Ro,Teff,chisq(4)
      real userff
      logical quiet
      common /verbosity/ quiet
      common /xmodage/ age,R_Ro,L_Lo,Teff,M_H,chisq

      external userff
      quiet = .TRUE.

c ---------------------------------------------
c identify this slave task
c ---------------------------------------------
      call mpi_comm_rank( MPI_COMM_WORLD, myid, ierr )

c ---------------------------------------------
c listen for a new job
c ---------------------------------------------
      master = 0
 25   msgtype = 1

c ---------------------------------------------
c receive data from master host
c ---------------------------------------------
      call mpi_recv( trial, 1, MPI_INTEGER, master,
     +               msgtype, MPI_COMM_WORLD, status, ierr )
      if (trial .EQ. -1) goto 99
      call mpi_recv( n, 1, MPI_INTEGER, master,
     +               msgtype, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( data, n, MPI_DOUBLE_PRECISION, master,
     +               msgtype, MPI_COMM_WORLD, status, ierr )

      do i=1,n
        data(i) = INT((100.*data(i))+0.5)/100.d0
      enddo
c ---------------------------------------------
c perform calculations with data
c ---------------------------------------------
      result = userff( n, data, myid )

c ---------------------------------------------
c send result to master host
c ---------------------------------------------      
      msgtype = 2 
      call mpi_send( myid, 1, MPI_INTEGER, master,
     +               msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( trial, 1, MPI_INTEGER, master,
     +               msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( result, 1, MPI_DOUBLE_PRECISION, 
     +       master, msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( age, 1, MPI_DOUBLE_PRECISION, 
     +       master, msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( R_Ro, 1, MPI_DOUBLE_PRECISION, 
     +       master, msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( L_Lo, 1, MPI_DOUBLE_PRECISION, 
     +       master, msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( Teff, 1, MPI_DOUBLE_PRECISION, 
     +       master, msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( M_H, 1, MPI_DOUBLE_PRECISION, 
     +       master, msgtype, MPI_COMM_WORLD, ierr )
      call mpi_send( chisq, 4, MPI_DOUBLE_PRECISION, 
     +       master, msgtype, MPI_COMM_WORLD, ierr )

c ---------------------------------------------
c go back for more work
c ---------------------------------------------
      goto 25

 99   return
      end

c*********************************************************************
