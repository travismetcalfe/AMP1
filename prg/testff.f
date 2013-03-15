      program testff

      implicit none

      double precision data(36),fitness,age
      real userff
      integer npar,IARGC,i,myid
      character*80 par_string
      logical quiet

      common /verbosity/ quiet
      common /xmodage/ age

      npar = IARGC()
      do i=1,npar
         CALL GETARG(i,par_string)
         READ(par_string,*) data(i)
      enddo

      if ((data(1).gt.1).or.(data(4).gt.1))
     + stop "please enter scaled units"

c  scaled units => fitting data
      myid=1
      quiet = .FALSE.

c  initialize routine using best model
      fitness = userff(npar,data,myid)

c  reproduce best model after initialization
      fitness = userff(npar,data,myid)

      end
