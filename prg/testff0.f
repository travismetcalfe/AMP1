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

      if ((data(1).lt.0.75).or.(data(2).gt.0.05).or.(data(4).lt.1))
     + stop "please enter physical units"

c  physical units => single model evaluation
      myid=0
      quiet = .FALSE.

      fitness = userff(npar,data,myid)

      end
