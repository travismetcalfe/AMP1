      program main
c
c  NOTE: all setups have been moved to s/r setups_main.
c
      implicit double precision (a-h,o-z)
c
      common/cversn/ ivrbcs,ivrrhs,ivreos,ivropc,ivreng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      call setups_main
c
c  for now, set i_paramset = 0, to flag for no repeated parameter
c  setting.
c
      i_paramset = 0
c
      call mnevol(i_paramset, ierr_param)
c
      stop
      end
