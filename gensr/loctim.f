      function loctim()
      character*24 loctim
      character*24 timest,fdate
c
c  returns date and time as a string. Note this routine
c  may depend on installation.
c
      timest= fdate()
c.. This needs update      loctim=timest
c..      loctim='Dummy_date_and_time'
      return
      end
