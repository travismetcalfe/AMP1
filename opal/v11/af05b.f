        program af05b
        implicit double precision (a-h,o-z)
c
c       purpose:
c               read inidividual alexander-ferguson (2005 format) 
c               opacity tables,  and create one ASCII and binary 
c               output file.
c
c       History: 
c
c       18/08/97:  creation from s/r alex2b.f
c
c       last modification: 22/02/06
c
c       Variables:
c       lun ....... logical unit number of opacity-table-file
c       ntab ...... nr. of different x-values
c       nval ...... max. nr. of op.values pro line (= # of rlg-values)
c       nlin ...... max. nr. of op.table lines (= # of log(T)=values)
c       nzva ...... nr. of different z-values
c                   (z = mass fraction of heavy elements)
c
c       rlg ....... array[1..nval,1..ntab], 
c                   decade log of r=density(gm/cm**3)/t6**3
c       tlg ....... array[1..nlin,1..ntab].
c                   decade log of temperature
c       opa ....... array[1..nvar,1..nlin,1..ntab],
c                   opacity values
c
c       setup array-dimensions
c
      parameter(ntab=8,nzva=13)
      parameter(nval=19,nlin=85)
      parameter(lun=21)
c
      dimension       rlg(nval)
      dimension       opa(nval)

      character*132   line
      character*8     tbname
      character*6     ofname
      character*6     ct6
      character*80    ifname
c
c     initialize input filenames
      data tbname  /'AF05TABS'/
c
c-----initialize output filename
      data ofname  /'af.bin'/
c
c     open outputfile
      lun1=lun+1
      close(lun1)
      open(lun1,file=ofname,status='unknown',
     .          form='unformatted',err=9011)
c
c     open file with names of tables (individual tables)
      lun2=lun+2
      close(lun2)
      open(lun2,file=tbname,status='old',err=9101)
c
c-----read tables
      do 1015 k=1,ntab
        do 3001 l=1,nzva
c         open inputfile (one single table)
          read(lun2,'(a)')ifname
          close(lun)
          open(lun,file=ifname,status='old',err=9001)
c         read 1st 3 lines 
          do 1001 i=1,3
           read(lun,'(a)') line
           if (i.eq.1) print '(a,$)',line(1:67)
 1001     continue      
c
c         print *,'read rlg header'
c         read/write string 'logT' + rlg-values
          read(lun,   '(a6,19f7.3)') ct6,(rlg(i),i=1,nval)
          write(lun1               )     (rlg(i),i=1,nval)
c
c         print *,'read opacity values'
c         read/write logT & opacity-values
          do 1010 j=1,nlin
            read(lun,'(f6.3,19f7.3)',err=9002)tlg,(opa(i),i=1,nval)
            if(j.eq.13) print '(a,f4.2)',' tlgmax=',tlg
            if(j.ge.13)       ! ignore first 12 line (start at logT=3.900)
     +      write(lun1,              err=9020)tlg,(opa(i),i=1,nval)
 1010     continue
c         print *,'af05b: ',j-1,' lines read in file ',ifname   

 2001       continue
 3001   continue
        close(lun)
 1015   continue
c
c
        close(lun1)
        close(lun2)
        stop
c
 9001   print *,'af05b: error in opening ',ifname,' lun= ',lun
        stop
 9002   print *,'af05b: error in reading: line: ',j
        stop
 9011   print *,'af05b: error in opening ',ofname,' lun= ',lun1
        stop
 9020   print *,'af05b: error in writing ',ofname,' lun= ',lun1 
        stop
 9101   print *,'af05b: error in opening ',tbname,' lun= ',lun2
        stop
c
        end
