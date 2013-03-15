        program af05b
        implicit double precision (a-h,o-z)
c
c       purpose:
c               read alexander-fergunson (1994) opacity tables extrapolate 
c               for log(R)=[-8,-7.5], add additional 54 (dummy) log(T) values
c               to be compatible with (2005)-tables, and write tables in
c               binary format to disk.
c
c       History: 
c
c       18/08/97:  creation from s/r alex2b.f
c
c       last modification: 28/02/06
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
      parameter(nval=19,nlin=23)
      parameter(lun=21)
c
      dimension       rlg(nval)
      dimension       opa(nval)
      dimension       di(2)

      character*132   line
      character*8     tbname
      character*6     ofname
      character*6     ct6
      character*80    ifname
c
c     initialize input filenames
      data tbname  /'AF94TABS'/
c
c-----initialize output filename
      data ofname  /'af.bin'/
c
c-----data for Akima extrapolation
      data rlg(1),rlg(2) /-8.0d0,-7.5d0/  ! requested data points
      data ni  /2/              ! number of requested data points
      data np  /3/              ! degree of polynomial for s/r uvip3d
c
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
c         read 1st 4 lines 
          do i=1,4
           read(lun,'(a)') line
           if (i.eq.1) print '(a,$)',line(1:67)
          enddo      
c
c         print *,'read rlg header'
c         read/write string 'logT' + rlg-values
          read(lun,   '(a6,17f7.2)') ct6,(rlg(i),i=3,nval)
          write(lun1               )     (rlg(i),i=1,nval)
c         print       '(a4,19f7.2)', ct6,(rlg(i),i=1,nval)
c
c         read next two blank lines
          do i=1,2
           read(lun,'(a)') line
          enddo      
c
c         print *,'read opacity values'
c         read/write logT & opacity-values
          ns=5               ! number of given table points for Akima
          do j=1,nlin
            read(lun,'(f4.2,17f7.3)',err=9002)tlg,(opa(i),i=3,nval)
            if(j.eq.5) print '(a,f4.2)',' tlgmax=',tlg
            if(j.ge.5) then    ! ignore first 4 lines (start at logT=3.90)
              if(j.gt.20)ns=2  ! makes extrapolation more `reasonable'
              call uvip3d(np,ns,rlg(3),opa(3),
     .                       ni,rlg(1),opa(1),di)
              write(lun1,              err=9020)tlg,(opa(i),i=1,nval)
c             print  '(f4.2,19f7.3)',           tlg,(opa(i),i=1,nval)
            endif
          enddo
c---------add additional dummy lines to be compatibLe with AF 2005 tables
          do j=1,54
            tlg=tlg-0.01d0
            write(lun1,              err=9020)tlg,(opa(i),i=1,nval)
c           print  '(f4.2,19f7.3)',           tlg,(opa(i),i=1,nval)
          enddo
 3001   continue
        close(lun)
 1015 continue
c
c
        close(lun1)
        close(lun2)
        stop
c
 9001   print *,'af94bext: error in opening ',ifname,' lun= ',lun
        stop
 9002   print *,'af94bext: error in reading: line: ',j
        stop
 9011   print *,'af94bext: error in opening ',ofname,' lun= ',lun1
        stop
 9020   print *,'af94bext: error in writing ',ofname,' lun= ',lun1 
        stop
 9101   print *,'af94bext: error in opening ',tbname,' lun= ',lun2
        stop
c
        end
