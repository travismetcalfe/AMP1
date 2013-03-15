        program af94X9intZ
        implicit double precision (a-h,o-z)
c
c       purpose:
c               read X=9, Z=(0.4, 0.7, 1.0) Alexander-Ferguson (1994 format) 
c               opacity tables,  create tables for Z=(0.6, 0.8) using
c               Akima interpolation, and write results into new 
c               new ASCII tables: g9.06.tron and g9.08.tron
c               
c
c       History: 
c
c       18/08/97:  creation from s/r alex2b.f
c
c       last modification: 27/02/06
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
      parameter(ntab=1,nzva=3)
      parameter(nval=17,nlin=23)
      parameter(lun=21)
c
      dimension       rlg(nval)
      dimension       opa(5,nval)
c
      dimension       zs(3)
      dimension       zi(2), di(2)
c
      character*125   line(3)
      character*60    tbname
      character*60    ofname(2)
      character*6     ct6
      character*80    ifname(3)
c
c     initialize input filenames
      data tbname  /'AF94X9TABS'/
c
c-----initialize output filenames
      data ofname  /'g9.06.tron','g9.08.tron'/
c-----degree of Akima polynomial 
      data np /3/
c-----number of table point for Akima interpolation (Z=0.4,0.7,1.0)
      data ns /3/
c-----number of output point for Akima interpolation (Z=0.6,0.8)
      data ni /2/
c-----table point for Akima interpolation
      data zs /0.04d0,0.07d0,0.1d0/    ! given points
      data zi /0.06d0,0.08d0/        ! requested points
c
      lun1=lun+1
      lun2=lun+2
      lun3=lun+3
      lun4=lun+4
      lun5=lun+5
      lun6=lun+6
      lun7=lun+7
c
c     open file with names of individual tables
      close(lun1)
      close(lun2)
      close(lun3)
      close(lun4)
      close(lun5)
      close(lun6)
      open(lun6,file=tbname,status='old',err=9101)
c
c     open input and output files
      read(lun6,'(a)')ifname(1)
      read(lun6,'(a)')ifname(2)
      read(lun6,'(a)')ifname(3)
      close(lun6)
c
      open(lun1,file=ifname(1),status='old',err=9001)
      open(lun2,file=ifname(2),status='old',err=9002)
      open(lun3,file=ifname(3),status='old',err=9003)
      open(lun4,file=ofname(1),status='unknown',err=9004)
      open(lun5,file=ofname(2),status='unknown',err=9005)
c
c-----read tables
c         read 1st 4 lines 
          do i=1,4
           read(lun1,'(a)') line(1)
           read(lun2,'(a)') line(2)
           read(lun3,'(a)') line(3)
           if (i.eq.1) then
              print '(a)',line(1)
              print '(a)',line(2)
              print '(a)',line(3)
              write(lun4,'(a26,f4.2)') line(2),zi(1)
              write(lun5,'(a26,f4.2)') line(2),zi(2)
           else
              write(lun4,'(a45)') line(2)
              write(lun5,'(a45)') line(2)
           endif
          enddo      
c
          print *,'reading log(R) values'
c         read string 'logT' + rlg-values
          read (lun1,'(a6,17f7.2)') ct6,(rlg(i),i=1,nval)
          read (lun2,'(a6,17f7.2)') ct6,(rlg(i),i=1,nval)
          read (lun3,'(a6,17f7.2)') ct6,(rlg(i),i=1,nval)
          write(lun4,'(a4,17f7.2)') ct6,(rlg(i),i=1,nval)
          write(lun5,'(a4,17f7.2)') ct6,(rlg(i),i=1,nval)
c
c         read next 2 lines 
          do i=1,2
           read (lun1,'(a)') line(1)
           read (lun2,'(a)') line(2)
           read (lun3,'(a)') line(3)
           write(lun4,'(a1)') line(1)
           write(lun5,'(a1)') line(2)
          enddo      
c
          print *,'reading opacity values'
c         read logT & opacity-values
          do j=1,nlin
            read(lun1,'(f4.2,17f7.3)',err=9006)tlg1,(opa(1,i),i=1,nval)
            read(lun2,'(f4.2,17f7.3)',err=9006)tlg2,(opa(2,i),i=1,nval)
            read(lun3,'(f4.2,17f7.3)',err=9006)tlg3,(opa(3,i),i=1,nval)
            if(tlg1.ne.tlg2)then
              print *,'tlg1.ne.tl2'
              stop
            endif
            if(tlg1.ne.tlg3)then
              print *,'tlg1.ne.tl3'
              stop
            endif
            do i=1,nval
               call uvip3d(np,ns,zs,opa(1,i),
     .                        ni,zi,opa(4,i),di)
            enddo
            write(lun4,'(f4.2,19f7.3)',err=9007)tlg1,(opa(4,i),i=1,nval)
            write(lun5,'(f4.2,19f7.3)',err=9007)tlg1,(opa(5,i),i=1,nval)
          enddo
c
c
        close(lun1)
        close(lun2)
        close(lun3)
        close(lun4)
        close(lun5)
        stop
c
 9001   print *,'af94X9intZ: error in opening ',ifname(1),' lun= ',lun
        stop
 9002   print *,'af94X9intZ: error in opening ',ifname(2),' lun= ',lun
        stop
 9003   print *,'af94X9intZ: error in opening ',ifname(3),' lun= ',lun
        stop
 9004   print *,'af94X9intZ: error in opening ',ofname(1),' lun= ',lun
        stop
 9005   print *,'af94X9intZ: error in opening ',ofname(2),' lun= ',lun
        stop
 9006   print *,'af94X9intZ: error in reading: line: ',nlin
        stop
 9007   print *,'af94X9intZ: error in writing: line: ',j
        stop
 9101   print *,'af94X9intZ: error in opening ',tbname
        stop
c
        end
