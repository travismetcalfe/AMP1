      program main
      character*80 file1, file2,ss
      dimension a(8),b(6),f(3,10)
      integer l1,l2,n1,n2,t
      real r,r0
c
      open(20,file='bcsum.0180.list',status='old')
      open(50,file='obs.0180.list',status='old')
      open(40,file='mode_f.0180.out',status='unknown')
      open(70,file='mode_r.0180.out',status='unknown')
     
c
      print *,'(write spherical function degree interval l1:l2 
     *egz.:001:002 or 1 2 or 1-2)'
      read(*,*) l1,l2
      write(6,'(''Choosen degree interval'',i3,'':'',i3)')l1,l2
      
      print *,'(write spherical function order interval n1:n2 
     *egz.:001:002 or 1 2 or 1-2)'
      read(*,*) n1,n2
      write(6,'(''Choosen order interval'',i3,'':'',i3)')n1,n2
     
      write(40,'(''   M/Msun   Age G_yr  R cm          Teff      '',
     *''L/Lsun   Xc        qc         logTeff  logL/Lsun  nu(n,l) 
     *(e-6 Hz), l, n ='',20i2)')
     *  ((i,j,j=n1,n2),i=l1,l2)
     
      write(70,'(''   M/Msun   Age G_yr  R cm          Teff      '',
     *''L/Lsun   Xc        qc         logTeff  logL/Lsun  nu(n,l) 
     *(e-6 Hz), l, n ='',20i2)')
     *  ((i,j,j=n1,n2),i=l1,l2)
      
      t=1
   10 read(20,*,end=20) file1
      open(30,file=file1,status='old')
      read(30,*) ss
      read(30,*) ss
      read(30,*) a
      read(50,*) file2
      open(35,file=file2,status='old')
      read(35,*) ss
      read(35,*) ss
      read(35,*) ss
      read(35,*) ss
      read(35,*) ss
      read(35,*) ss
      
      
      do i=1,3
        do j=1,10
	f(i,j)=0.0
	enddo
      enddo
      
      
      if (t.eq.1) then
         r0=a(4)	 
      endif
    
   15 read(35,*,end=16) b
      
            
      if (l1.le.b(1).and.b(1).le.l2.and.n1.le.b(2).and.b(2).le.n2) then
      f(1+b(1),b(2))=b(3)
      endif
      
      go to 15
      
   16 write(40,'(f10.5,f10.6,1pe13.5,0pf10.2,f10.4,f10.6,f10.6,
     *30f10.3)') 
     *(a(i),i=2,8),Log10(a(5)),Log10(a(6)),
     *((f(i,j),j=n1,n2),i=l1+1,l2+1)
      
      do i=1,3
        do j=1,10
      	f(i,j)=f(i,j)*((a(4)/r0)**(1.5))
      	enddo
      enddo
      
      
      write(70,'(f10.5,f10.6,1pe13.5,0pf10.2,f10.4,f10.6,f10.6,
     *30f10.3)')
     *(a(i),i=2,8),Log10(a(5)),Log10(a(6)),
     *((f(i,j),j=n1,n2),i=l1+1,l2+1)
      
      t=2
      go to 10
   20 continue
      stop
      end
