      Program exR5
      implicit double precision (a-h,o-z)
      common/cstdio/ istdin, istdou, istdpr, istder
c
c     extrapolate OPAL95 tables to log(R)=5
c     using linear extrapolation (Akima: uvip3p.f)
c
c     hg: 10/04/2003
c
      parameter(nx=7, nz=13)
      parameter(nr=19,nt=70)
      parameter(mr=27,mt=30)
      parameter(lun=31,nd=5,n8=8,nr1=nr+1)
      character*40 optabe,paderi,ivadat,optabee
c
      dimension rlg95(nr),rlg95e(mr),rlge(n8),tlg95(nt)
c
c     dimensions for OPAL95 tables
      dimension rlg(nr,nx,nz)
      dimension tlg(nt,nx,nz)
      dimension opa(nr,nt,nx,nz)
      dimension opae(mr,nt,nx,nz)
      dimension iva(nr,nt,nx,nz)
c
c     dimensions for Kurucz tables (kuru91_02X.bin)
      dimension tlk(mt)
      dimension rlk(nz,mr*mt)
      dimension opk(nz,mr*mt)
      dimension opkT4(8)
c
c     stupidly necessary for s/r rdi95
      common /tabdim/ nzvai,ntabi,nxiri,nyiri,nyisi,nyif,mdi,
     +                nti,iali
c
      nzvai=nz
      ntabi=nx
      nxiri=nr
      nyiri=nt
      mdi  =mr
c
      open(29,file='OPINTPATH_95XE',
     +     form='formatted',status='old')
      read(29,'(A)')optabe
      read(29,'(A)')paderi
      read(29,'(A)')ivadat
      read(29,'(A)')optabee
      close(29)
      write(istdpr,*)'Using tables:'
      print '(a,a)',optabe,'  (input)'
      print '(a,a)',optabee,' (output)'
c
c     read tables
      call rdi95(lun,optabe,nx,nr,nt,nz,rlg,tlg,opa,iva)
c
c     setup the opal95 T values
      do 3510 i=1,nt
       tlg95(i)=tlg(i,nx,nz)
3510  continue
c
c     setup extended opa95 R values
      do 3515 i=1,nr
       rlg95(i) =rlg(i,nx,nz)
       rlg95e(i)=rlg95(i)
3515  continue
      do 3516 i=1,n8
       rlge(i)     =1.0d0+float(i)/2.0d0
       rlg95e(nr+i)=rlge(i)
3516  enddo
c
c     do the extrapolation
      do iz=1,nz
        do ix=1,nx
          do it=1,nt
             call uvip3p(nd,nr,rlg95,opa (  1,it,ix,iz),
     +                      n8,rlge, opae(nr1,it,ix,iz))
             do ir=1,nr
                opae(ir,it,ix,iz)=opa(ir,it,ix,iz)  ! copy remaining values
             enddo ! ir
          enddo    ! it
        enddo      ! ix
      enddo        ! iz
c
c
c     output all extrapolated OPAL95 tables (format as in exakop95.f)
      open(lun+2,file=optabee,form='unformatted')
      do ix=1,nx
        do iz=1,nz
          write(lun+2)(rlg95e(ir),ir=1,mr)     ! mr=nr+n8
          do it=1,nt
            write(lun+2)tlg95(it),(opae(ir,it,ix,iz),ir=1,mr)
          enddo
        enddo ! iz
      enddo   ! ix
      close(lun+2)
c
c     output one extrapolated OPAL95 table
      ix=6        ! X=0.7
      iz=8        ! Z=0.02
      open( 8, form='unformatted')
      open( 9, form='unformatted')
      write(8) mr,(rlg95e(k),k=1,mr)
      write(9) nt,(tlg95(k), k=1,nt)
      close(8)
      close(9)
      open(80,form='unformatted')
      write(80)((opae(ir,it,ix,iz),ir=1,mr),it=1,nt)
      close(80)
c      
      stop
      end
