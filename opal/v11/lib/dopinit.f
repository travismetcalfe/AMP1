c      algorithm 677, collected algorithms from acm.
c      this work published in transactions on mathematical software,
c      vol. 15, no. 4, pp. 365-374.
c
      subroutine opinit(eps,iord,tabnam,imode)
      implicit double precision (a-h,o-z)
c
c  Initialize (read) Opazity-tables ,partial derivative values
c  and index-coordinates
c
c     History:
c     10.6.1992: creation
c     25.8.1992: reduce of iwk-      & wk-arrays to
c                          iwk(3*nt) & wk(2*n) respectively
c                where nt = # of triangles (1568) and n = # of datapoints
c                (nxi*nyi = 17*50 = 850)
c     30.9.1992: included 'status='old'' in open statement for VAX/VMS
c
c     20.12.1992 modified for z-domain interpolation
c                nzva different z-values
c                (z = mass fraction of heavy elements)
c
c                modified opacity values input from
c                'formatted' to 'unformatted'
c                first you have to run the program 'a2b'
c                (ascii-to-binary conversion of opacity-tables)
c
c      25.1.1993 changed variable t6 -> tlg = log10(t)
c                table-values changed to tlg in da2b
c
c       5.2.1993 modified pathname-treatment of inputfiles.
c                A single file will be used to define the
c                absolute pathnames of the 3 inputfiles
c                (optabe.bin, pderivs.dat, ival.dat) line by line
c                in the above order !
c                This file has to be placed in the
c                current working directory and must have
c                the filename OPINTPATH, or may be
c                assigned to a logical name: eg. 
c                under UNIX in a csh:
c              
c                setenv OPINTPATH /dir1/dir2/opacity.paths
c
c       5.3.1993 redefined pathname-treatment
c                defined new argument tabnam, which
c                defines pathname of file, which table
c                has to be used
c
c      11.3.1993 introduced parameter statement
c                for array-dimensions
c                removed ntab from common /tablex/
c                removed nzva from common /tablez/
c                add new common /arrdim/ for
c                consistency-check of the
c                common - array defs, because parameters
c                can not be passed through common-statement
c
c     12.3.1993  modified pderiv.dat read- handling
c                for use for the Kurucz-tables
c                modified common /jpoint/
c                         common /ipoint/
c
c     13.3.1993  use 'D' in 1-st column for
c                debugging commands
c
c     17.5.1993  use d(opa)/dX-extrapolated Kurucz-tables for all X-values
c                (s/r readkx, array opkx)
c                changed parameter md = 17
c
c     20.5.1993: iwk-array (triangle indices) will only be
c                read once (triangulation is the same for
c                all tables; independend from X and Z)
c
c     25.5.1993: included variable 'iorder' in argument list
c                and added iorder to common /tablex/
c
c      6.6.1993: change read-sequence of input-files (tabnam)
c
c      8.6.1993: mt=22 -> mt=31 (index of new Kurucz fittingpoint tk_fit=4.05)
c                (must be equal to ylo - output of kur2l2g)
c
c                introduced new parameter nyis, defining
c                the starting # of OPAL-temperature, for 
c                fit with Kurucz-table at tk_fit=4.05
c                (nyis=7 for t6=0.012 ->  tl_fit=4.079)
c
c                parameter nt   -> 2368 (from dopxext_s.f output)
c                          ndat -> 1275
c
c      7.11.1993: added H. Spath's birational spline algorithm
c                 on a regulary grid s/r rat2d & rbivpd (s/r dopints)
c                 !!!!!!! B U T  !!!!!!!!!
c                 the stored coefficients array ra will allocate
c                 > 10MBte ! (17*75*4*4*7*10)*8 byte
c
c                 new argument parameter imode
c
c                 if imode = 0 -> only minimum norm (677) algorithm (s/r opint{c,f})
c                    imode > 0 -> birational splines + 677 (s/r opints + opint{c,f})
c
c      8.11.1993: change drelpr->10.d0*drelpr for imode >0  prevent 'division by zero'
c                 for the analytically evaluated partial derivatives wrt x
c
c      28.1.1994:  fit with Kurucz-table at tk_fit=4.0 (mt=30)
c                 (nyis=6 for t6=0.011 ->  tl_fit=4.041)
c
c      13/11/95: modified for OPAL95 tables
c
c      23/05/97: new array ztabl = log10(ztab), used for 
c                logarithmic interpolation in Z
c
c     17/08/97: include electron conduction (new variable lec [logical])
c
c     19/08/97: include Alexander95 low-temperature tables
c
c     07/06/02: include 'external itohco' to initialize coefficients for
c               computing contribution due to electron conduction
c
c     04/04/03: check iord<=6
c
c     09/05/03: modified for using extrapolated OPAL95 and Kurucz tables
c
c     24Feb06 :modified for using extrapolated OPAL & "Alex-Ferg 2005" tables
c              (Kurucz tables are no longer supported)
c
c     16Aug10: Add flag for initialization of opacity tables in
c              common/c_opinit/
c
c     Last modification:
c     16Aug10
c
c  input parameters
c  
c  eps ............... relative machine precision
c  iord  ............. nr. of used tables for Akima-int. [2<=iorder<=6]
c  tabnam ............ filename, which tables to be used
c  imode ............. selection: if|imode|= 1 -> minimum norm (677) only
c                                   |imode|= 2 -> birational splines only
c                                   |imode|> 2 -> birational splines + 677
c                                    imode < 0 -> disable electron conduction
      external itohco
c
c     setup array-dimensions (max values for OPAL95 tables)
c
      parameter(ntab=8,nzva=13)
      parameter(nxir=27,nyir=70,nyis=5)     ! OPAL, tlg(nyis)=3.95
      parameter(ntal=73)                    ! Alexander-Ferguson 2005 (AF05)
      parameter(ndat=(nyir-nyis+1+ntal)*nxir)
c     nyii=nyir-nyis+1+ntal
c     nt=2*(nxir+nyii-2)+2*(nxir*nyii-2*(nxir+nyii-2)-1)
c     parameter(nt=3348,ial=2*ndat+5)       !ial=7511 for OPAL+AF05
      parameter(nt=7176)                    ! nt=7176 for OPAL+AF05
c     niwkl=3*nt,nwkl=2*nd=ial-5
      parameter(niwkl=3*nt,nwkl=2*ndat)
c
c---  dimension for s/2 rat2d
c 
      parameter(nyif=nyir+ntal-nyis+1)
c
      real ra
      logical*1 lec                  ! logical for electron conduction
      character*(*) tabnam
      character*256 optabe,pderiv,ivadat,lowtab
      dimension ival(nyir,ntab,nzva) ! for OPAL tables
c
c     dimension for Alexander-Ferguson-tables
c
      dimension opaal(nxir,ntal,ntab,nzva)
      dimension rlgal(nxir,ntab,nzva)
      dimension tlgal(ntal,ntab,nzva)
c
c --- definitions for s/r rat2d (rational splines)
c
      dimension p(nxir,nyif),
     +          q(nxir,nyif),
     +          r(nxir,nyif),
c    +          ra(nxir,nyif,4,4),
     +          rsdx(nxir),rsdy(nyif),
     +          rsax(nxir),rsbx(nxir),rscx(nxir),rsrx(nxir),
     +          rsay(nyif),rsby(nyif),rscy(nyif),rsry(nyif),
     +          ops(nxir,nyif,ntab,nzva)
c
c     partial derivatives and triangulation indices
      common /pderiv/ iwk(niwkl),wk(nwkl,ntab,nzva)
c
      common /opadat/ opa(nxir,nyir,ntab,nzva),
     +                rlg(nxir,ntab,nzva),
     +                tlg(nyir,ntab,nzva)
      common /valdat/ ivalo(nyir,ntab,nzva)
      common /tablex/ xtab(ntab),iorder,lec
      common /tablez/ ztab(nzva),ztabl(nzva)
      common /xyzdat/ xd(ndat,nzva),
     +                yd(ndat,nzva),
     +                zd(ndat,ntab,nzva)
      common /machin/ drelpr,toll,eps10
c
c
c---- birational spline coefficients for s/r rbivpd
c---- evaluated in s/r rat2d (s/r opinit)
c
      common /birasp/ pp,tls(nyif),
     +                ra(nxir,nyif,4,4,ntab,nzva)
c
c     common for array-consistency-check of common-arrays
c     between s/r opinit and s/r opintd
d     common /arrdim/ nxiv,nyiv,ndatv,ntabv,nzvav,niwkv,nwkv
c
c     table dimension for s/r opintc{f} and opints
      common /tabdim/ nzvai,ntabi,nxiri,nyiri,nyisi,iyifi,
     +                nti,iali
c
c  Flag for opacity initialization (reset in call of opintc or 
c  opints)
c
      common/c_opinit/ i_opinit
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      data xtab(1),xtab(2),xtab(3),xtab(4),
     +     xtab(5),xtab(6),xtab(7),xtab(8)
     +    /0.0d0,0.1d0,0.2d0,0.35d0,
     +     0.5d0,0.7d0,0.8d0,0.9d0/
c
      data ztab(1),ztab(2),ztab(3),ztab(4),ztab(5),
     +     ztab(6),ztab(7),ztab(8),ztab(9),ztab(10),
     +     ztab(11),ztab(12),ztab(13)
     +    /0.0d0,  1.0d-4 , 3.0d-4 , 1.0d-3 , 2.0d-3,
     +     4.0d-3, 1.0d-2 , 2.0d-2 , 3.0d-2 , 4.0d-2,
     +     6.0d-2, 8.0d-2 , 1.0d-1/
c
      data i_opinit /-1/
c
c  set flag for opacity initialization
c
      i_opinit = 1
c
c     compute log10(ztab)
      ztabl(1)=-7.0d0
      do 5 i=2,nzva
        ztabl(i)=dlog10(ztab(i))
    5 continue
c
c     tension parameter for rational splines
      pp = 0.d0
c
c     setup variables for common /arrdim/
d     nxiv=nxir
d     nyiv=nyir
d     ndatv=ndat
d     ntabv=ntab
d     nzvav=nzva
d     niwkv=niwkl
d     nwkv=nwkl
c 
c     setup iorder for common /tablex/
      iorder = iord
      if(iorder.lt.2)then             ! 040203
         iorder=2
         if(istdpr.gt.0) write(istdpr,*) ' opinit: WARNING set iord = 2'
      endif
      if(iorder.gt.6)then             !040203
         iorder=6
         if(istdpr.gt.0) write(istdpr,*) ' opinit: WARNING set iord = 6'
      endif
c---  relativemachine precision
c---  in case using linear interpolation for part. deriv.
c     drelpr = eps
c     toll   = 1.d0 - eps
c     eps10  = 10.d0 * eps
c---  in case of using analytical part.deriv.
      drelpr = 10.d0*eps
      toll   = 1.d0 - drelpr
      eps10  = drelpr
c
c     read the absolute pathnames of the
c     input files (optabe.bin,pderivs.dat,ival.dat,lowtab.dat)
c
      open(29,file=tabnam,form='formatted',status='old',err=9030)
      read(29,'(A)')optabe
      read(29,'(A)')lowtab
      read(29,'(A)')pderiv
      read(29,'(A)')ivadat
      close(29)
      if(istdpr.gt.0) then
        write(istdpr,'(/,a)') ' using tables:'
        write(istdpr,'(a)') optabe
        write(istdpr,'(a)') lowtab
        write(istdpr,'(a)') pderiv
        write(istdpr,'(a,/)') ivadat
      end if
c
c     read pointer and partial derivatives for subroutine masubi(e)
      open(33,file=pderiv,status='old',
     +        form='unformatted',err=9020)
c
c     # of different z-values (mass fraction of heavy elements)
      read(33)nzvai
      if(nzvai.gt.nzva)then
        write(istdou,*) 'opinit:error in Z-numbers'
        if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,*) 
     *    'opinit:error in Z-numbers'
        stop 'in dopinit'
      endif
c
c     # of opacity-tables
      read(33)ntabi
      if(ntabi.gt.ntab)then
        write(istdou,*) 'opinit:error in X-numbers'
        if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,*) 
     *    'opinit:error in X-numbers'
        stop 'in dopinit'
      endif
c
c
c     read # of triangles (nt), ial=2*nd+5 and triangulation-indices (iwk)
      read(33)nti
      if (abs(imode).eq.2) goto 4567    ! only rational splines
c
      if(istdpr.gt.0) write(istdpr,*) 
     *  'read partial derivative values and triangle indices...'
      if(nti.gt.nt)then
        write(istdou,*) 'error in triangle-numbers nt= ',nt,' nti= ',nti
        if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,*) 
     *    'error in triangle-numbers nt= ',nt,' nti= ',nti
        stop 'in dopinit'
      endif
c
      read(33)iali
      read(33)(iwk(n),n=1,3*nti)
c
c     read partial derivative values (wk ; n=2*nd=ial-5)
      do 2500 l=1,nzvai
       do 2301 m=1,ntabi
        read(33)(wk(n,m,l),n=1,iali-1-4)
 2301  continue
 2500 continue
c
 4567 close(33)
c
c     read ival-array from file 
      nxiri=nxir       ! extrapolated OPAL tables
      nyiri=nyir
      nyisi=nyis       ! nyis=5 : AF05<->OPAL fit at log(T)=3.90<->3.95
      nyifi=nyiri+ntal-nyisi+1
      if(istdpr.gt.0) write(istdpr,*) 
     *  'read ival.dat for extrapolation domain checking'
      open(32,file=ivadat,status='old',err=9010)
      do 3001 l=1,nzvai
         read(32,'(70i3)')((ivalo(m,k,l),m=1,nyiri),k=1,ntabi)
 3001 continue
      close(32)
c
c     read the opacity-tables
      lun =  31
      if(istdpr.gt.0) write(istdpr,*)  'read OPAL tables ....'
      call rdi95(lun,optabe,ntab,nxir,nyir,nzva,rlg,tlg,opa,ival)
      if(istdpr.gt.0) write(istdpr,*)  
     *  'read Alexander-Ferguson 2005 low-T tables ....'
c     write(istdou,*) 'AF05 log(T)-values: ',ntal
      call rdaf05(lun+2,lowtab,ntab,nxir,ntal,nzva,rlgal,tlgal,opaal)
      iyifi= nyifi                        ! total nr. of log(T) values
c     write(istdou,*) 'total log(T)-values: ',iyifi
c
c     do we use electron conduction ?
      if(imode.gt.0)then
         lec=.true.
         if(istdpr.gt.0) write(istdpr,'(a)') 
     *     ' electron conduction enabled'
      else
         lec=.false.
         if(istdpr.gt.0) write(istdpr,'(a)') 
     *     ' electron conduction disabled'
      endif
c
c     load arrays with opacity-values suitable for the s/r masub & rbival
c
      do 4010 l=1,nzvai
c       define xd- and yd- values (each of size nd=(nyir-nyis)*nxir)
        do 4001 itab=1,ntabi
          nd = 0
          do 1002 j=nyisi,nyiri
           jj = iyifi-nyiri+j
           do 1001 i=1,nxiri
               nd = nd + 1
               zd(nd,itab,l) = opa(i,j,itab,l)
               xd(nd,l) = rlg(i,itab,l)
               yd(nd,l) = tlg(j,itab,l)
               ops(i,jj,itab,l) = opa(i,j,itab,l)
 1001      continue
           tls(jj) = tlg(j,itab,l)
 1002     continue
c
c         load low-temperature opacity-values
c
            do js=1,ntal                    ! load AF05 tables
             do is=1,nxiri                  ! of size nd=27x27=729 ...OPAL
              nd=nd+1
              zd(nd,itab,l)     = opaal(is,js,itab,l)
              xd(nd,l)          = rlgal(is,itab,l)
              yd(nd,l)          = tlgal(js,itab,l)
              ops(is,ntal-js+1,itab,l) = opaal(is,js,itab,l)
             enddo
             tls(ntal-js+1)=tlgal(js,itab,l)
            enddo 
c
c      calculate birational spline coefficients ra if |imode|>1
c
        if (abs(imode).gt.1) then
          if((l.eq.1).and.(itab.eq.1).and.istdpr.gt.0)
     .      write(istdpr,*) 'compute birational spline coefficients...'
          ir=1             ! compute derivatives at boundaries
          call rat2d(nxiri,iyifi,nxir,nyif,
     .               rlg(1,itab,l),tls,ops(1,1,itab,l),
     .               ir,pp,drelpr,
     .               p,q,r,ra(1,1,1,1,itab,l),iflag,
     .               rsdx,rsax,rsbx,rscx,rsrx,
     .               rsdy,rsay,rsby,rscy,rsry)
          if(iflag.ne.0)then
             write(istdou,*) 'opinit: ERROR in s/r rat2d: iflag= ',iflag
             if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,*) 
     *         'opinit: ERROR in s/r rat2d: iflag= ',iflag
             stop 'in dopinit'
          endif
        endif
4001    continue
4010  continue      
c
      if(istdpr.gt.0) write(istdpr,'(a,i5,/)') 
     *  ' total data points nd=',nd
c
c
      return
c
9010  write(istdou,*) 'opinit: ERROR in opening ',ivadat
      if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,*) 
     *  'opinit: ERROR in opening ',ivadat
      stop 'in dopinit'
9020  write(istdou,*) 'opinit: ERROR in opening ',pderiv
      if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,*) 
     *  'opinit: ERROR in opening ',pderiv
      stop 'in dopinit'
9030  write(istdou,*) 'opinit: error in opening ',tabnam
      if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,*) 
     *  'opinit: error in opening ',tabnam
      stop 'in dopinit'
c
      end
