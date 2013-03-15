        subroutine rdaf05(lun,optabe,ntab,nval,nlin,nzva,rlg,
     +                    tlg,opa)
        implicit double precision (a-h,o-z)
        character*(*) optabe
c
c       purpose:
c               read extrapolated, binary ALexander-Furgenson (2005) tables
c
c       written: g.w. houdek
c
c       History: 13.6.1992 creation
c
c                19.12.1992 modified for z-domain interpolation
c                           nzva different z-values
c                           (z = mass fraction of heavy elements)
c
c                           modified opacity values input from
c                           'formatted' to 'unformatted'
c                           first you have to run the program 'a2b'
c                           (ascii-to-binary conversion of opacity-tables)
c
c       25.1.1993 changed variable t6 -> tlg = log10(t)
c                 table-values changed to tlg in da2b
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
c                setenv OPINTPATH /dir1/dir2/opacity.pathes
c
c                the character variable opinte is defined new
c                in the list of the subroutine arguments
c
c      13/11/95: modified for OPAL95 tables
c      05/10/02: modified for OPAL02 tables (Frank Pipers)
c      23/02/06: modified for AGS04 composition
c
c      last modification: 23/02/06
c
c       input values:
c       lun ....... logical unit number of opacity-table-file
c       optabe .... absolute filename path of the opacity table
c       ntab ...... nr. of different x-values [ntab=8]
c       nval ...... max. nr. of op.values per line (const. tlg)
c       nlin ...... max. nr. of op.table lines
c       nzva ...... max. nr. of different z-values [normaly z=13]
c                   (z = mass fraction of heavy elements)
c
c       output values:
c       rlg ....... array[1..nval,1..ntab,1..nzva], 
c                   decade log of r=density(gm/cm**3)/t6**3
c       tlg ....... array[1..nlin,1..ntab,1..nzva].
c                   temperature/10**6
c       opa ....... array[1..nvar,1..nlin,1..ntab,1..nzva],
c                   opacity values
c
        integer         lun,ntab,nval,nlin,nzva
        dimension       rlg(nval,ntab,nzva)
        dimension       tlg(nlin,ntab,nzva)
        dimension       opa(nval,nlin,ntab,nzva)
        dimension       ival(nlin,ntab,nzva)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c         open inputfile
          open(lun,file=optabe,form='unformatted',
     +             status='old',err=9001)
c
c       read ntabi tables
        do 3001 k=1,ntab
c
c         read nzva-times ntab tables
          do 1015 l=1,nzva
c
c          read line of rlg-values
           read(lun,err=9002) (rlg(i,k,l),i=1,nval)
c
c          read 'tlg' + opacity-values
           do 1010 j=1,nlin
            read(lun,end=1010,err=9002)
     +           tlg(j,k,l),(opa(i,j,k,l),i=1,nval)
c
            do 1005 m=nval,1,-1
               if(opa(m,j,k,l).ne.0.0d0)goto 1007
 1005       continue
 1007       ival(j,k,l) = m
 1010      continue
 1015     continue
 3001   continue
c       close input file lun
        close(lun)
c
        return
c
 9001   write(istdou,*) 'rdaf05: error in open statement'
        if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,*) 
     *    'rdaf05: error in open statement'
        stop 'in rdaf05'
 9002   write(istdou,*) 'rdaf05: error in read: line: ',j,k,l,i,nval
        if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,*) 
     *    'rdaf05: error in read: line: ',j,k,l,i,nval
        stop 'in rdaf05'
c
        end
