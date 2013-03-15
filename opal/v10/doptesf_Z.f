	program optesf_Z
	implicit double precision (a-h,o-z)
c
c	test-driver program for opacity
c       interpolation subprograms
c       opinit       (opacity initialisation)
c       opintc       (opacity interpolation 677)
c       opints       (opacity interpolation rational splines)
c
c       8.11.1993:
c               new argument imode to select algorithm:
c               if imode = 0 -> 677 only (minimum norm)  (s/r opint{c,f})
c                  imode > 0 -> 677 + birational splines (s/r opint{c,f} + opints)
c
c      13/11/95: modified for OPAL95 tables
c
        character*80 tabnam
        dimension zz(101),oplg(101)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c	get relative machine precision
	call maceps(eps)
c
        data tabnam /'OPINTPATH_02'/
c       data tabnam /'OPINTPATH_95'/
        data iorder /4/
        data imode  /-3/ !initialize minimum norm & rat. splines +EC
c
c       initialize opacity tables, use both algorithm (imode > 0)
        call opinit(eps,iorder,tabnam,imode)
c
c       define x   (hydrogen mass fraction)
c              z   (mass fraction of heavy elements)
c              tlg ( log10(temperature) )
c              rlg (rho/(t6**3))
c
        if(istdpr.gt.0) 
     *  write(istdpr,'(A)')' Exit with x < 0.'
1000    if(istdpr.gt.0) 
     *  write(istdpr,'(A,$)')' Enter x-value: '
        read(5,'(d12.5)')x
        if(x.lt.0.d0) goto 2001
        if(istdpr.gt.0) 
     *  write(istdpr,'(A,$)')' Enter tlg: '
        read(5,'(d12.5)')tlg
        if(istdpr.gt.0) 
     *  write(istdpr,'(A,$)')' Enter dlg: '
        read(5,'(d12.5)')dlg
        if(istdpr.gt.0) 
     *  write(istdpr,'(A,$)')' Enter first z-value: '
        read(5,'(d12.5)')z1
        if(istdpr.gt.0) 
     *  write(istdpr,'(A,$)')' Enter last z-value: '
        read(5,'(d12.5)')z2
c
        print *,''
        print *,'using opintf'
c
c       rlg (rho/(t6**3))
c       dlg ( log10(rho) )
c
        rlg = dlg - 3.d0*tlg +18.d0
c
c       677-interpolate to get opacity value opalg
        dz=(z2-z1)/1.0d2
        z=z1
        do i=1,101
          iexp=0
          ier=0
          call opintf(x,z,tlg,rlg,opalg,opr,opt,opx,opz,iexp,ier)
          if(ier.gt.0)if(istdpr.gt.0) 
     *  write(istdpr,'(a,i3)') 
     +               'optesf: error in interpolation s/r, ier=',ier
          if(iexp.gt.0)print *,'extrapolated points: ',iexp
c
c         print output
c         print *,'------------ 677 algorithm ---------'
c         print *,'log10 Opacity = ',opalg
c         print *,'      Opacity = ',10.d0**opalg
c         print *,'part. der. R  = ',opr
c         print *,'part. der. T  = ',opt
c         print *,'part. der. X  = ',opx
c         print *,'part. der. Z  = ',opz
c         print *
          print *,z,opalg
          zz(i)=z
          oplg(i)=opalg
          z=z+dz
        enddo
        write(13) (zz(i),i=1,101),(oplg(i),i=1,101)
c
        goto 2001
c
c	interpolate to get opacity value opalg
 	iexp=0
        ier=0
	call opints(x,z,tlg,rlg,opalg,opr,opt,opx,opz,iexp,ier)
        if(ier.gt.0)if(istdpr.gt.0) 
     *  write(istdpr,'(a,i3)') 
     +               'optesf: error in interpolation s/r, ier=',ier
        if(iexp.gt.0)print *,'extrapolated points: ',iexp
c
c	print output
        print *,'----------rational splines ---------'
	print *,'log10 Opacity = ',opalg
	print *,'      Opacity = ',10.d0**opalg
        print *,'part. der. R  = ',opr
        print *,'part. der. T  = ',opt
        print *,'part. der. X  = ',opx
        print *,'part. der. Z  = ',opz
        goto 1000
2001    continue
c
	end
