	program optesf
	implicit double precision (a-h,o-z)
      common/cstdio/ istdin, istdou, istdpr, istder
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
c      13/11/95: modified fort OPAL95 tables
c
        common /opatab/ioptab
        character*80 tabnam
c
c	get relative machine precision
	call maceps(eps)
c
        data tabnam /'OPINTPATH_AX'/
        data iorder /4/
        data imode  /3/ !initialize minimum norm & rat. splines +EC
c
c	initialize opacity tables, use both algorithm (imode > 0)
	call opinit(eps,iorder,tabnam,imode)
        write(istdpr,*)''
        write(istdpr,*)'TABLES: ',ioptab
        write(istdpr,*)''
c
c	define x   (hydrogen mass fraction)
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
     *  write(istdpr,'(A,$)')' Enter z-value: '
        read(5,'(d12.5)')z
        if(istdpr.gt.0) 
     *  write(istdpr,'(A,$)')' Enter tlg: '
        read(5,'(d12.5)')tlg
        if(istdpr.gt.0) 
     *  write(istdpr,'(A,$)')' Enter rlg: '
        read(5,'(d12.5)')rlg
c
c       rlg (rho/(t6**3))
c       dlg ( log10(rho) )
c
c       rlg = dlg - 3.d0*tlg +18.d0
c
c	677-interpolate to get opacity value opalg
 	iexp=0
        ier=0
	call opintf(x,z,tlg,rlg,opalg,opr,opt,opx,opz,iexp,ier)
        if(ier.gt.0)if(istdpr.gt.0) 
     *  write(istdpr,'(a,i3)') 
     +               'optesf: error in interpolation s/r, ier=',ier
        if(iexp.gt.0)write(istdpr,*)'extrapolated points: ',iexp
c
c	print output
        write(istdpr,*)'------------ 677 algorithm ---------'
	write(istdpr,*)'log10 Opacity = ',opalg
	write(istdpr,*)'      Opacity = ',10.d0**opalg
        write(istdpr,*)'part. der. R  = ',opr
        write(istdpr,*)'part. der. T  = ',opt
        write(istdpr,*)'part. der. X  = ',opx
        write(istdpr,*)'part. der. Z  = ',opz
        print *
c
c	interpolate to get opacity value opalg
 	iexp=0
        ier=0
	call opints(x,z,tlg,rlg,opalg,opr,opt,opx,opz,iexp,ier)
        if(ier.gt.0)if(istdpr.gt.0) 
     *  write(istdpr,'(a,i3)') 
     +               'optesf: error in interpolation s/r, ier=',ier
        if(iexp.gt.0)write(istdpr,*)'extrapolated points: ',iexp
c
c	print output
        write(istdpr,*)'----------rational splines ---------'
	write(istdpr,*)'log10 Opacity = ',opalg
	write(istdpr,*)'      Opacity = ',10.d0**opalg
        write(istdpr,*)'part. der. R  = ',opr
        write(istdpr,*)'part. der. T  = ',opt
        write(istdpr,*)'part. der. X  = ',opx
        write(istdpr,*)'part. der. Z  = ',opz
        goto 1000
2001    continue
c
	end
