      program main
c
c  print data from unformatted evolution model
c
c  ++++++++++++++++++++++++
c
c  Dated 13/3/90
c
c  Modified 13/8/91 to allow for both old and new model format.
c
      implicit double precision (a-h, o-z)
      include 'evolpr.incl'
c..      parameter (nnmax = 2411,ivarmx=12,nrdtmx=100,nidtmx=50,
c..     *  nbccmx = 100)
      character*60 fin
      dimension x(nnmax),y(iymax,nnmax),datmod(nrdtmx),
     *  ndtmod(nidtmx),bccoef(nbccmx)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      istdpr = 6
c
      write(istdpr,*) 'Enter file name'
      read(5,'(a)') fin
c
      ibcout=0
      write(istdpr,*) 'Enter 1 for bccoef output'
      write(istdpr,*) 'Enter 2 for long radius and luminosity output'
      read(5,*) ibcout
c
      open(2,file=fin,status='old',form='unformatted')
      write(istdpr,100) fin
c
      nr=0
      icase = 0
c
   10 call rdemdl(2,x,y,datmod,ndtmod,bccoef,iymax,iform,
     *  nn,nrdtmd,nidtmd,nbccf,nvar,icase,icry)
      if(icry.lt.0) go to 90
      icase = -1
      nr=nr+1
      if(ibcout.ne.2) then
	if(mod(nrdtmd,5).eq.0) then
	  krec = nrdtmd/5
        else
	  krec = nrdtmd/5+1
	end if
	write(istdpr,135) nr,nn
	do 30 k=1,krec
	i1=1+5*(k-1)
	i2=min0(i1+4,nrdtmd)
   30   write(istdpr,140) i1,(datmod(i),i=i1,i2)
c
	if(mod(nidtmd,5).eq.0) then
	  krec = nidtmd/5
        else
	  krec = nidtmd/5+1
	end if
	write(istdpr,145) 
	do 35 k=1,krec
	i1=1+5*(k-1)
	i2=min0(i1+4,nidtmd)
   35   write(istdpr,150) i1,(ndtmod(i),i=i1,i2)
        if(ibcout.eq.1) then
	  if(mod(nbccf,5).eq.0) then
	    krec = nbccf/5
          else
	    krec = nbccf/5+1
	  end if
          write(istdpr,160)
	  do 40 k=1,krec
	  i1=1+5*(k-1)
	  i2=min0(i1+4,nbccf)
   40     write(istdpr,140) i1,(bccoef(i),i=i1,i2)
        end if
c
c  output x0 and alfa with higher precision
c
        write(istdpr,175) datmod(2), datmod(3), datmod(14)
      else
        if(nr.eq.1) write(istdpr,180) nn,datmod(2), datmod(3)
        write(istdpr,185) nr,datmod(22),datmod(24),datmod(25)
      end if
      go to 10
c
   90 continue
      if(ibcout.eq.2.and.istdpr.gt.0) 
     *  write(istdpr,150) nn,datmod(2), datmod(3)
      stop
  100 format(' file name: ',a)
  135 format(/' record no',i4,'  number of mesh points =',i6,
     *   '  datmod:')
  140 format(i3,':',1p5e13.5)
  145 format(/' ndtmod:')
  150 format(i3,':',5i13)
  160 format(/' bccoef:')
  175 format(' X0 =',f15.11,'    alfa =',f15.11,'   sbcfct =',f15.11)
  180 format(/' nn =',i5,' X0 =',f15.11,'  alfa =',f15.11)
  185 format(' nr =',i4,' age =',1pe13.6,' Rs =',e15.8,' Ls =',e15.8)
      end
