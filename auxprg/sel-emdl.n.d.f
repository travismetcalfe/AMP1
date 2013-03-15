      program main
c
c  select single or several model(s) from file of 
c  unformatted evolution models
c
c  Double precision version
c  ++++++++++++++++++++++++
c
c  Dated 13/3/90
c
c  Modified 11/9/92 to allow for both old and new model format.
c
      implicit double precision (a-h, o-z)
      include 'evolpr.incl'
c..      parameter (nnmax = 2411,ivarmx=12,nrdtmx=100,nidtmx=50,
c..     *  nbccmx = 100)
      character*60 fin,fout
      dimension x(nnmax),y(iymax,nnmax),datmod(nrdtmx),
     *  ndtmod(nidtmx),bccoef(nbccmx)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      istdpr = 6
c
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'Enter input file name'
      read(5,'(a)') fin
c
      open(2,file=fin,status='old',form='unformatted')
c
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'Enter output file name'
      read(5,'(a)') fout
c
      open(3,file=fout,status='unknown',form='unformatted')
c
      nmod=0
      nmodw=0
c
    5 if(istdpr.gt.0) 
     *  write(istdpr,*) 'Enter model number'
      read(5,*,end=90,err=90) nmout
c
      if(istdpr.gt.0) 
     *  write(istdpr,100) fin
      icase = 0
c
   10 call rdemdl(2,x,y,datmod,ndtmod,bccoef,iymax,iform,
     *  nn,nrdtmd,nidtmd,nbccf,nvar,icase,icry)
      if(icry.lt.0) go to 30
      icase = -1
c
      nmod=nmod+1
      if(nmod.lt.nmout) go to 10
c
   30 call wdtmod(6,datmod,nrdtmd,ndtmod,nidtmd,0)
c
      if(nmod.gt.nmodw) then
        call wremdl(3,x,y,datmod,ndtmod,bccoef,iymax,iform,
     *    nn,nrdtmd,nidtmd,nbccf,nvar)
        nmodw=nmod
      end if
c
c  test for reading further models
c
      go to 5
c
   90 continue
      stop
  100 format(' file name: ',a60)
      end
