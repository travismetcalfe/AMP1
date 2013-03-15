      program main
c
c  scan file of evolution models using scnfil
c
c  Modified 14/8/91, to accomodate changed model format.
c
      implicit double precision (a-h, o-z)
      character*60 fin
      include 'evolpr.incl'
c..      parameter (nnmax = 4811,ivarmx=12,nrdtmx=100,nidtmx=50,
c..     *  nbccmx = 100)
      dimension x(nnmax),y(iymax,nnmax),datmod(nrdtmx),
     *  ndtmod(nidtmx),bc(nbccmx)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      istdpr = 6
c
      write(0,*) 'Enter input file'
      read(5,'(a)') fin
      open(2,file=fin,status='old',form='unformatted')
c
      icase = 0
c
      go to (1,2), ic
    1 write(istdpr,100) fin
      go to 5
    2 write(istdpr,102) fin
c
    5 nmod=0
c
   10 call rdemdl(2,x,y,datmod,ndtmod,bc,iymax,iform,
     *  nn,nrdtmd,nidtmd,nbccf,nvar,icase,icry)
      if(icry.lt.0) go to 90
      icase = -1
      nmod=nmod+1
      tef=10.**y(3,1)
      ams=datmod(23)/1.989e33
c
c  test for case
c
      go to (20,30), ic
c
   20 pc=bc(3)*1.e17
      tc=bc(4)*1.e7
      xc=bc(5)
      write(istdpr,110) nmod,ams,datmod(2),datmod(1),datmod(22),
     *  datmod(24),tef,datmod(25),pc,tc,xc,ndtmod(1)
      go to 10
c  test for setting of patmos
c
c  note that as a result of resetting  of datmod (23/11/82), datmod(21)
c  is non-zero and patmos is in datmod(20) instead of in datmod(19)
c
   30 iptmos=19
      if(datmod(21).ne.0) iptmos=20
      patmos=datmod(iptmos)
      write(istdpr,120) nmod,ams,datmod(22),datmod(2),datmod(1),
     *  datmod(3),patmos,datmod(24),datmod(25),ndtmod(1)
      go to 10
   90 continue
      ids=2
      if(nmod.eq.0) write(istdpr,190) ids
      stop
  100 format('# Evolution sequence on file ',a/'#'/
     *  '# nmod,mass,surface x,z,age,radius,',
     *  'eff. temperature,luminosity,pc,tc,xc,icase:'/'#')
  102 format('# Evolution sequence on file ',a/'#'/
     *  '# nmod,mass,age,surface x,z,alfa,patmos,radius,',
     *  'luminosity,icase:'/'#')
  110 format(i4,f11.6,2f10.6,1p6e15.8,0pf11.8,i12)
  120 format(i4,f12.6,1pe13.5,0p3f10.6,1p3e15.7,i12)
  190 format(//1x,10('*'),' no models on d/s',i3)
      end
