      program main
c
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      implicit double precision (a-h,o-z)
      include 'engenr.bz.d.incl'
c
      parameter (nwork = (iywstr+1)*nnmax, nsoon = 21*nnmax, 
     *  iy = 2*ivarmx, ivrmx1 = ivarmx+1, ivrmx4 = 4*ivarmx,
     *  nksidr = nbcprv + 5*idthm*idthm + 
     *  (nspcmx+2)*(7+nspcmx+(nspcmx+1)*(nspcmx+2)),
     *  nebcdr = (nspcmx+1)*(1+(nspcmx+2)*(nspcmx+3)),nspcm2=2*nspcmx)
c
c..      include '/usr/include/fortran/signal.h'
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
c  Modified 19/8/96, for introducing effects of centrifugal force
c
c  Corrected 5/8/97, by including engenr.nnz.d.incl, to set nspdmx etc
c
      integer sigvec
      common/ln10/ amm(3)
      common/xarr/ xxx(nnmax)
      common/yarr/ yyy(iy,nnmax)
      common/caddvr/ addvar(5,nnmax)
      common/yprtst/ yprt(ivrmx1,nnmax)
      common/cyfstr/ yzfstr(ivrmx4,nnmax)
      common/xnwvar/ xxnn(nnmax)
      common/heavy/ zatmos, zhc,zh(nnmax)
      common/compvr/ cvr(icvrmx, nnmax)
      common/convvr/ cvvarf(icvvar,6), cvvarl(icvvar,6),
     *  cvvars(icvvar,nnmax)
      common/comgrt/ isprot, irotcn, a5, riner0, riner, omgrt0,
     *  omgrot(nnmax)
c
c  used for storage for setting new oscillation variables, and in
c  s/r conmix and conmxc. 
c
      common/anwvar/ datdnn(8), aann(istrmx,nnmax)
      common/work/ www( nwork)  /ksider/ dkk(nksidr)
      common/bcprev/ dbcp(nbcprv)
      common/opccof/ iopccf,iopcdm,opcdat(40000)
      common/copcxz/ ixzopc,ixzdum,xzopc(2,20)
      common/excf/ dmm01(24)
      common/consts/ dm01(44)  /eqstd/ dm02(90)  /he3eql/ dm03(2)
      common/cnvout/ dddmmm(30)  /rnrout/ dddmm1(100)
     *  /sooner/ dddmm2(nsoon)  /convpt/ dddmm3(60)
     *  /convpp/ dddmm4(60)
      common/opctcl/ alamop,zatmop,tlmopc,dtmopc,fdgopl,iwdopc,
     *  iopacm,ifdgop
      common/bccn/ dmm03(6)   /bcatms/ dm04(7)  /opccnt/ dm06(10)
     *  /thetac/ dm07(ivarmx)  /opcxdr/ dm08(2)  /hviond/ dm09(135)
      common/cqhopf/ cqqhpf(7)
      common/crxstr/ rxstr(nspcm2,nnmax)
      common/rhcn/ dm10(13)
      common/clshft/ alshft
      common/cntmsh/ dm11(18)
      common/rbcder/ dm12(52)
      common/enggrv/ epsg(nnmax)
      common/ebcder/ dm13(nebcdr)
      common/cdgtnr/ idgtnr,idgtdm,dytnrk(ivarmx,nnmax)
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3,idgmxc
      common/cdiffu/ idiffus, itbdif, rhob, rmxdif, tlfdif, dtldif,
     *  ctdfmx, rctdif, rctdf1
c
      common/cversn/ ivrbcs,ivrrhs,ivreos,ivropc,ivreng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c..      external fpemes
c
c..      isig= sigvec(sigfpe,sig_call,fpemes,isig2,isig3)
c
      write(istdou,100) nnmax,nspcmx
      if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,100) 
     *  nnmax,nspcmx
c
c  in this version of the code, zero flag for diffusion
c
      idiffus=0
c
c  zero luminosity shift
c
      alshft=0.d0
c
c  zero flags for diagnostics
c
      idgbcs=0
      idgeng=0
      idgrhs=0
      idghe3=0
c
c  zero version numbers
c
      ivrbcs=0
      ivrrhs=0
      ivreos=0
      ivropc=0
      ivreng=0
c  reset tlmopc
      tlmopc=4.2
      dtmopc=0.2
c
      call mnevol
      stop
  100 format(//1x,70('*')//
     *  '  In this version of evolprg maximum number of mesh points is',
     *  i6//'  Maximum number of elements is ',i3//
     *  1x,70('*'))
      end
      subroutine fpemes(iarg1, iarg2)
c
c  warning message for floating point exception
c
      implicit double precision (a-h,o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.nnz.d.incl'
c
      logical norct
      common/rnratd/ al(10,krnrmx),norct
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      data icount /0/
c
      icount=icount+1
c
      if(istdpr.gt.0) write(istdpr,100) iarg1, iarg2
      if(icount.le.100) then
        return
      else
        if(istdpr.gt.0) write(istdpr,120) icount
        stop
      end if
  100 format(' ***** fpemes called with arguments ',2i10)
  110 format(//' al:',1p6e12.4)
  120 format(//' **** stopping after ',i4,' errors')
      end
      double precision function eamean(ea,idiffus)
c
c  finds average change of non-composition variables
c  Note that storage of variables depends on inclusion or not of diffusion,
c  as determined by idiffus
c
c  Original version: 15/12/92
c
      implicit double precision (a-h, o-z)
      dimension ea(1)
      dimension eas(4)
c
      do 10 i=1,3
   10 eas(i)=ea(i)
      if(idiffus.eq.0) then
	eas(4)=ea(4)
      else
	eas(4)=ea(5)
      end if
c
      eamean=amean(eas,4)
      return
      end
      double precision function amean(a,n)
c
c  finds mean of array a(i), i=1,...,n
c
      implicit double precision (a-h,o-z)
      dimension a(1)
c
      if(n.gt.0) go to 10
      amean=0
      return
c
   10 sum=0
      do 20 i=1,n
   20 sum=sum+a(i)
      amean=sum/n
      return
      end
