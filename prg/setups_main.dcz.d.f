      subroutine setups_main
c
c  Set up storage etc. for evolution code. This used to
c
c  be the main programme for the code.
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      implicit double precision (a-h,o-z)
      include 'engenr.bz.d.incl'
c
      parameter (nwork = (iywstr+1)*nnmax, nsoon = 50*nnmax, 
     *  nwork1 = 100*nnmax,
     *  iy = 2*(ivarmx+nspcmx)+1, ivrmx1 = ivarmx+1, ivrmx4 = 4*ivarmx,
     *  nksidr = nbcprv + 5*idthm*idthm + 
     *  (nspcmx+2)*(7+nspcmx+(nspcmx+1)*(nspcmx+2)),
     *  nebcdr = (nspcmx+1)*(1+(nspcmx+2)*(nspcmx+3)),nspcm2=2*nspcmx,
     *  nsdifr=2*nspdmx+2)
c
c  storage parameters for internal storage of evolution sequences
c  As set now, nstore_emdl is sufficient for, for example, 500 timesteps, 
c  601 mesh points and 8 variables per meshpoint
c
      parameter(nstore_emdl_p=2 800 000, 
     *  nstore_max=1000, nstore_datmod_p = nrdtmd*nstore_max,
     *  nstore_ndtmod_p = nidtmd*nstore_max, 
     *  nstore_bc_p = nbcprv*nstore_max)
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
c  Modified 14/10/04, to prepare for iteration using repeated fortran
c  calls. Note that mnevol now has argument i_paramset.
c
      integer sigvec
      dimension www1(nwork1)
      common/ln10/ amm(3)
      common/xarr/ xxx(nnmax)
      common/yarr/ yyy(iy,nnmax)
      common/xdarr/ xxxd(nnmax)
      common/ydarr/ yyyd(iy,nnmax)
      common/caddvr/ addvar(5,nnmax)
      common/cgngvr/ yvar(igvrmx,nnmax)
      common/yprtst/ yprt(ivrmx1,nnmax)
      common/cyfstr/ yzfstr(ivrmx4,nnmax)
      common/xnwvar/ xxnn(nnmax)
      common/heavy/ zatmos, zhc,zh(nnmax)
      common/compvr/ cvr(icvrmx, nnmax)
      common/convvr/ cvvarf(icvvar,6), cvvarl(icvvar,6),
     *  cvvars(icvvar,nnmax)
      common/comgrt/ isprot, irotcn, a5, riner0, riner, omgrt0,
     *  omgrot(nnmax)
      common/comgrp/ isprtp, irotcp, omgrtp(10001)
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
      common/ebcder/ dm13(nebcdr)
      common/enggrv/ epsg(nnmax)
      common/cdgtnr/ idgtnr,idgtdm,dytnrk(ivarmx,nnmax)
      common/cdifsv/ csdifr(nsdifr,nnmax)
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3,idgmxc
      common/cdiffu/ idiffus, itbdif, rhob, rmxdif, tlfdif, dtldif,
     *  ctdfmx, rctdif, rctdf1
c
c  commons for internal storage of model variables (total sizes set
c  in parameter statement above)
c
      common/cstore_emdl/ nstore_emdl, istore_emdl, 
     *  store_emdl(nstore_emdl_p)
      common/cstore_datmod/ nstore_datmod, istore_datmod, 
     *  store_datmod(nstore_datmod_p)
      common/cstore_ndtmod/ nstore_ndtmod, istore_ndtmod, 
     *  lstore_ndtmod(nstore_ndtmod_p)
      common/cstore_bc/ nstore_bc, istore_bc,store_bc(nstore_bc_p)
c
      common/cversn/ ivrbcs,ivrrhs,ivreos,ivropc,ivreng
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      equivalence(www(1),www1(1))
c
c..      external fpemes
c
c..      isig= sigvec(sigfpe,sig_call,fpemes,isig2,isig3)
c
      nstore_emdl   = nstore_emdl_p
      nstore_datmod = nstore_datmod_p
      nstore_ndtmod = nstore_ndtmod_p
      nstore_bc     = nstore_bc_p
      write(istdou,100) nnmax,nspcmx,ivarmx
      if(istdou.ne.istdpr.and.istdpr.gt.0)
     *  write(istdpr,100) nnmax,nspcmx,ivarmx
      write(istdpr,'(/a,i5)') ' Maximum size of ksider: ', nksidr
      write(istdou,110) nstore_emdl, nstore_datmod, nstore_ndtmod,
     *  nstore_bc
      if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,110) 
     *  nstore_emdl, nstore_datmod, nstore_ndatmd,
     *  nstore_ksidr
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
      return
  100 format(//1x,70('*')//
     *  '  In this version of evolprg maximum number of mesh points is',
     *  i6//'  Maximum number of elements is ',i3//
     *  '  Maximum number of variables is ',i3//
     *  1x,70('*'))
  110 format(/' Storage assigned for internal storage of models:'/
     *  ' nstore_emdl:  ', i8/
     *  ' nstore_datmod:', i8/
     *  ' nstore_ndtmod:', i8/
     *  ' nstore_bc:    ', i8/)
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
c..      write(6,110) (al(1,k),k=1,6)
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
c  If idiffus is entered as -1, simply take mean over all variables
c  If idiffus is entered as -2, take mean over all variables except
c  non-diffusing composition (e.g., assuming that diffusing composition
c  in convective core is treated by tnrkt).
c  If idiffus is entered as -3, take mean over all variables except
c  settling velocity and non-diffusing composition
c
c  Original version: 15/12/92
c
c  Modified 20/8/95 to allow for diffusion of several elements
c
c  Modified 3/10/03 to allow mean over all variables.
c
c  Modified 17/10/03 to allow mean including diffusing composition
c
      implicit double precision (a-h, o-z)
      dimension ea(1)
      dimension eas(20)
      common/rhcn/ a1,a2,a3,a4,zdummy,nvar,ifwrt,irhtst,inentc,
     *  ii1,ii2,ii3,icomp
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec, idcomp, iccomp
c
      if(idiffus.ge.0) then
c
        do 10 i=1,3
   10   eas(i)=ea(i)
        if(idiffus.eq.0) then
	  eas(4)=ea(4)
        else
	  eas(4)=ea(4+idcomp)
        end if
c
        eamean=amean(eas,4)
c
      else if(idiffus.eq.-1) then
c
        eamean=amean(ea,nvar)
c
      else if(idiffus.eq.-2) then
c
        eamean=amean(ea,4+2*idcomp)
c
      else 
c
        eamean=amean(ea,4+idcomp)
c
      end if
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
