      function fdrad(x,y,zh,ak,akr,akt,akx)
c
c  sets radiative gradient, by calling eqstf and opac
c  returns  ak = log10(kappa), akr = (dlog kappa/dlog rho)t,x,
c  akt = (dlog kappa/dlog t)rho,x and akx = (dlog kappa/dlog x)t,rho
c  in argument list.
c
c  also returns, in common /opcxdr/, akxa = (dlog kappa/d x) and,
c  when iwdopc = 1, 8 or 9, akz = (dlog kappa/dlog z)
c
c  Modified 5/6/03, including erroneously forgotten common/clshft/
c  in fdradp.
c
      implicit double precision (a-h,o-z)
      logical nosd,notd
c
c  include statement setting, inter alia, scratch directory
c
      include 'engenr.bz.d.incl'
      dimension y(*)
      common/ln10/ amm
      common/rhcn/ a1,a2,a3,a4,z,nvar,ifwrt,irhtst
      common/eqstd/ xii1(4),ane(10),rho(20),ht(20),p(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp
      common/noiter/ iter1, ntime
      common/clshft/ alshft
      common/cnvout/ ddrad,rrcnv,aacnv,ddacad,amach,
     *  drr(5), ddr(5), da(5),pturb,ycnv,dyda, dtxh
      common/cdiffu/ idiffus, itbdif, rhob, rmxdif, tlfdif, dtldif,
     *  ctdfmx, rctdif, rctdf1, idfxbc, ismdif, nsmdif, idiffc1
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3,idgmxc
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      nosd=.true.
      notd=.true.
c
      xh=y(5)
      yh=1-xh-zh
      fl=y(2)
      tl=y(3)
      call eqstf(fl,tl,xh,yh,zh,nosd,notd)
c
      pl=log10(p(1))
      rhl=log10(rho(1))
c
      call opact(rhl,tl,xh,zh,akl,akr,akt,akx)
      ak=10.d0**akl
c
c  set radiative gradient
c
      drad=(a3/a2)*10.d0**(akl-x+pl-4*tl)*(10.d0**y(4)-alshft)
      fdrad=drad
      if(drad.le.-1.and.istdpr.gt.0) then
        write(istdpr,*) 'a3, a2, akl, x, pl, tl, y(4), alshft in fdrad:'
        write(istdpr,*) a3, a2, akl, x, pl, tl, y(4), alshft
      end if
      ddrad=drad-dad(1)
      n=0
c
      return
  100 format(' fdrad called with x, fl, tl, drad ='/
     *  1pe13.5,0p2f10.5,1pe13.5)
  110 format('# Output from fdrad'/'#'/
     *  '# ntime, iter1, n, q, log rho, log T, X, log kappa,',
     *  ' ddrad, log L'/'#')
      end
      function fdradp(x,y,pl,zh,fl,ak,akr,akt,akx)
c
c  sets radiative gradient, by calling eqstf and opac
c  This version (unlike fdrad) evaluates gradient at
c  fixed log p and log T, iterating for equation of state
c  fl returns the resulting log f
c  ak, akr, akt, akx returns value of opacity and its derivatives
c
c  Modified 8/11/05, adding fl as returned quantity.
c
      implicit double precision (a-h,o-z)
      logical nosd,notd
c
c  include statement setting, inter alia, scratch directory
c
c  Note: engenr.n.d.incl replaced by engenr.nnz.d.incl, 5/6/02
      include 'engenr.nnz.d.incl'
      dimension y(*)
      common/ln10/ amm
      common/rhcn/ a1,a2,a3,a4,z,nvar,ifwrt,irhtst
      common/cxhder/ xhder(2,nnmax)
      common/eqstd/ xii1(4),ane(10),rho(20),ht(20),p(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp
      common/noiter/ iter1, ntime
      common/clshft/ alshft
      common/cnvout/ ddrad,rrcnv,aacnv,ddacad,amach,
     *  drr(5), ddr(5), da(5),pturb,ycnv,dyda, dtxh
      common/cdiffu/ idiffus, itbdif, rhob, rmxdif, tlfdif, dtldif,
     *  ctdfmx, rctdif, rctdf1, idfxbc, ismdif, nsmdif, idiffc1
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3,idgmxc
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      nosd=.true.
      notd=.true.
c
      xh=y(5)
      yh=1-xh-zh
      fl=y(2)
      tl=y(3)
      call eqstp(pl,tl,xh,yh,zh,nosd,notd,fl,nit)
c
      pl=log10(p(1))
      rhl=log10(rho(1))
c
      call opact(rhl,tl,xh,zh,akl,akr,akt,akx)
      ak=10.d0**akl
c
c  set radiative gradient
c
      drad=(a3/a2)*10.d0**(akl-x+pl-4*tl)*(10.d0**y(4)-alshft)
      if(istdpr.gt.0) write(istdpr,'(5f10.6,1p4e13.5,0pf10.6)')
     *  tl,pl,fl,xh,zh,10.d0**akl,10.d0**(-x+pl-4*tl),
     *  10.d0**y(4),alshft,drad
      fdradp=drad
      if(drad.le.-1.and.istdpr.gt.0) then
	write(istdpr,*) 
     *    'a3, a2, akl, x, pl, tl, y(4), alshft in fdradp:'
        write(istdpr,*) a3, a2, akl, x, pl, tl, y(4), alshft
      end if
      ddrad=drad-dad(1)
      n=0
c
      return
  100 format(' fdradp called with x, fl, tl, pl, drad ='/
     *  1pe13.5,0p3f10.5,1pe13.5)
  110 format('# Output from fdradp'/'#'/
     *  '# ntime, iter1, n, q, log rho, log T, X, log kappa,',
     *  ' ddrad, log L, log p'/'#')
      end
      function fdtxh(x,y,xh,dxh)
c
c  sets composition-gradient contribution to superadiabatic
c  temperature gradient in point characterized by x, y(.).
c  On input xh = X and dxh = d X/d x at that point.c
c  It is assumed that EOS quantities are already set in common/eqstd/
c  If idiffc1 = 0, gradient is set to zero.
c
c  Original version: 23/10/03
c
      implicit double precision(a-h, o-z)
      dimension y(*)
      common/ln10/ amm
      common/rhcn/ a1,a2,a3,a4,z,nvar,ifwrt,irhtst
      common/eqstd/ xii1(4),ane(10),rho(20),ht(20),p(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp
      common/cnvout/ ddrad,rrcnv,aacnv,ddacad,amach,
     *  drr(5), ddr(5), da(5),pturb,ycnv,dyda, dtxh
      common/cdiffu/ idiffus, itbdif, rhob, rmxdif, tlfdif, dtldif,
     *  ctdfmx, rctdif, rctdf1, idfxbc, ismdif, nsmdif, idiffc1
c
      pl=log10(p(1))
      if(idiffc1.gt.0) then
	dlpdx=-a2*10.d0**(2.d0*x-4.d0*y(1)-pl)
	dtxh=-trhp*rhxp*dxh/(amm*xh*dlpdx)
      else
	dtxh=0.d0
      end if
      fdtxh=dtxh
      return
      end
