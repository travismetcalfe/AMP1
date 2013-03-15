      subroutine tnrkt(x,y,zk,yp,zkp,ap,aq,rhs,bc,ii1,ii2,ii3,kk,
     .  ka,kb,ki,nn,id,ucy,ea,det,v,theta,dt,iter)
c     subroutine dtnrkt(same argument list)
c
c
c           time dependent newton-raphson-kantorovich programme
c           ***************************************************
c
c     if an error is detected, v(1) is set to zero before return
c
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
c  Modified 11/7/96, setting dimension parameters in parameter
c  statement.
c  ***** Note: the values still need checking for consistency
c
c  Modified 5/3/03, introducing common/cdgrhb/ kdgrhb, to allow
c  forcing of return if errors are encountered in s/r rhs or bc.
c
c  Modified 5/3/03, removing double-precision entry (which was not
c  effective anyway).
c
      implicit double precision (a-h,o-z)
      parameter(iid=30, ikd=35, igd=40, ihd=8, ialam=15)
      parameter(if1=2*iid, if2=3*iid)
      integer v(1)
      dimension f(if1),fd(iid,if2),alam(ialam,iid),
     .  alamd(ialam,iid,if1),g(ikd),
     .  gd(ikd,igd),gp(ikd,igd,2),gs(ikd,igd,2),
     .  a(iid,iid),d(iid,igd),h(ihd),
     .  hd(ihd,igd),q(iid,igd),p(1),theta(1),u(iid,iid),
     .  z(if1),dzdy(iid,if1),zp(if1)
      common /work/ q,f,h,fd,hd,g,err,gp,a,d,p
      common/ctdist/ dx,iik,in,in1,ir,k1,ik1,ik13,ii23,ii,i1,ik,
     .  iii1,iik1,iiik,ii11,iw
c
c  common for diagnostic from s/r rhs and bcs. If kdgrhb is set 
c  to -1, should exit from s/r tnrkt.
c
c  This is actually controlled in s/r ctdt75.
c
      common/cdgrhb/ kdgrhb
c
c  common for diagnostic output. if idiag .ge. 1, dyst contains
c  the corrections applied to the solution.
c
      common/cdgtnr/ idiag, idgdum, dyst(1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c     note that a, d and p must be stored consecutively and in order
c     err is placed after g to provide space for gst(ii+kk+1) if
c     necessary
c
c     note that the first dimensions of d,fd,hd must all be equal
c
c     note that with present dimensions the space occupied by variables
c     preceding p in common /work/ is 3046*4 bytes
c     dimension requirement for p is
c     p((ii1+ii2)*(ii1+ii2+kk-ka+1)*nn+ii3*(ii1+ii2+kk+1)*nn)
c     minimum space occupied by common /work/ is 3822*4 bytes
c
      dimension x(nn),y(id,nn),zk(kk),ea(id,3),ap(1),
     .          aq(1),yp(id,nn),zkp(kk)
      dimension dy(1,1),dyp(1,1),dzk(1),dzkp(1)
c
c  note that if double precision entry is enabled, this
c  dimensioning statement must be replaced by
c
c      dimension dy(id,nn),dyp(id,nn),dzk(kk),dzkp(kk)
c
c     logical equivalence
      equivalence (gd(1,1),gp(1,1,1))
c     space saving equivalence
      dimension gst(1),eb(1)
      equivalence (g(1),gst(1),eb(1)),(gs(1,1,1),p(1))
c
      data kdgrhb /0/
c
      external rhs,bc
c
      save
c
c     dimension information for leq, bc and rhs
c..    1 iid=12
c..      ikd=14
c..      igd=15
c..      ihd=8
c..      ialam=4
c     line printer dsrn
    1 iw=istdou
c
c     set counting limits
    2 ii=ii1+ii2
      ii23=ii2+ii3
      iii=ii+ii3
      iik=iii+kk
      ik=ii+kk
      kab=ka+kb
      kap=ik-ki
      ika=ii-ka
      ikka=ik-ka
      kp=kap-kab
      kip=ki+kp
c
    3 n1=nn-1
      i1=ii+1
      k1=kk+1
      kab1=kab+1
      kap1=kap+1
      ik1=ik+1
      ka1=ka+1
      ika1=ika+1
      ikka1=ikka+1
      iii1=iii+1
      iik1=iik+1
      ii11=ii1+1
      ik13=ik1*ii3
c
c     set counting indices for p and r (in ctdt75)
      ip=ikka1*ii
      ipu=ip*nn
      ipv=ipu+1-ik13
      ir=ii23+iii+ii2*iii
c
c
c     compatibility test
      if(ka.le.ii.and.kp.ge.0.and.ka.ge.0.and.kb.ge.0.and.ki.ge.0.and.
     .   kk.ge.0) go to 10
      write(iw,1000) ii,kk,ka,kb,ki
      v(1)=0
      return
   10 continue
c
c
c     set v if necessary
   20 iv=0
      do 21 i=1,ii
      if(v(i).gt.0) go to 21
      iv=1
      v(i)=i
      write(iw,1100) i,i
   21 continue
      if(ii3.eq.0) go to 24
      do 23 i=i1,iii
   23 v(i)=i
   24 if(iv.eq.1) write(iw,1101) (v(i), i=1,iii)
c
c     set range of independent variable
   25 r=x(nn)-x(1)
      if(r.eq.0.) go to 407
c
c
c     **********
   30 call zero(dzdy,288)
      ia=iid
      ib=ihd
      ic=ialam
      i=nn
      iter1=1
c     conditions solely at first boundary
c
c     empty derivative matrices for boundary conditions
   32 do 33 j=1,iik
      do 33 k=1,2
      do 33 i=1,kap
   33 gs(i,j,k)=0.d0
c
c
c     boundary conditions and equations at first boundary
   35 ia=ikd
      ib=igd
      in=2
      in1=1
      n=0
c     initial dx
      dx=(x(2)-x(1))/2.d0
      call ctdt75(x(1),y(1,1),zk,yp(1,1),zkp,z,zp,dzdy,theta,dt,ap,aq,
     .  rhs,f,fd,alam,alamd,a,d,h,hd,p(ipv),u,v,id,iid,ihd,ii1,ii2,
     .  ii3,iii,kk,ki,n,iter)
      call bc(x(1),x(nn),y(1,1),y(1,nn),zk,ap,aq,g,gs,ia,ib,nn)
      if(kdgrhb.eq.-1) go to 410
   40 if(v(1).eq.0) return
c     set boundary derivative matrix
c..      write(6,*) 'After statement 40'
c..      call wrleq(gd,gd(1,ka1),ka,ikka1,ikd,ikd,err)
      do 44 i=1,iik
      iv=v(i)
      if(i.gt.iii) iv=i
      do 41 j=1,kap
      do 41 k=1,2
   41 gp(j,i,k)=gs(j,iv,k)
c     (ensure that first and second b.c. matrices are loaded into gd)
      if(kab.eq.0) go to 43
      do 42 j=1,kab
   42 gd(j,i)=gd(j,i)+gp(j,i,2)
   43 if(i.gt.kap) go to 44
      gd(i,iik1)=g(i)
   44 continue
c..      write(6,*) 'In loop 44'
c..      call wrleq(gd,gd(1,ka1),ka,ikka1,ikd,ikd,err)
c     (ensure that eigenvalue dependence is loaded into gd)
      if(kp.eq.0.or.kk.eq.0) go to 46
      do 45 is=kab1,kap
      do 45 j=iii1,iik
   45 gd(is,j)=gd(is,j)+gp(is,j,2)
c     eliminate variables of type c
   46 if(ii3.eq.0) go to 50
      if(ka.eq.0) go to 47
      call elim75(gd,p(ipu+1),ikd,ii3,ka,ik1,ii,i1,iii,iii1,1)
   47 if(kp.eq.0) go to 50
      call elim75(gd(kab1,1),p(ipu+1),ikd,ii3,kp,ik1,ii,i1,iii,iii1,1)
c
c     first inversion (for equation 5)
   50 detsgn=1.d0
      if(ka.ne.0) go to 52
      det=0.d0
      err=1.d0
      go to 56
   52 continue
c..      write(6,*) 'In loop 52'
c..      call wrleq(gd,gd(1,ka1),ka,ikka1,ikd,ikd,err)
      call leq(gd,gd(1,ka1),ka,ikka1,ikd,ikd,err)
      if(err)  54,400,56
   54 detsgn=-detsgn
   56 det=log10( abs(err))
c     note that now gd(ial,nu)=-h(ial,nu,1) in equation 5
c
c     initial integration coefficients
      hn=(x(2)-x(1))/6.d0
      hn1=(x(3)-x(2))/6.d0
      if(ki.eq.0) go to 60
      wa=hn+hn1
      wb=hn-hn1
      wd=wa/hn
      ah=wd*(hn+wb)
      ah1=wa*wa*wd/hn1
      ah2=wa*(hn1-wb)/hn1
c
c
c     compute p and gamma (gd) at first mesh point
c     and store p in real*8 q for computation of e in loops 103 and 105
   60 if(ka.eq.0) go to 65
      do 64 nu=ka1,ik1
      j=(nu-ka-1)*ii+ika
      do 61 ial=1,ka
      wd=-gd(ial,nu)
      q(ika+ial,nu-ka)=wd
   61 p(ial+j)=wd
      if(kp.eq.0) go to 64
      do 63 is=kab1,kap
      wd=0.0d0
      do 62 ial=1,ka
   62 wd=wd-gd(is,ial)*gd(ial,nu)
   63 gd(is,nu)=wd+gd(is,nu)
   64 continue
c     contribution at first boundary to integral constraints
   65 if(ki.eq.0) go to 68
      do 67 nu=ka1,ik1
      do 67 ig=1,ki
      wd=0.0d0
      if(ka.eq.0) go to 67
      do 66 ial=1,ka
   66 wd=wd-hd(ig,ial)*gd(ial,nu)
   67 gd(ig+kap,nu)=ah*(wd+hd(ig,nu))
   68 continue
c
c
c
c     **********
c     beginning of preliminary outer loop
c
c     set storage indices
   70 in=1
      in1=2
c
   80 do 130 n=1,n1
      np1=n+1
      dx=3.d0*hn
c
c     check that independent variable is strictly monotonic
      if(r*dx.le.0.) go to 408
c
c     compute integration coefficients
      hn=hn1
      if(n+3.gt.nn) go to 81
      hn1=(x(n+3)-x(n+2))/6.d0
   81 if(ki.eq.0) go to 90
      if(n.eq.n1) go to 85
      if(in.eq.1) go to 84
c     check whether np1=nn-1 when nn is even
      if(np1.eq.n1) go to 84
      wa=hn+hn1
      wb=hn-hn1
      wd=wa/hn
      ah=wd*(hn+wb)+ah2
      ah2=wa/hn1
      ah1=wa*wd*ah2
      ah2=ah2*(hn1-wb)
      go to 90
c     integration coefficients at np1 when np1 is even
c     and nn-1 when nn is even
   84 ah=ah1
c     check whether np1=nn-2
      if(n.ne.nn-3) go to 90
c     integration coefficients at nn-2,nn-1 and nn when nn is even
      wa=hn+hn1
      wb=3.d0*hn+hn1
      ah=ah-hn1*hn1*hn1/(hn*wa)
      ah1=hn1*wb/hn+ah2
      ah2=hn1*(wb+hn1)/wa
      go to 90
c     integration coefficient at nn
   85 ah=ah2
c
c     equations at mesh point np1
   90 continue
      ipv=ipv+ik13
      call ctdt75(x(n),y(1,n),zk,yp(1,n),zkp,z,zp,dzdy,theta,dt,ap,aq,
     .  rhs,f,fd,alam,alamd,a,d,h,hd,p(ipv),u,v,id,iid,ihd,ii1,ii2,
     .  ii3,iii,kk,ki,n,iter)
c     test for singularity of a (in equation 2)
   96 if(v(1).eq.0) return
c     beginning of loop to construct e matrix
  100 do 109 i=1,ii
c
c     construct d matrix (b block is already loaded)
      if(ka.eq.0) go to 107
      do 104 m=i1,ik1
      wd=0.0d0
      do 103 ial=1,ka
  103 wd=wd+a(i,ial)*q(ika+ial,m-ka)
  104 d(i,m)=d(i,m)+wd
c
c     construct beginning of e matrix and set into end of a
      if(ika.eq.0) go to 107
      do 106 j=1,ika
      wd=0.0d0
      do 105 ial=1,ka
  105 wd=wd+a(i,ial)*q(ika+ial,j)
  106 a(i,j+ka)=wd+a(i,j+ka)
  107 continue
c
c     move d up to end of a to complete storage of e.
c     (note that if a were dimensioned exactly ii x ii or if object
c     time equivalence or storage allocation were permissible in fortran
c     this would be unnecessary. when a is exactly dimensioned, 108 loop
c     has no effect.)
      do 108 m=1,ik1
  108 a(i,ii+m)=d(i,m)
c     in 110 and 111,121 loops a(i,ii+m) will be used in place of d(i,m)
  109 continue
c     end of loop to construct e matrix
c
c
c     compute p matrix
c     and store current value in real*8 q for computation of e
c     in loops 103 and 105
c
  110 continue
      if(n.ge.500.and.n.le.550) then
c..        write(6,*) 'Loop 110, n, x =',n, x(n)
c..        call wrleq(a(1,ka1),a(1,ii+ka1),ii,ikka1,iid,iid,err)
      end if
      call leq(a(1,ka1),a(1,ii+ka1),ii,ikka1,iid,iid,err)
      if(err) 111,403,112
  111 detsgn=-detsgn
  112 k=n*ip
      det=det+log10( abs(err))
      do 115 nu=ka1,ik1
      j=(nu-ka-1)*ii+k
      do 115 i=1,ii
      wd=-a(i,ii+nu)
      q(i,nu-ka)=wd
  115 p(i+j)=wd
c
c     compute gamma matrix
  120 if(kip.eq.0.or.ika.eq.0) go to 125
      do 124 is=kab1,ik
      do 123 nu=ka1,ik1
      wd=0.0d0
      do 121 ib=ka1,ii
  121 wd=wd-gd(is,ib)*a(ib-ka,ii+nu)
      if(nu.gt.ii) go to 122
      gst(nu-ka)=wd
      go to 123
  122 gst(nu-ka)=wd+gd(is,nu)
  123 continue
      do 124 nu=ka1,ik1
  124 gd(is,nu)=gst(nu-ka)
c     integral constraints at point np1
  125 if(ki.eq.0) go to 129
      do 127 ig=1,ki
      do 127 nu=ka1,ik1
      wd=0.0d0
      if(ka.eq.0) go to 127
      do 126 ial=1,ka
  126 wd=wd-hd(ig,ial)*a(ika+ial,ii+nu)
  127 gd(ig+kap,nu)=gd(ig+kap,nu)+ah*(wd+hd(ig,nu))
c
c     reset storage indices
  129 i=in
      in=in1
      in1=i
  130 continue
c
c     end of preliminary outer loop
c     **********
c
c     eliminate variables of type c at final meshpoint
  140 if(ii3.eq.0) go to 200
      ipv=ipu+n1*ik13+1
      if(kb.eq.0) go to 145
      call elim75(gd(ka+1,1),p(ipv),ikd,ii3,kb,ik1,ii,i1,iii,iii1,1)
  145 if(kp.eq.0) go to 200
c     ensure that gp(is,iii+k,2) = 0
      do 146 is=kab1,kap
      do 146 k=iii1,iik1
  146 gp(is,k,2)=0.d0
c
      call elim75(gp(kab1,1,2),p(ipv),ikd,ii3,kp,ik1,ii,i1,iii,iii1,1)
c     reset gamma
      do 148 is=kab1,kap
      do 148 k=i1,ik1
  148 gd(is,k)=gd(is,k)+gp(is,k,2)
c
c     **********
c     remaining boundary conditions
c
  200 if(ka.eq.ik) go to 233
c     compute coefficients for equation 12a and load into end of gd
      if(kp.eq.0) go to 210
      do 205 is=kab1,kap
      if(ka.eq.0) go to 203
      do 202 nu=ka1,ik1
      wd=0.0d0
      do 201 ial=1,ka
  201 wd=wd-gp(is,ial,2)*a(ika+ial,ii+nu)
  202 gd(is,nu)=gd(is,nu)+wd
  203 do 204 ia=ka1,ii
  204 gd(is,ia)=gd(is,ia)+gp(is,ia,2)
  205 continue
  210 continue
c
c
c     compute coefficients for equation 13 and load into middle of gd
  220 if(ka.eq.0.or.kb.eq.0) go to 230
      do 223 m=ka1,kab
      do 222 nu=ka1,ik1
      wd=0.0d0
      do 221 ial=1,ka
  221 wd=wd-gd(m,ial)*a(ika+ial,ii+nu)
  222 gd(m,nu)=wd+gd(m,nu)
  223 continue
c
c     solve equations 12a and 13 (answer has wrong sign)
  230 if(ikka.eq.0) go to 233
  231 continue
c..      write(6,*) 'In loop 231'
c..      call wrleq(gd(ka1,ka1),gd(ka1,ik1),ikka,1,ikd,ikd,err)
      call leq(gd(ka1,ka1),gd(ka1,ik1),ikka,1,ikd,ikd,err)
      if(err) 233,404,234
  233 detsgn=-detsgn
  234 det=det+log10( abs(err))
      gd(ik1,ik1)=-1.0d0
      gd(ik1,1)=-1.0d0
c     (note that the first element of the second half of gp may now
c     have been overwritten)
c
c
c
c     **********
c     iterated solution
c
c     empty ea,eb
  300 do 301 i=1,iii
      eb(i)=0.0d0
      do 301 j=1,3
  301 ea(i,j)=0.0d0
c
c     solution at second boundary
  310 if(kk.eq.0) go to 313
      do 312 k=1,kk
      wd=gd(ii+k,ik1)
      gd(ii+k,1)=wd
      zk(k)=zk(k)-ucy*wd
  312 continue
  313 if(ka.ge.ii) go to 320
      do 318 ia=ka1,ii
      wa= abs(gd(ia,ik1))
      wd=y(v(ia),nn)
  315 wd=wd-ucy*gd(ia,ik1)
c
c  test for storing solution
c
      if(idiag.ge.1) then
	dyst(v(ia)+iii*(nn-1))=-ucy*gd(ia,ik1)
      end if
c
      wb= abs(wd)
      ea(v(ia),1)=wa
      eb(ia)=wb
      if(wb.lt.1.0e-10) go to 316
      ea(v(ia),2)=wa/wb
      ea(v(ia),3)= float(nn)
  316 continue
      y(v(ia),nn)=wd
  318 continue
c
c     remainder of solution   (beginning of second outer loop)
c
c     set storage indices
  320 in=1
      in1=ik1
      ipv=ipu+nn*ik13
c
      do 339 m=1,nn
      n=nn+1-m
      np1=ip*(n-1)
      ipv=ipv-ik13
      do 331 i=1,iii
      if(i.gt.ika) go to 321
      k=n-1
      ini=in
      if(k.lt.1) go to 331
      j=i+ka
      go to 322
  321 if(i.gt.ii) go to 324
      j=i-ika
      k=n
      ini=in1
  322 wd=0.0d0
      ia=i+np1
      do 323 nu=ka1,ik1
  323 wd=wd+p(ia+(nu-ka-1)*ii)*gd(nu,in1)
      go to 326
c     y's from equations of type c
  324 j=i
      k=n
      l=ipv+i-iii
      wd=0.d0
      do 325 iff=1,ik1
      l=l+ii3
  325 wd=wd+p(l)*gd(iff,in1)
  326 wa= abs(wd)
      wb=y(v(j),k)-ucy*wd
      y(v(j),k)=wb
  328 wb= abs(wb)
c
c  test for storing solution
c
      if(idiag.ge.1) then
	dyst(v(j)+iii*(k-1))=-ucy*wd
      end if
c
      ea(v(j),1)=ea(v(j),1)+wa
      eb(j)=eb(j)+wb
c..      if(v(j).eq.4) write(6,'(a,i5,1pe20.12,3e13.5)') 'COR: ',
c..     *  n, y(v(j),k), wa, ea(v(j),1), eb(j)
      if(wb.lt.1.0e-10) go to 329
      wb=wa/wb
  329 if(wb.lt.ea(v(j),2)) go to 330
      ea(v(j),2)=wb
      ea(v(j),3)= float(k)
  330 continue
      if(i.gt.ii) go to 331
      gd(j,ini)=wd
  331 continue
c
c     reset storage indices
  338 i=in
      in=in1
      in1=i
  339 continue
c
c     end of second outer loop
c
c     iteration complete
c     **********
c
c
c     set mean relative corrections
  340 do 341 i=1,ii
  341 ea(v(i),1)=ea(v(i),1)/eb(i)
c     set det
      if(detsgn) 345,350,350
  345 det=det*1.d10
c
c
  350 return
c
c
c     diagnostics
c     **********
  400 write(iw,1001)
  401 do 402 k=1,2
      do 402 j=1,ik
  402 write(iw,1002) (gs(i,j,k), i=1,kap)
      write(iw,1101) (v(i), i=1,ii)
      go to 500
  403 write(iw,1003)n
      go to 500
  404 write(iw,1004)
      do 405 i=1,kap
      do 405 j=1,ik
      do 405 k=1,2
  405 gs(i,j,k)=0.0d0
      call bc(x(1),x(nn),y,y(1,nn),zk,ap,aq,g,gs,ikd,igd,nn)
      go to 401
  407 write(iw,1006) nn,x(1)
      go to 500
  408 dx=2.0d0*dx
      write(iw,1007) n,x(n),np1,x(np1),dx,r,nn,(x(i), i=1,nn)
      go to 500
c
c  stop for error return from s/r rhs or bc
c
  410 write(iw,1008) kdgrhb
c
  500 v(1)=0
      return
c
c
 1000 format(//1x,10('*'),10x,'improper formulation detected by nrk',
     . 10x,10('*')/17x,'ii =',i3,4x,'kk =',i3,4x,'ka =',i3,4x,'kb =',i3,
     . 4x,'ki =',i3/)
 1001 format(//1x,10('*'),5x,'first boundary condition matrix singular i
     .n nrk',5x,10('*')//' gd ='/)
 1002 format(1x,1p14e9.1/)
 1003 format(//1x,10('*'),5x,'e matrix singular at n =',i4,' in nrk',
     . 5x,10('*'))
 1004 format(//1x,10('*'),5x,'second boundary matrix singular in nrk',
     .  5x,10('*')//' gd ='/)
 1005 format(' ')
 1006 format(//1x,10('*'),5x,'null range of independent variable in ',
     .       'nrk',5x,10('*')//21x,'x(1) = x(',i4,') =',1pe14.6/)
 1007 format(//1x,10('*'),5x,'independent variable not monotonic in',
     .       ' nrk',5x,10('*')//16x, 'x(',i4,') =',1pe12.4,
     .       ',    x(',i4,') =',1pe12.4/16x,'difference =',1pe12.4,
     .       ',    range =',1pe12.4//1x,'x(n), n=1,',i4,/
     .       (1x,1p10e13.5))
c
 1008 format(//' Error return from s/r bcs. kdgrhb =',i3)
c
c     warning message
 1100 format(1x,5('*'),3x,'v set in nrk        v(',i2,') =',i3)
 1101 format(/9x,'v =',7(1x,5i3))
c
c
      end
      subroutine ctdt75(x,y,zk,yp,zkp,z,zp,dzdy,theta,dt,ap,aq,rhs,f,
     .  fd,alam,alamd,a,d,h,hd,p,u,v,id,iid,ihd,ii1,ii2,ii3,iii,kk,ki,
     .  n,iter)
      implicit double precision (a-h,o-z)
      integer v
      logical bndry,notb
      dimension x(1),y(id,1),zk(1),yp(id,1),zkp(1),theta(1),ap(1),
     .  aq(1),f(1),fd(iii,1),alam(ii2,1),alamd(ii2,iii,1),a(iid,1),
     .  d(iid,1),h(1),hd(ihd,1),p(1),u(ii3,1),v(1),z(1),zp(1),
     .  dzdy(iid,1)
      dimension dy(1,1),dyp(1,1),dzk(1),dzkp(1)
c
c  note that if double precision entry is enabled, this
c  dimensioning statement must be replaced by
c
c      dimension dy(id,1),dyp(id,1),dzk(1),dzkp(1)
c
      common/ctdist/ dx,iik,in,in1,ir,k1,ik1,ik13,ii23,ii12,ii121,ik,
     .  iii1,iik1,iiik,ii11,iw
c
      common/cdgrhb/ kdgrhb
c
      external rhs
c
c
c
c  information at previous time level
      common/sooner/ r(1)
c
      save
c
c  dimension requirement for r is
c  r((ii2+ii3+(ii1+ii2+ii3)*(ii2+1))*nn)
c
   10 bndry=n.eq.0
      notb=ii2.eq.0
      np1=n+1
      ipp1=2
c  control for first call
      if(n.gt.0) go to 15
      ipp1=1
   15 iter1=iter-1
      inn=(in-1)*iik
      inn1=(in1-1)*iik
      inn11=inn1+1
      innp1=inn+1
      ini=(in-1)*iii
      ini1=(in1-1)*iii
      ini11=ini1+1
      iinn1=iii+inn1-1
      i=np1
      ia=iii
      ib=iid
      idz=iid
      ialam=ii2
      irn=n*ir
      irnz=irn+ii23
      irl=irnz+iii
c  values at previous time level
   20 if(iter.ne.1.or.ii23.eq.0) go to 30
      call rhs(x(ipp1),yp(1,ipp1),zkp,z(ini11),dzdy(1,inn11),ap,aq,
     .  f(ini11),fd(1,inn11),alam,alamd,h,d,idz,ia,ib,ialam,i,iter1)
      if(kdgrhb.eq.-1) go to 905
c  load z into r
   22 m=ini1
      do 23 j=1,iii
      m=m+1
   23 r(irnz+j)=z(m)
c  load f into r
      m=ini1+ii1
      if(ii23.eq.0) go to 30
      do 27 l=1,ii23
      m=m+1
   27 r(irn+l)=f(m)
c  load alam into r
      if(notb) go to 30
      j=irl
      do 28 l=1,iii
      do 28 k=1,ii2
      j=j+1
   28 r(j)=alam(k,l)
c  empty derivative matrices
   30 do 32 j=1,iii
      do 32 l=1,iii
   32 dzdy(l,j+inn1)=0.d0
      do 36 k=1,iik
      do 36 j=1,iii
   36 fd(j,k+inn1)=0.d0
      if(notb) go to 38
      do 37 k=1,iik
      do 37 j=1,iii
      do 37 l=1,iii
   37 alamd(l,j,k)=0.d0
   38 if(ki.eq.0) goto 40
      do 39 ig=1,ki
      do 39 j=1,ik
   39 d(ig,j)=0.d0
c
   40 iter1=1
      call rhs(x(ipp1),y(1,ipp1),zk,z(ini11),dzdy(1,inn11),ap,aq,
     .  f(ini11),fd(1,inn11),alam,alamd,h,d,idz,ia,ib,ialam,i,iter1)
      if(kdgrhb.eq.-1) go to 905
c  set hd
   42 if(ki.eq.0) go to 46
c
      do 43 i=1,iik
      iv=v(i)
      if(i.gt.iii) iv=i
      do 43 ig=1,ki
   43 hd(ig,i)=d(ig,iv)
      do 44 ig=1,ki
   44 hd(ig,iik1)=h(ig)
c  get z at previous time level from r
   46 if(ii23.eq.0) go to 50
      m=ini1
      do 47 i=1,iii
      m=m+1
   47 zp(m)=r(irnz+i)
c
c  equations of type c
c  -------------------
c
   50 if(ii3.eq.0) go to 60
c  set curly a, c and d
      do 59 iu=ii121,iii
      i=iu-ii12
      thetau=theta(iu)
      thetat=thetau*dt
      do 51 m=1,ii12
      iv=v(m)+inn1
   51 u(i,m)=thetat*fd(iu,iv)-dzdy(iu,iv)
      m=iu-ii1
      if(kk.eq.0) go to 55
   53 do 54 k=1,kk
      j=ii12+k
      l=iinn1+k
   54 u(i,j)=thetat*fd(iu,l)
   55 yy=z(iu+ini1)-zp(ini1+iu)
      u(i,ik1)=dt*(thetau*f(iu+ini1)+(1.0d0-thetau)*r(m+irn))
     .  -yy
      do 59 m=ii121,iii
      iv=m+inn1
c..      write(6,*) iu,i,m,iv,dzdy(iu,iv),fd(iu,iv),-thetat*fd(iu,iv)
      u(i,k1+m)=dzdy(iu,iv)-thetat*fd(iu,iv)
c..      write(6,*) i,k1+m,u(1,k1+m)
   59 continue
c  solve equation 2 to obtain u
c..      write(6,*) 'ii3 =',ii3, ' k1, ik1 =',k1, ik1
c..      do 59010 i=1,ii3
c..59010 write(6,*) i,(ik1+m,u(i,ik1+m),m=1,ii3)
c..      write(6,*) 'In loop 59'
c..      call wrleq(u(1,ik1+1),u,ii3,ik1,ii3,ii3,det)
      call leq(u(1,ik1+1),u,ii3,ik1,ii3,ii3,det)
      if(det) 60,901,60
c
c  equations of type b
c  -------------------
   60 if(ii2.eq.0) go to 80
      ddt=1.d0/dt
      do 72 ip=ii11,ii12
      ip1=ip-ii1
      thetap=theta(ip)
      thetau=1.d0-thetap
c  compute e at np1 and store in fd
      iir=ip1+irl-ii2
      do 64 i=1,iii
      iir=iir+ii2
      ee=0.d0
      do 62 j=1,iii
      yy=y(j,ipp1)-yp(j,ipp1)
   62 ee=ee+alamd(ip1,j,i)*yy
   64 fd(ip,i+inn1)=ddt*(thetap*(fd(ip,i+inn1)*dt+ee+alam(ip1,i))
     .  +thetau*r(iir))
c  compute g and load into fd
      if(kk.eq.0) go to 68
      do 67 k=1,kk
      j=iii+k
      l=j+inn1
      ee=0.d0
      do 66 i=1,iii
      yy=y(i,ipp1)-yp(i,ipp1)
   66 ee=ee+alamd(ip1,i,j)*yy
   67 fd(ip,l)=thetap*ddt*(fd(ip,l)*dt+ee)
c  compute h and store in f
   68 iir=ip1+irl-ii2
      ee=0.d0
      ff=0.d0
      do 71 i=1,iii
      iir=iir+ii2
      yy=y(i,ipp1)-yp(i,ipp1)
   70 ee=ee+alam(ip1,i)*yy
   71 ff=ff+r(iir)*yy
      f(ip+ini1)=ddt*((thetap*f(ip+ini1)+thetau*r(ip1+irn))*dt
     .  +thetap*ee+thetau*ff)
   72 continue
c  set a,b,c and d (b,c and d are loaded into d)
      if(bndry) go to 100
      do 78 ip=ii11,ii12
      ivp=ip
      thetap=theta(ip)
      thetau=1.d0-thetap
      do 73 i=1,iii
      iv=v(i)
      ivn=iv+inn
      ivn1=iv+inn1
      a(ip,i)=dx*fd(ivp,ivn)+thetap*dzdy(ivp,ivn)
   73 d(ip,i)=dx*fd(ivp,ivn1)-thetap*dzdy(ivp,ivn1)
c
      if(kk.eq.0) go to 75
      do 74 k=iii1,iik
   74 d(ip,k)=dx*(fd(ivp,k+inn)+fd(ivp,k+inn1))
c
   75 yy=z(ini1+ivp)-z(ini+ivp)
      yyp=zp(ini1+ivp)-zp(ini+ivp)
      d(ip,iik1)=dx*(f(ivp+ini)+f(ivp+ini1))
     .  -thetap*yy-thetau*yyp
c..      if(mod(n,10).eq.0.and.ip.eq.ii1+3)
c..     *  write(6,'(a,3i4,1p12e13.5)') 'In tnrkt:', 
c..     *  iter, ip, n, x, z(ini1+ivp), zp(ini1+ivp),
c..     *  yy, yyp, dx, ee*ddt, f(ivp+ini), f(ivp+ini1),
c..     *  dx*(f(ivp+ini)+f(ivp+ini1)), thetap,
c..     *  d(ip,iik1)
   77 continue
   78 continue
c
c
c  equations of type a
c  -------------------
c
   80 if(bndry) go to 100
      if(ii1.eq.0) go to 90
c  compute a,b,c,d (b,c,d are loaded into d)
   81 do 87 l=1,ii1
      ivl=l
      do 82 i=1,iii
      iv=v(i)
      ivn=iv+inn
      ivn1=iv+inn1
      a(l,i)=dx*fd(ivl,ivn)+dzdy(ivl,ivn)
   82 d(l,i)=dx*fd(ivl,ivn1)-dzdy(ivl,ivn1)
      if(kk.eq.0) go to 84
      do 83 k=iii1,iik
   83 d(l,k)=dx*(fd(ivl,k+inn)+fd(ivl,k+inn1))
   84 yy=z(ini1+ivl)-z(ini+ivl)
   86 d(l,iik1)=dx*(f(ivl+ini)+f(ivl+ini1))-yy
   87 continue
c..      write(6,*) 'After 87, a, d:'
c..      call wrleq(a,d,ii1,iik1,iid,iid,err)
c
c  computation of a,b,c and d thilde
c  ---------------------------------
c
   90 if(ii3.eq.0) return
      do 97 ia=1,ii12
c  a and b (b is loaded in end of d)
      l=-ii3
      do 92 ib=1,ii12
      l=l+ii3
      ee=0.d0
      ff=0.d0
      do 91 iu=ii121,iii
      i=iu-ii12
      ee=ee+a(ia,iu)*p(i+l)
   91 ff=ff+d(ia,iu)*u(i,ib)
      a(ia,ib)=a(ia,ib)+ee
   92 d(ia,ib)=d(ia,ib)+ff
c  c and d (loaded in d)
      do 94 k=1,k1
      j=ii12+k
      l=l+ii3
      m=iii+k
      ee=0.d0
      do 93 iu=ii121,iii
      i=iu-ii12
   93 ee=ee+a(ia,iu)*p(i+l)+d(ia,iu)*u(i,j)
   94 d(ia,m)=ee+d(ia,m)
c  move c and d to start at ii121
      j=ii12
      do 95 m=iii1,iik1
      j=j+1
   95 d(ia,j)=d(ia,m)
   97 continue
c
c  eliminate variables of type c from integral constraints
  100 if(ii3.eq.0) return
      if(ki.eq.0) go to 110
      call elim75(hd,u,ihd,ii3,ki,ik1,ii12,ii121,iii,iii1,2)
c  store u in p
  110 l=ik13-ii3
      do 111 iff=1,ik1
      l=l+ii3
      do 111 iv=1,ii3
  111 p(l+iv)=u(iv,iff)
c
  120 return
c
c  diagnostics
c  ***********
c
  901 write(iw,1001) dt, n, x
      do 902 i=1,iii
  902 write(iw,1002) i,theta(i),f(i+ini1),(fd(i,j+inn1),j=1,iik)
      v(1)=0
      return
c
c  stop for error return from s/r rhs or bc
c
  905 write(iw,1005) kdgrhb
      v(1)=0
      return
 1001 format(//,1x,5('*'),3x,'a matrix in equation 2 singular in tnrk',
     .  3x,5('*')//' time step =',1pe12.4,'  n =',i6,
     .  '  x =',1pe15.7//' i theta f fd ='/)
 1002 format(1x,i2,f6.3,1pe12.4,2x,9e12.4/,(23x,9e12.4))
c
 1005 format(//' Error return from s/r rhs. kdgrhb =',i3)
      end
      subroutine trns76(a,dzdy,ii,kk,ia,id,v)
c     transforms array a(k,i), i=1,ii, k=1,kk, with array dzdy
c      under permutation v
      implicit double precision (a-h,o-z)
      integer v
      dimension a(ia,1),dzdy(id,1),v(1)
      dimension a1(20)
c
      save
c
      if(kk.eq.0) return
c
   10 do 30 k=1,kk
      do 20 i=1,ii
      sum=0.d0
      do 15 j=1,ii
   15 sum=sum+a(k,j)*dzdy(j,i)
   20 a1(i)=sum
      do 30 i=1,ii
   30 a(k,i)=a1(v(i))
      return
      end
      subroutine elim75(ad,u,iad,ii3,nr,nc,ii,i1,iii,iii1,ipu)
c  eliminates variables of type c from array ad
c   nc: number of columns in the final array
c   nr: number of rows dealt with
c   iad: first dimension of ad
c   i1 = ii1 + ii2 + 1
c
      implicit double precision (a-h,o-z)
      dimension ad(iad,4),u(1)
c
      save
c
      iku=ipu*ii3
      ikuint=-iku-ipu+1
   10 do 20 i=1,nr
      j=0
      ikku=ikuint
      do 14 iff=1,nc
      ikku=ikku+iku
      ku=ikku
      j=j+1
      if(iff.eq.i1) j=iii1
      ee=0.d0
      do 12 iu=1,ii3
      ku=ku+ipu
   12 ee=ee+ad(i,ii+iu)*u(ku)
   14 ad(i,j)=ad(i,j)+ee
c  move the last part of ad up
      if(nc.le.ii) go to 20
      j=iii
      do 16 k=i1,nc
      j=j+1
   16 ad(i,k)=ad(i,j)
   20 continue
      return
      end
      subroutine wrleq(a,b,nn,mm,ia,ib,err)
      implicit double precision(a-h,o-z)
      dimension a(ia,1),b(ib,1)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
      data icount /20/
c
      if(istdpr.gt.0) then
        write(istdpr,100)
        do 10 i=1,nn
   10   write(istdpr,110) (a(i,j),j=1,nn)
c
        write(istdpr,120)
        do 20 i=1,nn
   20   write(istdpr,110) (b(i,j),j=1,mm)
      end if
      if(icount.gt.0) then
	icount=icount-1
        return
      else
c..	stop
      end if
  100 format(' In leq, a:')
  110 format(1p7e11.3)
  120 format(' In leq, b:')
      end
