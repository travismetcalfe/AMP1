      character*(*) function settrailer(am, z, icase, xseq, itcase)
c
c  returns model trailer corresponding to mass am, heavy-element
c  abundance Z and case number icase
c
c  If xseq .lt. 0 this is assumed to be an evolution sequence, and
c  the trailer ends in `.s' 
c  Otherwise xseq is the sequence number in the evolution sequence.
c  A real variable is used to allow interpolated models.
c
c  itcase flags for trailer case:
c  itcase .le. 0: trailer of form <mass>.Z<Z>.<case>.<s or xseq>
c  itcase    = 1: trailer of form <mass>.Z<Z>.X<X>.<case>.<s or xseq>
c
c  Variables other than mass and Z are taken from common/cvr_param/
c  (might need a check!)
c
c  Modified 4/10/06, to add choice of trailer format.
c
      implicit double precision(a-h, o-z)
      character*80 int2str, sam, sz, sxxh, scase, strcompr, sseq
c
      common/cvr_param/ par_am, par_z, par_agefin, par_rsfin, 
     *  par_alsfin, par_zxsfin, par_xmdtrl,
     *  par_xxh, par_fdgopl, par_tlopfg, par_dtlopf, par_tlopf1,
     *  par_fcno, par_agehe3, par_xrz12, par_xrz13, par_xrz14,
     *  par_xrz16, par_tprfct, par_alfa, par_etac, par_phc, par_alpove,
     *  par_alpovc, par_velrot
      common/cvi_param/ ipar_nt, ipar_icsove, ipar_icsovc, ipar_isprot,
     *  ipar_icmout, ipar_istosc, ipar_isetos
c
c  common defining standard input and output, standard error
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      iam=100*am+0.5
      iz=100000*z+0.5
      ixxh=10000*par_xxh+0.5
c
      sam = strcompr(int2str(iam))
      sz = strcompr(int2str(iz))
      sxxh = strcompr(int2str(ixxh))
      scase = strcompr(int2str(icase))
c
c  test for sequence
c
      if(xseq.lt.0) then
	sseq='s'
      else if(abs(xseq-intgpt(xseq)).le.1.d-7) then
	write(sseq,'(i4)') intgpt(xseq)
      else
	write(sseq,'(f10.3)') xseq
      end if
c
      sseq=strcompr(sseq)
c
      lam=length(sam)
      lz=length(sz)
      lxxh=length(sxxh)
      lcase=length(scase)
      lseq=length(sseq)
c
      if(itcase.le.0) then
        settrailer='.'//sam(1:lam)//'.Z'//sz(1:lz)//
     *  '.'//scase(1:lcase)//'.'//sseq(1:lseq)
      else
        settrailer='.'//sam(1:lam)//'.Z'//sz(1:lz)//
     *  '.X'//sxxh(1:lxxh)//'.'//scase(1:lcase)//'.'//sseq(1:lseq)
      end if
      return
      end
      character*(*) function int2str(i)
c
c sets string based on integer i, left-filled with zeros
c
      character*4 s, s1
      character*80 ss, strcompr
      i1=i
      if(i.lt.0) i1=-i
      write(s,'(i4)') i1
c
      if(i1.lt.10) then
	ss='000'//s
      else if(i1.lt.100) then
	ss='00'//s
      else if(i1.lt.1000) then
	ss='0'//s
      else
	ss=s
      end if
      s1=strcompr(ss)
      int2str = s1
      return
      end
