      subroutine mnevol
c
c  main programme for stellar evolution, with new equation of state
c  ****************************************************************
c
c  modified 20/8/1984 to rescale variables and coefficients in
c  central boundary condition routines, in order to avoid
c  underflow problems on recku univac.
c
c  the rescaling is flagged in models output to file by adding
c  10 000 000 to case number icase in datmod.
c
c  modified 29/12/1984, to make call of central boundary conditions
c  after iteration has converged, to get consistent set of
c  coefficients before output to file.
c  for the same reason, he3 abundance is set after convergence,
c  if he3 fudge is used, or he3 is in equilibrium.
c
c  modified 4/1/1985 to use log(r/1.d11) and log(l/1.d33) as
c  dependent variables, to reduce rounding error. these
c  variables are also used in the output to file of evolution
c  models.
c  the modification is flagged by setting icase/1.d7 to 2.
c
c  Modified 25/8/87, so that now uw(1) in common/rnrout/ returns
c  Uw. Previously it returned Uw/ln(10).
c
c  Modified 17/12/87, to output also extensive GONG set of variables
c
c  Modified 18/12/87, replacing common/time/ by common/cmtime/ (to
c  avoid conflict with system function time).
c
c  Modified 23/12/87, including factor albdbc in surface condition
c  on L when ntau .le. 1
c
c  Modified 18/1/88, changing definition of sbcfct and albdbc in
c  surface boundary conditions. Note that if sbcfct = 1 and 
c  albdbc = 1 (the previous and new default) the change has no
c  effect. Also change slightly the definition of output in datmod,
c  so that now sbcfct and, when ntaubc .le. 1 albdbc, are included.
c  This modification is flagged by setting icase/1.d7 to 3.
c
c  Also add Ne and rX to GONG model output.
c
c  Modified 19/1/88, updating derivatives in derlr based on
c  new evolution calculations.
c
c  Modified 15/3/88 to use consistent numerical constants
c  This modification is flagged by setting icase/1.d7 to 4.
c
c  Modified 20/2/89: iomfll added to input, to allow run-time
c     choice of correct statistical weight for heavy elements.
c     (see note in s/r eqstf.f of 5/1/84).
c
c  Modified 1/3/89: when resetting mesh in trial model with
c     itrial = 4, also reset radius at innermost mesh point.
c
c  Modified 16/3/89, to streamline setting of array sizes.
c     Should now be done in main programme.
c
c  Modified 12/4/89, rescaling also azt(3) and ax(3) when X is rescaled.
c     The omission to do so caused problems when resetting innermost
c     meshpoint in trial model (cf. modification 1/3/89)
c
c  Modified 11/5/89, adding option for ignoring entropy term in
c     energy equation.
c
c  Modified 5/6/89, adding option for smooth transition in temperature
c     gradient between specified T(tau) relation in atmosphere and
c     diffusion approximation in interior. This is flagged by icsrad.
c     This change is flagged in icase by increasing the leading digit
c     from 4 to 5. Also icase was changed to accomodate icsrad, 
c     by replacing 100*iopacm+1000*iopatm by 
c     100*(iopacm+2*iopatm)+1000*icsrad.
c
c  Modified 14/6/89, adding special output option flagged by parameter
c     ispcpr. This must be taken care of by user supplied routine
c     spcout. Also added flag icncbc to control central boundary
c     condition routine.
c     Furthermore age was added to common/cmtime/
c
c  Modified 23/6/89, re-installing the fudged mixing length theory
c     of C-D (Cambridge, 1986).
c
c  Modified 5/3/90, to increase size of arrays for atmospheric solution
c  from 51 to 201
c
c  Modified 7/3/90, adding input parameter icasex to choose 
c  optimal setting of r/R in adiabatic pulsation model. 
c  icasex = 1 corresponds to using optimal setting. 
c  Otherwise old setting is used.
c  Using icasex = 1 is flagged in adiabatic model by setting
c  data(8) = 1 (standard value is 0).
c
c  Modified 4/4/90, introducing different types of X-interpolation in
c  opacity, determined by the value of iwdopc, when iwdopc .ne. 1.
c
c  Modified 5/8/90, resetting X to 1.d-10 when it is less than the limit
c  xhzlm1 (currently hard-coded to be xhzlm1 = 1.d-6)
c
c  Modified 7/8/90, to reset equations for X in s/r rhs 
c  for n .ge. nxhzer where nxhzer is determined in s/r rhs 
c  based on the estimated (from the values at the previous time step) 
c  value of X is less c  than xhzlm1. - 
c  Also, in s/r engenr eps and rX are set to zero for
c  X .le. xhzlm2 (currently hard-coded to xhzlm2 = 1.d-7).
c
c  Modified 10/8/90, reducing timestep for evolution with 
c  convective core (see s/r setdt).
c
c  Modified 10/8/90, to reset outer boundary of a convective core
c  consistently with the mixed composition, by calling s/r conmxc
c  after iteration has converged.
c
c  Modified 12/8/90, to fix setting of coefficients al0 and al2 in
c  common/ksider/ (cf. s/r bciter).
c
c  Modified 12/8/90, including missing factor am in coefficient a3 for
c  s/r rhs. This has major effect on models of non-solar mass stars,
c  of course (but was probably ok before resetting coefficients
c  consistently in March 1988).
c
c  Modified 13/8/90, including missing term aml in constant b4 for 
c  s/r bcs. This error (like the error in a4) was probably introduced
c  in resetting coefficients "consistently" in March 1980.
c  Note, however, that it only has an effect for the (rarely used)
c  simple surface boundary condition without an atmosphere.
c
c  Modified 14/8/90, extending mixing in outermost convection zone
c  to the surface. Introduce re-initialization of opacity tables
c  if envelope composition has been changed due to dredge-up.
c
c  Modified 14/8/90, to fix up setting of ranges in s/r opact.
c
c  Modified 17/8/90, setting timestep and extrapolating before 
c  mesh stretching.
c
c  Modified 22/8/90, computing evolution of composition in convective
c  core by time-integrating (using essentially a predictor-corrector
c  method) with the mass-averaged reation rate. There are still 
c  remaining problems at the core boundaries which must be looked
c  into.
c  
c  Modified 28/8/90, correcting setting of Q-values, by reversing 
c  change in srncns which were introduced 20/3/89 to correct for what
c  was erroneously precieved as an error. New, corrected values
c  set in file srncns-3.d.f
c
c   Summary of versions of routine srncns:
c
c   srncns.d.f: Original setting of FCZ75 parameters.
c               Correct, after all.
c   srncns-1.d.f: Revised setting of FCZ75 parameters.
c                 *** Erroneous ****
c   srncns-2.d.f: Revised setting of Parker (1986) parameters.
c                 *** Erroneous ****
c   srncns-3.d.f: Original setting of Parker (1986) parameters.
c                 Correct (we hope).
c
c  Modified 29/8/90: Take out the resetting of convective core
c  parameters in s/r conmxc (no longer needed after introducing
c  consistent evolution of convective core composition)
c  Leave call of conmxc for checking, diagnostic output.
c
c  Modified 31/8/90, to correct error in temporary output of model
c  sequence to file during R, L iteration (due to error in
c  logics the ZAMS model was not output. Error 
c  introduced after August 13 1990).
c
c  Modified 16/10/90, to introduce treatment of convective
c  mixing for several elements (previously He3 was assumed
c  to be in equilibrium).
c
c  Modified 17/10/90, to introduce special treatment of He3 in 
c  convective core (analytical solution with averaged reaction rates
c  per pair of particles).
c
c  Modified 18/10/90, changing output on unit idscen to be formatted
c  output of central values.
c
c  ********************************************************************
c
c  Major modifications, starting July 91, to incorporate full treatment
c  of CNO cycle, and update output formats for more flexibility.
c
c  17/7/91 - 31/7/91: Change storage in s/r engenr etc., add 
c  first part of CNO cycle. Update, generalize storage in 
c  s/r bcs and bciter.
c
c  2/8/91: Incorporate rough treatment of convective overshoot
c  from convective envelope in standard version.
c
c  2/8/91 - 5/8/91: Set input parameters, define storage for abundance
c  variables, reset treatment of chemical evolution of convective core.
c
c  --------------------------------------------------------------------
c
c  Modified 13/5/92: Adding convection zone details to GONG output,
c  further modifications to output format. Flagged by changing
c  iform from -1 to -2
c
c  Modified 8/8/92, including option iwdopc = 9 to use WD/GH
c  routines for interpolation in OPAL and Kurucz tables
c
c  Modified 11/12/92, implementing halving of time step at 
c  convergence failure
c
c  Modified 20/2/93, allowing choice of nuclear parameters through
c  setting the new input parameter ivreng
c
c  Modified 10/9/93, to reset constants to double precision.
c
c  Modified 23/3/94, setting default iscren to 1 (this setting
c  had been inadvertently dropped)
c
c  Modified 19/7/95 to move heavy elements into array 
c  in common /heavy/, allowing Z varying with position.
c  This is in preparation for including heavy-element diffusion.
c
c  Modified 5/8/95, correcting erroneous resetting of X in outer
c  part of previous model when using OPAL opacities.
c
c  Modified 12/7/96, including JOP additions to energy generation
c  (principally to start pre-main-sequence evolution)
c
c  Modified 3/8/96, adding setting cqcp to zero if previous model
c  had no convective core, and fixing several instances of the
c  heavy-element abundance in s/r cmpcvc.
c
c  Modified 4/8/96, taking out option of writing data to disk with
c  ifwrt = 1 (but retaining ifwrt in common/rhcn/ for consistency)
c
c  Modified 5/8/96, adding more flexibility in control of 
c  mixed cores, with imixcr, and avoiding double setting
c  of mixed core with s/r mixcor for initial time step.
c
c  Modified 13/8/96, including new version number in ndtmod and
c  updating version number in icase to 9.
c
c  Modified 19/8/96, introducing effects of centrifugal force
c
c  Modified 24/8/96, converting Livermore routine to double precision
c
c  Modified 24/8/96, incorporating MJM convection treatment.
c
c  Modified 1/7/97, resetting lastmd to false when time step is
c  halved, to avoid problems in reaching last model
c
c  Modified 16/2/98, correcting setting of omgrot (previously too small
c  by factor 2*pi)
c
c  Modified 10/6/98, preparing for incorporating RT fits to T(tau) and
c  alpha (not yet done with full consistency!)
c
c  Modified 13/10/98, allowing extended iwdopc for use with GH
c  package version 9 (and later) (Note: implemented in evolmain.dnnz.d.f
c  on 1/4/98)
c
c  Modified 13/10/98, to treat He3 evolution in cases with a
c  convective core. Currently this is accomplished by resetting
c  the He3 abundance to the equilibrium value in the region
c  that was convective in the *previous* model. (Note: implemented 
c  in evolmain.dnnz.d.f on 21/4/95)
c
c  Modified 13/10/98, to establish consistency between this routine
c  and the diffusive version in evolmain.dnnz.d.f
c
c  Modified 13/11/98, implementing more flexible use of RT atmosphere
c  fit
c
c  Modified 16/11/98, extending datout to 150 variables, to store
c  full qqbc array
c
c  Modified 11/6/99, to allow setting Z from trial model in cases
c  with heavy-element settling
c
c  Modified 21/6/99, adding flsatm in datout(140), and using it for
c  ipartr = 1
c
c  Modified 16/11/99, adding storage of x at convection-zone boundaries
c  and of convection-zone parameters at previous time step
c
c  Modified 5/1/00, to add more general treatment of convective overshoot
c  from convective envelope and core.
c
c  Modified 15/5/00, introducing flag inmixc, such that inmixc .gt. 0
c  switches off various aspects of convective-core mixing.
c
c  Modified 15/5/00, redefining istrtc to flag also for resetting
c  of composition variables in cvr.
c
c  Modified 20/7/00, adding various options to emerge from 
c  convergence problems, particularly in case of growing convective
c  cores.
c
c  Modified 31/3/01, including new options for the treatment of the
c  convective core. This has been flagged by increasing invers to 1.
c  The treatment is controlled by the various components of imixcr.
c  Also, the source and executable names changed to, say, evolmain.az.d.f
c  (as in the present file), to get away from a steadily increasing number
c  of n-s.
c
c  Modified 2/8/02, including possibility of setting mesh point at
c  boundary of convection zone, flagged by istrtc.
c
c  Modified 5/8/02, including option for helium burning, flagged by
c  iheccs. Also, substantial changes in storage of abundances in 
c  common/compvr/ and output in gong files. As a result, inverb 
c  increased to 2.
c
c  Modified 29/8/02 in s/r mshstr, allowing effects of several discontinuities
c  in composition, etc., and modifying the definition of discontinuities
c  in s/r resxdr. Also reduced effects of very small convective cores
c  on mesh.
c
c  Modified 2/9/02, correcting resetting of central-boundary quantities
c  in case of failure to converge.
c
c  Modified 2/10/02, inserting possibility of shifting luminosity 
c  in case of very small, or negative, values. This is flagged, 
c  and at present switched on, by setting ishfal.
c  The value alshft of the shift is stored in datout(155).
c
c  Modified 14/10/02, resetting He abundance and Z (with 4He burning)
c  to ensure that Z is never smaller than initial surface value
c
c  Modified 5/6/03, correcting function fdradp by including 
c  erroneously forgotten common/clshft/ alshft. Note that this might
c  have had serious effects on the localization of the convective
c  core, in s/r mixcor
c
c  Modified 23/6/03, increasing cut-off in s/r shftal from 0.05 to 0.2
c  Use initial undercorrection after reduction in timestep when failure
c  to converge.
c
c  Modified 21/7/04, adding unformatted output of csum to unit iducen
c  if iducen .gt. 0 (default: -1)
c
c
c  ********************************************************************
c
c  Version numbers in icase and ndtmod:
c  ===================================
c
c  The version of the programme is characterized by the version number
c  in the case number icase of the models, given by the term in 10**7
c  in icase. As of 13/8/96, this has been set to 9, and a new version
c  number has been included in ndtmod.
c  The following changes of the version number have been made:
c
c  To flag for rescaling of central boundary condition coefficients,
c  and resetting of y(1,n) and y(4,n),
c  set version number to 2 (modification 4/1/1985).
c
c  To flag for redefinition of sbcfct and albdbc, and inclusion 
c  in datmod, set version number to 3 (modification 18/1/1988)
c
c  To flag for the use of consistent numerical constants,
c  set version number to 4 (modification 15/3/1988)
c
c  To flag for the possible use of different treatments of the
c  temperature gradient in the atmosphere (flagged by icsrad)
c  set version number to 5 (modification 5/6/1989)
c  
c  To flag for the changed definition of modeeq (to indicate use of
c  Coulomb terms in equation of state), and for the possible inclusion
c  of ivreng in icase, set version number to 6 (modification 14/5/90).
c
c  To flag for extensive modifications starting July 91, involving
c  CNO etc, set version number to 7 (modification 2/8/91).
c
c  Set version number in icase to 9, thus reaching the end of the
c  useful range. Add new version number invers in ndtmod(41),
c  defined as
c
c     invers = 100*inverd + inverb
c
c  Here inverb flags the basic (nondiffusive) structure of the code
c  Initial value of inverb = 0.
c  inverd flags the version of the diffusion treatment.
c  For the present (nondiffusive) version of the code, inverd = 0.
c  (Modification added 13/8/96).
c
c  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  Notes on storage of composition variables.
c
c  Variables considered and storage are controlled by the
c  parameters iche3 and icnocs determining the treatment
c  of He3 and the CNO cycle, respectively. In addition,
c  storage in central boundary condition array azt
c  depends on value of ihe3bc.
c
c  (As of 18/1/96) iche3 is given by iche3 = iche30+10*iche31
c  Here iche31 .gt. 0 is used to flag for using constant-coefficient
c  evolution of He3 abundance, while iche30 controls whether
c  He3 is taken to be in equilibrium. Storage, then depends on
c  iche30 (note that it is not meaningful to have iche30 = 2 for
c  iche31 .gt. 0)
c
c  For convenience in the description, here we introduce array
c  xcno containing CNO abundances. This is not in fact used
c  in the program.
c
c  For icnocs = 0: no CNO elements are considered.
c  For icnocs = 1 and icnocs = 2: xcno(1) = X(N14).
c  For icnocs = 3: xcno(1) = X(C13), xcno(2) = X(N14).
c  For icnocs = 4: xcno(1) = X(C12), xcno(2) = X(C13), 
c                  xcno(3) = X(N14).
c  The variable ispcno (in common/cmpstr/) gives the number of
c  CNO elements considered.
c
c  In array y(i,n): only store active variables, i.e. variables 
c  included in the solution. For simplicity, use y(i)
c  instead of y(i,n) in description. In all cases y(5) = X.
c  For iche03 = 2, icnocs = 0, no other variables considered.
c  For iche03 = 1, icnocs = 0: y(6) = X(He3).
c  For iche03 = 2, icnocs .gt. 0: y(5 + i) = xcno(i), i = 1, ispcno.
c  For iche03 = 1, icnocs .gt. 0: y(6) = X(He3),
c                                y(6 + i) = xcno(i), i = 1, ispcno.
c
c  In array cvr(i,n) in common/compvr/: store more extensive set of 
c  composition variables, for use in output. Storage 
c  only depends on icnocs. Variables are set in the array
c  cvr during initialization, depending on input model
c  and input parameters, and later in s/r rhs and s/r bciter.
c  For simplicity use cvr(i) instead of cvr(i,n) in description.
c  In all cases cvr(1) = X, cvr(2) = Y, cvr(3) = X(He3).
c  For icnocs = 1: cvr(4) = X(N14), cvr(5) = X(O16).
c  For icnocs = 2: cvr(4) = X(C12), cvr(5) = X(C13),
c                  cvr(6) = X(N14), cvr(7) = X(O16).
c  For icnocs = 3: cvr(4) = X(C12), cvr(5) = X(C13),
c                  cvr(6) = X(N14).
c  For icnocs = 4: cvr(4) = X(C12), cvr(5) = X(C13),
c                  cvr(6) = X(N14), cvr(7) = X(O16).
c
c  In arrays aztst and axst used for external storage of central
c  values and second derivatives. Storage only depends on icnocs.
c  For icnocs =    0: aztst(3) = X, aztst(4) = X(He3).
c  For icnocs .gt. 0: aztst(3) = X, aztst(4) = X(He3),
c                     aztst(k + 4) = xcno(k), k = 1, ispcno.
c
c  In arrays azt and ax used for internal storage of central
c  values and second derivatives in s/r bciter. 
c  Storage depends on ihe3bc and icnocs.
c  For ihe3bc = 1, icnocs =    0: azt(3) = X only.
c  For ihe3bc = 0, icnocs =    0: azt(3) = X, azt(4) = X(He3).
c  For icnocs .gt. 0: azt(3) = X, azt(4) = X(He3),
c  For ihe3bc = 1, icnocs .gt. 0: azt(3) = X,
c                     azt(k + 3) = xcno(k), k = 1, ispcno.
c  For ihe3bc = 0, icnocs .gt. 0: azt(3) = X, azt(4) = X(He3),
c                     azt(k + 4) = xcno(k), k = 1, ispcno.
c
c  Finally initial values of central CNO abundances, stored as 
c  in xcno, are in xcnoc in common/bccomp/.
c  
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
      implicit double precision (a-h,o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
c
      parameter(naztmx = nspcmx+3, iy = 2*ivarmx, 
     *  kaxst = 3+naztmx, nbmax = nspcmx+2)
      parameter(nitmax=50)
c
      character tail*5, cntrd*80, ctsxin*5, namez*5
      integer v
      logical time0,resrc,noder,nosd,notd,comout,dtest,skipt,norct,
     *  lastmd, concor, nscfil
      character*16 pmodu3
      character*80 file, filess
      dimension ea(iy,3),iea(iy),v(ivarmx),zk(1),ap(1),aq(1),
     .  thetad(ivarmx),dwf(40,1),
     .  axbc(naztmx),aztbc(naztmx),
     .  datmod(nrdtmd),ndtmod(nidtmd),datout(nrdtmd),ndtout(nidtmd),
     .  dtopfg(4),derlr(2,2,5),pmodu3(2),qqbcpr(5),xnum(10)
      common/yarr/ y(iy,1)
      common/xarr/ x(1)
      common/rhcn/ a1,a2,a3,a4,zdummy,nvar,ifwrt,irhtst,inentc,
     *  ii1,ii2,ii3,icomp
      common/cmtime/ age, time0, lastmd
      common/ln10/ amm,amm2,amm3
      common/bccn/ b1,b2,b3,b4,nb,iveb,icncbc,ihe3bc
      common/bcatms/ tmnbc,tmxbc,sbcfct,flsatm,ntaubc,iopatm,icsrad
      common/qatmos/ iqfit,iqqbc,qqbc(7)
      common/mxlcn/ c1,c2,etac,phc,alfa,tprfct,iconcs,imxlng,iturpr
      common/cnvovs/ icnvos, jcnvos, clcovs, cldovs, alphos, 
     *  rczl, rczlfx, rcnvos, qlcnos
      common/covsec/ icsove, icsovc, alpove, alpovc, 
     *  qmxove, qmxovc, rmxove, rmxovc, imxove, imxovc, 
     *  nmxove, nmxovc
      common/comgrt/ isprot, irotcn, a5, riner0, riner, omgrt0,
     *  omgrot(1)
      common/cdiffu/ idiffus, itbdif, rhob, rmxdif, tlfdif, dtldif,
     *  ctdfmx, rctdif, rctdf1,idfxbc, ismdif, nsmdif, idiffc1
c  size of work: 3046 + 22*nn, or if setting oscillation quantities
c     30*(nn + ntaubc) if this is larger.
      common/work/ dummy(20000)
c  size of sooner: 21*nn
      common/sooner/ soon(13000)
      common/noiter/ iter, ntime, eps, eam
      common/eqscnt/ anz0,anhe0,ihvz,iprrad,ihmins,igndeg
      common/eqstcl/ iradp,mode,dtest,skipt,noder
      common/eqstfl/ fleq
      common/hvabnd/ zab(10),iab
      common/hvcntl/ icnthv,iwrthv,dptst0,dptst1
      common/hvomcl/ iomfll
      common/clshft/ alshft
      common/heavy/ zatmos, zhc, zh(1)
      common/totmss/ am, rscmm
      common/potetc/ chi(125),amz(10),izz(10)
      common/hvname/ namez(10)
      common/opccnt/ xhsopc,tsmn,tstr,rhsmn,rhsmx,timx,rhimn,
     .  rhimx,sigstr,inopc,idmopc
      common/opctcl/ alamop,zatmop,tlmopc,dtmopc,fdgopl,iwdopc,
     *  iopcm1,ifdgop
      common/opctfg/ tlopfg,dtlopf,tlopf1,idopfd
      common/opccof/ iopccf,iopcdm,opcdat(10000)
      common/thetac/ theta(ivarmx)
      common/step/ dt
      common/cxhcnc/ compc(nspcmx)
      common/bccomp/ xhc,xhs,xcnoc(nspcmx)
      common/ksider/ bccoef(nbcprv)
      common/bcprev/ bccofp(nbcprv)
      common/caztax/ azt(nbmax),ax(nbmax)
      common/excf/ am0,al0,am2,al2
      common/logf/ flc /ebcder/ epsc
     *  /rbcder/ rhoc
      common/cntmsh/ wx,epsr,wr,wpt,wxh,wx3,wdgr,
     .  wmshos,drmsos,wmshcc,dqlmsc,qmscmn,qmscmx,
     .  intork,intvnt,iprmsh,icngrm,nnt,iwdgrd,istrtc,istrt1
      common/convpt/ dcs,dcc,nf(6),nl(6),inc,jnc,frcf(6),frcl(6),
     *  xrcf(6),xrcl(6)
      common/convpp/ dcsp,dccp,nfp(6),nlp(6),incp,jncp,
     *  frcfp(6),frclp(6),xrcfp(6),xrclp(6)
      common/cmxcor/ qmxcor,rmxcor,frmxc,nmxcor,nmxcp,qmxcp,rmxcp,
     *  qmxmin,qmxmax
      common/crxcor/ cqc, cqcp, crxmn(nspcmx), crxmnp(nspcmx)
      common/cymixc/ ymix(ivarmx,nnmax)
      common/cmxcnt/ ddrmix, inmixc, imixc0, imixc1, imixc2, imixc3
      common/eqstd/ xii(4),dum(50),pt(20),dum01(16)
      common/ccoulm/ epsdmu, icoulm, iclmit, iclmsd, nitdmu, epssdr
      common/nmbmsh/ nn
      common/rnratd/ al(10,krnrmx),norct
      common/he3fdg/ agesh,ifdhe3,iche30,iche31
      common/cnofrc/ fcno, xtlcno
      common/degfct/ thte,iscren
      common/he3eql/ xhe3eq,ieqhe3,iequs
      common/compos/ xseth, yset, xset3, xset12, xset13, xset14, xset16
      common/he3int/ xhe3(4)
      common/compvr/ cvr(icvrmx,1)
      common/compsz/ xzerh, yzer, xzer3, xrz12, xrz13, xrz14, xrz16
      common/cstcvr/ icvh1, icvhe4, icvhe3, icvc12, icvc13, icvn14, 
     *  icvo16
      common/cstyst/ iyh1, iyhe3, iyc12, iyc13, iyn14, iyo16, 
     *  iyche4, iycc12
      common/cstazt/ iah1, iahe3, iac12, iac13, ian14, iao16, 
     *  iache4, iacc12
      common/cmpstr/ ispxx3, ispcno, icvxx3, icvcno, ibcxx3, ibccno,
     *  isphec
      common/engcnt/ xhzlm1, xhzlm2, nxhzer, icnocs, nspec,
     *  iheccs, nspect
      common/engche/ xmxrhe
      common/engfdg/ epsfdg, qepsf1, qepsf2, ifdep1
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3,idgmxc,
     *  idgsdt,idgn75,idgn76,idgn77
      common/cversn/ ivrbcs,ivrrhs,ivreos,ivropc,ivreng
      common/ctabvr/ ivteos
      common/cdgtnr/ idgtnr,idgtdm,dytnrk(2000)
      common/cnvout/ dmm01(20)
      common/dmuder/ dmm03(11)
      common/rnrout/ dmm04(100)
      common/cofile/ nfiles, idsfil(20), file(20), iopen(20), filess(20)
c
c  commons flagging for errors in s/r rhs, bcs or physics routines
c
      common/cdgrhb/ kdgrhb
      common/cdgphs/ kdgeos, kdgopc, kdgeng
c
c  common giving fundamental constants
c
      common /fconst/ pi, ebase, amu, ame, clight, planck, boltzm, 
     *  cgrav, amsun, echar, ergev, syear
c
c  common giving derived constants, and constants used for
c  equation of state
c
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev1,exhm
c
c  dataset designators
c
      common/cdsdsg/ idstrl,idshvz,idszab,idsevl,idssum,idscen,idsnov,
     *  idsoov,idsgng,idsgsm,idstm1,idstm2,iducen
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
c  Note: if input is actually done with the namelist option,
c  some of the lines should probably be merged, to fit within
c  the maximum number of continuation lines.
c
c$nl      namelist/exec/ cntrd, idsexc, istart,
c$nl     *  istdpr, idsoov, idstrl, idsevl, idsefl, idssum,
c$nl     *  idscen, idsnov, idsgng, idsgng,
c$nl     *  idstm1, idstm2,
c$nl     *  am, z, nn,
c$nl     *  age0, nt, agefin, icntsl,
c$nl     *  iterlr, nitlr, rsfin, alsfin, epslr,
c$nl     *  itrial, nmdtrl, ix, xxh, ix3, istcno, isetzh, isthec, 
c$nl     *  nvartr, idatmd, itrdum, ipartr, iptr,
c$nl     *  ihvz, isethv, ianhe0, modeeq, iomfll,
c$nl     *  iwdopc, xhsopc, tsmn, tstr, rhsmn, rhsmx,
c$nl     *  timx, rhimn, rhimx, sigstr, inopc, iopccf,
c$nl     *  ifdgop, fdgopl, tlopfg, dtlopf, tlopf1,
c$nl     *  iopacm, iopatm, patmos, alamop, zatmop,
c$nl     *  fcno, iscren, thte, iche3, agehe3, icnocs, ivreng, iheccs,
c$nl     *  xzer3, xrz12, xrz13, xrz14, xrz16,
c$nl     *  ifdeps, epsfdg, qepsf1, qepsf2,
c$nl     *  iconcs, imxlng, iturpr, tprfct,
c$nl     *  alfa, etac, phc,
c$nl     *  imixcr, ddrmix,
c$nl     *  icnvos, jcnvos, clcovs, cldovs, alphos,
c$nl     *  ntaubc, sbcfct, albdbc, tmnbc, tmxbc, flsatm, iqfit,icsrad, 
c$nl     *  qqbc,
c$nl     *  ihe3bc, idfxbc,
c$nl     *  eps, nit, ucyt, nucy0, nucy, inentr, icncbc, 
c$nl     *  thetad,
c$nl     *  itndgn, itnprt,
c$nl     *  istrtc, epsr, wx, wpt, wr, wxh,wx3, wdgr, iwdgrd, 
c$nl     *  wmshos,drmsos,wmshcc,dqlmsc,qmscmn,qmscmx,
c$nl     *  nrpeat, iprmsh,
c$nl     *  dt0, dtmx, dtmn, dymx, xdtden, aldtrt, eta,
c$nl     *  icmout, istosc, ioscpr, ndsmth, icasex, intr, isetos,
c$nl     *  iastr, ifstr, ilstr, ispcpr,
c$nl     *  idgbcs,idgrhs,idgeos,idgopc,idgeng, idgtnr, idgmxc, idgsdt,
c$nl     *  itbcrh, irhtst
c
      equivalence(dummy(1),dwf(1,1)),(bccoef(3),aztbc(1)),
     .  (bccoef(kaxst),axbc(1))
c
      external bcs,rhs,dmpbcs,blstio,bleqst
c
      save
c
c  derivatives for L and R iteration. Written in order
c
c              d ln alpha   d ln X  d ln alpha  d ln X
c              ----------   ------  ----------  ------
c                d ln L     d ln L    d ln R    d ln R
c
c  or similarly, with alpha replaced by sbcfct
c
      data derlr /1.314,    -0.1501,   -4.918,   -0.04,
     *            0.809,    -0.229,    -4.529,   -0.043,
     *            2.723,    -0.152,    -10.97,   -0.046,
     *            2.075,    -0.235,    -10.95,   -0.052,
     *            6.640,    -0.152,    -19.25,   -0.044/
c
c  previous values for first line (from Eggleton models):
c
c..   data derlr /1.24,     -0.148,    -4.91,    -0.01,
c
      data pmodu3 /'  ','output on unit 3'/
      data tail /'    @'/
c
c  initialize data set controls
c
      data istdpp, idstrp, idshvp, idsevp, idsflp, idscnp, idssmp, 
     *     idsnvp, idsovp, idsgnp, idsgsp, idst1p, idst2p, idst3p, 
     *     idszap, inopcp, idopfp, iducnp
     *     /18*-1/
c
c  initialize centralization factors
c
      data thetad /ivarmx*0.5/
c
c  initialize control for opacity read
c
      data initop /0/
c
c  cut-off limits for zero hydrogen abundance
c  (when X .lt. xhzlm1, X is replaced by 1.d-10, when X .lt. xhzlm2
c   nuclear reactions are set to zero).
c
      xhzlm1=1.d-6
      xhzlm2=1.d-7
c
c  initialize flag for beginning of calculation (note that this
c  is reset to one after completion of parameter input)
c
      ibegin=0
c
c  flag for stop in case of no convergence
c
      istcon=0
c
c  maximum number of repeats at a given timestep (e.g., halving timestep)
c
      nrpmax=3
c
c  ******************************************************************
c
c  set version numbers
c
      iversn = 9
      invers=2
      write(istdou,101) iversn, invers
      if(istdpr.gt.0.and.istdou.ne.istdpr) 
     *  write(istdpr,101) iversn, invers
c
c  set format code for model output (changed to -2 on 13/5/92)
c
      iform = -2
c
c  initialize data files to zero
c
      call zero(datmod,nrdtmd)
      call zero(datout,nrdtmd)
      call izero(ndtmod,nidtmd)
      call izero(ndtout,nidtmd)
c
c  set scratch directory (initialized in common/cscrtc/)
c
      call stescr
c
c  flag for hydrogen test
c
      itsxin = 2
c
c  flags for 4He output
c
      iprhe4=0
      dprhe4=0.d0
      nphe41=580
      nphe42=601
      ndphe4=1
c
c  flag for shifting luminosity
c
      ishfal=1
c
c  ******************************************************************
c
c  defaults in exec
c  ****************
c
c  cntrd: string determining which control fields are read.
c      contains some or all of the following 
c      dsn: change dataset designator
c      mod: defines model
c      tri: determine trial model
c      eos: define equation of state
c      opa: define opacity
c      eng: define energy generation
c      con: define mixing length theory parameters, including
c           core mixing
c      ovs: define parameters related to convective overshoot.
c      rot: define parameters related to effects of rotation.
c      bcs: define boundary conditions, including atmosphere
c      int: read controls for integration procedure
c      msh: read controls for mesh stretching
c      tst: read controls for time step
c      out: read controls for output
c      dgn: read controls for diagnostics
c      fields may be separated by, e.g., "." (but  n o t comma
c      or blank).
      cntrd='dsn.mod.tri.eos.opa.eng.con.atm.int.msh.tst.out.dgn'
c
c  if idsexc .gt. 0 immediately read new exec from d/s idsexc
      idsexc=0
c  istart: if istart = 0, stop execution
      istart=1
c
c  dataset designators. implemented on 20/1/84. apart from
c  exceptions noted, old values have been used as defaults.
c
c  istdpr: unit for printed output. Initialized to (system dependent)
c     value in block data blstio
c  idstrl: input of trial model
      idstrl=2
c  idshvz: potential input of heavy-element abundance profile
      idshvz=12
c  idszab: potential input of relative composition of heavy elements
      idszab=14
c  idsevl: output evolution variables
      idsevl=3
c  idsefl: if ifstr = 1 (ilstr = 1) output first (last) model
c          on d/s idsefl
      idsefl=4
c  idssum: output formatted summary
      idssum=10
c  idscen: output values at centre
      idscen=11
c  idsnov: output new oscillation variables
      idsnov=15
c  idsoov: output of old (dogfish) oscillation variables
c          old value was idsoov=1.
      idsoov=16
c  idsgng: output GONG model variables
      idsgng=17
c  idsgsm: output GONG summary
      idsgsm=18
c  idstm1: temporary output during internal iteration for alfa
c          and x0. output first model in sequence
      idstm1=20
c  idstm2: temporary output during internal iteration for alfa
c          and x0. output all models in sequence
      idstm2=21
c  iducen: output of same variables as on idscen, but unformatted
c  (for full precision), if iducen .gt. 0
      iducen=-1
c
c  model mass, in solar units
      am=1.d0
c  heavy element abundance
      z=2.e-2
c  number of mesh points
      nn=201
c  initial age
      age0=0.d0
c  nt: number of time steps. if nt = 0, only initial static
c  model is computed.
      nt=0
c  final age (in years)
      agefin=4.75d9
c  icntsl: if icntsl ge 1, this is continuation of evolution sequence
c  Note: icntsl = 1 forces ipartr = 1.
      icntsl=0
c  iterlr: when iterlr.gt.0  iterate inside  programme for final 
c  rs and ls to get values in rsfin and alsfin, using derivatives in
c  derlr(.,.,iterlr). 
c  for iterlr = 1 and 2 iteration is the normal case, iterating on
c  X0 and alfa.
c  for iterlr = 3 and 4 iteration is for X0 and sbcfct (parameter in
c  surface boundary condition). This is intended for simplified model
c  with full ionization everywhere. Parameters set for albdbc = 6.
c  iterlr = 1 and 3 corresponds to iterating on a complete
c  evolution sequence, iterlr = 2 and 4 on a static model 
c  of the present sun.
c  For iterlr=5, use (old) CM formulation, as implemented by MJM,
c  with mixing length given by distance to surface
c
      iterlr=0
      nitlr=5
      rsfin=6.9599d10
      alsfin=3.8481d33
c  convergence limit
      epslr=1.d-5
c
c  itrial: flag for trial model
c  for itrial .gt. 0 use double precision trial models.
c  itrial = 1: use s/r trialp (probably obsolete)
c         = 2 - 4: read in previously computed model as trial
c           from d/s 2. when nmdtrl .gt. 0 use model no. nmdtrl.
c           if nmdtrl .ge. 1000 use last model on file.
c           for itrial = 3 reset number of mesh points to nn, using
c           same distribution of points.
c           for itrial = 4 redistribute mesh points.
c  for itrial .lt. 0 use single precision trial models, according to
c  abs(itrial) as above
c           
      itrial=2
      nmdtrl=0
c  ix: flag for resetting x in trial model.
c  ix = 1: do not reset x
c     = 2: reset x uniformly to xxh
c     = 3: scale x with constant factor, so that surface value is xxh.
      ix=2
      xxh=0.74564d0
c  when ix3 = 2, fudge he3 age agehe3 is reset to value from datmod,
c  if ipartr=1 or icntsl = 1 or isetos = 1
c  when ix3 = 3, fudge h3 age agehe3 is reset to value of age, provided
c  age .gt. 0.
      ix3=1
c  istcno: determines initial setting of CNO abundances
c  istcno = 1: set CNO abundances from trial model, if available
c  (otherwise use same setting as for istcno = 2)
c  istcno = 2: set CNO abundances from variables xrz12, ...
c  read in (see below, under energy generation section)
      istcno = 2
c  isetzh: controls setting of heavy-element abundance, as a function
c  possibly of position.
c  isetzh = 0: assume Z to be constant, given by input parameter z.
c  isetzh .gt. 0: set Z from call of s/r sethvz, reading data from
c  idshvz and interpolating to given mesh.
c  For isethz = 1: set Z directly from values read in
c  For isetzh = 2: rescale values read in such that surface
c  value is z, as given in parameters.
c  isetzh = -1: Set Z from trial model, either from varying
c  Z (with heavy-element settling) or from datmod
c  isetzh = -2: Set Z from trial model, either from Z as set in 
c  (with heavy-element settling) or from 1 - X - Y (in case with
c  4He burning)
c
c  In all cases, the heavy-element abundance is set into the
c  array zh in common/heavy/, at each meshpoint.
      isetzh = 0
c
c  isthec: determines initial setting of abundances for 4He burning.
c  Defined analogously to istcno.
      isthec = 2
c  nvartr: number of variables on trial model file, for itrial .ge. 2
      nvartr=6
c  idatmd: if idatmd = 1, trial model contains datmod (should be
c  the case for all presently conceivable models).
      idatmd=1
c  itrdum: Dummy input parameter, left for consistency with
c  old input files (may be of use later. Used to be itrleq).
      itrdum=0
c  ipartr: if ipartr = 1, reset all available parameters (including
c  age, mixing length etc) from values in datmod read in with trial
c  Note: when isetos = 1, or icntsl = 1, ipartr is forced to be 1.
      ipartr=0
c  if iptr .gt. 0, print trial solution using s/r prtsol, with
c  iptr steps.
      iptr=0
c
c  flag for heavy element treatment.
c  ihvz = 0: assume heavy elements to be fully ionized everywhere
c  ihvz = 1: complete treatment of ionization of C and O; first level
c            of ionization of Fe. Reset abundances to fit results of
c            full treatment of all 10 elements considered.
c  ihvz = 2: first level of ionization of all 10 elements.
c  ihvz = 3: complete treatment of C, N and O; first level of
c            ionization of remaining elements.
c  ihvz = 4: complete treatment of ionization of all 10 elements
c            included.
      ihvz=1
c  isethv: determines the setting up of abundances of heavy elements.
c  isethv = 0: read abundances from unit idszab
c  isethv = 1: use Ross & Aller C, O and Fe abundances
c  isethv = 2: use abundances set in block data bleqst (currently
c  Ross & Aller values for 10 elements).
c  isethv = 3: force setting of az = 16, anh0 = 0.5, anhe0 = 8.
c      these are the values to be used for the simple GONG models.
c  isethv = 4: read number densities, relative to hydrogen,
c              from unit idszab
c  note: this, together with the resetting of the heavy element
c  abundance in s/r hviona is rather a mess, and should be
c  cleared up.
c  see also resetting of heavy element parameters after statement
c  with label 33, below.
      isethv=2
c  ianhe0: when ianhe0 = 1 reset anhe0 to 6 (the value before
c          november 1982). otherwise keep value set in s/r setcns
      ianhe0=0
c  modeeq: From May 1990 used to flag for equation of state employed.
c  modeeq = 0: use EFF equation of state
c  modeeq = 1: use consistent Coulomb treatment in EFF equation of state
c  modeeq = 2: use inconsistent Coulomb treatment in EFF equation of 
c  state (including direct effect on pressure and enthalpy, but not in
c  ionization) 
c  Note: the last two options only has an effect when programme 
c  including Coulomb version eqstf is used.
      modeeq=0
c  iomfll: switch for statistical weight of fully ionized heavy elements
c  when iomfll = 0 set this to 15 (erroneously) to match situation in
c  old version of programme.
c  if iomfll = 1 use correct value.
      iomfll=0
c
c  opacity controls
c
c  When iwdopc = 1 use w. dappen opacities
c  When iwdopc .ge. 8 use WD/GH routines for interpolation 
c  in OPAL and Kurucz tables. Here
c  iwdopc = iwdop0 + 10*iwdop1, where
c    iwdop0 = 9: use minimum-norm interpolation (as in previous versions)
c    iwdop0 = 8: use birational splines
c
c    iwdop1 = 0: Do not include electron conduction
c    iwdop1 = 1: Include electron conduction
c  Note that iwdop1 .gt. 0 can only be used in versions 9 and later.
c  Otherwise the value of iwdopc determines the interpolation in X:
c  iwdopc = 0: use old scheme with parabolic interpolation.
c    Note: the cluster of points is chosen to be as close as possible
c    to target x. This leads to switch of cluster in the middle of
c    mesh interval, and hence to discontinuous interpolating function.
c    A little unfortunate!
c  iwdopc = 2: use linear interpolation
c  iwdopc = 3: use 4-point Lagrangian interpolation.
      iwdopc=0
c  xhsopc: value of X in outer region for spline opacity interpolation
      xhsopc=0.733078
c  tsmn, tstr, rhsmn, rhsmx, timx, rhimn: parameters controlling
c  regions in spline interpolation
      tsmn=3.5
      tstr=6
      rhsmn=-10
      rhsmx=0
      timx=7.3
      rhimn=-3
      rhimx=3
c  sigstr: stretching parameter in spline interpolation (note: must be
c  negative).
      sigstr=-5
c  inopc: unit number for opacity tables.
      inopc=13
c  size of internal storage for opacity table. This is set up in
c  common /opccof/
      iopccf=10000
c  ifdgop, fdgopl, tlopfg, dtlopf,tlopf1: controls for opacity fudge
c          in opacity routine opact.
c          if ifdgop = 1 fdgopl is added to log(opacity)
c          if ifdgop = 2 log(opacity) is modified with gaussian with
c          maximum fdgopl, centred at tlopfg and with width dtlopf.
c          if ifdgop = 3 log(opacity) is modified with gaussian with
c          maximum fdgopl and with width dtlopf, for log t .lt. tlopf1 
c          or log t .gt. tlopfg, and with constant fdgopl for
c          tlopf1 .le. log t .le. tlopfg.
c          if ifdgop = 4, read in opacity corrections from file,
c          unit idopfd, in the form (log T, delta log kappa)
      ifdgop=0
      fdgopl=0
      tlopfg=6.3
      dtlopf=0.3
      tlopf1=5
      idopfd=19
c  controls for opacity in atmosphere. to use special atmospheric
c  opacity iopacm must be 1 or 2, corresponding to w. d.
c  opacity at z = zatmop (in common/opctcl/; initialized to 0.02),
c  and auer et al opacity at z = zatmop, respectively.
c  when iopatm .ne. .1 make smooth matching between atmospheric and
c  interior opacity, at log t = tlmopc..
c  when iopatm = 1 use atmospheric opacity in s/r atmos only.
      iopacm=0
      iopatm=0
c  pressure at teff in present sun. if patmos .gt. 0 set opacity
c  fudge alamop in atmosphere so that pressure gets correct value.
      patmos=0
      alamop=1
c  zatmop: heavy element abundance in atmosphere, for opacity fudge.
      zatmop=0.02
c
c  energy generation parameters
c  fcno: fraction (by mass) of the heavy elements consisting of
c  C, N and O. Used only when icnocs = 0, for simple equilibrium
c  treatment of CNO cycle.
      fcno=0.28
c  iscren: controls electron screening 
c  iscren = 0: no electron screening 
c  Otherwise iscren = iscnuc+ 10*iscbe7.
c  iscnuc controls screening of nuclear reactions:
c  iscnuc = 1: use weak screening in original formulation, with
c     thte = thtec (as given in input)
c  iscnuc = 2: use weak screening, with consistent thte
c  iscnuc = 3: use intermediate screening, as done by Bahcall.
c  iscnuc = 4: use intermediate screening in Mitler formulation,
c     as used by Turck-Chieze.
c  iscbe7 controls screening of electron capture in Be7:
c  iscbe7 = 0: compute screening as Bahcall and Moeller
c     (ApJ, vol. 155, 511)
c  iscbe7 = 1: use approximate screening factor = 1.2, as
c     suggested by Bahcall
      iscren=1
c  thte: parameter in simple electron screening formulation
      thte=1
c  iche3: determines treatment of He3
c  iche3 = iche30 + 10*iche31
c  iche30 = 1: He3 is not in equilibrium
c  iche30 = 2: He3 is in equilibrium
c  iche31 = 0: use normal evolution equations
c  iche31 = 1: use backsubstitution with constant-coefficient
c     solution
      iche3=1
c  agehe3: age (in years) for initial He3 abundance fudge.
c  If agehe3 le 0, no fudge is used for initial He3
      agehe3=5.d7
c  icnocs: parameter determining treatment of CNO cycle.
c  icnocs = 0: Assume CNO cycle to be in equilibrium.
c  icnocs = 1: Follow conversion between O16 and N14, ignore
c  distribution amongst C12, C13 and N14.
c  icnocs = 2: Follow conversion between O16 and N14, include
c  distribution amongst C12, C13 and N14.
c  icnocs = 3: Ignore conversion between O16 and N14, include
c  conversion amongst C12, C13 and N14.
c  icnocs = 4: Include both conversion between O16 and N14 and
c  conversion amongst C12, C13 and N14.
      icnocs = 0
c  ivreng: parameter for choosing nuclear parameters
c  ivreng = 3: Use Parker (1986) parameters
c  ivreng = 4: Use Bahcall (1992) parameters 
c  ivreng = 5: Reaction constants from J.N. Bahcall routine energy,
c              from May 1992, with He3 + He4 and Be7 + H1 modified
c              to correspond to Dar and Shaviv
c  ivreng = 6: Original data from Gabriel (1990) (see file gab.coef.in)
c              have been replaced by values from
c              Bahcall & Pinsonneault (1994) where available.
c  ivreng = 7: Data from Adelberger et al.
c              (1998; Rev. Mod. Phys. 70, 1265)
c  ivreng = 8: Reaction rates computed with NACRE data
c              (Angulo et al. 1999, Nuclear Physics, A 656, 3-183)
c
      ivreng = 3
c  iheccs: flag for inclusion of 4He and 12C burning.
c  iheccs = 1: include both reactions
      iheccs = 0
c  xzer3: initial abundance of He3 
c  If xzer3 .lt. 0, use He3 abundance as read in from input model
c  unless fudge is used.
      xzer3 = 0
c  xrz12, xrz13, xrz14, xrz16: Initial abundances of C12, C13, N14
c  and O16, relative to total heavy element abundance Z,
c  for treatment of CNO cycle.
c#ai# Note: currently consistent treatment of thermodynamics
c  taking CNO conversion into account has not been implemented.
c  Note: in the special case of icnocs = 1, xrz14 is replaced
c  by fcno*A(N14)/acno, for consistency with old programme.
      xrz12 = 0.2254
      xrz13 = 0.
      xrz14 = 0.0549
      xrz16 = 0.4987
c  ifdeps, epsfdg, qepsf1, qepsf2: parameters for fudging energy
c  generation
c  ifdeps = 1: Increase epsilon by epsfdg for qepsf1 le (m/M) le qepsf2
c  ifdeps = 2: Increase epsilon by epsfdg*(T/Tc)
c              for qepsf1 le (m/M) le qepsf2
c  If ifdeps lt 0, the fudge corresponding to abs(ifdeps) is used
c  for the initial model, but not for the subsequent ones. This
c  is intended for starting pre-main-sequence evolution.
      ifdeps = 0
      epsfdg = 0
      qepsf1 = 0
      qepsf2 = 0
c
c  convection parameters
c  iconcs: determines treatment of convection
c  iconcs = 0: JC-D treatment of mixing-length theory, including
c  Cambridge fudge
c  iconcs = 1: MJM parametrized treatment
c  iconcs = 2: MJM version of old Canuto and Mazittelli formulation
c  iconcs = 3: MJM version of new Canuto and Mazittelli formulation
      iconcs=0
c  imxlng (relevant for MJM formulations): determines setting of 
c  mixing length
c  imxlng = 0: Mixing length proportional to pressure scale height
c  imxlng = 1: Mixing length proportional to depth
c  imxlng = 2: Mixing length varies according to f (???)
      imxlng=0
c  iturpr: controls treatment of turbulent pressure
c  iturpr = 0: Ignore turbulent pressure
c  iturpr = 1: Include turbulent pressure according to RT formulation.
      iturpr = 0
c  tprfct: Factor to be applied in turbulent-pressure calculation
      tprfct = 1
c  alfa: mixing length parameter
      alfa=1.549d0
c  etac: For iconcs = 0, parameter eta in mixing length theory
c        For iconcs = 1, parameter m in MJM parametrized formulation
      etac=0.15713484026d0
c  phc: For iconcs = 0 parameter PHIc in mixing length theory
c        For iconcs = 1 parameter beta in MJM parametrized formulation
      phc=2
c  Note: if etac .lt. 0 and iconcs = 0, use the fudge formulation of 
c  C-D (Cambridge, 1986)
c  Here the quantity alfa_C is set to alfa**(-0.8), and the
c  quantity beta is set to phc, where alfa and phc are the input
c  parameters. The peculiar definition of alfa_C in terms of alfa
c  is made in order to be able to use the same derivatives in
c  the iteration for the surface radius and luminosity.
C
c  imixcr: determines treatment of mixed core.
c  imixcr = imixc0 + 10*imixc1 + 100*imixc2 + 1000*imixc3
c  imixc0 = 1: set instability based on given log f (used before
c  22/5/93)
c  imixc0 = 2: set instability based on given log p 
c  imixc1 = 0: mixing in s/r rhs largely based on whether region is
c  convective or not
c  imixc1 = 1: base mixing in convective core in s/r rhs on 
c  q .le. qmxcor, as set in s/r mixcor, consistently with iteration 
c  for core composition. Also, test for oscillating core and
c  take appropriate action in that case.
c  imixc2 = 0: Use interpolation to unstable limit in s/r mixcor as
c  default (implementation before 28/12/00)
c  imixc2 = 1: Set unstable limit from extrapolation based on composition
c  outside core.
c  imixc3 = 0: No internal iteration for limit of convective core in
c  s/r cmpcvc (implementation before 28/12/00)
c  imixc3 = 1: Internal iteration using back-substitution
c  imixc3 = 2: Internal secant iteration
      imixcr = 1
c  ddrmix: critical difference between adiabatic and radiative
c  gradient for determining extent of mixed core
      ddrmix=0
c
c  Parameters for convective overshoot
c  icnvos: when 1, include overshoot from base of convective envelope
      icnvos = 0
c  jcnvos: additional parameter for overshoot. 
c  For icnvos = 1:
c    jcnvos .le. 1: Use original formulation
c    jcnvos  =   2: Set smooth transition to radiative gradient
c    jcnvos  =   3: Use MJM formulation
c  For icnvos = 2:
c    jcnvos = 1: Use simple tanh formulation, no subadiabatic region
c                in convection zone
c    jcnvos = 2: Use iterated solution that allows subadiabatic region
c                in convection zone (i.e., convectively unstable region)
      jcnvos = 0
c  clcovs: parameter determining width of transition at base of
c  overshoot region. For icnvos = 1 this is curly C in MJM formulation.
c  For icnvos = 2 (MR formulation) simply determines
c  width of tanh transition, in units of surface radius.
      clcovs = 1.d5
c  cldovs: parameter determining extent of overshoot region
c  For icnvos = 1 this is curly D in MJM formulation, i.e.,
c  maximum difference between radiative and adiabatic gradient
c  with overshoot. For icnvos = 2 it is simply the distance, in
c  units of R, below the unstable limit of the centre of the
c  transition.
      cldovs = 0.1
c  alphos: maximum allowed overshoot, in units of pressure scale
c  height at base of convection zone.
      alphos = 1
c  (Note: with overshoot, the parameter wmshos in the mesh stretching
c  group should also be set, to increase number of points in
c  overshoot region).
c  Overshoot parameters added 5/1/00.
c  icsove: flag for overshoot from convective envelope.
c  icsove = 1: Mix composition, set adiabatic temperature gradient
c  icsove = 2: Mix composition, do not change temperature gradient
      icsove = 0
c  icsovc: flag for overshoot from convective core, defined in
c  similar manner.
      icsovc = 0
c  alpove: Overshoot distance from envelope, in units of pressure
c  scale height at convection-zone boundary
      alpove = 0
c  alpovc: Overshoot distance from core, in units of pressure
c  scale height at convection-zone boundary or radius of core,
c  whichever is smaller.
      alpovc = 0
c
c  Parameters for rotation:
c  isprot: flag for inclusion of spherically symmetric effect of
c  rotation.
c  isprot = 0: No rotation
c  isprot = 1: Include rotation
      isprot=0
c  velrot: Rotational velocity. If velrot .le. 1, velrot is assumed to
c  give angular velocity, in sec**(-1). If velrot .gt. 1, velrot
c  is assumed to give equatorial linear velocity, in cm/sec.
      velrot=0
c
c  quantities for atmosphere
c  default are values used for calculating model 1 of c-d82.
c  ntaubc: number of mesh points in atmosphere
c  if ntaubc. le. 1, use simple pressure surface condition
      ntaubc=41
c  sbcfct: parameter for surface pressure condition 
c  for ntau .le. 1 use simple condition p*kappa = sbcfct*G*m/r**2
c  at point where integration stops (normally, for albdbc = 1,
c  at point where t = teff).
c  if ntau .gt. 1, pressure condition on top of atmosphere is 
c  p*kappa = 2*sbcfct*tau*G*M/R**2
      sbcfct=1
c  albdbc:  parameter in luminosity surface boundary condition
c  only applicable when ntau .le. 1. then condition on L
c  is written L = albdbc*4*pi*sigma*R**2*Teff**4
      albdbc=1
c  tmnbc, tmxbc:  minimum and maximum tau in atmosphere. Note:
c  tmxbc should correspond to where temperature is equal to the
c  effective temperature for the T(tau) relation used.
      tmnbc=1.d-4
      tmxbc=0.4116982d0
c  flsatm: trial log(f) for surface atmospheric iteration  
      flsatm=-10.  
c  iqfit: number of parameters for fit to T-tau relation.
c  For iqfit .lt. 0, use RT fitting formulas. In this case
c  iqfit = -(iqfit0 + 100*iqfit1) where iqfit1 describes various options
c  (see source code for qalphafit-art.d.f for details)
      iqfit=5
c  icsrad: parameter describing match between interior and atmosphere
c  (for ntaubc .gt. 0).
c  icsrad = 0: make no smooth transition to diffusion approximation
c  in interior.
c  icsrad = 1: make smooth transition in temperature gradient to
c  diffusion approximation in interior.
      icsrad=0
c  qqbc(1 - iqfit): parameters for fit to T-tau relation,
c  when iqfit .gt. 0.
      qqbc(1)=1.036
      qqbc(2)=-0.3134
      qqbc(3)=2.448
      qqbc(4)=-0.2959
      qqbc(5)=30.
      qqbc(6)=0.
      qqbc(7)=0.
c  ihe3bc: determines treatment of He3 in central boundary condition.
c  ihe3bc = 0: He3 not assumed to be in equilibrium at centre
c  ihe3bc = 1: He3 assumed to be in equilibrium at centre
      ihe3bc = 1
c  idfxbc: (only for calculations with diffusion) determines
c  central boundary condition on diffusion velocity of hydrogen
c  idfxbc = 0: set zero velocity condition
c  idfxbc = 1: set zero gradient condition
      idfxbc=0
c
c  integration controls.
c  eps: convergence criterion
      eps=1.d-6
c  nit: maximum number of iterations.
      nit=10
c  ucyt, nucy0, nucy: ucyt is undercorrection factor applied
c  for the first nucy0 iterations at the first time step, and
c  the first nucy iterations for subsequent time steps.
      ucyt=0.5
      nucy0=0
      nucy=0
c  inentr: when inentr = 1, ignore entropy terms in energy equation
      inentr = 0
c  icncbc: parameter controlling central boundary condition.
c  passed on to boundary condition routines.
c  (precise meaning has to be specified later).
      icncbc = 0
c  theta(1 - nvar): centralization parameters for time discretization
c  All variables (thetad(1 - ivarmx)) initialized to 0.5 in data
c  statement above.
      thetad(4)=1
      thetad(6)=1
c  itndgn: for itndgn = 1 print extended diagnostics from iteration,
c     otherwise only limited dagnostics.
      itndgn=0
c  itnprt: if itnprt .gt. 0, print solution with step itnprt
c  after each iteration.
      itnprt=0
c
c  parameters for stretching
c  istrtc: if istrtc .ge. 1, do stretching
c  note that no stretching is done if only initial model is
c  calculated.
c  To control aspects of mesh-setting define istrtc as
c  istrtc = istrt0 + 10*istrt1 +100*istrt2
c  If istrt0 = 1 and istrt1 .ge. 1, also reset heavy-element composition
c  in cvr in common/compvr/
c  If istrt2 .gt.0 also reset mesh near lower boundary of convective
c  envelope
      istrtc=0
c  epsr: parameter determining position of innermost mesh point
      epsr=1.5e-1
c  wx, wpt, wr, wxh, wx3, wdgr: parameters controlling mesh distribution
      wx=1
      wpt=3
      wr=3
      wx3=0.3
c  wxh: weight for gradient in hydrogen abundance. Added 8/6/92,
c  and initialized to zero for compatibility with previous versions
      wxh=0
      wdgr=1
c  iwdgrd: if iwdgrd = 1 use old weighting of term in second
c     derivative of superadiabatic gradient.
c          if iwdgrd = 2 use new formulation.  
c     if iwdgrd = 3, switch off term from second derivative.
      iwdgrd=1 
c  wmshos: weight for adding meshpoints in overshoot region at base
c  of convective envelope, or at base of convective envelope
      wmshos = 0
c  drmsos: extent, in r/R, of region where extra mesh points are
c  added with overshoot
      drmsos=0.05
c  wmshcc: weight for adding meshpoints at boundary of convective core
      wmshcc = 0.
c  dqlmsc: extent, in log10(m/M), of region where extra mesh points are
c  added at convective core
      dqlmsc=0.08
c  qmscmn, qmscmx: limits in m/M of region for adding extra points, 
c  in connection with mixed core. Note that these are normally set
c  in mshstret.
      qmscmn=0
      qmscmx=0
c  nrpeat: for first stretching in evolution sequence repeat 
c  nrpeat times unless this is later  iteration for rs and ls
      nrpeat=1
c  iprmsh: level of diagnostics in mesh stretching
c  iprmsh = 0: no diagnostics.
c  iprmsh .lt. 0: print stretch variable.
c  iprmsh .gt. 0: print stretch variable, and in addition details
c     at every iprmsh-th point.
      iprmsh=-1
c
c  parameters for time step determination
c  dt0: initial time step (seconds)
c  if dt0 is greater than the time remaining from the initial
c  age to agefin, the initial timestep is set to the remaining
c  time instead
      dt0=2.d14
c  dtmx: maximum time step (seconds)
      dtmx=2.d16
c  dtmn: minimum time step (seconds)
      dtmn=2.d13
c  dymx: allowed change in log p, log T, etc.
      dymx=1.d-1
c  xdtden: test on delta X/(X + xdtden)
      xdtden=0.1
c  aldtrt: for luminosity, limit (change in log l)*factl,  
c     where factl = (l/ls)*(1+aldtrt)/(aldtrt+l/ls)
      aldtrt=0 
c  eta: time extrapolation parameter, for s/r pushy
      eta=1.
c  dtfcrd, dtceps: parameters to reduce effect on setting of time step
c  of changes in narrow regions, for hydrogen abundance and luminocity.
c  Applied when dtceps .gt. 0.
c  dtfcrd defines the smallest value of the reduction factor.
c  See source for setdt for details
      dtfcrd=0.1
      dtceps=0.d0
c
c  in addition programme may input opacity data from d/s inopc
c  (default inopc=13) when using spline interpolation, or from
c  d/s 9 when using w. dappen opacity routines.
c  d/s idszab (default 14) may be used to read relative composition
c  of heavy elements
c
c  output control.
c  icmout: if icmout = 1, make complete printed output, with
c  s/r prtsol. Note: complete output is also obtained for
c  istosc = 1, or ioscpr = 1.
      icmout=0
c  when istosc = 1, set oscillation variables for all models (unless
c  iterating for rs and ls).
c  old oscillation variables (in dogfish format) output on d/s idsoov,
c  new oscillation variables on d/s idsnov.
      istosc=0
c  when ioscpr = 1 set oscillation variables as above, but only for
c  final model of present sun.
      ioscpr=1
c  when ndsmth .gt. 0 smooth transition between interior and atmosphere
c  for new oscillation variables, over range 2*ndsmth
      ndsmth=0
c  icasex: controls setting of fractional radius r/R in subroutine
c  prtsol and newvar. 
c  for icasex .ne. 1 use old (before 7/3/90) formulation, which leads
c  to serious accumulation of rounding errors near surface.
c  for icasex = 1 use r/R = 10.**(y(1,n) - y(1,1))
      icasex=0
c  intr: solution is printed with step intr.
      intr=5
c  for isetos = 1 set oscillation quantities from model on file,
c  read in as determined by itrial.
c  for isetos = -1, set complete output (including central values)
c  for model on file, but do not set oscillation quantities.
c  If nt .gt. 1, go through models nmdtrl, ..., nmdtrl + nt -1.
c  Note: forces ipartr = 1.
      isetos=0
c  iastr: for iastr = 1, store all models (except during iteration for
c  for surface R and L) on unit idsevl.
      iastr=1
c  ifstr: if ifstr = 1 output first model on d/s idsefl
      ifstr=0
c  ilstr: if ilstr = 1 output last model on d/s idsefl
      ilstr=0
c  ispcpr: for ispcpr .gt. 0 call user-supplied subroutine spcout
c     for special output.
      ispcpr=0
c
c  diagnostics.
c  idgbcs: for idgbcs = 1, print diagnostics from boundary conditions
      idgbcs=0
c  idgrhs: for idgrhs = 1, print diagnostics from right hand side
      idgrhs=0
c  idgeos: for idgeos = 1, print diagnostics from equation of state
c      (currently not implemented)
      idgeng=0
c  idgopc: for idgopc = 1, print diagnostics from opacity
      idgeng=0
c  idgeng: for idgeng = 1, print diagnostics from energy generation
      idgeng=0
c  idgtnr: for idgtnr .ge. 1, print all changes from tnrkt for the first
c     idgtnr iterations
      idgtnr=0
c  idgmxc: if idgmxc .ge. 1, print details from setting mixing in 
c     convective core
      idgmxc=0
c  idgsdt: if idgsdt .ge. 1, write details of time-step setting to 
c     unit 43
      idgsdt=0
c  idgn75, idgn76, idgn77: if non-zero, output mesh details,
c  models and right-hand-side results, respectively, during iteration,
c  such that latest set is always available.
c  If idgn75 .gt. 1 mesh details are kept for all time steps
c  If idgn77 .gt. 1 information is also stored about alam.
      idgn75=0
      idgn76=0
      idgn77=0
c  itbcrh: if itbcrh = 1 test right hand side routine for initial model
c     if itbcrh = 2, test in addition boundary condition routine.
      itbcrh=0
c  irhtst: if irhtst .gt. 0, print extensive set of data from
c  s/r rhs, etc.
      irhtst=0
c
c
c  end of defaults for /exec/
c
c  ************************************************************
c
c  set dataset assignments
c
      write(istdou,*)
     *   'Data sets needed:'
      write(istdou,*)
     *   'istdpr (if not ',istdou,'): printed output'
      write(istdou,*)
     *   'idstrl (default 2): input trial model'
      write(istdou,*) 
     *  'idshvz (default 12): input of heavy-element abundance profile'
      write(istdou,*) 
     *  'idszab (default 14): input of composition of heavy elements'
      write(istdou,*)
     *   'idsevl (default 3): output of evolution variables'
      write(istdou,*)
     *   '    if iastr = 1'
      write(istdou,*)
     *   'idsefl (default 4): output of evolution variables '
      write(istdou,*)
     *   '    for first or last model if ifstr = 1 or ilstr = 1'
      write(istdou,*)
     *   'inopc (default 13): opacity tables'
      write(istdou,*)
     *   'oscillation output, for isetos = 1, istosc = 1, or'
      write(istdou,*)
     *   '    ioscpr = 1:'
      write(istdou,*)
     *   'idsnov (default 15): new oscillation variables'
      write(istdou,*)
     *   'idsoov (default 16): old oscillation variables'
      write(istdou,*)
     *   'idsgng (default 17): GONG model variables'
      write(istdou,*)
     *   'idsgsm (default 18): GONG model summary'
      write(istdou,*)
     *   'For complete output (oscillation output or icmout = 1):'
      write(istdou,*)
     *   'idssum (default 10):   formatted central values'
      write(istdou,*)
     *   'idscen (default 11): more extensive formatted ',
     *   'central values (csum)'
      write(istdou,*)
     *   'idstm1 (default 20): temporary model file for iterlr = 1'
      write(istdou,*)
     *   'idstm2 (default 21): temporary model file for iterlr = 1'
      write(istdou,*)
     *   'iducen (default -1; 19 suggested): ',
     *   'more extensive unformatted central values, as idscen (cusum)'
c
      call ofiles
c
c                 **********************************
c
c  set constants for equation of state (further setting of nuclear
c  constants and setting up of opacity below, after read of input
c  parameters)
c
      call setcns
      call seteqs
c
c  initialize flag for first run
c
      ifsrun=1
c
c  open data set for scratch input, output
c
      open(91,status='scratch',form='unformatted')
c
      go to 22
c
c  entry for later pass. reset flag for first run
c
   20 ifsrun=0
c
c  read exec
c  *********
c
   22 continue
c
c$nl      read(5,exec,end=90)
c
      write(istdou,*) 'cntrd?'
      write(istdou,*) cntrd
      read(istdin,'(a)',end=90) cntrd
      write(istdou,*) '     -----------  New values  --------'
      write(istdou,*) cntrd
c
c  decode cntrd
c
      ircdsn = index(cntrd,'dsn')
      ircmod = index(cntrd,'mod')
      irctri = index(cntrd,'tri')
      irceos = index(cntrd,'eos')
      ircopa = index(cntrd,'opa')
      irceng = index(cntrd,'eng')
      irccon = index(cntrd,'con')
      ircovs = index(cntrd,'ovs')
      ircrot = index(cntrd,'rot')
      ircbcs = index(cntrd,'bcs')
      ircint = index(cntrd,'int')
      ircmsh = index(cntrd,'msh')
      irctst = index(cntrd,'tst')
      ircout = index(cntrd,'out')
      ircdgn = index(cntrd,'dgn')
c
      write(istdou,*) 'idsexc, istart?'
      write(istdou,*) idsexc, istart
      read(istdin,*,end=90) idsexc, istart
      write(istdou,*) '     -----------  New values  --------'
      write(istdou,*) idsexc, istart
c
c  test for stop
c
      if(istart.eq.0) go to 90
c
c  data set designators
c
      if(ircdsn.gt.0) then
        write(istdou,*) 
     *    'istdpr, idstrl, idshvz, idszab, idsevl, idsefl, idssum?'
        write(istdou,*) 
     *    istdpr, idstrl, idshvz, idszab, idsevl, idsefl, idssum
        read(istdin,*,end=90) 
     *    istdpr, idstrl, idshvz, idszab, idsevl, idsefl, idssum
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) 
     *    istdpr, idstrl, idshvz, idszab, idsevl, idsefl, idssum
        write(istdou,*) 'idscen, idsnov, idsoov, idsgng, idsgsm?'
        write(istdou,*) idscen, idsnov, idsoov, idsgng, idsgsm
        read(istdin,*,end=90) idscen, idsnov, idsoov, idsgng, idsgsm
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) idscen, idsnov, idsoov, idsgng, idsgsm
	write(istdou,*) 'idstm1, idstm2, iducen?'
        write(istdou,*) idstm1, idstm2, iducen
        read(istdin,*,end=90) idstm1, idstm2, iducen
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) idstm1, idstm2, iducen
      end if
c
c  model definition
c
      if(ircmod.gt.0) then
        write(istdou,*) 'am, z, nn?'
        write(istdou,*) am, z, nn
        read(istdin,*,end=90) am, z, nn
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) am, z, nn
        write(istdou,*) 'age0, nt, agefin, icntsl?'
        write(istdou,*) age0, nt, agefin, icntsl
        read(istdin,*,end=90) age0, nt, agefin, icntsl
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) age0, nt, agefin, icntsl
        write(istdou,*) 'iterlr, nitlr, rsfin, alsfin, epslr?'
        write(istdou,*) iterlr, nitlr, rsfin, alsfin, epslr
        read(istdin,*,end=90) iterlr, nitlr, rsfin, alsfin, epslr
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) iterlr, nitlr, rsfin, alsfin, epslr
      end if
c
c  trial model definition
c
      if(irctri.gt.0) then
        write(istdou,*) 
     *    'itrial, nmdtrl, ix, xxh, ix3, istcno, isetzh, isthec?'
        write(istdou,*) 
     *    itrial, nmdtrl, ix, xxh, ix3, istcno, isetzh, isthec
        read(istdin,*,end=90) 
     *    itrial, nmdtrl, ix, xxh, ix3, istcno, isetzh, isthec
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) 
     *    itrial, nmdtrl, ix, xxh, ix3, istcno, isetzh, isthec
        write(istdou,*) 'nvartr, idatmd, itrdum, ipartr, iptr?'
        write(istdou,*) nvartr, idatmd, itrdum, ipartr, iptr
        read(istdin,*,end=90) 
     *    nvartr, idatmd, itrdum, ipartr, iptr
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) nvartr, idatmd, itrdum, ipartr, iptr
      end if
c
c  equation of state definition
c
      if(irceos.gt.0) then
        write(istdou,*) 'ihvz, isethv, ianhe0, modeeq, iomfll?'
        write(istdou,*) ihvz, isethv, ianhe0, modeeq, iomfll
        read(istdin,*,end=90) ihvz, isethv, ianhe0, modeeq, iomfll
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) ihvz, isethv, ianhe0, modeeq, iomfll
      end if
c
c  opacity definition
c
      if(ircopa.gt.0) then
        write(istdou,*) 'iwdopc, xhsopc, tsmn, tstr, rhsmn, rhsmx?'
        write(istdou,*) iwdopc, xhsopc, tsmn, tstr, rhsmn, rhsmx
        read(istdin,*,end=90) iwdopc, xhsopc, tsmn, tstr, rhsmn, rhsmx
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) iwdopc, xhsopc, tsmn, tstr, rhsmn, rhsmx
        write(istdou,*) 'timx, rhimn, rhimx, sigstr, inopc, iopccf?'
        write(istdou,*) timx, rhimn, rhimx, sigstr, inopc, iopccf
        read(istdin,*,end=90) timx, rhimn, rhimx, sigstr, inopc, iopccf
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) timx, rhimn, rhimx, sigstr, inopc, iopccf
        write(istdou,*) 
     *    'ifdgop, fdgopl, tlopfg, dtlopf, tlopf1, idopfd?'
        write(istdou,*) ifdgop, fdgopl, tlopfg, dtlopf, tlopf1, idopfd
        read(istdin,*,end=90) 
     *    ifdgop, fdgopl, tlopfg, dtlopf, tlopf1, idopfd
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) ifdgop, fdgopl, tlopfg, dtlopf, tlopf1, idopfd
        write(istdou,*) 'iopacm, iopatm, patmos, alamop, zatmop?'
        write(istdou,*) iopacm, iopatm, patmos, alamop, zatmop
        read(istdin,*,end=90) iopacm, iopatm, patmos, alamop, zatmop
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) iopacm, iopatm, patmos, alamop, zatmop
      end if
c
c  energy generation definition
c
      if(irceng.gt.0) then
        write(istdou,*) 
     *    'fcno, iscren, thte, iche3, agehe3, icnocs, ivreng, iheccs?'
        write(istdou,*) 
     *    fcno, iscren, thte, iche3, agehe3, icnocs, ivreng, iheccs
        read(istdin,*,end=90) 
     *    fcno, iscren, thte, iche3, agehe3, icnocs, ivreng, iheccs
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) 
     *    fcno, iscren, thte, iche3, agehe3, icnocs, ivreng, iheccs
        write(istdou,*) 'xzer3, xrz12, xrz13, xrz14, xrz16?'
        write(istdou,*) xzer3, xrz12, xrz13, xrz14, xrz16
        read(istdin,*,end=90) xzer3, xrz12, xrz13, xrz14, xrz16
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) xzer3, xrz12, xrz13, xrz14, xrz16
	write(istdou,*) 'ifdeps, epsfdg, qepsf1, qepsf2?'
	write(istdou,*) ifdeps, epsfdg, qepsf1, qepsf2
	read(istdin,*) ifdeps, epsfdg, qepsf1, qepsf2
        write(istdou,*) '     -----------  New values  --------'
	write(istdou,*) ifdeps, epsfdg, qepsf1, qepsf2
      end if
c
c  mixing length theory definition
c
      if(irccon.gt.0) then
        write(istdou,*) 'iconcs, imxlng, iturpr, tprfct?'
        write(istdou,*)  iconcs, imxlng, iturpr, tprfct
        read(istdin,*,end=90) iconcs, imxlng, iturpr, tprfct
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*)  iconcs, imxlng, iturpr, tprfct
        write(istdou,*) 'alfa, etac, phc?'
        write(istdou,*) alfa, etac, phc
        read(istdin,*,end=90) alfa, etac, phc
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) alfa, etac, phc
        write(istdou,*) 'imixcr, ddrmix?'
        write(istdou,*) imixcr, ddrmix
        read(istdin,*,end=90) imixcr, ddrmix
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) imixcr, ddrmix
      end if
c
c  parameters for convective overshoot
c
      if(ircovs.gt.0) then
        write(istdou,*) 'icnvos, jcnvos, clcovs, cldovs, alphos?'
        write(istdou,*)  icnvos, jcnvos, clcovs, cldovs, alphos
        read(istdin,*,end=90) icnvos, jcnvos, clcovs, cldovs, alphos
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*)  icnvos, jcnvos, clcovs, cldovs, alphos
        write(istdou,*) 'icsove, icsovc, alpove, alpovc?'
        write(istdou,*)  icsove, icsovc, alpove, alpovc
        read(istdin,*,end=90) icsove, icsovc, alpove, alpovc
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) icsove, icsovc, alpove, alpovc
      end if
c
c  parameters for rotation
c
      if(ircrot.gt.0) then
        write(istdou,*) 'isprot, velrot?'
        write(istdou,*)  isprot, velrot
        read(istdin,*,end=90) isprot, velrot
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*)  isprot, velrot
      end if
c
c  boundary condition definition, including atmosphere
c
      if(ircbcs.gt.0) then
        write(istdou,*) 
     *    'ntaubc, sbcfct, albdbc, tmnbc, tmxbc, flsatm, iqfit,icsrad?'
        write(istdou,*) 
     *    ntaubc, sbcfct, albdbc, tmnbc, tmxbc, flsatm, iqfit,icsrad
        read(istdin,*,end=90) 
     *    ntaubc, sbcfct, albdbc, tmnbc, tmxbc, flsatm, iqfit,icsrad
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) 
     *    ntaubc, sbcfct, albdbc, tmnbc, tmxbc, flsatm, iqfit,icsrad
        write(istdou,*) 'qqbc?'
        write(istdou,*) qqbc
        read(istdin,*,end=90) qqbc
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) qqbc
        write(istdou,*) 'ihe3bc?'
        write(istdou,*) ihe3bc
        read(istdin,*,end=90) ihe3bc
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) ihe3bc
      end if
c
c  integration procedure definition
c
      if(ircint.gt.0) then
        write(istdou,*) 'eps, nit, ucyt, nucy0, nucy, inentr, icncbc?'
        write(istdou,*) eps, nit, ucyt, nucy0, nucy, inentr, icncbc
        read(istdin,*,end=90) eps, nit, ucyt, nucy0, nucy, inentr, 
     *    icncbc
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) eps, nit, ucyt, nucy0, nucy, inentr, icncbc
        write(istdou,*) 'thetad?'
        write(istdou,*) thetad
        read(istdin,*,end=90) thetad
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) thetad
        write(istdou,*) 'itndgn, itnprt?'
        write(istdou,*) itndgn, itnprt
        read(istdin,*,end=90) itndgn, itnprt
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) itndgn, itnprt
      end if
c
c  mesh stretching definition
c
      if(ircmsh.gt.0) then
        write(istdou,*) 
     *    'istrtc, epsr, wx, wpt, wr, wxh, wx3, wdgr, iwdgrd?'
        write(istdou,*) 
     *    istrtc, epsr, wx, wpt, wr, wxh, wx3, wdgr, iwdgrd
        read(istdin,*,end=90) 
     *    istrtc, epsr, wx, wpt, wr, wxh, wx3, wdgr, iwdgrd
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) 
     *    istrtc, epsr, wx, wpt, wr, wxh, wx3, wdgr, iwdgrd
        write(istdou,*) 'wmshos,drmsos,wmshcc,dqlmsc,qmscmn,qmscmx?'
        write(istdou,*)  wmshos,drmsos,wmshcc,dqlmsc,qmscmn,qmscmx
        read(istdin,*,end=90) wmshos,drmsos,wmshcc,dqlmsc,qmscmn,qmscmx
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*)  wmshos,drmsos,wmshcc,dqlmsc,qmscmn,qmscmx
        write(istdou,*) 'nrpeat, iprmsh?'
        write(istdou,*) nrpeat, iprmsh
        read(istdin,*,end=90) nrpeat, iprmsh
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) nrpeat, iprmsh
      end if
c
c  time step definition
c
      if(irctst.gt.0) then
        write(istdou,*) 
     *    'dt0, dtmx, dtmn, dymx, xdtden, aldtrt, eta, dtfcrd, dtceps?'
        write(istdou,*) 
     *    dt0, dtmx, dtmn, dymx, xdtden, aldtrt, eta, dtfcrd, dtceps
        read(istdin,*,end=90) 
     *    dt0, dtmx, dtmn, dymx, xdtden, aldtrt, eta, dtfcrd, dtceps
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) 
     *    dt0, dtmx, dtmn, dymx, xdtden, aldtrt, eta, dtfcrd, dtceps
      end if
c
c  output definition
c
      if(ircout.gt.0) then
        write(istdou,*) 
     *    'icmout, istosc, ioscpr, ndsmth, icasex, intr, isetos?'
        write(istdou,*) 
     *    icmout, istosc, ioscpr, ndsmth, icasex, intr, isetos
        read(istdin,*,end=90) 
     *    icmout, istosc, ioscpr, ndsmth, icasex, intr, isetos
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) 
     *    icmout, istosc, ioscpr, ndsmth, icasex, intr, isetos
        write(istdou,*) 'iastr, ifstr, ilstr, ispcpr?'
        write(istdou,*) iastr, ifstr, ilstr, ispcpr
        read(istdin,*,end=90) iastr, ifstr, ilstr, ispcpr
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) iastr, ifstr, ilstr, ispcpr
      end if
c
c  diagnostics definition
c
      if(ircdgn.gt.0) then
        write(istdou,*) 
     *    'idgbcs,idgrhs,idgeos,idgopc,idgeng, idgtnr, idgmxc, idgsdt?'
        write(istdou,*) 
     *    idgbcs,idgrhs,idgeos,idgopc,idgeng, idgtnr, idgmxc, idgsdt
        read(istdin,*) 
     *    idgbcs,idgrhs,idgeos,idgopc,idgeng, idgtnr, idgmxc, idgsdt
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) 
     *    idgbcs,idgrhs,idgeos,idgopc,idgeng, idgtnr, idgmxc, idgsdt
        write(istdou,*) 'idgn75, idgn76, idgn77?'
        write(istdou,*) idgn75, idgn76, idgn77
        read(istdin,*) idgn75, idgn76, idgn77
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) idgn75, idgn76, idgn77
        write(istdou,*) 'itbcrh, irhtst?'
        write(istdou,*) itbcrh, irhtst
        read(istdin,*,end=90) itbcrh, irhtst
        write(istdou,*) '     -----------  New values  --------'
        write(istdou,*) itbcrh, irhtst
      end if
c
c  succesfully completed input of parameters
c
      ibegin = 1
c
c
c      if(idsexc.gt.0) read(idsexc,exec,end=90)
c
      if(istart.le.0) stop 'mnevol'
c
c  set absolute value of itrial (note that sign is used to flag for
c  single precision input)
c
      itriab=iabs(itrial)
c
c  set dependencies and make consistency checks
c
c  for isetos = +-1 or icntsl = 1, force resetting of parameters from
c  trial model
c  also set itrial to 2 to avoid interpolation, and force using
c  CNO abundances from trial, if possible
c
      if(isetos.eq.1.or.isetos.eq.-1.or.icntsl.eq.1) then
	ipartr=1
	itriab=2
	itrial=2
	istcno=1
	isthec=1
      end if
c
c  for ipartr = 1, test for reading data
c
      if(ipartr.eq.1.and.(itriab.eq.1.or.idatmd.eq.0)) then
        write(istdou,190) 
        if(istdpr.gt.0.and.istdou.ne.istdpr) write(istdpr,190)
        go to 20
      end if
c
c  for ntaubc .le. 1 (i.e. no atmosphere match) force icsrad = 0
c
      if(ntaubc.le.1) icsrad = 0
c
c  for iopacm = 0 (no atmospheric opacity) force iopatm = 0
c
      if(iopacm.eq.0) iopatm = 0
c
c  test for equation of state version, depending on ivreos set
c  by setcns
c
      if(ivreos.ne.1) then
        if(modeeq.ne.0) then
          write(istdou,1001) ivreos, modeeq
          if(istdpr.gt.0.and.istdou.ne.istdpr) write(istdpr,1001) 
     *      ivreos, modeeq
          modeeq=0
        end if
c
      else
c
c  otherwise set icoulm appropriately, to match modeeq
c
        if(modeeq.eq.0) then
          icoulm=-1
        else if(modeeq.eq.1) then
          icoulm=0
        else if(modeeq.eq.2) then
          icoulm=2
        else
          write(istdou,1002) modeeq
          if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,1002) modeeq
          stop 'Stop 1 in evolmain'
        end if
      end if
c
c$nl      write(6,exec)
c
c  test for resetting printed output
c
      if(istdpr.gt.0.and.istdpr.ne.istdou) 
     *  call openfc(istdpr,istdpp,'u','f')
c
c  call routine for special (nonstandard) input
c  So far additional input must be set in commons
c
      call spcinp
c
c  end of control input
c
c              *************************************
c
c
c  set data sets
c
      call openfc(idstrl,idstrp,'o','u')
      if(isetzh.gt.0) call openfc(idshvz,idshvp,'u','f')
      if(iastr.eq.1) call openfc(idsevl,idsevp,'u','u')
      if(ifstr.eq.1.or.ilstr.eq.1) call openfc(idsefl,idsflp,'u','u')
c
c  test for oscillation output
c
      if(isetos.eq.1.or.istosc.eq.1.or.ioscpr.eq.1) then
        ist = 1
      else
        ist = 0
      end if
c
      if(ist.eq.1) then
        call openfc(idsoov,idsovp,'u','u')
        call openfc(idsnov,idsnvp,'u','u')
        call openfc(idsgng,idsgnp,'u','u')
        call openfc(idsgsm,idsgsp,'u','u')
      end if
c
      if(ist.eq.1.or.icmout.eq.1) then
        call openfc(idscen, idscnp, 'u','f')
        if(iducen.gt.0) call openfc(iducen, iducnp, 'u','u')
        call openfc(idssum, idssmp, 'u','f')
      end if
c
      if(iterlr.gt.0) then
        call openfc(idstm1,idst1p,'u','u')
        call openfc(idstm2,idst2p,'u','u')
      end if
c
      if(isethv.eq.0.or.isethv.eq.4) call openfc(idszab,idszap,'o','f')
      call openfc(inopc,inopcp,'o','u')
      if(ifdgop.eq.4) then
	call openfc(idopfd,idopfp,'o','f')
      end if
c
c  open log file to unit 99 (use default if unit has not been
c  assigned).
c
      idslog=99
      call stfile(idslog,ndslog)
      if(ndslog.eq.-1) then
        write(istdou,102)
        open(idslog,file='evol-file.log',status='unknown')
      else
        call openf(idslog,'u','f')
      end if
c
c  dump list of opened files to log file
c
      call dmpfil(idslog)
      close(idslog)
c
c  unpack control for stretching
c
      istrt0=mod(istrtc,10)
      istrt1=mod(istrtc/10,10)
      istrt2=mod(istrtc/100,10)
c
c  set maximum model no, for use if isetos = 1
c
      nmdmax=nmdtrl+nt-1
c
c  as fudge force tlmopc to 0 when iopatm = 1, thus forcing interior
c  opacity everywhere except in s/r atmos
c
      if(iopatm.eq.1) tlmopc=0
c
c  print information about equation of state
c
      if(istdpr.gt.0) then
        if(modeeq.eq.0) then
          write(istdpr,106)
        else if(modeeq.eq.1) then
          write(istdpr,107)
        else
          write(istdpr,108)
        end if
      end if
c
c  initialize previous time step
c
      dtp=dt
c
c  logics for  internal iteration for rs and ls
c
c  reset iterlr if .ne. 0 - 5
c
      if(iterlr.ne.0) iterlr=min0(5,max0(1,iterlr))
c
c  counting index for iteration
c
      itlr=1
c
c  counting indices for model output during iteration
c
      ipmlr=0
      ipmlr0=0
c
      go to 22000
c
c  entry point for continuing iteration
c
21001 itlr=itlr+1
      ipmlr0=ipmlr
c
c  as trial model set initial model from previous iteration
c
      if(nt.ne.0) then
        rewind idstm1
        read(idstm1) iform,nnt,nrdtmr,nidtmr,nvarrd,nbccf,
     *    (datmod(i),i=1,nrdtmr),(ndtmod(i),i=1,nidtmr),
     *    (x(n),(y(i,n),i=1,nvarrd),n=1,nnt),
     *    (bccoef(i),i=1,nbccf)
      else
        nnt=nn
      end if
c
c  reset xh
c
      fcxh=xxh/y(5,1)
      do 21005 n=1,nn
21005 y(5,n)=fcxh*y(5,n)
      aztbc(3)=fcxh*aztbc(3)
      axbc(3)=fcxh*axbc(3)
      go to 40
c
c
22000 agefns=syear*agefin
c
   33 if(ifsrun.eq.0) go to 40
c
c  relative composition of heavy elements
c  **************************************
c
      if(isethv.eq.0.or.isethv.eq.4) then
c
c  read relative abundances of heavy elements
c
	call skpcom(idszab)
        read(idszab,*) iab
	if(isethv.eq.0) then
c
c  read mass fractions
c
          do 34 i=1,iab
   34     read(idszab,*) zab(i)
	else
c
c  read number fractions, and reset with garbage element
c  (so far hard-coded to be Fe)
c
          do 35 i=1,iab
   35     read(idszab,*) xnum(i)
	  ngarb = 10
c
c  as a patch, use xxh in setting heavy-element abundance
c  *** Needs fixing
c
	  call sthvab(xxh, z, xnum, ah, amz, iab, ngarb, zab)
	end if
      else if(isethv.eq.1) then
c  set ross & aller c, o and fe abundances
        iab=10
        call zero(zab,10)
        zab(1)=0.2254
        zab(3)=0.4987
        zab(10)=0.0795
c
      end if
c
c  reset heavy element parameters. note that the whole treatment
c  of heavy element abundances needs to be cleared up
c
      if(isethv.eq.1) then
c  reset az
        az=16.389
      else if(isethv.ne.3) then
c
c  reset az,anz0
c
        az=0
        anz0=0
        do 36 i=1,iab
        az=az+zab(i)/amz(i)
   36   anz0=anz0+zab(i)*izz(i)/amz(i)
        az=1/az
        if(istdpr.gt.0) write(istdpr,136) az,anz0
c
      end if
c
c  test for resetting anhe0
c
   37 if(ianhe0.eq.1.or.isethv.eq.3) then
c
c  make dummy call of eqstf
c
        nosd=.true.
        notd=.true.
        xh=0.7
        yh=1-xh-z
        fl=-8.
        tl=4.5
        iwrthp=iwrthv
        iwrthv=1
        call eqstf(fl,tl,xh,yh,z,nosd,notd)
        iwrthv=iwrthp
        if(istdpr.gt.0.and.isethv.ne.3) write(istdpr,109) anhe0
        anhe0=6
      end if
c
c  test for setting GONG values of heavy element quantities
c
      if(isethv.eq.3) then
        az=16
        anz0=0.5
        anhe0=8
c
        if(istdpr.gt.0) write(istdpr,110) az, anz0, anhe0
c
      end if
c
c  storage indices
c
   40 in1=1
      in11=0
      in=ivarmx+1
c
c  test for setting of oscillation quantities from model on file
c
      if(isetos.eq.1.or.isetos.eq.-1) go to 44
c
c  parameters for tnrk to give initial solution
c
      ii1=4
      ii2=0
      ii3=0
      ii=ii1+ii2+ii3
      kk=0
      ka=2
      kb=2
      ki=0
      kt=0
c  set v and theta
   41 do 42 k=1,ivarmx
      v(k)=k
   42 theta(k)=thetad(k)
c
      if(icntsl.ge.1) go to 44
c  time = 0
      time0=.true.
c
c  initial age
c
      age=age0
c
c  initialize flag for final time steps
c
      ifindt=0
c
c  if iterating for rs and ls skip setting up of trial after first
c  iteration
c
      if(itlr.gt.1.or.itrial.eq.0) go to  48
c
   44 continue
c
c  trial solution
c  **************
c
      call settrl(idstrl,itrial,nmdtrl,nrdtrl,nvartr,idatmd,nspcmx,
     *  nbcprv,x,y,datmod,ndtmod,bccoef,iy,iformr,nnt,nrdtmr,nidtmr,
     *  nbccfr,nspcmr,icnotr,ihectr,icry)
c
      if(icry.lt.0) go to 20
c
c  test for continuing solution or setting oscillation quantities
c
      if(icntsl.ne.1.and.isetos.ne.1.and.isetos.ne.-1) then
c
c  reset xh?
c
        if(ix.eq.2) then
c  set xh uniformly to xxh
          do 235 n=1,nnt
  235     y(5,n)=xxh
          aztbc(3)=xxh
          axbc(3)=0
        else if(ix.eq.3) then
c  scale xh, to xxh on surface
          fcxh=xxh/y(5,1)
          do 238 n=1,nnt
  238     y(5,n)=fcxh*y(5,n)
          aztbc(3)=fcxh*aztbc(3)
          axbc(3)=fcxh*axbc(3)
c
        end if
      end if
c
c  test for resetting of parameters to values from file
c
      if(ipartr.eq.1) then
c
c  store datmod in datout, in case this is a run with isetos = 1
c
        call store(datmod,datout,nrdtmr)
        call istore(ndtmod,ndtout,nidtmr)
c
        icase=ndtmod(1)
c
c  note that icase was added at the time where qqbc was extended to
c  5 variables. thus we can test on icase to see how long qqbc is.
c  note, however, that icase may have been augmented by 1.d7 in
c  call of s/r sclcbc, to flag for rescaling.
c  hence base test on mod(icase,1e7).
c
        z=datmod(1)
        alfa=datmod(3)
        etac=datmod(4)
        phc=datmod(5)
c
	iche3=ndtmod(4)
        if(ix3.eq.2) agehe3=datmod(6)
        ntaubc=ndtmod(3)
        tmnbc=datmod(12)
        tmxbc=datmod(13)
        sbcfct=datmod(14)
	fcno=datmod(52)
	if(nrdtmr.ge.140.and.datmod(140).lt.0) flsatm=datmod(140)
c
        iqqbct=5
        iversr=mod(icase/10 000 000,10)
        if(iversr.eq.0) iqqbct=4
c
        if(ntaubc.le.1) then
          if(iversr.le.2) then
            albdbc=1
          else
            albdbc=datmod(15)
          end if
        else
          call store(datmod(15),qqbc,iqqbct)
        end if
c
        if(iversr.le.4) then
          icsrad=0
        else
          icsrad=mod(icase/1 000,10)
        end if
c
c
        age=datmod(22)*syear
        am=datmod(23)/amsun
c
        ifdgop=mod(icase/1 000 000,10)
        fdgopl=datmod(26)
c
c  reset ifdgop to 1 if fdgopl .ne. 0 and ifdgop = 0
c
        if(fdgopl.ne.0.and.ifdgop.eq.0) ifdgop=1
c
        if(ifdgop.gt.1) then
          if(datmod(27).gt.0) tlopfg=datmod(27)
          if(datmod(28).gt.0) dtlopf=datmod(28)
        end if
c
c  test for setting of patmos
c
c  note that as a result of resetting  of datmod (23/11/82), datmod(21)
c  is non-zero and patmos is in datmod(20) instead of in datmod(19)
c
        iptmos=19
        if(datmod(21).ne.0) iptmos=20
        if(datmod(iptmos).gt.0) patmos=datmod(iptmos)
c
	if(ndtmod(12).gt.100) iwdopc=ndtmod(12)/100
	ivreng=mod(icase/100 000 000,10)
c
c  test for setting rotation. Note that this only works for isprot = 1
c  (uniform rotation at all times). Otherwise, storage of varying rotation
c  in y must be organized, somehow.
c
	if(nidtmr.ge.35) then
	  isprot=ndtmod(35)
	  if(isprot.gt.0) velrot=datmod(136)
        end if
c
        if(nidtmr.ge.45) then
	  iconcs=ndtmod(31)
	  imxlng=ndtmod(32)
        end if
c
        if(nidtmr.ge.52) then
	  icsove=ndtmod(51)
	  icsovc=ndtmod(52)
	  if(icsove.ge.1) alpove=datmod(151)
	  if(icsovc.ge.1) alpovc=datmod(152)
        end if
c
	if(nidtmr.ge.61) iheccs=ndtmod(61)
c
        if(nrdtmr.ge.155.and.datmod(155).gt.0) alshft=datmod(155)
c
      end if
c
c                 -------------------------------
c
c  end setting trial solution
c
c  possibly output final set of control parameters
c
      if(istdpr.gt.0) then
	write(istdpr,
     *    '(// '' Input parameters, possibly reset from trial''/
     *         '' *******************************************'')')
        write(istdpr,*) 'cntrd'
        write(istdpr,*) cntrd,tail
        write(istdpr,*) 'idsexc, istart'
        write(istdpr,*) idsexc, istart,tail
        write(istdpr,*) 
     *    'istdpr, idstrl, idshvz, idszab, idsevl, idsefl, idssum'
        write(istdpr,*) 
     *    istdpr, idstrl, idshvz, idszab, idsevl, idsefl, idssum,tail
        write(istdpr,*) 'idscen, idsnov, idsoov, idsgng, idsgsm'
        write(istdpr,*) idscen, idsnov, idsoov, idsgng, idsgsm, tail
        write(istdou,*) 'idstm1, idstm2, iducen'
        write(istdou,*) idstm1, idstm2, iducen, tail
        write(istdpr,*) 'am, z, nn'
        write(istdpr,*) am, z, nn,tail
        write(istdpr,*) 'age0, nt, agefin, icntsl'
        write(istdpr,*) age0, nt, agefin, icntsl,tail
        write(istdpr,*) 'iterlr, nitlr, rsfin, alsfin, epslr'
        write(istdpr,*) iterlr, nitlr, rsfin, alsfin, epslr,tail
        write(istdpr,*) 
     *    'itrial, nmdtrl, ix, xxh, ix3, istcno, isetzh, isthec'
        write(istdpr,*) 
     *    itrial, nmdtrl, ix, xxh, ix3, istcno, isetzh, isthec, tail
        write(istdpr,*) 'nvartr, idatmd, itrdum, ipartr, iptr'
        write(istdpr,*) nvartr, idatmd, itrdum, ipartr, iptr,tail
        write(istdpr,*) 'ihvz, isethv, ianhe0, modeeq, iomfll'
        write(istdpr,*) ihvz, isethv, ianhe0, modeeq, iomfll,tail
        write(istdpr,*) 'iwdopc, xhsopc, tsmn, tstr, rhsmn, rhsmx'
        write(istdpr,*) iwdopc, xhsopc, tsmn, tstr, rhsmn, rhsmx,tail
        write(istdpr,*) 'timx, rhimn, rhimx, sigstr, inopc, iopccf'
        write(istdpr,*) timx, rhimn, rhimx, sigstr, inopc, iopccf,tail
        write(istdpr,*) 'ifdgop, fdgopl, tlopfg, dtlopf, tlopf1, idopfd'
        write(istdpr,*) 
     *    ifdgop, fdgopl, tlopfg, dtlopf, tlopf1, idopfd,tail
        write(istdpr,*) 'iopacm, iopatm, patmos, alamop, zatmop'
        write(istdpr,*) iopacm, iopatm, patmos, alamop, zatmop,tail
        write(istdpr,*) 
     *    'fcno, iscren, thte, iche3, agehe3, icnocs, ivreng, iheccs'
        write(istdpr,*) 
     *    fcno, iscren, thte, iche3, agehe3, icnocs, ivreng, iheccs,tail
        write(istdou,*) 'xzer3, xrz12, xrz13, xrz14, xrz16'
        write(istdou,*) xzer3, xrz12, xrz13, xrz14, xrz16,tail
        write(istdou,*) 'ifdeps, epsfdg, qepsf1, qepsf2'
        write(istdou,*) ifdeps, epsfdg, qepsf1, qepsf2
        write(istdou,*) 'iconcs, imxlng, iturpr, tprfct'
        write(istdou,*)  iconcs, imxlng, iturpr, tprfct, tail
        write(istdpr,*) 'alfa, etac, phc'
        write(istdpr,*) alfa, etac, phc,tail
        write(istdpr,*) 'imixcr, ddrmix'
        write(istdpr,*) imixcr, ddrmix,tail
        write(istdpr,*) 'icnvos, jcnvos, clcovs, cldovs, alphos'
        write(istdpr,*)  icnvos, jcnvos, clcovs, cldovs, alphos, tail
        write(istdpr,*) 'icsove, icsovc, alpove, alpovc'
        write(istdpr,*)  icsove, icsovc, alpove, alpovc, tail
        write(istdpr,*) 'isprot, velrot'
        write(istdpr,*)  isprot, velrot, tail
        write(istdpr,*) 
     *    'ntaubc, sbcfct, albdbc, tmnbc, tmxbc, flsatm, iqfit,icsrad'
        write(istdpr,*) 
     *    ntaubc,sbcfct,albdbc,tmnbc,tmxbc,flsatm,iqfit,icsrad,tail
        write(istdpr,*) 'qqbc'
        write(istdpr,*) qqbc,tail
        write(istdou,*) 'ihe3bc'
        write(istdou,*) ihe3bc
        write(istdpr,*) 'eps, nit, ucyt, nucy0, nucy, inentr, icncbc'
        write(istdpr,*) eps, nit, ucyt, nucy0, nucy, inentr, icncbc,tail
        write(istdpr,*) 'thetad'
        write(istdpr,*) thetad,tail
        write(istdpr,*) 'itndgn, itnprt'
        write(istdpr,*) itndgn, itnprt,tail
        write(istdpr,*) 
     *    'istrtc, epsr, wx, wpt, wr, wxh, wx3, wdgr, iwdgrd'
        write(istdpr,*) 
     *    istrtc, epsr, wx, wpt, wr, wxh, wx3, wdgr, iwdgrd, tail
        write(istdou,*) 'wmshos,drmsos,wmshcc,dqlmsc,qmscmn,qmscmx'
        write(istdou,*)  wmshos,drmsos,wmshcc,dqlmsc,qmscmn,qmscmx
        write(istdpr,*) 'nrpeat, iprmsh'
        write(istdpr,*) nrpeat, iprmsh,tail
        write(istdpr,*) 
     *    'dt0, dtmx, dtmn, dymx, xdtden, aldtrt, eta, dtfcrd, dtceps'
        write(istdpr,*) 
     *    dt0,dtmx,dtmn,dymx,xdtden,aldtrt,eta,dtfcrd,dtceps,tail
        write(istdpr,*) 
     *    'icmout, istosc, ioscpr, ndsmth, icasex, intr, isetos'
        write(istdpr,*) 
     *    icmout, istosc, ioscpr, ndsmth, icasex, intr, isetos,tail
        write(istdpr,*) 'iastr, ifstr, ilstr, ispcpr'
        write(istdpr,*) iastr, ifstr, ilstr, ispcpr,tail
        write(istdpr,*) 
     *    'idgbcs,idgrhs,idgeos,idgopc,idgeng, idgtnr, idgmxc, idgsdt'
        write(istdpr,*) 
     *    idgbcs,idgrhs,idgeos,idgopc,idgeng,idgtnr,idgmxc,idgsdt,tail
        write(istdpr,*) 'idgn75, idgn76, idgn77'
        write(istdpr,*) idgn75, idgn76, idgn77, tail
        write(istdpr,*) 'itbcrh, irhtst'
        write(istdpr,*) itbcrh, irhtst,tail
      end if
c
c  decode flag for selecting opacity
c
      iwdop0=mod(iwdopc,10)
      iwdop1=iwdopc/10
c
c unpack flag for iche3
c
      iche30=mod(iche3,10)
      if(iche30.eq.2) then
	iche31=0
	iche3 =iche30
      else
        iche31=iche3/10
      end if
c
c unpack flag for core mixing
c
      imixc0=mod(imixcr,10)
      imixc1=mod(imixcr/10,10)
      imixc2=mod(imixcr/100,10)
      imixc3=mod(imixcr/1000,10)
c
c  for consistency with diffusive version of programme,
c  introduce flag (set to zero) for diffusion
c
      idiffus=0
      idiffc1=0
c
c  test for case of fudge of energy generation
c
      if(ifdeps.gt.0) then
	ifdep1=ifdeps
      else if(ifdeps.lt.0) then
	ifdep1=-ifdeps
      else
	ifdep1=0
      end if
c
c  set number of terms in T(tau) fit, from RT parameter iqfit
c
      iqqbc=mod(iabs(iqfit),100)
c
c                   ****************************
c
c  set constants for energy generation
c
      call srncns
c
c  initialize index that may be reset in rhs for
c  zeroing equations for hydrogen abundance, nuclear energy
c  generation rate.
c
      nxhzer=nn+1
c
   48 continue
c
c  set Z into common/heavy/, depending on isetzh, in initial
c  model (with possible resetting later, after stretch)
c
      call sethvz(x,y,nnt,iy,datmod,ndtmod,z,isetzh,idshvz)
c
c  set xhsopc
c
      xhsopc=y(5,1)
c
c  initialize opacity table, except for subsequent  iterations for
c  rs and ls when using w.d. opacity
c
      if(initop.eq.0.or.
     *  (iwdop0.ne.1.and.iwdop0.ne.8.and.iwdop0.ne.9)) then
        rewind inopc
        call opcscf(iwdopc,iopacm)
	initop=1
      end if
c
c  initial time step
c
      delage=agefns-age
      if(agefin.gt.0.and.dt0.gt.delage) then
        dt=delage
        ifindt=1
      else if(agefin.gt.0.and.2*dt0.gt.delage) then
        dt=delage/2
        ifindt=1
      else
        dt=dt0
      end if
c
c  constants for rhs and bcs
c  *************************
c
c  note that am is mass in units of 1.989d33 g.
c  the constant of gravity = 6.6732e-8 cgs
c  the radiation density constant is 7.5647e-15 cgs.
c
c  changed 18/1/84, with corresponding changes in
c  s/r rhs, to express equations in terms of log r - 11.
c  this is required to avoid floating-point overflow on univac.
c
c  changed 4/1/85, to use log(r/1.d11) and log(l/1.d33) as
c  dependent variables.
c
c  numerical constants reset consistently 15/3/88
c  (but incorrectly). Factor am was missing in a3.
c  This was corrected 12/8/90.
c
c++      a1=1.5827959e-1*am
c++      a2=2.1008442d14*am**2
c++      a3=4.1654684d23*am
c++      a4=1.989*am
c
      amsunr=amsun/1.d33
      a1=amsunr*am/(4*pi)
      a2=1.d22*cgrav*amsunr*amsunr*am*am/(4*pi)
      a3=3*amsun*am/(64.0d11*pi*pi*clight*car)
      a4=amsunr*am
c
c  constant for including centrifugal acceleration in hydrostatic
c  support
c
      a5=0.6666666666667d22*amsunr*am/(4.d0*pi)
c
c..      write(6,*) 'consistent and actual a1', a1c, a1
c..      write(6,*) 'consistent and actual a2', a2c, a2
c..      write(6,*) 'consistent and actual a3', a3c, a3
c..      write(6,*) 'consistent and actual a4', a4c, a4
c
c  constants for bcs
c  note: b2 is never used, and hence not set.
c
c  b1 changed 4/1/85 to log(m/1.d33)
c  b3 and b4 changed to account for change of dependent variables.
c  also small changes of numerical values, to make b3 and b4
c  consistent with a1 - a4 (previous values were b1 = 33.29864+aml),
c  b3 = -3.14727 and b4 = 26.1226+aml+log10(sbcfct))
c
c  Note: b4 fixed on 13/9/90, by addition of aml = log(M/Msun)
c
      aml=log10(am)
c++      b1=0.29863478+aml
c++      b3=-14.1472377+log10(albdbc)
c++      b4=4.1229689+aml+log10(sbcfct)
c
      b1=log10(amsunr)+aml
      b3=log10(pi*clight*car)-11+log10(albdbc)
      b4=log10(cgrav*amsunr)+11+log10(sbcfct)+aml
c
c..      write(6,*) 'consistent and actual b1', b1c, b1
c..      write(6,*) 'consistent and actual b3', b3c, b3
c..      write(6,*) 'consistent and actual b4', b4c, b4
c
c  central and surface xh, for bcs at time=0
c
      xhc=y(5+in11,nnt)
      xhs=y(5+in11,1)
c
c  constants for mixlng
c  note that c1 is apparently not used
c
c++      c1=am*1.32730d26
c++      c2=3.02379e-4
      c1=am*cgrav*amsun
      c2=4*car*clight/3
c..      write(6,*) 'consistent and actual c1', c1c, c1
c..      write(6,*) 'consistent and actual c2', c2c, c2
c
c  set storage indices for abundances, including CNO elements.
c  Initialize abundance variable array cvr in common/compvr/
c
      call incomp(y,nnt,iy,aztbc,z,iche30,icnocs,iheccs,ihe3bc,
     *  istcno,isthec,icnotr,ihectr, icntsl,itlr,nt)
c
c  treatment of he3
c
      if(iche30.eq.1) then
c
c  he3 not in equilibrium
c
        ieqhe3=0
      else
c
c  he3 in equilibrium
c
        ieqhe3=1
      end if
c
c  initial he3 fudge, when he3 is not assumed to be in equilibrium
c
      if(iche30.ne.2.and.agehe3.gt.0) then
        ifdhe3=1
c
c  test for resetting fudge he3 age to actual model age
c  (only if age is .gt. 0)
c
        if(ix3.eq.3.and.age.gt.0) agehe3=age/syear
        agesh=agehe3*syear
        call zero(xhe3,4)
c
      else
        ifdhe3=0
        agesh=-1.
      end if
c
c  total number of variables for time-dependent solution
c
      nvar = 4+ispxx3+ispcno+isphec
c
      if(istdpr.gt.0) then
        write(istdpr,141) am,alfa,ix,ix3,y(5,1),z,agehe3
        if(ntaubc.le.1) then
          write(istdpr,142) sbcfct, albdbc
        else
          write(istdpr,143) ntaubc,tmnbc,tmxbc,sbcfct,qqbc
        end if
c
        write(istdpr,144) (k,theta(k),k=1,nvar)
c
      end if
c
c  initialize parameters for convective overshoot
c  Note that this may need to be moved (or repeated) after setting
c  of mesh in trial model.
c 
      incovs = 1
      call setovs(x,y,nnt,iy,incovs)
      incovs = 0
c
c  zero dtopfg(1-3)
c
      call zero(dtopfg,4)
c
c  test for diagnostics for opacity fudge
c
      if(ifdgop.gt.0) then
        fdgopc=10.**fdgopl
	if(istdpr.gt.0) then
          if(ifdgop.eq.1) write(istdpr,145) fdgopl,fdgopc
          if(ifdgop.eq.2) write(istdpr,146) fdgopl,tlopfg,dtlopf
          if(ifdgop.eq.3) write(istdpr,147) fdgopl,dtlopf,tlopf1,tlopfg
	end if
        if(ifdgop.eq.4) then
	  call stfile(idopfd,nopfd)
	  if(istdpr.gt.0) write(istdpr,148) file(nopfd)(1:75)
	else
c
          dtopfg(1)=fdgopl
c
c  note: this storage was implemented for ifdgop = 3 at 09:10 
c  on may 16 1983.
c  models computed and stored before this time with ifdgop = 3 
c  do not have tlopfg and dtlopf set into datmod.
c
          if(ifdgop.gt.1) then
            dtopfg(2)=tlopfg
            dtopfg(3)=dtlopf
            dtopfg(4)=tlopf1
          end if
        end if
c
      end if
c
c  test for setting of alamop to give fixed surface pressure
c  in model of present sun
c
      if(patmos.gt.0) then
c  assumed values for rs and ls for present sun
        arsatm=6.9599d10
        alsatm=3.8481d33
        amsatm=amsun*am
        if(iterlr.ne.0) then
          arsatm=rsfin
          alsatm=alsfin
        end if
c
        call patmit(patmos,alsatm,arsatm,amsatm,xhs,10,1.d-6,0,icry)
        if(istdpr.gt.0) write(istdpr,104) patmos,alamop
c
      end if
c
c  skip test and strecth for itlr.gt. 1
c
      if(itlr.gt.1) go to 49
c
c  test of bcs and rhs
c
      if(itbcrh.ge.2) then
        call tstbcs(bcs,x(1),x(nnt),y(1,1),y(1,nnt),zk,ap,aq,ii,5,
     *    kk,ka,kb,nnt,1.d-4)
      end if
      if(itbcrh.ge.1) then
c
        call tstrhs(rhs,x,y,zk,ap,aq,ii,5,kk,nnt,iy,1.d-4,20)
c
      end if
c
c  change number of mesh points
c
      if(itriab.eq.3) then
        icngrm=1
        resrc=.false.
c
c  set al0, al2, and zero am0, am2, for setting radius at innermost
c  meshpoint
c
        al0=bccoef(1)
        al2=bccoef(2)
        am0=0
        am2=0
c
        intork=2
        intvnt=2
c
	nvarev=nvar
	nvar=nvartr
c
        call mshstr(x,y,in1,in,iy,nnt,nn,rhs,resrc)
c
	nvar=nvarev
c
      end if
c
      icngrm=0
c
c  change mesh in trial solution?
c
      if(itriab.eq.4) then
        resrc=.true.
c
c  set al0, al2, and zero am0, am2, for setting radius at innermost
c  meshpoint
c
        al0=bccoef(1)
        al2=bccoef(2)
        am0=0
        am2=0
c
        intork=1
        intvnt=1
c
	nvarev=nvar
	nvar=nvartr
c
        call mshstr(x,y,in1,in,iy,nnt,nn,rhs,resrc)
c
	nvar=nvarev
      end if
c
c  test for setting angular velocity 
c
      if(isprot.gt.0) then
c
c  initialize angular velocity and moment of inertia
c
	introt=-1
	call setrot(x,y(in1,1),y(in,1),iy,nn,velrot,introt)
      else
	omgrot(1)=0
      end if
c
c  test if this is only setting oscillation quantities 
c  or printing for given model
c
      if(isetos.eq.1.or.isetos.eq.-1) then
        nn=nnt
        ifdhe3=0
        time0=age.eq.0
        ifprt=-2
        xhc=y(5,nn)
        xhs=y(5,1)
        xhe3(1)=aztbc(4)
        agey=age/syear
        comout=.true.
        if(isetos.eq.1) then
          ist=1
        else
          ist=-1
        end if
      end if
c
c  set quantities at convection-zone boundaries and possibly
c  prepare for convective overshoot
c
      icscbr=0
      call setcbr(x,y(in1,1),rhs,nn,iy,icscbr,iqcres)
c
      in11=in1-1
c
c  skip to printing when printing or setting oscillation quantities
c
      if(isetos.eq.1.or.isetos.eq.-1) then
c
c  set inc to -1 to flag for setting convection zone boundaries 
c  in s/r rhs
c     
        inc=-1
c
        if(istdpr.gt.0) write(istdpr,132) nrdtrl,age,agey
c
        go to 705
      end if
c
      if(iptr.gt.0) then
c
c  print trial solution
c
        if(istdpr.gt.0) write(istdpr,100)
c  ensure that y(1,.) is monotonic
        nn1=nn-1
        do 26 n=2,nn1
        if(y(1,n).le.y(1,n+1)) y(1,n)=(y(1,n-1)+y(1,n+1))/2
   26   continue
        call prtsol(x,y,y,zk,ap,aq,rhs,nn,iy,iptr,-1,-1,icasex,iformr,
     *    agey,datmod,ndtmod,ndsmth)
c
      end if
c
c  copy solution to y(in,.) for case of convergence failure
c
      do n=1,nn
        do i=1,ivarmx
	  y(in+i-1,n)=y(i,n)
	end do
      end do
c
c              -----------------------------------
c
c  end of section for trial solution and setting up
c
c  initial parameters for stretch
c
   49 resrc=.true.
      intork=1
      intvnt=1
      ntime=0
c
c  reset Z into common/heavy/, depending on isetzh.
c
      call sethvz(x,y,nn,iy,datmod,ndtmod,z,isetzh,idshvz)
c
c  set initial Z-value
c
      if(icntsl.eq.1) then
        zhinit=datmod(1)
      else
        zhinit=zh(1)
      end if
c
c  initialize nmxcor
c
      nmxcor=0
c
c  test whether this is the continuation of a previously computed
c  sequence
c
      if(icntsl.ge.1) then
c
c  possibly set boundary of convective core and mixed region
c
        inc=0
        jc=1
        ireset=0
        icomp0=1
        if(istdpr.gt.0) write(istdpr,127)

        call conmxc(x,y(1,1),iy,jc,nn,icomp0,ireset)
c
c  test for setting mixed-core boundary
c
        if(inc.gt.0.and.nl(inc).eq.nn) then
          if(ddrmix.eq.0) then
	    iextrp=3
          else
	    iextrp=10
          end if
	  iter=1
          call mixcor(x,y(1,1),iy,nn,compc,iextrp,0)
	  iter=0
	else
	  nmxcor=0
	  qmxcor=0
        end if
        go to 76
      end if
c
c  new model sequence
c
      ntime=-1
      age=age-dt
c
c  counting index for first stretching
c
c  for first stretching in evolution sequence repeat nrpeat times
c  unless this is later iteration for rs and ls
c
      if(istrtc.ge.1.and.itlr.le.1) then
        nmsh=0
      else
        nmsh=nrpeat
      end if
c
      nucy1=nucy0+1
      irnucy=0
c
c  stepping in time
c  ----------------
c
   50 ntime=ntime+1
c
c  possibly reset parameters for 4He output
c
      if(ntime.eq.-800) then
	iprhe4=0
	nphe41=280
	nphe42=310
	ndphe4=1
      end if
c
c  possibly reset parameters for diagnostics
c
      if(ntime.eq.-50) then
	idgn77=2
      end if
      if(istdpr.gt.0.and.idgbcs.ge.1) write(istdpr,445) bccoef
  445 format(///' bccoef:'/(1p10e13.5))
c
c  initialization for test of integrated hydrogen abundance
c
      xtimep=age
      xhintp=xhint
c
      if(itsxin.ge.1) call tstxin(x,y(in1,1),nn,iy,1,xhint,xhintp,
     *  age,xtimep,icry,'pos.0')
c
      in11=in1-1
      age=age+dt
      agey=age/syear
      dty=dt/syear
      write(istdou,105) ntime,age,agey,dt,dty
      if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,105) 
     *  ntime,age,agey,dt,dty
c
c  store parameters for convective borders at previous model
c
      dcsp=dcs
      dccp=dcp
      call istore(nf,nfp,6)
      call istore(nl,nlp,6)
      incp=inc
      call store(frcf,frcfp,6)
      call store(frcl,frclp,6)
      call store(xrcf,xrcfp,6)
      call store(xrcl,xrclp,6)
c
c  set flag for last model in sequence
c
      lastmd=ntime.eq.nt.or.age.ge.agefns*(1.d0-1.d-7)
c
c  the iteration
c  *************
c
   53 ucy=ucyt
c
      xtst1=0.
      xtst2=0.8
      if(ntime.ge.75) then
	xtst1=0.07
	xtst2=0.03
      end if
      ndxtst=2
      if(xtst1.gt.xtst2) then
        do n=1,nn
	  qq=10.**x(n)
	  if(qq.ge.xtst1) nxtst1=n
	  if(qq.ge.xtst2) nxtst2=n
        end do
      else
	nxtst1=0
	nxtst2=0
      end if
c
      if(nxtst1.lt.nxtst2.and.istdpr.gt.0) then
        write(istdpr,57099) 'y(in,.)',0,(n, 10.**x(n),y(in+4,n),
     *    y(in+iyche4-1,n),n=nxtst1,nxtst2,ndxtst)
        write(istdpr,57099) 'y(in1,.)',0,(n, 10.**x(n),y(in1+4,n),
     *    y(in1+iyche4-1,n),n=nxtst1,nxtst2,ndxtst)
      end if
c
c  initialize flag for including entropy term (may be reset as a last
c  resort, if iteration fails to converge)
c
       inentc=inentr
c
c  test for setting composition in convective core
c  Note: change test to ignore very small radiative core
c  Also, use a crude fudge to eliminate completely convective models
c  (with "convective cores" extending to very low temperature)
c
      if(ntime.gt.0.and.inc.gt.0) then
	concor=nl(inc).ge.nn-2.and.y(3,nf(inc)).ge.5
      else
        concor=.false.
      end if
c
      if(concor) then
c
c  make predictor step for treatment of convective-core composition
c
        iprdcr=1
	inmixc=0
        it=-1
        if(istdpr.gt.0.and.iprhe4.ge.1) write(istdpr,
     *    '(a/(i5,1p5e15.7))') 
     *    '4He abundances before predictor call of cmpcvc:',
     *    (n,x(n),y(in11+iyche4,n)-dprhe4,y(in11+iycc12,n),
     *    cvr(icvhe4,n)-dprhe4,cvr(icvc12,n),n=nphe41,nphe42,ndphe4)
        call cmpcvc(x,y,in,in1,iy,nn,compc,icomp,iche30,icnocs,
     *    iheccs,dt,it,iprdcr,icrycm)
        if(istdpr.gt.0.and.iprhe4.ge.1) write(istdpr,
     *    '(a/(i5,1p5e15.7))') 
     *    '4He abundances after predictor call of cmpcvc:',
     *    (n,x(n),y(in11+iyche4,n)-dprhe4,y(in11+iycc12,n),
     *    cvr(icvhe4,n)-dprhe4,cvr(icvc12,n),n=nphe41,nphe42,ndphe4)
c
      else
c
c  set both versions of the convective-core mass to zero
c  (Note: usage of these masses badly needs clearing up)
c
	nmxcp=0
	inmixc=1
	qmxcp=0
        cqcp=0
      end if
c
c  test for including possible luminosity shift
c
      if(ishfal.gt.0) then
        initsh=1
        call shftal(x,y(in,1),y(in1,1),nn,iy,initsh)
        initsh=0
      end if
c
c  set counter for possibly repeating solution this time step
c
      nrpsol=0
c
   57 if(idgn76.ne.0) then
        open(76,file='ttt.emdl',status='unknown',form='unformatted')
	write(76) iform,nn,nrdtmd,nidtmd,nvar,nbcprv,
     *    (datout(i),i=1,nrdtmd),(ndtout(i),i=1,nidtmd),
     *    (x(n),(y(in11+i,n),i=1,nvar),n=1,nn),
     *    (bccoef(i),i=1,nbcprv)
      end if
c
      if(idgn77.ne.0) then
        open(77,file='ttt.rhs',status='unknown',form='unformatted')
        write(77) nn,nvar,dt
      end if
c
c  augment counter for repeating solution
c
      nrpsol=nrpsol+1
c
      irczms=0
      iterlm=4
      irscvg=0
c
c..      if(ntime.ge.78) idgbcs = 2
c
      do 60 it=1,nit
c
      if(nxtst1.lt.nxtst2.and.istdpr.gt.0) then
        write(istdpr,57099) 'y(in,.)',1,(n, 10.**x(n),y(in+4,n),
     *    y(in+iyche4-1,n),n=nxtst1,nxtst2,ndxtst)
        write(istdpr,57099) 'y(in1,.)',1,(n, 10.**x(n),y(in1+4,n),
     *    y(in1+iyche4-1,n),n=nxtst1,nxtst2,ndxtst)
      end if
57099 format(//'X,Y-test for ',a10,' location ',i3/'n, q, X, Y:'/
     *  (i5,1p3e15.7))
c
c  test for consistency of measures of convective-core masses
c
      if(abs(qmxcor-cqc).gt.1.e-7.and.istdpr.gt.0) 
     *  write(istdpr,'(/'' **** Warning: qmxcor = '',1pe15.8,
     * '' .ne. cqc ='',e15.8,'' in evolmain'')') qmxcor, cqc
      if(abs(qmxcp-cqcp).gt.1.e-7.and.istdpr.gt.0) 
     *  write(istdpr,'(/'' **** Warning: qmxcp = '',1pe15.8,
     * '' .ne. cqcp ='',e15.8,'' in evolmain'')') qmxcp, cqcp
      if(istdpr.gt.0) write(istdpr,'(a,i5,1pe13.5,i7,e13.5)') 
     *  ' Before entering tnrkt, nmxcor, qmxcor, nmxcp, qmxcp =',
     *  nmxcor, qmxcor, nmxcp, qmxcp
c
      iter=it
c
      if(it.ge.nucy1) ucy=1.d0
      if(istdpr.gt.0) then
        write(istdpr,*) 'Before call of tnrkt, iopatm =',iopatm
c
        if(iprhe4.ge.1) write(istdpr,'(a/(i5,1p5e15.7))') 
     *    '4He abundances before call of tnrkt:',
     *    (n,x(n),y(in11+iyche4,n)-dprhe4,y(in11+iycc12,n),
     *    cvr(icvhe4,n)-dprhe4,cvr(icvc12,n),n=nphe41,nphe42,ndphe4)
c
      end if
c..      if(ishfal.gt.0) call shftal(x,y(in1,1),nn,iy,1,initsh)
      call tnrkt(x,y(in1,1),zk,y(in,1),zk,ap,aq,rhs,bcs,ii1,ii2,ii3,
     .  kk,ka,kb,ki,nn,iy,ucy,ea,det,v,theta,dt,iter)
      icrczm=1
      if(irczms.eq.1.and.mod(iter-itczms,2).eq.0) then
        if(istdpr.gt.0) write(istdpr,'(/a/)') ' Before call of rczmsh'
        j=1
        xxf=frcf(j)*x(nf(j))+(1-frcf(j))*x(nf(j)-1)
        xxl=frcl(j)*x(nl(j))+(1-frcl(j))*x(nl(j)+1)
        if(istdpr.gt.0) write(istdpr,116)
     *    j,nf(j),10.**x(nf(j)),nl(j),10.**x(nl(j)),
     *    x(nl(j)), frcl(j),
     *    xxf, 10.d0**xxf, xxl, 10.d0**xxl
        inrczm=0
        call rczmsh(x,y(in1,1),y(in,1),nn,iy,3,icrczm,inrczm)
      end if
c..      if(ishfal.gt.0) call shftal(x,y(in1,1),nn,iy,2,initsh)
c
c  test for errors in tnrkt
c
      if(v(1).le.0) then
	eam=1.d0
	go to 60010
      end if
c
c  As part of the fudge of the He3 abundance in convective cores,
c  reset He3 abundance in that part of the model which was convective
c  in the previous model
c
      if(cqcp.gt.0.and.iche30.eq.1) then
	do 57091 n=1,nn
	if(10.d0**x(n).lt.cqcp) y(in1+5,n)=cvr(3,n)
57091   continue
      end if
c
      if(istdpr.gt.0.and.itnprt.gt.0) write(istdpr,351) 
     *  (n,x(n),(y(in11+j,n),j=1,6), n=1,nn,itnprt)
  351 format(//' n, x, y(1-6, n):'/(i4,1pe13.5,0p5f10.5,1pe13.5))
c
       if(idgn76.ne.0) write(76) iform,nn,nrdtmd,nidtmd,nvar,nbcprv,
     *    (datout(i),i=1,nrdtmd),(ndtout(i),i=1,nidtmd),
     *    (x(n),(y(in11+i,n),i=1,nvar),n=1,nn),
     *    (bccoef(i),i=1,nbcprv)
c
      eam=eamean(ea,idiffus)
c
c  test for undefined solution (not a number)
c
      if(infnte(eam).ne.0) then
        if(istdpr.gt.0) write(istdpr,112) iter,(ea(k,1),k=1,ii),eam
	call store(bccofp,bccoef,nbcprv)
	eam=1.
	go to 60010
      end if
c
      if(idgbcs.ge.1.and.istdpr.gt.0) call dmpbcs
c
      if(itndgn.ne.1) then
        if(istdpr.gt.0) write(istdpr,112) iter,(ea(k,1),k=1,ii),eam
      else
c
        if(istdpr.gt.0) then
	  write(istdpr,111) iter,eam
          do 58 j=1,3
   58     write(istdpr,115) j,(ea(k,j),k=1,ii)
          do 58010 k=1,ii
58010     iea(k)=ea(k,3)+0.5
          write(istdpr,58091) (y(in1,iea(k))-y(in,iea(k)),k=1,ii)
58091     format(2x,1p7e11.3)
          do 58092 j=1,inc
          xxf=frcf(j)*x(nf(j))+(1-frcf(j))*x(nf(j)-1)
          xxl=frcl(j)*x(nl(j))+(1-frcl(j))*x(nl(j)+1)
58092     write(istdpr,116)
     *      j,nf(j),10.**x(nf(j)),nl(j),10.**x(nl(j)),
     *      x(nl(j)), frcl(j),
     *      xxf, 10.d0**xxf, xxl, 10.d0**xxl
        end if
c
      end if
c
c  test for properties of convergence, with a possible reset of 
c  solution. If so, continue iteration with no further reset.
c
      if(ntime.gt.0) then
        call tstcvg(x,y(in1,1),y(in,1),nn,iy,ea,dt,irscvg,icry)
        if(icry.ne.0) then
	  ucy=0.2d0
	  nucy1=iter+5
	  if(istdpr.gt.0) write(istdpr,'(/ '' nucy, ucy reset to'',
     *      i3,f10.5)') nucy1, ucy
          irscvg=0
	  go to 60
        end if
      end if
c
c  as a crude check for discontinuously growing convective core,
c  switch off mixing etc if separate convection zone in core has
c  developed. NEEDS LATER REFINEMENT, ALMOST SURELY!
c
      iccc=0
      do 58093 j=1,inc
      if(x(nf(j)).le.-2) iccc=iccc+1
58093 continue
c..      if(iccc.ge.2) then
      if(iccc.le.-2) then
	if(istdpr.gt.0) write(istdpr,58900) iccc
	inmixc=1
      end if
58900 format(//' ***** Warning: ',i2,' convective regions in core.'/
     *         '       Switching off core-composition iteration.')
c
c  test for output of changes
c
      if(istdpr.gt.0.and.idgtnr.ge.1.and.iter.le.idgtnr) then
        write(istdpr,352)
        do 59 n=1,nn
   59   write(istdpr,353) n,x(n),(dytnrk(j+ii*(n-1)),j=1,ii)
      end if
  352 format(//' n, x, changes in y(i):'/)
  353 format(i5,1pe13.5,10e11.3)
c
c  reset parameters for convective overshoot
c  In an attempt to stabilize solution, freeze zones for slow
c  convergence
c
      if(iter.le.10) then
        call setovs(x,y(in1,1),nn,iy,incovs)
      end if
c
c  test for breaking iteration if correction is too large, or X has
c  become substantially negative. In such cases, jump to halve
c  timestep
c
      if(ntime.ge.2.and.eam.ge.0.3) then
	if(istdpr.gt.0) write(istdpr,59081) eam
	go to 60010
      else if(iter.gt.iterlm.and.eam.gt.2*eamp) then
	itwarn=itwarn+1
	if(istdpr.gt.0) write(istdpr,59080) eam,eamp
	if(itwarn.eq.3) go to 60010
      else
	itwarn=0
      end if
      eamp=eam
c
c..      do 59010 n=1,nn
c..      if(y(in1+4,n).lt.-0.1.or.y(in1+4,n).gt.1.1) then
c..	write(istdpr,59090) y(in1+4,n), n, x(n)
c..	eam=1.d0
c..	go to 60010
c..      end if
c
      do 59010 n=1,nn
      if(y(in1+4,n).lt.-0.1.or.y(in1+4,n).gt.1.1) then
	n1=max0(1,n-2)
	n2=min0(n+2,nn)
	xnew=0
	do 59008 i=n1,n2
59008   xnew=xnew+y(in1+4,i)
        xnew=xnew/(n2-n1+1)
	if(istdpr.gt.0) write(istdpr,59091) y(in1+4,n), n, x(n), xnew
	if(xnew.le.-1.or.xnew.gt.1.1) then
	  eam=1.d0
	  go to 60010
        else
	  y(in1+4,n)=xnew
        end if
      end if
59010 continue
59080 format(/' ***** Warning in evolmain. Iteration not converging.',
     *  ' eam, eamp =',1p2e13.5)
59081 format(/' ***** Warning in evolmain. Iteration not converging.',
     *        '       eam =',1pe13.5)
59090 format(/' ***** Warning in evolmain. X =',1pe13.5,
     *  ' at n, x =',i5,e13.5)
59091 format(/' ***** Warning in evolmain. X =',1pe13.5,
     *  ' at n, x =',i5,e13.5/
     *        '       Replaced by ',e13.5)
c
c  test for setting composition in convective core
c  Note: change test to ignore very small radiative core
c
c  Also, for the time being do not call if previous model
c  had no mixed core. This deserves to be looked into,
c  for the few cases where a mixed core appears during evolution.
c
      if(istdpr.gt.0) write(istdpr,*) 
     *  'in test in evolmain, inmixc, inentc =', inmixc, inentc
      if(ntime.gt.0.and.inc.gt.0.and.nl(inc).ge.nn-2
     *  .and.nmxcor.ne.-1.and.inmixc.eq.0) then
c
c  make corrector step for treatment of convective-core composition
c
c  Note: this now includes call of mixcor, since we wish to
c  iterate repeatedly
c
        iprdcr=2
        call cmpcvc(x,y,in,in1,iy,nn,compc,icomp,iche30,icnocs,
     *    iheccs,dt,it,iprdcr,icrycm)
c
c  test for error in cmpcvc
c
	if(icrycm.lt.0) then
	  eam=1.d0
	  go to 60010
        end if
c
      end if
c
      if(itsxin.ge.2) then
	xhintp = xhint
	xtimep = age
	write(ctsxin,'(''it.'',i2)') iter
	call tstxin(x,y(in1,1),nn,iy,1,xhint,xhintp,age,xtimep,
     *    icry,ctsxin)
      end if
c
      if(nxtst1.lt.nxtst2.and.istdpr.gt.0) then
        write(istdpr,57099) 'y(in,.)',2,(n, 10.**x(n),y(in+4,n),
     *    y(in+iyche4-1,n),n=nxtst1,nxtst2,ndxtst)
        write(istdpr,57099) 'y(in1,.)',2,(n, 10.**x(n),y(in1+4,n),
     *    y(in1+iyche4-1,n),n=nxtst1,nxtst2,ndxtst)
      end if
c
      if(eam.lt.eps) then
	if(irczms.eq.0.and.istrt2.gt.0) then
	  irczms = 1
	  itczms = iter
	  eamp = 100*eam
          inrczm=1
	  call rczmsh(x,y(in1,1),y(in,1),nn,iy,3,icrczm,inrczm)
          iterlm=iter+5
        else
	  go to 61
        end if
      end if
c
c  end of iteration loop
c  *********************
c
   60 continue
c
c  test for convergence
c
60010 if(eam.gt.eps) then
        write(istdou,150) eam, eps
        if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,150) eam, eps
c..      	if(ntime.eq.86) stop 'Convergence problems at timestep 86'
	if(kdgrhb.lt.0) stop 'Convergence problems with physics errors'
c
c  test for maximum number of repeats reached
c
        if(nrpsol.gt.nrpmax) then
	  write(istdou,'(//'' Number of repeats ='',i2,
     *      '' exceeds maximum. Stop'')') nrpsol
	  if(istdou.ne.istdpr.and.istdpr.gt.0)
     *      write(istdpr,'(//'' Number of repeats ='',i2,
     *      '' exceeds maximum. Stop'')') nrpsol
	  stop 'Too many repeats with convergence problems'
	end if
	if(time0) then
c
c  test for increasing number of steps for initial iteration with
c  undercorrection
c
	  if(irnucy.eq.1) then
            write(istdou,154)
            if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,154)
            stop 'Stop time0, no convergence in evolmain'
	  else
	    irnucy=1
            nucy1=1+nit/2
	    ucy=ucyt
	    write(istdou,151) nucy1-1
	    if(istdpr.gt.0.and.istdpr.ne.istdou) 
     *        write(istdpr,151) nucy1-1
            write(istdou,105) ntime,age,agey,dt,dty
            if(istdpr.gt.0.and.istdpr.ne.istdou) 
     *        write(istdpr,105) ntime,age,agey,dt,dty
            if(idgn76.ne.0) close(76)
            if(idgn77.ne.0) close(77)
	    i=in
	    in=in1
	    in1=i
            in11=in1-1
            if(istcon.gt.0) stop 'Stop 2a, no convergence in evolmain'
	    go to 57
	  end if
c
c  tests for non-zero time
c
	else if(dt.lt.2*dtmn.and.inmixc.eq.1.and.inentc.eq.1) then
          write(istdou,154) 
          if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,154) 
          stop 'Stop 2, no convergence in evolmain'
	else if(dt.lt.2*dtmn.and.inentc.eq.0) then
	  iqcrsp = 0
          write(istdou,153) iqcrsp
          if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,153) iqcrsp
	  inentc = 1
	  inmixc=1
	  redfac=1.
          if(istcon.gt.0) stop 'Stop 2d, no convergence in evolmain'
	  if(iheccs.ne.0) stop 'Convergence problems'
	else 
c
c  test for problems with small convective core (the alpha Cen A syndrom)
c
	  icscbr=1
          iqcres=0
	  call setcbr(x,y(in,1),rhs,nn,iy,icscbr,iqcres)
	  if(iqcres.eq.1) then
	    inmixc = -1
	    redfac=1.
          else if(iqcres.eq.2) then
	    inmixc = 1
	    inentc = 1
	    redfac=1.
            write(istdou,153) iqcres
            if(istdpr.gt.0.and.istdpr.ne.istdou) 
     *        write(istdpr,153) iqcres
            if(istcon.gt.0) stop 'Stop 2c, no convergence in evolmain'
	    if(iheccs.ne.0) stop 'Convergence problems'
	  else if(inmixc.eq.0.and.dt.lt.2*dtmn) then
            write(istdou,153) iqcres
            if(istdpr.gt.0.and.istdpr.ne.istdou) 
     *        write(istdpr,153) iqcres
	    inmixc=1
	    inentc=1
	    redfac=1.
            if(istcon.gt.0) stop 'Stop 2b, no convergence in evolmain'
	    if(iheccs.ne.0) stop 'Convergence problems'
          else
	    redfac=0.5
	  end if
        end if
	age=age-dt
	call redtst(y(in,1),zk,y(in1,1),zk,dt,redfac,eam,ii,kk,nn,iy)
        write(istdou,158) redfac
        if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,158) redfac
	age=age+dt
        agey=age/syear
        dty=dt/syear
        write(istdou,105) ntime,age,agey,dt,dty
        if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,105) 
     *    ntime,age,agey,dt,dty
c
c  reset central boundary condition quantities from previous 
c  time step in bccoef
c
        call store(bccofp, bccoef, nbcprv)
c
c  set convective-core parameters based on convection zone in last
c  converged models (a somewhat flakey fix)
c
	qmxcor=qmxcp
	cqc=qmxcor
c  
c  reset flag for last model in sequence
c
	lastmd=ntime.eq.nt.or.age.ge.agefns*(1.d0-1.d-7)
c
        if(idgn76.ne.0) close(76)
        if(idgn77.ne.0) close(77)
        if(istcon.gt.0) stop 'Stop 2e, no convergence in evolmain'
c
c  reset undercorrection properties to initial values
c
        nucy1=nucy0
	ucy=ucyt
	if(istdpr.gt.0) write(istdpr,'(/ '' nucy, ucy reset to'',
     *    i3,f10.5)') nucy1, ucy
	go to 57
      end if
c
c  close diagnostic files
c
   61 if(idgn76.ne.0) close(76)
      if(idgn77.ne.0) close(77)
c
      nmxl=-20
      if(nmxl.gt.0) call tstmxl(x(nmxl),y(in1,nmxl),nmxl)
c
      if(iprhe4.ge.1.and.istdpr.gt.0) write(istdpr,'(a/(i5,1p5e15.7))') 
     *  '4He abundances after end of tnrkt iteration:',
     *  (n,x(n),y(in11+iyche4,n)-dprhe4,y(in11+iycc12,n),
     *  cvr(icvhe4,n)-dprhe4,cvr(icvc12,n),n=nphe41,nphe42,ndphe4)
c  set number of iterations
c
      itmax=it
c
      if(itsxin.ge.1) call tstxin(x,y(in1,1),nn,iy,1,xhint,xhintp,
     *  age,xtimep,icry,'pos.1')
c
c  reset composition from ymix when core composition has been 
c  modified as a result of convective core
c
      if(ntime.gt.0.and.inc.gt.0.and.nl(inc).ge.nn-2
     *  .and.nmxcor.ne.-1.and.inmixc.eq.0) then
	do 61005 n=1,nn
	do 61005 i=1,icomp
61005   y(in1+3+i,n)=ymix(4+i,n)
      end if
c
c  reset X to 1.d-10 where less than xhzlim
c  Also possibly reset Y (but only once) in case of He burning
c
c  Note that 4He composition is only reset if there are no reactions,
c  as determined by xmxrhe
c
      if(istdpr.gt.0) 
     *  write(istdpr,*) 'In evolmain, xmxrhe =',xmxrhe
      do 61010 n=1,nn
      if(y(in1+4,n).lt.xhzlm1.or.n.ge.nxhzer) then
	y(in1+4,n)=1.d-10
	if(x(n).gt.xmxrhe.and.iheccs.ne.0) then
	  if(n.eq.600.and.istdpr.gt.0) 
     *      write(istdpr,*) 'In evolmain, at n =',n,
     *      ' reset Y from ',y(in11+iyche4,n),' to',1.d0-1.d-10-zh(n)
	  y(in11+iyche4,n)=1.d0-1.d-10-zh(n)
        end if
      end if
c
c  in region of hydrogen exhaustion, with 4He burning reset zh(n)
c
      if(y(in1+4,n).le.2.d-10.and.iheccs.ne.0) then
        zh(n)=1.d0-y(in11+iyche4,n)-1.d-10
c..	if(n.eq.-600.or.abs(zh(n)-0.02d0).gt.1.d-8) then
c..	  write(istdpr,*) 'At n =',n,' set Z to',zh(n)
c..          write(istdpr,*) 'X, Y =', 1.d-10, y(in11+iyche4,n)
c..        end if
      end if
61010 continue
c
      istmxc=0
c
c  set nmxcor etc to zero when there is no convective core
c
      if(inc.eq.0.or.nl(inc).lt.nn-2.or.imixcr.eq.-1) then
	nmxcor=0
	qmxcor=0
	rmxcor=0
	cqc=0
c
c  else test for setting mixed-core boundary, if this was not set 
c  during iteration, and set flag istmxc for core mixing
c  with call of conmxc later.
c
      else if(inc.gt.0.and.nl(inc).eq.nn.and.nmxcor.eq.0) then
        if(ddrmix.eq.0) then
          iextrp=3
        else
          iextrp=10
        end if
	if(istdpr.gt.0) write(istdpr,*) 
     *    'Call mixcor after completion of iteration'
        call mixcor(x,y(in1,1),iy,nn,compc,iextrp,0)
	istmxc=1
	cqc=qmxcor
c
c  set extent of mixed region at core for initial time step
c
      else if(ntime.eq.0) then
        if(ddrmix.eq.0) then
c
c  modified from iextrp=1 on 13/10/98, for consistency 
c  with evolmain.dnnz.d.f
c
	  iextrp=3
        else
	  iextrp=10
        end if
        call mixcor(x,y(in1,1),iy,nn,compc,iextrp,0)
      end if
c
c  zero iter (for reasons none to clear)
c
      iter_conv=iter
      iter=0
c
c  test for setting he3 abundance in common/conpvr/
c#ai#  Note: this should not be necessary as cvr should
c#ai# have been set already through calling s/r rhs.
c#ai# However, keep for the moment.
c
      if(ifdhe3.ne.0.or.iche30.eq.2) then
        call setxh3(x,y(in1,1),cvr(3,1),nn,iy,icvrmx,iche30)
      end if
c
c  reset coefficients for central expansion
c
      call bcscen(x(nn),y(in1,nn),nn)
c
c  in case of a convective core mix material, except when computing
c  static models
c
      if(inc.gt.0.and.nt.gt.0) then
        do 66 j=1,inc
        nmx=nl(j)
        if(nmx.eq.nn) then
c
c  reset composition only for initial model, when He3
c  fudge has been used
c
          if(ifdhe3.ne.0) then
            ireset=1
            icomp0=2
	  else if(istmxc.eq.1) then
            ireset=1
            icomp0=icomp
          else
            ireset=0
            icomp0=1
          end if
      if(iprhe4.ge.1.and.istdpr.gt.0) write(istdpr,'(a/(i5,1p5e15.7))') 
     *  '4He abundances, before conmxc:',
     *  (n,x(n),y(in11+iyche4,n)-dprhe4,y(in11+iycc12,n),
     *  cvr(icvhe4,n)-dprhe4,cvr(icvc12,n),n=nphe41,nphe42,ndphe4)
          call conmxc(x,y(in1,1),iy,j,nn,icomp0,ireset)
c
	end if
c
c  mixing in convection zones that do not include the
c  centre, or if inmixc .ne. 0, indicating that mixing is
c  not dealt with in s/r rhs. Make separate call to mix variables in cvr
c  when treating CNO cycle.
c
        if((nmx.lt.nn.and.abs(y(in11+5,nmx)-xhs).gt.1.d-7)
     *       .or.inmixc.ne.0) then
          call conmix(x,y(in1,1),iy,j,nn,azt(3),icomp,1)
c..          if(icnocs.ge.1) then
c..            call conmix(x,cvr(4,1),icvrmx,j,nn,cvr(4,nn+1),icvcno,0)
c..          end if
        end if
   66   continue
      end if
c
      if(iprhe4.ge.1.and.istdpr.gt.0) write(istdpr,'(a/(i5,1p5e15.7))') 
     *  '4He abundances, after conmix:',
     *  (n,x(n),y(in11+iyche4,n)-dprhe4,y(in11+iycc12,n),
     *  cvr(icvhe4,n)-dprhe4,cvr(icvc12,n),n=nphe41,nphe42,ndphe4)
c
      if(nxtst1.lt.nxtst2.and.istdpr.gt.0) then
        write(istdpr,57099) 'y(in,.)',3,(n, 10.**x(n),y(in+4,n),
     *    y(in+iyche4-1,n),n=nxtst1,nxtst2,ndxtst)
        write(istdpr,57099) 'y(in1,.)',3,(n, 10.**x(n),y(in1+4,n),
     *    y(in1+iyche4-1,n),n=nxtst1,nxtst2,ndxtst)
      end if
c
      if(itsxin.ge.1) then
	xhintp = xhint
	xtimep = age
	call tstxin(x,y(in1,1),nn,iy,1,xhint,xhintp,age,xtimep,
     *    icry,'pos.2')
      end if
c
c  test whether envelope composition has changed and, if so, 
c  reset opacity tables if spline tables are used
c  Note that this must be changed at this point, to 
c  prepare for calls of mshstr and prtsol.
c
      xhs=y(in1+4,1)
      if(abs(xhs-xhsopc).ge.1.d-5.and.
     *  iwdop0.ne.1.and.iwdop0.ne.8.and.iwdop0.ne.9) then
        if(istdpr.gt.0) write(istdpr,167) xhsopc, xhs
        xhsopc=xhs
        rewind inopc
        call opcscf(iwdopc, iopacm)
      end if
c
c  test for resetting angular velocity 
c
      if(isprot.gt.0) then
c
c  initialize angular velocity and moment of inertia
c
	if(ntime.eq.0) then
	  introt=1
        else
	  introt=0
	end if
	call setrot(x,y(in1,1),y(in,1),iy,nn,velrot,introt)
      end if
c
c  output on terminal
c
      if(istdpr.ne.istdou) then
        y1=1.d11*10.d0**y(in11+1,1)
        y4=1.d33*10.d0**y(in11+4,1)
        write(istdou,113) iter_conv,eam,y1,y4
      end if
c
c  print the solution, economy output
c#ai# For the moment, and for testing, print same variables as
c#ai# previously. Later to be changed to just active variables.
c
      if(istdpr.gt.0) then
        write(istdpr,120)
        write(istdpr,122)
        do 67 n=1,nn,intr
        y1=1.d11*10.d0**y(in11+1,n)
        y3=10.d0**y(in11+3,n)
        y4=1.d33*(10.d0**y(in11+4,n)-alshft)
   67   write(istdpr,123) 
     *    n,x(n),y1,y(in11+2,n),y3,y4,y(in11+5,n),cvr(3,n)
      end if
c
c  test for boundary condition diagnostics
c
      if(idgbcs.ge.1.and.istdpr.gt.0) call dmpbcs
c
      ifprt=1
c
c  test for initial stretching and repeated solution on new mesh
c  for initial timestep
c
      if(ntime.eq.0.and.istrtc.ge.1.and.nmsh.lt.nrpeat) then
        nmsh=nmsh+1
        if(nmsh.eq.nrpeat) intork=2
        if(istdpr.gt.0) then
          write(istdpr,*)
     *      'wx,epsr,wr,wpt,wxh, wx3,wdgr,wmshos,intork,intvnt,iprmsh,',
     .      'icngrm,nnt,iwdgrd'
          write(istdpr,*)
     *      wx,epsr,wr,wpt,wxh, wx3,wdgr,wmshos,intork,intvnt,iprmsh,
     .      icngrm,nnt,iwdgrd
	end if
        call mshstr(x,y,in1,in,iy,nn,nn,rhs,resrc)
        if(istdpr.gt.0) write(istdpr,130)
        nucy1=1
        go to 53
      end if
c
c  output model to file
c  set parameters for output
c
      amms=amsun*am
      arrs=1.d11*10.d0**y(in1,1)
      alls=1.d33*(10.d0**y(in1+3,1)-alshft)
      sigstd=sigstr
c
c  test for consistency of modeeq and ivreos.
c
      if(ivreos.ne.modeeq.and.ivreos.ne.9) then
        write(istdou,125) ivreos, modeeq
        if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,125) 
     *    ivreos, modeeq
      end if
c
c  case number, added 23/11/1982.
c
c  note also that before this date qqbc has length 4, after this date
c  length 5. thus the position of patmos and alamop in datmod has been
c  shifted by 1. the position of the remaining variables is unchanged.
c  *******************************************************************
c
      icase=iqqbc+10*iwdop0+100*(iopacm+2*iopatm)
     * +1 000*icsrad+10 000*ihvz
     * +100 000*ivreos+1 000 000*ifdgop+100 000 000*ivreng
c
c  include version number in icase
c
      icase=icase+iversn*10 000 000
c
c  set datout by writing to and reading from scratch file
c  unless this is run with isetos = 1
c
      if(isetos.ne.1.and.isetos.ne.-1) then
c
c  set qqbcpr, depending on ntaubc
c
        if(ntaubc.le.1) then
          qqbcpr(1)=albdbc
          call zero(qqbcpr(2),4)
        else
          call store(qqbc,qqbcpr,5)
        end if
c
c  Use write and read to set first part of datout
c
        rewind 91
        write(91) zhinit,y(in1+4,1),alfa,etac,phc,
     .    agehe3,sigstd,wx,wpt,wr,wx3,tmnbc,tmxbc,sbcfct,qqbcpr,
     .    patmos,alamop,agey,amms,arrs,alls,(dtopfg(i),i=1,4)
        rewind 91
        read(91) (datout(i),i=1,29)
c
c  for remaining parameters in datout, use more explicit setting
c  (note that by a stupid mistake, datmod(34) was previously set
c  to the non-existing parameter timn)
c
        datout(30)=tsmn
        datout(31)=tstr
        datout(32)=rhsmn
        datout(33)=rhsmx
        datout(34)=0
        datout(35)=timx
        datout(36)=rhimn
        datout(37)=rhimx
        datout(38)=sigstr
        datout(39)=zatmop
c
        call store(zab,datout(41),10)
c
        datout(51)=thte
        datout(52)=fcno
        datout(53)=xzer3
        datout(54)=xrz12
        datout(55)=xrz13
        datout(56)=xrz14
        datout(57)=xrz16
c
        datout(61)=clcovs
        datout(62)=cldovs
        datout(63)=alphos
        datout(65)=albdbc
        datout(66)=eps
c
        if(nvar.le.10) then
          call store(thetad,datout(71),nvar)
        else
          if(istdpr.gt.0) write(istdpr,135) nvar
          call store(thetad,datout(71),10)
        end if
c
        datout(81)=epsr
        datout(82)=wdgr
        datout(83)=wmshos
        datout(84)=drmsos
        datout(85)=wmshcc
        datout(86)=dqlmsc
        datout(87)=wxh
c
        datout(91)=dt0
        datout(92)=dtmx
        datout(93)=dtmn
        datout(94)=dymx
        datout(95)=xdtden
        datout(96)=aldtrt
	datout(97)=dtfcrd
	datout(98)=dtceps
c
	if(ifdeps.ne.0) then
	  datout(111)=epsfdg
	  datout(112)=qepsf1
	  datout(113)=qepsf2
        end if
c
	if(iturpr.gt.0) then
	  datout(115)=tprfct
        end if
c
        datout(121)=amu
        datout(122)=ame
        datout(123)=clight
        datout(124)=planck
        datout(125)=boltzm
        datout(126)=cgrav
        datout(127)=amsun
        datout(128)=echar
        datout(129)=ergev
        datout(130)=syear
c
	if(isprot.gt.0) then
	  datout(135)=velrot
	  datout(136)=omgrot(1)
        end if
c
	datout(140)=flsatm
c
	if(iqfit.lt.0) call store(qqbc,datout(141),iqqbc)
c
c  parameters for new overshoot formulation (5/1/00)
c
        datout(151)=alpove
        datout(152)=alpovc
c
c  shift in luminosity (added 2/10/02)
c
	datout(155)=alshft
c
c  set ndtout explicitly
c
        call izero(ndtout,nidtmd)
        ndtout(1)=icase
        ndtout(2)=istrtc
        ndtout(3)=ntaubc
c
        ndtout(4)=iche3
        ndtout(5)=icnocs
        ndtout(6)=ihe3bc
        ndtout(7)=nspcmx
        ndtout(8)=idthm
        ndtout(9)=isethv
        ndtout(10)=ianhe0
        ndtout(11)=iomfll
	ndtout(12)=ivropc+100*iwdopc
        ndtout(13)=iscren
        ndtout(14)=icnvos
        ndtout(15)=jcnvos
        ndtout(16)=inentc
        ndtout(17)=icncbc
        ndtout(18)=iwdgrd
        ndtout(19)=icvcno
c
        ndtout(24)=nn
        ndtout(25)=ntime
c
	ndtout(26)=inc
	ndtout(27)=imixcr
	ndtout(30)=ifdeps
	ndtout(31)=iconcs
	ndtout(32)=imxlng
	ndtout(33)=iturpr
	ndtout(35)=isprot
	ndtout(41)=invers
	ndtout(42)=ivteos
	ndtout(45)=iqfit
c
c  parameters for new overshoot formulation (5/1/00)
c
        ndtout(51)=icsove
        ndtout(52)=icsovc
c
	if(qmxcor.gt.0) ndtout(55)=inmixc
c
c  total number of iterations
c
	ndtout(56)=itmax
c
	ndtout(61) = iheccs
c
      end if
c
      if(iterlr.eq.0) then
c
        if(iastr.eq.1) 
     *    write(idsevl) iform,nn,nrdtmd,nidtmd,nvar,nbcprv,
     *    (datout(i),i=1,nrdtmd),(ndtout(i),i=1,nidtmd),
     *    (x(n),(y(in11+i,n),i=1,nvar),n=1,nn),
     *    (bccoef(i),i=1,nbcprv)
	call flush(idsevl)
      else
c
c  temporary output during iteration for alfa and xh
c
c  first model in sequence
c
        if(nt.ne.0.and.ntime.eq.0)  then
          rewind idstm1
          write(idstm1) iform,nn,nrdtmd,nidtmd,nvar,nbcprv,
     *      (datout(i),i=1,nrdtmd),(ndtout(i),i=1,nidtmd),
     *      (x(n),(y(in11+i,n),i=1,nvar),n=1,nn),
     *      (bccoef(i),i=1,nbcprv)
        end if
c
c  current model in sequence
c
        write(idstm2) iform,nn,nrdtmd,nidtmd,nvar,nbcprv,
     *    (datout(i),i=1,nrdtmd),(ndtout(i),i=1,nidtmd),
     *    (x(n),(y(in11+i,n),i=1,nvar),n=1,nn),
     *    (bccoef(i),i=1,nbcprv)
c
        ipmlr=ipmlr+1
c
      end if
c
c  flag for setting oscillation variables
c
      ist=0
      if(iterlr.eq.0.and.(istosc.eq.1.or.
     *  (ioscpr.eq.1.and.lastmd))) ist=1
c
c  test for complete output
c
      comout=icmout.eq.1.or.ist.eq.1
c
  705 if(comout)
     .  call prtsol(x,y(in1,1),y(in,1),zk,ap,aq,rhs,nn,iy,intr,ifprt,
     .  ist,icasex,iform,agey,datout,ndtout,ndsmth)
c
      if(iprhe4.ge.1.and.istdpr.gt.0) write(istdpr,'(a/(i5,1p5e15.7))') 
     *  '4He abundances, after prtsol:',
     *  (n,x(n),y(in11+iyche4,n)-dprhe4,y(in11+iycc12,n),
     *  cvr(icvhe4,n)-dprhe4,cvr(icvc12,n),n=nphe41,nphe42,ndphe4)
c
c
c  output boundaries of convection zones
c
      if(inc.gt.0.and.istdpr.gt.0) then
c
        write(istdpr,117)
c
        write(istdpr,*) 'inc, nf,nl,frcf,frcl', inc, nf,nl,frcf,frcl
c
        do 713 j=1,inc
        nff=nf(j)
        nll=nl(j)
c
c  upper edge
c
        frc=frcf(j)
        frc1=1-frc
        write(istdpr,*) ' Convection, j, nff, nll, frc, frc1, in11 =',
     *    j,nff,nll, frc, frc1, in11
c
        if(nff.gt.1) then
          y1=1.d11*10.d0**(frc*y(in11+1,nff)+frc1*y(in11+1,nff-1))
          y2=              frc*y(in11+2,nff)+frc1*y(in11+2,nff-1)
          y3=      10.d0**(frc*y(in11+3,nff)+frc1*y(in11+3,nff-1))
          y4=1.d33*10.d0**(frc*y(in11+4,nff)+frc1*y(in11+4,nff-1))
          y5=              frc*y(in11+5,nff)+frc1*y(in11+5,nff-1)
          y6=              frc*y(in11+6,nff)+frc1*y(in11+6,nff-1)
c
          xx=        frc*x(nff)+frc1*x(nff-1)
c
	else
c
          y1=1.d11*10.d0**y(in11+1,nff)
          y2=             y(in11+2,nff)
          y3=      10.d0**y(in11+3,nff)
          y4=1.d33*10.d0**y(in11+4,nff)
          y5=             y(in11+5,nff)
          y6=             y(in11+6,nff)
c
          xx=             x(nff)
c
	end if
c
        write(istdpr,118) j,nff,xx,y1,y2,y3,y4,y5,y6
c
	xrcf(j)=xx
c
c  lower boundary
c
        frc=frcl(j)
        frc1=1-frc
c
        if(nll.lt.nn) then
          y1=1.d11*10.d0**(frc*y(in11+1,nll)+frc1*y(in11+1,nll+1))
          y2=              frc*y(in11+2,nll)+frc1*y(in11+2,nll+1)
          y3=      10.d0**(frc*y(in11+3,nll)+frc1*y(in11+3,nll+1))
          y4=1.d33*10.d0**(frc*y(in11+4,nll)+frc1*y(in11+4,nll+1))
          y5=              frc*y(in11+5,nll)+frc1*y(in11+5,nll+1)
          y6=              frc*y(in11+6,nll)+frc1*y(in11+6,nll+1)
c
          xx=        frc*x(nll)+frc1*x(nll+1)
c
	else
c
          y1=1.d11*10.d0**y(in11+1,nll)
          y2=             y(in11+2,nll)
          y3=      10.d0**y(in11+3,nll)
          y4=1.d33*10.d0**y(in11+4,nll)
          y5=             y(in11+5,nll)
          y6=             y(in11+6,nll)
c
          xx=             x(nll)
c
	end if
c
	write(istdpr,119) j,nll,xx,y1,y2,y3,y4,y5,y6
c
	xrcl(j)=xx
  713   continue
c
      end if
c
c  test for special output from user-supplied routine spcout
c
      if(ispcpr.gt.0) then
        call spcout(x,y(in1,1),y(in,1),nn,iy,intr,agey,datout,ispcpr)
      end if
c
      if(isetos.eq.1.or.isetos.eq.-1) then
c
c  test for reading a sequence of models, for nt .gt. 1
c
        if(nmdtrl.lt.nmdmax) then
          nmdtrl=nmdtrl+1
          go to 44
        else
          go to 20
        end if
      end if
c
c  check whether this is last  model
c  *********************************
c
      if(lastmd) go to 85
c
      if(ntime.lt.0) then
        write(istdou,*) ' Error. ntime .lt.0'
        if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,*) 
     *    ' Error. ntime .lt.0'
        stop 'Stop 3 in evolmain'
      else if(ntime.gt.0) then
c
c  subsequent time steps. reset dt
c
        call setdt(x,y(in1,1),y(in,1),dymx,dtmx,dtmn,age,agefns,dt,dtp,
     *    xdtden,aldtrt,dtfcrd,dtceps,iy,nn)
c
c  set new trial solution by extrapolation
c
c  to stabilize solution, do not allow extrapolation for the first
c  few time steps
c
	if(ntime.le.3) then
	  etares=0.
        else
	  etares=eta
        end if
c
        call pushy(x,y(in1,1),zk,y(in,1),zk,dt,dtp,etares,ii,kk,nn,iy)
c
c..     write(6,72091) ' After pushy',in1,in,((y(in1+i-1,n),i=1,2),
c..     *    (y(in+i-1,n),i=1,2),n=1,20)
c..72091   format(/a/' in1 =',i4,'  in =',i4/
c..     *    ' y(in1), y(in):'//(1p4e15.7))
      end if
c
c  change mesh
c  ***********
c
      if(istrtc.ge.1) 
     *  call mshstr(x,y,in1,in,iy,nn,nn,rhs,resrc)
c
      if(iprhe4.ge.1.and.istdpr.gt.0) write(istdpr,'(a/(i5,1p5e15.7))') 
     *  '4He abundances after call of mshstr:',
     *  (n,x(n),y(in11+iyche4,n)-dprhe4,y(in11+iycc12,n),
     *  cvr(icvhe4,n)-dprhe4,cvr(icvc12,n),n=nphe41,nphe42,ndphe4)
      if(xtst1.gt.xtst2) then
        do n=1,nn
	  qq=10.**x(n)
	  if(qq.ge.xtst1) nxtst1=n
	  if(qq.ge.xtst2) nxtst2=n
        end do
      else
	nxtst1=0
	nxtst2=0
      end if
c
      if(itsxin.ge.1) then
	xhintp = xhint
	call tstxin(x,y(in1,1),nn,iy,1,xhint,xhintp,age,xtimep,
     *    icry,'pos.3')
      end if
c
      if(nxtst1.lt.nxtst2.and.istdpr.gt.0) then
        write(istdpr,57099) 'y(in,.)',4,(n, 10.**x(n),y(in+4,n),
     *    y(in+iyche4-1,n),n=nxtst1,nxtst2,ndxtst)
        write(istdpr,57099) 'y(in1,.)',4,(n, 10.**x(n),y(in1+4,n),
     *    y(in1+iyche4-1,n),n=nxtst1,nxtst2,ndxtst)
      end if
c
      if(ntime.gt.0) go to 80
c
c  ntime = 0.
c
c  output first model to file
c
      if(ifstr.eq.1) write(idsefl) iform,nn,nrdtmd,nidtmd,nvar,nbcprv,
     *    (datout(i),i=1,nrdtmd),(ndtout(i),i=1,nidtmd),
     *    (x(n),(y(in11+i,n),i=1,nvar),n=1,nn),
     *    (bccoef(i),i=1,nbcprv)
c
c  reset ii1,ii2,ii3,y(in+.,.) and time0 after initial solution
c
c  Note: this is also entry point for icntsl ge 1
c  this is clearly bad style, and should be fixed up.
c
   76 time0=.false.
      intork=2
      ifdhe3=0
c
      nucy1=nucy+1
      ii1=3
      ii2=1
      ii3=ispxx3 + ispcno + isphec
      ii=ii1+ii2+ii3
      icomp=ii3
      jn=in+nvar-1
      k=in11
      do 78 j=in,jn
      k=k+1
      do 78 n=1,nn
   78 y(j,n)=y(k,n)
c
c  test for resetting case of fudge of energy generation
c  after initial solution
c
      if(ifdeps.lt.0) then
	ifdep1=0
	if(istdpr.gt.0) write(istdpr,*) 'Reset ifdep1 to 0'
      end if
c
c  ensure that envelope composition is uniform and unchanged 
c  after pushy and mshstr
c  Precisely how this is best done is not clear, but at least
c  region corresponding to envelope opacity plus a bit should
c  have the envelope opacity composition.
c
c  Note that extrapolation is problematic in convective core,
c  where we extrapolate uniform composition. This can hardly
c  be helped, however.
c
c  ***** Note: until 5/8/95, this was (erroneously) done in all
c        cases, compromising calculations with diffusion
c
   80 if(iwdop0.ne.1.and.iwdop0.ne.8.and.iwdop0.ne.9) then
	xhsprv = y(in+4,1)
        do 82 n=1,nn
        if(y(in1+2,n).gt.tstr+0.2) go to 83
        nres = n
   82   y(in+4,n)=xhsopc
   83   continue
c
        if(istdpr.gt.0) write(istdpr,*) 'X reset from', 
     *    xhsprv,' to', xhsopc, ' outside n =',nres
      end if
c
c  ensure that X is never negative, nor exceeds surface value
c  In case of 4He burning, also check that Z is never smaller 
c  than initial surface value
c
      do 84 n=1,nn
      if(y(in1+4,n).le.0) then
	if(istdpr.gt.0) write(istdpr,168) n, y(in1+4,n)
	y(in1+4,n)=1.d-10
      end if
      if(y(in+4,n).le.0) then
	if(istdpr.gt.0) write(istdpr,168) n, y(in+4,n)
	y(in+4,n)=1.d-10
      end if
      if(y(in1+4,n).gt.y(in1+4,1)+1.d-13) then
	if(istdpr.gt.0) write(istdpr,169) n, y(in1+4,n),y(in1+4,1)
	y(in1+4,n)=y(in1+4,1)
	if(iheccs.ne.0.and.x(n).gt.xmxrhe) 
     *	  y(in1+iyche4-1,n)=y(in1+iyche4-1,1)
      end if
      if(y(in+4,n).gt.y(in+4,1)+1.d-13) then
	if(istdpr.gt.0) write(istdpr,169) n, y(in+4,n),y(in+4,1)
	y(in+4,n)=y(in+4,1)
	if(iheccs.ne.0.and.x(n).gt.xmxrhe) 
     *	  y(in+iyche4-1,n)=y(in+iyche4-1,1)
      end if
c
      if(iheccs.ne.0) then
        if(1.d0-y(in1+4,n)-y(in1+iyche4-1,n).lt.zhinit-1.d-12) then
	  if(istdpr.gt.0) write(istdpr,*) 'Reset Y. Y_old, Y_new:',
     *      y(in1+iyche4-1,n),1.d0-y(in1+4,n)-zhinit
	  y(in1+iyche4-1,n)=1.d0-y(in1+4,n)-zhinit
	  zh(n)=zhinit
        end if
        if(1.d0-y(in+4,n)-y(in+iyche4-1,n).lt.zhinit) then
	  y(in+iyche4-1,n)=1.d0-y(in+4,n)-zhinit
        end if
      end if
c
   84 continue
c
c  reset storage indices
c
      i=in1
      in1=in
      in=i
c
      if(nxtst1.lt.nxtst2.and.istdpr.gt.0) then
        write(istdpr,57099) 'y(in,.)',5,(n, 10.**x(n),y(in+4,n),
     *    y(in+iyche4-1,n),n=nxtst1,nxtst2,ndxtst)
        write(istdpr,57099) 'y(in1,.)',5,(n, 10.**x(n),y(in1+4,n),
     *    y(in1+iyche4-1,n),n=nxtst1,nxtst2,ndxtst)
      end if
c
      go to 50
c
c  now model sequence is  finished
c  *******************************
c
c  test for internal iteration for rs and ls
c
   85 if(iterlr.eq.0) go to 89
c
      if(age.le.0.99*agefns-1) then
c
c  diagnostics for final model age too small
c
        write(istdou,171) agey,agefin
        if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,171) 
     *    agey,agefin
        go to 87
      end if
c
c  find errors in ls and rs
c
      dls=log(alsfin/alls)
      drs=log(rsfin/arrs)
      write(istdou,170) itlr,alfa,sbcfct,y(in1+4,1),alls,arrs,dls,drs
      if(istdpr.gt.0.and.istdpr.ne.istdou) 
     *  write(istdpr,170) itlr,alfa,sbcfct,y(in1+4,1),alls,arrs,dls,drs
c
c  test for convergence
c
      if(abs(dls).ge.epslr.or.abs(drs).ge.epslr) then
c  test for excessive number of iterations
        if(itlr.ge.nitlr) then
          write(istdou,172) nitlr
          if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,172) nitlr
        else
c  find new alfa or sbcfct and xxh
          dalfa=derlr(1,1,iterlr)*dls+derlr(1,2,iterlr)*drs
          dxxh=derlr(2,1,iterlr)*dls+derlr(2,2,iterlr)*drs
c
          if(iterlr.le.2.or.iterlr.ge.5) then
            alfa=alfa*exp(dalfa)
          else
            sbcfct=sbcfct*exp(dalfa)
          end if
c
          xxh=y(in1+4,1)*exp(dxxh)
          if(iterlr.le.2.or.iterlr.ge.5) then
            write(istdou,175) dalfa,dxxh,alfa,xxh
            if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,175) 
     *        dalfa,dxxh,alfa,xxh
          else
            write(istdou,176) dalfa,dxxh,sbcfct,xxh
            if(istdpr.gt.0.and.istdpr.ne.istdou) write(istdpr,176) 
     *        dalfa,dxxh,sbcfct,xxh
          end if
          if(istdpr.gt.0) write(istdpr,177)
          go to 21001
c
        end if
c
      end if
c
c  end iteration for alfa, xxh. output final model on unit idsevl if
c  iastr = 1.
c
   87 if(istdpr.gt.0) write(istdpr,177)
      if(iastr.eq.1) then
c
        if(nt.eq.0) then
c  when nt  = 0 the model is already in y.
          write(idsevl) iform,nn,nrdtmd,nidtmd,nvar,nbcprv,
     *      (datout(i),i=1,nrdtmd),(ndtout(i),i=1,nidtmd),
     *      (x(n),(y(in11+i,n),i=1,nvar),n=1,nn),
     *      (bccoef(i),i=1,nbcprv)
c
c  store datout in datmod for later call of prtsol
c
          call  store(datout,datmod,nrdtmd)
          call istore(ndtout,ndtmod,nidtmd)
c
        else
c
c  otherwise read in last model sequence from  temporary storage
c  on unit idstm2
c
          rewind idstm2
c  last sequence is models ipmlr0+1, ..., ipmlr
          if(istdpr.gt.0) write(istdpr,180)
          k=1
          do 88001 i=1,ipmlr
          read(idstm2,end=88002) iform,nn,nrdtmr,nidtmr,nvarrd,nbccf,
     *      (datmod(j),j=1,nrdtmr),(ndtmod(j),j=1,nidtmr),
     *      (x(n),(y(j,n),j=1,nvarrd),n=1,nn),
     *      (bccoef(j),j=1,nbccf)
          if(i.gt.ipmlr0) k=2
          if(istdpr.gt.0) write(istdpr,182) 
     *      i,datmod(22),datmod(24),datmod(25),pmodu3(k)
          if(k.eq.2) write(idsevl) iform,nn,nrdtmr,nidtmr,nvar,nbccf,
     *      (datmod(j),j=1,nrdtmd),(ndtmod(j),j=1,nidtmd),
     *      (x(n),(y(j,n),j=1,nvar),n=1,nn),
     *      (bccoef(j),j=1,nbccf)
88001     continue
c
88002     in11=0
c
        end if
c
      end if
c
c  test for setting of oscillation variables
c
      if(ioscpr.eq.1) then
        in1=in11+1
        in=ivarmx+2-in1
        agey=datmod(22)
c
c  reset inc to -1, to force setting of convection zone boundaries 
c  in s/r rhs
c
        inc=-1
c
        call prtsol(x,y(in1,1),y(in,1),zk,ap,aq,rhs,nn,iy,intr,
     *    1,1,icasex,iform,agey,datmod,ndtmod,ndsmth)
c
        comout=.true.
c
      end if
c
c  output last model
c
   89 if(ilstr.eq.1) write(idsefl) iform,nn,nrdtmd,nidtmd,nvar,nbcprv,
     *    (datout(i),i=1,nrdtmd),(ndtout(i),i=1,nidtmd),
     *    (x(n),(y(in11+i,n),i=1,nvar),n=1,nn),
     *    (bccoef(i),i=1,nbcprv)
      go to 20
c
c  end of input
c  ************
c
   90 continue
      if(ibegin.eq.0) then
	write(istdou,159)
	if(istdpr.gt.0.and.istdou.ne.istdpr) write(istdpr,159)
        stop 'Stop 4 in evolmain'
      end if
      if(comout.and.nscfil(idssum).and.istdpr.gt.0) then
c  print global values from d/s idssum
        close(idssum) 
        call openf(idssum,'o','f')
        write(istdpr,160)
        call cpychr(idssum,istdpr)
c
      end if
c  scan modes on file
      if(iastr.eq.1.and.istdpr.gt.0) then
        close(idsevl) 
        call openf(idsevl,'o','u')
        call scnfil(idsevl,istdpr,1)
      end if
c
      stop
 1001 format(//' **** error in evolution programme.'/
     *         '      Coulomb effects not implemented in this version,',
     *              ' ivreos =',i3/
     *         '      modeeq reset from',i3,'  to 0')
 1002 format(//' **** error in evolution programme.'/
     *         '      modeeq =',i3,'  not allowed.')
  100 format(//' the trial solution:'/)
  101 format(//1x,75('*')//'    Old version number =',i4,
     *                     '    New version number =',i5//1x,75('*')//)
  102 format(//' Log in standard name: evol-file.log')
  104 format(//' iterate for atmospheric pressure. patmos =',
     *  1pe13.5,'   new alamop =',0pf10.5)
  105 format(///40x,50(1h*)//
     .  ' the ',i4,'-th timestep. the age = ',1pe13.5,
     .  ' seconds = ',e13.5,' years'/22x,' timestep =',e13.5,
     .  ' seconds = ',e13.5,' years.'/)
  106 format(///' use basic EFF equation of state')
  107 format(///
     *  ' use EFF equation of state with consistent Coulomb effect')
  108 format(///
     *  ' use EFF equation of state with inconsistent Coulomb effect')
  109 format(///' anhe0 reset from',f10.5,'  to 6.')
  110 format(//' set heavy element quantities suitable for simple ',
     *  'GONG models'/
     *  ' az =',f10.5,'  anh0 =',f10.5,'   anhe0 =',f10.5)
  111 format(///' iteration no',i3,' eam =',1pe13.3/1x,34(1h+)/
     .  /' ea:')
  112 format(' iteration no',i3,' eam =',1pe11.3/
     *  ' ea(1-ii):',7e11.3/(10x,7e11.3))
  113 format(/' end of iteration after ',i3,' steps. eam =',1pe11.3/
     *  ' rs =',e13.5,' cm.  Ls =',e13.5,' ergs/sec')
  115 format(i2,1p7e11.3/(2x,7e11.3))
  116 format(/' the ',i2,'. convective zone goes from n =',i4,
     .  ',q =',1pe12.5,' to n =',i4,', q =',e12.5/
     .  ' Last x =',e15.7,'  frcl =',0pf12.7/
     .  ' Actual range: x =',1pe13.5,' q =',e12.5,
     .  ' to x =',e13.5,' q =',e12.5)
  117 format(//' variables at edges of convection zones.'//
     *  25x,'n',7x,'log q',12x,'r',11x,'log f',10x,'t',14x,'l',
     *  13x,'x',10x,'x3'/)
  118 format(/' upper edge of zone',i3,i5,1p2e15.7,0pf12.7,1p2e15.7,
     *  0pf12.7,1pe13.5)
  119 format(' lower edge of zone',i3,i5,1p2e15.7,0pf12.7,1p2e15.7,
     *  0pf12.7,1pe13.5)
  120 format(//' the solution'/1x,12('*')/)
  122 format(///' shortened output.'/
     *  '   n',7x,'log q',12x,'r',11x,'log f',10x,'T',14x,'L',
     *  13x,'X',10x,'X3'/)
  123 format(i5,1p2e15.7,0pf12.7,1p2e15.7,0pf12.7,1pe13.5)
  125 format(//' *** error in setting case number. ivreos =',i3,
     *  ' differs from modeeq =',i3)
  127 format(//
     *  ' Setting convective core boundary for continuing solution'/)
  130 format(//' repeat solution at time 0 on new mesh'/  
     *         ' -------------------------------------'//)
  132 format(///40x,50(1h*)//
     .  ' Print solution from file at the ',i4,'-th timestep.'/
     .  ' The age = ',1pe13.5,' seconds = ',e13.5,' years')
  135 format(/' ***** Error in setting datout. nvar = ',i3,
     *  ' exceeds 10.'/
     *  '             Only thetad(1 - 10) stored in datout')
  136 format(//' az = ',f10.3,' anz0 =',f10.5)
  141 format(//' mass =',f10.5,'  solar masses'//
     *  ' alfa =',f12.6//' ix =',i3,' ix3 =',i3//
     *  ' surface x =',f10.7,'    z =',f10.7,'  agehe3 =',1pe13.5,
     *  ' years')
  142 format(/' simple surface conditions.'/
     *  ' sbcfct =',f10.5,'   albdbc =',f10.5)
  143 format(/' parameters in atmospheric solution:'/
     *  ' ntau =',i4,' taumin =',e13.5,'  taumax =',e13.5,
     *  ' sbcfct =',0pf10.5/
     *  ' q(1-7) =',5f12.5/(9x,5f12.7))
  144 format(//(' theta(',i2,') =',f10.6))
  145 format(/' log(opacity) artificially increased by',f10.5,
     *  ' .  corresponding opacity factor =',f10.5)
  146 format(//' log(opacity) fudged by adding ',f8.5,
     *  '*exp(-((log t -',f8.5,')/',f8.5,')**2)')
  147 format(//' log(opacity) fudged by adding ',f8.5,
     *  '*exp(-(dtl/',f8.5,')**2)  where dtl = min(log t - ',f8.5,
     *  ', 0,',f8.5,' - log t)')
  148 format(//' log(opacity) fudged by adding ',
     * 'delta log kappa(log T) interpolated from file'/1x,a/)
  150 format(//' ***** Error. Iteration not converged.'/
     *       '       eam = ',1pe11.3,'  .gt. eps =',e11.3)
  151 format(/' Change nucy0 to',i3,' and retry')
  153 format(/ '       As last desperate attempt, ',
     *  'try without mixing in convective regions,'/
     *         '       and without entropy terms in energy equation'/
     *         '       case: iqcres =',i3)
  154 format(/ '       Computation terminated.'/)
  155 format(/ '       As last even more desperate attempt, ',
     *         'try without entropy terms in energy equation'/)
  158 format(//'       Try repeating solution,',
     *  ' reducing timestep by factor',f8.3//)
  159 format(///' ***** Error in parameter inputs')
  160 format(//1x,78(1h*)/)
  167 format(//' envelope opacity X has been reset from',f12.7,
     *  '  to',f12.7/)
  168 format(' *** Warning. At n =',i5,' X =',1pe11.3,
     *  ' X reset to 1.d-10')
  169 format(' *** Warning. At n =',i5,' X =',1pe13.5,
     *  ' exceeds surface value =',e13.5/
     *       '     X reset to surface value')
  170 format(///1x,80('*')//
     *  ' iteration for alfa or sbcfct and xh to get',
     *  ' ls and rs'//' iteration no',i3,'  alfa =',f12.7,
     *  '  sbcfct =',f12.5,'  xhs =',f12.7/
     *  '  ls =',1pe15.7,'  rs =',e15.7//
     *  ' relative error in ls =',e13.5,'  in rs =',e13.5)
  171 format(///1x,80('*')//' in iteration for alfa and xh ',
     *  ' final age =',1pe13.5,' years is smaller than agefin =',e13.5,
     *  ' years.'/' iteration abandoned')
  172 format(//1x,10('*'),' iteration failed to converge in less than',
     *  i4,' iterations')
  175 format(//' relative change in alfa =',1pe13.5,'  and in xhs =',
     *  e13.5//' new alfa =',0pf12.7,'  new xhs =',f12.7)
  176 format(//' relative change in sbcfct =',1pe13.5,'  and in xhs =',
     *  e13.5//' new sbcfct =',0pf12.7,'  new xhs =',f12.7)
  177 format(//1x,80('*'))
  180 format(///' transfer final  model sequence from unit 21 to ',
     *  'unit 3'//' models on unit 21 (i, age, radius, luminosity):'/)
  182 format(i4,1p3e15.7,5x,a)
  190 format(//' ***** error. ipartr = 1, but trial model provides',
     *  ' no data')
      end
