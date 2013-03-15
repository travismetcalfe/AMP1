      subroutine srncns
c
c  ************************************************************************
c
c  Sets constants for reaction rate calculation.
c
c  The constants are set in the form defined by the 
c  Fowler, Caughlan and Zimmerman reaction rate expressions.
c  To set these coefficients from a more transparent table
c  of cross-section factors and their energy derivatives,
c  the programme set-fczcoef.f can be used.
c
c  ************************************************************************
c
c  Expanded to include parameters for reactions in CNO cycle.
c  For convenience the same numbering of the first 5 nuclear reactions,
c  and the first 6 Q-values, have
c  been maintained, from the previous version.
c  Also, these parameters are the same as in srncns-3.d.f
c  (i.e. corrected Parker86 parameters). The remaining parameters
c  are, except where otherwise noted, from Gabriels list of
c  July 1990.
c
c  Note on storage assignment: storage is set up in common/rcncns/
c  It is assumed that the amount of storage required for atomic
c  weights is no greater than the number knucst of reactions,
c  and that the number of Q-values required is no greater than
c  knucst + 2.
c
c  Numbering of reactions in arrays storing parameters for reaction
c  rates (i.e. a, b, ss and zz):
c
c   1: H1(H1,e+ nu)H2
c   2: He3(He3,2H1)He4
c   3: He3(He4,gamma)Be7
c   4: Be7(H1,gamma)B8
c   5: N14(H1,gamma)O15
c   6: C12(H1,gamma)N13
c   7: C13(H1,gamma)N14
c   8: N15(H1,He4)C12
c   9: N15(H1,gamma)O16
c  10: O16(H1,gamma)F17
c  11: O17(H1,He4)N14
c
c  Numbering of atomic weights:
c
c  awght(1): He3
c  awght(2): C12
c  awght(3): C13
c  awght(4): N14
c  awght(5): N15
c  awght(6): O16
c  awght(7): O17
c
c  Numbering of Q-values:
c
c   q(1): H1(H1, e+ nu) D2(H1,gam)He3
c   q(2): He3(He3, 2 H1) He4
c   q(3): He3(He4, gam) Be7
c   q(4): Be7(e-, nu) Li7(H1, He4) He4
c   q(5): Be7(H1, gam) B8(, e+ nu) Be8(, He4) He4
c   q(6): Combined Q-value for CN cycle operating in equilibrium.
c   q(7): C12(H1,gamma)N13(,e+ nu)C13
c   q(8): C13(H1,gamma)N14
c   q(9): N14(H1,gamma)O15(,e+ nu)N15
c  q(10): N15(H1,He4)C12
c  q(11): N15(H1,gamma)O16
c  q(12): O16(H1,gamma)F17(,e+ nu)O17
c  q(13): O17(H1,He4)N14
c
c  Note that the coefficients are for computing Na*lambda(i,j),
c  where Na is Avogadro's number and lambda(i,j) is the rate
c  of reactions per pair (i,j)
c
c  Sets of coefficients available (as determined by the parameter ivreng
c  in common/cversn/):
c
c  ivreng = 1: Parameters from Fowler, Caughlan & Zimmerman (1975; 
c              ARAA 13, 69)
c
c  ivreng = 3: Default parameters, essentially from Parker 
c              (Physics of the Sun, Notes from Yale 1987 meeting)
c  ivreng = 4: Reaction constants from J.N. Bahcall routine energy,
c              from May 1992
c  ivreng = 5: Reaction constants from J.N. Bahcall routine energy,
c              from May 1992, with He3 + He4 and Be7 + H1 modified
c              to correspond to Dar and Shaviv
c  ivreng = 6: Original data from Gabriel (1990) (see file gab.coef.in)
c              have been replaced by values from 
c              Bahcall & Pinsonneault (1995; RMP 67, 781) where available.
c  ivreng = 7: Data from Adelberger et al. 
c              (1998; Rev. Mod. Phys. 70, 1265).
c
c
c  Corrected 20/3/89: before, q(4), q(5) and q(6) had been permuted.
c  This was discovered by Yveline Lebreton. In this version, the
c  order has been corrected, and the values have been changed slightly.
c  **** Note: this correction was in fact erroneous. See below. ****
c
c  Modified 20/3/89: Now uses cross-sections for nuclear reactions from
c  Parker (Physics of the Sun, Notes from Yale 1987 meeting) for
c  the pp-reactions.
c
c  Modified 28/8/90:
c  ****************
c  Correct error introduced in Q-values 20/3/89 (note that the
c  "correction introduced there was incorrect).
c
c  *****  Note (30/8/90): it was found that the "correction" applied
c         20/3/89 was in fact
c         incorrect, by checking in detail how the Q-values
c         are used in s/r engenr. In the present version the
c         original order has been restored, but using the Parker
c         parameters.
c
c  Modified 19/2/93: Set up array zz12 containing individual charges
c  for each reaction, to allow for intermediate screening
c
c  Modified 20/2/93: Allow several different sets of parameters,
c  selected by value of ivreng. Default, ivreng = 3, corresponds to
c  parameters used before this date and is set in data statement.
c
c  Modified 22/12/94: Include parameters from Bahcall & Pinsonneault (1994)
c  corresponding to ivreng=6
c
c  Modified 22/7/96: Include parameters from Fowler et al. (1975),
c  corresponding to ivreng=1
c
c  Modified 17/4/98: Include Adelberger et al (1998) parameters
c  corresponding to ivreng=7
c  
c  Modified 06/08/05: Include NACRE reaction rates for ivreng=8
c  
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 11/3/90
c
c..      parameter (knucst = 11, isnuc = 5, nspcmx = 7)
c..      parameter (knucwg = knucst, knucq = knucst + 3, 
c..     *  krnrmx = knucst + 2)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      implicit double precision (a-h,o-z)
c  Note: engenr.n.d.incl replaced by engenr.bz.d.incl, 3/10/05
      include 'engenr.bz.d.incl'
c
      common/rcncns/ a(knucst),b(knucst),s(knucst,isnuc),
     *  zz(knucst),zz12(2,knucst),zsmh,acno,awght(knucwg),q(knucq), 
     *  knucpm, kqvals, kawght
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3,idgmxc,
     *  idgsdt,idgn75,idgn76,idgn77
      common/cstdio/ istdin, istdou, istdpr, istder
c
      common/cversn/ ivrbcs,ivrrhs,ivreos,ivropc,ivreng
c
c  *****************************************************************
c
c  set default reaction rate version to 3, to flag for corrected Parker
c  rates.
c
      data ivreng /3/
c
c  *****************************************************************
c
      save
c
      if(istdpr.gt.0) then
        write(istdpr,'(//'' Setting nuclear parameters''/)')
c
        write(istdpr,*) 'istdpr, istdou, ivreng, idgeng', 
     *    istdpr, istdou, ivreng, idgeng
      end if
c
      call zero(s,knucst*isnuc)
c
c  test for parameter set
c
      if(ivreng.eq.1) then
c
c  Set Fowler, Caughlan and Zimmerman (1975) parameter set
c
        if(istdpr.gt.0) write(istdpr,105)
        a(1)=4.21e-15
        b(1)=3.38
        s(1,1)=0.123
        s(1,2)=1.09
        s(1,3)=0.938
c
        a(2)=5.96e10
        b(2)=12.276
        s(2,1)=0.034
        s(2,2)=-0.199
        s(2,3)=-0.047
        s(2,4)=0.032
        s(2,5)=0.019
c
        a(3)=6.33e6
        b(3)=12.826
        s(3,1)=0.033
        s(3,2)=-0.350
        s(3,3)=-0.080
        s(3,4)=0.056
        s(3,5)=0.033
c
        a(4)=4.35e5
        b(4)=10.262
c
        a(5)=5.08e7
        b(5)=15.228
        s(5,1)=0.027
        s(5,2)=-0.778
        s(5,3)=-0.149
        s(5,4)=0.261
        s(5,5)=0.127
c
c  Z1 =   1.0  Z2 =   6.0 A1 =   1.00782  A2 =  12.00000
c
        a(6)=  2.04E+07
        b(6)=   13.690
        s(6,1)=    0.030
        s(6,2)=    1.19
        s(6,3)=    0.254
        s(6,4)=    2.06
        s(6,5)=    1.12
c
c  Z1 =   1.0  Z2 =   6.0 A1 =   1.00782  A2 =  13.00335
c
        a(7)=  8.01E+07
        b(7)=   13.717
        s(7,1)=    0.030
        s(7,2)=    0.958
        s(7,3)=    0.204
        s(7,4)=    1.39
        s(7,5)=    0.753
c
c  Z1 =   1.0  Z2 =   7.0 A1 =   1.00782  A2 =  15.00011
c
        a(8)=  8.16E+11
        b(8)=   15.251
        s(8,1)=    0.027
        s(8,2)=    6.74
        s(8,3)=    1.29
c
c  Z1 =   1.0  Z2 =   7.0 A1 =   1.00782  A2 =  15.00011
c
        a(9)=  9.78E+08
        b(9)=   15.251
        s(9,1)=    0.027
        s(9,2)=    0.219
        s(9,3)=    0.042
        s(9,4)=    6.83
        s(9,5)=    3.32
c
c  Z1 =   1.0  Z2 =   8.0 A1 =   1.00782  A2 =  15.99492
c
        a(10)=  1.50E+08
        b(10)=   16.692
c
c  O17(H1,He4)N14
c
        a(11)=  1.53d7
        b(11)=   16.712
        s(11,1)=    0.0250
        s(11,2)=    5.39
        s(11,3)=    0.94
        s(11,4)=   13.5
        s(11,5)=    5.98
c
      else if(ivreng.eq.3) then
c
c  Parker (1986) set. Default
c
        if(istdpr.gt.0) write(istdpr,110)
        a(1)=4.006e-15
        b(1)=3.38
        s(1,1)=0.123
        s(1,2)=1.09
        s(1,3)=0.938
c
        a(2)=5.593d10
        b(2)=12.276
        s(2,1)=0.034
        s(2,2)=-0.198
        s(2,3)=-0.047
        s(2,4)=0.032
        s(2,5)=0.019
c
        a(3)=5.610d6
        b(3)=12.826
        s(3,1)=0.033
        s(3,2)=-0.211
        s(3,3)=-0.048
        s(3,4)=0.056
        s(3,5)=0.033
c
        a(4)=3.120d5
        b(4)=10.262
c
c  N14+H1: still from Fowler et al (1975)
c
        a(5)=5.08d7
        b(5)=15.228
        s(5,1)=0.027
        s(5,2)=-0.778
        s(5,3)=-0.149
        s(5,4)=0.261
        s(5,5)=0.127
c
c  Z1 =   1.0  Z2 =   6.0 A1 =   1.00782  A2 =  12.00000
c
        a(6)=  2.11433E+07
        b(6)=   13.6923
        s(6,1)=    0.0304
        s(6,2)=    0.6646
        s(6,3)=    0.1416
        s(6,4)=    0.3627
        s(6,5)=    0.1965
c
c  Z1 =   1.0  Z2 =   6.0 A1 =   1.00782  A2 =  13.00335
c
        a(7)=  8.00387E+07
        b(7)=   13.7197
        s(7,1)=    0.0304
        s(7,2)=    0.9602
        s(7,3)=    0.2041
        s(7,4)=    1.3936
        s(7,5)=    0.7533
c
c  Z1 =   1.0  Z2 =   7.0 A1 =   1.00782  A2 =  15.00011
c
        a(8)=  1.19112E+12
        b(8)=   15.2535
        s(8,1)=    0.0273
        s(8,2)=    1.9717
        s(8,3)=    0.3770
        s(8,4)=   13.6599
        s(8,5)=    6.6418
c
c  Z1 =   1.0  Z2 =   7.0 A1 =   1.00782  A2 =  15.00011
c
        a(9)=  9.77327E+08
        b(9)=   15.2535
        s(9,1)=    0.0273
        s(9,2)=    0.2054
        s(9,3)=    0.0393
        s(9,4)=    5.9993
        s(9,5)=    2.9170
c
c  Z1 =   1.0  Z2 =   8.0 A1 =   1.00782  A2 =  15.99492
c
        a(10)=  1.49882E+08
        b(10)=   16.6955
        s(10,1)=    0.0250
        s(10,2)=   -1.1734
        s(10,3)=   -0.2050
        s(10,4)=    0.7340
        s(10,5)=    0.3261
c
c  O17(H1,He4)N14
c  Coefficients taken from Caughlan & Fowler (1988), 
c  excluding extra terms that do not conform with standard form
c
        a(11)=  1.53d7
        b(11)=   16.712
        s(11,1)=    0.0250
        s(11,2)=    5.39
        s(11,3)=    0.94
        s(11,4)=   13.5
        s(11,5)=    5.98
c
      else if(ivreng.eq.4) then
c
c  Reaction constants from J.N. Bahcall routine energy, from May 1992
c
c  Reaction: H1 + H1
c
	if(istdpr.gt.0) write(istdpr,120)
c
        a(1) =  3.93678E-15
        b(1) =  3.37960E+00
        s(1,1) =  1.23100E-01
        s(1,2) =  1.07900E+00
        s(1,3) =  9.30400E-01
        s(1,4) =  0.00000E+00
        s(1,5) =  0.00000E+00
c
c  Reaction: He3 + He3
c
        a(2) =  5.42084E+10
        b(2) =  1.22720E+01
        s(2,1) =  3.39000E-02
        s(2,2) = -6.16000E-02
        s(2,3) = -1.46400E-02
        s(2,4) =  0.00000E+00
        s(2,5) =  0.00000E+00
c
c  Reaction: He3 + He4
c
        a(3) =  5.53004E+06
        b(3) =  1.28220E+01
        s(3,1) =  3.26000E-02
        s(3,2) = -2.11500E-01
        s(3,3) = -4.80900E-02
        s(3,4) =  0.00000E+00
        s(3,5) =  0.00000E+00
c
c  Reaction: Be7 + H
c
        a(4) =  2.90463E+05
        b(4) =  1.02620E+01
        s(4,1) =  4.05650E-02
        s(4,2) = -3.64000E-01
        s(4,3) = -1.03400E-01
        s(4,4) =  0.00000E+00
        s(4,5) =  0.00000E+00
c
c  Reaction: N14 + H
c
        a(5) =  5.07728E+07
        b(5) =  1.52247E+01
        s(5,1) =  2.74000E-02
        s(5,2) = -7.78800E-01
        s(5,3) = -1.49100E-01
        s(5,4) =  2.61200E-01
        s(5,5) =  1.27170E-01
c
c  Reaction: C12 + H
c
        a(6) =  2.11426E+07
        b(6) =  1.36869E+01
        s(6,1) =  3.04000E-02
        s(6,2) =  6.64500E-01
        s(6,3) =  1.41550E-01
        s(6,4) =  3.62700E+00
        s(6,5) =  1.96470E+00
c
c  Reaction: C13 + H
c
        a(7) =  8.00371E+07
        b(7) =  1.37142E+01
        s(7,1) =  3.03000E-02
        s(7,2) =  9.60100E-01
        s(7,3) =  2.04100E-01
        s(7,4) =  1.39400E+00
        s(7,5) =  7.53000E-01
c
c  Z1 =   1.0  Z2 =   7.0 A1 =   1.00782  A2 =  15.00011
c  (same as for ivreng = 3)
c
        a(8)=  1.19112E+12
        b(8)=   15.2535
        s(8,1)=    0.0273
        s(8,2)=    1.9717
        s(8,3)=    0.3770
        s(8,4)=   13.6599
        s(8,5)=    6.6418
c
c  Z1 =   1.0  Z2 =   7.0 A1 =   1.00782  A2 =  15.00011
c  (same as for ivreng = 3)
c
        a(9)=  9.77327E+08
        b(9)=   15.2535
        s(9,1)=    0.0273
        s(9,2)=    0.2054
        s(9,3)=    0.0393
        s(9,4)=    5.9993
        s(9,5)=    2.9170
c
c  Reaction: O16 + H
c
        a(10) =  1.49883E+08
        b(10) =  1.66888E+01
        s(10,1) =  2.49400E-02
        s(10,2) = -1.17300E+00
        s(10,3) = -2.05000E-01
        s(10,4) =  7.34000E-01
        s(10,5) =  3.26100E-01
c
c  O17(H1,He4)N14
c  Coefficients taken from Caughlan & Fowler (1988), 
c  excluding extra terms that do not conform with standard form
c  (same as for ivreng = 3)
c
        a(11)=  1.53d7
        b(11)=   16.712
        s(11,1)=    0.0250
        s(11,2)=    5.39
        s(11,3)=    0.94
        s(11,4)=   13.5
        s(11,5)=    5.98
c
      else if(ivreng.eq.5) then
c
c  Reaction constants from J.N. Bahcall routine energy, from May 1992
c  with He3 + He4 and Be7 + H1 modified to correspond to Dar and Shaviv
c
c  Reaction: H1 + H1
c
	if(istdpr.gt.0) write(istdpr,125)
c
        a(1) =  3.93678E-15
        b(1) =  3.37960E+00
        s(1,1) =  1.23100E-01
        s(1,2) =  1.07900E+00
        s(1,3) =  9.30400E-01
        s(1,4) =  0.00000E+00
        s(1,5) =  0.00000E+00
c
c  Reaction: He3 + He3
c
        a(2) =  5.42084E+10
        b(2) =  1.22720E+01
        s(2,1) =  3.39000E-02
        s(2,2) = -6.16000E-02
        s(2,3) = -1.46400E-02
        s(2,4) =  0.00000E+00
        s(2,5) =  0.00000E+00
c
c  Reaction: He3 + He4
c
        a(3) =  4.66947E+06
        b(3) =  1.28220E+01
        s(3,1) =  3.26000E-02
        s(3,2) = -2.11500E-01
        s(3,3) = -4.80900E-02
        s(3,4) =  0.00000E+00
        s(3,5) =  0.00000E+00
c
c  Reaction: Be7 + H
c
        a(4) =  2.20451E+05
        b(4) =  1.02620E+01
        s(4,1) =  4.05650E-02
        s(4,2) = -3.64000E-01
        s(4,3) = -1.03400E-01
        s(4,4) =  0.00000E+00
        s(4,5) =  0.00000E+00
c
c  Reaction: N14 + H
c
        a(5) =  5.07728E+07
        b(5) =  1.52247E+01
        s(5,1) =  2.74000E-02
        s(5,2) = -7.78800E-01
        s(5,3) = -1.49100E-01
        s(5,4) =  2.61200E-01
        s(5,5) =  1.27170E-01
c
c  Reaction: C12 + H
c
        a(6) =  2.11426E+07
        b(6) =  1.36869E+01
        s(6,1) =  3.04000E-02
        s(6,2) =  6.64500E-01
        s(6,3) =  1.41550E-01
        s(6,4) =  3.62700E+00
        s(6,5) =  1.96470E+00
c
c  Reaction: C13 + H
c
        a(7) =  8.00371E+07
        b(7) =  1.37142E+01
        s(7,1) =  3.03000E-02
        s(7,2) =  9.60100E-01
        s(7,3) =  2.04100E-01
        s(7,4) =  1.39400E+00
        s(7,5) =  7.53000E-01
c
c  Z1 =   1.0  Z2 =   7.0 A1 =   1.00782  A2 =  15.00011
c  (same as for ivreng = 3)
c
        a(8)=  1.19112E+12
        b(8)=   15.2535
        s(8,1)=    0.0273
        s(8,2)=    1.9717
        s(8,3)=    0.3770
        s(8,4)=   13.6599
        s(8,5)=    6.6418
c
c  Z1 =   1.0  Z2 =   7.0 A1 =   1.00782  A2 =  15.00011
c  (same as for ivreng = 3)
c
        a(9)=  9.77327E+08
        b(9)=   15.2535
        s(9,1)=    0.0273
        s(9,2)=    0.2054
        s(9,3)=    0.0393
        s(9,4)=    5.9993
        s(9,5)=    2.9170
c
c  Reaction: O16 + H
c
        a(10) =  1.49883E+08
        b(10) =  1.66888E+01
        s(10,1) =  2.49400E-02
        s(10,2) = -1.17300E+00
        s(10,3) = -2.05000E-01
        s(10,4) =  7.34000E-01
        s(10,5) =  3.26100E-01
c
c  O17(H1,He4)N14
c  Coefficients taken from Caughlan & Fowler (1988), 
c  excluding extra terms that do not conform with standard form
c  (same as for ivreng = 3)
c
        a(11)=  1.53d7
        b(11)=   16.712
        s(11,1)=    0.0250
        s(11,2)=    5.39
        s(11,3)=    0.94
        s(11,4)=   13.5
        s(11,5)=    5.98
c
      else if(ivreng.eq.6) then
c
	if(istdpr.gt.0) write(istdpr,130)
c
c Reaction parameters. Original data from Gabriel (1990)
c (see file gab.coef.in) have been replaced by values from
c Bahcall & Pinsonneault (1994) where available.
c
c Z1 =   1.0  Z2 =   1.0 A1 =   1.00782  A2 =   1.00782
c
        n=1
        a(n)=  3.93686E-15
        b(n)=    3.3810
        s(n,1)=    0.1232
        s(n,2)=    1.0974
        s(n,3)=    0.9467
        s(n,4)=    0.0000
        s(n,5)=    0.0000
c
c Z1 =   2.0  Z2 =   2.0 A1 =   3.01603  A2 =   3.01603
c
        n=2
        a(n)=  5.42080E+10
        b(n)=   12.2772
        s(n,1)=    0.0339
        s(n,2)=   -0.2186
        s(n,3)=   -0.0519
        s(n,4)=    0.0348
        s(n,5)=    0.0210
c
c Z1 =   2.0  Z2 =   2.0 A1 =   3.01603  A2 =   4.00260
c
        n=3
        a(n)=  5.53071E+06
        b(n)=   12.8274
        s(n,1)=    0.0325
        s(n,2)=   -0.2143
        s(n,3)=   -0.0487
        s(n,4)=    0.0000
        s(n,5)=    0.0000
c
c Z1 =   1.0  Z2 =   4.0 A1 =   1.00782  A2 =   7.01693
c
        n=4
        a(n)=  2.90476E+05
        b(n)=   10.2643
        s(n,1)=    0.0406
        s(n,2)=   -0.3949
        s(n,3)=   -0.1122
        s(n,4)=    0.0000
        s(n,5)=    0.0000
c
c Z1 =   1.0  Z2 =   7.0 A1 =   1.00782  A2 =  14.00307
c
        n=5
        a(n)=  5.07745E+07
        b(n)=   15.2308
        s(n,1)=    0.0274
        s(n,2)=   -0.7788
        s(n,3)=   -0.1491
        s(n,4)=    0.2612
        s(n,5)=    0.1272
c
c Z1 =   1.0  Z2 =   6.0 A1 =   1.00782  A2 =  12.00000
c
        n=6
        a(n)=  2.11433E+07
        b(n)=   13.6923
        s(n,1)=    0.0304
        s(n,2)=    0.6646
        s(n,3)=    0.1416
        s(n,4)=    0.3627
        s(n,5)=    0.1965
c
c Z1 =   1.0  Z2 =   6.0 A1 =   1.00782  A2 =  13.00335
c
        n=7
        a(n)=  8.00387E+07
        b(n)=   13.7197
        s(n,1)=    0.0304
        s(n,2)=    0.9602
        s(n,3)=    0.2041
        s(n,4)=    1.3936
        s(n,5)=    0.7533
c
c Z1 =   1.0  Z2 =   7.0 A1 =   1.00782  A2 =  15.00011
c
        n=8
        a(n)=  1.19112E+12
        b(n)=   15.2535
        s(n,1)=    0.0273
        s(n,2)=    1.9717
        s(n,3)=    0.3770
        s(n,4)=   13.6599
        s(n,5)=    6.6418
c
c Z1 =   1.0  Z2 =   7.0 A1 =   1.00782  A2 =  15.00011
c
        n=9
        a(n)=  9.77327E+08
        b(n)=   15.2535
        s(n,1)=    0.0273
        s(n,2)=    0.2054
        s(n,3)=    0.0393
        s(n,4)=    5.9993
        s(n,5)=    2.9170
c
c Z1 =   1.0  Z2 =   8.0 A1 =   1.00782  A2 =  15.99492
c
        n=10
        a(n)=  1.49882E+08
        b(n)=   16.6955
        s(n,1)=    0.0250
        s(n,2)=   -1.1734
        s(n,3)=   -0.2050
        s(n,4)=    0.7340
        s(n,5)=    0.3261
c
c  O17(H1,He4)N14
c  Coefficients taken from Caughlan & Fowler (1988), 
c  excluding extra terms that do not conform with standard form
c  (same as for ivreng = 3)
c
        a(11)=  1.53d7
        b(11)=   16.712
        s(11,1)=    0.0250
        s(11,2)=    5.39
        s(11,3)=    0.94
        s(11,4)=   13.5
        s(11,5)=    5.98
c
      else if(ivreng.eq.7) then
c
	if(istdpr.gt.0) write(istdpr,140)
c
c  Results from Adelberger et al (1998)
c
c Z1 =   1.0  Z2 =   1.0 A1 =   1.00782  A2 =   1.00782
c
        n=1
        a(n)=  3.93686E-15
        b(n)=    3.3810
        s(n,1)=    0.1232
        s(n,2)=    1.0877
        s(n,3)=    0.9383
        s(n,4)=    0.0000
        s(n,5)=    0.0000
c
c Z1 =   2.0  Z2 =   2.0 A1 =   3.01603  A2 =   3.01603
c
        n=2
        a(n)=  5.85447E+10
        b(n)=   12.2772
        s(n,1)=    0.0339
        s(n,2)=   -0.2678
        s(n,3)=   -0.0636
        s(n,4)=    0.0530
        s(n,5)=    0.0320
c
c Z1 =   2.0  Z2 =   2.0 A1 =   3.01603  A2 =   4.00260
c
        n=3
        a(n)=  5.49958E+06
        b(n)=   12.8274
        s(n,1)=    0.0325
        s(n,2)=   -0.2086
        s(n,3)=   -0.0474
        s(n,4)=    0.0000
        s(n,5)=    0.0000
c
c Z1 =   1.0  Z2 =   4.0 A1 =   1.00782  A2 =   7.01693
c
        n=4
        a(n)=  2.46386E+05
        b(n)=   10.2643
        s(n,1)=    0.0406
        s(n,2)=   -0.2095
        s(n,3)=   -0.0595
        s(n,4)=    0.1677
        s(n,5)=    0.1212
c
c Z1 =   1.0  Z2 =   7.0 A1 =   1.00782  A2 =  14.00307
c
        n=5
        a(n)=  5.35273E+07
        b(n)=   15.2308
        s(n,1)=    0.0274
        s(n,2)=   -1.6000
        s(n,3)=   -0.3064
        s(n,4)=    0.0000
        s(n,5)=    0.0000
c
c Z1 =   1.0  Z2 =   6.0 A1 =   1.00782  A2 =  12.00000
c
        n=6
        a(n)=  1.95393E+07
        b(n)=   13.6923
        s(n,1)=    0.0304
        s(n,2)=    0.7631
        s(n,3)=    0.1626
        s(n,4)=    4.7908
        s(n,5)=    2.5950
c
c Z1 =   1.0  Z2 =   6.0 A1 =   1.00782  A2 =  13.00335
c
        n=7
        a(n)=  1.10599E+08
        b(n)=   13.7197
        s(n,1)=    0.0304
        s(n,2)=   -0.4045
        s(n,3)=   -0.0860
        s(n,4)=    7.4590
        s(n,5)=    4.0322
c
c Z1 =   1.0  Z2 =   7.0 A1 =   1.00782  A2 =  15.00011
c
        n=8
        a(n)=  1.03077E+12
        b(n)=   15.2535
        s(n,1)=    0.0273
        s(n,2)=    2.0123
        s(n,3)=    0.3848
        s(n,4)=   17.0646
        s(n,5)=    8.2973
c
c Z1 =   1.0  Z2 =   7.0 A1 =   1.00782  A2 =  15.00011
c
        n=9
        a(n)=  9.77327E+08
        b(n)=   15.2535
        s(n,1)=    0.0273
        s(n,2)=    0.1438
        s(n,3)=    0.0275
        s(n,4)=    6.1492
        s(n,5)=    2.9899
c
c Z1 =   1.0  Z2 =   8.0 A1 =   1.00782  A2 =  15.99492
c
        n=10
        a(n)=  1.49882E+08
        b(n)=   16.6955
        s(n,1)=    0.0250
        s(n,2)=   -1.2244
        s(n,3)=   -0.2139
        s(n,4)=    0.6973
        s(n,5)=    0.3098
c
c  O17(H1,He4)N14
c  Coefficients taken from Caughlan & Fowler (1988), 
c  excluding extra terms that do not conform with standard form
c  (same as for ivreng = 3)
c
        a(11)=  1.53d7
        b(11)=   16.712
        s(11,1)=    0.0250
        s(11,2)=    5.39
        s(11,3)=    0.94
        s(11,4)=   13.5
        s(11,5)=    5.98
      else if(ivreng.eq.8) then
c
c NACRE reaction rates
c
c Test
c
        if(idgeng.gt.0.and.istdpr.gt.0)
     *    write(istdpr,*) "OK for NACRE reaction rates setting"
c
c End test
c
        if(istdpr.gt.0) write(istdpr,150)
        call snrnacre
        if(idgeng.ge.2.and.istdpr.gt.0)
     *    write(istdpr,155) s(1,3),a(4),a(10),b(10)
c
      else
	write(istdou,190) ivreng
	if(istdou.ne.istdpr.and.istdpr.gt.0) write(istdpr,190) ivreng
	stop
      end if
c
c  charges and their products
c
      zz(1)=1
      zz12(1,1)=1
      zz12(2,1)=1
c
      zz(2)=4
      zz12(1,2)=2
      zz12(2,2)=2
c
      zz(3)=4
      zz12(1,3)=2
      zz12(2,3)=2
c
      zz(4)=4
      zz12(1,4)=1
      zz12(2,4)=4
c
      zz(5)=7
      zz12(1,5)=1
      zz12(2,5)=7
c
      zz(6)= 6
      zz12(1,6)=1
      zz12(2,6)=6
c
      zz(7)= 6
      zz12(1,7)=1
      zz12(2,7)=6
c
      zz(8)= 7
      zz12(1,8)=1
      zz12(2,8)=7
c
      zz(9)= 7
      zz12(1,9)=1
      zz12(2,9)=7
c
      zz(10)= 8
      zz12(1,10)=1
      zz12(2,10)=8
c
      zz(11)= 8
      zz12(1,11)=1
      zz12(2,11)=8
      zsmh=5.02
c
c  average atomic weight for CNO elements
c
      acno=14.92
c
      awght(1) = 3.016
      awght(2) = 12.0
      awght(3) = 13.00335
      awght(4) = 14.00307
      awght(5) = 15.00011
      awght(6) = 15.99492
      awght(7) = 16.99913
c
c  Q-values, in MeV, corrected for nu-energy
c
      q(1)=6.686
      q(2)=12.859
      q(3)=1.587
      q(4)=17.397
c
c  PPIII effective Q from Bahcall & Ulrich 1988
c
      q(5)=11.499
      q(6)=24.97
c
c  partial CNO Q-values, largely from Gabriel (1990)
c
      q(7)=4.0943
      q(8)=7.551
      q(9)=9.057
      q(10)=4.966
      q(11)=12.128
c
c  q(12) and q(13) are set from tables of mass defects in Clayton,
c  using neutrino energy from Gabriel. The sum of the two values
c  is 0.002 MeV below value for combined reaction given by Gabriel.
c
      q(12)=2.366
      q(13)=1.193
c
c  set numbers of parameters that have been set
c
      knucpm = 11
      kqvals = 13
      kawght = 7
c
      return
  105 format(/
     *   ' **** setting Fowler et al. (1975) reaction parameters'/)
  110 format(/
     *   ' **** setting Parker 1987 reaction parameters'/)
  120 format(/
     *   ' **** setting Bahcall 1992 reaction parameters'/)
  125 format(/
     *   ' **** setting Bahcall 1992 reaction parameters'/
     *   '      with Dar and Shaviv He3 + He4 and Be7 + H1'/)
  130 format(/
     *   ' **** setting Bahcall & Pinsonneault 1994',
     *   ' reaction parameters')
  140 format(/
     *   ' **** setting Adelberger et al. 1998',
     *   ' reaction parameters')
  150 format(/
     *   ' **** setting NACRE reaction parameters')
  155   format('Some coefficients, s(1,3)',e13.5,' a(4) :',e13.5,
     &  ' a(10) :',e13.5,' b(10) :',e13.5,' ...Seems to work')
  190 format(//' ***** Error in s/r srncns. ivreng = ',i4,
     *  ' not implemented')
      end
      block data blengr
c
c  initialize controls for reaction rate and energy generation
c  routines
c
      implicit double precision (a-h,o-z)
      common/he3eql/ xhe3eq,ieqhe3,iequsd
      common/engcnt/ xhzlm1, xhzlm2, nxhzer, icnocs, nspec,
     *  iheccs, nspect
      common/cnofrc/ fcno, xtlcno
      common/degfct/ thte,iscren
c
      data ieqhe3 /0/
      data fcno /0.28/
      data thte,iscren /1.,1/
      data icnocs /0/
c
      end
