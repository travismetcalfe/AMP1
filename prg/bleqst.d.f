      block data bleqst
c  initialize data for equation of state
c
c  Modified 22/5/90 to include all levels of ionization of Ar and Fe.
c  Previously only the first 15 levels were included, the remainder
c  being forced to be unionized
c
      implicit double precision (a-h,o-z)
      character name*5
      common/eqscnt/ anh0,anhe0,ihvz,iprrad,ihmin,igndeg
      common/hvabnd/ ab(10),iab
      common/eqdpco/ frhi,bdcoh,bdcoz,idpco
      common/eqphcs/ ic,c(48)
      common/potetc/ chi(125),am(10),iz(10)
      common/hvname/ name(10)
      common /hvcntl/ icnthv,iwrthv,dptst0,dptst1
      common/hvomeg/ iom(26),iom1(20)
      common/hvomcl/ iomfll
      common/ccoulm/ epsdmu, icoulm, iclmit, iclmsd, nitdmu, epssdr
c
      data anh0,anhe0,ihvz,iprrad,ihmin /0.5,6.,1,1,0/
      data ab,iab /0.2254,0.0549,0.4987,0.0335,0.00197,0.0436,
     *  0.00403,0.0565,0.00180,0.0795,10/
      data frhi,bdcoh,bdcoz,idpco /1.,1.,1.,0/
c  coefficients for s/r phder.
      data ic,c/4,2.315472,7.128660,7.504998,2.665350,7.837752,23.507934
     .,23.311317,7.987465,9.215560,26.834068,25.082745,8.020509,3.693280
     .,10.333176,9.168960,2.668248,2.315472,6.748104,6.564912,2.132280,
     .7.837752,21.439740,19.080088,5.478100,9.215560,23.551504,19.015888
     .,4.679944,3.693280,8.859868,6.500712,1.334124,1.157736,3.770676,
     .4.015224,1.402284,8.283420,26.184486,28.211372,10.310306,14.755480
     .,45.031658,46.909420,16.633242,7.386560,22.159680,22.438048,
     .7.664928/
      data iz/6,7,8,10,11,12,13,14,18,26/
      data name/'   C','   N','   O','  Ne','  Na',
     .  '  Mg','  Al','  Si','  Ar','  Fe'/
      data am/12.00,14.01,16.00,20.17,22.99,24.31,26.98,28.08,
     .  39.94,55.84/
      data chi/11.26,24.38,47.86,64.48,391.99,489.84,
     .  14.54,29.60,47.43,77.45,97.86,551.92,666.83,
     .  13.61,35.15,54.93,77.39,113.87,138.08,739.11,871.12,
     .  21.56,41.07,63.5,97.16,126.4,157.91,207.3,239.,1196.,1360.,
     .  5.14,47.29,71.65,98.88,138.60,172.36,208.44,264.15,299.78,
     .  1465.,1646.,
     .  7.64,15.03,80.12,109.29,141.23,186.86,225.31,265.96,327.90,
     .  367.36,1761.2,2085.,
     .  5.98,18.82,28.44,119.96,153.77,190.42,241.93,285.13,330.1,
     .  398.5,441.9,2085.5,2299.,
     .  8.15,16.34,33.46,45.13,166.73,205.11,246.41,303.87,
     .  351.83,401.3,476.0,523.2,2436.,2666.,
     .  15.75,27.62,40.90,59.79,75.0,91.3,124.,143.46,421.,480.,
     .  539.5,621.1,688.5,755.5,854.4,918.,4121.,4426.,
     .  7.90,16.18,30.64,56.,79.,105.,133.,151.,235.,262.,290.,321.,
     .  355.,390.,457.,489.,1266.,1358.,1456.,1582.,1689.,1799.,1950.,
     .  2045.,8828.,9278./
      data icnthv,iwrthv /0,0/
c
c  limiting arguments in exponential in test for full or no
c  ionization in s/r hviona.
c  dptst0 is used for first level and element to ionize
c  dptst1 is used for remaining levels
c
      data dptst0,dptst1 /85.,19./
c
c  quantities for calculating statistical weights in
c  function omega.
c
      data iom/2,1,2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,10,21,28,25,6,25,
     .  28,21/
      data iom1/2,1,1,2,10,15,21,28,28,25,7,6,6,7,25,30,28,21,21,10/
c  flag for statistical weight omega for fully ionized atom.
c  before 5/1/84 omega was always set to 15, for some bizarre reason.
c  for transition introduce flag iomfll so that iomfll = 0 corresponds
c  to old situation, and iomfll .ne. 0 gives correct value omega = 1.
c  Change default to 1, 30/9/91
      data iomfll /1/
c
c  controls for Coulomb effect.
c
      data epsdmu, icoulm, iclmit, iclmsd
     *    / 1.e-12,   0,       1,      0 /
      end
