      subroutine sthvab(x,z,xnum,aharg,amarg,nspec,ngarb,zab)
c
c  set heavy element abundances, assuming fixed X, Z and 
c  number fractions except for garbage.
c  xnum must contain number fractions on input, am the atomic masses,
c  and ngarb the index of the garbage element.
c  returns mass fractions, relative to Z, in zab
c
c  Original version: 25/9/91
c
      implicit double precision (a-h, o-z)
      dimension xnum(1), amarg(1), zab(1)
      dimension am(10)
c
      ah=aharg
      do 10 n=1,nspec
   10 am(n)=amarg(n)
c
c  reset relevant atomic weights to WD values
c
      ah    =1.008
      am(1) =12.01150
      am(2) =14.00700
      am(3) =15.99900
      am(10)=55.84700
c
      sum = 0
      do 15 n=1,nspec
      if(n.ne.ngarb) then
	zab(n) = xnum(n)*am(n)*x/(ah*z)
	sum = sum+zab(n)
      end if
   15 continue
c
      zab(ngarb) = 1 - sum
      return
      end
