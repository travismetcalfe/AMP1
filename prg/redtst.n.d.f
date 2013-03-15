      subroutine redtst(y,zk,yp,zkp,dt,redfac,eam,ii,kk,nn,id)
c
c  Reduces timestep dt by factor redfac, and sets new trial
c  solution into yp. If eam is small (with hard-coded condition)
c  interpolates between previous timestep and unconverged solution.
c  Otherwise use solution at previous time step as trial.
c
c  Original version: 11/12/92
c
      implicit double precision (a-h,o-z)
      dimension y(id,*),zk(*),yp(id,*),zkp(*)
      common/noiter/ iter, ntime, eps, eamcon, it_force
c
c  common defining standard input and output, standard error
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      save
c
c  set new timestep
c
      dt=redfac*dt
c
c  test for interpolation or not
c
      if(ntime.eq.-244) then
        eamlim=1.d-5
      else
c
c  reset criterion to force use of previous solution.
c
        write(istder,
     *    '(/'' ***** Force previous solution in s/r redtst'')')
        eamlim=1.d-50
      end if
c
      if(eam.gt.eamlim) then
        do 15 n=1,nn
        do 15 i=1,ii
   15   yp(i,n)= y(i,n)
c
      else
	fct1=1-redfac
	fct2=redfac
        do 20 n=1,nn
        do 20 i=1,ii
   20   yp(i,n)= fct1*y(i,n)+fct2*yp(i,n)
c
      end if
c
      if(kk.eq.0) return
      do 30 k=1,kk
   30 zkp(k)=fct1*zk(k)+fct2*zkp(k)
      return
      end
