      subroutine setcvr(n)
c
c  stores CNO composition variables from values in common/compos/
c  into cvr(.,n) in common/compvr/, if values are non-negative
c
c  Original version: 10/8/02
c
      implicit double precision (a-h,o-z)
c
c  include statement setting storage indices for reaction rate
c  and energy generation parameters
c
      include 'engenr.bz.d.incl'
c
      common/compos/ xseth, yset, xset3, xset12, xset13, xset14, xset16
      common/compvr/ cvr(icvrmx, 1)
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
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng,idghe3,idgmxc
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr, istder
c
      if(icnocs.eq.1) then
        if(xset14.ge.0) cvr(icvn14,n)=xset14
        if(xset16.ge.0) cvr(icvo16,n)=xset16
      else if(icnocs.eq.2.or.icnocs.eq.4) then
        if(xset12.ge.0) cvr(icvc12,n)=xset12
        if(xset13.ge.0) cvr(icvc13,n)=xset13
        if(xset14.ge.0) cvr(icvn14,n)=xset14
        if(xset16.ge.0) cvr(icvo16,n)=xset16
      else if(icnocs.eq.3) then
        if(xset12.ge.0) cvr(icvc12,n)=xset12
        if(xset13.ge.0) cvr(icvc13,n)=xset13
        if(xset14.ge.0) cvr(icvn14,n)=xset14
      end if
c
      return
      end
