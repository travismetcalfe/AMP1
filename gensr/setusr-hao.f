      subroutine setusr(user)
c
c  returns user name of user running the programme.
c  Note: this is currently hardcoded to return jcd
c  for simple use at HAO
c
c  Original version: 3/12/92
c
      character*(*) user
      integer getuid
c
      user = 'jcd'
c
c..      id=getuid()
c..c
c..      if(id.eq.534) then
c..        user='aake'
c..      else if(id.eq.535) then
c..        user='brandenb'
c..      else if(id.eq.503) then
c..        user='bt'
c..      else if(id.eq.519) then
c..        user='fgj'
c..      else if(id.eq.509) then
c..        user='hans'
c..      else if(id.eq.501) then
c..        user='jcd'
c..      else if(id.eq.530) then
c..        user='jens'
c..      else if(id.eq.515) then
c..        user='jesm'
c..      else if(id.eq.554) then
c..        user='jones'
c..      else if(id.eq.542) then
c..        user='kj'
c..      else if(id.eq.544) then
c..        user='michel'
c..      else if(id.eq.505) then
c..        user='mjt'
c..      else if(id.eq.514) then
c..        user='mlo'
c..      else if(id.eq.537) then
c..        user='ms'
c..      else if(id.eq.547) then
c..        user='mv'
c..      else if(id.eq.541) then
c..        user='nha'
c..      else if(id.eq.556) then
c..        user='nuspl'
c..      else if(id.eq.513) then
c..        user='pen'
c..      else if(id.eq.512) then
c..        user='pg'
c..      else if(id.eq.543) then
c..        user='pot'
c..      else if(id.eq.510) then
c..        user='schou'
c..      else if(id.eq.500) then
c..        user='srf'
c..      else if(id.eq.529) then
c..        user='tang'
c..      else if(id.eq.536) then
c..        user='th'
c..      else
c..	write(6,100) id
c..	stop
c..      endif
      return
  100 format(//' ***** Error in setusr. ID = ',i5,' not in table')
      end
