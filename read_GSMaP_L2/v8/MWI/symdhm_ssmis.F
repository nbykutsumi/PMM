      subroutine symdhm( iyear, month, iday, ihour, minute, iseq)
c--------- Rev  2012.06.12 tashima@restec -s ---------
c      integer*2 iyear, month, iday,idayall, ihour, minute
c      integer*4 iseq,iseq2,iseq3,iseq4
      integer iyear, month, iday,idayall, ihour, minute
      integer iseq,iseq2,iseq3,iseq4
c--------- Rev  2012.06.12 tashima@restec -e ---------
      integer*2 mend(12),mend2(12)
      data mend/31,59,90,120,151,181,212,243,273,304,334,365/
      data mend2/31,60,91,121,152,182,213,244,274,305,335,366/

      minute = mod( iseq, 60 )
      iseq2 = (iseq -minute)/60
      ihour = mod(iseq2,24)
      iseq3 = (iseq2-ihour)/24
      idayall = mod(iseq3,366)
      iseq4 = (iseq3-idayall)/366
      iyear = 1900+iseq4

      if ( ((mod(iyear,4).eq.1).and.(idayall.eq.0)).or.
     &      (mod(iyear,4).eq.0) ) then

         if (idayall.eq.0) then
            iyear=1900+iseq4-1
            month=12
            iday=31
         else if ( (idayall.gt.0).and.(idayall.le.mend2(1)) ) then
            iyear=1900+iseq4
            month=1
            iday=idayall
         else if (idayall.gt.mend2(1)) then
            do m=2,12               
               do md=mend2(m-1),mend2(m)
                  if (md.eq.idayall) then
                     month=m
                     iday=idayall-mend2(m-1)
                     GO TO 100
                  endif
               enddo
            enddo
 100        iyear=1900+iseq4
         endif

      else
         
         if(idayall.le.mend(1))then
            iyear=1900+iseq4
            month=1
            iday=idayall
         else
            do m=2,12
               do md=mend(m-1),mend(m)
                  if (md.eq.idayall) then
                     month=m
                     iday=idayall-mend(m-1)
                      GO TO 200
                   end if
                enddo
             enddo
 200         iyear=1900+iseq4
          endif
       endif

      return
      end
