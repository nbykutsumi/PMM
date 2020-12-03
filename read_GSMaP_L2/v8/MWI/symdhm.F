cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         calculate sequnence no( iseq)
c           for iyear, month, iday, ihour and minute)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine symdhm( iyear, month, iday, ihour, minute, iseq)
      minute = mod( iseq, 60 )
      call stymdh( iyear, month, iday, ihour, iseq/ 60)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         determine year, month, day, hour from record sequence no
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine stymdh( iy, im, id, ih, irecsq)
      is =irecsq / 24 + 1
      ih = mod( irecsq, 24)
      call stoymd( is, iy, im, id)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculating date of gregorian calendar               c
c     from sequence day number since 1801.01.01            c
c     xx  ;  sequence number  (integer type)               c
c     y   ;  year (ex.1983)   (integer type)               c
c     m   ;  month            (integer type)               c
c     d   ;  day              (integer type)               c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine stoymd(xx,y,m,d)
      implicit integer(a-z)
      x=xx
      y=int((x+306-0.1)/365.25)+1802
 1020 k=int(365.25*(y-1801)-306)
      k=k-((((y-1501)/100)*3/4)-2)
      if(x-k)2020,2020,2021
 2020 y=y-1
      goto 1020
 2021 x=x-k
      m=int((x-0.5)/30.6)+3
      k=int(30.6*(m-3)+0.5)
      d=x-k
      if(m-12)2040,2040,2041
 2041 m=m-12
      y=y+1
 2040 y=y-1
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         calculate sequence no ( iseq)
c           for iyear, month, iday, ihour, minute
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ymdhms( iyear, month, iday, ihour, minute, iseq)
      call ymdhts( iyear, month, iday, ihour, iseq)
      iseq = iseq * 60 + minute
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculate record sequence no for year, month, day, hour
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ymdhts( iy, im, id, ih, irecsq)
      call ymdtos( is, iy, im, id)
      irecsq = ( is - 1) * 24 + ih
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculating sequence day number since 1801.01.01     c
c     from date of gregorian calendar                      c
c     x   ;   sequence number  (integer type)              c
c     yy  ;   year (ex.1983)   (integer type)              c
c     mm  ;   month            (integer type)              c
c     dd  ;   day              (integer type)              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ymdtos(x,yy,mm,dd)
      implicit integer(a-z)
      y=yy
      m=mm
      d=dd
      if(m-3)2020,2021,2021
 2020 y=y-1
      m=m+12
 2021 y=y+1
      x=int(365.25*(y-1801)-306)
      x=x-((((y-1501)/100)*3/4)-2)
      x=x+int(30.6*(m-3)+0.5)+d
      return
      end
