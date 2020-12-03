      IMPLICIT  NONE

      integer NPIX
      integer MAX_SCANS
      parameter(NPIX=128)
      parameter(MAX_SCANS=9000)

C-----------
      integer iyy,imm,idd,izz,minute,idseq
      INTEGER*2 year(MAX_SCANS)
      INTEGER*2 month(MAX_SCANS),day(MAX_SCANS)
      INTEGER*2 hour(MAX_SCANS)
      INTEGER*2 minu(MAX_SCANS)

      character*300 input_seq
      character*13  etxt1
      character*8   etxt2
      character*5   etxt3

      integer i, j, k
      integer iscan, ipix, icount, ncount, nwrite

C----------- rain_spc variables
      real lon(NPIX, MAX_SCANS)
      real lat(NPIX, MAX_SCANS)
      real rain_status2(NPIX, MAX_SCANS)
      real sfctype2(NPIX, MAX_SCANS)
      real sfcrain(NPIX, MAX_SCANS)

      integer*4  iseqm(MAX_SCANS)
      real zlon(npix),zlat(npix),zitoil(npix)
      real zirflg(npix), rainpct3785, rain_SSpol(npix)

C----------- Function Declarations
      integer iargc,n_arguments

C----------------------------------------
      etxt2 = 'Scan No:'
      etxt3 = 'ipix:'

      n_arguments = iargc()
      if (n_arguments .eq. 1) then
        call getarg(1, input_seq)
      else
        stop
      endif

C----------- Open input_seq
      iscan=0
      
      open(10, file=input_seq, access='SEQUENTIAL', 
     *     form='UNFORMATTED', status='OLD', err=103)

      read(10, err=105) ncount     
      write(*,*)'Scan No. : ',ncount

      if ( ( ncount.le.0 ).or.( ncount.gt.MAX_SCANS ) ) go to 202

      do icount=1,ncount

         do ipix=1,NPIX
            read(10, err=105) 
     &           zlon(ipix),zlat(ipix),iseqm(icount),zitoil(ipix),
     &           zirflg(ipix), rainpct3785, rain_SSpol(ipix)
            
            if ((zlon(ipix).lt.-180.0).or.(zlon(ipix).gt.360.0)) then
               etxt1 = '[Error: lon]'
               write(*,600) etxt1,etxt2,icount,etxt3,ipix,zlon(ipix)
            endif

            if ((zlat(ipix).lt.-90.0).or.(zlat(ipix).gt.90.0)) then
               etxt1 = '[Error: lat]'
               write(*,600) etxt1,etxt2,icount,etxt3,ipix,zlat(ipix)
            endif                                 
         end do

         iscan=iscan+1
         idseq=iseqm(icount)
         
         call symdhm(iyy,imm,idd,izz,minute,idseq)
c         write(*,*) 'date check ',iyy,imm,idd,izz,minute

         year(iscan)=iyy
         month(iscan)=imm
         day(iscan)=idd
         hour(iscan)=izz           
         minu(iscan)=minute  
         
         do ipix=1,NPIX
            lon(ipix,iscan)=zlon(ipix)
            lat(ipix,iscan)=zlat(ipix)            
            rain_status2(ipix,iscan)=ZIRFLG(ipix)
            sfctype2(ipix,iscan)=ZITOIL(ipix)
            sfcrain(ipix,iscan)=rain_SSpol(ipix)
         enddo
         
      enddo

      close(10)

C------------------------------------------      
      nwrite=iscan
      write(*,*)'Scan No.(nwrite) : ',nwrite

      do j=1,nwrite
         do i=1,NPIX
            if (sfcrain(i,j).gt.0.0) then
               write(*,*) lon(i,j),lat(i,j),sfcrain(i,j)
            endif
         enddo
      enddo

      call exit(0)

C-------- (exit code) ---------------------
 103  write(*,*) " [E] 103: open error (input) "//input_seq
      close(10)
      call exit(103)

 105  write(*,*) " [E] 105: read error (input) "//input_seq
      close(10)
      call exit(105)

 202  write(*,*) " [E] 202: scan number error" 
      close(10)
      call exit(202)

 600  format(A13,1X,A8,1X,I5,3X,A5,1X,I3,5X,F9.4)

      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine stoymd(xx,y,m,d)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculating date of gregorian calendar               c
c     from sequence day number since 1801.01.01            c
c     xx  ;  sequence number  (integer type)               c
c     y   ;  year (ex.1983)   (integer type)               c
c     m   ;  month            (integer type)               c
c     d   ;  day              (integer type)               c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
      subroutine ymdtos(x,yy,mm,dd)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculating sequence day number since 1801.01.01     c
c     from date of gregorian calendar                      c
c     x   ;   sequence number  (integer type)              c
c     yy  ;   year (ex.1983)   (integer type)              c
c     mm  ;   month            (integer type)              c
c     dd  ;   day              (integer type)              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
      subroutine    ymdhts( iy, im, id, ih, irecsq)
c
cccccccc  document  cccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculate record sequence no for year, month, day, hour
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call ymdtos( is, iy, im, id)
c
      irecsq = ( is - 1) * 24 + ih
c
      return
      end
      subroutine    stymdh( iy, im, id, ih, irecsq)
c
cccccccc  document  cccccccccccccccccccccccccccccccccccccccccccccccccccc
c         determine year, month, day, hour from record sequence no
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      is =irecsq / 24 + 1
      ih = mod( irecsq, 24)
c
      call stoymd( is, iy, im, id)
c
      return
      end
      subroutine symdhm( iyear, month, iday, ihour, minute, iseq)
c
cccccccc  document  cccccccccccccccccccccccccccccccccccccccccccccccccccc
c         calculate sequnence no( iseq)
c           for iyear, month, iday, ihour and minute)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      minute = mod( iseq, 60 )
      call stymdh( iyear, month, iday, ihour, iseq/ 60)
      return
      end
      subroutine ymdhms( iyear, month, iday, ihour, minute, iseq)
c
cccccccc  document  cccccccccccccccccccccccccccccccccccccccccccccccccccc
c         calculate sequence no ( iseq)
c           for iyear, month, iday, ihour, minute
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call ymdhts( iyear, month, iday, ihour, iseq)
      iseq = iseq * 60 + minute
      return
      end
