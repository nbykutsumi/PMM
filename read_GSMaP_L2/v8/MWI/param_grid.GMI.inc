cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c param_grid_GMI.inc
c 
c 2019.10.07   NC (GSMaP v8)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c--- missing value
      real UNDEF
      PARAMETER(UNDEF=-999.9)
      integer IUNDEF
      PARAMETER(IUNDEF=-999)

c--- set path parameter
      integer NPIX, MAX_SCANS
c     PARAMETER(NPIX=221,MAX_SCANS=3000)   ! Rev 2013.10.09
      PARAMETER(NPIX=221,MAX_SCANS=3100)   ! Rev 2014.06.18

      integer IM, JM
c     PARAMETER(IM=221,JM=3000)            ! Rev 2014.06.18
      PARAMETER(IM=NPIX,JM=MAX_SCANS)

c---  set grid parameter
      integer IDIM, JDIM
c     PARAMETER(IDIM=3600,JDIM=1200)       ! Rev 2019.10.03
      PARAMETER(IDIM=3600,JDIM=1800)

C--- Consideration of EFOV
c 考慮するフットプリントの範囲は、85GHzの半径(7.2x4.4)のX倍に
c 0.1度格子の半分の距離を加えた距離。
      real radius1, radius2, mult
c      PARAMETER(radius1=3.6, radius2=2.3, mult=2.0)   ! Rev 2012.10.22
c      PARAMETER(radius1=3.6, radius2=2.3, mult=1.0)   ! Rev 2014.05.23
      PARAMETER(radius1=3.6, radius2=2.2, mult=1.0)
