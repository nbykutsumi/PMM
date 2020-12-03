cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c param_grid_MIS.inc
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
c     PARAMETER(NPIX=180,MAX_SCANS=9000)   ! Rev 2014.06.18
      PARAMETER(NPIX=180,MAX_SCANS=20000)

      integer IM, JM
c     PARAMETER(IM=180,JM=6000)            ! Rev 2013.04.15
c     PARAMETER(IM=180,JM=9000)            ! Rev 2014.06.18
      PARAMETER(IM=NPIX,JM=MAX_SCANS)

c---  set grid parameter
      integer IDIM, JDIM
c     PARAMETER(IDIM=3600,JDIM=1200)       ! Rev 2019.10.03
      PARAMETER(IDIM=3600,JDIM=1800)

C--- Consideration of EFOV
c 考慮するフットプリントの範囲は、85GHzの半径(7.2x4.4)のX倍に
c 0.1度格子の半分の距離を加えた距離。
      real radius1, radius2, mult
c     PARAMETER(radius1=7.0, radius2=6.5, mult=2.5)   ! Rev 2012.10.22
      PARAMETER(radius1=7.0, radius2=6.5, mult=1.0)
