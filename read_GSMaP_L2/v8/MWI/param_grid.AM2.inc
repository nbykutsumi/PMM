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
      PARAMETER(NPIX=486,MAX_SCANS=18000)

      integer IM, JM
      PARAMETER(IM=NPIX,JM=MAX_SCANS)

c---  set grid parameter
      integer IDIM, JDIM
c     PARAMETER(IDIM=3600,JDIM=1200)       ! Rev 2019.10.03
      PARAMETER(IDIM=3600,JDIM=1800)

C--- Consideration of EFOV
c 考慮するフットプリントの範囲は、85GHzの半径(7.2x4.4)のX倍に
c 0.1度格子の半分の距離を加えた距離。
      real radius1, radius2, mult
      PARAMETER(radius1=2.5, radius2=1.5, mult=1.0)
