module pmm_fsub

!---------------------------------
CONTAINS
!*************************************************
!*********************************************************
SUBROUTINE saone_ra2org_ra(a2in, a2out)
implicit none
!---- input --------
real,dimension(32,28)         :: a2in
!f2py intent(in)                   a2in
!---- out ----------
real,dimension(3200,2800)       :: a2out
!f2py intent(out)                  a2out
!---- calc ---------
integer            ix,iy, iix,iiy
!-------------------
do iy = 1,28
  do ix = 1,32
    do iiy = (iy-1)*100+1, iy*100
      do iix = (ix-1)*100+1, ix*100
        a2out(iix,iiy) = a2in(ix,iy)
      end do
    end do
  end do
end do
!-------------------
return
END SUBROUTINE saone_ra2org_ra
!*************************************************
SUBROUTINE obt_match_map(a2map, lat_first, lon_first, dlat, dlon, a1lat, a1lon, nl, ny, nx, a1out)
implicit none
!-----------------------------------
integer                               nl, ny, nx
!---- in --------
double precision,dimension(nx,ny)  :: a2map
!real,dimension(nx,ny)  :: a2map
!f2py intent(in)                      a2map
double precision                      lat_first, lon_first
!f2py intent(in)                      lat_first, lon_first
double precision                      dlat, dlon
!f2py intent(in)                      dlat, dlon


double precision,dimension(nl)     :: a1lat, a1lon
!f2py intent(in)                      a1lat, a1lon  ! lon: 0-360
!---- out --------
double precision,dimension(nl)     :: a1out
!f2py intent(out)                     a1out
!---- calc -------
integer                               i, y, x
double precision                      lat, lon
double precision                      LBLAT, LBLON  ! lower Boundary
!-----------------
!a1out = 0d0
LBLAT = lat_first - 0.5d0*dlat
LBLON = lon_first - 0.5d0*dlon
do i = 1,nl
  lat = a1lat(i)
  lon = a1lon(i)
  y   = floor((lat - LBLAT)/dlat) +1
  x   = floor((lon - LBLON)/dlon) +1
  a1out(i) = a2map(x,y)
end do

return
END SUBROUTINE obt_match_map

!*************************************************
SUBROUTINE obt2wnpac261x265(a2dat, a2lon, a2lat, nw, nl, a2sum, a2num)
implicit none
!--------------------------------------
! On Western North Pacific Cloud Data grid boxes
! Output data
! BBox = [[-0.1, 113.875],[52.1, 180.125]]
! with 0.20 (lat) x 0.25 (lon) degree resolution Y x X =(261 x 265)
!--------------------------------------
integer                          nl, nw
!---- in -----------
real,dimension(nw,nl)         :: a2dat, a2lat, a2lon
!f2py intent(in)                 a2dat, a2lat, a2lon
!---- out -----------
real,dimension(265,261)       :: a2sum, a2num
!f2py intent(out)                a2sum, a2num
!---- calc ----------
integer                          il,iw
integer                          ix,iy
real                             dat,lon,lat
!---- para ----------
real,parameter                :: dlat  = 0.20
real,parameter                :: dlon  = 0.25
real,parameter                :: lat_first = -0.1
real,parameter                :: lat_last  = 52.1
real,parameter                :: lon_first = 113.875
real,parameter                :: lon_last  = 180.125

!---------------------
a2sum    = 0.0
a2num    = 0.0
do il = 1,nl
  do iw = 1,nw
    dat = a2dat(iw,il)
    !------------
    if (dat .lt. 0.0)then
      cycle
    end if
    !------------
    lon = a2lon(iw,il)
    lat = a2lat(iw,il)
    if ((lat .ge. lat_first).and.(lat .le. lat_last)&
         .and.(lon .ge. lon_first).and.(lon .le. lon_last))then
      ix  = floor( (mod(lon+360., 360.) - lon_first) / dlon) +1
      iy  = floor((lat - lat_first)/dlat) + 1
    else
      cycle
    end if
    a2sum(ix,iy) = a2sum(ix,iy) + dat
    a2num(ix,iy) = a2num(ix,iy) + 1.0
  end do
end do

return
END SUBROUTINE obt2wnpac261x265

!*************************************************
SUBROUTINE obt2jp2800x3200(a2dat, a2lon, a2lat, nw, nl, a2sum, a2num)
implicit none
!--------------------------------------
! Output data
! BBox = [[20., 118.],[48., 150.]]
! with 0.01 degree resolution Y x X =(2800x3200)
!--------------------------------------
integer                          nl, nw
!---- in -----------
real,dimension(nw,nl)         :: a2dat, a2lat, a2lon
!f2py intent(in)                 a2dat, a2lat, a2lon
!---- out -----------
real,dimension(3200,2800)     :: a2sum, a2num
!f2py intent(out)                a2sum, a2num
!---- calc ----------
integer                          il,iw
integer                          ix,iy
real                             dat,lon,lat
!---- para ----------
real,parameter                :: dlat  = 0.01
real,parameter                :: dlon  = 0.01
real,parameter                :: lat_first = 20.0
real,parameter                :: lat_last  = 48.0
real,parameter                :: lon_first = 118.0
real,parameter                :: lon_last  = 150.0

!---------------------
a2sum    = 0.0
a2num    = 0.0
do il = 1,nl
  do iw = 1,nw
    dat = a2dat(iw,il)
    !------------
    if (dat .lt. 0.0)then
      cycle
    end if
    !------------
    lon = a2lon(iw,il)
    lat = a2lat(iw,il)
    if ((lat .ge. lat_first).and.(lat .le. lat_last)&
         .and.(lon .ge. lon_first).and.(lon .le. lon_last))then
      ix  = floor( (mod(lon+360., 360.) - lon_first) / dlon) +1
      iy  = floor((lat - lat_first)/dlat) + 1
    else
      cycle
    end if
    a2sum(ix,iy) = a2sum(ix,iy) + dat
    a2num(ix,iy) = a2num(ix,iy) + 1.0
  end do
end do

return
END SUBROUTINE obt2jp2800x3200

!*************************************************
SUBROUTINE obt2jp280x320(a2dat, a2lon, a2lat, nw, nl, a2sum, a2num)
implicit none
!--------------------------------------
! Output data
! BBox = [[20., 118.],[48., 150.]]
! with 0.01 degree resolution Y x X =(280x320)
!--------------------------------------
integer                          nl, nw
!---- in -----------
real,dimension(nw,nl)         :: a2dat, a2lat, a2lon
!f2py intent(in)                 a2dat, a2lat, a2lon
!---- out -----------
real,dimension(320,280)       :: a2sum, a2num
!f2py intent(out)                a2sum, a2num
!---- calc ----------
integer                          il,iw
integer                          ix,iy
real                             dat,lon,lat
!---- para ----------
real,parameter                :: dlat  = 0.1
real,parameter                :: dlon  = 0.1
real,parameter                :: lat_first = 20.0
real,parameter                :: lat_last  = 48.0
real,parameter                :: lon_first = 118.0
real,parameter                :: lon_last  = 150.0

!---------------------
a2sum    = 0.0
a2num    = 0.0
do il = 1,nl
  do iw = 1,nw
    dat = a2dat(iw,il)
    !------------
    if (dat .lt. 0.0)then
      cycle
    end if
    !------------
    lon = a2lon(iw,il)
    lat = a2lat(iw,il)
    if ((lat .ge. lat_first).and.(lat .le. lat_last)&
         .and.(lon .ge. lon_first).and.(lon .le. lon_last))then
      ix  = floor( (mod(lon+360., 360.) - lon_first) / dlon) +1
      iy  = floor((lat - lat_first)/dlat) + 1
    else
      cycle
    end if
    a2sum(ix,iy) = a2sum(ix,iy) + dat
    a2num(ix,iy) = a2num(ix,iy) + 1.0
  end do
end do

return
END SUBROUTINE obt2jp280x320




!*************************************************
SUBROUTINE obt2dec(a2dat, a2lon, a2lat, nw, nl, a2sum, a2num)
implicit none
!--------------------------------------
integer                          nl, nw
!---- in -----------
real,dimension(nw,nl)         :: a2dat, a2lat, a2lon
!f2py intent(in)                 a2dat, a2lat, a2lon
!---- out -----------
real,dimension(3600,1800)     :: a2sum, a2num
!f2py intent(out)                a2sum, a2num
!---- calc ----------
integer                          il,iw
integer                          ix,iy
real                             dat,lon,lat
!---- para ----------
real,parameter                :: dlat  = 0.1
real,parameter                :: dlon  = 0.1
real,parameter                :: lat_first = -90.0
!---------------------
a2sum    = 0.0
a2num    = 0.0
do il = 1,nl
  do iw = 1,nw
    dat = a2dat(iw,il)
    !------------
    if (dat .lt. 0.0)then
      cycle
    end if
    !------------
    lon = a2lon(iw,il)
    lat = a2lat(iw,il)
    
    ix  = floor( mod(lon+360., 360.)*100. / (dlon*100.)) +1
    iy  = floor((lat - lat_first)/dlat) + 1

    a2sum(ix,iy) = a2sum(ix,iy) + dat
    a2num(ix,iy) = a2num(ix,iy) + 1.0
  end do
end do

return
END SUBROUTINE obt2dec
!*************************************************
SUBROUTINE obt2saone(a2dat, a2lon, a2lat, nw, nl, a2sum, a2num)
implicit none
!--------------------------------------
integer                          nl, nw
!---- in -----------
real,dimension(nw,nl)         :: a2dat, a2lat, a2lon
!f2py intent(in)                 a2dat, a2lat, a2lon
!---- out -----------
real,dimension(360,180)       :: a2sum, a2num
!f2py intent(out)                a2sum, a2num
!---- calc ----------
integer                          il,iw
integer                          ix,iy
real                             dat,lon,lat
!---- para ----------
real,parameter                :: dlat  = 1.0
real,parameter                :: dlon  = 1.0
real,parameter                :: lat_first = -90.0
!---------------------
a2sum    = 0.0
a2num    = 0.0
do il = 1,nl
  do iw = 1,nw
    dat = a2dat(iw,il)
    !------------
    if (dat .lt. 0.0)then
      cycle
    end if
    !------------
    lon = a2lon(iw,il)
    lat = a2lat(iw,il)
    
    ix  = mod( floor(lon + 360.0), 360) + 1
    iy  = floor((lat - lat_first)/dlat) + 1
    a2sum(ix,iy) = a2sum(ix,iy) + dat
    a2num(ix,iy) = a2num(ix,iy) + 1.0
  end do
end do

return
END SUBROUTINE obt2saone
!******************************************
end module pmm_fsub
