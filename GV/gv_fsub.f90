module gv_fsub


!-------------------
CONTAINS
!**********************************************************
! SUBROUTINE & FUNCTION
!**********************************************************
SUBROUTINE fill_esurf(a2joinprof, miss, nl, nh, a2out)
  implicit none
  !-------------
  integer                       nl, nh
  !--- input ---
  real,dimension(nh,nl)      :: a2joinprof
  !f2py intent(in)              a2joinprof

  real                          miss
  !f2py intent(in)              miss

  !--- output --
  real,dimension(nh,nl)      :: a2out
  !f2py intent(out)             a2out
  !--- calc ----
  real                          esurf
  integer                       il,ih, doneflag
  !-------------
  a2out= a2joinprof
  do il = 1,nl
    esurf    = a2joinprof(1,il)
    if (esurf .eq. miss)then
      esurf=0
    end if

    !-- search nsurf bin  --
    doneflag = 0
    do ih = 2,nh
      if (doneflag .ne. 0) exit
      if (a2joinprof(ih,il).eq.miss)then
        a2out(ih,il) = esurf
      else
        doneflag = 1
      end if
    end do
  end do
  !-------------
  return
END SUBROUTINE

!**********************************************************
SUBROUTINE fill_esurf_interp(a2joinprof, miss, nl, nh, a2out)
  implicit none
  !-------------
  integer                       nl, nh
  !--- input ---
  real,dimension(nh,nl)      :: a2joinprof
  !f2py intent(in)              a2joinprof

  real                          miss
  !f2py intent(in)              miss

  !--- output --
  real,dimension(nh,nl)      :: a2out
  !f2py intent(out)             a2out
  !--- calc ----
  real                          esurf, nsurf
  integer                       il,ih
  integer                       insurf
  !-------------
  a2out= a2joinprof
  do il = 1,nl
    esurf    = a2joinprof(1,il)
    if (esurf .eq. miss)then
      esurf=0
    end if

    !-- search nsurf bin  --
    insurf   = 0
    do ih = 2,nh
      if (insurf .ne. 0) exit
      if (a2joinprof(ih,il).ne.miss)then
        insurf = ih
      end if
    end do
    !-- linear interp --
    if ((insurf.ne.0).and.(insurf.ne.nh))then
      nsurf = a2joinprof(insurf,il)
      do ih = 2,nh
        a2out(ih,il) =  ((ih-1)*nsurf + (insurf-ih)*esurf )/(insurf-1)
      end do
    end if
  end do
  !-------------
  return
END SUBROUTINE



!**********************************************************
SUBROUTINE extract_slice_clusterprof(a2in, a1iidxpy, nxout, nx, ny, a2out)
  implicit none
  !-------------
  integer                       nx, ny
  !--- input ---
  real,dimension(nx,ny)      :: a2in
  !f2py intent(in)              a2in
  integer,dimension(ny)      :: a1iidxpy
  !f2py intent(in)              a1iidxpy
  integer                       nxout
  !f2py intent(in)              nxout
  !--- output --
  real,dimension(nxout,ny)   :: a2out
  !f2py intent(out)             a2out
  !--- calc ----
  integer                       x,y, x0
  integer                       ix,ex,num
  real,parameter             :: miss=-9999
  !-------------
  a2out= miss
  do y = 1,ny
    num  = 0
    ix   = a1iidxpy(y)+1
    ex   = ix +  nxout -1
    if (ex.gt.nx)then
      ex = nxout
    end if

    x0 = 0
    do x = ix,ex
      x0 = x0+1
      a2out(x0,y) = a2in(x,y)
    end do
  end do
  !-------------
  return
END SUBROUTINE
!**********************************************************
SUBROUTINE mean_slice_negativemask(a2in, a1iidxpy, nxave, nx, ny, a1out)
  implicit none
  !-------------
  integer                 nx, ny
  !--- input ---
  double precision,dimension(nx,ny)  :: a2in
  !f2py intent(in)                      a2in
  integer,dimension(ny)              :: a1iidxpy
  !f2py intent(in)                      a1iidxpy
  integer                               nxave
  !f2py intent(in)                      nxave
  !--- output --
  double precision,dimension(ny)     :: a1out
  !f2py intent(out)                     a1out
  !--- calc ----
  integer                               x,y
  integer                               ix,ex,num
  double precision                      vsum
  double precision,parameter         :: miss=-9999d0 
  !-------------
  a1out= miss
  do y = 1,ny
    vsum = 0d0
    num  = 0
    ix   = a1iidxpy(y)+1
    ex   = ix +  nxave -1
    if (ex.gt.nx)then
      ex = nx
    end if

    do x = ix,ex
      if (a2in(x,y) .ge. 0d0)then
        vsum = vsum + a2in(x,y)
        num  = num + 1
      end if
    end do
    if (num.gt.0)then
      a1out(y) = vsum/num
    else
      a1out(y) = miss
    end if
  end do
  !-------------
  return
END SUBROUTINE
!*****************************************************************
SUBROUTINE gauge_match_pyxy(a2satelon, a2satelat, glon, glat, thdist, nx, ny, a1x, a1y)
  implicit none

  !-----------
  integer                               nx, ny
  !-- input --
  double precision,dimension(nx,ny)  :: a2satelon, a2satelat
  !f2py intent(in)                      a2satelon, a2satelat
  double precision                      glon, glat, thdist
  !f2py intent(in)                      glon, glat, thdist

  !-- output--
  integer,dimension(nx*ny)           :: a1x, a1y
  !f2py intent(out)                     a1x, a1y

  !-- calc --
  integer                               x, y, nout
  double precision                      lon, lat, dist
  !--

  a1x  = -9999
  a1y  = -9999
 
  nout = 0 
  do y = 1,ny
    do x = 1, nx
      lon = a2satelon(x,y)
      lat = a2satelat(x,y)
      dist = hubeny(glat, glon, lat, lon) /1000d0 ! km
      if (dist .le. thdist) then
        nout = nout +1
        a1x(nout) = x -1   ! index for python
        a1y(nout) = y -1   ! index for python
      end if
    end do
  end do 

RETURN
END SUBROUTINE gauge_match_pyxy

!*****************************************************************
SUBROUTINE point2map_distmask(a1dat, a1lon_obs, a1lat_obs, a1lon_map, a1lat_map, power, thdist, nobs, nx, ny, a2out)

  implicit none
  !------------------------------------
  integer                               ny, nx, nobs
  
  !-- input --
  double precision,dimension(nobs)   :: a1dat, a1lon_obs, a1lat_obs
  !f2py intent(in)                      a1dat, a1lon_obs, a1lat_obs
  double precision,dimension(nx)     :: a1lon_map
  !f2py intent(in)                      a1lon_map

  double precision,dimension(ny)     :: a1lat_map
  !f2py intent(in)                      a1lat_map

  double precision                      power
  !f2py intent(in)                      power
  double precision                      thdist  ! km
  !f2py intent(in)                      thdist

  !-- calc --
  double precision                      v, w, w_sum, wv_sum, dist, dist_min
  double precision                      lon0, lat0, lon1, lat1
  integer                               ix, iy, iobs

  !-- out --
  double precision,dimension(nx,ny)  :: a2out
  !f2py intent(out)                     a2out

  !--------------
  do ix = 1,nx
    do iy = 1,ny
      lon0 = a1lon_map(ix)
      lat0 = a1lat_map(iy) 

      wv_sum =0d0
      w_sum  =0d0
      dist_min= 1d5
      do iobs = 1,nobs
        lon1 = a1lon_obs(iobs)
        lat1 = a1lat_obs(iobs)
        v    = a1dat(iobs)
        dist = hubeny(lat0, lon0, lat1, lon1) /1000d0 ! km
        w    = 1d0 / (dist**power + 1d-3)   ! 1d-3 is to avoid dividing by zero
        w_sum = w_sum + w
        wv_sum= wv_sum + w * v
        
        if (dist.lt.dist_min)then
          dist_min = dist
        end if

      end do
      if (dist_min .le. thdist)then
        a2out(ix,iy) = wv_sum / w_sum
      else
        a2out(ix,iy) = -9999d0
      end if

    end do

  end do


RETURN
END SUBROUTINE point2map_distmask
!*****************************************************************


SUBROUTINE point2map(a1dat, a1lon_obs, a1lat_obs, a1lon_map, a1lat_map, power, nobs, nx, ny, a2out)

  implicit none
  !------------------------------------
  integer                               ny, nx, nobs
  
  !-- input --
  double precision,dimension(nobs)   :: a1dat, a1lon_obs, a1lat_obs
  !f2py intent(in)                      a1dat, a1lon_obs, a1lat_obs
  double precision,dimension(nx)     :: a1lon_map
  !f2py intent(in)                      a1lon_map

  double precision,dimension(ny)     :: a1lat_map
  !f2py intent(in)                      a1lat_map

  double precision                      power
  !f2py intent(in)                      power

  !-- calc --
  double precision                      v, w, w_sum, wv_sum, dist
  double precision                      lon0, lat0, lon1, lat1
  integer                               ix, iy, iobs

  !-- out --
  double precision,dimension(nx,ny)  :: a2out
  !f2py intent(out)                     a2out

  !--------------
  do ix = 1,nx
    do iy = 1,ny
      lon0 = a1lon_map(ix)
      lat0 = a1lat_map(iy) 

      wv_sum =0d0
      w_sum  =0d0
      do iobs = 1,nobs
        lon1 = a1lon_obs(iobs)
        lat1 = a1lat_obs(iobs)
        v    = a1dat(iobs)

        dist = hubeny(lat0, lon0, lat1, lon1) /1000d0 ! km
        w    = 1d0 / (dist**power + 1d-3)   ! 1d-3 is to avoid dividing by zero
        w_sum = w_sum + w
        wv_sum= wv_sum + w * v

      end do
      a2out(ix,iy) = wv_sum / w_sum


    end do

  end do


RETURN
END SUBROUTINE point2map
!*****************************************************************
SUBROUTINE point2mindistmap(a1lon_obs, a1lat_obs, a1lon_map, a1lat_map, nobs, nx, ny, a2out)

  implicit none
  !------------------------------------
  integer                               ny, nx, nobs
  
  !-- input --
  double precision,dimension(nobs)   :: a1lon_obs, a1lat_obs
  !f2py intent(in)                      a1lon_obs, a1lat_obs
  double precision,dimension(nx)     :: a1lon_map
  !f2py intent(in)                      a1lon_map

  double precision,dimension(ny)     :: a1lat_map
  !f2py intent(in)                      a1lat_map

  !-- calc --
  double precision                      dist, dist_min
  double precision                      lon0, lat0, lon1, lat1
  integer                               ix, iy, iobs

  !-- out --
  double precision,dimension(nx,ny)  :: a2out
  !f2py intent(out)                     a2out

  !--------------
  do ix = 1,nx
    do iy = 1,ny
      lon0 = a1lon_map(ix)
      lat0 = a1lat_map(iy) 

      dist_min = 1d5
      do iobs = 1,nobs
        lon1 = a1lon_obs(iobs)
        lat1 = a1lat_obs(iobs)

        dist = hubeny(lat0, lon0, lat1, lon1) /1000d0 ! km
        if (dist .lt. dist_min) then
            dist_min = dist
        end if
      end do

      a2out(ix,iy) = dist_min

    end do
  end do


RETURN
END SUBROUTINE point2mindistmap
!*****************************************************************



FUNCTION hubeny(lat1, lon1, lat2, lon2)
  implicit none
  !-- for input -----------
  double precision                      lat1, lon1, lat2, lon2
!f2py intent(in)                        lat1, lon1, lat2, lon2
  !-- for output-----------
  double precision                      hubeny   ! (m)
!f2py intent(out)                       hubeny
  !-- for calc ------------
  double precision,parameter         :: pi = atan(1.0d0)*4.0d0
  double precision,parameter         :: a  = 6378137d0
  double precision,parameter         :: b  = 6356752.314140d0
  double precision,parameter         :: e2 = 0.00669438002301188d0
  double precision,parameter         :: a_1_e2 = 6335439.32708317d0
  double precision                      M, N, W
  double precision                      latrad1, latrad2, lonrad1, lonrad2
  double precision                      latave, dlat, dlon
  double precision                      dlondeg
  !------------------------
  latrad1   = lat1 * pi / 180.0d0
  latrad2   = lat2 * pi / 180.0d0
  lonrad1   = lon1 * pi / 180.0d0
  lonrad2   = lon2 * pi / 180.0d0
  !
  latave    = (latrad1 + latrad2)/2.0d0
  dlat      = latrad2 - latrad1
  dlon      = lonrad2 - lonrad1
  !
  dlondeg   = lon2 - lon1
  if ( abs(dlondeg) .gt. 180.0d0) then
    dlondeg = 180.0d0 - mod(abs(dlondeg), 180.0d0)
    dlon    = dlondeg * pi / 180.0d0
  end if
  !-------
  W  = sqrt(1.0d0 - e2 * sin(latave)**2.0d0 )
  M  =  a_1_e2 / (W**3.0d0)
  N  =  a / W
  hubeny  = sqrt( (dlat * M)**2.0d0 + (dlon * N * cos(latave))**2.0d0 )
  !print *, hubeny * 0.001d0
  !print *, "a_1_e2=", a_1_e2
  !print *, "calc",a*(1.0d0 - e2)

  !M  = 6334834.0d0 / sqrt( (1.0d0 - 0.006674d0 * sin(latave) **2.0d0)**2.0d0 )
  !N  = 6377397.0d0 / sqrt( 1.0d0 - 0.006674d0 * sin(latave) **2.0d0)
  !hubeny = sqrt( (M * dlat )**2.0d0 + ( N * cos(latave) * dlon)**2.0d0 )
RETURN
END FUNCTION hubeny
!*****************************************************************

end module gv_fsub
