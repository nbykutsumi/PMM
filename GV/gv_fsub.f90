module gv_fsub


!-------------------
CONTAINS
!**********************************************************
! SUBROUTINE & FUNCTION
!**********************************************************
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
  double precision,dimension(ny,nx)  :: a2out
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

        dist = hubeny(lat0, lon0, lat1, lon1) /1000d0  ! km
        w    = 1d0 / dist**power
        w_sum = w_sum + w
        wv_sum= wv_sum + w * v

      end do
      a2out(ix,iy) = wv_sum / w_sum


    end do
  end do

RETURN
END SUBROUTINE point2map
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
