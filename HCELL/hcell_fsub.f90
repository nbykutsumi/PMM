module hcell_fsub

!------------------------------
CONTAINS

!***************************************
SUBROUTINE obt2grid_sumnum(a1in, a1lon, a1lat, lon_first, lat_first, dlon, dlat, nx, ny, nl, a2sum, a2num)
implicit none
!---- input ------
integer                           nl,nx,ny
!f2py intent(in)                  nl,nx,ny
double precision,dimension(nl) :: a1in, a1lat, a1lon
!f2py intent(in)                  a1in, a1lat, a1lon
double precision                  dlon, dlat
!f2py intent(in)                  dlon, dlat
double precision                  lon_first, lat_first
!f2py intent(in)                  lon_first, lat_first
!---- out --------
integer,dimension(nx,ny)       :: a2sum,a2num
!f2py intent(out)                 a2sum,a2num
!---- calc -------
integer                           il, iy, ix
double precision                  lon, lat
!-----------------
a2sum   = 0
a2num   = 0
do il = 1,nl
    lon = a1lon(il)
    lat = a1lat(il)
    ix  = floor((mod(lon+360., 360.) - lon_first) / dlon) + 1
    iy  = floor((lat - lat_first)/dlat) + 1
    a2num(ix,iy) = a2num(ix,iy) + 1
    a2sum(ix,iy) = a2sum(ix,iy) + a1in(il)

end do 
return
END SUBROUTINE

!***************************************
SUBROUTINE obt2grid_count(a1lon, a1lat, lon_first, lat_first, dlon, dlat, nx, ny, nl, a2out)
implicit none
!---- input ------
integer                           nl,nx,ny
!f2py intent(in)                  nl,nx,ny
double precision,dimension(nl) :: a1lat, a1lon
!f2py intent(in)                  a1lat, a1lon
double precision                  dlon, dlat
!f2py intent(in)                  dlon, dlat
double precision                  lon_first, lat_first
!f2py intent(in)                  lon_first, lat_first
!---- out --------
integer,dimension(nx,ny)       :: a2out
!f2py intent(out)                 a2out
!---- calc -------
integer                           il, iy, ix
double precision                  lon, lat
!-----------------
a2out   = 0
do il = 1,nl
    lon = a1lon(il)
    lat = a1lat(il)
    ix  = floor((mod(lon+360., 360.) - lon_first) / dlon) + 1
    iy  = floor((lat - lat_first)/dlat) + 1
    a2out(ix,iy) = a2out(ix,iy) + 1
end do 
return
END SUBROUTINE
!***************************************
!***************************************
SUBROUTINE obt2grid_hist(a1in, a1lon, a1lat, a1bin, lon_first, lat_first, dlon, dlat, nx, ny, nl, nbin, a3sum, a3num)
implicit none
!---- input ------
integer                           nl,nbin,nx,ny
!f2py intent(in)                  nl,nbin,nx,ny
double precision,dimension(nl) :: a1in, a1lat, a1lon
!f2py intent(in)                  a1in, a1lat, a1lon
double precision,dimension(nbin) :: a1bin
!f2py intent(in)                    a1bin
double precision                  dlon, dlat
!f2py intent(in)                  dlon, dlat
double precision                  lon_first, lat_first
!f2py intent(in)                  lon_first, lat_first
!---- out --------
integer,dimension(nx,ny,nbin+1)   :: a3sum, a3num
!f2py intent(out)                    a3sum, a3num


!---- calc -------
integer                           il, ibin, iy, ix
double precision                  lon, lat
!-----------------
a3sum   = 0
a3num   = 0
do il = 1,nl
    lon = a1lon(il)
    lat = a1lat(il)

    ix  = floor((mod(lon+360., 360.) - lon_first) / dlon) + 1
    iy  = floor((lat - lat_first)/dlat) + 1

    do ibin = 1,nbin
        if (a1in(il) .lt. a1bin(ibin))then
            a3sum(ix,iy,ibin) = a3sum(ix,iy,ibin)+a1in(il)
            a3num(ix,iy,ibin) = a3num(ix,iy,ibin)+1
            exit
        else if (ibin.eq.nbin) then
            a3sum(ix,iy,ibin+1) = a3sum(ix,iy,ibin+1)+a1in(il)
            a3num(ix,iy,ibin+1) = a3num(ix,iy,ibin+1)+1
            exit
        else
            cycle
        end if
    end do 

end do 
return
END SUBROUTINE
!***************************************



!***************************************
SUBROUTINE obt2grid_hist_old(a1in, a1lon, a1lat, a1bin, lon_first, lat_first, dlon, dlat, nx, ny, nl, nbin, a3out)
implicit none
!---- input ------
integer                           nl,nbin,nx,ny
!f2py intent(in)                  nl,nbin,nx,ny
double precision,dimension(nl) :: a1in, a1lat, a1lon
!f2py intent(in)                  a1in, a1lat, a1lon
double precision,dimension(nbin) :: a1bin
!f2py intent(in)                    a1bin
double precision                  dlon, dlat
!f2py intent(in)                  dlon, dlat
double precision                  lon_first, lat_first
!f2py intent(in)                  lon_first, lat_first
!---- out --------
integer,dimension(nx,ny,nbin+1)   :: a3out
!f2py intent(out)                    a3out
!---- calc -------
integer                           il, ibin, iy, ix
double precision                  lon, lat
!-----------------
a3out   = 0
do il = 1,nl
    lon = a1lon(il)
    lat = a1lat(il)


    ix  = floor((mod(lon+360., 360.) - lon_first) / dlon) + 1
    iy  = floor((lat - lat_first)/dlat) + 1

    do ibin = 1,nbin
        if (a1in(il) .lt. a1bin(ibin))then
            a3out(ix,iy,ibin) = a3out(ix,iy,ibin)+1
            exit
        else if (ibin.eq.nbin) then
            a3out(ix,iy,ibin+1) = a3out(ix,iy,ibin+1)+1
            exit
        else
            cycle
        end if
    end do 

end do 
return
END SUBROUTINE
!***************************************

end module hcell_fsub 
