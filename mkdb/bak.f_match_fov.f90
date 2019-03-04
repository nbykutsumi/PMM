module f_match_fov
CONTAINS
!*********************************************************
SUBROUTINE extract_3d(a3in, a2x, a2y, miss_idx, miss_out, nxin, nyin, nzin, nxidx, nyidx, a3out)
!-----------------------
! accepts python type index 0,1,2,...
! input (y,x,z).T = (z,x,y)
!-----------------------
implicit none
!-------------------------------------------------
integer             nxin, nyin, nzin, nxidx, nyidx

!---- in ----
double precision,dimension(nzin,nxin,nyin)   :: a3in
!f2py intent(in)                                a3in
integer,dimension(nxidx,nyidx)          :: a2x, a2y  ! python type index 0,1,2..
!f2py intent(in)                           a2x, a2y
integer                                    miss_idx
!f2py intent(in)                           miss_idx
double precision                           miss_out
!f2py intent(in)                           miss_out

!---- out ---
double precision,dimension(nzin,nxidx,nyidx) :: a3out
!f2py intent(out)                               a3out

!---- cal ---
integer                 iyidx,ixidx,x,y
!------------
do iyidx = 1,nyidx
  do ixidx = 1,nxidx
    x = a2x(ixidx,iyidx)
    y = a2y(ixidx,iyidx)
    if (x.eq.miss_idx)then
      a3out(:,ixidx,iyidx)=miss_out
    else
      a3out(:,ixidx,iyidx)=a3in(:,x+1,y+1)
    end if
  end do
end do
return
END SUBROUTINE
!*****************************************************************



!*********************************************************
SUBROUTINE extract_2d(a2in, a2x, a2y, miss_idx, miss_out, nxin, nyin, nxidx, nyidx, a2out)
!-----------------------
! accepts python type index 0,1,2,...
!-----------------------
implicit none
!-------------------------------------------------
integer             nxin, nyin, nxidx, nyidx

!---- in ----
double precision,dimension(nxin,nyin)   :: a2in
!f2py intent(in)                           a2in
integer,dimension(nxidx,nyidx)          :: a2x, a2y  ! python type index 0,1,2..
!f2py intent(in)                           a2x, a2y
integer                                    miss_idx
!f2py intent(in)                           miss_idx
double precision                           miss_out
!f2py intent(in)                           miss_out

!---- out ---
double precision,dimension(nxidx,nyidx) :: a2out
!f2py intent(out)                          a2out

!---- cal ---
integer                 iyidx,ixidx,x,y
!------------
do iyidx = 1,nyidx
  do ixidx = 1,nxidx
    x = a2x(ixidx,iyidx)
    y = a2y(ixidx,iyidx)
    if (x.eq.miss_idx)then
      a2out(ixidx,iyidx)=miss_out
    else
      a2out(ixidx,iyidx)=a2in(x+1,y+1)
    end if
  end do
end do
return
END SUBROUTINE


!*****************************************************************
SUBROUTINE match_gmi_dpr(a2lon1, a2lat1, a2lon2, a2lat2, nx1, ny1, nx2, ny2, a2x1, a2x2, a2x3, a2x4, a2y1, a2y2, a2y3, a2y4)
implicit none
!-------------------------------------------------
! return X and Y in fortran indexing: 1,2,3,4...., NOT 0,1,2, ..
!-------------------------------------------------
integer         nx1, ny1, nx2, ny2
!---- in ------
double precision,dimension(nx1,ny1)   :: a2lon1, a2lat1
!f2py intent(in)                         a2lon1, a2lat1

double precision,dimension(nx2,ny2)   :: a2lon2, a2lat2
!f2py intent(in)                         a2lon2, a2lat2
!---- out ------
integer,dimension(nx1,ny1)            :: a2x1, a2x2, a2x3, a2x4
!f2py intent(out)                        a2x1, a2x2, a2x3, a2x4
integer,dimension(nx1,ny1)            :: a2y1, a2y2, a2y3, a2y4
!f2py intent(out)                        a2y1, a2y2, a2y3, a2y4


!---- calc -----
integer                                  x1, x2, y1, y2
integer                                  ystart2, yend2, xstart2, xend2
integer                                  ymem2
integer                                  ypre1
integer                                  xpre2, ypre2
integer                                  xtemp2, ytemp2
integer                                  xtemp2_2, ytemp2_2
integer                                  xtemp2_3, ytemp2_3
integer                                  xtemp2_4, ytemp2_4
integer                                  flagfound
double precision                         lat1, lat2, lon1, lon2
double precision                         dist
double precision                         mindist, mindist2, mindist3, mindist4


integer,parameter                     :: miss_int=-9999
integer,parameter                     :: dx2 = 5
integer,parameter                     :: dy2 = 5
integer,parameter                     :: dycurve2 = 9
integer,parameter                     :: dyforward2=300
double precision,parameter            :: thdist = 7   ! km
!-------------------------------
flagfound = 0
ypre1 =miss_int
xpre2 =miss_int
ypre2 =miss_int
ymem2 = miss_int
a2x1   = miss_int
a2y1   = miss_int
a2x2   = miss_int
a2y2   = miss_int
a2x3   = miss_int
a2y3   = miss_int
a2x4   = miss_int
a2y4   = miss_int


do y1 = 1,ny1
!do y1 = 780,785
!do y1 = 30,1508
  !write(*,*) y1

  do x1 = 1, nx1
  !do x1 = 44,44
    lon1 = a2lon1(x1, y1)
    lat1 = a2lat1(x1, y1)

    !if ((lat1.gt.-1).and.(lat1.lt.1))then
    !  write(*,*) y1,x1,lat1,lon1
    !end if

    !-- DPR loop --
    mindist =1d+6   ! [km]
    mindist2=1d+6
    mindist3=1d+6
    mindist4=1d+6   ! [km]

    ! set start & end of y and x 
    if (flagfound.eq.0)then
      !write(*,*) 'AAAAAAAAAAAAAAAAAAAAAAAAAAAaa'
      ystart2 = 1
      xstart2 = 1
      !yend2   = ny2
      !xend2   = nx2
      yend2   = dyforward2
      xend2   = nx2
      !write(*,*) 'a',y1,x1,ystart2,ymem2
    else if (ypre1.ne.y1)then
      !write(*,*) 'BBBBBBBBBBBBBB'
      ystart2 = 1
      ystart2 = max(1, ymem2 - dycurve2)
      xstart2 = 1
      yend2   = min(ny2, ymem2 + dycurve2 + dy2) 
      xend2   = nx2
      !write(*,*) 'a',y1,x1,ystart2,ymem2
    else if (ypre2.eq.miss_int)then
      !write(*,*) 'CCCCCCCCCCCCCC'
      ystart2 = max(1,ymem2 - dycurve2)
      xstart2 = 1
      yend2   = min(ny2, ymem2 + dycurve2 + dy2)
      xend2   = nx2
      !write(*,*) 'b',y1,x1,ystart2,ymem2
    else
      !write(*,*) 'DDDDDDDDDDDDDD'
      ystart2 = max(1, ypre2 - dy2)
      xstart2 = max(1, xpre2 - dx2)
      yend2   = min(ny2, ypre2 + dy2)
      xend2   = min(nx2, xpre2 + dx2)
      !write(*,*) 'c',y1,x1,ystart2,ymem2
    end if 
    !write(*,*) y1,x1, ystart2,yend2,xstart2,xend2 
    !---- start nearest DPR searching loop ----
    ytemp2 = miss_int
    xtemp2 = miss_int
    ytemp2_2=miss_int
    ytemp2_3=miss_int
    ytemp2_4=miss_int
    xtemp2_2=miss_int
    xtemp2_3=miss_int
    xtemp2_4=miss_int

    do y2 = ystart2, yend2
      do x2 = xstart2, xend2
        lon2 = a2lon2(x2, y2)
        lat2 = a2lat2(x2, y2)
        dist = hubeny(lat1, lon1, lat2, lon2)*0.001d0
        !if ((y2.ge.720).and.(y2.le.722))then
        !  if ((x2.ge.3).and.(x2.le.5))then
        !    write(*,*) y2,x2,lat2,lon2,dist
        !  end if
        !end if
        if (dist.gt.thdist)then
          cycle
        end if

        flagfound = flagfound + 1

        if (dist.lt.mindist)then
          mindist4 = mindist3
          xtemp2_4 = xtemp2_3
          ytemp2_4 = ytemp2_3

          mindist3 = mindist2
          xtemp2_3 = xtemp2_2
          ytemp2_3 = ytemp2_2

          mindist2 = mindist
          xtemp2_2 = xtemp2
          ytemp2_2 = ytemp2

          mindist = dist
          xtemp2  = x2
          ytemp2  = y2
          ymem2   = y2
        else if (dist.lt.mindist2)then
          mindist4 = mindist3
          xtemp2_4 = xtemp2_3
          ytemp2_4 = ytemp2_3

          mindist3 = mindist2
          xtemp2_3 = xtemp2_2
          ytemp2_3 = ytemp2_2

          mindist2 = dist
          xtemp2_2 = x2
          ytemp2_2 = y2
        else if (dist.lt.mindist3)then
          mindist4 = mindist3
          xtemp2_4 = xtemp2_3
          ytemp2_4 = ytemp2_3

          mindist3 = dist
          xtemp2_3 = x2
          ytemp2_3 = y2
        else if (dist.lt.mindist4)then
          mindist4 = dist
          xtemp2_4 = x2
          ytemp2_4 = y2
        end if

     end do
    end do
    a2x1(x1,y1) = xtemp2
    a2y1(x1,y1) = ytemp2

    a2x2(x1,y1) = xtemp2_2
    a2y2(x1,y1) = ytemp2_2

    a2x3(x1,y1) = xtemp2_3
    a2y3(x1,y1) = ytemp2_3

    a2x4(x1,y1) = xtemp2_4
    a2y4(x1,y1) = ytemp2_4

    xpre2      = xtemp2
    ypre2      = ytemp2
    ypre1      = y1

    !if ((x1.ge.11).and.(x1.le.50))then
    !  write(*,*) y1,x1,ytemp2,ytemp2_2,ytemp2_3,ytemp2_4,xtemp2,xtemp2_2,xtemp2_3,xtemp2_4
    !end if
  end do
end do 
!-------------------------------
END SUBROUTINE


!!!*****************************************************************
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

!*********************************************************
end module f_match_fov
