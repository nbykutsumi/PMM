module temp_fsub
contains

function f(lon, dlon)
  double precision        lon, dlon
  !f2py intent(in)
  integer                 f
  !f2py intent(out)       f
  f = floor( mod( lon + 360.0d0, 360.0d0)/dlon) + 1
 
end function f

end module temp_fsub
