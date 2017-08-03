from numpy import *
from PMM.pmm_fsub import *
import myfunc.IO.CLOUDTYPE as CLOUDTYPE

ver     = ".M.3"
cl      = CLOUDTYPE.MyCloudWNP(ver=ver)
[[lllat,lllon],[urlat,urlon]] = cl.BBox
dlat = cl.dLat
dlon = cl.dLon
miss_out = -9999.

a1lat= arange(-3.0,3.0,1.0)
a1lon= arange(140,146,1.0)

a1Lat= arange(lllat,urlat,dlat)
a1Lon= arange(lllon,urlon,dlon)

Lon,Lat = meshgrid(a1Lon,a1Lat)
a2dat= Lat
shapeOut = [10]+list(a2dat.shape)
a3dat= resize(a2dat, shapeOut)

out = pmm_fsub.pickup_data(a3dat.T, lllat, lllon, urlat, urlon, dlat, dlon, a1lon, a1lat, miss_out).T
print "*"*50
print a1lat
print "*"*50
print a3dat
print "*"*50
print out
