from numpy import *
from detect.detect_fsub import *
import myfunc.IO.RadarAMeDAS as RadarAMeDAS
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
from bisect import bisect_left 

prj  = "ra_0.01"
ra   = RadarAMeDAS.RadarAMeDAS(prj=prj)
Lat  = ra.Lat
Lon  = ra.Lon
ny   = ra.ny
nx   = ra.nx

prj  = "ra_0.1"
ra   = RadarAMeDAS.RadarAMeDAS(prj=prj)
Lat  = ra.Lat
Lon  = ra.Lon
ny   = ra.ny
nx   = ra.nx



#prj  = "0.20x0.25WNP"
#cl   = CLOUDTYPE.CloudWNP()
#Lat  = cl.Lat
#Lon  = cl.Lon
#ny   = cl.ny
#nx   = cl.nx

miss = -9999.

srcDir  = "/tank/utsumi/data/RadarAMeDAS/mask"
srcPath = srcDir + "/RadarSites.csv"
f = open(srcPath,"r")
lines = f.readlines()
f.close()

dlatlon = {}
for sline in lines[1:]:
  line = sline.split(",")
  name = line[0]
  latDeg= float(line[1]) 
  latMin= float(line[2]) 
  lonDeg= float(line[3]) 
  lonMin= float(line[4]) 

  lat   = latDeg + latMin/60.0
  lon   = lonDeg + lonMin/60.0
  dlatlon[name] = [lat,lon]

print dlatlon
a2mask = zeros([ny,nx],float32)
for site in dlatlon.keys():
  print site
  lat, lon = dlatlon[site]

  ytmp  = bisect_left(Lat, lat)
  if (lat - Lat[ytmp-1]) < (Lat[ytmp]-lat):
    y = ytmp-1
  else:
    y = ytmp

  xtmp  = bisect_left(Lon, lon)
  if (lon - Lon[xtmp-1]) < (Lon[xtmp]-lon):
    x = xtmp-1
  else:
    x = xtmp

  a2site =ones([ny,nx],float32)*miss
  a2site[y,x] = 1.0

  if site in ["Naze","Okinawa","Ishigakijima"]:
    radkm = 150.0
  else:
    radkm = 200.0
  a2tmp = detect_fsub.mk_territory_reg(a2site.T, Lon, Lat, radkm*1000., miss, 0.0).T
  a2mask = a2mask + a2tmp

a2mask = ma.masked_equal(a2mask,0.0).filled(miss)
a2mask = ma.masked_greater(a2mask,0.0).filled(1.0)
#oDir   = srcDir
oDir   = "/tank/utsumi/PMM/WNP.261x265/MASK"
oPath  = srcDir + "/RAmask.kubota.%s.%dx%d"%(prj,ny,nx)
a2mask.astype(float32).tofile(oPath)
print oPath
