from numpy import *
import numpy as np
import os, sys, socket
import myfunc.util as util
import myfunc.regrid.Regrid as Regrid

#-- Read GTOPO ----
lullat = [90,40,-10]
lullon = [-180,-140,-100,-60,-20,20,60,100,140]
llatlon0 = [(ullat,ullon) for ullat in lullat for ullon in lullon]

lullon = [-180,-120,-60,0,60,120]
llatlon1 = [(-60,ullon) for ullon in lullon]

llatlon = llatlon0 + llatlon1

#lullat = [90,40,-60]
#lullon = [-140,-100]
#llatlon = [(ullat,ullon) for ullat in lullat for ullon in lullon]
#print llatlon
#sys.exit()

a2out = zeros([180,360],float32)

us = Regrid.UpScale()
for (ullat,ullon) in llatlon:
    lat1 = ullat
    lon0 = ullon

    if lat1==-60:
        lat0 = -90
        lon1 = lon0+60
        nyOrg,nxOrg= 3600, 7200
    else:
        lat0 = lat1-50
        lon1 = lon0+40
        nyOrg,nxOrg= 6000,4800

    hostname = socket.gethostname()
    if   hostname in ['shui']:
        orogDir  = "/data1/hjkim/GTOPO30"
        tankbaseDir='/tank'
    elif hostname in ['well']:
        orogDir  = "/media/disk2/share/data/GTOPO30"
        tankbaseDir='/home/utsumi/mnt/lab_tank'
    else:
        print 'check hostname',hostname
        sys.exit()

    if ullat >0:
        SN = "N"
    else:
        SN = "S"

    if ullon >180:
        WE = "W"
    elif (-180<=ullon)&(ullon<=0):
        WE = "W"
    else:
        WE = "E"

    #orogPath = orogDir + "/E060N40.DEM"
    orogPath = orogDir + "/%s%03d%s%02d.DEM"%(WE, abs(ullon), SN, abs(ullat))
    dmwPath  = orogDir + "/%s%03d%s%02d.DMW"%(WE, abs(ullon), SN, abs(ullat))

    if not os.path.exists(orogPath):
        continue

    print orogPath

    """
    BYTEORDER      M
    LAYOUT       BIL
    NROWS         6000
    NCOLS         4800
    NBANDS        1
    NBITS         16
    BANDROWBYTES         9600
    TOTALROWBYTES        9600
    BANDGAPBYTES         0
    NODATA        -9999
    ULXMAP        60.00416666666667
    ULYMAP        39.99583333333333
    XDIM          0.00833333333333
    YDIM          0.00833333333333
    """

    #ddem[(ullat,ullon)] = flipud(fromfile(orogPath, "int16").byteswap().reshape(-1,nx))
    #ddem[(ullat,ullon)] = ma.masked_equal(flipud(fromfile(orogPath, "int16").byteswap().reshape(-1,nxOrg)),-9999).filled(0)
    a2fin = ma.masked_equal(flipud(fromfile(orogPath, "int16").byteswap().reshape(-1,nxOrg)),-9999).filled(0)

    print a2fin.shape

    dlat = 0.00833333333333
    dlon = 0.00833333333333

    #*** Upscale **********
    LatOrg = arange(lat0+0.5*dlat, lat1-0.5*dlat+0.0001,dlat)
    LonOrg = arange(lon0+0.5*dlon, lon1-0.5*dlon+0.0001,dlon)
    LatUp  = arange(lat0+0.5,lat1-0.5+0.01,1.0)
    LonUp  = arange(lon0+0.5,lon1-0.5+0.01,1.0)

    print LatOrg.shape, a2fin.shape

    us(LatOrg, LonOrg, LatUp, LonUp, globflag=False)

    a2up = us.upscale(a2fin, pergrid=False, miss_in=-9999., miss_out=-9999.) 

    print a2up.shape

    y0 = int(lat0-(-90))
    x0 = int(lon0-(-180))
    nyUp,nxUp = a2up.shape

    a2out[y0:y0+nyUp, x0:x0+nxUp] = a2up

#-- Save ---
outDir = tankbaseDir + '/utsumi/validprof/const'
outPath= outDir + '/orog.meter.sp.one.180x360.npy'
np.save(outPath, a2out)
print outPath

