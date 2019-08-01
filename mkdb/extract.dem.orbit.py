from numpy import *
import h5py
import sys, os
import numpy as np
from netCDF4 import *
from bisect import bisect_left
from datetime import datetime, timedelta
import myfunc.util as util
import glob
import socket

iDTime = datetime(2017,1,1)
eDTime = datetime(2017,12,31)
lDTimeDay = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

gmibaseDir = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05'

#-- Read GTOPO ----
lullat = [90,40,-10,-60]
lullon = [-180,-140,-100,-60,-20,20,60,100,140]
llatlon = [(ullat,ullon) for ullat in lullat for ullon in lullon]


ddem    = {}
for (ullat,ullon) in llatlon:
    hostname = socket.gethostname()
    if   hostname in ['shui']:
        orogDir  = "/data1/hjkim/GTOPO30"
    elif hostname in ['well']:
        orogDir  = "/media/disk2/share/data/GTOPO30"
    else:
        print 'check hostname',hostname
        sys.exit()


    #ullat = int( (lat - (-60))/50. )*50. + 50 -60.
    #ullon = int( (lon - (-180))/40.)*40. -180.

    if ullat >0:
        SN = "N"
    else:
        SN = "S"

    if ullon >180:
        WE = "W"
    elif (-180<=ullon)&(ullon<0):
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

    #ny = 6000
    nx = 4800
    #ddem[(ullat,ullon)] = flipud(fromfile(orogPath, "int16").byteswap().reshape(-1,nx))
    ddem[(ullat,ullon)] = ma.masked_equal(flipud(fromfile(orogPath, "int16").byteswap().reshape(-1,nx)),-9999).filled(0)

    ## load DMW file
    #f=open(dmwPath, "r"); lines=f.readlines(); f.close()
    #lonmin = float(lines[4])
    #latmax = float(lines[5])
    #lonmax = lonmin + 0.00833333333333*(nx-1)
    #latmin = latmax - 0.00833333333333*(ny-1)

    #dlat = 0.00833333333333
    #dlon = 0.00833333333333
    #para = {}
    #para["ny"] = ny
    #para["nx"] = nx
    #para["lllat"]=latmin
    #para["lllon"]=lonmin
    #para["urlat"]=latmax
    #para["urlon"]=lonmax
    #para["dlat" ]=dlat
    #para["dlon" ]=dlon

    #para["Lat"]  = arange(latmin,latmax+0.5*dlat, dlat)
    #para["Lon"]  = arange(lonmin,lonmax+0.5*dlon, dlon)


dlat = 0.00833333333333
dlon = 0.00833333333333
#------------------
for DTimeDay in lDTimeDay:
    print DTimeDay
    YearDir,MonDir, DayDir = DTimeDay.timetuple()[:3] 
    gmiDir   = gmibaseDir + '/%04d/%02d/%02d'%(YearDir,MonDir,DayDir)

    ssearch  = gmiDir + '/1C.GPM.GMI.*.HDF5'
    lsrcPath = sort(glob.glob(ssearch))
    #print ssearch
    #print lsrcPath
    for srcPath in lsrcPath:
        #srcPath = '/home/utsumi/temp/1C.GPM.GMI.XCAL2016-C.20171206-S141617-E154850.021437.V05A.HDF5'
        print srcPath
        gNum  = srcPath.split('.')[-3]
       
        with h5py.File(srcPath,'r') as h5:
            a2lat = h5['/S1/Latitude'][:]
            a2lon = h5['/S1/Longitude'][:]
        
            #a1year  = h5['S1/ScanTime/Year'        ][:]
            #a1mon   = h5['S1/ScanTime/Month'       ][:]
            #a1day   = h5['S1/ScanTime/DayOfMonth'  ][:]
            #a1hour  = h5['S1/ScanTime/Hour'        ][:]
            #a1mnt   = h5['S1/ScanTime/Minute'      ][:]

        ny,nx  = a2lat.shape
        X,Y    = np.meshgrid(range(nx),range(ny))
        a2elev = zeros([ny,nx],int16)

        for (lat1,lon0) in llatlon:
            if (lat1,lon0) not in ddem.keys():
                #print '(lat1,lon0)'
                continue

            if lat1==-60:
                lat0 = -90
                lon1 = lon0+40
            else:
                lat0 = lat1-50
                lon1 = lon0+40 
 

 
            a2mask1 = ma.masked_greater_equal(a2lat,lat1).mask
            a2mask2 = ma.masked_less(a2lat,lat0).mask
            a2mask3 = ma.masked_greater_equal(a2lon,lon1).mask
            a2mask4 = ma.masked_less(a2lon,lon0).mask
            a2mask  = a2mask1+a2mask2+a2mask3+a2mask4

            a1satx  = ma.masked_where(a2mask, X).compressed()
            a1saty  = ma.masked_where(a2mask, Y).compressed()

            a1satlat = a2lat[a1saty,a1satx]
            a1satlon = a2lon[a1saty,a1satx]

            if len(a1satlat)==0: continue

            #print ''
            #print 'lon'
            #print lon0,lon1
            #print a1satlon.min(), a1satlon.max()

            a1demy  = ((a1satlat - lat0)/dlat).astype(int32)
            a1demx  = ((a1satlon - lon0)/dlon).astype(int32)

            #print 'max demx'
            #print a1demx.max()
            #print 'max satx'
            #print a1satx.max()
            a2elev[a1saty,a1satx] = ddem[(lat1,lon0)][a1demy,a1demx]

        #------------ 
        outbaseDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.gtopo'
        outDir     = outbaseDir + '/%04d/%02d/%02d'%(YearDir,MonDir,DayDir)
        outPath    = outDir + '/%s.%s.npy'%('gtopo', gNum)
        util.mk_dir(outDir)
        np.save(outPath, a2elev.astype(int16))
        print outPath
        
        #sys.exit() 
    
    
