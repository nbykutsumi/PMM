from numpy import *
import h5py
from datetime import datetime, timedelta
import glob
from math import cos, sin, acos 
import numpy as np
#lori = [0,180]
#lori = [180]
lori = [0]

drev = {0:'016608',180:'016171'}  # ori=0: 016608: 2017-1-30,  ori=180: 016171: 2017-1-2
dictDTime= {0:datetime(2017,1,30), 180:datetime(2017,1,2)}


nx = 221
verGMI = '05'
subverGMI = 'A'
fullverGMI = '%s%s'%(verGMI,subverGMI)
Year = 2017
Mon  = 1
thdist = 15
liy  = range(0,2900+1,100)  # GMI: ny=~2963
#liy  = range(1000,2900+1,100)  # GMI: ny=~2963




#********************************************
def hubeny(lat1, lon1, lat2, lon2):
    #-- for calc ------------
    a  = 6378137
    b  = 6356752.314140
    e2 = 0.00669438002301188
    a_1_e2 = 6335439.32708317
    #------------------------
    latrad1   = lat1 * pi / 180
    latrad2   = lat2 * pi / 180
    lonrad1   = lon1 * pi / 180
    lonrad2   = lon2 * pi / 180
    
    latave    = (latrad1 + latrad2)/2.0
    dlat      = latrad2 - latrad1
    dlon      = lonrad2 - lonrad1
    
    dlondeg   = lon2 - lon1
    if ( abs(dlondeg) > 180.0):
        dlondeg = 180.0 - mod(abs(dlondeg), 180.0)
        dlon    = dlondeg * pi / 180.0
    
    W  = sqrt(1.0 - e2 * sin(latave)**2.0 )
    M  =  a_1_e2 / (W**3.0)
    N  =  a / W
    hubeny  = sqrt( (dlat * M)**2.0 + (dlon * N * cos(latave))**2.0 )
    return hubeny*0.001

def dist_simple(lat1, lon1, lat2, lon2):
    RADEARTH = 6371
    RTD = 57.29578
    DTR = 0.017453
    
    dist = RADEARTH*acos(cos(DTR*lon1-DTR*lon2)*cos(DTR*lat1)*cos(DTR*lat2) + sin(DTR*lat1)*sin(DTR*lat2))
    return dist


#********************************************
for ori in lori:
    rev   = drev[ori]
    DTime = dictDTime[ori]
    Year,Mon,Day = DTime.timetuple()[:3]

    baseDirGMI = '/work/hk01/PMM/NASA/GPM.GMI/1C/V%s'%(verGMI)
    srcDirGMI  = baseDirGMI+ '/%04d/%02d/%02d'%(Year,Mon,Day)
    #ssearchGMI = srcDirGMI + '/1C.GPM.GMI.XCAL2016-C.20170102-S011732-E025005.016171.V05A.HDF5'
    ssearchGMI = srcDirGMI + '/1C.GPM.GMI.*.%s.*.HDF5'%(rev)
    print ssearchGMI
    lsrcPathGMI = glob.glob(ssearchGMI)
    srcPathGMI  = lsrcPathGMI[0]
    print srcPathGMI

    #-- Read HDF file --
    with h5py.File(srcPathGMI) as h:
        a2lat1 = h['S1/Latitude'][:]
        a2lon1 = h['S1/Longitude'][:]
        a2lat2 = h['S2/Latitude'][:]
        a2lon2 = h['S2/Longitude'][:]
    
    #-- dydx file --
    dyPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.dydx/dy.%03d.npy'%(ori)
    dxPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.dydx/dx.%03d.npy'%(ori)

    a2dy = np.load(dyPath)
    a2dx = np.load(dxPath)

    ny,nx = a2dx.shape
    for iy1 in range(ny):
        for ix1 in range(nx):

            dy = a2dy[iy1,ix1]
            dx = a2dx[iy1,ix1]


            if dy !=-9999:
                iy2 = iy1 + dy
                ix2 = ix1 + dx

                lat1 = a2lat1[iy1,ix1]
                lon1 = a2lon1[iy1,ix1]
                lat2 = a2lat2[iy2,ix2]
                lon2 = a2lon2[iy2,ix2]
                km = hubeny(lat1,lon1,lat2,lon2)
            else:
                iy2 = -9999
                ix2 = -9999
                km = -9999.

            print 'iy1,ix1=',iy1,ix1,'dy,dx=',dy,dx,'km=',km


