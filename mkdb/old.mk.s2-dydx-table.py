from numpy import *
import h5py
from datetime import datetime, timedelta
import glob
from math import cos, sin, acos 
#lori = [0,180]
lori = [180]

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
    h = h5py.File(srcPathGMI)
    a2lat1 = h['S1/Latitude'][:]
    a2lon1 = h['S1/Longitude'][:]
    a2lat2 = h['S2/Latitude'][:]
    a2lon2 = h['S2/Longitude'][:]


    ldy2 = arange(0,20).astype(int32)
    ldy3 = arange(0,5).astype(int32)
    if ori==180:
        ldy2 = -ldy2
        ldy3 = -ldy3
    else:
        ldy2 = ldy2
        ldy3 = ldy3

    for iy1 in liy:
        for ix1 in range(0,nx):
        #for ix1 in range(103,105):

            lat1 = a2lat1[iy1,ix1]
            lon1 = a2lon1[iy1,ix1]

            found = 0
            distmin = 99999
            dyopt   = -9999
            dxopt   = -9999

            for dy2 in ldy2:
                iy2 = iy1 + dy2
                #print iy1, ix1, iy2
                for ix2 in range(0,nx):
                    dx2   = ix2-ix1
                    try:
                        lat2 = a2lat2[iy2,ix2]
                        lon2 = a2lon2[iy2,ix2]
                    except IndexError:
                        continue
                    km = hubeny(lat1,lon1,lat2,lon2)
                    #km= dist_simple(lat1,lon1,lat2,lon2)
                    #print iy1,ix1,'--',iy2,ix2,'lat,lon=', lat1,lon1,'--',lat2,lon2, km
                    if km < thdist:
                        kmmin = km
                        dyopt = iy2 - iy1
                        dxopt = ix2 - ix1
                        for dx3 in range(-5,5):
                            for dy3 in ldy3:
                                iy3 = iy2 + dy3
                                ix3 = ix2 + dx3
                                try:
                                    lat3 = a2lat2[iy3,ix3]
                                    lon3 = a2lon2[iy3,ix3]
                                except IndexError:
                                    continue
                                
                                km = hubeny(lat1,lon1,lat3,lon3)
                                if km <kmmin:
                                    kmmin = km
                                    dyopt = iy2 - iy1
                                    dxopt = ix2 - ix1
                                    

                        fount = 1
                        print '**********************************************'
                        print 'iy1, ix1=',iy1, ix1, '  dy, dx=', dyopt, dxopt, '  dist=',kmmin
                        print '**********************************************'
                        break
                if found ==1:
                    break
