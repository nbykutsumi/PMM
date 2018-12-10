import matplotlib
matplotlib.use('Agg')
from numpy import *
from datetime import datetime, timedelta
import numpy as np
import myfunc.util as util
import GPMGV
import os, sys
from decimal import Decimal, ROUND_HALF_UP, ROUND_HALF_DOWN
import myfunc.IO.MERRA2 as MERRA2
from collections import deque
import socket

#iYM = [1998,4]
#eYM = [2004,10]
iYM = [2005,4]
eYM = [2014,10]

lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
lYM = lYM[::-1]
print lYM

gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()
dgName  = gv.ret_ddomYM2gName()

merra = MERRA2.M2I1NXASM()

#ldomain = gv.domains
#ldomain = [domain for domain in ldomain if not domain in ['FLORIDA-KSC','FLORIDA-SFL-N']]
ldomain = ['FLORIDA-KSC', 'KWAJALEIN-KWA', 'TEXAS-HAR', 'FLORIDA-SFL-N', 'N.Carolina-IPHEx_Duke', 'FLORIDA-STJ', 'N.Carolina-IPHEx_NASA']


#ldomain = ['BRAZIL-INP', 'BRAZIL-LBA', 'CALIFORNIA-ERK', 'DARWIN-CSC', 'FLORIDA-KAM-E', 'FLORIDA-KAM-W', 'FLORIDA-KAP', 'FLORIDA-KP2', 'FLORIDA-KSC', 'FLORIDA-NNN', 'FLORIDA-SFL-N', 'FLORIDA-SFL-S', 'FLORIDA-STJ', 'FLORIDA-TFB', 'FRANCE-HyMeX-E', 'FRANCE-HyMeX-W', 'IOWA-IFloodS', 'IOWA-IFloodS_APU_Gauges', 'IOWA-SFB', 'KWAJALEIN-KWA', 'KWAJALEIN-RMI', 'MARYLAND-GSFC', 'MARYLAND-PCMK-N', 'MARYLAND-PCMK-S', 'N.Carolina-IPHEx_Duke', 'N.Carolina-IPHEx_NASA', 'OKLAHOMA-MC3E', 'TEXAS-HAR', 'VIRGINIA-HFD', 'VIRGINIA-NASA-C', 'VIRGINIA-NASA-NE', 'VIRGINIA-NASA-SE', 'VIRGINIA-NASA-W', 'VIRGINIA-NSWD-N', 'VIRGINIA-NSWD-S', 'VIRGINIA-WFF', 'WASHINGTON-OLYMPEX_NASA', 'WASHINGTON-OLYMPEX_STDALN']
#ldomain = ['VIRGINIA-HFD', 'VIRGINIA-NASA-C', 'VIRGINIA-NASA-NE', 'VIRGINIA-NASA-SE', 'VIRGINIA-NASA-W', 'VIRGINIA-NSWD-N', 'VIRGINIA-NSWD-S', 'VIRGINIA-WFF', 'WASHINGTON-OLYMPEX_NASA', 'WASHINGTON-OLYMPEX_STDALN']

#ldomain = ['IOWA-SFB']

miss= -9999.

#thdist  = 2.5 # km
thdist  = 5.0 # km
#cbins   = 11
cbins   = 49

ddtime_1min = timedelta(seconds=60)


hostname  = socket.gethostname()
if   hostname in ['shui']:
    baseDir     = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.L2A25'
elif hostname in ['well']:
    baseDir     = '/media/disk2/share/GPMGV/MATCH.L2A25'

def load_gtopo(lat,lon):
    if   hostname in ['shui']:
        orogDir  = "/data1/hjkim/GTOPO30"
    elif hostname in ['well']:
        orogDir  = "/media/disk2/share/data/GTOPO30"
    else:
        print 'check hostname',hostname
        sys.exit()

    ullat = int( (lat - (-60))/50. )*50. + 50 -60.
    ullon = int( (lon - (-180))/40.)*40. -180.

    if ullat >0:
        SN = "N"
    else:
        SN = "S"

    if ullon >180:
        WE = "W"
    elif (-180<ullon)&(ullon<0):
        WE = "W"
    else:
        WE = "E"
    #orogPath = orogDir + "/E060N40.DEM"
    orogPath = orogDir + "/%s%03d%s%02d.DEM"%(WE, abs(ullon), SN, abs(ullat))
    dmwPath  = orogDir + "/%s%03d%s%02d.DMW"%(WE, abs(ullon), SN, abs(ullat))
    print dmwPath
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
    ny, nx = 6000, 4800
    a = flipud(fromfile(orogPath, "int16").byteswap().reshape(ny,nx))

    # load DMW file
    f=open(dmwPath, "r"); lines=f.readlines(); f.close()
    lonmin = float(lines[4])
    latmax = float(lines[5])
    lonmax = lonmin + 0.00833333333333*(nx-1)
    latmin = latmax - 0.00833333333333*(ny-1)

    dlat = 0.00833333333333
    dlon = 0.00833333333333
    para = {}
    para["ny"] = ny
    para["nx"] = nx
    para["lllat"]=latmin
    para["lllon"]=lonmin
    para["urlat"]=latmax
    para["urlon"]=lonmax
    para["dlat" ]=dlat
    para["dlon" ]=dlon

    para["Lat"]  = arange(latmin,latmax+0.5*dlat, dlat)
    para["Lon"]  = arange(lonmin,lonmax+0.5*dlon, dlon)

    return a,para


def dtime_round_Mnt(DTime):
    year = DTime.year
    mon  = DTime.month
    day  = DTime.day
    hour = DTime.hour
    mnt  = DTime.minute
    return datetime(year,mon,day,hour,mnt)


def ret_elevation(lat,lon):
    a2elev, para = load_gtopo(lat, lon)
    lllat= para['lllat']
    lllon= para['lllon']
    dlat = para['dlat']
    dlon = para['dlon']
    y    = int((lat - lllat)/dlat)
    x    = int((lon - lllon)/dlon)
    return a2elev[y,x]

def calc_merra_rh_watar(DTime):
    #a2t  = merra.load_var('T10M' ,DTime)
    #a2q  = merra.load_var('QV10M',DTime)
    #a2ps = merra.load_var('PS',   DTime)
    a2t  = merra.load_var('t10m' ,DTime)
    a2q  = merra.load_var('qv10m',DTime)
    a2ps = merra.load_var('ps',   DTime)

    # saturation vapor pressure
    # Tetens equation
    # from Sugawara and Kondo 1994, JSHWR
    # 'Errors in Calculation of Satulation Vapor Pressure'
    K = 7.5* (a2t - 273.15)/(a2t-273.15 + 237.3)
    a2esat= 6.1078 * 10**K  # hPa

    # varpor pressure
    #a2e = a2q * a2ps*0.01/(0.622-0.622*a2q+a2q)
    a2e = a2q * a2ps*0.01/0.622
    a2rh= a2e / a2esat
    return a2rh

def latlon2yxpy_merra(lat,lon):
    dlat = merra.dLat
    dlon = merra.dLon
    y = int((lat - (-90-0.5*dlat))/dlat)
    x = int((lon - (-180-0.5*dlon))/dlon)
    return y,x

#*******************************************************
for domain in ldomain:
    readelev = 1

    for YM in lYM:
        Year, Mon = YM

        if (domain,Year,Mon) not in dgName.keys():
            continue

        # load satellite lat, lon, time
        datDir   = baseDir + '/cbin.%d/%.1fkm/%s/%04d%02d'%(cbins, thdist, domain, Year,Mon)

        print datDir
        if not os.path.exists(datDir):
            print 'no data',domain,YM
            continue
 
        satelatPath  = datDir + '/p_sateLat.npy'
        satelonPath  = datDir + '/p_sateLon.npy'

        a1lat = np.load(satelatPath)
        a1lon = np.load(satelonPath)

        if readelev==1:
            readelev=0
            tlat = a1lat[0]
            tlon = a1lon[0]
            a2elev,para = load_gtopo(tlat,tlon)
            lllat= para['lllat']
            lllon= para['lllon']
            dlat = para['dlat']
            dlon = para['dlon']
            nyelev,nxelev =a2elev.shape

        # empty container
        a1selev = deque([])

        for idat in range(len(a1lat)):
            lat = a1lat[idat]
            lon = a1lon[idat]
            y   = int((lat - lllat)/dlat)
            x   = int((lon - lllon)/dlon)
            try:
                elev =a2elev[y,x]
            except:
                tlat = lat
                tlon = lon
                a2elev = load_gtopo(tlat,tlon)
                lllat= para['lllat']
                lllon= para['lllon']
                dlat = para['dlat']
                dlon = para['dlon']
                nyelev,nxelev =a2elev.shape
                elev =a2elev[y,x]

            a1selev.append(elev)

        # save data                 
        a1selev  = array(a1selev).astype(float32)
        elevPath = datDir + '/p_sateElev.npy'
        np.save(elevPath, a1selev)
        print elevPath


