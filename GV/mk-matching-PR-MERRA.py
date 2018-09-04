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


#iYM = [2010,4]
#eYM = [2012,10]
#iYM = [2010,4]
iYM = [2008,4]
eYM = [2014,10]

lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
lYM = lYM[::-1]
print lYM

gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()
dgName  = gv.ret_ddomYM2gName()

merra = MERRA2.M2I1NXASM()

ldomain = gv.domains
#ldomain = [domain for domain in ldomain if not domain in ['FLORIDA-KSC','FLORIDA-SFL-N']]
ldomain = ['FLORIDA-KSC','FLORIDA-SFL-N']

#ldomain = ['BRAZIL-INP', 'BRAZIL-LBA', 'CALIFORNIA-ERK', 'DARWIN-CSC', 'FLORIDA-KAM-E', 'FLORIDA-KAM-W', 'FLORIDA-KAP', 'FLORIDA-KP2', 'FLORIDA-KSC', 'FLORIDA-NNN', 'FLORIDA-SFL-N', 'FLORIDA-SFL-S', 'FLORIDA-STJ', 'FLORIDA-TFB', 'FRANCE-HyMeX-E', 'FRANCE-HyMeX-W', 'IOWA-IFloodS', 'IOWA-IFloodS_APU_Gauges', 'IOWA-SFB', 'KWAJALEIN-KWA', 'KWAJALEIN-RMI', 'MARYLAND-GSFC', 'MARYLAND-PCMK-N', 'MARYLAND-PCMK-S', 'N.Carolina-IPHEx_Duke', 'N.Carolina-IPHEx_NASA', 'OKLAHOMA-MC3E', 'TEXAS-HAR', 'VIRGINIA-HFD', 'VIRGINIA-NASA-C', 'VIRGINIA-NASA-NE', 'VIRGINIA-NASA-SE', 'VIRGINIA-NASA-W', 'VIRGINIA-NSWD-N', 'VIRGINIA-NSWD-S', 'VIRGINIA-WFF', 'WASHINGTON-OLYMPEX_NASA', 'WASHINGTON-OLYMPEX_STDALN']
#ldomain = ['VIRGINIA-HFD', 'VIRGINIA-NASA-C', 'VIRGINIA-NASA-NE', 'VIRGINIA-NASA-SE', 'VIRGINIA-NASA-W', 'VIRGINIA-NSWD-N', 'VIRGINIA-NSWD-S', 'VIRGINIA-WFF', 'WASHINGTON-OLYMPEX_NASA', 'WASHINGTON-OLYMPEX_STDALN']

#ldomain = ['IOWA-SFB']

miss= -9999.

#thdist  = 2.5 # km
thdist  = 5.0 # km

ddtime_1min = timedelta(seconds=60)

satebaseDir = '/home/utsumi/mnt/wellshare/GPMGV/L2A25'
gvmapbaseDir= '/home/utsumi/mnt/wellshare/GPMGV/GVMAP'

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
    a2t  = merra.load_var('T10M' ,DTime)
    a2q  = merra.load_var('QV10M',DTime)
    a2ps = merra.load_var('PS',   DTime)

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
for YM in lYM:
    Year, Mon = YM
    for domain in ldomain:
        if (domain,Year,Mon) not in dgName.keys():
            continue

        # load satellite lat, lon, time
        baseDir = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.L2A25'
        datDir   = baseDir + '/%.1fkm/%s/%04d%02d'%(thdist, domain, Year,Mon)

        if not os.path.exists(datDir):
            print 'no data',domain,YM
            continue
 
        satelatPath  = datDir + '/p_sateLat.npy'
        satelonPath  = datDir + '/p_sateLon.npy'
        dtimePath    = datDir + '/p_dtime.npy'

        a1lat = np.load(satelatPath)
        a1lon = np.load(satelonPath)
        a1dtime= np.load(dtimePath)

        dtime_bef = datetime(1900,1,1,0)

        # empty container
        a1rh = deque([])

        for idat, dtime in enumerate(a1dtime):
            if dtime != dtime_bef:
                dtime_bef = dtime

                #-- load merra2 --
                a2rh= calc_merra_rh_watar(dtime)
                lat = a1lat[idat]
                lon = a1lon[idat]
                y,x = latlon2yxpy_merra(lat,lon)
                rh  = a2rh[y,x]


            a1rh.append(rh)
        
        # save data                 
        a1rh   = array(a1rh).astype(float32)
        rhPath = datDir + '/p_rh.npy'
        np.save(rhPath, a1rh)
        print rhPath


