import matplotlib
matplotlib.use('Agg')
from numpy import *
from datetime import datetime, timedelta
import numpy as np
import myfunc.util as util
import GPMGV
import os, sys, socket
from decimal import Decimal, ROUND_HALF_UP, ROUND_HALF_DOWN
import myfunc.IO.MERRA2 as MERRA2
from collections import deque


#iYM = [2002,4]
#eYM = [2014,10]
#iYM = [1998,4]
#eYM = [2004,10]
iYM = [1998,4]
eYM = [2004,10]

lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
#lYM = lYM[::-1]
print lYM

#gv = GPMGV.GPMGV()
#gv.load_sitelist_reclassified()
#dgName  = gv.ret_ddomYM2gName()
prdName   = 'L2A25'
merra_asm = MERRA2.M2I1NXASM()
merra_slv = MERRA2.M2T1NXSLV()
a1lat_merra= merra_slv.Lat
a1lon_merra= merra_slv.Lon
a2lon_merra, a2lat_merra = meshgrid(a1lon_merra, a1lat_merra)


#ldomain = gv.domains
ldomain = ['FLORIDA-KSC', 'KWAJALEIN-KWA', 'TEXAS-HAR', 'FLORIDA-SFL-N', 'N.Carolina-IPHEx_Duke', 'FLORIDA-STJ', 'N.Carolina-IPHEx_NASA']
#ldomain = ['FLORIDA-KSC']

#ldomain = ['BRAZIL-INP', 'BRAZIL-LBA', 'CALIFORNIA-ERK', 'DARWIN-CSC', 'FLORIDA-KAM-E', 'FLORIDA-KAM-W', 'FLORIDA-KAP', 'FLORIDA-KP2', 'FLORIDA-KSC', 'FLORIDA-NNN', 'FLORIDA-SFL-N', 'FLORIDA-SFL-S', 'FLORIDA-STJ', 'FLORIDA-TFB', 'FRANCE-HyMeX-E', 'FRANCE-HyMeX-W', 'IOWA-IFloodS', 'IOWA-IFloodS_APU_Gauges', 'IOWA-SFB', 'KWAJALEIN-KWA', 'KWAJALEIN-RMI', 'MARYLAND-GSFC', 'MARYLAND-PCMK-N', 'MARYLAND-PCMK-S', 'N.Carolina-IPHEx_Duke', 'N.Carolina-IPHEx_NASA', 'OKLAHOMA-MC3E', 'TEXAS-HAR', 'VIRGINIA-HFD', 'VIRGINIA-NASA-C', 'VIRGINIA-NASA-NE', 'VIRGINIA-NASA-SE', 'VIRGINIA-NASA-W', 'VIRGINIA-NSWD-N', 'VIRGINIA-NSWD-S', 'VIRGINIA-WFF', 'WASHINGTON-OLYMPEX_NASA', 'WASHINGTON-OLYMPEX_STDALN']
#ldomain = ['VIRGINIA-HFD', 'VIRGINIA-NASA-C', 'VIRGINIA-NASA-NE', 'VIRGINIA-NASA-SE', 'VIRGINIA-NASA-W', 'VIRGINIA-NSWD-N', 'VIRGINIA-NSWD-S', 'VIRGINIA-WFF', 'WASHINGTON-OLYMPEX_NASA', 'WASHINGTON-OLYMPEX_STDALN']

#ldomain = ['IOWA-SFB']

miss= -9999.

cbins = 11
ddtime_1min = timedelta(seconds=60)

hostname = socket.gethostname()
if hostname == 'mizu':
    baseDir  = "/home/utsumi/mnt/wellshare/GPMGV/DOM.%s"%(prdName)

elif hostname=='well':
    baseDir = "/media/disk2/share/GPMGV/DOM.%s"%(prdName)


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
    a2t  = merra_asm.load_var('t10m' ,DTime)
    a2q  = merra_asm.load_var('qv10m',DTime)
    a2ps = merra_asm.load_var('ps',   DTime)

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


def calc_cdiff_x(a2dat, miss=None):
    ''' central difference '''
    ny,nx    = a2dat.shape
    a2cdif_x = zeros([ny,nx],float32)

    if miss ==None:
        a2cdif_x[:,1:-1] = a2dat[:,2:] - a2dat[:,:-2]
        a2cdif_x[:,0]    = a2dat[:,1] - a2dat[:,-1]
        a2cdif_x[:,-1]   = a2dat[:,0] - a2dat[:,-2]
    else:
        a2dat = ma.masked_equal(a2dat, miss)
        a2cdif_x[:,1:-1] = (a2dat[:,2:] - a2dat[:,:-2]).filled(miss)
        a2cdif_x[:,0]    = (a2dat[:,1] - a2dat[:,-1]).filled(miss)
        a2cdif_x[:,-1]   = (a2dat[:,0] - a2dat[:,-2]).filled(miss)

    return a2cdif_x

def calc_cdiff_lon(a2dat, miss=None):
    ''' central difference '''
    ny,nx    = a2dat.shape
    a2cdif_x = zeros([ny,nx],float32)

    if miss ==None:
        a2cdif_x[:,1:-1] = a2dat[:,2:] - a2dat[:,:-2]
        a2cdif_x[:,0]    = a2dat[:,1] - a2dat[:,-1]+360
        a2cdif_x[:,-1]   = a2dat[:,0] - a2dat[:,-2]+360
    else:
        a2dat = ma.masked_equal(a2dat, miss)
        a2cdif_x[:,1:-1] = (a2dat[:,2:] - a2dat[:,:-2]).filled(miss)
        a2cdif_x[:,0]    = (a2dat[:,1] - a2dat[:,-1]+360).filled(miss)
        a2cdif_x[:,-1]   = (a2dat[:,0] - a2dat[:,-2]+360).filled(miss)

    return a2cdif_x


def calc_cdiff_y(a2dat,miss=None):
    ''' central difference '''
    ny,nx    = a2dat.shape
    a2cdif_y = zeros([ny,nx],float32)

    if miss==None:
        a2cdif_y[1:-1,:] = a2dat[2:,:] - a2dat[:-2,:]
        a2cdif_y[0,:]    = -9999.
        a2cdif_y[-1,:]   = -9999.
    else:
        a2dat = ma.masked_equal(a2dat,miss)
        a2cdif_y[1:-1,:] = (a2dat[2:,:] - a2dat[:-2,:]).filled(miss)
        a2cdif_y[0,:]    = -9999.
        a2cdif_y[-1,:]   = -9999.

    return a2cdif_y

def calc_merra_div(DTime, lev='850'):

    if lev=='850':
        a2u  = merra_slv.load_var('u850' ,DTime).filled(miss)
        a2v  = merra_slv.load_var('v850' ,DTime).filled(miss)
    elif lev=='10m':
        a2u  = merra_slv.load_var('u10m' ,DTime)
        a2v  = merra_slv.load_var('v10m' ,DTime)

    ny,nx= a2u.shape

    #-- dif --
    a2dlat = calc_cdiff_y(a2lat_merra)
    a2dlon = calc_cdiff_lon(a2lon_merra)
    rad    = 3.14159/180
    a2dy   = a2dlat * rad * 6.37e6
    a2dx   = a2dlon * rad * 6.37e6 * cos(a2lat_merra*rad)

    a2du   = calc_cdiff_x(a2u, miss)
    a2dv   = calc_cdiff_y(a2v, miss)

    a2div  = a2du/a2dx + a2dv/a2dy
    a2div  = ma.masked_where(a2du==-9999., a2div)
    a2div  = ma.masked_where(a2dv==-9999., a2div)

    return a2div.filled(miss)

def latlon2yxpy_merra(lat,lon):
    dlat = merra_slv.dLat
    dlon = merra_slv.dLon
    y = int((lat - (-90-0.5*dlat))/dlat)
    x = int((lon - (-180-0.5*dlon))/dlon)
    return y,x

#*******************************************************
for YM in lYM:
    Year, Mon = YM
    for domain in ldomain:
        #if (domain,Year,Mon) not in dgName.keys():
        #    continue

        # load satellite lat, lon, time
        #baseDir = '/home/utsumi/mnt/wellshare/GPMGV/GLOC.L2A25'
        #baseDir = '/home/utsumi/mnt/wellshare/GPMGV/DOM.L2A25'
        datDir   = baseDir + '/cbin.%d/%s/%04d%02d'%(cbins, domain, Year,Mon)

        if not os.path.exists(datDir):
            print 'no data',domain,YM
            print datDir
            continue
 
        satelatPath  = datDir + '/sateLat.npy'
        satelonPath  = datDir + '/sateLon.npy'
        dtimePath    = datDir + '/dtime.npy'

        a1lat = np.load(satelatPath)
        a1lon = np.load(satelonPath)
        a1dtime= np.load(dtimePath)

        dtime_bef = datetime(1900,1,1,0)

        # empty container
        a1rh     = deque([])
        a1div850 = deque([])
        a1div10m = deque([])

        for idat, dtime in enumerate(a1dtime):
            dtimeTmp = dtime_round_Mnt(dtime)
            if dtimeTmp != dtime_bef:
                dtime_bef = dtimeTmp

                Y,M = dtimeTmp.timetuple()[:2]
                if [Y,M] in lYM:
                    #-- load merra2 --
                    a2rh = calc_merra_rh_watar(dtimeTmp)
                    a2div850= calc_merra_div(dtimeTmp,lev='850')
                    a2div10m= calc_merra_div(dtimeTmp,lev='10m')
    
                    lat = a1lat[idat]
                    lon = a1lon[idat]
                    y,x = latlon2yxpy_merra(lat,lon)
                    rh  = a2rh[y,x]
                    div850 = a2div850[y,x]
                    div10m = a2div10m[y,x]
                else:
                    print 'not',YM
                    rh  = miss
                    div850 = miss
                    div10m = miss
    

            a1rh.append(rh)
            a1div850.append(div850)
            a1div10m.append(div10m)


        # save data                 
        a1rh    = array(a1rh).astype(float32)
        a1div850= array(a1div850).astype(float32)
        rhPath  = datDir + '/rh.npy'
        div850Path = datDir + '/div850.npy'
        div10mPath = datDir + '/div10m.npy'
        np.save(rhPath, a1rh)
        np.save(div850Path, a1div850)
        np.save(div10mPath, a1div10m)
        print rhPath


