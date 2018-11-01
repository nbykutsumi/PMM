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
iYM = [2001,4]
eYM = [2004,10]
#eYM = [2005,3]

lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
#lYM = lYM[::-1]
print lYM

gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()
dgName  = gv.ret_ddomYM2gName()

merra_asm = MERRA2.M2I1NXASM()
merra_slv = MERRA2.M2T1NXSLV()
a1lat_merra= merra_slv.Lat
a1lon_merra= merra_slv.Lon
a2lon_merra, a2lat_merra = meshgrid(a1lon_merra, a1lat_merra)




ldomain = gv.domains
#ldomain = ['FLORIDA-KSC', 'KWAJALEIN-KWA', 'TEXAS-HAR', 'FLORIDA-SFL-N', 'N.Carolina-IPHEx_Duke', 'FLORIDA-STJ', 'N.Carolina-IPHEx_NASA']
ldomain = ['N.Carolina-IPHEx_Duke']

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

def calc_merra_div(DTime):

    a2u  = merra_slv.load_var('u850' ,DTime).filled(miss)
    a2v  = merra_slv.load_var('v850' ,DTime).filled(miss)

    ny,nx= a2u.shape

    #-- div --
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


    iy, ix = 100, 160
    du = a2u[iy,ix+1] - a2u[iy,ix-1]
    dv = a2v[iy+1,ix] - a2v[iy-1,ix]
   
    dlat = a2lat_merra[iy+1,ix] - a2lat_merra[iy-1,ix]
    dlon = a2lon_merra[iy,ix+1] - a2lon_merra[iy,ix-1]
    lat  = a2lat_merra[iy,ix]
    dy   = dlat*rad * 6.37e6
    dx   = dlon*rad * 6.37e6 * cos(lat*rad)

    div  = du/dx + dv/dy

    print a2dy[iy,ix], dy
    print a2dx[iy,ix], dx
    print a2div[iy,ix], div 


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
        baseDir = '/home/utsumi/mnt/wellshare/GPMGV/GLOC.L2A25'
        datDir   = baseDir + '/%.1fkm/%s/%04d%02d'%(thdist, domain, Year,Mon)

        if not os.path.exists(datDir):
            print 'no data',domain,YM
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

        for idat, dtime in enumerate(a1dtime):
            dtimeTmp = dtime_round_Mnt(dtime)
            if dtimeTmp != dtime_bef:
                dtime_bef = dtimeTmp

                Y,M = dtimeTmp.timetuple()[:2]
                if [Y,M] in lYM:
                    #-- load merra2 --
                    a2rh = calc_merra_rh_watar(dtimeTmp)
                    a2div= calc_merra_div(dtimeTmp)
    
                    lat = a1lat[idat]
                    lon = a1lon[idat]
                    y,x = latlon2yxpy_merra(lat,lon)
                    rh  = a2rh[y,x]
                    div = a2div[y,x]



                    sys.exit()

                else:
                    print 'not',YM
                    rh  = miss
                    div = miss
    

            a1rh.append(rh)
            a1div850.append(div)
        
        # save data                 
        a1rh    = array(a1rh).astype(float32)
        a1div850= array(a1div850).astype(float32)
        rhPath  = datDir + '/rh.npy'
        div850Path = datDir + '/div850.npy'


