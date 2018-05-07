from numpy import *
from datetime import datetime, timedelta
import numpy as np
import myfunc.util as util
import GPMGV
import os, sys
import glob
from decimal import Decimal, ROUND_HALF_UP, ROUND_HALF_DOWN


gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()

ldomain = gv.domains
ldomain = ['FLORIDA-STJ']

iYM = [2014,1]
eYM = [2014,1]
lYM = util.ret_lYM(iYM, eYM)
res = 0.02
miss= -9999.
offset_bef = 5    # minutes
offset_aft = 30   # minutes

ddtime_1min = timedelta(seconds=60)

satebaseDir = '/home/utsumi/mnt/wellshare/GPMGV/L2A25'
gvmapbaseDir= '/home/utsumi/mnt/wellshare/GPMGV/GVMAP'

def ret_idx_array(Lat, Lon, lllatBnd, lllonBnd, nymap, nxmap):
    Xmap = ( (Lon - lllonBnd)/res ).astype(int32)
    Ymap = ( (Lat - lllatBnd)/res ).astype(int32)
    msk1 = ma.masked_outiside(Xmap, 0, nx-1)
    msk2 = ma.masked_outiside(Ymap, 0, ny-1)

    nysate, nxsate = Lat.shape
    Xsate, Ysate = meshgrid(range(nysate), range(nxsate))

    msk  = ~(msk1 + msk2)

    a1xmap = Xmap[msk] # automatically flatten
    a1ymap = Ymap[msk]
    a1xsate= Xsate[msk]
    a1ysate= Ysate[msk]
 
    return a1xmap, a1ymap, a1xsate, a1ysate    

def ret_dBBoxes(res):
    #-- read sitelist_summary ---
    listDir  = "/work/a01/utsumi/data/GPMGV/sitelist"
    listPath = listDir + "/sitelist_summary.csv"
    f = open(listPath, "r"); lines=f.readlines(); f.close()
    dnynx       = {}
    dBBox       = {}
    dBBoxBnd    = {}
    for line in lines[1:]:
        line = line.strip().split(",")
        domain = line[3]
        lllat  = float(line[4])
        urlat  = float(line[5])
        lllon  = float(line[6])
        urlon  = float(line[7])
   
        #-- pixel centers -- 
        lat0 = float( Decimal(lllat).quantize(Decimal('0.01'), rounding=ROUND_HALF_DOWN) )
        lon0 = float( Decimal(lllon).quantize(Decimal('0.01'), rounding=ROUND_HALF_DOWN) )
        lat1 = float( Decimal(urlat).quantize(Decimal('0.01'), rounding=ROUND_HALF_UP) )
        lon1 = float( Decimal(urlon).quantize(Decimal('0.01'), rounding=ROUND_HALF_UP) )

        a1lat_map = arange(lat0,lat1+res*0.1, res)
        a1lon_map = arange(lon0,lon1+res*0.1, res)

        lllat = a1lat_map.min()
        lllon = a1lon_map.min()
        urlat = a1lat_map.max()
        urlon = a1lon_map.max()

        ny = len(a1lat_map)
        nx = len(a1lon_map)
        #-- pixcel edge --
        lllatBnd = lllat - 0.5*res
        lllonBnd = lllon - 0.5*res
        urlatBnd = urlat + 0.5*res
        urlonBnd = urlon + 0.5*res

        dBBox[domain]    = [[lllat, lllon],[urlat,urlon]]
        dBBoxBnd[domain] = [[lllatBnd,lllonBnd],[urlatBnd,urlonBnd]]
        dnynx[domain]    = [ny,nx] 
    
        return dnynx, dBBox, dBBoxBnd


dnynx, dBBox, dBBoxBnd = ret_dBBoxes(res)

def dtime_round_Mnt(DTime):
    year = DTime.year
    mon  = DTime.month
    day  = DTime.day
    hour = DTime.hour
    mnt  = DTime.minute
    return datetime(year,mon,day,hour,mnt)


for domain in ldomain:
    for YM in lYM:

        Year, Mon = YM
        #-- satellite overpass
        sateDir = satebaseDir + '/%s/%04d%02d'%(domain, Year, Mon)
        if not os.path.exists(sateDir):
            print 'No directory', sateDir
            continue
        ssearch = sateDir + '/prcp.*.npy'
        lsatePath = glob.glob(ssearch)
        lsatePath = sorted(lsatePath)

        #for satePath in lsatePath:
        for satePath in lsatePath[:2]:
            fileName = os.path.basename(satePath)
            gNum    = fileName.split('.')[2]
            ietime  = fileName.split('.')[1]
            #itime   = ietime.split('-')[0]
            #etime   = ietime.split('-')[1]
            print ''
            print ietime            

            satelatPath = sateDir + '/lat.%s.%s.npy'%(ietime, gNum)
            satelonPath = sateDir + '/lon.%s.%s.npy'%(ietime, gNum)
            satetimePath= sateDir + '/dtime.%s.%s.npy'%(ietime, gNum)
            # load sateprcp
            a3sateprcp = np.load(satePath)
            
            a2satelat  = np.load(satelatPath)
            a2satelon  = np.load(satelonPath)
            a1satetime = np.load(satetimePath)

            print a3sateprcp.shape, a2satelat.shape, a2satelon.shape, a1satetime.shape

            print a1satetime.min(), a1satetime.max()
            dtimeA = dtime_round_Mnt(a1satetime.min())
            dtimeB = dtime_round_Mnt(a1satetime.max()) 
            ldtime0 = util.ret_lDTime(dtimeA, dtimeB, ddtime_1min)

            print ldtime0
            for dtime0 in ldtime0:
                print ietime, dtime0
                ldtimeTmp = util.ret_lDTime(dtime0-ddtime_1min*offset_bef, dtime0 + ddtime_1min*offset_aft, ddtime_1min)

                print '-'*50
                print dtime0
                print ldtimeTmp
                print '-'*50

                # mask satellite data
                mskTmp = ma.masked_outside(a1satetime, dtime0, dtime0+ddtime_1min - timedelta(microseconds=1)).mask
                a1timeidx  = ma.masked_where(mskTmp, arange(len(a1satetime)).astype(int32)).compressed()

                print ietime

                a3sateprcpTmp = a3sateprcp[a1timeidx,:,:]
                a2satelatTmp  = a2satelat[a1timeidx,:]
                a2satelonTmp  = a2satelon[a1timeidx,:]

                print a3sateprcpTmp.shape, a2satelatTmp.shape, a2satelonTmp.shape
                for dtimeTmp in ldtimeTmp:
                    dtimeGV = dtimeTmp + ddtime_1min
                    YearGV = dtimeGV.year
                    MonGV  = dtimeGV.month
                    DayGV  = dtimeGV.day
                    HourGV = dtimeGV.hour
                    MntGV  = dtimeGV.minute

                    # load gvmap
                    gvmapDir = gvmapbaseDir + '/%s/%04d%02d'%(domain,YearGV,MonGV)

                    gvmapPath = gvmapDir + '/map.%s.%04d%02d%02d.%02d%02d.npy'%(domain,YearGV,MonGV,DayGV,HourGV,MntGV)
                    a2gvmap = np.load(gvmapPath)

                     
 
