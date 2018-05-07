import matplotlib
matplotlib.use('Agg')
from numpy import *
from datetime import datetime, timedelta
from gv_fsub import *
import numpy as np
import GPMGV
import myfunc.util as util
import os, sys
import glob
from decimal import Decimal, ROUND_HALF_UP, ROUND_HALF_DOWN
import matplotlib.pyplot as plt



gv = GPMGV.GPMGV()

gv.load_sitelist_reclassified()
satebaseDir = '/home/utsumi/mnt/wellshare/GPMGV/L2A25'
gaugebaseDir= '/work/a01/utsumi/data/GPMGV/2A56'

gvmapbaseDir= '/home/utsumi/mnt/wellshare/GPMGV/GVMAP'

ldomain = ['FLORIDA-STJ']
#ldomain = ['MARYLAND-GSFC']
#ldomain = gv.domains
print ldomain
iYM = [2014,1]
eYM = [2014,1]
#eYM = [2014,9]
lYM = util.ret_lYM(iYM,eYM)
res = 0.02
miss= -9999.

dgName = gv.ret_ddomYM2gName()


def load_gauge(domain, gName, Year,Mon):
    if len(domain.split('-'))==2:
        region, nwName= domain.split('-')
    else:
        region, nwName, tmp = domain.split('-')

    nwCode   = gv.dnwCode[domain]
    gaugeDir = gaugebaseDir + '/%s/%s/%04d'%(region,nwName,Year)

    #-- asc type --
    ssearch  = gaugeDir + '/2A56_%s_%s_%04d%02d_?.asc'%(nwCode, gName, Year, Mon)
    lgaugePath= glob.glob(ssearch)
    
    #-- 2a56 type --
    if len(lgaugePath)==0:
        ssearch  = gaugeDir + '/%s-%s-%04d%02d.2a56'%(nwName, gName, Year, Mon)
        lgaugePath= glob.glob(ssearch)
    if len(lgaugePath) !=0:
        gaugePath = lgaugePath[0] 
    else:
        print 'no file',domain,gName,YM
        sys.exit()

    # load
    sfx = gaugePath.split('.')[-1]
    if sfx =='asc':
        head, aDTime, aPrcp = gv.load_obsfile_asc(gaugePath)
    elif sfx =='2a26':
        head, aDTime, aPrcp = gv.load_obsfile_asc(gaugePath)
    else:
        print 'Check gauge file',gaugePath
        sys.exit()

    return aPrcp



#-- read sitelist_summary ---
listDir  = "/work/a01/utsumi/data/GPMGV/sitelist"
listPath = listDir + "/sitelist_summary.csv"
f = open(listPath, "r"); lines=f.readlines(); f.close()

dBBox    = {}
for line in lines[1:]:
    line = line.strip().split(",")
    region = line[0]
    nwName = line[1]
    nwCode = line[2]
    domain = line[3]
    lllat  = float(line[4])
    urlat  = float(line[5])
    lllon  = float(line[6])
    urlon  = float(line[7])
    sYearDat  = int(line[8])
    sMonDat   = int(line[9])
    eYearDat  = int(line[10])
    eMonDat   = int(line[11])
    lyyyymm   = line[12:]

    key         = domain
    dBBox[key]  = [[lllat,lllon],[urlat,urlon]]


for domain in ldomain:
    #-- set map coordination
    BBox = dBBox[domain]
    [[lat0,lon0],[lat1,lon1]] = BBox

    lllat = float( Decimal(lat0).quantize(Decimal('0.01'), rounding=ROUND_HALF_DOWN) )
    lllon = float( Decimal(lon0).quantize(Decimal('0.01'), rounding=ROUND_HALF_DOWN) )
    urlat = float( Decimal(lat1).quantize(Decimal('0.01'), rounding=ROUND_HALF_UP) )
    urlon = float( Decimal(lon1).quantize(Decimal('0.01'), rounding=ROUND_HALF_UP) )

    a1lat_map = arange(lllat,urlat+res*0.1, res)
    a1lon_map = arange(lllon,urlon+res*0.1, res)

    for YM in lYM:
        Year   = YM[0]
        Mon    = YM[1]


        #-- check satellite overpass time
        sateDir = satebaseDir + '/%s/%04d%02d'%(domain, Year, Mon) 
        if not os.path.exists(sateDir):
            print 'No directory', sateDir
            continue
        ssearch = sateDir + '/prcp.*.npy'
        lsatePath = glob.glob(ssearch)
        lsatePath = sorted(lsatePath)

        ltime = []
        for satePath in lsatePath:
            fileName = os.path.basename(satePath)
            ietime  = fileName.split('.')[1]
            itime   = ietime.split('-')[0]
            etime   = ietime.split('-')[1] 
            iYear,iMon,iDay,iHour,iMnt,iSec = int(itime[:4]), int(itime[4:6]), int(itime[6:8]), int(itime[8:10]), int(itime[10:12]), int(itime[12:14])
            eYear,eMon,eDay,eHour,eMnt,eSec = int(etime[:4]), int(etime[4:6]), int(etime[6:8]), int(etime[8:10]), int(etime[10:12]), int(etime[12:14])
        
            iDTime = datetime(iYear,iMon,iDay,iHour,iMnt)
            eDTime = datetime(eYear,eMon,eDay,eHour,eMnt)
            lDTimeTmp = util.ret_lDTime(iDTime,eDTime,timedelta(seconds=60))
            for DTimeTmp in lDTimeTmp:
                if DTimeTmp not in ltime:
                    ltime.append(DTimeTmp)


        if len(ltime)==0:
            print 'No satellite overpass for',domain,YM
            continue

        #-- add offset to ltime ----
        offset_bef = 4   # minutes
        offset_aft = 31  # minutes


        ltime_tmp =[]
        for DTime0 in ltime:
            for imin in arange(-offset_bef, offset_aft +1):
                DTime = DTime0 + timedelta(minutes=imin)
                if DTime not in ltime_tmp:
                    ltime_tmp.append(DTime)

        ltime = sorted(ltime_tmp)




        #-- load gauge data
        lgName = dgName[domain, Year, Mon]  
        daPrcp = {}
        dLat   = {}
        dLon   = {} 
        for gName in lgName:
            dLat[gName], dLon[gName] = gv.dlatlon[domain,gName]

            daPrcp[gName] = load_gauge(domain, gName, Year, Mon)

        #----
        DTime0 = datetime(Year,Mon,1,0,0)
        for DTime in ltime:
            if DTime.month != Mon:
                continue

            imin = int( (DTime - DTime0).total_seconds() /60 )

            a1lat = []
            a1lon = []
            a1dat = []
            for gName in lgName:
                a1lat.append(dLat[gName])
                a1lon.append(dLon[gName])
                a1dat.append(daPrcp[gName][imin])


            a1lat = array(a1lat)
            a1lon = array(a1lon)
            a1dat = array(a1dat)
   
            print DTime, a1dat.max() 
            if a1dat.max() ==0.0:
                a2gv = array([miss]).astype(float32)
            else:
                a2gv = gv_fsub.point2map(a1dat, a1lon, a1lat, a1lon_map, a1lat_map, 3).T.astype(float32)
    

            gvmapDir  = gvmapbaseDir + '/%s/%04d%02d'%(domain,Year,Mon)
            util.mk_dir(gvmapDir)

            Year = DTime.year
            Mon  = DTime.month
            Day  = DTime.day
            Hour = DTime.hour
            Mnt  = DTime.minute

            gvmapPath = gvmapDir + '/map.%s.%04d%02d%02d.%02d%02d.npy'%(domain,Year,Mon,Day,Hour,Mnt)
            np.save(gvmapPath, a2gv)
            print gvmapPath




