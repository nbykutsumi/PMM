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



iYM = [2014,4]
eYM = [2014,9]
lYM = util.ret_lYM(iYM,eYM)


#ldomain = ['BRAZIL-INP', 'BRAZIL-LBA', 'CALIFORNIA-ERK', 'DARWIN-CSC', 'FLORIDA-KAM-E', 'FLORIDA-KAM-W', 'FLORIDA-KAP', 'FLORIDA-KP2', 'FLORIDA-KSC', 'FLORIDA-NNN', 'FLORIDA-SFL-N', 'FLORIDA-SFL-S', 'FLORIDA-STJ', 'FLORIDA-TFB', 'FRANCE-HyMeX-E', 'FRANCE-HyMeX-W', 'IOWA-IFloodS', 'IOWA-IFloodS_APU_Gauges', 'IOWA-SFB', 'KWAJALEIN-KWA', 'KWAJALEIN-RMI', 'MARYLAND-GSFC', 'MARYLAND-PCMK-N', 'MARYLAND-PCMK-S', 'N.Carolina-IPHEx_Duke', 'N.Carolina-IPHEx_NASA', 'OKLAHOMA-MC3E', 'TEXAS-HAR', 'VIRGINIA-HFD', 'VIRGINIA-NASA-C', 'VIRGINIA-NASA-NE', 'VIRGINIA-NASA-SE', 'VIRGINIA-NASA-W', 'VIRGINIA-NSWD-N', 'VIRGINIA-NSWD-S', 'VIRGINIA-WFF', 'WASHINGTON-OLYMPEX_NASA', 'WASHINGTON-OLYMPEX_STDALN']

ldomain = ['FLORIDA-SFL-S', 'FLORIDA-STJ', 'FLORIDA-TFB', 'FRANCE-HyMeX-E', 'FRANCE-HyMeX-W', 'IOWA-IFloodS', 'IOWA-IFloodS_APU_Gauges', 'IOWA-SFB', 'KWAJALEIN-KWA', 'KWAJALEIN-RMI', 'MARYLAND-GSFC', 'MARYLAND-PCMK-N', 'MARYLAND-PCMK-S', 'N.Carolina-IPHEx_Duke', 'N.Carolina-IPHEx_NASA', 'OKLAHOMA-MC3E', 'TEXAS-HAR', 'VIRGINIA-HFD', 'VIRGINIA-NASA-C', 'VIRGINIA-NASA-NE', 'VIRGINIA-NASA-SE', 'VIRGINIA-NASA-W', 'VIRGINIA-NSWD-N', 'VIRGINIA-NSWD-S', 'VIRGINIA-WFF', 'WASHINGTON-OLYMPEX_NASA', 'WASHINGTON-OLYMPEX_STDALN']

#ldomain = ['FLORIDA-SFL-N']
#ldomain = ['FLORIDA-STJ']
#ldomain = ['N.Carolina-IPHEx_Duke']
#ldomain = gv.domains
ldomain = sorted(ldomain)
print ldomain

res = 0.01
miss= -9999.
thdist = 2.5  # km

gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()
satebaseDir = '/home/utsumi/mnt/wellshare/GPMGV/L2A25'
gaugebaseDir= '/work/a01/utsumi/data/GPMGV/2A56'

gvmapbaseDir= '/home/utsumi/mnt/wellshare/GPMGV/GVMAP'

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
        ssearch  = gaugeDir + '/%s-%s-%04d%02d.2a56'%(nwName.split('_')[0], gName, Year, Mon)
        lgaugePath= glob.glob(ssearch)
    if len(lgaugePath) !=0:
        gaugePath = lgaugePath[0] 
    else:
        print 'no file',domain,gName,YM
        print 'search=',ssearch
        print 'in',gaugeDir
        sys.exit()

    # load
    sfx = gaugePath.split('.')[-1]
    if sfx =='asc':
        head, aDTime, aPrcp = gv.load_obsfile_asc(gaugePath)
    elif sfx =='2a56':
        head, aDTime, aPrcp = gv.load_obsfile_2a56(gaugePath)
    else:
        print 'Check gauge file',gaugePath
        sys.exit()

    return aPrcp


def ret_a2mindist(domain, Year, Mon):
    lgName = dgName[domain, Year, Mon]
    a1lat  = []
    a1lon  = []
    for gName in lgName:
        lat, lon = gv.dlatlon[domain,gName]
        a1lat.append(lat) 
        a1lon.append(lon) 
    a2mindist = gv_fsub.point2mindistmap(a1lon, a1lat, a1lon_map, a1lat_map).T.astype(float32)
    return a2mindist


def dtime_round_Mnt(DTime):
    year = DTime.year
    mon  = DTime.month
    day  = DTime.day
    hour = DTime.hour
    mnt  = DTime.minute
    return datetime(year,mon,day,hour,mnt)



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

        if (domain,Year,Mon) not in dgName.keys():
            print 'no obs',domain,Year,Mon
            continue
        #-- distance from nearest gauge -
        a2mindist = ret_a2mindist(domain, Year, Mon)
        gvmapDir  = gvmapbaseDir + '/%s/%04d%02d'%(domain,Year,Mon)
        util.mk_dir(gvmapDir)
        distPath  = gvmapDir + '/map.mindist.%s.%04d%02d.npy'%(domain,Year,Mon)
        np.save(distPath, a2mindist)
        print 'make minimum distance file'
        print distPath

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

            ##--- test -----
            #if ietime != '20140102172040-20140102172120': continue
            #print fileName
            ##--------------


            itime   = ietime.split('-')[0]
            etime   = ietime.split('-')[1] 
            iYear,iMon,iDay,iHour,iMnt,iSec = int(itime[:4]), int(itime[4:6]), int(itime[6:8]), int(itime[8:10]), int(itime[10:12]), int(itime[12:14])
            eYear,eMon,eDay,eHour,eMnt,eSec = int(etime[:4]), int(etime[4:6]), int(etime[6:8]), int(etime[8:10]), int(etime[10:12]), int(etime[12:14])
        
            iDTime = dtime_round_Mnt(datetime(iYear,iMon,iDay,iHour,iMnt))
            eDTime = dtime_round_Mnt(datetime(eYear,eMon,eDay,eHour,eMnt))
            lDTimeTmp = util.ret_lDTime(iDTime,eDTime,timedelta(seconds=60))
            for DTimeTmp in lDTimeTmp:
                if DTimeTmp not in ltime:
                    ltime.append(DTimeTmp)


        if len(ltime)==0:
            print 'No satellite overpass for',domain,YM
            continue

        #-- add offset to ltime ----
        offset_bef = 9   # minutes
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

            a1latGd = []  # only dood data
            a1lonGd = []
            a1datGd = []

            a1latAll = []  # All data (inc. low quality flagged)
            a1lonAll = []
            a1datAll = []   


            for gName in lgName:
                a1latAll.append(dLat[gName])
                a1lonAll.append(dLon[gName])
                a1datAll.append(abs(daPrcp[gName][imin]))

                if daPrcp[gName][imin] >=0:
                    a1latGd.append(dLat[gName])
                    a1lonGd.append(dLon[gName])
                    a1datGd.append(daPrcp[gName][imin])

            if max(a1datAll) ==0.0:
                a2gvAll = array([miss]).astype(float32)
            else:
                a2gvAll = gv_fsub.point2map_distmask(a1datAll, a1lonAll, a1latAll, a1lon_map, a1lat_map, 3, thdist).T.astype(float32)


            if max(a1datGd) ==0.0:
                a2gvGd = array([miss]).astype(float32)
            else:
                a2gvGd = gv_fsub.point2map_distmask(a1datGd, a1lonGd, a1latGd, a1lon_map, a1lat_map, 3, thdist).T.astype(float32)
 

            ##--- test -----
            #a2gvAll  = ma.masked_less(a2gvAll.astype(int32),0)*0 + DTime.hour*100 + DTime.minute
            #a2gvGd   = ma.masked_less(a2gvGd.astype(int32),0)*0 + DTime.hour*100 + DTime.minute

            #a2gvAll  = a2gvAll.data
            #a2gvGd   = a2gvGd.data
            ##--------------


    

            gvmapDir  = gvmapbaseDir + '/%s/%04d%02d'%(domain,Year,Mon)
            util.mk_dir(gvmapDir)

            Year = DTime.year
            Mon  = DTime.month
            Day  = DTime.day
            Hour = DTime.hour
            Mnt  = DTime.minute

            gvmapPathAll = gvmapDir + '/map.All.%s.%04d%02d%02d.%02d%02d.npy'%(domain,Year,Mon,Day,Hour,Mnt)

            gvmapPathGd = gvmapDir + '/map.Gd.%s.%04d%02d%02d.%02d%02d.npy'%(domain,Year,Mon,Day,Hour,Mnt)

            np.save(gvmapPathAll, a2gvAll)
            np.save(gvmapPathGd,  a2gvGd)
            print gvmapPathAll




