from numpy import *
from datetime import datetime, timedelta
import GPMGV
import myfunc.util as util
import os, sys
import glob

gv = GPMGV.GPMGV()

gv.load_sitelist_reclassified()
satebaseDir = '/home/utsumi/mnt/wellshare/GPMGV/L2A25'
gaugebaseDir= '/work/a01/utsumi/data/GPMGV/2A56'
#ldomain = ['FLORIDA-STJ']
ldomain = ['MARYLAND-GSFC']

iYM = [2014,8]
eYM = [2014,8]
lYM = util.ret_lYM(iYM,eYM)


dgName = gv.ret_ddomYM2gName()

print gv.dnwCode.keys()

for domain in ldomain:
    for YM in lYM:
        Year   = YM[0]
        Mon    = YM[1]


        #-- check satellite overpass time
        sateDir = satebaseDir + '/%s/%04d%02d'%(domain, Year, Mon) 
        if not os.path.exists(sateDir): continue
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
        
            iDTime = datetime(iYear,iMon,iDay,iHour,iMnt,iSec)
            eDTime = datetime(eYear,eMon,eDay,eHour,eMnt,eSec)
            dsec   = (eDTime - iDTime).total_seconds()
            mDTime = iDTime + timedelta(seconds= int(dsec*0.5))

            ltime.append(mDTime)
            
        #-- load gauge data
        lgName = dgName[domain, Year, Mon]  
        for gName in lgName:
            region, nwName= domain.split('-')
            nwCode   = gv.dnwCode[domain]
            gaugeDir = gaugebaseDir + '/%s/%s/%04d'%(region,nwName,Year)

            #-- asc type --
            ssearch  = gaugeDir + '/2A56_%s_%s_%04d%02d_?.asc'%(nwCode, gName, Year, Mon)
            lgaugePath= glob.glob(ssearch)
            
            elif len(lgaugePath) >0:
                print 'more than one data'
                print lgaugePath

            #-- 2a56 type --
            else:
                ssearch  = gaugeDir + '/2A56_%s_%s_%04d%02d_?.asc'%(nwCode, gName, Year, Mon)
                lgaugePath= glob.glob(ssearch)
            
            if len(lgaugePath) !=0:
                gaugePath = lgaugePath[0] 
                print gaugePath





#d = gv.ret_ddomYM2gName()
#print d.keys()
#key = d.keys()[0]
#print ''
#print key, d[key]
