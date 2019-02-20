from numpy import *
import myfunc.util as util
from datetime import datetime, timedelta
import myfunc.IO.GRDC as GRDC
from collections import deque
import numpy as np

grdc = GRDC.GRDC()

iYear = 1998
eYear = 1999

prcp  = 'GPCC'

#grdcDir     = '/home/utsumi/work/temp/GRDC/dat_day'
#simbaseDir  = '/home/utsumi/work/temp/outflw'
#loclistDir  = '/home/utsumi/work/temp/global_15min'
#loclistPath = loclistDir + '/grdc_loc.txt'
#outDir      = '/home/utsumi/work/temp/timeseries'
#figDir      = '/home/utsumi/work/temp/fig'

grdcDir     = '/work/data2/GRDC/dat_day'
simbaseDir  = '/work/hk01/utsumi/PMM/hydro/outflw.trip'
#loclistPath = '/work/hk01/utsumi/PMM/hydro/list/grdc_loc_rev_1998_2009.txt'
outDir      = '/work/hk01/utsumi/PMM/hydro/timeseries.trip'
figDir      = '/work/hk01/utsumi/PMM/hydro/fig'




util.mk_dir(outDir)
util.mk_dir(figDir)
#*** Functions *********************
def load_sim(prcp, stnid, iYear, eYear):
    nYear = eYear - iYear + 1
    lYear = range(iYear,eYear+1)
    a1out = deque([])
    for Year in lYear:
        simDir  = simbaseDir + '/%s/%04d'%(prcp, Year)
        srcPath = simDir + '/flwout.%07d.npy'%(stnid)
        a1tmp   = np.load(srcPath)
        a1out.extend(a1tmp) 

    return array(a1out)

def get_timerange_grdc(srcPath):
    f=open(srcPath,'r'); lines=f.readlines(); f.close()
    for line in lines:
        line = line.split()
        if len(line) ==6:
            if line[2]=='series:':
                iym = map(int, line[3].split('-'))
                eym = map(int, line[5].split('-'))
                break
    return iym, eym
#***********************************

iDTime = datetime(iYear,1,1,0)
eDTime = datetime(eYear,12,31,0)
dDTime = timedelta(days=1)
lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)

#-- Read GRDC location info for TRIP ---
#loclistDir  = '/work/data1/hjkim/cf.mizu'
#loclistPath = loclistDir + '/siteinfo_1.0.csv'
loclistDir  = '/work/hk01/utsumi/PMM/hydro/list'
loclistPath = loclistDir + '/siteinfo_1.0_utm_v02.csv'

f=open(loclistPath,'r'); lines=f.readlines(); f.close()

lstnid = []
lxy    = []
lrivName=[]
lstnName=[]
luparea =[]
for line in lines[1:]:
    line  = line.split(',')
    stnid = int(line[0])
    rivName=line[1]
    rivName= '_'.join(rivName.split())
    stnName=line[2]
    uparea =line[7]
    if uparea < 1000:
        continue
    lstnid.append(stnid)
    lrivName.append(rivName)
    lstnName.append(stnName)

#--------------------
#lstnid = [3629000]

for istn, stnid in enumerate(lstnid):
    grdcPath = grdcDir + '/%07d.day'%(stnid)
    
    try:
        iYMgrdc,eYMgrdc = get_timerange_grdc(grdcPath)
    except IOError:
        continue
    iDTimegrdc = datetime(iYMgrdc[0],iYMgrdc[1],1,0)
    eDTimegrdc = datetime(eYMgrdc[0],eYMgrdc[1],1,0)
     
    if (eDTimegrdc <iDTime) or (eDTime < iDTimegrdc):
        continue
    
   
    # Read GRDC and Simulated data ---- 
    a1grdc = grdc.mk_grdcDailyArray(grdcPath, iDTime, eDTime, dattype='calc')
    #print a1grdc
    #a1sim  = load_sim(prcp, stnid, iYear,eYear) * 1000.
    a1sim  = load_sim(prcp, stnid, iYear,eYear) * 0.001

    a1grdc = ma.masked_less(a1grdc, 0)
    a1sim  = ma.masked_less(a1sim,  0)

    #----------- 
    lout = []
    for i in range(len(a1grdc)):
        obs    = a1grdc[i]
        sim    = a1sim[i]
        DTime  = lDTime[i]
        Year,Mon,Day = DTime.timetuple()[:3]
        stmp   = '%04d,%02d,%02d,%4f,%4f'%(Year,Mon,Day,obs,sim)
        lout.append(stmp)

    #-- Save -- 
    rivName  = lrivName[istn]
    stnName  = lstnName[istn]
    sout     = '%s,%s\n'%(rivName,stnName)
    sout     = sout + 'Year,Mon,Day,Obs,Sim\n'
    sout     = sout + '\n'.join(lout).strip()
    outPath  = outDir + '/flwout.%s.%07d.%s.%s.csv'%(prcp, stnid, rivName, stnName)

    f = open(outPath, 'w'); f.write(sout); f.close()
    print outPath
