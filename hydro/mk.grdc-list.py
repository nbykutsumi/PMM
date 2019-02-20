from numpy import *
from datetime import datetime, timedelta
import glob
import myfunc.IO.GRDC
import numpy as np
import myfunc.util as util
import sys

grdcDir = '/work/data2/GRDC/dat_day'
#iDTime = datetime(1900,1,1,0)
#eDTime = datetime(2019,12,31,0)

iDTime = datetime(1998,1,1,0)
eDTime = datetime(2009,12,31,0)

res    = '1deg'
res    = '15min'

if   res=='1deg':
    dlatlon=1.0
elif res=='15min':
    dlatlon=0.25
else:
    print 'check res',res
    sys.exit()


#iDTime = datetime(2000,1,1,0)
#eDTime = datetime(2009,12,31,0)

#-- Functions --
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

def readHeaderDaily(srcPath):
    f=open(srcPath,'r'); lines=f.readlines(); f.close()
    d = {}
    #-- River ------
    #d['River']   = lines[9].split(':')[1].strip()
    river = lines[9].split(':')[1].strip()
    d['River'] = '_'.join(river.split())

    #-- Station ----
    #d['Station'] = lines[10].split(':')[1].strip()
    station = lines[10].split(':')[1].strip()
    d['Station'] = '_'.join(station.split())

    #-- Lat & Lon--
    d['Lat']  = float(lines[12].split()[4])
    d['Lon']  = float(lines[13].split()[4])

    #-- Catchment area [km2] --
    d['uparea']=float(lines[14].split()[4])    

    #-- Altitude [m.a.s.l] ---
    d['Altitude'] = float(lines[15].split()[3])   

    #-- timerange --
    line = lines[23].split()
    iym  = map(int, line[3].split('-'))
    eym  = map(int, line[5].split('-'))
    d['timerange'] = iym, eym

    return d

#---------------
#-- Read uparea --
upareaDir  = '/work/hk01/hjkim/CaMa-Flood/CaMa-Flood_v3.6.2_20140909/map/global_15min'
upareaPath = upareaDir + '/uparea.bin'
a2uparea   = fromfile(upareaPath, float32).reshape(720,1440)

#-- Read basin number --
basinDir   = '/work/hk01/hjkim/CaMa-Flood/CaMa-Flood_v3.6.2_20140909/map/global_15min'
basinPath  = basinDir + '/basin.bin'
a2basin    = fromfile(basinPath, int32).reshape(720,1440)

#-----------------
ssearch = grdcDir + '/*.day'
lsrcPath= sort(glob.glob(ssearch))

sout = ''
n = 0
for srcPath in lsrcPath:
    iym,eym = get_timerange_grdc(srcPath)
    iDTimeTmp = datetime(iym[0],iym[1],1,0)
    eDTimeTmp = datetime(eym[0],eym[1],1,0)

    stnid   = int(srcPath.split('/')[-1].split('.')[0])

    #if stnid !=1134900: continue

    if (eDTimeTmp<iDTime)or(eDTime<iDTimeTmp):
        continue

    n = n+1 
    head = readHeaderDaily(srcPath)
    lat  = head['Lat']
    lon  = head['Lon']

    y0   = int((90-lat)/0.25)
    x0   = int((lon+180)/0.25)
    uparea0 = head['uparea']
    ldydx = [[dy,dx] for dy in [-2,-1,0,1,2] for dx in [-2,-1,0,1,2]]
    lbias = []
    for dy,dx in ldydx:
        y=y0+dy
        x=x0+dx
        uparea1 = a2uparea[y,x]*10**(-6)
        bias = abs(uparea1-uparea0)/uparea0
        #print stnid, n,x0+dx, y0+dy, dx,dy,'%.2f'%(bias)
        lbias.append(bias)
    imin    = np.argmin(lbias)
    dy,dx   = ldydx[imin]
    y,x     = y0+dy, x0+dx
    uparea1 = a2uparea[y,x]*1e-6
    bias_area= (uparea1-uparea0)/uparea0

    iym,eym = head['timerange']
    rivnum  = a2basin[y,x]
    iym     = '%04d %02d'%(iym[0],iym[1])
    eym     = '%04d %02d'%(eym[0],eym[1])

    rivName = head['River'].ljust(37)
    stnName = head['Station'].ljust(36)

    stmp = '%07d %s %s x'%(stnid, rivName, stnName)
    stmp = stmp + ' %7.2f %6.2f %4d %4d %7d %7d %5.2f'%(lon, lat, x, y, uparea0, uparea1, bias_area)
    stmp = stmp + ' %7d %s %s'%(rivnum, iym, eym)
    sout = sout + '\n' + stmp
    print stmp

#-- Save --
listDir = '/work/hk01/utsumi/PMM/hydro/list'
util.mk_dir(listDir)
listPath= listDir + '/grdc_loc_rev_%04d_%04d.txt'%(iDTime.year, eDTime.year)
f=open(listPath,'w'); f.write(sout.strip()); f.close()
print listPath
