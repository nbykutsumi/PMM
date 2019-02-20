import matplotlib
matplotlib.use('Agg')
from numpy import *
import myfunc.util as util
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from datetime import datetime, timedelta
import myfunc.IO.GRDC as GRDC
import numpy as np
from collections import deque

grdc = GRDC.GRDC()

iYear = 2000
eYear = 2009

lprcp  =['JRA55']
lmetName=['CC','NSE','RBias']
#calc  = True
calc  = False

grdcDir     = '/home/utsumi/work/temp/GRDC/dat_day'
simbaseDir  = '/home/utsumi/work/temp/outflw'
loclistDir  = '/home/utsumi/work/temp/global_15min'
loclistPath = loclistDir + '/grdc_loc.txt'
outDir      = '/home/utsumi/work/temp/metric'
figDir      = '/home/utsumi/work/temp' + '/fig'
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

#--- Read loc list --
f=open(loclistPath,'r'); lines=f.readlines(); f.close()
lstnid = []
dyx    = {}
dlatlon= {}
for line in lines:
    line  = line.split()
    stnid = int(line[0])
    x       = int(line[6])  # 0,1,2..
    y       = int(line[7])  # 0,1,2..
    lon     = float(line[4])
    lat     = float(line[5])
    dyx[stnid]     = [y,x]
    dlatlon[stnid] = [lat,lon]

    lstnid.append(stnid)
#--------------------
#lstnid = lstnid[:5]  # test

for metName in lmetName:
    for prcp in lprcp:
        if calc != True: break
    
        lout = []
        for stnid in lstnid:
            grdcPath = grdcDir + '/%07d.day'%(stnid)
            
            iYMgrdc,eYMgrdc = get_timerange_grdc(grdcPath)
            iDTimegrdc = datetime(iYMgrdc[0],iYMgrdc[1],1,0)
            eDTimegrdc = datetime(eYMgrdc[0],eYMgrdc[1],1,0)
            
            if (eDTimegrdc <iDTime) or (eDTime < iDTimegrdc):
                continue
            
            # Read GRDC and Simulated data ---- 
            a1grdc = grdc.mk_grdcDailyArray(grdcPath, iDTime, eDTime, dattype='calc')
            #print a1grdc
            #a1sim  = load_sim(prcp, stnid, iYear,eYear) * 1000.
            a1sim  = load_sim(prcp, stnid, iYear,eYear)
        
            a1grdc = ma.masked_less(a1grdc, 0)
            a1sim  = ma.masked_less(a1sim,  0)
            a1grdc = ma.masked_invalid(a1grdc)
            a1sim  = ma.masked_invalid(a1sim)
            # Metrics ----------- 
            if metName == 'CC':
                var = np.ma.corrcoef(a1grdc, a1sim, allow_masked=True)[0,1]
            elif metName =='NSE':
                var = ((a1grdc - a1sim)**2).sum() \
                    / ((a1grdc - a1grdc.mean())**2).sum()
                var = 1.0 - var

            elif metName =='RBias':
                var = (a1sim.mean()/a1grdc.mean())

            stmp   = '%07d  %.4f'%(stnid, var)
            lout.append(stmp)
        
        #-- Save metric -- 
        sout     = '\n'.join(lout).strip()
        corrPath = outDir + '/%s.%s.txt'%(metName, prcp)
        f = open(corrPath, 'w'); f.write(sout); f.close()
        print corrPath
        #for i in range(len(a1grdc)):
        #    print a1grdc[i], a1sim[i]
#***************************************************
# Draw map 
#***************************************************
BBox = [[-60,-180],[90,180]]
[[lllat,lllon],[urlat,urlon]] = BBox

for metName in lmetName:
    for prcp in lprcp:
        #-- Read CC file --
        corrPath = outDir + '/%s.%s.txt'%(metName, prcp)
        f = open(corrPath,'r'); lines = f.readlines(); f.close()
        lvar  = []
        llat  = []
        llon  = []
        for line in lines:
            line = line.split()
            stnid = int(line[0])
            var   = float(line[1])
            lat,lon  = dlatlon[stnid]
            lvar.append(var)
            llat.append(lat)
            llon.append(lon)
        
        #-- Draw map -
        if   metName=='CC':
            vmin,vmax = [0,1]
            cmap      = 'jet'
        elif metName=='NSE':
            vmin,vmax = [0,1]
            cmap      = 'jet'
        elif metName=='RBias':
            vmin,vmax = [0,2] 
            cmap      = 'Spectral_r'

        fig  = plt.figure(figsize=(8,4))
        fig  = plt.figure()
        ax   = fig.add_subplot(1,1,1)
        M    = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)
        im  = M.scatter(llon, llat, c=lvar, s=10, marker='s', cmap=cmap, vmin=vmin, vmax=vmax)
        plt.colorbar(im, orientation='horizontal')
        
        meridians = arange(-180,180,30)
        parallels = arange(-90,90,30)
        M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=8, linewidth=0.7)
        M.drawparallels(parallels, labels=[1,0,0,0], fontsize=8, linewidth=0.7)
        M.drawcoastlines()
        stitle  = '%s %s'%(metName, prcp)
        plt.title(stitle)
        #-- save --
        figPath = figDir + '/map.%s.%s.png'%(metName, prcp)
        plt.savefig(figPath)
        print figPath


