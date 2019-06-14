import numpy as np
import pylab as pl
import matplotlib.gridspec as gridspec
from glob import glob
import numpy.ma as ma
import sys,os
from datetime import datetime, timedelta
import myfunc.util as util
import calendar
from collections import deque

#***********************************************************
# Functions
#***********************************************************

def read_Tc(lDTime=None, ldydx=None, isurf=None):
    a2tc = deque([])
    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]
        a2tcTmp = None
        for idydx,(dy,dx) in enumerate(ldydx):
            srcDir = '/work/hk01/utsumi/PMM/stop/data/Tc/%04d/%02d/%02d'%(Year,Mon,Day)
            #srcDir = '/mnt/j/PMM/stop/data/Tc/%04d/%02d/%02d'%(Year,Mon,Day)
            srcPath1=srcDir + '/Tc1.%ddy.%ddx.%02dsurf.npy'%(dy,dx,isurf)
            srcPath2=srcDir + '/Tc2.%ddy.%ddx.%02dsurf.npy'%(dy,dx,isurf)
            if not os.path.exists(srcPath1):
                print 'No file:',srcPath1
                continue
            atc1 = np.load(srcPath1)
            atc2 = np.load(srcPath2)
            atc  = np.c_[atc1, atc2]

            if atc.shape[0]==0:  # new
                continue

            if a2tcTmp is None:  # new
                a2tcTmp = atc
            else:
                a2tcTmp = np.c_[a2tcTmp, atc]


        if a2tcTmp is None:
            continue
        else:
            a2tcTmp = np.array(a2tcTmp)

        a2tc.extend(a2tcTmp)
    return np.array(a2tc)

def read_var_collect(varName=None, lDTime=None, ldydx=None, isurf=None):
    a2var = deque([])
    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]
        a2varTmp = None
        for idydx,(dy,dx) in enumerate(ldydx):
            srcDir = '/work/hk01/utsumi/PMM/stop/data/%s/%04d/%02d/%02d'%(varName,Year,Mon,Day)
            #srcDir = '/mnt/j/PMM/stop/data/%s/%04d/%02d/%02d'%(varName,Year,Mon,Day)
            srcPath=srcDir + '/%s.%ddy.%ddx.%02dsurf.npy'%(varName,dy,dx,isurf)
            if not os.path.exists(srcPath): continue
            avar = np.load(srcPath)

            if avar.shape[0]==0:  # new
                continue

            if a2varTmp is None:
                a2varTmp = avar
            else:
                a2varTmp = np.c_[a2varTmp, avar]

        if a2varTmp is None:
            continue
        else:
            a2varTmp = np.array(a2varTmp)
        #**********************
        a2var.extend(a2varTmp)
    return np.array(a2var)

print 'Define functions'
#***********************************************************
# Main loop
#***********************************************************
Year   = 2017
lMon   = range(1,12+1)
#lMon   = [6]
#lisurf = [8]
lisurf = range(1,14+1)
#lisurf = [8]
#llatminmax = [[-90,90],[-90,-20],[-30,30],[20,90]]
llatminmax = [[20,90]]
#llatminmax = [[-90,90]]

dnum   = {}
for isurf in lisurf:
    for Mon in lMon:
        ldy   = [-1,0,1]
        ldx   = [-3,-2,-1,0,1,2,3]
        #ldx   = [-2,-1,0,1,2]
        imid  = int((len(ldy)*len(ldx)-1)/2)
        ldydx = [[dy,dx] for dy in ldy for dx in ldx]
        eDay  = calendar.monthrange(Year,Mon)[1]
        lDTime = util.ret_lDTime(datetime(Year,Mon,1), datetime(Year,Mon,eDay),timedelta(days=1))

        
        #****************************************************
        # Read data
        #****************************************************
        trainTc   = read_Tc(lDTime, ldydx, isurf)
        trainStop = read_var_collect('stop', lDTime, [[0,0]], isurf)
        trainLat  = read_var_collect('Latitude', lDTime, [[0,0]], isurf)
        

        #****************************************************
        # Screen invalid Tc
        #****************************************************
        a1flagtc = ma.masked_inside(trainTc, 50, 350).all(axis=1).mask
        trainTc  = trainTc  [a1flagtc]
        trainStop= trainStop[a1flagtc]
        trainLat = trainLat [a1flagtc]


        for latmin,latmax in llatminmax:
            print isurf,Mon,latmin,latmax
            #****************************************************
            # Screen latitude
            #****************************************************
            index_keep = []
            for i in range(trainTc.shape[0]):
                lat = trainLat[i]
                if (latmin<=lat)and(lat<=latmax):
                    index_keep.append(i)
           
            if len(index_keep)==0:
                dnum[isurf, Mon, latmin, latmax] = 0
            else:
                tmpTc   = trainTc[index_keep]
                tmpStop = trainStop[index_keep]
                tmpLat  = trainLat[index_keep]
                dnum[isurf, Mon, latmin, latmax] = tmpTc.shape[0]


for (latmin,latmax) in llatminmax:
    lout   = []
    llabel = [''] + lMon
    lout.append(llabel)

    for isurf in lisurf:
        ltmp = [dnum[isurf, Mon, latmin, latmax] for Mon in lMon] 
        ltmp = [isurf] + ltmp
        lout.append(ltmp)

    sout = util.list2csv(lout)
    csvPath = '/home/utsumi/temp/stop/num.min.%d.max.%d.csv'%(latmin,latmax)
    #csvPath = '/mnt/c/ubuntu/temp/num.min.%d.max.%d.csv'%(latmin,latmax)
    f=open(csvPath, 'w'); f.write(sout); f.close()
    
    print csvPath
