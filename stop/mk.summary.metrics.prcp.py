import numpy as np
from numpy import *
import h5py
from bisect import bisect_left
import sys, os, glob, socket
import myfunc.util as util
import calendar
import pickle
import random

myhost = socket.gethostname()
if myhost == 'shui':
    listDir    = '/work/hk01/utsumi/PMM/US/obtlist'
    workbaseDir= '/work'
    tankbaseDir= '/tank'
    epcbaseDir = '/tank/utsumi/PMM/retepc'
    gprofbaseDir = '/work/hk01/PMM/NASA/GPM.GMI/2A/V05'

    figDir   = '/home/utsumi/temp/stop'

elif myhost == 'well':
    listDir    = '/home/utsumi/mnt/lab_tank/utsumi/PMM/US/obtlist'
    workbaseDir= '/home/utsumi/mnt/lab_work'
    tankbaseDir= '/home/utsumi/mnt/lab_tank'
    epcbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc'
    gprofbaseDir = '/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.GMI/2A/V05'
    figDir   = '/home/utsumi/temp/stop'

else:
    print 'check hostname',myhost
    sys.exit()

lsurftype = ['ocean','vegetation','coast','snow']
#lsurftype = ['ocean']

dsurflabel={ 'ocean':'Class1 (Ocean)'
            ,'vegetation':'Classes 3-7 (Vegetation)'
            ,'snow':'Classes 8-11 (Snow)'
            ,'coast':'Class 13 (Coast)'
            }

region = 'GLB'
#season = 'JJA'
season = 6
DB_MAXREC = 10000
DB_MINREC = 1000

lexpr = ['glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC),
        #'glb.stop-wgt-obs-01.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC),
        #'glb.stop-rng-obs-01.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC),

        'glb.stop-wgt-obs-01.minrec%d.maxrec%d-th0.0mm'%(DB_MINREC,DB_MAXREC),
        'glb.stop-rng-obs-01.minrec%d.maxrec%d-th0.0mm'%(DB_MINREC,DB_MAXREC),
        'glb.stop-wgt-ret-01.minrec%d.maxrec%d-th1.0mm'%(DB_MINREC,DB_MAXREC),
        'glb.stop-wgt-ret-01.minrec%d.maxrec%d-th5.0mm'%(DB_MINREC,DB_MAXREC),
        'glb.stop-wgt-ret-01.minrec%d.maxrec%d-th7.0mm'%(DB_MINREC,DB_MAXREC),

        ]

dexprshort = {}
for expr in lexpr:
    tmp = expr.split('.')[1]
    if tmp[:4]=='stop':
        if expr[-2:] =='mm':
            thpr = expr.split('-')[-1]
            exprshort = '-'.join(tmp.split('-')[1:3])+ '-' + thpr
        else:
            exprshort = '-'.join(tmp.split('-')[1:3])
    else:
        exprshort = 'org'
    dexprshort[expr] = exprshort
#*************************************************
lvar = ['num','cc','nbias','rmse']
dat = {}
for expr in lexpr:
    for surftype in lsurftype:
        #*******************************
        # Metrics
        #-------------------------------
        stampOut = '%s.%s.%s.%s'%(expr, region, surftype, season)
        iPath = figDir + '/metrics.prcp.%s.csv'%(stampOut)
        f=open(iPath,'r'); lines =f.readlines(); f.close()

        line = lines[1].strip().split(',')
        for i,var in enumerate(lvar):
            dat[expr,surftype,var] = line[i] 

nexpr = len(lexpr)
llabel1 = [''] + ['num']*nexpr + ['cc']*nexpr + ['nbias']*nexpr + ['rmse']*nexpr
llabel2 = [''] + [dexprshort[expr] for expr in lexpr]*4

lout = []
lout.append(llabel1)
lout.append(llabel2)
for surftype in lsurftype:
    line = [surftype]
    for var in lvar:
        for expr in lexpr:
            key = (expr, surftype,var)
            line.append(dat[key])

    lout.append(line)

sout = util.list2csv(lout)
oPath = figDir + '/summary.prcp.%s.%s.csv'%(region,season)
f=open(oPath,'w'); f.write(sout); f.close()
print oPath
