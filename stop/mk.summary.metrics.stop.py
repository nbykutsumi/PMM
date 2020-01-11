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
season = 8
DB_MAXREC = 10000
DB_MINREC = 1000

lexpr = ['NScmb.glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC),
        'ml.best01-HTQZ'
        ]

dexprshort = {}
for expr in lexpr:
    tmp = expr.split('.')[0]
    if tmp=='ml':
        exprshort = 'ml'
    else:
        exprshort = 'org'
    dexprshort[expr] = exprshort
#*************************************************
lvar = ['cc','nbias','rmse']
dat = {}
for expr in lexpr:
    for surftype in lsurftype:
        #*******************************
        # Metrics
        #-------------------------------
        stampOut = '%s.%s.%s.%s'%(expr, region, surftype, season)
        iPath = figDir + '/metrics.stop.%s.csv'%(stampOut)
        f=open(iPath,'r'); lines =f.readlines(); f.close()

        line = lines[1].strip().split(',')
        for i,var in enumerate(lvar):
            dat[expr,surftype,var] = line[i] 

llabel1 = [''] + ['cc']*2 + ['nbias']*2 + ['rmse']*2
llabel2 = [''] + [dexprshort[expr] for expr in lexpr]*3
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
oPath = figDir + '/summary.stop.%s.%s.csv'%(region,season)
f=open(oPath,'w'); f.write(sout); f.close()
print oPath
