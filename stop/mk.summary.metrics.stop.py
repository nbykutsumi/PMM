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

#lsurftype = ['ocean','vege','coast','snow']
lsurftype = ['ocean','seaice','vege','snow','swater','coast','siedge']

#lsurftype = ['ocean']


region = 'GLB'
#season = 'ALL'
season = 'JJA'
DB_MAXREC = 10000
DB_MINREC = 1000

lexpr = ['NScmb.glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC),
        'ml.best01-HTQZ',
        #'ml-2017.best01-HTQZ',
        #'ml-2017.best01cr-HTQZ'
        ]

dexprshort = {}
for expr in lexpr:
    tmp = expr.split('.')[0]
    if tmp in ['ml','ml-2017']:
        exprshort = expr
    else:
        exprshort = 'org'
    dexprshort[expr] = exprshort
#*************************************************
lvar = ['cc','nbias','rmse','slope','intercept']
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

nexpr = len(lexpr)
llabel1 = [''] + ['cc']*nexpr + ['nbias']*nexpr + ['rmse']*nexpr + ['slope']*nexpr + ['intercept']*nexpr
llabel2 = [''] + [dexprshort[expr] for expr in lexpr]*5
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
