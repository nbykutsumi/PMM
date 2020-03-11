from numpy import *
import myfunc.util as util
from datetime import datetime, timedelta
import random
import numpy as np
import os, sys
import socket

hostname = socket.gethostname()
if hostname == 'shui':
    tankbaseDir = '/tank'
    workbaseDir = '/work'
    figDir = '/home/utsumi/temp/stop'
elif hostname == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    workbaseDir = '/home/utsumi/mnt/lab_work'
    figDir = '/home/utsumi/temp/stop'

else:
    print 'check hostname',hostname
    sys.exit()

#lvarName = ['epc']
lvarName = ['Ku_NS_heightStormTop']


lepcid = range(0,29*29*29)

lYear = [2017]
lMon  = range(1,12+1)
lYM   = [[Year,Mon] for Year in lYear for Mon in lMon]
iMon,eMon = lMon[0],lMon[-1]


#-- Read nrec files --
drec = {}
for (Year,Mon) in lYM:
    srcPath = tankbaseDir + '/utsumi/PMM/EPCDB/list/nrec-wetcase.%04d%02d.csv'%(Year,Mon)
    f=open(srcPath,'r'); lines=f.readlines(); f.close()
    lnrecTmp = []
    for line in lines:
        line  = line.strip().split(',')
        nrec  = int(line[1])
        lnrecTmp.append(nrec)

    drec[Year,Mon] = array(lnrecTmp)

drec[0] = zeros(len(drec[lYear[0],lMon[0]]),int32)
for (Year,Mon) in lYM:
    drec[0] = drec[0] + drec[Year,Mon]


lout = []
for i in range(len(drec[0])):
    lout.append([i,drec[0][i]])

sout = util.list2csv(lout)
outPath = '/home/utsumi/temp/stop/temp.nrec-wet.csv'
f=open(outPath,'w'); f.write(sout); f.close()
print outPath


