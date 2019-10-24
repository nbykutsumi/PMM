from numpy import *
from scipy.stats import gaussian_kde
import socket
import random
import os, sys, glob
import numpy as np
import epcfunc
import EPCDB
import myfunc.util as util
from datetime import datetime, timedelta
import calendar

iYM = [2017,1]
eYM = [2017,12]
lYM = util.ret_lYM(iYM, eYM)

DB_MAXREC  = 10000
DB_MINREC  = 1000

sensor = 'GMI'
myhost = socket.gethostname()
if myhost =="shui":
    tankbaseDir = '/tank'
    outDir  = '/work/hk01/utsumi/PMM/EPCDB/list'
elif myhost =="well":
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    outDir  = '/media/disk2/share/PMM/EPCDB/list'


NEM_USE    = 3
NPCHIST    = 29

#epcbaseDir = workbaseDir + '/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp103-117.GMI.epcid-s1'
precbaseDir = tankbaseDir + '/utsumi/PMM/EPCDB/GMI.V05A.S1.ABp103-117/DPRGMI_NS_surfPrecipTotRate'
#************************
# Initialize
#------------------------
lthpr = [0,0.1,1,5,10,30]
dnum = {}
lepcid = range(NPCHIST**NEM_USE)
for i in lepcid:
    dnum[i] = []
    for thpr in lthpr:
        dnum[i].append(0)


#************************
# Count
#------------------------
for (Year,Mon) in lYM:
    lsrcPath = sort(glob.glob(precbaseDir + '/%04d%02d/DPRGMI_NS_surfPrecipTotRate.?????.npy'%(Year,Mon)))
    for srcPath in lsrcPath:
        epcid = int(srcPath.split('.')[-2])
        a1prec = ma.masked_invalid(np.load(srcPath))

        for ithpr,thpr in enumerate(lthpr): 
            num = ma.masked_less(a1prec, thpr).count()
            dnum[epcid][ithpr] = dnum[epcid][ithpr] + num
        print Year,Mon,epcid, dnum[epcid]

#************************
# Save
#------------------------
sout = 'epcid' +',' +  ','.join(map(str, lthpr)) + '\n'
for epcid in lepcid:
    sout = sout + '%d'%(epcid) + ',' + ','.join(['%d'%x for x in dnum[epcid]]) + '\n'

outPath = outDir + '/count.epc.csv'
f=open(outPath,'w'); f.write(sout); f.close()
print outPath
