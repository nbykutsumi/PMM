from numpy import *
from scipy.stats import gaussian_kde
import socket
import random
import os, sys, glob
import numpy as np
import epcfunc
import JPLDB
import myfunc.util as util
from datetime import datetime, timedelta
import calendar

iYM = [2017,1]
eYM = [2017,12]
lYM = util.ret_lYM(iYM, eYM)

DB_MAXREC  = 10000
DB_MINREC  = 1000

#sensor = 'GMI'
#sensor = 'AMSR2'
#sensor = 'SSMIS'
#lsensor = ['AMSR2','SSMIS','ATMS','MHS']
lsensor = ['SSMIS','ATMS','MHS']


myhost = socket.gethostname()
if myhost =="shui":
    tankbaseDir = '/tank'
    outDir  = tankbaseDir + '/utsumi/PMM/JPLDB/list'
elif myhost =="well":
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    #outDir  = '/media/disk2/share/PMM/EPCDB/list'
    outDir  = tankbaseDir + '/utsumi/PMM/JPLDB/list'


NEM_USE    = 3
NPCHIST    = 29

dlsatid = {'AMSR2':[30], 'SSMIS':[16,17,18,19], 'ATMS':[100,101], 'MHS':[201,202,318,319]}

for sensor in lsensor:
    dbDir = tankbaseDir + '/utsumi/PMM/JPLDB/EPC_DB/%s_EPC_DATABASE_TEST29'%(sensor)
    #************************
    # Initialize
    #------------------------
    lsatid = dlsatid[sensor]
    lthpr = [0,0.1,1,5,10,30]
    
    lepcid = range(NPCHIST**NEM_USE)
    danum = {}
    for satid in [999] + lsatid:
        danum[satid] = np.zeros([len(lepcid), len(lthpr)]).astype(int32)
    
    
    db  = JPLDB.JPLDB(sensor)
    #************************
    # Count
    #------------------------
    for epcid in lepcid:
        srcPath =  dbDir + '/db_%05d.bin'%(epcid)
    
        if not os.path.exists(srcPath):
            continue
    
        print sensor,epcid
        db.set_file(srcPath)
        a1prec = ma.masked_invalid(db.get_var('precip_NS_cmb'))
        a1satid= db.get_var('satid')
    
    
        #-- for satid=999 --
        satid = 999
        a1precTmp = a1prec.filled(-9999)
        a1hist =  np.histogram(a1precTmp, lthpr + [9999])[0]
        a1num  =  (a1hist[::-1].cumsum())[::-1]
        danum[satid][epcid] = a1num
     
        for satid in lsatid:
            a1precTmp = ma.masked_where(a1satid != satid, a1prec).filled(-9999)
    
            a1hist =  np.histogram(a1precTmp, lthpr + [9999])[0]
            a1num  =  (a1hist[::-1].cumsum())[::-1]
            danum[satid][epcid] = a1num
    
    ##************************
    ## Save
    ##------------------------
    for satid in [999] + lsatid:
        sout = 'epcid' +',' +  ','.join(map(str, lthpr)) + '\n'
        for epcid in lepcid:
            sout = sout + '%d'%(epcid) + ',' + ','.join(['%d'%x for x in danum[satid][epcid]]) + '\n'
        
        outPath = outDir + '/count.jpldb.%s.%d.csv'%(sensor, satid)
        util.mk_dir(outDir)
        f=open(outPath,'w'); f.write(sout); f.close()
        print outPath
