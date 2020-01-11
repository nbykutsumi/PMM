#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt

import sys, os, shutil, socket
import myfunc.util as util
import glob
import subprocess
from datetime import datetime, timedelta
import numpy as np
import shutil
from numpy import *
import EPCDB

iDTime = datetime(2014,6,1)
eDTime = datetime(2014,7,31)
dDTime = timedelta(days=1)
lDTime = util.ret_lDTime(iDTime,eDTime,dDTime)
#lvar = ['zmNS','prwatprofNS']
#lvar = ['zmNS']
#lvar = ['prwatprofNS']
#lvar = ['tbNS']
lvar = ['heightStormTopNS']


#** Constants ******
sensor  = 'GMI'
iscan = -9999
escan = -9999
#iscan = 917
#escan = 1117
target_oid = None
#target_oid = 3556


DB_MAXREC = 10000
DB_MINREC = 1000
NLEV_DPR = 50    # extract this number of layers
NLEV_PRECIP = 50

expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
#expr = 'glb.v04.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)


myhost  = socket.gethostname()
if myhost =="shui":
    gmibaseDir  = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05'
    matchbaseDir= '/tank/utsumi/PMM/MATCH.GMI.V05A'
    coefDir = '/tank/utsumi/PMM/EPCDB/EPC_COEF/%s'%(sensor)
    dbDir   = '/tank/utsumi/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
    outbaseDir = '/tank/utsumi/PMM/retepc/%s'%(expr)
elif myhost =="well":
    #gmibaseDir  = '/media/disk2/share/data/PMM/NASA/GPM.GMI/1C/V05'
    gmibaseDir  = '/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.GMI/1C/V05'
    matchbaseDir= '/media/disk2/share/PMM/MATCH.GMI.V05A'
    coefDir = '/media/disk2/share/PMM/EPCDB/EPC_COEF/%s'%(sensor)
    dbDir   = '/media/disk2/share/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
    #outbaseDir = '/media/disk2/share/PMM/retepc/%s'%(expr)
    outbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/%s'%(expr)

stampshort = 'y%04d-%04d.nrec%d'%(iscan,escan,DB_MAXREC)


db = EPCDB.EPCDB()
for DTime in lDTime:
    for varName in lvar:
        scan = varName[-2:]
        #---------------------
        if varName in ['tbNS','tbMS']:
            dbvarName = 'tb'
            nz = 13
        elif varName=='zmNS':
            dbvarName = 'z_ku'
            nz = NLEV_DPR
        elif varName=='zmMS':
            dbvarName = 'z_ka'
            nz = NLEV_DPR
        elif varName in ['prwatprofNS','prwatprofMS']:
            dbvarName = 'precip_water_prof_%s'%(scan)
            nz = NLEV_PRECIP
        elif varName =='heightStormTopNS':
            dbvarName = 'storm_height_ku'
            nz = 1

        else:
            print 'check varName',varName
            sys.exit()
        #---------------------
        Year,Mon,Day = DTime.timetuple()[:3]
        retDir = outbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        ssearch = retDir + '/top-idxdb%s.*.%s.npy'%(scan,stampshort)
        lidxdbPath = np.sort(glob.glob(ssearch))
        if len(lidxdbPath)==0:
            print 'No idx_db file'
            print ssearch
            continue

        for idxdbPath in lidxdbPath:
            oid = int(os.path.basename(idxdbPath).split('.')[1])
            if (target_oid is not None) and (oid != target_oid):
                continue

            stamp= '%06d.y%04d-%04d.nrec%d'%(oid, iscan, escan,DB_MAXREC)
           
            idxPath  = retDir + '/top-idxdb%s.%s.npy'%(scan, stamp)
            irecPath = retDir + '/top-irec%s.%s.npy'%(scan, stamp)

            a2idxdb= np.load(idxdbPath).astype(int32)
            a2irec = np.load(irecPath ).astype(int32)

            a1idxdb = a2idxdb.flatten()
            a1irec  = a2irec.flatten()
            lidx_db = np.sort(list(set(a1idxdb)))

            ny,nx  = a2idxdb.shape

            X,Y = np.meshgrid(np.arange(nx), np.arange(ny))
            a1x = X.flatten()
            a1y = Y.flatten()
            aout = None

            lidx_db = lidx_db  # test
            for idx_db in lidx_db:
                if idx_db ==-9999: continue

                #print idx_db
                a1flag = ma.masked_equal(a1idxdb, idx_db).mask
                a1xTmp = a1x[a1flag]
                a1yTmp = a1y[a1flag]
                a1irecTmp  = a1irec[a1flag]

                db.set_idx_db(dbDir, idx_db)

                if varName in ['zmNS','zmMS','prwatprofNS','prwatprofMS']:
                    dat = db.get_var(dbvarName)[a1irecTmp][:,-nz:]
                else:
                    dat = db.get_var(dbvarName)[a1irecTmp]
 
                if aout is None:
                    if nz ==1:
                        aout = np.ones([ny,nx], dat.dtype)*(-9999)
                    else:
                        aout = np.ones([ny,nx,nz], dat.dtype)*(-9999)

                aout[a1yTmp, a1xTmp] = dat

            #--- Save -------------
            outDir = retDir
            #outDir  = '/home/utsumi/temp/ret'
            outPath = outDir + '/top-%s.%s.npy'%(varName, stamp)
            np.save(outPath, aout)
            print outPath
