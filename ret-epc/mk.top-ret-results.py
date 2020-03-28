# %%
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

#onlyshift = True
onlyshift = False

iDTime = datetime(2014,10,14)
eDTime = datetime(2014,10,14)
dDTime = timedelta(days=1)
lDTime = util.ret_lDTime(iDTime,eDTime,dDTime)
#lvar = ['zmNS','prwatprofNS']
#lvar = ['zmNS']
#lvar = ['prwatprofNS']
lvar = ['tbNS']
#lvar = ['heightStormTopNS']
#lvar = ['elevNS']

#** Constants ******
sensor  = 'GMI'
iscan = -9999
escan = -9999
#iscan = 917
#escan = 1117
#target_oid = None
target_oid = 3556


DB_MAXREC = 10000
DB_MINREC = 1000
NLEV_DPR = 50    # extract this number of layers
NLEV_PRECIP = 50

#expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
expr = 'glb.relsurf01.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)


myhost  = socket.gethostname()
if myhost =="shui":
    tankDir    = '/tank'
    gmibaseDir  = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05'
    matchbaseDir= '/tank/utsumi/PMM/MATCH.GMI.V05A'
    coefDir = '/tank/utsumi/PMM/EPCDB/EPC_COEF/%s'%(sensor)
    dbDir   = '/tank/utsumi/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
    outbaseDir = '/tank/utsumi/PMM/retepc/%s'%(expr)
elif myhost =="well":
    tankDir    = '/home/utsumi/mnt/lab_tank'
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
        elif varName in ['elevNS']:
            dbvarName = varName[:-2]
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
            if onlyshift is True: continue

            oid = int(os.path.basename(idxdbPath).split('.')[1])
            if (target_oid is not None) and (oid != target_oid):
                continue
            print varName, Year,Mon,Day,oid

            stamp= '%06d.y%04d-%04d.nrec%d'%(oid, iscan, escan,DB_MAXREC)
           
            idxPath  = retDir + '/top-idxdb%s.%s.npy'%(scan, stamp)
            irecPath = retDir + '/top-irec%s.%s.npy'%(scan, stamp)

            a2idxdb= np.load(idxdbPath).astype(int32)
            a2irec = np.load(irecPath ).astype(int32)

            a1idxdb = a2idxdb.flatten()
            a1irec  = a2irec.flatten()
            lidx_db = np.sort(list(set(a1idxdb)))

            #----------------

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
            if varName in ['zmNS','zmMS','prwatprofNS']:
                outvarName = varName + '-rs'
            else:
                outvarName = varName

            outDir = retDir
            #outDir  = '/home/utsumi/temp/ret'
            outPath = outDir + '/top-%s.%s.npy'%(outvarName, stamp)
            np.save(outPath, aout)
            print outPath

        #*********************************************
        # Shift top-profiles
        #*********************************************
        for idxdbPath in lidxdbPath:
            oid = int(os.path.basename(idxdbPath).split('.')[1])
            if (target_oid is not None) and (oid != target_oid):
                continue
            print varName, Year,Mon,Day,oid

            if varName in ['zmNS','zmMS','prwatprofNS']:
                scan = varName[-2:]
                stamp= '%06d.y%04d-%04d.nrec%d'%(oid, iscan, escan,DB_MAXREC)

                #-- Elevation (Obs) ---
                elevDir = tankDir + '/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.gtopo/%04d/%02d/%02d'%(Year,Mon,Day) 
                a2elevobs = np.load(elevDir + '/gtopo.%06d.npy'%(oid))

                #-- Elevation (Top) ---
                a2elevtop = np.load(retDir + '/top-elev%s.%s.npy'%(scan,stamp))

                #-- Profile data to be shifted ---
                a3in = np.load(retDir + '/top-%s-rs.%s.npy'%(varName, stamp))   # From top to bottom
                 
                #--- Shift the profile so that top-elev == obs-elev ---
                vres = 500. # m
                nl,nx,nz = a3in.shape
                tmpnz = 3*nz
                miss_int  = -9999

                a2surfbinshift = (a2elevobs/vres).astype('int16') - (a2elevtop/vres).astype('int16')

                a3tmp = np.ones([nl,nx,tmpnz], int32)* miss_int

                X,Y = np.meshgrid( np.arange(nx), np.arange(nl) )
                a1x = X.flatten().astype('int32')
                a1y = Y.flatten().astype('int32')
                a1surfbinshift = a2surfbinshift.flatten()

                for iz in range(nz):   # iz=0 (top) to iz=nz-1 (bottom)
                    a1k = nz + iz - a1surfbinshift   # if iz=0 (top) --> a1k =nz-a1surfbinshift;   if iz=nz-1 (bottom) --> a1k = 2*nz-1 - a1surfbinshift

                    a3tmp[a1y,a1x,a1k] = a3in[a1y,a1x,iz]

                a3out = a3tmp[:,:,nz:2*nz]
                #-- Save -----
                opath = retDir + '/top-%s.%s.npy'%(varName, stamp)   # From top to bottom
                np.save(opath, a3out)
                print opath

                ## %%
                #a2max = a3in.max(axis=2)
                #a2mask = ma.masked_less(a2max,0)
                #a2mask = ma.masked_where(a2surfbinshift !=2, a2mask)
                #a2mask = a2mask.mask

                #axtmp = ma.masked_where(a2mask, X).compressed()
                #aytmp = ma.masked_where(a2mask, Y).compressed()

                #i = 10
                #x,y = axtmp[i],aytmp[i]
                #print ''
                #print a3in[y,x,:]
                #print a3out[y,x,:]
                #print ''
                #print 'elev-top',a2elevtop[y,x]
                #print 'elev-obs',a2elevobs[y,x]
