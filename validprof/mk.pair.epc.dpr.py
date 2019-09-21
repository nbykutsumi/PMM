from numpy import *
import myfunc.util as util
import os, sys
import glob
import h5py
import numpy as np
from datetime import datetime, timedelta
import socket

#*******************************
myhost = socket.gethostname()
if myhost =='shui':
    tankbaseDir = '/tank'
    workbaseDir = '/work'
    epcbaseDir  = '/tank/utsumi/PMM/retepc'
elif myhost == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    workbaseDir = '/home/utsumi/mnt/lab_work'
    epcbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc'

else:
    print 'check myhost'
    sys.exit()
#*******************************
iDTime = datetime(2014,12,1)
eDTime = datetime(2015,2,28)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
nrec   = 20000
miss_out= -9999.
#varName = 'prwatprofNS'
varName = 'top-prwatprofNS'
expr = 'glb.wprof.org'
#------------------------------------------------

def prof250mTo500m(a3prof, miss_out):
    ny,nx,nz = a3prof.shape
    #a3out = array([ma.masked_less(a3prof[:,:,i:i+2],0).mean(2) for i in range(0,nz,2)])
    a3out = ((ma.masked_less(a3prof[:,:,0::2],0) + ma.masked_less(a3prof[:,:,1::2],0) )*0.5).filled(miss_out)
    #a3out = ma.masked_less( concatenate([a3prof[:,:,0::2].reshape(ny,nx,-1,1), a3prof[:,:,1::2].reshape(ny,nx,-1,1)], axis=3), 0).mean(axis=3).filled(miss_out)
    
    return a3out
#------------------------------------------------
for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]
    pmwDir = epcbaseDir + '/%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)
    ssearch  = pmwDir + '/%s.??????.y-9999--9999.nrec%05d.npy'%(varName,nrec)
    lpmwPath = sort(glob.glob(ssearch))
   
    for pmwPath in lpmwPath: 
        oid = int(pmwPath.split('/')[-1].split('.')[1])
        if not os.path.exists(pmwPath):
            continue
    
        print 1
        #-- Read PMW retrieval data (profile) ----------------------
        a3profp = np.load(pmwPath)[:,83:137+1,:]  # 50 (250m) layers: 0-12.5 km.
        a3profp = prof250mTo500m(a3profp, miss_out=-9999.)  # 25 (500m) layers

        #-- Read PMW retrieval data (surface precip) ---------------
        a2sfcprecp = np.load(pmwDir + '/nsurfMScmb.%06d.y-9999--9999.nrec%05d.npy'%(oid, nrec))[:,83:137+1]
        #-- Reshape PMW --
        a1sfcprecp=a2sfcprecp.flatten()
        a2profp = a3profp.reshape(-1,25)  # 25 (500m) layers
        a2profp = a2profp[:,::-1]  # Top to bottom --> Bottom to top 
        #-- Read DPR -----------------------------------------------
        dprbaseDir = workbaseDir + '/hk01/PMM/NASA/GPM.DPRGMI/2B/V06'
        dprDir     = dprbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        ssearch = dprDir + '/2B.GPM.DPRGMI.*.%06d.V???.HDF5'%(oid)
        try:
            dprPath = glob.glob(ssearch)[0]
        except:
            print 'No DPR file for oid=',oid
            continue

        print 2
        with h5py.File(dprPath, 'r') as h:
            a2sfcprecd = h['NS/surfPrecipTotRate'][:]
            a3profd = h['NS/precipTotWaterCont'][:,:,16:]     # g/m3  (0-18 g/m3), 250m layers, missing=-9999.9,  Cut-off first 16 layers (22km-18km)
    
    
        a3profd = prof250mTo500m(a3profd, miss_out=miss_out)[:,:,::-1] # convert to 500m layers up to 18km (36-layers): From bottom to top. 
    
        #-- Read GMI-DPR matching index file ----------------------
        xyDir = tankbaseDir + '/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/%04d/%02d/%02d'%(Year,Mon,Day)
        xPath = xyDir + '/Xpy.1.%06d.npy'%(oid)
        yPath = xyDir + '/Ypy.1.%06d.npy'%(oid)

        print xPath
        if not os.path.exists(xPath):
            continue    
        a1x = np.load(xPath).flatten()
        a1y = np.load(yPath).flatten()
    
        print 3
        #-- Extract matching pixels from DPR array ---------------- 
        a1mask1 = ma.masked_less(a1x,0).mask
        a1mask2 = ma.masked_less(a1y,0).mask
        a1mask  = a1mask1 + a1mask2
    
        a1x = ma.masked_where(a1mask, a1x).filled(0)
        a1y = ma.masked_where(a1mask, a1y).filled(0)
    
        a1sfcprecd = a2sfcprecd[a1y,a1x]
        a2profd    = a3profd[a1y,a1x,:]
    
        a1sfcprecd[a1mask] = miss_out
        a2profd[a1mask,:] = miss_out
        #-- Screen no precipitation cases -----------------------
        a1flag1 = ma.masked_greater(a1sfcprecd, 0).mask
        a1flag2 = ma.masked_greater(a1sfcprecp, 0).mask
        a1flag  = a1flag1 + a1flag2 
       
        # screen a1sfcprecd==-9999. --
        a1flag3 = ma.masked_not_equal(a1sfcprecd, -9999.).mask
        a1flag4 = ma.masked_not_equal(a1sfcprecp, -9999.).mask
        a1flag  = a1flag * a1flag3 * a1flag4
 
        a1sfcprecd = a1sfcprecd[a1flag] 
        a1sfcprecp = a1sfcprecp[a1flag] 
        a2profd    = a2profd[a1flag]   # Bottom to top
        a2profp    = a2profp[a1flag]   # Bottom to top
   
        outDir     = tankbaseDir + '/utsumi/validprof/pair.epc/%04d/%02d/%02d'%(Year,Mon,Day)
        util.mk_dir(outDir)
        if varName =='prwatprofNS':
            np.save(outDir + '/profpmw.%06d.npy'%(oid), a2profp)
        elif varName=='top-prwatprofNS':
            np.save(outDir + '/top-profpmw.%06d.npy'%(oid), a2profp)
        np.save(outDir + '/profrad.%06d.npy'%(oid), a2profd)
        np.save(outDir + '/precpmw.%06d.npy'%(oid), a1sfcprecp)    
        np.save(outDir + '/precrad.%06d.npy'%(oid), a1sfcprecd)
        print outDir, oid
