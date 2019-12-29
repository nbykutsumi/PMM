from numpy import *
import myfunc.util as util
import os, sys
import glob
import h5py
import numpy as np
from datetime import datetime, timedelta
import socket
import EPCDB
#*******************************
iDTime = datetime(2015,5,1)
eDTime = datetime(2015,5,31)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

myhost = socket.gethostname()
if myhost =='shui':
    tankbaseDir = '/tank'
    workbaseDir = '/work'
    epcbaseDir  = '/tank/utsumi/PMM/retepc'
    dbDir   = '/tank/utsumi/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)

elif myhost == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    workbaseDir = '/home/utsumi/mnt/lab_work'
    epcbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc'
    dbDir   = '/media/disk2/share/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)

else:
    print 'check myhost'
    sys.exit()
#*******************************
thpr    = 0.1
NLEV_PREC = 50
miss_out= -9999.
db = EPCDB.EPCDB()

lvar = ['heightStormTop']
#lvar = ['tqv']
#lvar = ['vfracConv']

ddbvar = {'heightStormTop':'storm_height_ku'
        }

#------------------------------------------------

def prof250mTo500m(a3prof, miss_out):
    ny,nx,nz = a3prof.shape
    #a3out = array([ma.masked_less(a3prof[:,:,i:i+2],0).mean(2) for i in range(0,nz,2)])
    #a3out = ((ma.masked_less(a3prof[:,:,0::2],0) + ma.masked_less(a3prof[:,:,1::2],0) )*0.5).filled(miss_out)
    #a3out = ma.masked_less( concatenate([a3prof[:,:,0::2].reshape(ny,nx,-1,1), a3prof[:,:,1::2].reshape(ny,nx,-1,1)], axis=3), 0).mean(axis=3).filled(miss_out)

    a4out = np.empty([2,ny,nx,nz/2])  # 2019/11/12
    a4out[0]=a3prof[:,:,0::2]
    a4out[1]=a3prof[:,:,1::2]
    a3out = ma.masked_less(a4out ,0).mean(axis=0).filled(miss_out)

 
    return a3out

def average_2ranges_3d(a3in,miss=None,dtype=float32, fill=True):
    '''
    a3in: (ny, nx, nz) --> a2out: (ny, nx, nz/2)
    nz should be an even number
    '''
    ny,nx,nz = a3in.shape
    a4out = empty([2,ny,nx,nz/2], dtype)
    a4out[0] = a3in[:,:,::2]
    a4out[1] = a3in[:,:,1::2]
    if (miss is not None) and (fill==True):
        a3out = ma.masked_equal(a4out,miss).mean(axis=0).filled(miss).astype(dtype)
    elif (miss is not None) and (fill==False):
        a3out = ma.masked_equal(a4out,miss).mean(axis=0).astype(dtype)
    else:
        a3out = a4out.mean(axis=0)
    return a3out

def ave_9grids_3d(a3in, a1y, a1x, miss):
    '''
    returns 2-d array with the size of (nl,nz)
    a3in: (ny,nx,nz)
    nl = len(a1y)=len(a1x)
    output: (nl, nz)
    '''

    if ma.is_masked(a3in):
        a3in = a3in.filled(miss)  # 2019/12/02
    #-- Average 9 grids (over Linearlized Z)--
    nydpr,nxdpr,nzdpr= a3in.shape
    ldydx = [[dy,dx] for dy in [-1,0,1] for dx in [-1,0,1]]

    a1dprmask   = False

    a3datTmp    = empty([9,len(a1y),nzdpr], float32)

    for itmp, [dy,dx] in enumerate(ldydx):
        a1yTmp  = ma.masked_outside(a1y + dy, 0,nydpr-1)
        a1xTmp  = ma.masked_outside(a1x + dx, 0,nxdpr-1)
        a1dprmask= a1dprmask + a1yTmp.mask + a1xTmp.mask

        a2datTmp= a3in[a1yTmp.filled(0),a1xTmp.filled(0),:]

        a3datTmp[itmp,:] = a2datTmp


    a2datTmp = ma.masked_equal(a3datTmp,miss).mean(axis=0)
    a2datTmp[a1dprmask,:] = miss
    return a2datTmp


def ave_9grids_2d(a2in, a1y, a1x, miss):
    '''
    returns 1-d array with the size of (nl)
    a2in: (ny,nx)
    nl = len(a1y)=len(a1x)
    output: (nl)
    '''

    if ma.is_masked(a2in):
        a2in = a2in.filled(miss)   # 2019/12/02
    #-- Average 9 grids --
    nydpr,nxdpr = a2in.shape
    ldydx = [[dy,dx] for dy in [-1,0,1] for dx in [-1,0,1]]

    a1dprmask   = False

    a2datTmp    = empty([9,len(a1y)], float32)

    for itmp, [dy,dx] in enumerate(ldydx):
        a1yTmp  = ma.masked_outside(a1y + dy, 0,nydpr-1)
        a1xTmp  = ma.masked_outside(a1x + dx, 0,nxdpr-1)
        a1dprmask= a1dprmask + a1yTmp.mask + a1xTmp.mask

        a1datTmp= a2in[a1yTmp.filled(0),a1xTmp.filled(0)]

        a2datTmp[itmp,:] = a1datTmp

    #------------
    a1datTmp = ma.masked_equal(a2datTmp,miss).mean(axis=0)
    a1datTmp[a1dprmask] = miss


    return a1datTmp

def sum_9grids_2d(a2in, a1y, a1x, miss):
    '''
    returns 1-d array with the size of (nl)
    a2in: (ny,nx)
    nl = len(a1y)=len(a1x)
    output: (nl)
    '''

    if ma.is_masked(a2in):
        a2in = a2in.filled(miss)   # 2019/12/02
    #-- Average 9 grids (over Linearlized Z)--
    nydpr,nxdpr = a2in.shape
    ldydx = [[dy,dx] for dy in [-1,0,1] for dx in [-1,0,1]]

    a1dprmask   = False

    a2datTmp    = empty([9,len(a1y)], float32)

    for itmp, [dy,dx] in enumerate(ldydx):
        a1yTmp  = ma.masked_outside(a1y + dy, 0,nydpr-1)
        a1xTmp  = ma.masked_outside(a1x + dx, 0,nxdpr-1)
        a1dprmask= a1dprmask + a1yTmp.mask + a1xTmp.mask

        a1datTmp= a2in[a1yTmp.filled(0),a1xTmp.filled(0)]

        a2datTmp[itmp,:] = a1datTmp

    a1datTmp = ma.masked_equal(a2datTmp,miss).sum(axis=0)
    a1datTmp[a1dprmask] = miss
    return a1datTmp


#------------------------------------------------
for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]
    pmwDir = epcbaseDir + '/%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)
    ssearch  = pmwDir + '/%s.??????.y-9999--9999.nrec%05d.npy'%('prwatprofNS',DB_MAXREC)
    lpmwPath = sort(glob.glob(ssearch))
   
    for pmwPath in lpmwPath: 
        oid = int(pmwPath.split('/')[-1].split('.')[1])
        if not os.path.exists(pmwPath):
            continue
    
        ##-- Read PMW retrieval data (profile) ----------------------
        #a3profp = np.load(pmwPath)[:,83:137+1,:]  # 50 (250m) layers: 0-12.5 km.
        #a3profp = prof250mTo500m(a3profp, miss_out=-9999.)  # 25 (500m) layers

        #-- Read PMW retrieval data (surface precip) ---------------
        a2sfcprecp = np.load(pmwDir + '/nsurfMScmb.%06d.y-9999--9999.nrec%05d.npy'%(oid, DB_MAXREC))[:,83:137+1]

        #-- Read top-rank idxdb and irec
        a2topirec = np.load(pmwDir + '/top-irecNS.%06d.y-9999--9999.nrec%05d.npy'%(oid,DB_MAXREC))[:,83:137+1].astype('int32')
        a2topidxdb= np.load(pmwDir + '/top-idxdbNS.%06d.y-9999--9999.nrec%05d.npy'%(oid,DB_MAXREC))[:,83:137+1]

        #-- Reshape PMW --
        a1sfcprecp= a2sfcprecp.flatten()
        a1topirec = a2topirec.flatten()
        a1topidxdb= a2topidxdb.flatten()

        #a2profp = a3profp.reshape(-1,25)  # 25 (500m) layers
        #a2profp = a2profp[:,::-1]  # Top to bottom --> Bottom to top 
        #-- Read DPR-Combined -----------------------------------------------
        dprbaseDir = workbaseDir + '/hk01/PMM/NASA/GPM.DPRGMI/2B/V06'
        dprDir     = dprbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        ssearch = dprDir + '/2B.GPM.DPRGMI.*.%06d.V???.HDF5'%(oid)
        try:
            dprPath = glob.glob(ssearch)[0]
        except:
            print 'No DPR file for oid=',oid
            continue

        with h5py.File(dprPath, 'r') as h:
            a2sfcprecd = h['NS/surfPrecipTotRate'][:]
            a3profd = h['NS/precipTotWaterCont'][:,:,16:]     # g/m3  (0-18 g/m3), 250m layers, missing=-9999.9,  Cut-off first 16 layers (22km-18km)
    
    
        a3profd = prof250mTo500m(a3profd, miss_out=miss_out)[:,:,::-1] # convert to 500m layers up to 18km (36-layers): From bottom to top. 
   
        #-- Read DPR-Ku  ----------------------------------------
        dprbaseDir = tankbaseDir + '/utsumi/data/PMM/NASA/GPM.Ku/2A/V06'
        dprDir     = dprbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        ssearch = dprDir + '/2A.GPM.Ku.*.%06d.V06A.HDF5'%(oid)
        try:
            dprPath = glob.glob(ssearch)[0]
        except:
            print 'No Ku file for oid=',oid
            continue

        #with h5py.File(dprPath, 'r') as h:
        #    a2prtyped0= h['NS/CSF/typePrecip'][:]

        #-- Read GMI-DPR matching index file ----------------------
        xyDir = tankbaseDir + '/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/%04d/%02d/%02d'%(Year,Mon,Day)
        xPath = xyDir + '/Xpy.1.%06d.npy'%(oid)
        yPath = xyDir + '/Ypy.1.%06d.npy'%(oid)

        if not os.path.exists(xPath):
            continue    
        a1x = np.load(xPath).flatten()
        a1y = np.load(yPath).flatten()
    
        #-- Extract matching pixels from DPR array ---------------- 
        a1mask1 = ma.masked_less(a1x,0).mask
        a1mask2 = ma.masked_less(a1y,0).mask
        a1mask  = a1mask1 + a1mask2
    
        #************************
        #*** surface precip *****
        a2sfcprecd = ma.masked_less(a2sfcprecd,0)  
        a1sfcprecd = ave_9grids_2d(a2sfcprecd, a1y, a1x, miss=-9999.).filled(-9999.)


        ##*** Precipitation type *****
        #a2prtyped0 = (a2prtyped0/10000000).astype('int16')
        #DatStrat = ma.masked_not_equal(a2prtyped0,1).filled(0)
        #DatConv  = ma.masked_not_equal(a2prtyped0,2).filled(0)/2
        #DatOther = ma.masked_not_equal(a2prtyped0,3).filled(0)/3

        #DatStrat = sum_9grids_2d(DatStrat, a1y, a1x, miss=0).filled(0)
        #DatConv  = sum_9grids_2d(DatConv,  a1y, a1x, miss=0).filled(0)
        #DatOther = sum_9grids_2d(DatOther, a1y, a1x, miss=0).filled(0)
        #a1prtyped = (DatStrat + DatConv*10 + DatOther*100)

        ##*** water content profile *****
        #a2profd    = ave_9grids_3d(a3profd, a1y, a1x, miss=-9999.).filled(-9999.)  # use a1y and a1x, not a1yTmp and a1xTmp.

        a1sfcprecd[a1mask] = miss_out
        #a2profd[a1mask,:] = miss_out

        #--------------------


        #-- Screen no precipitation cases -----------------------
        a1flag1 = ma.masked_greater(a1sfcprecd, thpr).mask
        a1flag2 = ma.masked_greater(a1sfcprecp, thpr).mask
        a1flag  = a1flag1 + a1flag2 
       
        # screen a1sfcprecd==-9999. --
        a1flag3 = ma.masked_not_equal(a1sfcprecd, -9999.).mask
        a1flag4 = ma.masked_not_equal(a1sfcprecp, -9999.).mask

        a1flag  = a1flag * a1flag3 * a1flag4

 
        #a1sfcprecd = a1sfcprecd[a1flag] 
        #a1sfcprecp = a1sfcprecp[a1flag] 

        #a1prtyped  = a1prtyped[a1flag]

        #a2profd    = a2profd[a1flag]   # Bottom to top
        #a2profp    = a2profp[a1flag]   # Bottom to top

        a1topirec  = a1topirec[a1flag]
        a1topidxdb = a1topidxdb[a1flag]  
        ##*********************************
        # Read top ranked variables
        #print 'START Top'
        lidx_db = np.sort(list(set(a1topidxdb)))

        ny  = len(a1topirec)
        a1y = np.arange(ny).astype('int32')


        dout = {}
        for var in lvar:
            if var in ['heightStormTop','tqv','vfracConv']:
                dout[var] = np.ones([ny],float32)*-9999.

            else:
                print 'check var',var
                sys.exit()

    
        for idx_db in lidx_db:
            a1flagdb = ma.masked_equal(a1topidxdb, idx_db).mask
            a1irecTmp = a1topirec[a1flagdb]
            a1yTmp    = a1y[a1flagdb]

            db.set_idx_db(dbDir, idx_db)

            for var in lvar:
                if var in ddbvar.keys():
                    dbvar = ddbvar[var]
                else:
                    dbvar = var

                if len(a1yTmp)<10:
                    atmp = db.get_var(dbvar, idx=a1irecTmp)
                else:
                    atmp = db.get_var(dbvar)[a1irecTmp]

                dout[var][a1yTmp] = atmp
        #--------------------------------- 
        outDir     = tankbaseDir + '/utsumi/PMM/validprof/pair/epc.%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)
        util.mk_dir(outDir)
        #np.save(outDir + '/profpmw.%06d.npy'%(oid), a2profp.astype('float32'))
        #np.save(outDir + '/profrad.%06d.npy'%(oid), a2profd.astype('float32'))
        #np.save(outDir + '/precpmw.%06d.npy'%(oid), a1sfcprecp.astype('float32'))    
        #np.save(outDir + '/precrad.%06d.npy'%(oid), a1sfcprecd.astype('float32'))

        #np.save(outDir + '/typePreciprad.%06d.npy'%(oid), a1prtyped.astype('int16'))
        #np.save(outDir + '/top-profpmw.%06d.npy'%(oid), a2topprofp.astype('float32'))
        #np.save(outDir + '/top-typePrecippmw.%06d.npy'%(oid), a1topprtypep.astype('int16'))
        #np.save(outDir + '/top-idxdb.%06d.npy'%(oid), a1topidxdb.astype('int32'))
        #np.save(outDir + '/top-irec.%06d.npy'%(oid), a1topirec.astype('int32'))
        for var in lvar:
            outPath = outDir + '/top-%spmw.%06d.npy'%(var,oid)
            np.save(outPath, dout[var])


        print outPath, oid
