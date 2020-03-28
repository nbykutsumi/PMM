# %%
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
iDTime = datetime(2014,11,1)
eDTime = datetime(2015,5,31)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
lskipdates = [[2014,12,9],[2014,12,10]]
useorblist = True

DB_MAXREC = 10000
DB_MINREC = 1000
#expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
expr = 'glb.relsurf01.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
#expr = 'glb.relsurf02.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

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

#lvar = [['DPRGMI','NS/Input/zeroDegAltitude']]
#lvar = [['DPRGMI','NS/Input/surfaceElevation'],['DPRGMI','NS/Input/zeroDegAltitude'],['Ku','NS/PRE/heightStormTop'],['Ku','dprx']]
#lvar = [['Ku','NS/PRE/heightStormTop']]
#lvar = [['Ku','dpry'],['Ku','dprx']]
lvar = [['DPRGMI','NS/vfracConv']] 

db = EPCDB.EPCDB()

#*******************
# oid list
#*******************
if useorblist is True:
    orblistPath= tankbaseDir + '/utsumi/PMM/retepc/list-orbit/list.GMI.well.201406-201505.1000obts.txt'
    #orblistPath= tankbaseDir + '/utsumi/PMM/retepc/list-orbit/list.GMI.well.201406-201505.2obts.txt'

    f=open(orblistPath,'r'); lines=f.readlines(); f.close()
    loid = []
    for gmiPath in lines:
        gmiPath = gmiPath.strip()
        oid = int(os.path.basename(gmiPath).split('.')[-3])
        Year,Mon,Day= map(int, os.path.dirname(gmiPath).split('/')[-3:])
        if [Year,Mon,Day] in lskipdates:
            continue
        DTimeTmp = datetime(Year,Mon,Day)
        if DTimeTmp < iDTime: continue
        if DTimeTmp > eDTime: continue
        loid.append([oid,Year,Mon,Day])

else:
    loid = []
    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]
        if [Year,Mon,Day] in lskipdates:
            continue
        pmwDir = epcbaseDir + '/%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)
        ssearch  = pmwDir + '/%s.??????.y-9999--9999.nrec%05d.npy'%('prwatprofNS',DB_MAXREC)
        lpath = sorted(glob.glob(ssearch))
        if len(lpath) == 0:
            continue
        for tmppath in lpath:
            oid = int(os.path.basename(int(tmppath)).split('.')[-4])
            loid.append([oid,Year,Mon,Day])
#------------------------------------------------


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
for oid,Year,Mon,Day in loid:
    pmwDir = epcbaseDir + '/%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)
    profPath = pmwDir + '/%s.%06d.y-9999--9999.nrec%05d.npy'%('prwatprofNS',oid,DB_MAXREC)
    if not os.path.exists(profPath):
        continue

    #-- Read PMW retrieval data (surface precip) ---------------
    a2sfcprecp = np.load(pmwDir + '/nsurfMScmb.%06d.y-9999--9999.nrec%05d.npy'%(oid, DB_MAXREC))[:,83:137+1]

    #-- Read top-rank idxdb and irec
    a2topirec = np.load(pmwDir + '/top-irecNS.%06d.y-9999--9999.nrec%05d.npy'%(oid,DB_MAXREC))[:,83:137+1].astype('int32')
    a2topidxdb= np.load(pmwDir + '/top-idxdbNS.%06d.y-9999--9999.nrec%05d.npy'%(oid,DB_MAXREC))[:,83:137+1]

    #-- Reshape PMW --
    a1sfcprecp= a2sfcprecp.flatten()
    a1topirec = a2topirec.flatten()
    a1topidxdb= a2topidxdb.flatten()

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
   
    #- Read target variables --------------------- 
    #davarorg = {}
    davarorgrad = {}
    davarorgpmw = {}
    for [prod,varName] in lvar:
        if varName in ['dpry','dprx']: continue

        #-- Read ----------------------------------------
        if prod=='DPRGMI':
            dprbaseDir = workbaseDir + '/hk01/PMM/NASA/GPM.DPRGMI/2B/V06'
            dprDir     = dprbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
            ssearch = dprDir + '/2B.GPM.DPRGMI.*.%06d.V???.HDF5'%(oid)
        #-- Read DPR-Ku  ----------------------------------------
        elif prod=='Ku':
            dprbaseDir = tankbaseDir + '/utsumi/data/PMM/NASA/GPM.Ku/2A/V06'
            dprDir     = dprbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
            ssearch = dprDir + '/2A.GPM.Ku.*.%06d.V06A.HDF5'%(oid)


        try:
            dprPath = glob.glob(ssearch)[0]
        except:
            print 'No DPR file for oid=',oid
            continue
    

        if varName=='NS/vfracConv':
            #-- DPR ---
            with h5py.File(dprPath, 'r') as h:
                davarorgrad[prod,varName] = h['/NS/Input/precipitationType'][:]
                aprecrad = h['NS/surfPrecipTotRate'][:]

            #-- PMW ---
            davarorgpmw[prod,varName] = np.load(pmwDir + '/nsurfConvNScmb.%06d.y-9999--9999.nrec%05d.npy'%(oid, DB_MAXREC))[:,83:137+1]
            aprecpmw = np.load(pmwDir + '/nsurfNScmb.%06d.y-9999--9999.nrec%05d.npy'%(oid, DB_MAXREC))[:,83:137+1]

            davarorgpmw[prod,varName] = (ma.masked_where(aprecpmw==0, davarorgpmw[prod,varName])/aprecpmw).filled(0.0)

    


        else: 
            with h5py.File(dprPath, 'r') as h:
                davarorg[prod,varName] = h[varName][:]
    
    #-- Read GMI-DPR matching index file ----------------------
    xyDir = tankbaseDir + '/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/%04d/%02d/%02d'%(Year,Mon,Day)
    xPath = xyDir + '/Xpy.1.%06d.npy'%(oid)
    yPath = xyDir + '/Ypy.1.%06d.npy'%(oid)

    if not os.path.exists(xPath):
        continue    
    a1x = np.load(xPath).flatten()
    a1y = np.load(yPath).flatten()


    for [prod,varName] in lvar:
        if varName =='dpry':
            davarorg[prod,varName] = a1y
        elif varName == 'dprx': 
            davarorg[prod,varName] = a1x

    #-- Extract matching pixels from DPR array ---------------- 
    a1mask1 = ma.masked_less(a1x,0).mask
    a1mask2 = ma.masked_less(a1y,0).mask
    a1mask  = a1mask1 + a1mask2
    
    #************************
    #*** surface precip *****
    a2sfcprecd = ma.masked_less(a2sfcprecd,0)  
    a1sfcprecd = ave_9grids_2d(a2sfcprecd, a1y, a1x, miss=-9999.).filled(-9999.)
    #-- Screen no precipitation cases -----------------------
    a1flag1 = ma.masked_greater(a1sfcprecd, thpr).mask
    a1flag2 = ma.masked_greater(a1sfcprecp, thpr).mask
    a1flag  = a1flag1 + a1flag2 
    
    # screen a1sfcprecd==-9999. --
    a1flag3 = ma.masked_not_equal(a1sfcprecd, -9999.).mask
    a1flag4 = ma.masked_not_equal(a1sfcprecp, -9999.).mask

    a1flag  = a1flag * a1flag3 * a1flag4

    davarrad = {}
    davarpmw = {}
    for [prod,varName] in lvar:
        if (prod=='DPRGMI')and(varName in ['NS/Input/zeroDegAltitude'
                                          ,'NS/Input/surfaceElevation']):
            avarorg = davarorgrad[prod,varName]
            avar = ave_9grids_2d(avarorg, a1y, a1x, miss=-9999.9).filled(-9999.)
            avar[a1mask] = -9999
            davarrad[prod,varName] = avar[a1flag]

        elif (varName =='NS/vfracConv'):
            #-- DPR ---
            avarorgrad = davarorgrad[prod,varName]
            avar0= ma.masked_less(aprecrad,0).filled(0)
            avar1= (avarorgrad/10000000).astype('int32')

            avar1= ma.masked_where(avar1 !=2, aprecrad).filled(0.0)

            avar0= sum_9grids_2d(avar0, a1y, a1x, miss=0).filled(0)
            avar1= sum_9grids_2d(avar1, a1y, a1x, miss=0).filled(0)
            avar = (ma.masked_where(avar0==0, avar1)/avar0).filled(0.0)
            avar[a1mask] = -9999.
            davarrad[prod,varName] = avar[a1flag]

            #-- PMW ---
            avar = davarorgpmw[prod,varName].flatten() # volume fraction
            avar[a1mask] = -9999.
            davarpmw[prod,varName] = avar[a1flag]  

        elif (varName in ['NS/PRE/heightStormTop']):
            avarorg = davarorgrad[prod,varName]
            avar = ma.masked_greater(avarorg,32767).filled(32767).astype('int16')  # miss is converted from -9999.9 --> -9999
            #a2stopd = ma.masked_less(a2stopd,0)
            avar = ma.masked_less_equal(avar,0) # mask zero and less
            avar = ave_9grids_2d(avar, a1y, a1x, miss=-9999).filled(-9999).astype('int16')
            avar[a1mask] = -9999
            davar[prod,varName] = avar[a1flag]

        else:
            avarorg = davarorgrad[prod,varName]
            avarorg[a1mask] = -9999
            davar[prod,varName] = avarorg[a1flag]
    #--------------------------------- 
    for [prod,varName] in lvar:
        if varName in ['NS/vfracConv']:
            outradName = varName.split('/')[-1] + 'rad'
            outpmwName = varName.split('/')[-1] + 'pmw'

        elif len(varName.split('/'))>1:
            outradName = varName.split('/')[-1] + 'rad'

        elif varName in ['dprx','dpry']:
            outradName = varName

        else:
            outradName = varName + 'rad'
        
        if varName in ['NS/Input/zeroDegAltitude'
                      ,'NS/Input/surfaceElevation'
                      ]:
            dtype='float32'

        elif varName in ['dpry','dprx']:
            dtype='int16'

        else: 
            dtype='float32'

        outDir     = tankbaseDir + '/utsumi/PMM/validprof/pair/epc.%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)
        util.mk_dir(outDir)

        if varName in ['NS/vfracConv']:
            outradPath = outDir + '/%s.%06d.npy'%(outradName,oid)
            np.save(outradPath, davarrad[prod,varName].astype(dtype)) 

            outpmwPath = outDir + '/%s.%06d.npy'%(outpmwName,oid)
            np.save(outpmwPath, davarpmw[prod,varName].astype(dtype)) 

        else:
            outradPath = outDir + '/%s.%06d.npy'%(outradName,oid)
            np.save(outradPath, davarrad[prod,varName].astype(dtype)) 


        print oid, outradPath


# %%
