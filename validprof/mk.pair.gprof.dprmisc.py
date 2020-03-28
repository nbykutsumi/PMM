# %%
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
elif myhost == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    workbaseDir = '/home/utsumi/mnt/lab_work'
else:
    print 'check myhost'
    sys.exit()
#*******************************
iDTime = datetime(2014,6,1)
eDTime = datetime(2015,5,31)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
lskipdates = [[2014,12,9],[2014,12,10]]
thpr   = 0.1
miss_out= -9999.
rettype = 'gprof-shift'
useorblist = True

#lvar = [['DPRGMI','NS/Input/zeroDegAltitude']]
#lvar = [['Ku','NS/CSF/typePrecip']]
#lvar = [['Ku','NS/CSF/typePrecip']]
#lvar = [['Ku','NS/PRE/heightStormTop'],['Ku','NS/CSF/typePrecip'],['DPRGMI','NS/Input/zeroDegAltitude'],['DPRGMI','NS/vfracConv'],['DPRGMI','NS/Input/surfaceElevation']]
#lvar = [['Ku','NS/CSF/typePrecip'],['DPRGMI','NS/Input/zeroDegAltitude'],['Ku','dprx']]
#lvar = [['DPRGMI','NS/Input/surfaceElevation']]
lvar = [['DPRGMI','NS/vfracConv']]
#lvar = [['Ku','NS/PRE/heightStormTop']]
#lvar = [['Ku','dprx']]

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

        gprofbaseDir = workbaseDir + '/hk01/PMM/NASA/GPM.GMI/2A/V05'
        gprofDir = gprofbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        ssearch  = gprofDir + '/2A.GPM.GMI.GPROF2017v1.*.??????.V05A.HDF5'
        lgprofPath = sorted(glob.glob(ssearch))
        if len(lgprofPath)==0:
            continue

        for gprofPath in lgprofPath: 
            oid = int(gprofPath.split('/')[-1].split('.')[-3])
            loid.append([oid,Year,Mon,Day])

#------------------------------------------------
def ret_aprof(a4clusterProf, a2tIndex, a3profNum, a3profScale, lspecies=[0,2,3,4]):
    nh = 28
    ny,nx = a2tIndex.shape
    a4out = empty([len(lspecies),ny, nx, nh], dtype='float32')

    for i,species in enumerate(lspecies):
        a1profScale= a3profScale[:,:,species].flatten()
        a1profNum  = a3profNum[:,:,species].flatten()
        a1tIndex   = a2tIndex.flatten()

        #-- Handle non-precipitation pixels --
        a1flag = ma.masked_equal(a1profNum, 0).mask
        a1profNum[a1flag] = 1
        a1tIndex[a1flag] = 1

        a2prof = a1profScale.reshape(-1,1) * a4clusterProf[a1profNum-1,:, a1tIndex-1, species]
        a2prof[a1flag,:] = 0.0 
        a4out[i] = a2prof.reshape(ny,nx,nh)

    return a4out


def gprofLayerconversion(a3prof): 
    ny,nx,nz = a3prof.shape
    a3outTop = zeros([ny,nx, (nz-20)*2],float32)
    a3outTop[:,:,0::2] = a3prof[:,:,20:]
    a3outTop[:,:,1::2] = a3prof[:,:,20:]
    return concatenate([a3prof[:,:,:20], a3outTop], axis=2)


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


#***************************************
for (oid,Year,Mon,Day) in loid:
    gprofbaseDir = workbaseDir + '/hk01/PMM/NASA/GPM.GMI/2A/V05'
    gprofDir  = gprofbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
    ssearch  = gprofDir + '/2A.GPM.GMI.GPROF2017v1.*.%06d.V05A.HDF5'%(oid)
    lgprofPath = sort(glob.glob(ssearch))
    if len(lgprofPath)==0:
        print 'no flies'
        print ssearch

    gprofPath = lgprofPath[0] 
    print gprofPath
 
    #-- Read GPROF ---------------------------------------------
    with h5py.File(gprofPath,'r') as h: 
        a2qFlag    = h['S1/qualityFlag'][:,83:137+1]  # (Y,X)
        a2sfcprecg = h['S1/surfacePrecipitation'][:,83:137+1]  # (Y,X)

    #-- Reshape GPROF --
    a1sfcprecg=a2sfcprecg.flatten()
    
    #-- Read DPR -----------------------------------------------
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
    davarorgrad = {}
    davarorgpmw = {}
    for [prod,varName] in lvar:

        if varName in ['dpry','dprx']: continue

        #-- Read DPR-Ku  ----------------------------------------

        elif prod=='DPRGMI':
            dprbaseDir = workbaseDir + '/hk01/PMM/NASA/GPM.DPRGMI/2B/V06'
            dprDir     = dprbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
            ssearch = dprDir + '/2B.GPM.DPRGMI.*.%06d.V???.HDF5'%(oid)
            try:
                dprPath = glob.glob(ssearch)[0]
            except:
                print 'No DPR file for oid=',oid
                continue

            if varName == 'NS/vfracConv':
                with h5py.File(dprPath, 'r') as h:
                    davarorgrad[prod,varName] = h['/NS/Input/precipitationType'][:]
                    aprecrad = h['NS/surfPrecipTotRate'][:]

                with h5py.File(gprofPath, 'r') as h:
                    aconvprectmp = h['S1/convectivePrecipitation'][:,83:137+1]
                    davarorgpmw[prod,varName] = (ma.masked_where( a2sfcprecg<=0, aconvprectmp) / a2sfcprecg).filled(0.0)


            else:
                with h5py.File(dprPath, 'r') as h:
                    davarorgrad[prod,varName] = h[varName][:]

        elif prod=='Ku':
            dprbaseDir = tankbaseDir + '/utsumi/data/PMM/NASA/GPM.Ku/2A/V06'
            dprDir     = dprbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
            ssearch    = dprDir + '/*.%06d.V06A.HDF5'%(oid)
            try:
                dprPath = glob.glob(ssearch)[0]
            except:
                print 'No DPR file for oid=',oid
                continue

            with h5py.File(dprPath, 'r') as h:
                davarorgrad[prod,varName] = h[varName][:]

    #-- Read GMI-DPR matching index file ----------------------
    xyDir = tankbaseDir + '/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/%04d/%02d/%02d'%(Year,Mon,Day)
    xPath = xyDir + '/Xpy.1.%06d.npy'%(oid)
    yPath = xyDir + '/Ypy.1.%06d.npy'%(oid)

    if not os.path.exists(xPath):
        continue    
    a1x = np.load(xPath).flatten()
    a1y = np.load(yPath).flatten()

    for [prod,varName] in lvar:
        if varName =='dprx':
            davarorgrad[prod,varName] = a1x
        elif varName =='dpry':
            davarorgrad[prod,varName] = a1y
    
    #-- Extract matching pixels from DPR array ---------------- 
    a1mask1 = ma.masked_less(a1x,0).mask
    a1mask2 = ma.masked_less(a1y,0).mask
    a1mask  = a1mask1 + a1mask2
    
    nyg,nxg = a2sfcprecg.shape

    #*** surface precip *****
    a2sfcprecd = ma.masked_less(a2sfcprecd,0)
    a1sfcprecd = ave_9grids_2d(a2sfcprecd, a1y, a1x, miss=-9999.).filled(-9999.)

    a1sfcprecd[a1mask] = miss_out
    #-- Screen no precipitation cases -----------------------
    a1flag1 = ma.masked_greater(a1sfcprecd, thpr).mask
    a1flag2 = ma.masked_greater(a1sfcprecg, thpr).mask
    a1flag  = a1flag1 + a1flag2 
    
    # screen a1sfcprecd==-9999. --
    a1flag3 = ma.masked_not_equal(a1sfcprecd, -9999.).mask
    a1flag  = a1flag * a1flag3

    davarrad = {}
    davarpmw = {}
    for [prod,varName] in lvar:
        if (prod=='DPRGMI')and(varName in ['NS/Input/zeroDegAltitude'
                                          ,'NS/Input/surfaceElevation']):
            avarorg = davarorgrad[prod,varName]
            avar = ave_9grids_2d(avarorg, a1y, a1x, miss=-9999.9).filled(-9999.)
            avar[a1mask] = -9999
            davarrad[prod,varName] = avar

        elif varName in ['NS/PRE/heightStormTop']:
            avarorg = davarorgrad[prod,varName]
            Dat0 = ma.masked_greater(avarorg,32767).filled(32767).astype('int16')  # miss is converted from -9999.9 --> -9999
            #Dat0 = ma.masked_less(Dat0,0)
            Dat0 = ma.masked_less_equal(Dat0,0) # 2019/12/02
            davarrad[prod,varName] = ave_9grids_2d(Dat0, a1y, a1x, miss=-9999).filled(-9999).astype('int16')

        elif varName in ['NS/vfracConv']:
            #-- DPR ------
            avarorg = davarorgrad[prod,varName]
            avar0= ma.masked_less(aprecrad,0).filled(0)
            avar1= (avarorg/10000000).astype('int32')

            avar1= ma.masked_where(avar1 !=2, aprecrad).filled(0.0)

            avar0= sum_9grids_2d(avar0, a1y, a1x, miss=0).filled(0)
            avar1= sum_9grids_2d(avar1, a1y, a1x, miss=0).filled(0)
            avar = (ma.masked_where(avar0==0, avar1)/avar0).filled(0.0)
            avar[a1mask] = -9999.
            davarrad[prod,varName] = avar

            #-- PMW ------
            avar = davarorgpmw[prod,varName].flatten()
            avar[a1mask] = -9999.
            davarpmw[prod,varName] = avar


        elif varName in ['NS/CSF/typePrecip']:
            Dat0 = (davarorgrad[prod,varName]/10000000).astype('int16')
            DatStrat = ma.masked_not_equal(Dat0,1).filled(0)
            DatConv  = ma.masked_not_equal(Dat0,2).filled(0)/2
            DatOther = ma.masked_not_equal(Dat0,3).filled(0)/3

            DatStrat = sum_9grids_2d(DatStrat, a1y, a1x, miss=0).filled(0)
            DatConv  = sum_9grids_2d(DatConv,  a1y, a1x, miss=0).filled(0)
            DatOther = sum_9grids_2d(DatOther, a1y, a1x, miss=0).filled(0)
            davarrad[prod,varName] = (DatStrat + DatConv*10 + DatOther*100)

        elif varName in ['dprx','dpry']:
            davarrad[prod,varName] = davarorgrad[prod,varName]

        else:
            print 'check varName',varName
            sys.exit()
    #---------------------------------
    for [prod,varName] in lvar:
        if len(varName.split('/'))>1:
            outradName = varName.split('/')[-1] + 'rad'
            outpmwName = varName.split('/')[-1] + 'pmw'

        elif varName in ['dprx','dpry']:
            outradName = varName
 
        else:
            outradName = varName + 'rad'

        if varName in ['NS/Input/zeroDegAltitude'
                      ,'NS/Input/surfaceElevation']:
            dtype='float32'

        elif varName in ['NS/PRE/heightStormTop']:
            dtype='int32'

        elif varName in ['dpry','dprx']:
            dtype='int16'

        else:
            dtype='float32'

        outDir     = tankbaseDir + '/utsumi/PMM/validprof/pair/%s/%04d/%02d/%02d'%(rettype,Year,Mon,Day)
        util.mk_dir(outDir)

        if varName in ['NS/vfracConv']:
            outradPath = outDir + '/%s.%06d.npy'%(outradName,oid)
            np.save(outradPath, davarrad[prod,varName][a1flag].astype(dtype))
            outpmwPath = outDir + '/%s.%06d.npy'%(outpmwName,oid)
            np.save(outpmwPath, davarpmw[prod,varName][a1flag].astype(dtype))

        else:
            outPath = outDir + '/%s.%06d.npy'%(outradName,oid)
            np.save(outPath, davarrad[prod,varName][a1flag].astype(dtype))

        print outradPath

# %%
