import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import myfunc.util as util
#import myfunc.IO.GPM.l2_dpr as l2_dpr
import h5py
import myfunc.IO.GPM.l1_gmi as l1_gmi
import glob
from datetime import datetime, timedelta
import numpy as np
import sys
from f_match_fov import *
import os

gmi  = l1_gmi.L1_GMI()
#dpr  = l2_dpr.L2_DPR()
mwscan= 'S1'
radar = 'Ku'

#iDTime = datetime(2017,7,1)
#eDTime = datetime(2018,1,1)
#iDTime = datetime(2017,1,1)
#eDTime = datetime(2018,1,1)
iDTime = datetime(2014,9,2)
eDTime = datetime(2015,5,31)
#iDTime = datetime(2014,10,14)
#eDTime = datetime(2014,10,14)


#-- exclude missing files --
lDTimeTmp = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
lDTime = []
for DTime in lDTimeTmp:
    ##if ((datetime(2014,8,25)<=DTime) & (DTime <= datetime(2014,8,25))): continue
    #if ((datetime(2014,8,29)<=DTime) & (DTime <= datetime(2014,8,29))): continue
    #if ((datetime(2014,9,4)<=DTime) & (DTime <= datetime(2014,9,4))): continue
    #if ((datetime(2014,9,16)<=DTime) & (DTime <= datetime(2014,9,16))): continue
    #if ((datetime(2014,9,22)<=DTime) & (DTime <= datetime(2014,9,22))): continue
    #if ((datetime(2017,9,26)<=DTime) & (DTime <= datetime(2017,9,29))): continue
    #if ((datetime(2014,10,1)<=DTime) & (DTime <= datetime(2014,10,2))): continue
    #if ((datetime(2014,10,7)<=DTime) & (DTime <= datetime(2014,10,7))): continue
    #if ((datetime(2014,10,22)<=DTime) & (DTime <= datetime(2014,10,24))): continue
    #if ((datetime(2014,11,5)<=DTime) & (DTime <= datetime(2014,11,6))): continue
    #if ((datetime(2014,11,20)<=DTime) & (DTime <= datetime(2014,11,20))): continue
    #if ((datetime(2014,12,1)<=DTime) & (DTime <= datetime(2014,12,1))): continue
    #if ((datetime(2014,12,1)<=DTime) & (DTime <= datetime(2014,12,1))): continue

    lDTime.append(DTime)
#---------------------------

ix0 = 83   # in python indexing. GMI angle bins= 0, 1, 2, ..., 220 : in total=221
ex0 = 137  # in python indexing. GMI angle bins= 0, 1, 2, ..., 220 : in total=221

cx  = 110  # GMI center angle bin (py-idx)
cw  = 55    # extract this width around center
#cw  = 15    # extract this width around center
w   = int(cw/2)


verGMI = '05'
subverGMI = 'A'
fullverGMI = '%s%s'%(verGMI,subverGMI)

verDPR = '06'
subverDPR = 'A'
fullverDPR = '%s%s'%(verDPR,subverDPR)

baseDirGMI = '/work/hk02/PMM/NASA/GPM.GMI/1C/V%s'%(verGMI)
baseDirDPR = '/work/hk02/PMM/NASA/GPM.Ku/2A/V%s'%(verDPR)
idxbaseDir = '/tank/utsumi/PMM/MATCH.GMI.V%s/%s.ABp%03d-%03d.%s.V%s.IDX'%(fullverGMI, mwscan, ix0, ex0, radar, fullverDPR)

outrootDir = '/tank/utsumi/PMM/MATCH.GMI.V%s'%(fullverGMI)

#lvar = ['/NS/SLV/precipRate']
#lvar = ['/NS/CSF/typePrecip']
lvar = ['NS/PRE/heightStormTop']
#lvar = ['/NS/CSF/typePrecip','NS/PRE/heightStormTop','NS/CSF/flagAnvil','/NS/SLV/precipRate']
#lvar = ['/NS/CSF/typePrecip','NS/PRE/heightStormTop','NS/CSF/flagAnvil','/NS/SLV/precipRate']
#lvar = ['/NS/CSF/typePrecip','NS/PRE/heightStormTop','NS/CSF/flagAnvil','/NS/VER/heightZeroDeg']
#lvar = ['/NS/Latitude','/NS/Longitude']
#lvar = ['/NS/VER/heightZeroDeg']
#lvar = ['NS/SLV/precipRateESurface']

#------------------------------------------
def ave_9grids_3d(a3in, a1y, a1x, imiss, omiss, omiss_nomatch):
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

    a1matchmask   = False

    a3datTmp    = empty([9,len(a1y),nzdpr], float32)

    for itmp, [dy,dx] in enumerate(ldydx):
        a1yTmp  = ma.masked_outside(a1y + dy, 0,nydpr-1)
        a1xTmp  = ma.masked_outside(a1x + dx, 0,nxdpr-1)
        a1matchmask= a1matchmask + a1yTmp.mask + a1xTmp.mask

        a2datTmp= a3in[a1yTmp.filled(0),a1xTmp.filled(0),:]

        a3datTmp[itmp,:] = a2datTmp


    a2datTmp = ma.masked_equal(a3datTmp,imiss).mean(axis=0).filled(omiss)
    a2datTmp[a1matchmask,:] = omiss_nomatch
    return a2datTmp

#------------------------------------------
def ave_9grids_2d(a2in, a1y, a1x, imiss, omiss, omiss_nomatch):
    '''
    returns 1-d array with the size of (nl)
    a2in: (ny,nx)
    nl = len(a1y)=len(a1x)
    output: (nl)
    '''

    if ma.is_masked(a2in):
        a2in = a2in.filled(imiss)   # 2019/12/02
    #-- Average 9 grids (over Linearlized Z)--
    nydpr,nxdpr = a2in.shape
    ldydx = [[dy,dx] for dy in [-1,0,1] for dx in [-1,0,1]]

    a1matchmask   = False

    a2datTmp    = empty([9,len(a1y)], float32)


    for itmp, [dy,dx] in enumerate(ldydx):
        a1yTmp  = ma.masked_outside(a1y + dy, 0,nydpr-1)
        a1xTmp  = ma.masked_outside(a1x + dx, 0,nxdpr-1)
        a1matchmask= a1matchmask + a1yTmp.mask + a1xTmp.mask

        a1datTmp= a2in[a1yTmp.filled(0),a1xTmp.filled(0)]
        
        a2datTmp[itmp,:] = a1datTmp

    a1datTmp = (ma.masked_equal(a2datTmp,imiss).mean(axis=0)).filled(omiss)
    a1datTmp[a1matchmask] = omiss_nomatch
    return a1datTmp

#------------------------------------------
for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]

    srcDirDPR   = baseDirDPR + '/%04d/%02d/%02d'%(Year,Mon,Day)
    ssearch     = srcDirDPR  + '/2A.GPM.%s.*.V%s.HDF5'%(radar,fullverDPR)
    lsrcPathDPR = sort(glob.glob(ssearch))


    if len(lsrcPathDPR)==0:
        print 'No DPR file',Year,Mon,Day
        print ssearch
        #sys.exit()
        continue

    for srcPathDPR in lsrcPathDPR:
        oid = srcPathDPR.split('.')[-3]

        idxDir      = idxbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        idxPathX    = idxDir + '/Xpy.1.%s.npy'%(oid)
        idxPathY    = idxDir + '/Ypy.1.%s.npy'%(oid)
   
        if not os.path.exists(idxPathX):
            print 'No file'
            print idxPathX
            continue
        X    = np.load(idxPathX)
        Y    = np.load(idxPathY)
    
        a2x  = X[:,cx-w-ix0:cx+w+1-ix0]
        a2y  = Y[:,cx-w-ix0:cx+w+1-ix0]
        nygmi, nxgmi = a2x.shape
        for var in lvar:
            varName = var.split('/')[-1]
            #DatDPR = dpr.load_var_granule(srcPathDPR, var)
            hdpr   = h5py.File(srcPathDPR)
            DatDPR = hdpr[var][:]
            hdpr.close() 

            datatype   = DatDPR.dtype


            if   datatype =='float32':
                imiss = -9999.9
                omiss = -9999.
                omiss_nomatch = -9999.
            elif datatype =='int32':
                miss = -9999
                imiss = -9999
                omiss = -9999
                omiss_nomatch = -9999

            if var=='NS/PRE/heightStormTop':
                DatDPR = ma.masked_greater(DatDPR,32767).filled(32767)  # miss is converted from -9999.9 --> -9999
                DatDPR = ma.masked_less_equal(DatDPR,0) # mask zero and less
                imiss = -9999
                omiss = -9999
                datatype= 'int16'
            if   len(DatDPR.shape)==2:
                datout = ave_9grids_2d(DatDPR, a2y.flatten(), a2x.flatten(), imiss, omiss, omiss_nomatch)
                datout = datout.reshape(nygmi,nxgmi)
 
            elif len(DatDPR.shape)==3:
                datout = ave_9grids_3d(DatDPR, a2y.flatten(), a2x.flatten(), miss)
                datout = datout.reshape(nygmi,nxgmi,-1)

            if varName in ['precipRate']:
                datatype = 'int16'
                miss     = -9999
                datout   = (ma.masked_less(datout,0)*100).astype(datatype)
                datout   = ma.masked_less(datout,0).filled(miss)

            else:
                datout     = datout.astype(datatype)

            ##-- test --
            #plt.imshow(ma.masked_less_equal(datout[:150],0), origin='lower')
            #plt.colorbar()
            #plt.savefig('/home/utsumi/temp/temp.png')
            #sys.exit()
            ##-- test --

            varName    = var.split('/')[-1] 
            outbaseDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V%s/%s.ABp%03d-%03d.%s.V%s.9ave.%s'%(fullverGMI, mwscan, cx-w, cx+w, radar, fullverDPR, varName)
            outDir     = outbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
            outPath    = outDir + '/%s.%s.npy'%(varName, oid)

            util.mk_dir(outDir)
            np.save(outPath, datout)
            print outPath
