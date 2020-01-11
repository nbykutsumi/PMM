import numpy as np
import pylab as pl
#import matplotlib.gridspec as gridspec
from glob import glob
import numpy.ma as ma
import sys,os, glob
from datetime import datetime, timedelta
import myfunc.util as util
from collections import deque
import socket
import h5py
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.models import model_from_json
from tensorflow.keras import layers, initializers

ldy   = [0]
#ldy   = [-1,0,1]
#ldx   = [-3,-2,-1,0,1,2,3]
ldx   = [0]
ldydx = [[dy,dx] for dy in ldy for dx in ldx]
ntc1  = 9
ntc2  = 4
ncomb = (ntc1+ntc2)* len(ldydx)
imid  = int((len(ldy)*len(ldx)-1)/2)
lsurf = ['ocean','seaice','vege','snow','swater','coast','siedge']
nsample = 10000
stopexpr = 'best01'
#stopact  = 'LTQZ'
stopact  = 'HTQZ'

#iDTime = datetime(2015,1,1)
#eDTime = datetime(2015,5,31)
iDTime = datetime(2014,9,1)
eDTime = datetime(2014,11,30)
#iDTime = datetime(2014,12,16)
#eDTime = datetime(2014,12,31)


lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
lskipdates = [[2014,10,22],[2014,10,23],[2014,10,24],[2014,11,25],[2014,12,9],[2014,12,10]]


myhost = socket.gethostname()
if myhost =='shui':
    pmmbaseDir ='/tank/utsumi/PMM'
    hdfbaseDir ='/work/hk02/PMM/NASA'
    stopbaseDir= '/tank/utsumi/PMM/stop'

    figDir = '/home.rainbow/utsumi/public_html/tempfig/stop'
elif myhost =='well':
    pmmbaseDir ='/home/utsumi/mnt/lab_tank/utsumi/PMM'
    hdfbaseDir ='/home/utsumi/mnt/lab_work/hk02/PMM/NASA'
    stopbaseDir='/home/utsumi/mnt/lab_tank/utsumi/PMM/stop'

else:
    print 'check myhost',myhost


tcmin, tcmax = 50, 350
dvarminmax = {'T':[200,350],'Q':[0,120],'Z':[-100,8800]}

#get_ipython().magic(u'matplotlib inline')
#******************************************
# Functions
#******************************************
def shift_array(ain=None, dy=None,dx=None,miss=-9999):
    ny,nx,nz = ain.shape
    aout = np.ones([ny,nx,nz]).astype(ain.dtype)*miss
    if   dy<=0: iy0=0; ey0=ny-abs(dy); iy1=abs(dy); ey1=ny
    elif dy> 0: iy0=abs(dy); ey0=ny; iy1=0; ey1=ny-abs(dy)
    if   dx<=0: ix0=0; ex0=nx-abs(dx); ix1=abs(dx); ex1=nx
    elif dx> 0: ix0=abs(dx); ex0=nx; ix1=0; ex1=nx-abs(dx)

    aout[iy0:ey0,ix0:ex0] = ain[iy1:ey1,ix1:ex1]
    return aout


def my_unit(x,Min,Max):
    return (x-Min)/(Max-Min)



def ret_lisurf(surf):
    if surf in ['ocean','seaice','vege','snow','swater','coast','siedge']:
        if   surf=='ocean':  lisurf=[1]
        elif surf=='seaice': lisurf=[2]
        elif surf=='vege':   lisurf=[3,4,5,6,7]
        elif surf=='snow':   lisurf=[8,9,10,11]
        elif surf=='swater': lisurf=[12]
        elif surf=='coast':  lisurf=[13]
        elif surf=='siedge': lisurf=[14]
    elif surf in range(1,14+1):
        lisurf = [surf]
    else:
        print 'check surf',surf
        sys.exit()

    return lisurf

#******************************************
#**************************************************************
# Read parameters for ML Stop estimation
#--------------------------------------------------------------
cpDir     = stopbaseDir + '/cp/%s-%s-ssn%s'%(stopexpr,stopact, 0)
outbaseDir= stopbaseDir + '/orbit/%s-%s-ssn%s'%(stopexpr,stopact, 0)

dmodel = {}
for surf in lsurf:
    modelPath = cpDir + '/model-s-%s.json'%(surf)
    #modelPath = cpDir + '/model.json'
    with open(modelPath,'r') as f:
        json_string = f.read()
    dmodel[surf] = model_from_json(json_string)

    cpPath = cpDir + '/cp-s-%s.ckpt'%(surf)
    dmodel[surf].load_weights(cpPath)

    print dmodel[surf] 
#**************************************************************
# Main loop start
#******************************************
for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]
    if [Year,Mon,Day] in lskipdates: continue

    #** Read HDF GMI Tb ****
    gmiDir = hdfbaseDir + '/GPM.GMI/1C/V05/%04d/%02d/%02d'%(Year,Mon,Day)

    ssearch = gmiDir + '/1C.GPM.GMI.*.??????.????.HDF5'

    lgmiPath= glob.glob(ssearch)
    lgmiPath = np.sort(lgmiPath)

    for gmiPath in lgmiPath:
        oid = int(gmiPath.split('.')[-3])

        #if oid <=4564: continue  # test

        with h5py.File(gmiPath,'r') as h:
            a2latgmi = h['/S1/Latitude'][:]
            a2longmi = h['/S1/Longitude'][:]
            atc1     = h['/S1/Tc'][:]
            if 'H' in stopact:
                atc2tmp  = h['/S2/Tc'][:]

        #** Read S2 index ******
        idxxPath = pmmbaseDir + '/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d/Xpy.1.%06d.npy'%(Year,Mon,Day,oid)
        idxyPath = pmmbaseDir + '/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d/Ypy.1.%06d.npy'%(Year,Mon,Day,oid)
    
        a2x = np.load(idxxPath)
        a2y = np.load(idxyPath)
    
        #***********************

        if 'H' in stopact:
            a2maskx = ma.masked_less(a2x,0).mask
            a2masky = ma.masked_less(a2y,0).mask
            a2maskidx  = a2maskx + a2masky
            a2x = ma.masked_where(a2maskidx, a2x).filled(0)
            a2y = ma.masked_where(a2maskidx, a2y).filled(0)
    
            atc2 = atc2tmp[a2y,a2x,:]
            for i in range(4):
                atc2[:,:,i] = ma.masked_where(a2maskidx, atc2[:,:,i]).filled(-9999.)

 
    
        #** Read HDF GPROF surface type ****
        gprofDir = hdfbaseDir + '/GPM.GMI/2A/V05/%04d/%02d/%02d'%(Year,Mon,Day)
        gprofPath= glob.glob(gprofDir + '/2A.GPM.GMI.GPROF*.%06d.????.HDF5'%(oid))[0]
    
        with h5py.File(gprofPath,'r') as h:
            a2surf = h['/S1/surfaceTypeIndex'][:]
       
        #** Read other parameters **
        dvarOrg = {}
        lvar = []
        if 'T' in stopact:
            tmpPath = pmmbaseDir + '/MATCH.GMI.V05A/S1.ABp000-220.MERRA2.t2m/%04d/%02d/%02d/t2m.%06d.npy'%(Year,Mon,Day,oid)
            dvarOrg['T'] = np.load(tmpPath).flatten()
            lvar.append('T')
    
        if 'Q' in stopact:
            tmpPath = pmmbaseDir + '/MATCH.GMI.V05A/S1.ABp000-220.MERRA2.tqv/%04d/%02d/%02d/tqv.%06d.npy'%(Year,Mon,Day,oid)
            dvarOrg['Q'] = np.load(tmpPath).flatten()
            lvar.append('Q')
    
    
        if 'Z' in stopact:
            tmpPath = pmmbaseDir + '/MATCH.GMI.V05A/S1.ABp000-220.gtopo/%04d/%02d/%02d/gtopo.%06d.npy'%(Year,Mon,Day,oid)
            dvarOrg['Z'] = np.load(tmpPath).flatten()
            lvar.append('Z')
    
    
        #** Tc vector ************

        if 'H' in stopact:
            atc  = np.concatenate([atc1,atc2],axis=2)
        elif 'L' in stopact:
            atc  = atc1
    
        aTcOrg  = np.concatenate([ shift_array(atc, dy, dx) for dy in ldy for dx in ldx], axis=2)
    
        #******************************************
        # Prep prediction
        #******************************************
        #print 'predTc.shape',aTcOrg.shape
        nz = aTcOrg.shape[2]
        aTcOrg     = aTcOrg.reshape(-1, nz)
        asurfOrg   = a2surf.flatten()
    
        #******************************************
        ny,nx = atc.shape[:2]
        predStop = np.ones(ny*nx)*(-9999.)  # 1-D
    
        #******************************************
        # Screen surf type
        #******************************************
        for surf in lsurf:
            if   surf == 'ocean':
                a1flagsurf  = ma.masked_equal(asurfOrg, 1).mask
            elif surf == 'seaice':
                a1flagsurf  = ma.masked_equal(asurfOrg, 2).mask
            elif surf == 'vege':
                a1flagsurf  = ma.masked_inside(asurfOrg, 3,7).mask
            elif surf == 'snow':
                a1flagsurf  = ma.masked_inside(asurfOrg, 8,11).mask
            elif surf == 'swater':
                a1flagsurf  = ma.masked_equal(asurfOrg, 12).mask
            elif surf == 'coast':
                a1flagsurf  = ma.masked_equal(asurfOrg, 13).mask
            elif surf == 'siedge':
                a1flagsurf  = ma.masked_equal(asurfOrg, 14).mask
            else:
                print 'check surf=',surf
                sys.exit()

            if a1flagsurf.sum() ==0:
                continue
            print 'surf',surf
            #******************************************
            # Screen invalid Tc
            #******************************************
            a1flagtc  = ma.masked_inside(aTcOrg, 50, 350).all(axis=1).mask
    
            #******************************************
            # Screen invalid data
            #******************************************
            a1flagvar = np.array([True])
            for var in lvar:
                varmin, varmax = dvarminmax[var]
                a1flagvar = a1flagvar * ma.masked_inside(dvarOrg[var],varmin, varmax).mask
                a1flagvar = a1flagvar * ~ma.masked_invalid(dvarOrg[var]).mask
    
    
            #******************************************
            # Extract
            #******************************************
            a1flag = a1flagsurf * a1flagtc * a1flagvar
            
            a1idx = np.arange(aTcOrg.shape[0]).astype('int32')
            a1idx = a1idx[a1flag]
            #print 'bef',aTcOrg.shape, a1idx.shape
             
            aTc   = aTcOrg[a1idx]
            dvar = {}
            for var in lvar:
                dvar[var] = dvarOrg[var][a1idx]
    
            #******************************************
            # Preprocess
            #******************************************
            aTc = my_unit(aTc, tcmin, tcmax)
            for var in lvar:
                #print 'dvar.shape',var,dvar[var].shape
                varmin, varmax = dvarminmax[var]
                avar= my_unit(dvar[var], varmin, varmax) 
                aTc = np.hstack([aTc, avar.reshape(-1,1)])
    
    
            #******************************************
            # Prediction
            #******************************************
            model = dmodel[surf]
            pred  = model.predict(aTc).flatten()
    
            #print ''
            #print ''
            #print 'surf=',surf
            #print 'pred',pred.shape
            #print 'predStop',predStop.shape
            #print 'a1idx',a1idx.shape
            #print 'pred.max()',pred.max()
            #print ''
            #print ''

            predStop[a1idx] = pred
             
           
        #******************************************
        # Reshape
        #******************************************
        predStop = predStop.reshape(ny,nx).astype('float32')
        print 'predStop.max=',predStop.max()

        #******************************************
        # Save
        #******************************************
        outDir = outbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        util.mk_dir(outDir)
        outPath= outDir + '/stop.%06d.npy'%(oid)
        np.save(outPath, predStop)
        print outPath
