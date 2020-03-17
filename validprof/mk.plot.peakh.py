# %%
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import ListedColormap
import seaborn as sns
#%matplotlib inline

from numpy import *
import myfunc.util as util
from datetime import datetime, timedelta
import glob, os, sys
import numpy as np
import calendar
import random
import pandas as pd
from collections import deque

tankbaseDir = '/home/utsumi/mnt/lab_tank'
figDir  = '/home/utsumi/temp/ret'

calcflag = True
#calcflag = False
figflag  = True
#figflag  = False

iYM = [2014,6]
eYM = [2015,5]
#eYM = [2015,5]
#nsample = 30
nsample = 1000
lYM = util.ret_lYM(iYM,eYM)
lat0   = -60
lon0   = -180
dlatlon= 2.5
#nz     = 25  # 500m layers
nz     = 20  # 500m layers
ny,nx  = 120,360
miss   = -9999.
thpr = 0.5  # mm/h
thelev = 3000  # m
#rel  = 'FL'
rel  = 'SF'
dexpr = {'epc':'glb.relsurf01.minrec1000.maxrec10000',
         'gprof-shift':'v01'
        }

lrettype= ['epc','gprof-shift']

lfzrange = [(0,1),(1,2),(2,3),(3,4),(4,5)]
lstype = ['sea','veg']
lregion= ['tro','mid']
lptype = ['conv','stra']

#lfzrange = [(0,1)]
#lstype = ['veg']
#lregion= ['mid']
#lptype = ['conv']

#***************************************

for irettype,rettype in enumerate(lrettype):
    if calcflag is False: continue

    expr = dexpr[rettype]
    if rettype=='epc':
        pairDir = tankbaseDir + '/utsumi/PMM/validprof/pair/%s.%s'%(rettype, expr)
    elif rettype in ['gprof-shift']:
        pairDir = tankbaseDir + '/utsumi/PMM/validprof/pair/%s'%(rettype)
    else:
        print 'check rettype',rettype
    
    llatPath = []
    for Year,Mon in lYM:
        print Year,Mon
        srcDir = pairDir + '/%04d/%02d/??'%(Year,Mon)
        #srcDir = pairDir + '/%04d/%02d/01'%(Year,Mon)  # test
        ssearch= srcDir + '/Latitude.*.npy'
        llatPathTmp= glob.glob(ssearch)
        llatPath = llatPath + llatPathTmp

    random.seed(0)
    llatPath = np.sort(random.sample(llatPath, nsample))


    #-- Initialize --
    o1lat  = deque([])
    o1lon  = deque([])
    o1vfracconv= deque([])
    o1stype= deque([])
    o1fzhgt= deque([])
    o1elev = deque([])
    o1peakhrad = deque([])
    o1peakhpmw = deque([])
    #o2profrad = deque([])
    #o2profpmw = deque([])

    for ipath,latPath in enumerate(llatPath):
        print ipath
        srcDir = os.path.dirname(latPath)
        oid    = int(latPath.split('.')[-2])
        a1lat  = np.load(srcDir + '/Latitude.%06d.npy'%(oid))
        a1lon  = np.load(srcDir + '/Longitude.%06d.npy'%(oid))
        a1elev    = np.load(srcDir + '/surfaceElevationrad.%06d.npy'%(oid))
        a1flev    = np.load(srcDir + '/zeroDegAltituderad.%06d.npy'%(oid))
        a1vfracconv= np.load(srcDir + '/vfracConvrad.%06d.npy'%(oid))
        a1stype   = np.load(srcDir + '/surfaceTypeIndex.%06d.npy'%(oid))
        a1dprx    = np.load(srcDir + '/dprx.%06d.npy'%(oid))
        a1precpmw = np.load(srcDir + '/precpmw.%06d.npy'%(oid))
        a1precrad = np.load(srcDir + '/precrad.%06d.npy'%(oid))
        a2profpmw = np.load(srcDir + '/profpmw.%06d.npy'%(oid))[:,:nz]  # nz(500m) layers, bottom to top
        a2profrad = np.load(srcDir + '/profrad.%06d.npy'%(oid))[:,:nz]  # nz(500m) layers, bottom to top

        if rettype in ['gprof-shift','gprof']:
            a1qflag = np.load(srcDir + '/qualityFlag.%06d.npy'%(oid))


        #-------------------------------------------
        a1idx = range(a2profrad.shape[0])
        if len(a1idx)==0: continue

        a1idx = ma.masked_where(a1precpmw<thpr, a1idx)
        a1idx = ma.masked_where(a1precrad<thpr, a1idx)
        if rettype in ['gprof','gprof-shift']:
            a1idx = ma.masked_where(a1qflag !=0, a1idx)

        a1idx = ma.masked_where(a1elev>thelev, a1idx)   # mask elevation>thelev

        #a1idx = ma.masked_where(a1flev-a1elev<=1000, a1idx)   # mask freezing level <= surface elevation + 1000m

        a1idx = ma.masked_where(ma.masked_outside(a1dprx, 24-6,24+6).mask, a1idx)  # mask outer DPR swath (center=24 out of 0-48)
        a1idx = a1idx.compressed()

        #-- Screen ----------------
        a2profrad = ma.masked_less(a2profrad[a1idx],0)
        a2profpmw = ma.masked_less(a2profpmw[a1idx],0)
        a1lat     = a1lat[a1idx]
        a1lon     = a1lon[a1idx]
        a1elev    = a1elev[a1idx]
        a1flev    = a1flev[a1idx]
        a1vfracconv= a1vfracconv[a1idx]
        a1stype   = a1stype[a1idx]
        a1dprx    = a1dprx[a1idx]

        #-------------------------------------------
        #-------------------------------------------
        #-- Shift profiles relative to surface --
        dz = 500. # m
        nl = a2profrad.shape[0]
        a1idx = np.arange(nl).astype('int16')
        a2profradtmp = np.ones([nl,nz*2],'float32')*miss  # double height
        a2profpmwtmp = np.ones([nl,nz*2],'float32')*miss  # double height
        a1sfcbin = (a1elev/dz).astype('int16')  # freezing level bin, bottom to top order
        setsfcbin = set(list(a1sfcbin))
        for sfcbin in setsfcbin:
            a1idxtmp = a1idx[a1sfcbin==sfcbin]
            a2profradtmp[a1idxtmp, nz-sfcbin:2*nz-sfcbin] = a2profrad[a1idxtmp,:]  # if sfcbin=0, replace nz:2*nz
            a2profpmwtmp[a1idxtmp, nz-sfcbin:2*nz-sfcbin] = a2profpmw[a1idxtmp,:]  # if sfcbin=0, replace nz:2*nz

        a2profrad = a2profradtmp[:,nz:2*nz]  # nz:2*nz, nz layers in total
        a2profpmw = a2profpmwtmp[:,nz:2*nz]  # nz:2*nz, nz layers in total


        #-- Freezing height from surface ---
        a1fzhgt = ((ma.masked_less(a1flev,0) - a1elev)*0.001).filled(miss)

        #-- mask missing values ----
        a2mask = ma.masked_less(a2profrad,0).mask
        a2mask = a2mask + ma.masked_less(a2profpmw,0).mask
        a2profrad = ma.masked_where(a2mask, a2profrad).filled(miss)
        a2profpmw = ma.masked_where(a2mask, a2profpmw).filled(miss)

        #-- Peak height ----------
        a1peakhrad = a2profrad.argmax(axis=1)*0.5 + 0.25
        a1peakhpmw = a2profpmw.argmax(axis=1)*0.5 + 0.25

        #-- Stack ---
        o1lat.extend(list(a1lat)) 
        o1lon.extend(list(a1lon)) 
        o1vfracconv.extend(list(a1vfracconv))
        o1stype.extend(list(a1stype))
        o1fzhgt.extend(list(a1fzhgt))
        o1elev.extend(list(a1elev))
        o1peakhrad.extend(list(a1peakhrad))
        o1peakhpmw.extend(list(a1peakhpmw))
        #o2profrad.extend(a2profrad.tolist())
        #o2profpmw.extend(a2profpmw.tolist())


    #*********************************
    # Make mean profiles
    #---------------------------------
    a1lat   = np.array(o1lat)
    a1lon   = np.array(o1lon)
    a1vfracconv= np.array(o1vfracconv)
    a1stype = np.array(o1stype)
    a1fzhgt = np.array(o1fzhgt)
    a1elev  = np.array(o1elev)
    a1peakhrad=np.array(o1peakhrad)
    a1peakhpmw=np.array(o1peakhpmw)
    #a2profrad=np.array(o2profrad)
    #a2profpmw=np.array(o2profpmw)

    #-- Peak height ----------
    p1peakhradlabel = pd.Series(a1peakhrad.astype('int16'))

    #-- Region ---------------
    a1flagtro  = ma.masked_inside(a1lat,-15,15).mask
    a1flagmid  = ma.masked_inside(np.abs(a1lat),35,60).mask

    p1region  = pd.Series([np.nan]*a1lat.shape[0])
    p1region[a1flagtro] = 'tro'
    p1region[a1flagmid] = 'mid'

    #-- Precip type ----------
    a1flagconv = ma.masked_greater_equal(a1vfracconv,0.6).mask
    a1flagstra = ma.masked_less(a1vfracconv,0.4).mask

    p1ptype  = pd.Series([np.nan]*a1lat.shape[0])
    p1ptype[a1flagconv] = 'conv'
    p1ptype[a1flagstra] = 'stra'

    #-- Surface type --
    a1flagsea  = ma.masked_equal(a1stype,1).mask
    a1flagveg  = ma.masked_inside(a1stype,3,7).mask
    a1flagsnow = ma.masked_inside(a1stype,8,11).mask
    a1flagcoast= ma.masked_equal(a1stype,13).mask
    a1flagland = a1flagveg + a1flagsnow

    p1stype  = pd.Series([np.nan]*a1lat.shape[0])
    p1stype[a1flagsea] = 'sea'
    p1stype[a1flagveg] = 'veg'
    p1stype[a1flagsnow] = 'snow'
    p1stype[a1flagcoast] = 'coast'

    #-- Ret type ----
    p1rettype= pd.Series([rettype]*a1lat.shape[0])

    #-- DataFrame ---
    if irettype==0:
        df = pd.DataFrame({
                           'region':p1region,
                           'rettype':p1rettype,
                           'ptype':p1ptype,
                           'stype':p1stype,
                           'peakhradlabel':p1peakhradlabel,
                           'peakhrad':a1peakhrad,
                           'peakhpmw':a1peakhpmw,
        })
    else:
        df = pd.concat([
            df,
            pd.DataFrame({
                           'region':p1region,
                           'rettype':p1rettype,
                           'ptype':p1ptype,
                           'stype':p1stype,
                           'peakhradlabel':p1peakhradlabel,
                           'peakhrad':a1peakhrad,
                           'peakhpmw':a1peakhpmw,
            })
        ],axis=0)

#-- Screen -----
df = df[(df['region']=='mid')]
df = df[(df['stype']=='sea')|(df['stype']=='veg')]
df = df[(df['ptype']=='conv')|(df['ptype']=='stra')]
df = df[(df['peakhrad']>=0)&(df['peakhrad']<5)]
#df = df[(df['peakhrad']>=3)&(df['peakhrad']<4)]
#---------------

ncol = 2
nrow = 4
fig, axes = plt.subplots(nrow,ncol,figsize=(10,14))
lkey = [(rettype,ptype,stype) for rettype in lrettype for ptype in lptype for stype in lstype]
for iax,(ax, key) in enumerate(zip(axes.flat, lkey)):
    rettype,ptype,stype = key
    dftmp = df[(df['rettype']==rettype)&(df['ptype']==ptype)&(df['stype']==stype)]

    sns.boxplot(x='peakhradlabel', y='peakhpmw', ax=ax, data=dftmp)

    ax.set_xticklabels(['0-1','1-2','2-3','3-4','4-5'])
    ax.set_ylim([-2,10])
    ax.set_title('%s %s %s'%(rettype,ptype,stype))

    #-- Average peakhrad ---
    lx = [0,1,2,3,4]
    for x in lx:
        ax.plot([x-0.4,x+0.4],[x+0.5, x+0.5],'-',color='r')

fig.tight_layout()
figpath = '/home/utsumi/temp/ret/temp.peakh.violin.png'
plt.savefig(figpath)
print figpath


ncol = 2
nrow = 4
fig, axes = plt.subplots(nrow,ncol,figsize=(10,14))
lkey = [(rettype,ptype,stype) for rettype in lrettype for ptype in lptype for stype in lstype]
for iax,(ax, key) in enumerate(zip(axes.flat, lkey)):
    rettype,ptype,stype = key
    dftmp = df[(df['rettype']==rettype)&(df['ptype']==ptype)&(df['stype']==stype)]
    dcount= dftmp['peakhradlabel'].value_counts()

    sns.barplot(x=dcount.index, y=dcount.to_numpy())
    ax.set_xticklabels(['0-1','1-2','2-3','3-4','4-5'])
    ax.set_title('%s %s %s'%(rettype,ptype,stype))

    #-- Average peakhrad ---
    lx = [0,1,2,3,4]
    for x in lx:
        ax.plot([x-0.4,x+0.4],[x+0.5, x+0.5],'-',color='r')

fig.tight_layout()
figpath = '/home/utsumi/temp/ret/temp.peakh.num.png'
plt.savefig(figpath)
print figpath




# %%
