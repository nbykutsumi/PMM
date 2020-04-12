# %%
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import ListedColormap
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
import scipy.stats
tankbaseDir = '/home/utsumi/mnt/lab_tank'
figDir  = '/home/utsumi/temp/ret'

calcflag = True
#calcflag = False
figflag  = True
#figflag  = False

iYM = [2014,6]
eYM = [2015,5]
#nsample = 100
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
for rettype in lrettype:
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
    o2profrad = deque([])
    o2profpmw = deque([])

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
        #-- Shift profiles relative to freezing level --
        if rel=='FL':
            dz = 500. # m
            nl = a2profrad.shape[0]
            a1idx = np.arange(nl).astype('int16')

            a2profradtmp = np.ones([nl,nz*2],'float32')*miss  # double height
            a2profpmwtmp = np.ones([nl,nz*2],'float32')*miss  # double height
            a1fzbin = (a1flev/dz).astype('int16')  # freezing level bin, bottom to top order
            setfzbin = set(list(a1fzbin))
            for fzbin in setfzbin:
                a1idxtmp = a1idx[a1fzbin==fzbin]
                a2profradtmp[a1idxtmp, nz-fzbin:2*nz-fzbin] = a2profrad[a1idxtmp,:]  # if fzbin=0, replace nz:2*nz
                a2profpmwtmp[a1idxtmp, nz-fzbin:2*nz-fzbin] = a2profpmw[a1idxtmp,:]  # if fzbin=0, replace nz:2*nz

            a2profrad = a2profradtmp[:,nz-int(nz*0.5):nz+int(nz*0.5)+1]  # nz-10:nz+11, 21 layers in total
            a2profpmw = a2profpmwtmp[:,nz-int(nz*0.5):nz+int(nz*0.5)+1]  # nz-10:nz+11, 21 layers in total

        #-------------------------------------------
        #-- Shift profiles relative to surface --
        if rel=='SF':
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

                #if sfcbin==2:
                #    idx=a1idxtmp[0]
                #    print ''
                #    print idx,a2profrad[idx]
                #    print ''
                #    print idx,a2profradtmp[idx]
                #    print ''
                #    print idx,a2profpmw[idx]
                #    print ''
                #    print idx,a2profpmwtmp[idx]
                #    sys.exit()
            a2profrad = a2profradtmp[:,nz:2*nz]  # nz:2*nz, nz layers in total
            a2profpmw = a2profpmwtmp[:,nz:2*nz]  # nz:2*nz, nz layers in total


        #-- Freezing height from surface ---
        a1fzhgt = ((ma.masked_less(a1flev,0) - a1elev)*0.001).filled(miss)

        #-- mask missing values ----
        a2mask = ma.masked_less(a2profrad,0).mask
        a2mask = a2mask + ma.masked_less(a2profpmw,0).mask
        a2profrad = ma.masked_where(a2mask, a2profrad).filled(miss)
        a2profpmw = ma.masked_where(a2mask, a2profpmw).filled(miss)

        #-- Stack ---
        o1lat.extend(list(a1lat)) 
        o1lon.extend(list(a1lon)) 
        o1vfracconv.extend(list(a1vfracconv))
        o1stype.extend(list(a1stype))
        o1fzhgt.extend(list(a1fzhgt))
        o2profrad.extend(a2profrad.tolist())
        o2profpmw.extend(a2profpmw.tolist())

    #*********************************
    # Make mean profiles
    #---------------------------------
    a1lat   = np.array(o1lat)
    a1lon   = np.array(o1lon)
    a1vfracconv= np.array(o1vfracconv)
    a1stype = np.array(o1stype)
    a1fzhgt = np.array(o1fzhgt)
    a2profrad=np.array(o2profrad)
    a2profpmw=np.array(o2profpmw)

    #-- Region ---------------
    a1flagtro  = ma.masked_inside(a1lat,-15,15).mask
    a1flagmid  = ma.masked_inside(np.abs(a1lat),35,60).mask

    d1flagreg = {}
    d1flagreg['tro'] = a1flagtro
    d1flagreg['mid'] = a1flagmid

    #-- Precip type ----------
    a1flagconv = ma.masked_greater_equal(a1vfracconv,0.6).mask
    a1flagstra = ma.masked_less(a1vfracconv,0.4).mask

    d1flagptype= {}
    d1flagptype['conv']=a1flagconv
    d1flagptype['stra']=a1flagstra

    #-- Surface type --
    a1flagsea  = ma.masked_equal(a1stype,1).mask
    a1flagveg  = ma.masked_inside(a1stype,3,7).mask
    a1flagsnow = ma.masked_inside(a1stype,8,11).mask
    a1flagcoast= ma.masked_equal(a1stype,13).mask
    a1flagland = a1flagveg + a1flagsnow

    d1flagstype= {}
    d1flagstype['sea'] = a1flagsea
    d1flagstype['veg'] = a1flagveg
    d1flagstype['snow'] = a1flagsnow
    d1flagstype['coast'] = a1flagcoast
    d1flagstype['land'] = a1flagland

    d1flagfz = {}
    for fzrange in lfzrange:
        fz0,fz1 = fzrange
        d1flagfz[fzrange] = ma.masked_inside(a1fzhgt, fz0, fz1).mask

    lkey = [(region,stype,ptype,fzrange)
            for region in lregion
            for stype  in lstype
            for ptype  in lptype
            for fzrange in lfzrange
            ]

    for (region,stype,ptype,fzrange) in lkey:
        a1flag = d1flagreg[region] * d1flagstype[stype] * d1flagptype[ptype] * d1flagfz[fzrange]

        a2profradtmp = ma.masked_less(a2profrad[a1flag],0)
        a2profpmwtmp = ma.masked_less(a2profpmw[a1flag],0)

        a1averad = a2profradtmp.mean(axis=0).filled(miss)
        a1avepmw = a2profpmwtmp.mean(axis=0).filled(miss)

        a1stdrad = a2profradtmp.std(axis=0).filled(miss)
        a1stdpmw = a2profpmwtmp.std(axis=0).filled(miss)

        a1numrad = a2profradtmp.count(axis=0)
        a1numpmw = a2profpmwtmp.count(axis=0)

        a1p05rad = np.nanpercentile(a2profradtmp.filled(nan), 5, axis=0)
        a1p05pmw = np.nanpercentile(a2profpmwtmp.filled(nan), 5, axis=0)
        a1p95rad = np.nanpercentile(a2profradtmp.filled(nan), 95, axis=0)
        a1p95pmw = np.nanpercentile(a2profpmwtmp.filled(nan), 95, axis=0)

        #sys.exit()

        odir = tankbaseDir + '/utsumi/PMM/validprof/meanprof-orbit/%s.%s'%(rettype,expr)
        util.mk_dir(odir)

        fz0, fz1 = fzrange
        stamp ='rel%s.s-%s.p-%s.r-%s.fz-%d-%d'%(rel,stype,ptype,region,fz0,fz1)

        np.save(odir + '/prof.%s.ave.rad.npy'%(stamp), a1averad)
        np.save(odir + '/prof.%s.ave.pmw.npy'%(stamp), a1avepmw)
        np.save(odir + '/prof.%s.std.rad.npy'%(stamp), a1stdrad)
        np.save(odir + '/prof.%s.std.pmw.npy'%(stamp), a1stdpmw)
        np.save(odir + '/prof.%s.num.rad.npy'%(stamp), a1numrad)
        np.save(odir + '/prof.%s.num.pmw.npy'%(stamp), a1numpmw)

        np.save(odir + '/prof.%s.p05.rad.npy'%(stamp), a1p05rad)
        np.save(odir + '/prof.%s.p05.pmw.npy'%(stamp), a1p05pmw)
        np.save(odir + '/prof.%s.p95.rad.npy'%(stamp), a1p95rad)
        np.save(odir + '/prof.%s.p95.pmw.npy'%(stamp), a1p95pmw)

    if a2profrad.shape[0]==0:
        continue

    nzout  = a2profrad.shape[1]
    if rel=='FL':
        a1hgt = np.arange(-int(nz*0.5),int(nz*0.5)+1)*0.5 # km
    elif rel=='SF':
        a1hgt = np.arange(nz)*0.5  # km

    print a1hgt
    np.save(odir + '/prof.%s.rel-height.npy'%(rel), a1hgt)
    print odir


#********************************************** 
# Figure
#********************************************** 
lkey = [(region,stype,ptype,fzrange)
        for region in lregion
        for stype  in lstype
        for ptype  in lptype
        for fzrange in lfzrange
        ]


for key in lkey:
    region,stype,ptype,fzrange = key
    fz0,fz1 = fzrange
    stamp ='rel%s.s-%s.p-%s.r-%s.fz-%d-%d'%(rel,stype,ptype,region,fz0,fz1)
    for rettype in lrettype:
        expr = dexpr[rettype]
        odir = tankbaseDir + '/utsumi/PMM/validprof/meanprof-orbit/%s.%s'%(rettype,expr)

        if rettype=='epc': 
            a1rad   = np.load(odir + '/prof.%s.ave.rad.npy'%(stamp))
            a1epc   = np.load(odir + '/prof.%s.ave.pmw.npy'%(stamp))
            a1rad05 = np.load(odir + '/prof.%s.p05.rad.npy'%(stamp))
            a1epc05 = np.load(odir + '/prof.%s.p05.pmw.npy'%(stamp))
            a1rad95 = np.load(odir + '/prof.%s.p95.rad.npy'%(stamp))
            a1epc95 = np.load(odir + '/prof.%s.p95.pmw.npy'%(stamp))

            a1numrad= np.load(odir + '/prof.%s.num.rad.npy'%(stamp))
            a1numepc= np.load(odir + '/prof.%s.num.pmw.npy'%(stamp))

            a1hgt  = np.load(odir + '/prof.%s.rel-height.npy'%(rel))

        elif rettype in ['gprof','gprof-shift']:
            a1gpr   = np.load(odir + '/prof.%s.ave.pmw.npy'%(stamp))
            a1gpr05 = np.load(odir + '/prof.%s.p05.pmw.npy'%(stamp))
            a1gpr95 = np.load(odir + '/prof.%s.p95.pmw.npy'%(stamp))
            a1numgpr= np.load(odir + '/prof.%s.num.pmw.npy'%(stamp))

    a1rad = ma.masked_less(a1rad,0)
    a1epc = ma.masked_less(a1epc,0)
    a1gpr = ma.masked_less(a1gpr,0)
    a1rad05 = ma.masked_less(a1rad05,0)
    a1epc05 = ma.masked_less(a1epc05,0)
    a1gpr05 = ma.masked_less(a1gpr05,0)
    a1rad95 = ma.masked_less(a1rad95,0)
    a1epc95 = ma.masked_less(a1epc95,0)
    a1gpr95 = ma.masked_less(a1gpr95,0)

    #-- mask with height --
    a1rad[a1hgt<=1] = np.nan
    a1epc[a1hgt<=1] = np.nan
    a1gpr[a1hgt<=1] = np.nan
    a1rad05[a1hgt<=1] = np.nan
    a1epc05[a1hgt<=1] = np.nan
    a1gpr05[a1hgt<=1] = np.nan
    a1rad95[a1hgt<=1] = np.nan
    a1epc95[a1hgt<=1] = np.nan
    a1gpr95[a1hgt<=1] = np.nan
  

    lvar = ['ave','num']
    for var in lvar:
        if var=='ave':
            a1varrad = a1rad
            a1varepc = a1epc
            a1vargpr = a1gpr
        elif var=='num':
            a1varrad = a1numrad
            a1varepc = a1numepc
            a1vargpr = a1numgpr

        if a1varrad.shape[0]==0:
            print 'No profile',key
            continue
            #sys.exit()

        fig = plt.figure(figsize=(2.5,3.2))
        #fig = plt.figure(figsize=(3,4))
        ax  = fig.add_axes([0.2,0.15,0.65,0.7])
    
        a1y = a1hgt
    
        ax.plot( a1varrad, a1y, '-',  c='k', linewidth=2, label='CMB', clip_on=False) 
        ax.plot( a1varepc, a1y, '-',  c='k', linewidth=1, label='EPC', clip_on=False)
        ax.plot( a1vargpr, a1y, '--', c='k', linewidth=1.3, label='GPROF', clip_on=False)

        #if var =='ave': 
        #    ax.plot( a1rad05, a1y, '-',  c='k', linewidth=2) 
        #    ax.plot( a1epc05, a1y, '-',  c='k', linewidth=1)
        #    ax.plot( a1gpr05, a1y, '--', c='k', linewidth=1.3)

        #    ax.plot( a1rad95, a1y, '-',  c='k', linewidth=2) 
        #    ax.plot( a1epc95, a1y, '-',  c='k', linewidth=1)
        #    ax.plot( a1gpr95, a1y, '--', c='k', linewidth=1.3)
    
        xmax = {'ave':0.65, 'num':None}[var]
        ax.set_ylim([0.9,9.9])
        ax.set_xlim([0,xmax])
        plt.xlabel('(g/m3)',fontsize=12)
        plt.ylabel('(km)',fontsize=12)
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)


        stitle = '%s %s %s %s'%(region, stype, ptype, var)
        stitle = stitle +'\n' +'FL=%d-%dkm'%(fz0,fz1)
        plt.title(stitle, fontsize=11)
        #plt.legend()   
    
        #---------------------------------
        figPath = figDir + '/prof.%s.%s.png'%(stamp,var)
        plt.savefig(figPath)
        print figPath



# %%
