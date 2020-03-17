# %%
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import ListedColormap
%matplotlib inline

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

#calcflag = True
calcflag = False
figflag  = True
#figflag  = False

iYM = [2014,6]
eYM = [2015,5]
#nsample = 10
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
expr = 'glb.relsurf01.minrec1000.maxrec10000'

lseason= ['JJA','DJF']
lstype = ['sea','veg']
lregion= ['tro','mid']

#***************************************

    pairDir = tankbaseDir + '/utsumi/PMM/validprof/pair/%s.%s'%(rettype, expr)

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

        a1idx = ma.masked_where(a1flev-a1elev<=1000, a1idx)   # mask freezing level <= surface elevation + 1000m

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
        o2profrad.extend(a2profrad.tolist())
        o2profpmw.extend(a2profpmw.tolist())

    #*********************************
    # Make mean profiles
    #---------------------------------
    a1lat   = np.array(o1lat)
    a1lon   = np.array(o1lon)
    a1vfracconv= np.array(o1vfracconv)
    a1stype = np.array(o1stype)
    a2profrad=np.array(o2profrad)
    a2profpmw=np.array(o2profpmw)

    #-- Region ---------------
    a1flagtro  = ma.masked_inside(a1lat,-15,15).mask
    a1flagmid  = ma.masked_inside(np.abs(a1lat),35,50).mask

    d1flagreg = {}
    d1flagreg['tro'] = a1flagtro
    d1flagreg['mid'] = a1flagmid

    #-- Precip type ----------
    a1flagconv = ma.masked_greater_equal(a1vfracconv,0.5).mask
    a1flagstra = ma.masked_less(a1vfracconv,0.5).mask

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

    lkey = [(region,stype,ptype) for region in lregion
                                 for stype  in lstype
                                 for ptype  in lptype]

    for (region,stype,ptype) in lkey:
        a1flag = d1flagreg[region] * d1flagstype[stype] * d1flagptype[ptype]

        a2profradtmp = ma.masked_less(a2profrad[a1flag],0)
        a2profpmwtmp = ma.masked_less(a2profpmw[a1flag],0)

        a1averad = a2profradtmp.mean(axis=0).filled(miss)
        a1avepmw = a2profpmwtmp.mean(axis=0).filled(miss)

        a1stdrad = a2profradtmp.std(axis=0).filled(miss)
        a1stdpmw = a2profpmwtmp.std(axis=0).filled(miss)

        a1numrad = a2profradtmp.count(axis=0)
        a1numpmw = a2profpmwtmp.count(axis=0)

        a1p25rad = np.nanpercentile(a2profradtmp.filled(nan), 5, axis=0)
        a1p25pmw = np.nanpercentile(a2profpmwtmp.filled(nan), 5, axis=0)
        a1p75rad = np.nanpercentile(a2profradtmp.filled(nan), 95, axis=0)
        a1p75pmw = np.nanpercentile(a2profpmwtmp.filled(nan), 95, axis=0)

        #sys.exit()

        odir = tankbaseDir + '/utsumi/PMM/validprof/meanprof-orbit/%s.%s'%(rettype,expr)
        util.mk_dir(odir)

        stamp ='s-%s.p-%s.r-%s'%(stype,ptype,region)

        np.save(odir + '/prof.relFL.%s.ave.rad.npy'%(stamp), a1averad)
        np.save(odir + '/prof.relFL.%s.ave.pmw.npy'%(stamp), a1avepmw)
        np.save(odir + '/prof.relFL.%s.std.rad.npy'%(stamp), a1stdrad)
        np.save(odir + '/prof.relFL.%s.std.pmw.npy'%(stamp), a1stdpmw)
        np.save(odir + '/prof.relFL.%s.num.rad.npy'%(stamp), a1numrad)
        np.save(odir + '/prof.relFL.%s.num.pmw.npy'%(stamp), a1numpmw)

        np.save(odir + '/prof.relFL.%s.p25.rad.npy'%(stamp), a1p25rad)
        np.save(odir + '/prof.relFL.%s.p25.pmw.npy'%(stamp), a1p25pmw)
        np.save(odir + '/prof.relFL.%s.p75.rad.npy'%(stamp), a1p75rad)
        np.save(odir + '/prof.relFL.%s.p75.pmw.npy'%(stamp), a1p75pmw)

    nzout  = a2profrad.shape[1]
    a1relh = np.arange(-int(nz*0.5),int(nz*0.5)+1)*0.5 # km
    print a1relh
    np.save(odir + '/prof.relFL.rel-height.npy', a1relh)
    print odir



