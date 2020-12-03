# %%
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import ListedColormap
#%matplotlib inline
from numpy import *
import myfunc.util as util
import os, sys
import numpy as np
import random
from collections import deque
import JPLDB, EPCDB

tankDir = '/home/utsumi/mnt/lab_tank'
figDir  = '/home/utsumi/temp/ret'

calcflag = True
#calcflag = False
figflag  = True
#figflag  = False

#lsensor = ['AMSR2','SSMIS','MHS','ATMS']
lsensor = ['GMI']
#lsensor = ['AMSR2']
lidx_db = range(29*29*29)  # 0 - 24388
#lidx_db = range(3000)
lat0   = -60
lon0   = -180
dlatlon= 2.5
#nz     = 25  # 500m layers
nz     = 20  # 500m layers
ny,nx  = 120,360
miss   = -9999.
thpr = 0.5  # mm/h
thelev = 3000  # m
expr = 'prf'


def prof250mTo500m(a2prof, miss_out):
    ''' 2D version '''
    ny,nz = a2prof.shape

    a3out = np.empty([2,ny,nz/2])  # 2019/11/12
    a3out[0]=a2prof[:,0::2]
    a3out[1]=a2prof[:,1::2]
    a2out = ma.masked_less(a3out ,0).mean(axis=0).filled(miss_out)
 
    return a2out


def gprofLayerconversion(a3prof): 
    ''' Convert GPROF 1km layers at the high altitude to 500m layers '''
    ny,nx,nz = a3prof.shape
    a3outTop = zeros([ny,nx, (nz-20)*2],float32)
    a3outTop[:,:,0::2] = a3prof[:,:,20:]
    a3outTop[:,:,1::2] = a3prof[:,:,20:]
    return concatenate([a3prof[:,:,:20], a3outTop], axis=2)

#***************************************
for sensor in lsensor:

    if calcflag is False: continue

    if sensor == 'GMI':
        db = EPCDB.EPCDB()
    else:
        db = JPLDB.JPLDB(sensor)

    #-- Initialize --
    o1lat  = deque([])
    o1lon  = deque([])
    o1vfracconv= deque([])
    o1stype= deque([])
    o1t2m  = deque([])
    o1inc  = deque([])
    o1idx_db  = deque([]) 
    o1precrad = deque([])
    o1precepc = deque([])
    o1precgpr = deque([])
    o2profrad = deque([])
    o2profepc = deque([])
    o2profgpr = deque([])

    for idx_db in lidx_db:
        tankDir = '/home/utsumi/mnt/lab_tank'
        retDir = tankDir + '/utsumi/PMM/retsynt/prf.AMSR2.smp1000/%05d'%(idx_db)

        if sensor == 'GMI':
            dbDir  = '/media/disk2/share/PMM/EPCDB/samp.10000.GMI.V05A.S1.ABp103-117.01-12'
            db.set_idx_db(dbDir, idx_db)

            retDir = tankDir + '/utsumi/PMM/retsynt/%s.%s.smp1000/%05d'%(expr,sensor,idx_db)
            if not os.path.exists(retDir):
                #print 'No file', retDir
                continue
            a1precepc = np.load(retDir + '/nsurfNScmb.est.%05d.npy'%(idx_db))
            a1precrad = np.load(retDir + '/nsurfNScmb.obs.%05d.npy'%(idx_db))
            a2profepc = np.load(retDir + '/precip_water_prof_NS_relsurf.est.%05d.npy'%(idx_db))[:,-nz*2:][:,::-1] # Top to bottom --> Bottom to top
            a1irec    = np.load(retDir + '/irec.obs.%05d.npy'%(idx_db))

            a1lat  = db.get_var('glat')[a1irec]
            a1lon  = db.get_var('glon')[a1irec]
            a1elev = db.get_var('elev')[a1irec]
            a1t2m  = db.get_var('t2m') [a1irec]

            a1inc  = np.zeros([a1lat.shape[0]],'float32')

            a1vfracconv = db.get_var('vfrac_conv_NS_cmb')[a1irec]
            a1stype     = db.get_var('sfc_class')[a1irec]
            a1precgpr   = db.get_var('precip_GPROF')[a1irec]
            a1qflag     = db.get_var('qual_flag_GPROF')[a1irec]

            a2profrad = db.get_var('precip_water_prof_NS_relsurf')[:,-nz*2:][:,::-1] # Top to bottom --> Bottom to top

            a2profgpr   = db.get_var('prof_GPROF')[a1irec,-nz:][:,::-1] # Top to bottom --> Bottom to top  

        else:
            dbDir  = '/media/disk2/share/PMM/JPLDB/EPC_DB/%s_EPC_DATABASE'%(sensor)
            dbPath = dbDir + '/db_%05d.bin'%(idx_db)
            db.set_file(dbPath)

            retDir = tankDir + '/utsumi/PMM/retsynt/%s.%s.smp1000/%05d'%(expr,sensor,idx_db)
            if not os.path.exists(retDir): continue

            a1precepc = np.load(retDir + '/nsurfNScmb.est.%05d.npy'%(idx_db))
            a1precrad = np.load(retDir + '/nsurfNScmb.obs.%05d.npy'%(idx_db))
            a2profepc = np.load(retDir + '/precip_water_prof_NS.est.%05d.npy'%(idx_db))[:,-nz*2:][:,::-1] # Top to bottom --> Bottom to top
            a2profrad = np.load(retDir + '/precip_water_prof_NS.obs.%05d.npy'%(idx_db))[:,-nz*2:][:,::-1] # Top to bottom --> Bottom to top
            a1irec    = np.load(retDir + '/irec.obs.%05d.npy'%(idx_db))

            a1lat  = db.get_var('glat')[a1irec]
            a1lon  = db.get_var('glon')[a1irec]
            a1elev = db.get_var('elev')[a1irec]
            a1t2m  = db.get_var('t2m') [a1irec]
            a1inc  = db.get_var('inc_S1') [a1irec]
            a1vfracconv = db.get_var('vfrac_conv_NS_cmb')[a1irec]
            a1stype     = db.get_var('sfc_class')[a1irec]
            a1precgpr   = db.get_var('precip_GPROF')[a1irec]
            a1qflag     = db.get_var('qual_flag_GPROF')[a1irec]
            a2profgpr   = db.get_var('prof_GPROF')[a1irec,-nz:][:,::-1] # Top to bottom --> Bottom to top  

        #-------------------------------------------
        if len(a1irec)==0: continue

        print idx_db
        a1idx = np.arange(len(a1irec))

        #-- Screen by precipitation --
        a1idx = ma.masked_where(a1precepc<thpr, a1idx)
        a1idx = ma.masked_where(a1precrad<thpr, a1idx)
        a1idx = ma.masked_where(a1precgpr<thpr, a1idx)

        #-- Quality flag -------
        a1idx = ma.masked_where(a1qflag !=0, a1idx)

        #-----------------------
        a1idx = a1idx.compressed()
        if len(a1idx)==0: continue


        #-- Screen ----------------
        a2profrad   = ma.masked_less(a2profrad[a1idx],0)
        a2profepc   = ma.masked_less(a2profepc[a1idx],0)
        a2profgpr   = ma.masked_less(a2profgpr[a1idx],0) *0.001

        a1precrad = a1precrad[a1idx]
        a1precepc = a1precepc[a1idx]
        a1precgpr = a1precgpr[a1idx]
        a1lat     = a1lat[a1idx]
        a1lon     = a1lon[a1idx]
        a1elev    = a1elev[a1idx]
        a1t2m     = a1t2m[a1idx]
        a1inc     = a1inc[a1idx]
        a1vfracconv= a1vfracconv[a1idx]
        a1stype   = a1stype[a1idx]


        #-- 250m --> 500m (EPC & CMB) ---
        a2profepc = prof250mTo500m(a2profepc, miss_out=miss)
        a2profrad = prof250mTo500m(a2profrad, miss_out=miss)

        #-- Shift profiles relative to surface --
        ''' Only for EPC and DPR '''
        dz = 500. # m
        nl = a2profrad.shape[0]
        a1idx = np.arange(nl).astype('int16')
        a2profradtmp = np.ones([nl,nz*2],'float32')*miss  # double height
        a2profepctmp = np.ones([nl,nz*2],'float32')*miss  # double height
        a1sfcbin = (a1elev/dz).astype('int16')  # freezing level bin, bottom to top order
        setsfcbin = set(list(a1sfcbin))
        for sfcbin in setsfcbin:
            a1idxtmp = a1idx[a1sfcbin==sfcbin]
            a2profradtmp[a1idxtmp, nz-sfcbin:2*nz-sfcbin] = a2profrad[a1idxtmp,:]  # if sfcbin=0, replace nz:2*nz
            a2profepctmp[a1idxtmp, nz-sfcbin:2*nz-sfcbin] = a2profepc[a1idxtmp,:]  # if sfcbin=0, replace nz:2*nz

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
        a2profepc = a2profepctmp[:,nz:2*nz]  # nz:2*nz, nz layers in total

        ##-- mask missing values ----
        #a2mask = ma.masked_less(a2profrad,0).mask
        #a2mask = a2mask + ma.masked_less(a2profepc,0).mask
        #a2mask = a2mask + ma.masked_less(a2profgpr,0).mask

        #a2profrad = ma.masked_where(a2mask, a2profrad).filled(miss)
        #a2profepc = ma.masked_where(a2mask, a2profepc).filled(miss)
        #a2profgpr = ma.masked_where(a2mask, a2profgpr).filled(miss)


        #-- Stack ---
        o1lat.extend(list(a1lat)) 
        o1lon.extend(list(a1lon)) 
        o1vfracconv.extend(list(a1vfracconv))
        o1stype.extend(list(a1stype))
        o1t2m.extend(list(a1t2m))
        o1inc.extend(list(a1inc))
        o1precrad.extend(a1precrad)
        o1precepc.extend(a1precepc)
        o1precgpr.extend(a1precgpr)
        o2profrad.extend(a2profrad.tolist())
        o2profepc.extend(a2profepc.tolist())
        o2profgpr.extend(a2profgpr.tolist())

        nrectemp = a1precrad.shape[0]
        o1idx_db.extend(np.array([idx_db]*nrectemp, 'int32'))

        print sensor,idx_db,len(o1lat)
    #*********************************
    # Make mean profiles
    #---------------------------------
    a1lat   = np.array(o1lat).astype('float32')
    a1lon   = np.array(o1lon).astype('float32')
    a1vfracconv= np.array(o1vfracconv).astype('float32')
    a1stype = np.array(o1stype).astype('int16')
    a1t2m   = np.array(o1t2m)  .astype('float32')
    a1inc   = np.array(o1inc)  .astype('float32')
    a1idx_db= np.array(o1idx_db).astype('int32')
    a1precrad=np.array(o1precrad).astype('float32')
    a1precepc=np.array(o1precepc).astype('float32')
    a1precgpr=np.array(o1precgpr).astype('float32')
    a2profrad=np.array(o2profrad).astype('float32')
    a2profepc=np.array(o2profepc).astype('float32')
    a2profgpr=np.array(o2profgpr).astype('float32')

    #-- Save data -----
    outDir = tankDir + '/utsumi/PMM/multi/pair-prof-rs/%s.%s'%(expr, sensor)
    #outDir = tankDir + '/utsumi/PMM/multi/pair-prof-rs/test-%s.%s'%(expr, sensor)
    util.mk_dir(outDir)

    np.save(outDir + '/lat.npy', a1lat)
    np.save(outDir + '/lon.npy', a1lon)
    np.save(outDir + '/vfracconv.npy', a1vfracconv)
    np.save(outDir + '/stype.npy', a1stype)
    np.save(outDir + '/t2m.npy', a1t2m)
    np.save(outDir + '/inc.npy', a1inc)
    np.save(outDir + '/idx_db.npy', a1idx_db)
    np.save(outDir + '/precrad.npy', a1precrad)
    np.save(outDir + '/precepc.npy', a1precepc)
    np.save(outDir + '/precgpr.npy', a1precgpr)
    np.save(outDir + '/profrad.npy', a2profrad)
    np.save(outDir + '/profepc.npy', a2profepc)
    np.save(outDir + '/profgpr.npy', a2profgpr)

    print outDir


# %%
