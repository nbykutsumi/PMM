# %%
import matplotlib
matplotlib.use('Agg')
#%matplotlib inline

from numpy import *
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import h5py
import numpy as np
import glob
from mpl_toolkits.basemap import Basemap
import sys, os
import myfunc.util as util

DTime = datetime(2018,1,1)
Year,Mon,Day = DTime.timetuple()[:3]

amsr2     = ["GCOMW1","AMSR2","1C","1C","V05","A"]
ssmis_f16 = ["F16","SSMIS","1C","1C","V05","B"]
ssmis_f17 = ["F17","SSMIS","1C","1C","V05","A"]
ssmis_f18 = ["F18","SSMIS","1C","1C","V05","A"]
atms_npp  = ["NPP","ATMS","1C","1C","V05","A"]
atms_noaa20= ["NOAA20","ATMS","1C","1C","V05","A"]

mhs_metopa= ["METOPA","MHS","1C","1C","V05","A"]
mhs_metopb= ["METOPB","MHS","1C","1C","V05","A"]
mhs_noaa18= ["NOAA18","MHS","1C","1C","V05","A"]
mhs_noaa19= ["NOAA19","MHS","1C","1C","V05","A"]

#lspec = [amsr2, ssmis_f16,ssmis_f17,ssmis_f18,atms_npp, atms_noaa20, mhs_metopa, mhs_metopb, mhs_noaa18, mhs_noaa19]
#lspec = [amsr2, ssmis_f16,ssmis_f17,ssmis_f18]
#lspec = [amsr2, ssmis_f16, atms_noaa20]
#lspec = [ssmis_f16]
lspec = [atms_noaa20]

dlscan = {
        'AMSR2':[1,2,3,4,5],
        'SSMIS':[1,2,3,4],
        'ATMS':[1,2,3,4],
        'MHS':[1]
        }

for spec in lspec:
    sate      = spec[0]
    sensor    = spec[1]
    prdName   = spec[2]
    prj       = spec[3]
    ver       = spec[4]
    minorver  = spec[5]

    basedir  = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/%s.%s/%s/%s'%(sate,sensor,prdName,ver)

    srcdir   = basedir + '/%04d/%02d/%02d'%(Year,Mon,Day)
    #ssearchGMI = srcDirGMI + '/1C.GPM.GMI.XCAL2016-C.20170102-S011732-E025005.016171.V05A.HDF5'
    ssearch   = srcdir + '/1C.NOAA19.MHS.XCAL2016-V.20180101-S031436-E045636.045859.V05A.HDF5'
    ssearch   = srcdir + '/*.%s.%s.*.HDF5'%(sate,sensor)
    print ssearch
    lsrcpath  = glob.glob(ssearch)
    srcpath   = lsrcpath[0]

    print srcpath

    lscan = dlscan[sensor]
    dlat = {}
    dlon = {}
    with h5py.File(srcpath,'r') as h:
        for scan in lscan:
            dlat[scan] = h['S%d/Latitude'%(scan)][:]
            dlon[scan] = h['S%d/Longitude'%(scan)][:]
            a1scori      = h['S%d/SCstatus/SCorientation'%(scan)][:]

    #--- GPROF -----
    oid = srcpath.split('.')[-3]
    gprofdir = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/%s.%s/2A-CLIM/V05/%04d/%02d/%02d'%(sate,sensor,Year,Mon,Day)
    ssearch = gprofdir + '/2A-CLIM.*.%s.????.HDF5'%(oid)
    print ssearch
    gprofpath= glob.glob(ssearch)[0]
    with h5py.File(gprofpath, 'r') as h:
        a2latgpr = h['S1/Latitude'][:]
        a2longpr = h['S1/Longitude'][:]
    #---------------

    ny,nx0 = dlat[1].shape
    #y0 = int(0.9*0.5*ny)

    nrow, ncol = 2, 2
    fig, axs = plt.subplots(nrow, ncol, figsize=(9,9))
    ipos = -1
    for latpos in ['trop','south']:
        if latpos =='south':
            y0 = 10
        else:
            y0 = int(0.8*0.25*ny)

        scori = a1scori[y0]
        for xpos in ['side','cent']:
            ipos = ipos + 1
            if xpos == 'side':
                x0 = 1
                #x0 = nx0-1
            else:
                x0 = int(0.9*nx0/2)

            #-- S1 --
            ax = axs.flatten()[ipos] 

            for scan in lscan:
                a2lat = dlat[scan]
                a2lon = dlon[scan]
                nxtmp = dlat[scan].shape[1]

                if nxtmp ==nx0:
                    ax.scatter(a2lon[y0-1:y0+2,x0-1:x0+2], a2lat[y0-1:y0+2,x0-1:x0+2], marker='$%d$'%(scan),edgecolor='k',facecolor='none',s=50)

                    ax.scatter(a2lon[y0-1,x0-1], a2lat[y0-1,x0-1], marker='o', edgecolor='k',facecolor='none',s=200)

                else:
                    y1 = y0
                    x1 = x0*2
                    ax.scatter(a2lon[y1-1:y1+2,x1-2:x1+3], a2lat[y1-1:y1+2,x1-2:x1+3], marker='$%d$'%(scan),edgecolor='r',facecolor='none',s=50)

                    ax.scatter(a2lon[y1-1,x1-2], a2lat[y1-1,x1-2], marker='o', edgecolor='r',facecolor='none',s=200)

            #- GPROF ----
            if sensor not in ['ATMS']:
                if dlat[lscan[-1]].shape[1]==nx0:
                    x1 = x0
                    ax.scatter(a2longpr[y0-1:y0+2,x1-1:x1+2], a2latgpr[y0-1:y0+2,x1-1:x1+2], marker='o', edgecolor='b',facecolor='none',s=150)
                else:
                    x1 = x0*2
                    ax.scatter(a2longpr[y0-1:y0+2,x1-2:x1+3], a2latgpr[y0-1:y0+2,x1-2:x1+3], marker='o', edgecolor='b',facecolor='none',s=150)


            #------------

            ori = {0:'ascending', 1:'descending'}[scori]
            ax.set_title('%s %s %s %s (#%s)\niscan=%d'%(sensor, latpos, xpos, ori, oid, y0-1))

            # X tick labels 
            ax.tick_params(axis='x', labelrotation=30)

    plt.tight_layout()
    figpath = '/home/utsumi/temp/multi/footprints.%s.png'%(sensor)
    plt.savefig(figpath)
    print figpath

    plt.show()

# %%
