import matplotlib
matplotlib.use('Agg')
from numpy import *
import netCDF4
import numpy as np
from mpl_toolkits.basemap import Basemap
import h5py
import glob
import matplotlib.pyplot as plt
import socket
import sys, os

myhost =socket.gethostname()

if myhost=='well':
    dprbaseDir   = '/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.Ku/2A/%s'%('V06')
    exvarbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/env/pair'
else:
    print 'check myhost',myhost




#lvar = ['dtvdz_low']
#lvar = ['deptdz_low']
lvar = ['deptdz_up']
oid = 19015
Year,Mon,Day,Hour = 2017,7,3,22
exvarDir = exvarbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)

for var in lvar:
    if var in ['cape']:
        exvarPath = exvarDir +'/full.%s.%04.1fkm.%06d.npy'%(var,0,oid)

        a2var = np.load(exvarPath)
        a2var = ma.masked_less(a2exvar,0)

    elif var in ['dtvdz_low']:
        dz = 4.5 - 1.5
        a2var15 = np.load(exvarDir +'/full.%s.%04.1fkm.%06d.npy'%('tv',1.5,oid))
        a2var45 = np.load(exvarDir +'/full.%s.%04.1fkm.%06d.npy'%('tv',4.5,oid))
        a2var    = -(ma.masked_less(a2var45,0) - ma.masked_less(a2var15,0)) /dz

    elif var in ['deptdz_low']:
        dz = 4.5 - 1.5
        a2var15 = np.load(exvarDir +'/full.%s.%04.1fkm.%06d.npy'%('ept',1.5,oid))
        a2var45 = np.load(exvarDir +'/full.%s.%04.1fkm.%06d.npy'%('ept',4.5,oid))
        a2var    = -(ma.masked_less(a2var45,0) - ma.masked_less(a2var15,0)) /dz

    elif var in ['deptdz_up']:
        dz = 7.5 - 4.5
        a2var45 = np.load(exvarDir +'/full.%s.%04.1fkm.%06d.npy'%('ept',4.5,oid))
        a2var75 = np.load(exvarDir +'/full.%s.%04.1fkm.%06d.npy'%('ept',7.5,oid))
        a2var    = -(ma.masked_less(a2var75,0) - ma.masked_less(a2var45,0)) /dz

    else:
        print 'check var',var
        sys.exit()

    #---- Read lat & lon data fron dpr --
    dprDir   = dprbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
    ssearch  = dprDir + '/2A.GPM.Ku.*.%06d.%s.HDF5'%(oid,'V06A')
    dprPath  = glob.glob(ssearch)[0]

    with h5py.File(dprPath,'r') as h5:
        a2lat = h5['/NS/Latitude'][:]
        a2lon = h5['/NS/Longitude'][:]

    #---- draw -------
    lllat = 25
    lllon = 120
    urlat = 45
    urlon = 140
    print a2var 
    print a2var.min(), a2var.mean(),a2var.max()
    if var in ['cape']:
        vmin, vmax = 200
    elif var in ['dtvdz_low']:
        vmin, vmax = 3,7
    else:
        vmin, vmax= None, None

    fig = plt.figure(figsize=(4,4))
    
    ax  = fig.add_axes([0.1,0.1,0.8,0.8])
    M = Basemap(llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax, resolution='l')
    #im = M.scatter(a2lon, a2lat, c=a2var,vmin=vmin, vmax=vmax)
    im = M.pcolormesh(a2lon, a2lat, a2var,vmin=vmin, vmax=vmax)
    plt.colorbar(im, orientation='horizontal')
    M.drawcoastlines()
    stitle = '%s\n%04d/%02d/%02d oid=%s'%(var,Year,Mon,Day,oid)
    plt.title(stitle)

    figPath = '/home/utsumi/temp/env/orbit.era.%s.%s.png'%(var,oid)
    plt.savefig(figPath)
    print figPath

