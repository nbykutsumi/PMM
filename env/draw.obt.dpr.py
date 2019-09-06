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
else:
    print 'check myhost',myhost




#lvar = ['dtvdz_low']
#lvar = ['deptdz_low']
#lvar = ['/NS/CSF/typePrecip']
lvar = ['/NS/CSF/typePrecip','/NS/SLV/precipRateNearSurface']
oid = 19015
Year,Mon,Day,Hour = 2017,7,3,22

for var in lvar:
    varName = var.split('/')[-1]
    #---- Read DPR var and lat & lon data fron dpr --
    dprDir   = dprbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
    ssearch  = dprDir + '/2A.GPM.Ku.*.%06d.%s.HDF5'%(oid,'V06A')
    dprPath  = glob.glob(ssearch)[0]

    with h5py.File(dprPath,'r') as h5:
        a2var = h5[var][:]
        a2lat = h5['/NS/Latitude'][:]
        a2lon = h5['/NS/Longitude'][:]
    #---- draw -------
    lllat = 25
    lllon = 120
    urlat = 45
    urlon = 140
    print a2var 
    print a2var.min(), a2var.mean(),a2var.max()

    if varName =='typePrecip':
        a2var = (ma.masked_less(a2var,0)/10000000).astype(int32)
        vmin,vmax=None,None
        cbarlbl = 'mm/h'

    elif varName =='precipRateNearSurface':
        a2var = ma.masked_less_equal(a2var, 0)
        vmin,vmax=0,10
        cbarlbl = '1:Strat 2:Conv 3:Others'

    else:
        print 'check var',varName
        sys.exit()

    fig = plt.figure(figsize=(8,8))
    ax  = fig.add_axes([0.1,0.1,0.75,0.8])
    cax = fig.add_axes([0.87,0.11,0.03,0.79])
    M = Basemap(llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax, resolution='l')
    #im = M.scatter(a2lon, a2lat, c=a2var,vmin=vmin, vmax=vmax)
    im = M.pcolormesh(a2lon, a2lat, a2var,vmin=vmin, vmax=vmax)
    M.drawcoastlines()
    stitle = '%s\n%04d/%02d/%02d oid=%s'%(varName,Year,Mon,Day,oid)
    ax.set_title(stitle, fontsize=15)

    M.drawparallels(arange(-90,90,5), labels=[1,0,0,0], fontsize=12, linewidth=0.5, fmt='%d')
    M.drawmeridians(arange(-180,180,5), labels=[0,0,0,1], fontsize=12, linewidth=0.5, fmt='%d')

    #-- colorbar --
    cbar = plt.colorbar(im, cax=cax, orientation='vertical')
    cbar.set_label(cbarlbl,fontsize=13)
    cbar.ax.tick_params(labelsize=13)
    #-- Save ----
    figPath = '/home/utsumi/temp/env/orbit.dpr.%s.%s.png'%(varName,oid)
    plt.savefig(figPath)
    print figPath
