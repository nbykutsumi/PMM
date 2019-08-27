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



#lvar = ['w_low','w_mid','w_up']
#lvar = ['tcwv','mvimd']
#lvar = ['tp','cape','cin']
#lvar = ['deptdz_low']
#lvar = ['deptdz_low','deptdz_up','dtvdz_low']
lvar = ['pre.deptdz_low','pre.deptdz_up','pre.dtvdz_low']
oid = 19015
Year,Mon,Day,Hour = 2017,7,3,22
exvarDir = exvarbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
miss = -9999.

for var in lvar:
    if var in ['cape','cin','tcwv']:
        exvarPath = exvarDir +'/full.%s.%04.1fkm.%06d.npy'%(var,0,oid)

        a2var = np.load(exvarPath)
        a2var = ma.masked_less(a2var,0)
        if   var in ['cape','cin']: cbarlbl  = 'J/kg'
        else:
            print 'check var in unit',var
            cbarlbl=''
            #sys.exit()

    elif var in ['mvimd']:
        exvarPath = exvarDir +'/full.%s.%04.1fkm.%06d.npy'%(var,0,oid)

        a2var = -np.load(exvarPath) # divergence-->convergence
        #a2var = ma.masked_less(a2var,0)
        cbarlbl='kg/m2/s'

    elif var in ['tp']:
        exvarPath = exvarDir +'/full.%s.%04.1fkm.%06d.npy'%(var,0,oid)

        a2var = np.load(exvarPath)
        a2var = ma.masked_less_equal(a2var,0)*1000  # m/hour -> mm/hour
        if   var == 'tp': cbarlbl  = 'mm/h'
        else:
            print 'check var in unit',var
            sys.exit()


    elif var in ['dtvdz_low','pre.dtvdz_low']:
        dz = 4.5 - 1.5
        a2var15 = np.load(exvarDir +'/full.%s.%04.1fkm.%06d.npy'%('tv',1.5,oid))
        a2var45 = np.load(exvarDir +'/full.%s.%04.1fkm.%06d.npy'%('tv',4.5,oid))
        a2var    = -(ma.masked_less(a2var45,0) - ma.masked_less(a2var15,0)) /dz
        cbarlbl  = 'K/km'

    elif var in ['deptdz_low','pre.deptdz_low']:
        dz = 4.5 - 1.5
        a2var15 = np.load(exvarDir +'/full.%s.%04.1fkm.%06d.npy'%('ept',1.5,oid))
        a2var45 = np.load(exvarDir +'/full.%s.%04.1fkm.%06d.npy'%('ept',4.5,oid))
        a2var    = -(ma.masked_less(a2var45,0) - ma.masked_less(a2var15,0)) /dz
        cbarlbl  = 'K/km'

    elif var in ['deptdz_up','pre.deptdz_up']:
        dz = 7.5 - 4.5
        a2var45 = np.load(exvarDir +'/full.%s.%04.1fkm.%06d.npy'%('ept',4.5,oid))
        a2var75 = np.load(exvarDir +'/full.%s.%04.1fkm.%06d.npy'%('ept',7.5,oid))
        a2var    = -(ma.masked_less(a2var75,0) - ma.masked_less(a2var45,0)) /dz
        cbarlbl  = 'K/km'

    elif var in ['w_low','w_mid','w_up']:
        varTmp,lev = var.split('_')
        levkm = {'low':1.5, 'mid':4.5, 'up':7.5}[lev]
        a2var   = np.load(exvarDir +'/full.%s.%04.1fkm.%06d.npy'%(varTmp,levkm,oid))
        cbarlbl = 'm/s'

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

    #--- vmin, vmax --
    a2mask1 = ma.masked_outside(a2lon,lllon,urlon).mask
    a2mask2 = ma.masked_outside(a2lat,lllat,urlat).mask
    a2mask  = a2mask1 * a2mask2
    a2varTmp = ma.masked_where(a2mask, a2var)
    vmin    = np.percentile(a2varTmp,10)
    vmax    = np.percentile(a2varTmp,90)
    print vmin,vmax
    print a2varTmp.min(), a2varTmp.max()
    #if var in ['cape']:
    #    vmin, vmax = None,None
    #if var in ['cin']:
    #    vmin,vmax=25,45
    #elif var in ['dtvdz_low']:
    #    vmin, vmax = 3,7
    #else:
    #    vmin, vmax= None, None

    fig = plt.figure(figsize=(8,8))
    
    ax  = fig.add_axes([0.1,0.1,0.75,0.8])
    cax = fig.add_axes([0.87,0.11,0.03,0.79]) 
    M = Basemap(llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax, resolution='l')
    #im = M.scatter(a2lon, a2lat, c=a2var,vmin=vmin, vmax=vmax)
    im = M.pcolormesh(a2lon, a2lat, a2var,vmin=vmin, vmax=vmax)
    M.drawcoastlines()
    stitle = '%s\n%04d/%02d/%02d oid=%s'%(var,Year,Mon,Day,oid)
    ax.set_title(stitle, fontsize=15)

    M.drawparallels(arange(-90,90,5), labels=[1,0,0,0], fontsize=12, linewidth=0.5, fmt='%d')
    M.drawmeridians(arange(-180,180,5), labels=[0,0,0,1], fontsize=12, linewidth=0.5, fmt='%d')

    #-- colorbar --
    cbar = plt.colorbar(im, cax=cax, orientation='vertical')
    cbar.set_label(cbarlbl, fontsize=13)
    cbar.ax.tick_params(labelsize=13)


    figPath = '/home/utsumi/temp/env/orbit.era.%s.%s.png'%(var,oid)
    plt.savefig(figPath)
    print figPath

