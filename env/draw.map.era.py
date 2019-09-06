import matplotlib
matplotlib.use('Agg')
from numpy import *
import netCDF4
import numpy as np
from mpl_toolkits.basemap import Basemap
import glob
import matplotlib.pyplot as plt
import socket
import sys, os
from datetime import datetime, timedelta
import myfunc.util as util

myhost =socket.gethostname()

if myhost =='shui':
    erabaseDir = '/tank/utsumi/era5'
    figDir = '/home.rainbow/utsumi/public_html/tempfig'
elif myhost == 'well':
    erabaseDir = '/media/disk2/share/data/era5'
    figDir = '/home/utsumi/mnt/lab_home_rainbow/public_html/tempfig'
else:
    print 'check myhost'
    sys.exit()


iDTime = datetime(2017,7,1,0)
eDTime = datetime(2017,7,1,0)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(hours=1))

lplev = [975,900,825,750,600,450,300,200]
lvar  = ['2t','skt']
dim   = 'single'
#dim   = 'plev'
ilev  = 0

def read_var_3d_hour(var,Year,Mon,Day,Hour):
    srcDir = erabaseDir + '/%s/%04d%02d'%(var,Year,Mon)
    srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(var,Year,Mon,Day)
    with netCDF4.Dataset(srcPath) as np:
        a3var = np.variables[var][Hour]
    return a3var


def read_var_2d_hour(var,Year,Mon,Day,Hour):
    srcDir = erabaseDir + '/%s/%04d%02d'%(var,Year,Mon)
    srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(var,Year,Mon,Day)
    ncvar = {'2t':'t2m', '2d':'d2m', 'sp':'sp','skt':'skt'}[var]

    with netCDF4.Dataset(srcPath) as np:
        a2var = np.variables[ncvar][Hour]
    return a2var


for DTime in lDTime:
    Year,Mon,Day,Hour = DTime.timetuple()[:4]
    for var in lvar:
        if dim == 'single':
            a2var = read_var_2d_hour(var, Year,Mon,Day,Hour)
        elif dim == 'plev':
            a3var = read_var_3d_hour(var, Year,Mon,Day,Hour)[::-1,:,:]
            a2var = a3var[ilev]
        else:
            print 'check dim',dim
            sys.exit()

        #*** Figure ***********
        latRA0= 90   # from North to South
        lonRA0= 0
        dlatRA = 0.25
        dlonRA = 0.25
        [[lllat,lllon],[urlat,urlon]] = [[-90,0],[90,360-0.25]]
        a1lon  = np.arange(lonRA0,360-dlonRA+0.001, dlonRA)
        a1lat  = np.arange(latRA0,-90-0.001,-dlatRA)

        a2lon,a2lat = np.meshgrid(a1lon, a1lat)


        if var in ['2t','skt']:
            vmin,vmax = 250,315
        else:
            vmin,vmax = None,None

        fig = plt.figure(figsize=(6,4))
        
        ax  = fig.add_axes([0.1,0.1,0.75,0.8])
        cax = fig.add_axes([0.87,0.11,0.03,0.79]) 
        M = Basemap(llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax, resolution='l')
        im = M.pcolormesh(a2lon, a2lat, a2var, cmap='jet',vmin=vmin, vmax=vmax)
        M.drawcoastlines()
        stitle = '%s\n%04d/%02d/%02d UTC%02d'%(var,Year,Mon,Day,Hour)
        ax.set_title(stitle, fontsize=15)

        M.drawparallels(arange(-90,90,30), labels=[1,0,0,0], fontsize=12, linewidth=0.5, fmt='%d')
        M.drawmeridians(arange(-180,180,30), labels=[0,0,0,1], fontsize=12, linewidth=0.5, fmt='%d')

        #-- colorbar --
        cbar = plt.colorbar(im, cax=cax, orientation='vertical')
        #cbar.set_label(cbarlbl, fontsize=13)
        cbar.ax.tick_params(labelsize=13)


        figPath = figDir + '/%s.%04d.%02d.%02d.%02d.%04dhPa.png'%(var,Year,Mon,Day,Hour,lplev[ilev])
        plt.savefig(figPath)
        print figPath

