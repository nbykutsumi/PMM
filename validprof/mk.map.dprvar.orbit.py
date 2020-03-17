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

#calcflag = True
calcflag = False
figflag  = True
#figflag  = False

lseason = ['JJA','DJF']
iYM = [2014,6]
eYM = [2015,5]
nsample = 300
#nsample = 10
lYM = util.ret_lYM(iYM,eYM)
lat0   = -60
lon0   = -180
dlatlon= 2.5
#nz     = 25  # 500m layers
nz     = 20  # 500m layers
ny,nx  = 120,360
miss   = -9999.
thpr = 0.5  # mm/h
lvar = ['fzlev']
expr = 'glb.relsurf01.minrec1000.maxrec10000'
rettype='epc'
def draw_map(a2dat, a2cont=None, a1contlev=None, figPath=None, textcbartop=None, bounds=None, centers=None, extend=None, mincolor=None, maxcolor=None, ):
    res = 2.5
    fig = plt.figure(figsize=(6,3))
    #fig = plt.figure(figsize=(5,2.5))
    ax  = fig.add_axes([0.08,0.1,0.8,0.8])
    a1lat = np.arange(-60+res*0.5,60-res*0.5+0.1,res)
    a1lon = np.arange(-180+res*0.5,180-res*0.5+0.01,res)
    X,Y = np.meshgrid(a1lon,a1lat)
    M = Basemap(resolution='l', llcrnrlat=-60, llcrnrlon=-180, urcrnrlat=60, urcrnrlon=180, ax=ax)

    #-- Discrete colormap ---
    if (centers is not None)&(bounds is not None):
        print 'bounds and centers cannot be used at the same time'
        print 'bounds',bounds
        print 'centers',centers
        sys.exit()

    elif bounds is not None:
        cmap = plt.cm.get_cmap(mycm, len(bounds)+1)  # define the colormap
        cmaplist = [cmap(i) for i in range(cmap.N)]
        if mincolor is not None:
            cmaplist[0] = mincolor # grey 
        if maxcolor is not None:
            cmaplist[1] = maxcolor # grey 

        cmap = matplotlib.colors.ListedColormap(cmaplist)

        norm = matplotlib.colors.BoundaryNorm(bounds, ncolors=cmap.N, clip=False)  # normalize

    elif centers is not None:
        cmap = plt.cm.get_cmap(mycm, len(centers))  # define the colormap
        cmaplist = [cmap(i) for i in range(cmap.N)]
        if mincolor is not None:
            cmaplist[0] = mincolor # grey 
        if maxcolor is not None:
            cmaplist[1] = maxcolor # grey 

        cmap = matplotlib.colors.ListedColormap(cmaplist)
        centers= np.array(centers)
        bndmin = centers[0]-0.5*(centers[1]-centers[0])
        bndmax = centers[-1]+0.5*(centers[-1]-centers[-2])
        boundstmp = [bndmin] + list(0.5*(centers[:-1]+centers[1:])) + [bndmax]
        print centers
        print boundstmp
        norm = matplotlib.colors.BoundaryNorm(boundstmp, ncolors=cmap.N, clip=False)  # normalize


    else:
        cmap = mycm
        norm = None

    cmap.set_bad(color='silver')
    #------------------------

    im  = M.pcolormesh(X, Y, a2dat, vmin=vmin, vmax=vmax, cmap=cmap, norm=norm)

    plt.title(stitle, fontsize=15, pad=11)
    M.drawcoastlines(linewidth=1)

    M.drawmeridians(np.arange(-180,180+1,30), labels=[0,0,0,1], fontsize=10, linewidth=0.5, fmt='%d',rotation=50, yoffset=7)
    M.drawparallels(np.arange(-60,60+1,30), labels=[1,0,0,0], fontsize=10, linewidth=0.5, fmt='%d')

    cax = fig.add_axes([0.89,0.25,0.02,0.50])
    cbar= plt.colorbar(im, orientation='vertical', cax=cax)

    if bounds is not None:
        if extend=='both':
            cbar.ax.set_yticklabels([''] + list(bounds[1:-1]) + [''])
        if extend=='min':
            cbar.ax.set_yticklabels([''] + list(bounds[1:]))
        if extend=='max':
            cbar.ax.set_yticklabels(list(bounds[:-1]) + [''])
    elif centers is not None:
        cbar.set_ticks(centers)
        cbar.set_ticklabels(centers)


    cax.tick_params(labelsize=13)

    if textcbartop is not None:
        cax.text(0.1,1.05,textcbartop,fontsize=13)

    #-- Contour -----------
    if a2cont is not None:
        M.contour(X,Y, a2cont, levels=a1contlev)
    #----------------------

    plt.savefig(figPath)
    print figPath
    plt.show()
    plt.clf()  

#***************************************
lkey = [(season,var) for season in lseason
                  for var    in lvar
       ]
for key in lkey:
    season,var = key

    if calcflag is not True: continue
    #**** File list ********************
    pairDir = tankbaseDir + '/utsumi/PMM/validprof/pair/%s.%s'%(rettype, expr)

    llatPath = []
    for Year,Mon in lYM:
        if (season =='JJA'):
            if Mon not in [6,7,8]: continue
        elif season =='DJF': 
            if Mon not in [12,1,2]: continue
        else:
            print 'check season',season
            sys.exit()

        print season,Year,Mon
        srcDir = pairDir + '/%04d/%02d/??'%(Year,Mon)
        ssearch= srcDir + '/Latitude.*.npy'
        llatPathTmp= glob.glob(ssearch)
        llatPath = llatPath + llatPathTmp

    random.seed(0)
    llatPath = np.sort(random.sample(llatPath, nsample))


    #-- Initialize --
    a1lat  = deque([])
    a1lon  = deque([])
    a1var  = deque([])

    for ipath,latPath in enumerate(llatPath):
        srcDir = os.path.dirname(latPath)
        oid    = int(latPath.split('.')[-2])
        a1lattmp  = np.load(srcDir + '/Latitude.%06d.npy'%(oid))
        a1lontmp  = np.load(srcDir + '/Longitude.%06d.npy'%(oid))
        a1elev    = np.load(srcDir + '/surfaceElevationrad.%06d.npy'%(oid))

        a1prec  = np.load(srcDir + '/precrad.%06d.npy'%(oid))

        if var=='fzlev':
            a1vartmp = np.load(srcDir + '/zeroDegAltituderad.%06d.npy'%(oid))
            a1vartmp = ((ma.masked_less(a1vartmp,0) - a1elev)*0.001).filled(miss)
        #-------------------------------------------
        a1idx = range(len(a1prec))
        a1idx = ma.masked_where(a1prec<thpr, a1idx)
        a1idx = a1idx.compressed()

        print ipath,len(a1idx), a1prec.shape[0]
        print latPath
        if a1idx.sum()==0: continue

        a1lattmp = a1lattmp[a1idx]
        a1lontmp = a1lontmp[a1idx]
        a1vartmp = a1vartmp[a1idx]

        a1lat.extend(list(a1lattmp)) 
        a1lon.extend(list(a1lontmp)) 
        a1var.extend(list(a1vartmp))

    #-- Project over map ----
    lonbnd = np.arange(-180,180+0.1, dlatlon)
    latbnd = np.arange(-60,60+0.1, dlatlon)

    print len(a1lat), len(a1lon), len(a1var)

    if var == 'fzlev':
        a1var = ma.masked_equal(a1var, miss).filled(-1.0)

    a2ave = scipy.stats.binned_statistic_2d(a1lat,a1lon,a1var,statistic='mean',bins=[latbnd,lonbnd]).statistic
    a2std = scipy.stats.binned_statistic_2d(a1lat,a1lon,a1var,statistic='std',bins=[latbnd,lonbnd]).statistic
    a2num = scipy.stats.binned_statistic_2d(a1lat,a1lon,a1var,statistic='count',bins=[latbnd,lonbnd]).statistic

    odir = tankbaseDir + '/utsumi/PMM/validprof/map-orbit/%s.%s'%(rettype,expr)
    util.mk_dir(odir)
    avepth = odir + '/ave.%s.rad.%s.npy'%(var,season)
    stdpth = odir + '/std.%s.rad.%s.npy'%(var,season)
    numpth = odir + '/num.%s.rad.%s.npy'%(var,season)
    np.save(avepth, a2ave.astype('float32'))
    np.save(stdpth, a2std.astype('float32'))
    np.save(numpth, a2num.astype('int32'))
    print rettype,a2ave.shape
    print avepth

    #plt.imshow(ma.masked_invalid(a2ave), origin='lower'); plt.colorbar()
    #plt.show() 

#**************************************
# Figure
#--------------------------------------
for season in lseason:
    for var in lvar:
        if figflag is not True: continue

        odir = tankbaseDir + '/utsumi/PMM/validprof/map-orbit/%s.%s'%(rettype,expr)

        dbndave = {
                    #'fzlev':np.arange(-1,6+0.1,1),
                    'fzlev':np.arange(-0.5,6.5+0.1,1),
                   }
        dbndstd = {'fzlev':np.arange(0,2+0.1,0.5),
                   }
        dcntave = {'dpeakh':np.arange(-1.5,1.5+0.1,0.5),
                   'dstop':np.arange(-1.5,1.5+0.1,0.5),
                   }
        dcntstd = {'dpeakh':np.arange(0,3+0.1,0.5),
                   'dstop':np.arange(0,3+0.1,0.5),
                   }
        dcmave = {'fzlev':'rainbow'}
        dcmstd = {'fzlev':'rainbow'}

        dvarname = {'fzlev':'Freezing height (from surface) [km]',
                    }

        if var in ['fzlev']:
            bndave = dbndave[var]
            bndstd = dbndstd[var]
            cntave = None
            cntstd = None
            avemin,avemax=bndave[0],bndave[-1]
            stdmin,stdmax=bndstd[0],bndstd[-1]
        elif var in ['dpeakh','dstop']:
            bndave = None
            bndstd = None
            cntave = dcntave[var]
            cntstd = dcntstd[var]
            avemin,avemax=cntave[0],cntave[-1]
            stdmin,stdmax=cntstd[0],cntstd[-1]

        a2ave = np.load(odir + '/ave.%s.rad.%s.npy'%(var,season))
        a2std = np.load(odir + '/std.%s.rad.%s.npy'%(var,season))
        a2num = np.load(odir + '/num.%s.rad.%s.npy'%(var,season))

        a2ave = ma.masked_invalid(a2ave)
        a2std = ma.masked_invalid(a2std)
        print season,a2ave.min(), a2ave.max()
        print season,a2std.min(), a2std.max()
        #-- Ave ---
        figPath = figDir + '/map.%s.rad.ave.%s.png'%(var,season)
        mycm = dcmave[var]
        vmin,vmax = avemin,avemax
        stitle = '%s %s'%(dvarname[var], season)
        draw_map(a2dat=a2ave, figPath=figPath, bounds=bndave, centers=cntave)


        figPath = figDir + '/map.%s.rad.std.%s.png'%(var,season)
        mycm = dcmstd[var]
        vmin,vmax = stdmin,stdmax
        stitle = '%s STD %s'%(dvarname[var], season)
        draw_map(a2dat=a2std, figPath=figPath, bounds=bndstd, centers=cntstd)


# %%
