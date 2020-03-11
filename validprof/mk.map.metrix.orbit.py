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

iYM = [2014,6]
eYM = [2015,5]
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
#thpr = 3.0  # mm/h
thwat = 0.033  # g/m3 for storm top
lrettype= ['epc']
#lrettype= ['gprof-shift']
lvar = ['cc','dpeakh','dpeakv','dstop','dcond','rmse']+['peakhrad','peakhpmw','peakvrad','peakvpmw','stoprad','stoppmw','condrad','condpmw']
#lvar = ['dcond']
dexpr = {'epc':'glb.relsurf01.minrec1000.maxrec10000','gprof-shift':'v01'}


def calc_cc(x,y,axis):
    ''' masked elements are not used '''

    if len(x.shape)==3:
        ny,nx,nz = x.shape
        xm = x.mean(axis=axis).reshape(ny,nx,1)
        ym = y.mean(axis=axis).reshape(ny,nx,1)
    elif len(x.shape)==2:
        ny,nz = x.shape
        xm = x.mean(axis=axis).reshape(ny,1)
        ym = y.mean(axis=axis).reshape(ny,1)

    A  = ((x-xm)*(y-ym)).sum(axis=axis)
    B  = ((x-xm)**2).sum(axis=axis)
    C  = ((y-ym)**2).sum(axis=axis)
    return A/( np.sqrt(B*C) )

def ret_stormtop(a2prof, miss_out=-9999.):
    thmin = 0.1 # g/m3
    ny,nz = a2prof.shape
    a1h = arange(nz).reshape(1,-1)*0.5 + 0.25
    a2h = np.repeat(a1h, ny, axis=0)
    a1stop = ma.masked_where(a2prof< thmin, a2h).max(axis=1).filled(miss_out)
    return a1stop 

def ret_peakheight(a2prof, miss_out=-9999.):
    thmin = 0.1 # g/m3
    a2prof = ma.masked_less(a2prof, thmin)
    a1mask = a2prof.mask.all(axis=1)
    a1h = a2prof.argmax(axis=1)*0.5 + 0.25
    a1h = ma.masked_where(a1mask, a1h).filled(miss_out)
    return a1h

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

for rettype in lrettype:
    if calcflag is not True: continue
    #**** File list ********************
    expr = dexpr[rettype]
    if rettype=='epc': 
        pairDir = tankbaseDir + '/utsumi/PMM/validprof/pair/%s.%s'%(rettype, expr)
    elif rettype=='gprof-shift':
        pairDir = tankbaseDir + '/utsumi/PMM/validprof/pair/%s'%(rettype)
    else:
        print 'check rettype',rettype
        sys.exit()

    llatPath = []
    for Year,Mon in lYM:
        print Year,Mon
        srcDir = pairDir + '/%04d/%02d/??'%(Year,Mon)
        ssearch= srcDir + '/Latitude.*.npy'
        llatPathTmp= glob.glob(ssearch)
        llatPath = llatPath + llatPathTmp

    random.seed(0)
    llatPath = np.sort(random.sample(llatPath, nsample))


    #-- Initialize --
    a1lat  = deque([])
    a1lon  = deque([])
    a1cc   = deque([])
    a1dpeakh = deque([])
    a1dpeakv = deque([])
    a1dstop  = deque([])
    a1dcond  = deque([])
    a1rmse   = deque([])

    a1peakhrad = deque([])
    a1peakvrad = deque([])
    a1stoprad  = deque([])
    a1condrad  = deque([])

    a1peakhpmw = deque([])
    a1peakvpmw = deque([])
    a1stoppmw  = deque([])
    a1condpmw  = deque([])


    for ipath,latPath in enumerate(llatPath):
        srcDir = os.path.dirname(latPath)
        oid    = int(latPath.split('.')[-2])
        a1lattmp  = np.load(srcDir + '/Latitude.%06d.npy'%(oid))
        a1lontmp  = np.load(srcDir + '/Longitude.%06d.npy'%(oid))
        a1elev    = np.load(srcDir + '/surfaceElevationrad.%06d.npy'%(oid))

        a1precpmw = np.load(srcDir + '/precpmw.%06d.npy'%(oid))
        a1precrad = np.load(srcDir + '/precrad.%06d.npy'%(oid))
        a2profpmw = np.load(srcDir + '/profpmw.%06d.npy'%(oid))[:,:nz]  # nz(500m) layers, bottom to top
        a2profrad = np.load(srcDir + '/profrad.%06d.npy'%(oid))[:,:nz]  # nz(500m) layers, bottom to top

        if rettype in ['gprof-shift','gprof']:
            a1qflag = np.load(srcDir + '/qualityFlag.%06d.npy'%(oid))


        #-- Set missing to the lower levels -------
        a1idx = np.arange(len(a1precpmw))
        vres = 500  # m
        a1sfcbin = (a1elev /vres).astype('int16')
        a1cltbin = a1sfcbin + 2  # surface + 1km
        a1cltbin = ma.masked_less(a1cltbin,0).filled(0)
        a1cltbin = ma.masked_greater(a1cltbin,24).filled(24)
        setcltbin= set(list(a1cltbin))

        for cltbin in setcltbin:
            a1idxtmp = ma.masked_where(a1cltbin !=cltbin, a1idx).compressed()
            a2profpmw[a1idxtmp, :cltbin] = miss
            a2profrad[a1idxtmp, :cltbin] = miss
        #-------------------------------------------

        a1idx = ma.masked_where(a1precpmw<thpr, a1idx)
        a1idx = ma.masked_where(a1precrad<thpr, a1idx)
        if rettype in ['gprof','gprof-shift']:
            a1idx = ma.masked_where(a1qflag !=0, a1idx)

        a1idx = a1idx.compressed()


        print ipath,len(a1idx), a1precpmw.shape[0]
        print latPath
        #a1idx = random.sample(a1idx,30)  # test

        a2profradtmp = ma.masked_less(a2profrad[a1idx],0)
        a2profpmwtmp = ma.masked_less(a2profpmw[a1idx],0)
        a1lattmp     = a1lattmp[a1idx]
        a1lontmp     = a1lontmp[a1idx]

        a1lat.extend(list(a1lattmp)) 
        a1lon.extend(list(a1lontmp)) 

        #-- mask missing values ----
        a2mask = ma.masked_less(a2profradtmp,0).mask
        a2mask = a2mask + ma.masked_less(a2profpmwtmp,0).mask
        a2profradtmp = ma.masked_where(a2mask, a2profradtmp)
        a2profpmwtmp = ma.masked_where(a2mask, a2profpmwtmp)
        nltmp = len(a1idx)

        #-- CC -----
        a1cctmp = calc_cc(a2profradtmp[:,:15], a2profpmwtmp[:,:15], axis=1)  # to 7.5km
        a1cc.extend(list(a1cctmp))

        #-- peakh difference -----
        a1peakhradtmp = a2profradtmp.argmax(axis=1)
        a1peakhpmwtmp = a2profpmwtmp.argmax(axis=1)
        a1dpeakhtmp   = (a1peakhpmwtmp - a1peakhradtmp)*0.5  # km
        a1dpeakh.extend(list(a1dpeakhtmp))
        a1peakhrad.extend(list(a1peakhradtmp*0.5))
        a1peakhpmw.extend(list(a1peakhpmwtmp*0.5))

        #-- peak condensed water content --
        a1peakvradtmp = a2profradtmp[range(nltmp),a1peakhradtmp]
        a1peakvpmwtmp = a2profpmwtmp[range(nltmp),a1peakhpmwtmp]
        a1dpeakvtmp= (a1peakvpmwtmp - a1peakvradtmp)/a1peakvradtmp
        a1dpeakv.extend(list(a1dpeakvtmp))
        a1peakvrad.extend(list(a1peakvradtmp))
        a1peakvpmw.extend(list(a1peakvpmwtmp))


        #-- storm top height ------
        a2iz  = np.array(range(nz)*nltmp).reshape(-1,nz)
        a1stopradtmp = ma.masked_where(a2profradtmp<thwat, a2iz).argmax(axis=1)*0.5
        a1stoppmwtmp = ma.masked_where(a2profpmwtmp<thwat, a2iz).argmax(axis=1)*0.5
        a1dstoptmp= a1stoppmwtmp - a1stopradtmp
        a1dstop.extend(list(a1dstoptmp))
        a1stoprad.extend(list(a1stopradtmp))
        a1stoppmw.extend(list(a1stoppmwtmp))

        #-- Mean condensed water content ----
        a1condradtmp = a2profradtmp.mean(axis=1)
        a1condpmwtmp = a2profpmwtmp.mean(axis=1)
        a1dcondtmp= (a1condpmwtmp - a1condradtmp)/a1condradtmp
        a1dcond.extend(list(a1dcondtmp))
        a1condrad.extend(list(a1condradtmp))
        a1condpmw.extend(list(a1condpmwtmp))

        #-- Normalized RMSE -----
        a1rmsetmp = np.sqrt((np.square(a2profradtmp-a2profpmwtmp)).mean(axis=1)) / a1condradtmp
        a1rmse.extend(list(a1rmsetmp))

    #-- Project over map ----
    lonbnd = np.arange(-180,180+0.1, dlatlon)
    latbnd = np.arange(-60,60+0.1, dlatlon)

    for var in lvar:
        if var=='cc': a1var = a1cc     
        if var=='dpeakh': a1var = a1dpeakh 
        if var=='dpeakv': a1var = a1dpeakv 
        if var=='dstop': a1var = a1dstop  
        if var=='dcond': a1var = a1dcond  
        if var=='rmse' : a1var = a1rmse   
        if var=='peakhrad':a1var=a1peakhrad
        if var=='peakhpmw':a1var=a1peakhpmw
        if var=='peakvrad':a1var=a1peakvrad
        if var=='peakvpmw':a1var=a1peakvpmw
        if var=='stoprad' :a1var=a1stoprad 
        if var=='stoppmw' :a1var=a1stoppmw
        if var=='condrad' :a1var=a1condrad
        if var=='condpmw' :a1var=a1condpmw

        print len(a1lat), len(a1lon), len(a1var)
        a2ave = scipy.stats.binned_statistic_2d(a1lat,a1lon,a1var,statistic='mean',bins=[latbnd,lonbnd]).statistic
        a2std = scipy.stats.binned_statistic_2d(a1lat,a1lon,a1var,statistic='std',bins=[latbnd,lonbnd]).statistic
        a2num = scipy.stats.binned_statistic_2d(a1lat,a1lon,a1var,statistic='count',bins=[latbnd,lonbnd]).statistic

        odir = tankbaseDir + '/utsumi/PMM/validprof/map-orbit/%s.%s'%(rettype,expr)
        util.mk_dir(odir)
        avepth = odir + '/ave.%s.npy'%(var)
        stdpth = odir + '/std.%s.npy'%(var)
        numpth = odir + '/num.%s.npy'%(var)
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
lkey = [(rettype,var) for rettype in lrettype
                      for var     in lvar]

for key in lkey:
    if figflag is not True: continue

    rettype,var = key
    expr = dexpr[rettype]
    odir = tankbaseDir + '/utsumi/PMM/validprof/map-orbit/%s.%s'%(rettype,expr)

    dbndave = {'cc':np.arange(0,1+0.01,0.1),
               'dpeakv':np.arange(-0.5,0.5+0.01,0.2),
               'dcond':np.arange(-0.5,0.5+0.01,0.2),
               'rmse':np.arange(0,1.6+0.01,0.4),
               'peakhrad':np.arange(0,8+0.1,1),
               'peakhpmw':np.arange(0,8+0.1,1),
               'peakvrad':np.arange(0,0.8+0.01,0.1),
               'peakvpmw':np.arange(0,0.8+0.01,0.1),
               'stoprad' :np.arange(0,12+0.1,2),
               'stoppmw' :np.arange(0,12+0.1,2),
               'condrad' :np.arange(0,0.3+0.01, 0.05),
               'condpmw' :np.arange(0,0.3+0.01, 0.05),

               }
    dbndstd = {'cc':np.arange(0,1+0.01,0.1),
               'dpeakv':np.arange(0,1+0.01,0.2),
               'dcond':np.arange(0,1+0.01,0.2),
               'rmse':np.arange(0,1+0.01,0.2),
               'peakhrad':np.arange(0,4+0.1,1),
               'peakhpmw':np.arange(0,4+0.1,1),
               'peakvrad':np.arange(0,0.8+0.01,0.1),
               'peakvpmw':np.arange(0,0.8+0.01,0.1),
               'stoprad' :np.arange(0,5+0.1,1),
               'stoppmw' :np.arange(0,5+0.1,1),
               'condrad' :np.arange(0,0.3+0.01, 0.05),
               'condpmw' :np.arange(0,0.3+0.01, 0.05),
               }
    dcntave = {'dpeakh':np.arange(-1.5,1.5+0.1,0.5),
               'dstop':np.arange(-1.5,1.5+0.1,0.5),
               }
    dcntstd = {'dpeakh':np.arange(0,3+0.1,0.5),
               'dstop':np.arange(0,3+0.1,0.5),
               }
    dcmave = {'cc':'rainbow','dpeakh':'RdBu_r','dpeakv':'RdBu_r','dstop':'RdBu_r','dcond':'RdBu_r','rmse':'rainbow'}
    dcmstd = {'cc':'rainbow','dpeakh':'rainbow' ,'dpeakv':'rainbow' ,'dstop':'rainbow' ,'dcond':'rainbow','rmse':'rainbow'}

    if var in ['peakhrad','peakhpmw','stoprad','stoppmw']:
        dcmave[var] = 'rainbow'
        dcmstd[var] = 'rainbow'
    elif var in ['peakvrad','peakvpmw','condrad','condpmw']:
        dcmave[var] = 'gist_stern_r'
        dcmstd[var] = 'gist_stern_r'

    dretname = {'epc':'EPC','gprof-shift':'GPROF','gprof':'GPROF'}
    dvarname = {'cc':'corr. coef',
                'rmse':'Normalized RMSE',
                'dpeakh':'Peak condensed water height diff.[km]',
                'peakhrad':'Peak condensed water height [km]',
                'peakhpmw':'Peak condensed water height [km]',
                'dpeakv':'Peak condened waater cont. diff. [g/m3]',
                'peakvrad':'Peak condensed water cont. [g/m3]',
                'peakvpmw':'Peak condensed water cont. [g/m3]',
                'dstop':'Storm top height diff. [km]',
                'stoprad':'Storm top height [km]',
                'stoppmw':'Storm top height [km]',
                'dcond':'Mean condensed water diff (Normed)',
                'condrad':'Mean condensed water cont. [g/m3]',
                'condpmw':'Mean condensed water cont. [g/m3]',
                }

    if var in ['cc','dpeakv','dcond','rmse']+['peakhrad','peakhpmw','peakvrad','peakvpmw','stoprad','stoppmw','condrad','condpmw']:
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

    a2ave = np.load(odir + '/ave.%s.npy'%(var))
    a2std = np.load(odir + '/std.%s.npy'%(var))
    a2num = np.load(odir + '/num.%s.npy'%(var))

    a2ave = ma.masked_invalid(a2ave)
    a2std = ma.masked_invalid(a2std)

    #-- Ave ---
    figPath = figDir + '/map.%s.%s.ave.png'%(rettype,var)
    mycm = dcmave[var]
    vmin,vmax = avemin,avemax
    stitle = '%s %s'%(dretname[rettype],dvarname[var])
    draw_map(a2dat=a2ave, figPath=figPath, bounds=bndave, centers=cntave)


    figPath = figDir + '/map.%s.%s.std.png'%(rettype,var)
    mycm = dcmstd[var]
    vmin,vmax = stdmin,stdmax
    stitle = '%s %s STD'%(dretname[rettype],dvarname[var])
    draw_map(a2dat=a2std, figPath=figPath, bounds=bndstd, centers=cntstd)


# %%
