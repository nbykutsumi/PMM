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
workbaseDir = '/home/utsumi/mnt/lab_work'
tankbaseDir = '/home/utsumi/mnt/lab_tank'
figDir  = '/home/utsumi/temp/ret'

useorblist = True

calcflag = True
#calcflag = False
figflag  = True
#figflag  = False

iDTime = datetime(2014,6,1)
eDTime = datetime(2015,5,30)
#eDTime = datetime(2014,6,10)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
lskipdates = [[2014,8,29],[2014,9,16],[2014,10,1],[2014,10,2],[2014,11,5],[2014,12,8],[2014,12,9],[2014,12,10]]
nsample = 1000
#nsample = 20
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
#lrettype= ['epc']
lrettype= ['epc','gprof-shift']
#lvar = ['cc','dpeakh','dpeakv','dstop','dcond','rmse']+['peakhrad','peakhpmw','peakvrad','peakvpmw','stoprad','stoppmw','condrad','condpmw'] + ['dvfracconv','vfracconvrad','vfracconvpmw']

#lvar = ['condpmw','condrad','dcond','stoppmw','stoprad','dstop','cc']
lvar = ['dstop','stoprad','stoppmw']
dexpr = {'epc':'glb.relsurf01.minrec1000.maxrec10000','gprof-shift':'v01'}

lvar_pmwptype =['stoppmw','condpmw']
lvar_radptype =['stoprad','condrad','dstop','dcond','cc']

#*******************
# oid list
#*******************
if useorblist is True:
    orblistPath= tankbaseDir + '/utsumi/PMM/retepc/list-orbit/list.GMI.well.201406-201505.1000obts.txt'
    #orblistPath= tankbaseDir + '/utsumi/PMM/retepc/list-orbit/list.GMI.well.201406-201505.2obts.txt'

    f=open(orblistPath,'r'); lines=f.readlines(); f.close()
    loid = []
    for gmiPath in lines:
        gmiPath = gmiPath.strip()
        oid = int(os.path.basename(gmiPath).split('.')[-3])
        Year,Mon,Day= map(int, os.path.dirname(gmiPath).split('/')[-3:])
        if [Year,Mon,Day] in lskipdates:
            continue
        DTimeTmp = datetime(Year,Mon,Day)
        if DTimeTmp < iDTime: continue
        if DTimeTmp > eDTime: continue
        loid.append([oid,Year,Mon,Day])


else:
    loid = []
    lYM = util.ret_lYM(iYM,eYM)
    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]
        if [Year,Mon,Day] in lskipdates:
            continue
        tempbaseDir  = workbaseDir + '/hk02/PMM/NASA/GPM.GMI/1C/V05'
        lpath = sorted(glob.glob(tempbaseDir + '/%04d/%02d/%02d/1C.GPM.GMI.XCAL2016-C.*.??????.????.HDF5'%(Year,Mon,Day)))

        if len(lpath) == 0:
            continue
        for tmppath in lpath:
            oid = int(os.path.basename(int(tmppath)).split('.')[-4])
            loid.append([oid,Year,Mon,Day])

random.seed(0)
print loid
print len(loid)
a1idxtmp = range(len(loid))
a1idxtmp = sorted(random.sample(a1idxtmp, min(nsample, len(loid))))
loid = (np.array(loid)[a1idxtmp]).tolist()
#*******************
# Functions
#*******************
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
    ax  = fig.add_axes([0.08,0.1,0.8,0.7])
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

    plt.title(stitle, fontsize=13, pad=11)
    M.drawcoastlines(linewidth=1)

    M.drawmeridians(np.arange(-180,180+1,30), labels=[0,0,0,1], fontsize=10, linewidth=0.5, fmt='%d',rotation=50, yoffset=7)
    M.drawparallels(np.arange(-60,60+1,30), labels=[1,0,0,0], fontsize=10, linewidth=0.5, fmt='%d')

    cax = fig.add_axes([0.89,0.2,0.02,0.50])
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


    #-- Initialize --
    o1lat  = deque([])
    o1lon  = deque([])

    o1vfracrad = deque([])
    o1vfracpmw = deque([])

    for (oid,Year,Mon,Day) in loid:
        print Year,Mon,Day,oid
        srcDir = pairDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        a1lattmp  = np.load(srcDir + '/Latitude.%06d.npy'%(oid))
        a1lontmp  = np.load(srcDir + '/Longitude.%06d.npy'%(oid))
        a1elev    = np.load(srcDir + '/surfaceElevationrad.%06d.npy'%(oid))

        a1precpmw = np.load(srcDir + '/precpmw.%06d.npy'%(oid))
        a1precrad = np.load(srcDir + '/precrad.%06d.npy'%(oid))
        a2profpmw = np.load(srcDir + '/profpmw.%06d.npy'%(oid))[:,:nz]  # nz(500m) layers, bottom to top
        a2profrad = np.load(srcDir + '/profrad.%06d.npy'%(oid))[:,:nz]  # nz(500m) layers, bottom to top

        a1vfracrad= np.load(srcDir + '/vfracConvrad.%06d.npy'%(oid))
        a1vfracpmw= np.load(srcDir + '/vfracConvpmw.%06d.npy'%(oid))

        if rettype in ['gprof-shift','gprof']:
            a1qflag = np.load(srcDir + '/qualityFlag.%06d.npy'%(oid))


        #-- Set missing to the lower levels -------
        a1idx = np.arange(len(a1precpmw))
        vres = 500  # m
        a1sfcbin = (a1elev /vres).astype('int16')
        a1cltbin = a1sfcbin + 2  # surface + 1km
        a1cltbin = ma.masked_less(a1cltbin,0).filled(0)
        a1cltbin = ma.masked_greater(a1cltbin,nz-1).filled(nz-1)
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


        a2profradtmp = ma.masked_less(a2profrad[a1idx],0)
        a2profpmwtmp = ma.masked_less(a2profpmw[a1idx],0)
        a1lattmp     = a1lattmp[a1idx]
        a1lontmp     = a1lontmp[a1idx]
        a1elevtmp    = a1elev[a1idx]

        a1vfracradtmp= a1vfracrad[a1idx]
        a1vfracpmwtmp= a1vfracpmw[a1idx]

        o1lat.extend(list(a1lattmp)) 
        o1lon.extend(list(a1lontmp)) 

        #-- mask missing values ----
        a2mask = ma.masked_less(a2profradtmp,0).mask
        a2mask = a2mask + ma.masked_less(a2profpmwtmp,0).mask
        a2profradtmp = ma.masked_where(a2mask, a2profradtmp)
        a2profpmwtmp = ma.masked_where(a2mask, a2profpmwtmp)
        nltmp = len(a1idx)

        #-- convective volume fraction --
        o1vfracrad.extend(list(a1vfracradtmp))
        o1vfracpmw.extend(list(a1vfracpmwtmp))

    #-- Project over map ----
    lonbnd = np.arange(-180,180+0.1, dlatlon)
    latbnd = np.arange(-60,60+0.1, dlatlon)

    a1convrad = ma.masked_greater(o1vfracrad,0.5).mask.astype('int32')
    a1convpmw = ma.masked_greater(o1vfracpmw,0.5).mask.astype('int32')

    a1lat = np.array(o1lat)
    a1lon = np.array(o1lon)

    a2convrad = scipy.stats.binned_statistic_2d(a1lat,a1lon,a1convrad,statistic='sum',bins=[latbnd,lonbnd]).statistic
    a2totrad  = scipy.stats.binned_statistic_2d(a1lat,a1lon,a1convrad,statistic='count',bins=[latbnd,lonbnd]).statistic
    a2fracrad= ma.masked_where(a2totrad==0, a2convrad).astype('float32')/a2totrad.astype('float32')

    a2convpmw = scipy.stats.binned_statistic_2d(a1lat,a1lon,a1convpmw,statistic='sum',bins=[latbnd,lonbnd]).statistic
    a2totpmw  = scipy.stats.binned_statistic_2d(a1lat,a1lon,a1convpmw,statistic='count',bins=[latbnd,lonbnd]).statistic
    a2fracpmw= ma.masked_where(a2totpmw==0, a2convpmw).astype('float32')/a2totpmw.astype('float32')

    a2fracrad = a2fracrad.filled(-9999.)
    a2fracpmw = a2fracpmw.filled(-9999.)

    odir = tankbaseDir + '/utsumi/PMM/validprof/map-orbit/%s.%s'%(rettype,expr)
    util.mk_dir(odir)
    radpth = odir + '/numconvfrac.rad.npy'
    pmwpth = odir + '/numconvfrac.pmw.npy'
    np.save(radpth, a2fracrad.astype('float32'))
    np.save(pmwpth, a2fracpmw.astype('float32'))

    print radpth

#**************************************
# Figure
#--------------------------------------
for rettype in lrettype:
    if figflag is not True: continue

    expr = dexpr[rettype]
    srcdir = tankbaseDir + '/utsumi/PMM/validprof/map-orbit/%s.%s'%(rettype,expr)

    bnd = np.arange(0,1.0+0.01,0.2)

    radpth = odir + '/numconvfrac.rad.npy'
    pmwpth = odir + '/numconvfrac.pmw.npy'

    a2rad = np.load(radpth)
    a2pmw = np.load(pmwpth)

    a2rad = ma.masked_invalid(ma.masked_equal(a2rad,-9999.))
    a2pmw = ma.masked_invalid(ma.masked_equal(a2pmw,-9999.))

    #-- CMB ---
    a2dat   = a2rad
    figPath = figDir + '/map.numconvfrac.%s.rad.png'%(rettype)
    stitle = 'convective event fraction (CMB)'
    mycm = 'gist_stern_r'
    vmin,vmax = 0, 1
    draw_map(a2dat=a2dat, figPath=figPath, bounds=bnd, centers=None, extend=None)

    #-- PMW ---
    a2dat   = a2pmw
    figPath = figDir + '/map.numconvfrac.%s.pmw.png'%(rettype)
    stitle = 'convective event fraction (%s)'%{'epc':'EPC', 'gprof-shift':'GPROF'}[rettype]
    mycm = 'gist_stern_r'
    vmin,vmax = 0, 1
    draw_map(a2dat=a2dat, figPath=figPath, bounds=bnd, centers=None, extend=None)



# %%
