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
#thwat = 0.033  # g/m3 for storm top
lthwat = [0.05, 0.2]  # g/m3
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

for (thwat,rettype) in [(thwat,rettype)
                    for thwat in lthwat
                    for rettype in lrettype]:

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
    o1cc   = deque([])
    o1dpeakh = deque([])
    o1dpeakv = deque([])
    o1dstop  = deque([])
    o1dcond  = deque([])
    o1dvfrac = deque([])
    o1rmse   = deque([])

    o1peakhrad = deque([])
    o1peakvrad = deque([])
    o1stoprad  = deque([])
    o1condrad  = deque([])
    o1vfracrad = deque([])

    o1peakhpmw = deque([])
    o1peakvpmw = deque([])
    o1stoppmw  = deque([])
    o1condpmw  = deque([])
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

        #-- Mask profiles that has too small magnitude --
        a1idx = ma.masked_where(a2profpmw.max(axis=1) <thwat, a1idx)
        a1idx = ma.masked_where(a2profrad.max(axis=1) <thwat, a1idx)

        #------------
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

        #-- CC -----
        a1cctmp = calc_cc(a2profradtmp[:,:15], a2profpmwtmp[:,:15], axis=1)  # to 7.5km
        o1cc.extend(list(a1cctmp))

        #-- convective volume fraction --
        a1dvfractmp = ma.masked_less(a1vfracpmwtmp,0) - ma.masked_less(a1vfracradtmp,0)
        a1dvfractmp = a1dvfractmp.filled(miss)

        o1dvfrac.extend(list(a1dvfractmp))
        o1vfracrad.extend(list(a1vfracradtmp))
        o1vfracpmw.extend(list(a1vfracpmwtmp))

        #-- peakh difference -----
        a1peakhradtmp = a2profradtmp.argmax(axis=1)*0.5 - a1elevtmp*0.001 + 0.5
        a1peakhpmwtmp = a2profpmwtmp.argmax(axis=1)*0.5 - a1elevtmp*0.001 + 0.5
        a1dpeakhtmp   = (a1peakhpmwtmp - a1peakhradtmp)  # km
        o1dpeakh.extend(list(a1dpeakhtmp))
        o1peakhrad.extend(list(a1peakhradtmp))
        o1peakhpmw.extend(list(a1peakhpmwtmp))

        #-- peak condensed water content --
        a1peakhbinradtmp = a2profradtmp.argmax(axis=1)
        a1peakhbinpmwtmp = a2profpmwtmp.argmax(axis=1)

        a1peakvradtmp = a2profradtmp[range(nltmp),a1peakhbinradtmp]
        a1peakvpmwtmp = a2profpmwtmp[range(nltmp),a1peakhbinpmwtmp]
        a1dpeakvtmp= (a1peakvpmwtmp - a1peakvradtmp)/a1peakvradtmp
        o1dpeakv.extend(list(a1dpeakvtmp))
        o1peakvrad.extend(list(a1peakvradtmp))
        o1peakvpmw.extend(list(a1peakvpmwtmp))


        #-- storm top height ------
        a2iz  = np.array(range(nz)*nltmp).reshape(-1,nz)
        a1stopradtmp = ma.masked_where(a2profradtmp<thwat, a2iz).argmax(axis=1)*0.5 - a1elevtmp*0.001 + 0.5
        a1stoppmwtmp = ma.masked_where(a2profpmwtmp<thwat, a2iz).argmax(axis=1)*0.5 - a1elevtmp*0.001 + 0.5


        #a1tmp = ma.masked_where(a2profradtmp<thwat, a2iz).argmax(axis=1)*0.5 + 0.5
        #sys.exit()

        a1dstoptmp= a1stoppmwtmp - a1stopradtmp
        o1dstop.extend(list(a1dstoptmp))
        o1stoprad.extend(list(a1stopradtmp))
        o1stoppmw.extend(list(a1stoppmwtmp))

        #-- Mean condensed water content ----
        a1condradtmp = a2profradtmp.mean(axis=1)
        a1condpmwtmp = a2profpmwtmp.mean(axis=1)
        a1dcondtmp= (a1condpmwtmp - a1condradtmp)/a1condradtmp
        o1dcond.extend(list(a1dcondtmp))
        o1condrad.extend(list(a1condradtmp))
        o1condpmw.extend(list(a1condpmwtmp))

        #-- Normalized RMSE -----
        a1rmsetmp = np.sqrt((np.square(a2profradtmp-a2profpmwtmp)).mean(axis=1)) / a1condradtmp
        o1rmse.extend(list(a1rmsetmp))

    #-- Project over map ----
    lonbnd = np.arange(-180,180+0.1, dlatlon)
    latbnd = np.arange(-60,60+0.1, dlatlon)

    for var in lvar:
        if var=='cc': a1var = o1cc     
        if var=='dpeakh': a1var = o1dpeakh 
        if var=='dpeakv': a1var = o1dpeakv 
        if var=='dstop': a1var = o1dstop  
        if var=='dcond': a1var = o1dcond  
        if var=='dvfracconv': a1var = o1dvfrac
        if var=='rmse' : a1var = o1rmse   
        if var=='peakhrad':a1var=o1peakhrad
        if var=='peakhpmw':a1var=o1peakhpmw
        if var=='peakvrad':a1var=o1peakvrad
        if var=='peakvpmw':a1var=o1peakvpmw
        if var=='stoprad' :a1var=o1stoprad 
        if var=='stoppmw' :a1var=o1stoppmw
        if var=='condrad' :a1var=o1condrad
        if var=='condpmw' :a1var=o1condpmw
        if var=='vfracconvrad':a1var=o1vfracrad
        if var=='vfracconvpmw':a1var=o1vfracpmw

        a1lat = np.array(o1lat)
        a1lon = np.array(o1lon)

        a2ave = scipy.stats.binned_statistic_2d(a1lat,a1lon,a1var,statistic='mean',bins=[latbnd,lonbnd]).statistic
        a2std = scipy.stats.binned_statistic_2d(a1lat,a1lon,a1var,statistic='std',bins=[latbnd,lonbnd]).statistic
        a2num = scipy.stats.binned_statistic_2d(a1lat,a1lon,a1var,statistic='count',bins=[latbnd,lonbnd]).statistic

        odir = tankbaseDir + '/utsumi/PMM/validprof/map-orbit/%s.%s'%(rettype,expr)
        util.mk_dir(odir)
        avepth = odir + '/ave.%s.th-%.3f.npy'%(var,thwat)
        stdpth = odir + '/std.%s.th-%.3f.npy'%(var,thwat)
        numpth = odir + '/num.%s.th-%.3f.npy'%(var,thwat)
        np.save(avepth, a2ave.astype('float32'))
        np.save(stdpth, a2std.astype('float32'))
        np.save(numpth, a2num.astype('int32'))
        print rettype,a2ave.shape
        print avepth

        ##-- Classified by self-precipitation type ----
        #for bywhat in ['rad','pmw']:
        #    if (bywhat =='rad')&(var in lvar_radptype):
        #        a1vfrac = o1vfracrad

        #    elif (bywhat =='pmw')&(var in lvar_pmwptype):
        #        a1vfrac = o1vfracpmw

        #    else:
        #        continue

        #    for ptype in ['conv','stra']:
        #        if ptype=='conv':
        #            a1flag = ma.masked_greater(a1vfrac, 0.5).mask
        #        else:
        #            a1flag = ma.masked_less(a1vfrac, 0.5).mask

        #        a1lattmp = np.array(o1lat)[a1flag]
        #        a1lontmp = np.array(o1lon)[a1flag]
        #        a1vartmp = np.array(a1var)[a1flag]

        #        a2ave = scipy.stats.binned_statistic_2d(a1lattmp,a1lontmp,a1vartmp,statistic='mean',bins=[latbnd,lonbnd]).statistic
        #        a2std = scipy.stats.binned_statistic_2d(a1lattmp,a1lontmp,a1vartmp,statistic='std',bins=[latbnd,lonbnd]).statistic
        #        a2num = scipy.stats.binned_statistic_2d(a1lattmp,a1lontmp,a1vartmp,statistic='count',bins=[latbnd,lonbnd]).statistic

        #        odir = tankbaseDir + '/utsumi/PMM/validprof/map-orbit/%s.%s'%(rettype,expr)
        #        util.mk_dir(odir)
        #        avepth = odir + '/ave.%s.by-%s-%s.th-%.3f.npy'%(var,bywhat,ptype,thwat)
        #        stdpth = odir + '/std.%s.by-%s-%s.th-%.3f.npy'%(var,bywhat,ptype,thwat)
        #        numpth = odir + '/num.%s.by-%s-%s.th-%.3f.npy'%(var,bywhat,ptype,thwat)
        #        np.save(avepth, a2ave.astype('float32'))
        #        np.save(stdpth, a2std.astype('float32'))
        #        np.save(numpth, a2num.astype('int32'))
        #        print rettype,a2ave.shape
        #        print avepth

        #----------------------------------------


        #plt.imshow(ma.masked_invalid(a2ave), origin='lower'); plt.colorbar()
        #plt.show() 

#**************************************
# Figure
#--------------------------------------
#lkey = [(rettype,var) for rettype in lrettype
#                      for var     in lvar]

lkey = [(thwat,rettype,var,bywhat,ptype)
                    for thwat   in lthwat
                    for rettype in lrettype
                    for var     in lvar
                    #for bywhat  in ['','rad','pmw']
                    #for ptype   in ['','conv','stra']
                    for bywhat  in ['']
                    for ptype   in ['']
                    ]


for key in lkey:
    thwat,rettype,var,bywhat,ptype = key
    if figflag is not True: continue

    if bywhat=='':
        if ptype !='': continue
    if bywhat=='rad':
        if var not in lvar_radptype: continue
    if bywhat=='pmw':
        if var not in lvar_pmwptype: continue

    expr = dexpr[rettype]
    srcdir = tankbaseDir + '/utsumi/PMM/validprof/map-orbit/%s.%s'%(rettype,expr)

    dbndave = { 'cc':np.arange(0,1+0.01,0.1),
                'dpeakh': [-3,-2,-1,-0.5,0.5,1,2,3],
                'dpeakv': [-9999,-0.5,-0.2,0.2,0.5,9999],
                'dcond': [-9999,-0.5,-0.2,0.2,0.5,9999],
                'dstop': [-2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2],
                'dvfracconv':[-9999,-0.5,-0.3,-0.1,0.1,0.3,0.5,9999],
                'rmse':np.arange(0,1.6+0.01,0.4),
                'peakhrad':[0,2,3,4,5,6,7,8],
                'peakhpmw':[0,2,3,4,5,6,7,8],
                'peakvrad':np.arange(0,0.8+0.01,0.1),
                'peakvpmw':np.arange(0,0.8+0.01,0.1),
                'stoprad' :np.arange(0,12+0.1,2),
                'stoppmw' :np.arange(0,12+0.1,2),
                'condrad' :np.arange(0,0.3+0.01, 0.05),
                'condpmw' :np.arange(0,0.3+0.01, 0.05),
                'vfracconvrad':np.arange(0,1.0+0.01,0.2),
                'vfracconvpmw':np.arange(0,1.0+0.01,0.2),
               }
    dbndstd = { 'cc':np.arange(0,1+0.01,0.1),
                'dpeakh':np.arange(0,2.5+0.01,0.5),
                'dpeakv':np.arange(0,1+0.01,0.2),
                'dcond':np.arange(0,1+0.01,0.2),
                'dstop':np.arange(0,2.5+0.01,0.5),
                'dvfracconv':np.arange(0,0.5+0.01,0.1),
                'rmse':np.arange(0,1+0.01,0.2),
                'peakhrad':np.arange(0,2.5+0.01,0.5),
                'peakhpmw':np.arange(0,2.5+0.01,0.5),
                'peakvrad':np.arange(0,0.8+0.01,0.1),
                'peakvpmw':np.arange(0,0.8+0.01,0.1),
                'stoprad' :np.arange(0,5+0.1,1),
                'stoppmw' :np.arange(0,5+0.1,1),
                'condrad' :np.arange(0,0.3+0.01, 0.05),
                'condpmw' :np.arange(0,0.3+0.01, 0.05),
                'vfracconvrad':np.arange(0,0.5+0.01,0.1),
                'vfracconvpmw':np.arange(0,0.5+0.01,0.1),
                }
    #dcntave = { 'dpeakh':np.arange(-1.5,1.5+0.1,0.5),
    #            'dstop':np.arange(-1.5,1.5+0.1,0.5),
    #            }
    #dcntstd = { 'dpeakh':np.arange(0,3+0.1,0.5),
    #            'dstop':np.arange(0,3+0.1,0.5),
    #            }

    dcmave = {}
    dcmstd = {}
    if var in ['cc','rmse']:
        dcmave[var] = 'rainbow'
        dcmstd[var] = 'hot_r'
    elif var in ['dpeakh','dpeakv','dstop','dcond','dvfracconv']:
        dcmave[var] = 'RdBu_r'
        dcmstd[var] = 'hot_r'
    elif var in ['peakhrad','peakhpmw','stoprad','stoppmw']:
        dcmave[var] = 'rainbow'
        dcmstd[var] = 'hot_r'
    elif var in ['peakvrad','peakvpmw','condrad','condpmw','vfracconvrad','vfracconvpmw']:
        dcmave[var] = 'gist_stern_r'
        dcmstd[var] = 'hot_r'

    dretname = {'epc':'EPC','gprof-shift':'GPROF','gprof':'GPROF'}
    dvarname = {'cc':'corr. coef',
                'rmse':'Normalized RMSE',
                'dpeakh':'Diff. of the height of the peak [km]',
                'peakhrad':'Height of the condensed water peak [km]',
                'peakhpmw':'Height of the condensed water peak [km]',
                'dpeakv':'Peak condened water cont. diff. (Normed)',
                'peakvrad':'Peak condensed water cont. [g/m3]',
                'peakvpmw':'Peak condensed water cont. [g/m3]',
                'dstop':'Storm top height diff. [km]',
                'stoprad':'Storm top height [km]',
                'stoppmw':'Storm top height [km]',
                'dcond':'Mean condensed water diff (Normed)',
                'condrad':'Mean condensed water cont. [g/m3]',
                'condpmw':'Mean condensed water cont. [g/m3]',
                'dvfracconv':'Convective vol. frac. diff',
                'vfracconvrad':'Mean convective vol. frac.',
                'vfracconvpmw':'Mean convective vol. frac.',
                }

    if var in ['cc','dstop','dpeakh','dpeakv','dcond','rmse','dvfracconv']+['peakhrad','peakhpmw','peakvrad','peakvpmw','stoprad','stoppmw','condrad','condpmw','vfracconvrad','vfracconvpmw']:
        bndave = dbndave[var]
        bndstd = dbndstd[var]
        cntave = None
        cntstd = None
        avemin,avemax=bndave[0],bndave[-1]
        stdmin,stdmax=bndstd[0],bndstd[-1]
    #elif var in ['dpeakh','dstop']:
    #    bndave = None
    #    bndstd = None
    #    cntave = dcntave[var]
    #    cntstd = dcntstd[var]
    #    avemin,avemax=cntave[0],cntave[-1]
    #    stdmin,stdmax=cntstd[0],cntstd[-1]

    if var in ['dpeakv','dcond','dvfracconv']:
        extend = 'both'
    elif var in ['peakhrad','peakhpmw']:
        extend = 'min'
    else:
        extend = None


    if ptype=='':
        a2ave = np.load(srcdir + '/ave.%s.th-%.3f.npy'%(var,thwat))
        a2std = np.load(srcdir + '/std.%s.th-%.3f.npy'%(var,thwat))
        a2num = np.load(srcdir + '/num.%s.th-%.3f.npy'%(var,thwat))
    else:
        a2ave = np.load(srcdir + '/ave.%s.by-%s-%s.th-%.3f.npy'%(var,bywhat,ptype,thwat))
        a2std = np.load(srcdir + '/std.%s.by-%s-%s.th-%.3f.npy'%(var,bywhat,ptype,thwat))
        a2num = np.load(srcdir + '/num.%s.by-%s-%s.th-%.3f.npy'%(var,bywhat,ptype,thwat))


    a2ave = ma.masked_invalid(a2ave)
    a2std = ma.masked_invalid(a2std)

    if var[-3:]=='rad':
        datname='CMB'
    elif rettype=='epc':
        datname='EPC'
    elif rettype in ['gprof','gprof-shift']:
        datname='GPROF'
    else:
        datname='XXX'

    #-- Ave ---
    if ptype=='':
        figPath = figDir + '/map.%s.%s.th-%.3f.ave.png'%(rettype,var,thwat)
        stitle = '%s %s %.3fg/m3'%(datname,dvarname[var],thwat)
    else:
        figPath = figDir + '/map.%s.%s.by-%s-%s.th-%.3f.ave.png'%(rettype,var,bywhat,ptype,thwat)
        stitle = '%s %s (%s-by-%s) %.3fg/m3'%(datname,dvarname[var], str.upper(ptype),str.upper(bywhat), thwat)

    mycm = dcmave[var]
    vmin,vmax = avemin,avemax
    draw_map(a2dat=a2ave, figPath=figPath, bounds=bndave, centers=cntave, extend=extend)


    #-- STD ---
    if ptype !='': continue
    if ptype =='':
        figPath = figDir + '/map.%s.%s.th-%.3f.std.png'%(rettype,var,thwat)
        stitle = '%s %s %.3fg/m3 STD'%(datname,dvarname[var],thwat)
    else:
        figPath = figDir + '/map.%s.%s.by-%s-%s.th-%.3f.std.png'%(rettype,var,bywhat,ptype,thwat)
        stitle = '%s %s (%s-by-%s) %.3fg/m3 STD'%(datname,dvarname[var], str.upper(ptype),str.upper(bywhat),thwat)

    mycm = dcmstd[var]
    vmin,vmax = stdmin,stdmax
    draw_map(a2dat=a2std, figPath=figPath, bounds=bndstd, centers=cntstd)


# %%
