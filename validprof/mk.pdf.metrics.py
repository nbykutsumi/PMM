# %%
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import ListedColormap
#%matplotlib inline
from scipy.stats import binned_statistic
from numpy import *
import myfunc.util as util
from datetime import datetime, timedelta
import glob, os, sys
import numpy as np
import calendar
import random
from collections import deque
import pickle
workbaseDir = '/home/utsumi/mnt/lab_work'
tankbaseDir = '/home/utsumi/mnt/lab_tank'
figDir  = '/home/utsumi/temp/ret'

useorblist = True

#calcflag = True
calcflag = False
figflag  = True
#figflag  = False

iDTime = datetime(2014,6,1)
eDTime = datetime(2015,5,31)
#eDTime = datetime(2014,6,5)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
lskipdates = [[2014,8,29],[2014,9,16],[2014,10,1],[2014,10,2],[2014,11,5],[2014,12,8],[2014,12,9],[2014,12,10]]
nsample = 1000
#nsample = 100
#nsample = 20
lat0   = -60
lon0   = -180
dlatlon= 2.5
#nz     = 25  # 500m layers
nz     = 20  # 500m layers
ny,nx  = 120,360
miss   = -9999.
thpr = 0.5  # mm/h
thwat = 0.033  # g/m3 for storm top
#lrettype= ['epc']
lrettype= ['epc','gprof-shift']
#lvar = ['cc','dpeakh','dpeakv','dstop','dcond','rmse']+['peakhrad','peakhpmw','peakvrad','peakvpmw','stoprad','stoppmw','condrad','condpmw'] + ['dvfracconv','vfracconvrad','vfracconvpmw']

#lvar = ['cond','stop','vfracconv']
lvar = ['cond']
dexpr = {'epc':'glb.relsurf01.minrec1000.maxrec10000','gprof-shift':'v01'}


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
print len(loid)
a1idxtmp = range(len(loid))
a1idxtmp = sorted(random.sample(a1idxtmp, min(nsample, len(loid))))
loid = (np.array(loid)[a1idxtmp]).tolist()
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
    o1dpeakh = deque([])
    o1dpeakv = deque([])
    o1dstop  = deque([])
    o1dcond  = deque([])
    o1dvfrac = deque([])

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

    #-- Calc ave and variablity ----
    lonbnd = np.arange(-180,180+0.1, dlatlon)
    latbnd = np.arange(-60,60+0.1, dlatlon)

    dvar = {}
    for var in lvar:
        if var=='stop':
            a1varrad=o1stoprad
            a1varpmw=o1stoppmw
        elif var=='cond':
            a1varrad=o1condrad
            a1varpmw=o1condpmw
        elif var=='vfracconv':
            a1varrad=o1vfracrad
            a1varpmw=o1vfracpmw

        a1lat = np.array(o1lat)
        a1lon = np.array(o1lon)

        #-- Classified by precipitation type ----
        lkey = [(dattype,bywhat,ptype,region)
                for dattype in ['rad','pmw']
                for bywhat  in ['rad']
                for ptype   in ['all','conv','stra']
                for region  in ['mid','tro','all']
        ]
        for (dattype,bywhat,ptype,region) in lkey:
            if (dattype=='rad')&(bywhat=='pmw'):
                continue
            print dattype,bywhat,ptype,region

            if dattype=='rad':
                a1var = a1varrad
            elif dattype=='pmw':
                a1var = a1varpmw
            else:
                print 'check dattype',dattype

            if (bywhat =='rad'):
                a1vfrac = o1vfracrad
            elif (bywhat =='pmw'):
                a1vfrac = o1vfracpmw
            else:
                continue

            if   ptype=='all':
                a1flagptype = np.array([True])
            elif ptype=='conv':
                a1flagptype = ma.masked_greater(a1vfrac, 0.5).mask
            elif ptype=='stra':
                a1flagptype = ma.masked_less(a1vfrac, 0.5).mask
            else:
                print 'check ptype',ptype

            if   region=='all':
                a1flaglat = np.array([True])
            elif region=='mid':
                a1flaglat = ma.masked_inside(np.abs(a1lat), 35,60).mask
            elif region=='tro':
                a1flaglat = ma.masked_inside(a1lat, -15,15).mask

            #--- Screening ----
            a1flag = a1flagptype * a1flaglat

            if (len(a1flag)==1)&(a1flag.sum()==1):
                a1vartmp = np.array(a1var)
            else:
                a1vartmp = np.array(a1var)[a1flag]

            a1vartmp = ma.masked_invalid(a1vartmp)
            a1vartmp = ma.masked_equal(a1vartmp,miss)

            #--- Histogram ----------
            a1bin = np.power(10, np.arange(-3,2.6+0.01,0.1))
            #a1bin = np.concatenate(np.arange(0,1,0.01)), 
            a1num, a1bin_edges = binned_statistic(x=a1vartmp, values=a1vartmp, statistic='count', bins=a1bin)[:2]
            a1bin_center = (a1bin_edges[1:] + a1bin_edges[:-1])*0.5
            a1bin_width  = (a1bin_edges[1:] - a1bin_edges[:-1])
            a1pdf        = a1num / a1num.sum() / a1bin_width
            dvar[var,dattype,bywhat,ptype,region,'num'] = a1num
            dvar[var,dattype,bywhat,ptype,region,'pdf'] = a1pdf
            dvar[var,dattype,bywhat,ptype,region,'bin_edges'] = a1bin_edges
            dvar[var,dattype,bywhat,ptype,region,'bin_center'] = a1bin_center
            dvar[var,dattype,bywhat,ptype,region,'bin_width'] = a1bin_width


    odir = tankbaseDir + '/utsumi/PMM/validprof/metrix-orbit/%s.%s'%(rettype,expr)
    util.mk_dir(odir)
    opath = odir + '/pdf.%s.pickle'%(rettype)
    with open(opath,'w') as f:
        pickle.dump(dvar, f)
    print opath


##***********************
## Figures: cPDF
##***********************
if figflag is False:
    print 'No figures'
    sys.exit()

for rettype in lrettype:
    expr = dexpr[rettype]
    srcdir = tankbaseDir + '/utsumi/PMM/validprof/metrix-orbit/%s.%s'%(rettype,expr)
    srcpath = srcdir + '/pdf.%s.pickle'%(rettype)
    with open(srcpath,'rb') as f:
        dvar = pickle.load(f) 
    if rettype =='epc':
        dvarepc = dvar
    elif rettype in ['gprof','gprof-shift']:
        dvargpr = dvar

#for var in lvar:
ptype = 'all'
bywhat= 'rad'
region= 'all'
var   = 'cond'
a1bin_center = dvarepc[var,'rad','rad' ,ptype,region,'bin_center']
a1pdfcmb = dvarepc[var,'rad','rad' ,ptype,region,'pdf']
a1pdfepc = dvarepc[var,'pmw',bywhat,ptype,region,'pdf']
a1pdfgpr = dvargpr[var,'pmw',bywhat,ptype,region,'pdf']

a1numcmb = dvarepc[var,'rad','rad' ,ptype,region,'num']
a1numepc = dvarepc[var,'pmw',bywhat,ptype,region,'num']
a1numgpr = dvargpr[var,'pmw',bywhat,ptype,region,'num']


print a1pdfgpr
print ''
print a1numgpr
print ''
print a1bin_center
a1x = a1bin_center

nrow, ncol = 1, 1
fig = plt.figure(figsize=(6,4))
ax  = fig.add_axes([0.1,0.1,0.8,0.8])
ax.plot(a1x, a1pdfcmb, '-'  ,linewidth=2,   color='k')
ax.plot(a1x, a1pdfepc, '-'  ,linewidth=1,   color='k')
ax.plot(a1x, a1pdfgpr, '--' ,linewidth=1.5,  color='k')

#plt.yscale('log')
plt.xscale('log')
#-- title --
ax.set_title('PDF of condensed water content')

#-- X,Y-limit
#ymax = 1.1
#ax.set_ylim([0,ymax])
ax.set_xlim([0.01,20])
#-- Y-axis label --
#ax.set(ylabel='probability density')


plt.tight_layout(rect=[0, 0, 1, 0.9])
plt.show()

figpath = '/home/utsumi/temp/ret/pdf.%s.png'%(var)
plt.savefig(figpath)
print figpath

