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

lvar = ['cond','stop','vfracconv']
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
#*******************
# Functions
#*******************
def calc_cc(x,y,axis):
    ''' masked elements are not used '''
    x = ma.masked_invalid(ma.masked_less(x,0))
    y = ma.masked_invalid(ma.masked_less(y,0))
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
                for bywhat  in ['rad','pmw']
                for ptype   in ['all','conv','stra']
                for region  in ['mid','tro']
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
            a1vartmp = np.array(a1var)[a1flag]


            a1vartmp = ma.masked_invalid(a1vartmp)
            a1vartmp = ma.masked_equal(a1vartmp,miss)

            #--- Calc ----------
            dvar[var,dattype,bywhat,ptype,region,'ave'] = a1vartmp.mean()
            dvar[var,dattype,bywhat,ptype,region,'num'] = a1vartmp.count()
            dvar[var,dattype,bywhat,ptype,region,'p05'] = np.nanpercentile(a1vartmp.filled(np.nan),5)
            dvar[var,dattype,bywhat,ptype,region,'p25'] = np.nanpercentile(a1vartmp.filled(np.nan),25)
            dvar[var,dattype,bywhat,ptype,region,'p50'] = np.nanpercentile(a1vartmp.filled(np.nan),50)
            dvar[var,dattype,bywhat,ptype,region,'p75'] = np.nanpercentile(a1vartmp.filled(np.nan),75)
            dvar[var,dattype,bywhat,ptype,region,'p95'] = np.nanpercentile(a1vartmp.filled(np.nan),95)

            #if (var=='cond')&(ptype=='stra')&(region=='mid'):
            #    print a1vartmp.min()
            #    sys.exit()



    odir = tankbaseDir + '/utsumi/PMM/validprof/metrix-orbit/%s.%s'%(rettype,expr)
    util.mk_dir(odir)
    opath = odir + '/metrix.%s.pickle'%(rettype)
    with open(opath,'w') as f:
        pickle.dump(dvar, f)
    print opath


    #*********************************************************
    # Calc convective count fraction (with Bootstrap)
    #*********************************************************
    a1vfracrad = np.array(o1vfracrad)
    a1vfracpmw = np.array(o1vfracpmw)

    dcfrac = {}
    lkey = [(dattype,radcondition,region)
            for dattype       in ['rad','pmw']
            for radcondition  in ['all','conv','stra']
            for region        in ['all','mid','tro']
            ]
    for (dattype,radcondition,region) in lkey:
        if (dattype=='rad')&(radcondition=='conv'): continue

        if dattype=='rad':
            a1vfrac = a1vfracrad
        elif dattype=='pmw':
            a1vfrac = a1vfracpmw
        else:
            print 'check dattype',dattype
            sys.exit()

        #-- Screen ----
        if   radcondition=='all':
            a1flagradcondition= np.array([True])
        elif radcondition=='conv':
            a1flagradcondition= ma.masked_greater(a1vfracrad,0.5).mask
        elif radcondition=='stra':
            a1flagradcondition= ma.masked_less_equal(a1vfracrad,0.5).mask
        else:
            print 'check radcondition',radcondition

        if   region=='all':
            a1flaglat = np.array([True]*a1lat.shape[0])
        elif region=='mid':
            a1flaglat = ma.masked_inside(np.abs(a1lat), 35,60).mask
        elif region=='tro':
            a1flaglat = ma.masked_inside(a1lat, -15,15).mask

        a1flag = a1flagradcondition * a1flaglat

        a1conv = (ma.masked_greater(a1vfrac, 0.5).mask).astype('float32')[a1flag]

        a1conv = ma.masked_less(a1conv,0)
        a1conv = ma.masked_invalid(a1conv).compressed()
        #-- Bootstrap---
        nchoice = min(20000, a1conv.shape[0])
        nloop   = 500 
        print 'Bootstrap',dattype,region
        a1vartmp = [np.random.choice(a1conv, nchoice, replace=True).sum()/float(nchoice) for i in range(nloop)]
        print 'Bootstrap Done'

        a1vartmp = np.array(a1vartmp)
        dcfrac[dattype,radcondition,region,'ave'] = a1vartmp.mean()
        dcfrac[dattype,radcondition,region,'p05'] = np.percentile(a1vartmp,5)
        dcfrac[dattype,radcondition,region,'p25'] = np.percentile(a1vartmp,25)
        dcfrac[dattype,radcondition,region,'p50'] = np.percentile(a1vartmp,50)
        dcfrac[dattype,radcondition,region,'p75'] = np.percentile(a1vartmp,75)
        dcfrac[dattype,radcondition,region,'p95'] = np.percentile(a1vartmp,95)

    opath = odir + '/convCountFrac.%s.pickle'%(rettype)
    with open(opath,'w') as f:
        pickle.dump(dcfrac, f)
    print opath




##***********************
## Figures: condensed water content
##***********************
if figflag is False:
    print 'No figures'
    sys.exit()

for rettype in lrettype:
    expr = dexpr[rettype]
    srcdir = tankbaseDir + '/utsumi/PMM/validprof/metrix-orbit/%s.%s'%(rettype,expr)
    srcpath = srcdir + '/metrix.%s.pickle'%(rettype)
    with open(srcpath,'rb') as f:
        dvar = pickle.load(f) 
    if rettype =='epc':
        dvarepc = dvar
    elif rettype in ['gprof','gprof-shift']:
        dvargpr = dvar

#for var in lvar:
for var in ['cond']:
    for bywhat in ['rad','pmw']:
        nrow, ncol = 1, 4
        fig, axs = plt.subplots(nrows=nrow, ncols=ncol, figsize=(8,1.6))
        #lkey = [(ptype,region)
        #                    for ptype  in ['stra','conv']
        #                    for region in ['tro','mid']  # 'all','mid','tro'
        #]

        lkey = [(region,ptype)
                            for region in ['mid','tro']  # 'all','mid','tro'
                            for ptype  in ['conv','stra']
        ]


        #for ikey, (ptype,region) in enumerate(lkey):
        for ikey, (region,ptype) in enumerate(lkey):
            ax = axs[ikey]
            avecmb = dvarepc[var,'rad','rad' ,ptype,region,'ave']
            aveepc = dvarepc[var,'pmw',bywhat,ptype,region,'ave']
            avegpr = dvargpr[var,'pmw',bywhat,ptype,region,'ave']
            p05cmb = dvarepc[var,'rad','rad' ,ptype,region,'p05']
            p05epc = dvarepc[var,'pmw',bywhat,ptype,region,'p05']
            p05gpr = dvargpr[var,'pmw',bywhat,ptype,region,'p05']
            p25cmb = dvarepc[var,'rad','rad' ,ptype,region,'p25']
            p25epc = dvarepc[var,'pmw',bywhat,ptype,region,'p25']
            p25gpr = dvargpr[var,'pmw',bywhat,ptype,region,'p25']
            p75cmb = dvarepc[var,'rad','rad' ,ptype,region,'p75']
            p75epc = dvarepc[var,'pmw',bywhat,ptype,region,'p75']
            p75gpr = dvargpr[var,'pmw',bywhat,ptype,region,'p75']
            p95cmb = dvarepc[var,'rad','rad' ,ptype,region,'p95']
            p95epc = dvarepc[var,'pmw',bywhat,ptype,region,'p95']
            p95gpr = dvargpr[var,'pmw',bywhat,ptype,region,'p95']

            aerrmin= -np.array([p25cmb,p25epc,p25gpr]) + np.array([avecmb,aveepc,avegpr])
            aerrmax=  np.array([p75cmb,p75epc,p75gpr]) - np.array([avecmb,aveepc,avegpr])

            labels = ['CMB','EPC','GPR']
            rects  = ax.bar(labels, [avecmb,aveepc,avegpr], yerr=[aerrmin,aerrmax], color='gray')

            #-- title --
            ptypename = {'stra':'Strat','conv':'Conv'}[ptype]
            regionname= {'mid':'Mid-high','tro':'Tropics'}[region]
            #ax.set_title('%s (%s)'%(ptypename, regionname))
            ax.set_title('%s (%s)'%(regionname, ptypename))

            #-- Y-limit
            #ymax = {('stra','mid'):0.35,
            #        ('stra','tro'):0.35,
            #        ('conv','mid'):0.35,
            #        ('conv','tro'):0.35,
            #    }[(ptype,region)]
            ymax  = 0.27  # 0.12, 0.27 ?
            ax.set_ylim([0,ymax])

            #-- Y-axis label --
            if ikey==0:
                ax.set(ylabel='[g/m3]')

            #-- text annotation --
            for irect,rect in enumerate(rects):
                height = rect.get_height()
                #if irect ==0:                
                #    xy = (rect.get_x()*0.93, height + ymax*0.1)
                #else:
                #    xy = (rect.get_x()*0.93, height + ymax*0.3)
                xy = (rect.get_x()*0.93, height + ymax*0.1)
                ax.annotate('%.2f'%(height), xy=xy)

            ##-- text annotation(percentatge) --
            #print ''
            #print ptype,region,bywhat
            #for irect,rect in enumerate(rects):
            #    height = rect.get_height()
            #    if irect==0:
            #        height0 = height
            #        continue
            #    else:
            #        ratio = height/height0
            #    xy = (rect.get_x()*0.91, height + ymax*0.1)
            #    ax.annotate('(%.2f)'%(ratio), xy=xy)
            #    print height0, height, height/height0, ratio


        varname={'cond':'Mean condensed water content'}['cond']
        bywhatname={'rad':"Conditional on CMB's precip type",
                    'pmw':"Conditional on own precip type"}[bywhat]
        plt.suptitle('%s (%s)'%(varname, bywhatname))
        plt.tight_layout(rect=[0, 0, 1, 0.9])
        plt.show()

        bywhatfilename = {'rad':'rad','pmw':'self'}[bywhat]
        figpath = '/home/utsumi/temp/ret/bar.%s.by-%s.png'%(var,bywhatfilename)
        plt.savefig(figpath)
        print figpath


##***********************
## Figures: condensed water content (unconditional)
##***********************
for var in ['cond']:
    for bywhat in ['rad','pmw']:
        nrow, ncol = 1, 2
        ptype= 'all'

        fig, axs = plt.subplots(nrows=nrow, ncols=ncol, figsize=(4.1,1.6))
        lkey = [region for region in ['tro','mid']]
        for ikey, (region) in enumerate(lkey):
            ax = axs[ikey]
            avecmb = dvarepc[var,'rad','rad' ,ptype,region,'ave']
            aveepc = dvarepc[var,'pmw',bywhat,ptype,region,'ave']
            avegpr = dvargpr[var,'pmw',bywhat,ptype,region,'ave']
            p25cmb = dvarepc[var,'rad','rad' ,ptype,region,'p25']
            p25epc = dvarepc[var,'pmw',bywhat,ptype,region,'p25']
            p25gpr = dvargpr[var,'pmw',bywhat,ptype,region,'p25']
            p75cmb = dvarepc[var,'rad','rad' ,ptype,region,'p75']
            p75epc = dvarepc[var,'pmw',bywhat,ptype,region,'p75']
            p75gpr = dvargpr[var,'pmw',bywhat,ptype,region,'p75']

            aerrmin= -np.array([p25cmb,p25epc,p25gpr]) + np.array([avecmb,aveepc,avegpr])
            aerrmax=  np.array([p75cmb,p75epc,p75gpr]) - np.array([avecmb,aveepc,avegpr])

            labels = ['CMB','EPC','GPR']
            rects  = ax.bar(labels, [avecmb,aveepc,avegpr], yerr=[aerrmin,aerrmax], color='gray')

            #-- title --
            regionname= {'mid':'Mid-high','tro':'Tropics'}[region]
            ax.set_title('%s '%(regionname))

            #-- Y-limit
            #ymax  = 0.35  # 0.12, 0.27 ?
            ymax  = 0.27  # 0.12, 0.27 ?
            ax.set_ylim([0,ymax])

            #-- Y-axis label --
            if ikey==0:
                ax.set(ylabel='[g/m3]')

            #-- text annotation --
            for irect,rect in enumerate(rects):
                height = rect.get_height()
                #if irect ==0:                
                #    xy = (rect.get_x()*0.93, height + ymax*0.1)
                #else:
                #    xy = (rect.get_x()*0.93, height + ymax*0.3)
                xy = (rect.get_x()*0.93, height + ymax*0.1)
                ax.annotate('%.2f'%(height), xy=xy)

            ##-- text annotation(percentatge) --
            #for irect,rect in enumerate(rects):
            #    height = rect.get_height()
            #    if irect==0:
            #        height0 = height
            #        continue
            #    else:
            #        ratio = height/height0
            #    xy = (rect.get_x()*0.91, height + ymax*0.1)
            #    ax.annotate('(%.2f)'%(ratio), xy=xy)



        varname={'cond':'Mean condensed water content'}['cond']
        plt.suptitle('%s'%(varname), x=0.55, y=0.98 )
        plt.tight_layout(rect=[0, 0, 1, 0.9])
        plt.show()

        figpath = '/home/utsumi/temp/ret/bar.%s.all-type.png'%(var)
        plt.savefig(figpath)
        print figpath







#********************************************
#** Figure convective event fraction ****
#********************************************
for rettype in lrettype:
    expr = dexpr[rettype]
    srcdir = tankbaseDir + '/utsumi/PMM/validprof/metrix-orbit/%s.%s'%(rettype,expr)
    srcpath = srcdir + '/convCountFrac.%s.pickle'%(rettype)
    with open(srcpath,'rb') as f:
        dcfrac = pickle.load(f) 
    if rettype =='epc':
        dcfracepc = dcfrac
    elif rettype in ['gprof','gprof-shift']:
        dcfracgpr = dcfrac

#-----------------------
nrow, ncol = 1,4
fig, axs = plt.subplots(nrows=nrow, ncols=ncol, figsize=(8,1.6))
#lkey = [(radcondition,region)
#            for radcondition in ['stra','conv']
#            for region in ['tro','mid']
#            ]

lkey = [(region,ptype)
                    for region in ['mid','tro']  # 'all','mid','tro'
                    for ptype  in ['conv','stra']
]




#for iax,(radcondition,region) in enumerate(lkey):
for iax,(region,radcondition) in enumerate(lkey):
    ax = axs[iax]

    #avecmb = dcfracepc['rad',radcondition,region,'ave']
    aveepc = dcfracepc['pmw',radcondition,region,'ave']
    avegpr = dcfracgpr['pmw',radcondition,region,'ave']

    #p05cmb = dcfracepc['rad',radcondition,region,'p05']
    p05epc = dcfracepc['pmw',radcondition,region,'p05']
    p05gpr = dcfracgpr['pmw',radcondition,region,'p05']
    #p25cmb = dcfracepc['rad',radcondition,region,'p25']
    p25epc = dcfracepc['pmw',radcondition,region,'p25']
    p25gpr = dcfracgpr['pmw',radcondition,region,'p25']
    #p75cmb = dcfracepc['rad',radcondition,region,'p75']
    p75epc = dcfracepc['pmw',radcondition,region,'p75']
    p75gpr = dcfracgpr['pmw',radcondition,region,'p75']
    #p95cmb = dcfracepc['rad',radcondition,region,'p95']
    p95epc = dcfracepc['pmw',radcondition,region,'p95']
    p95gpr = dcfracgpr['pmw',radcondition,region,'p95']

    if radcondition=='stra':
        avecmb=0
    elif radcondition=='conv':
        avecmb=1
    else:
        print 'check radcondition',radcondition; sys.exit()

    #aerrmin= -np.array([p25epc,p25gpr]) + np.array([aveepc,avegpr])
    #aerrmax=  np.array([p75epc,p75gpr]) - np.array([aveepc,avegpr])
    aerrmin= -np.array([0,p25epc,p25gpr]) + np.array([0,aveepc,avegpr])
    aerrmax=  np.array([0,p75epc,p75gpr]) - np.array([0,aveepc,avegpr])



    labels = ['CMB','EPC','GPR']
    #ax.bar(labels, [aveepc,avegpr], yerr=[aerrmin,aerrmax], color='gray')
    rects  = ax.bar(labels, [avecmb,aveepc,avegpr], yerr=[aerrmin,aerrmax], color='gray')
    ymax = {'stra':0.55,'conv':1.2}[radcondition]
    ax.set_ylim([0,ymax])

    if iax==0:
        ax.set(ylabel='[frac. count]')

    #-- title --
    ptypename = {'stra':'Strat','conv':'Conv'}[radcondition]
    regionname= {'mid':'Mid-high','tro':'Tropics'}[region]
    #ax.set_title('%s (%s)'%(ptypename, regionname))
    ax.set_title('%s (%s)'%(regionname, ptypename))

    #-- text annotation --
    for rect in rects:
        height = rect.get_height()
        xy = (rect.get_x()*0.93, height + ymax*0.05)
        ax.annotate('%.2f'%(height), xy=xy)


plt.suptitle("Convective pixel fraction (Conditional on CMB's precip type)")
plt.tight_layout(rect=[0, 0, 1, 0.9])
plt.show()
figpath = '/home/utsumi/temp/ret/bar.convCountFrac.png'
plt.savefig(figpath)
print figpath


#********************************************
#-- Figure: Unconditional convective event fraction ---
#********************************************
nrow, ncol = 1,2
fig, axs = plt.subplots(nrows=nrow, ncols=ncol, figsize=(4.1,1.6))
radcondition = 'all'
for iregion,region in enumerate(['mid','tro']):
    ax = axs[iregion]
    avecmb = dcfracepc['rad',radcondition,region,'ave']
    aveepc = dcfracepc['pmw',radcondition,region,'ave']
    avegpr = dcfracgpr['pmw',radcondition,region,'ave']

    p25cmb = dcfracepc['rad',radcondition,region,'p25']
    p25epc = dcfracepc['pmw',radcondition,region,'p25']
    p25gpr = dcfracgpr['pmw',radcondition,region,'p25']
    p75cmb = dcfracepc['rad',radcondition,region,'p75']
    p75epc = dcfracepc['pmw',radcondition,region,'p75']
    p75gpr = dcfracgpr['pmw',radcondition,region,'p75']


    aerrmin= -np.array([0,p25epc,p25gpr]) + np.array([0,aveepc,avegpr])
    aerrmax=  np.array([0,p75epc,p75gpr]) - np.array([0,aveepc,avegpr])

    labels = ['CMB','EPC','GPR']
    rects  = ax.bar(labels, [avecmb,aveepc,avegpr], yerr=[aerrmin,aerrmax], color='gray')
    ymax = 1.2
    ax.set_ylim([0,ymax])

    if iregion==0:
        ax.set(ylabel='[frac. count]')

    #-- title --
    regionname= {'mid':'Mid-high','tro':'Tropics'}[region]
    ax.set_title('%s'%(regionname))

    #-- text annotation --
    for rect in rects:
        height = rect.get_height()
        xy = (rect.get_x()*0.93, height + ymax*0.05)
        ax.annotate('%.2f'%(height), xy=xy)


plt.suptitle("Convective pixel fraction")
plt.tight_layout(rect=[0, 0, 1, 0.9])
plt.show()
figpath = '/home/utsumi/temp/ret/bar.convCountFrac-unconditional.png'
plt.savefig(figpath)
print figpath


#********************************************
#-- Figure: Unconditional convective event fraction (Only CMB)
#********************************************
fig = plt.figure(figsize=(2.2,1.6))
ax  = fig.add_axes([0.2,0.25,0.6,0.45])
radcondition = 'all'
avetro = dcfracepc['rad',radcondition,'tro','ave']
avemid = dcfracepc['rad',radcondition,'mid','ave']

p25tro = dcfracepc['rad',radcondition,'tro','p25']
p25mid = dcfracepc['rad',radcondition,'mid','p25']
p75tro = dcfracepc['rad',radcondition,'tro','p75']
p75mid = dcfracepc['rad',radcondition,'mid','p75']

aerrmin= -np.array([p25tro,p25mid]) + np.array([avetro,avemid])
aerrmax=  np.array([p75tro,p75mid]) - np.array([avetro,avemid])

labels = ['Trop','Mid-high']
rects  = ax.bar(labels, [avetro,avemid], yerr=[aerrmin,aerrmax], color='gray')
ymax = 1.2
ax.set_ylim([0,ymax])

if iregion==0:
    ax.set(ylabel='[frac. count]')

#-- title --
#ax.set_title('CMB')

#-- text annotation --
for rect in rects:
    height = rect.get_height()
    xy = (rect.get_x()+0.2, height + ymax*0.05)
    print 'xy=',xy
    ax.annotate('%.2f'%(height), xy=xy)


plt.suptitle("Convective pixel\n fraction (CMB)")
plt.tight_layout(rect=[0, 0, 1, 0.9])
plt.show()
figpath = '/home/utsumi/temp/ret/bar.convCountFrac-CMB.png'
plt.savefig(figpath)
print figpath






# %%
