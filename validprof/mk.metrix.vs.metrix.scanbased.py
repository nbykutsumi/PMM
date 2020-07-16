# %%
import matplotlib
matplotlib.use('Agg')
#%matplotlib inline
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.ticker import ScalarFormatter
import numpy as np
import os, sys
from numpy import ma
import myfunc.util as util
from numpy import *
from collections import deque
import calendar, glob
from datetime import datetime, timedelta
import random
import socket
import scipy.stats as stats
import pickle
import string

useorblist = True
calcflag = True
#calcflag = False
#lrettype = ['epc']
lrettype = ['epc','gprof-shift']
nsample = 1000
#nsample = 10

iDTime = datetime(2014,6,1)
eDTime = datetime(2015,5,31)
#eDTime = datetime(2014,6,30)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
lskipdates = [[2014,8,29],[2014,9,16],[2014,10,1],[2014,10,2],[2014,11,5],[2014,12,8],[2014,12,9],[2014,12,10]]

DB_MAXREC = 10000
DB_MINREC = 1000
dexpr = {'epc':'glb.relsurf01.minrec1000.maxrec10000','gprof-shift':'v01'}

thpr = 0.5 # mm/h
thwat = 0.033  # g/m3 for storm top
nz    = 20
miss  = -9999.
#xymetric = ['dvfracConv','ndprec']
#xymetric = ['dwatNorm','ndprec']
#xymetric = ['cc','ndprec']
#xymetric = ['dstop-prof','ndprec']
#xymetric = ['dstop','ndprec']

xymetric = ['ndprec','cc']
#xymetric = ['ndprec','dwatNorm']
#xymetric = ['ndprec','dstop-prof']
#xymetric = ['ndprec','dstop']
#xymetric = ['ndprec','dvfracConv']
#xymetric = ['ndprec','rmse']


#xymetric = ['dtqv','dwatNorm']
#xymetric = ['dwatNorm','dprec']
#xymetric = ['dconvfrac','dwatNorm']
#xymetric = ['dstop','dprec']
#xymetric = ['dstop','dwatNorm']
#xymetric = ['dvfracConv','dwatNorm']
#xymetric = ['dvfracConv','dprec']
#xymetric = ['dprec','dstop']
#xymetric = ['dprec','dvfracConv']
#xymetric = ['dprec','dwatNorm']
#xymetric = ['dprec','rmse']
#xymetric = ['dprec','cc']
#xymetric = ['cc','dprec']
#xymetric = ['rmse','dprec']

#lsurf = ['all','sea','veg']
lsurf = ['all']

myhost = socket.gethostname()
if myhost =='shui':
    tankbaseDir = '/tank'
    workbaseDir = '/work'
    epcbaseDir  = '/tank/utsumi/PMM/retepc'

elif myhost == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    workbaseDir = '/home/utsumi/mnt/lab_work'
    epcbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc'

else:
    print 'check myhost'
    sys.exit()


#lprec = arange(0.25, 20+0.01, 0.5)
#lbiaslev = [0,1,2,3]
dprecrange = {-1:[0.5, 999],0:[1,3], 1:[5,7], 2:[10,12]}
#dprecrange = {2:[10,12]}
aprecrange = dprecrange.values()
#lbias = np.arange(0.,20+0.1, 1.0)
dxbnd = {
         'dprec':np.arange(-10,10+0.01,0.5)
        ,'ndprec':np.arange(-3,3,0.1)
        ,'dwatNorm':np.arange(-1,5+0.01,0.2)
        ,'dtqv': np.arange(-20,20+0.1, 2.0)  # kg/m2
        ,'dconvfrac': np.arange(-1,1+0.01, 0.1)  # kg/m2
        ,'dstop':np.arange(-4-0.25, 4+0.25, 0.5) # km
        ,'dstop-prof':np.arange(-4-0.25, 4+0.25, 0.5) # km
        ,'dvfracConv':np.arange(-1-0.05, 1+0.05, 0.1) # km
        #,'cc':np.arange(-0.2-0.05, 1+0.05, 0.05) # km
        ,'cc':np.arange(-0.2, 1+0.01, 0.05) # km
        ,'rmse':np.arange(0, 1.0+0.01, 0.05) # km
        }
#lbias = np.arange(-10,10+0.1,1.0)


dminmax = {
         'dprec':[-10,10]
        ,'ndprec':[-1,2]
        ,'dwatNorm':[-1,5]
        ,'dtqv': [None,None]
        ,'dconvfrac': [None,None]
        ,'dstop':[None,None]
        ,'dstop-prof':[None,None]
        ,'dvfracConv': [-1,1]
        ,'cc': [0, 1]
        ,'rmse':[-1,1]
        }


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


#*****************************************
def calc_rmse(a2ref, a2dat):
    a2ref = ma.masked_less(a2ref,0)
    a2dat = ma.masked_less(a2dat,0)
    return np.sqrt( ((a2dat-a2ref)**2).mean(axis=1) )    

def calc_cc(x,y,axis):
    x = ma.masked_less(x,0)
    y = ma.masked_less(y,0)
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

def calc_dwatNorm(a2ref, a2dat):
    a1ref = a2ref.sum(axis=1)
    a1dat = a2dat.sum(axis=1)
    return ma.masked_where(a1ref==0, a1dat - a1ref)/a1ref

def calc_ptypefrac(avar, ptype='conv'):
    avar = ma.masked_less(avar,0).astype('int16')
    strat = (avar %10).astype('int16')
    conv  = (avar%100-strat).astype('int16')/10
    other = (avar/100).astype('int16')
    ntot  = (strat + conv + other).astype('float32')
    
    if ptype=='conv':
        return conv / ntot
    elif ptype=='strat':
        return strat / ntot
    elif ptype=='other':
        return other / ntot
    else:
        print 'check',ptype
        sys.exit()

def calc_dptypefrac(a1ref, a1dat, ptype='conv'):
    a1reffrac = calc_ptypefrac(a1ref,ptype=ptype)
    a1datfrac = calc_ptypefrac(a1dat,ptype=ptype)
    return ma.masked_less(a1datfrac,0) - ma.masked_less(a1reffrac,0)



def get_mycm(cmap=plt.cm.jet):
    cm = cmap
    cm_list = cm(np.arange(cm.N))[:-20]
    return ListedColormap(cm_list)

#*****************************************
xmetric, ymetric = xymetric
for rettype in lrettype:
    if calcflag is not True: continue

    #**** File list ********************
    expr = dexpr[rettype]
    if rettype=='epc': 
        pairbaseDir = tankbaseDir + '/utsumi/PMM/validprof/pair/%s.%s'%(rettype, expr)
    elif rettype=='gprof-shift':
        pairbaseDir = tankbaseDir + '/utsumi/PMM/validprof/pair/%s'%(rettype)
    else:
        print 'check rettype',rettype
        sys.exit()



    ##-- Initialize ------
    o1precrad = deque([])
    o1precpmw = deque([])
    o1xmetric = deque([])
    o1ymetric = deque([])
    o1lat     = deque([])
    o1lon     = deque([])
    dout      = {}

    #--------------------
    for (oid,Year,Mon,Day) in loid:
        pairDir = pairbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)

        a1lat = np.load(pairDir + '/Latitude.%06d.npy'%(oid))
        a1lon = np.load(pairDir + '/Longitude.%06d.npy'%(oid))
    
        a1precpmw = np.load(pairDir+ '/precpmw.%06d.npy'%(oid))
        a1precrad = np.load(pairDir+ '/precrad.%06d.npy'%(oid))
        a1surftype= np.load(pairDir + '/surfaceTypeIndex.%06d.npy'%(oid))
        a1elev    = np.load(pairDir + '/surfaceElevationrad.%06d.npy'%(oid))

        a2pmw = np.load(pairDir+ '/profpmw.%06d.npy'%(oid))[:,:nz]  # nz(500m) layers, bottom to top
        a2rad = np.load(pairDir+ '/profrad.%06d.npy'%(oid))[:,:nz]  # nz(500m) layers, bottom to top

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
            a2pmw[a1idxtmp, :cltbin] = miss
            a2rad[a1idxtmp, :cltbin] = miss

        #*****************************************
        # Screen input data
        #*****************************************
        if rettype == 'gprof':
            a1qflag = np.load(pairDir + '/qualityFlag.%06d.npy'%(oid))
            a1flagQ = ma.masked_equal(a1qflag,0).mask   # Good quality (flag=0)
        else:
            a1flagQ = np.array([True])

        #--- Screen Sea-ice(2), Snow(8-11), Sea-ice edge(14) ---
        a1flagS = ma.masked_equal(np.load(pairDir + '/surfaceTypeIndex.%06d.npy'%(oid)) ,2)
        a1flagS = ma.masked_inside(a1flagS,8,11)
        a1flagS = ma.masked_equal(a1flagS,14)
        a1flagS = ~(a1flagS.mask)
    
        #--- Screen missing surface precipitation ---
        a1flagP1 = ma.masked_invalid(ma.masked_greater_equal(a1precpmw,thpr)).mask
        a1flagP2 = ma.masked_invalid(ma.masked_greater_equal(a1precrad,thpr)).mask
        a1flagP  = a1flagP1 * a1flagP2 

        ##--- Screen high mountains ------------------
        a1flagE = np.array([True])
        #a1flagE  = ma.masked_less(a1elev, 1000).mask

        #--- Check quality flag --
        a1flag  = a1flagS * a1flagP *a1flagQ * a1flagE
        
        ##---------------------
        if a1flag.sum()==0:
            continue

        a1precpmw = a1precpmw[a1flag]
        a1precrad = a1precrad[a1flag]
        a1surftype= a1surftype[a1flag]
        a1elev    = a1elev[a1flag]
        a1lat     = a1lat[a1flag]
        a1lon     = a1lon[a1flag]
        a2pmw     = a2pmw[a1flag]
        a2rad     = a2rad[a1flag]

        #-- mask missing values ----
        a2mask1 = ma.masked_less(a2pmw,0).mask
        a2mask2 = ma.masked_invalid(a2pmw).mask
        a2mask3 = ma.masked_less(a2rad,0).mask
        a2mask4 = ma.masked_invalid(a2rad).mask
        a2mask  = a2mask1 + a2mask2 + a2mask3 + a2mask4

        a2pmw = ma.masked_where(a2mask, a2pmw)
        a2rad = ma.masked_where(a2mask, a2rad)


        a1flagM   = np.array([True])
        for imetric,metric in enumerate(xymetric):

            if metric=='dprec':
                a1metric = ma.masked_less(a1precpmw,0) - ma.masked_less(a1precrad,0)
                a1flagM = ~ma.masked_invalid(a1metric).mask

            elif metric=='ndprec':
                a1metric = ma.masked_less(a1precpmw,0) - ma.masked_less(a1precrad,0)
                a1metric = ma.masked_where(a1precrad==0, a1metric)/a1precrad
                a1flagM = ~ma.masked_invalid(a1metric).mask

            elif metric=='rmse':
                a1metric = calc_rmse(a2rad, a2pmw)

            elif metric=='cc':
                a1metric = calc_cc(a2rad[:,:15], a2pmw[:,:15], axis=1)  # to 7.5km
                a1flagM  = a1flagM * ~ma.masked_invalid(a1metric).mask

            elif metric=='dwatNorm':
                a1metric = calc_dwatNorm(a2ref=a2rad, a2dat=a2pmw)
                a1flagM  = a1flagM * ~ma.masked_invalid(a1metric).mask

            elif metric=='dconvfrac':
                a1pmw=np.load(pairDir + '/top-typePrecippmw.%06d.npy'%(oid))
                a1rad=np.load(pairDir + '/typePreciprad.%06d.npy'%(oid))
 
                a1metric = calc_dptypefrac(a1rad, a1pmw, ptype='conv')


            elif metric=='dvfracConv':
                if rettype=='epc':
                    a1pmw=np.load(pairDir + '/top-vfracConvpmw.%06d.npy'%(oid))
                elif rettype=='gprof':
                    a1pmw=np.load(pairDir + '/vfracConvpmw.%06d.npy'%(oid))

                a1rad=np.load(pairDir + '/vfracConvrad.%06d.npy'%(oid))
                a1metric = ma.masked_less(a1pmw,0) - ma.masked_less(a1rad,0)

            elif metric=='dstop-prof':
                #a1metric = (ma.masked_less_equal(a1pmw,0) - ma.masked_less_equal(a1rad,0)) * 0.001  # [km]
                nltmp = a2rad.shape[0]
                a2iz  = np.array(range(nz)*nltmp).reshape(-1,nz)
                a1stoprad = ma.masked_where(a2rad < thwat, a2iz).argmax(axis=1)*0.5 - a1elev*0.001 + 0.5
                a1stoppmw = ma.masked_where(a2pmw < thwat, a2iz).argmax(axis=1)*0.5 - a1elev*0.001 + 0.5
                a1metric = a1stoppmw - a1stoprad
                a1flagM  = a1flagM * ~ma.masked_invalid(a1metric).mask


            elif metric=='dstop':
                a1pmw = np.load(pairDir + '/top-heightStormToppmw.%06d.npy'%(oid))
                a1rad = np.load(pairDir + '/heightStormToprad.%06d.npy'%(oid))

                a1metric = (ma.masked_less_equal(a1pmw,0) - ma.masked_less_equal(a1rad,0)) * 0.001  # [km]



            elif metric=='dtqv':
                a1pmw =np.load(pairDir + '/top-tqvpmw.%06d.npy'%(oid))
                a1rad =np.load(pairDir + '/tqv.%06d.npy'%(oid))
                a1metric = ma.masked_less(a1pmw,0) - ma.masked_less(a1rad,0)
        

            else:
                print 'check metric',metric
                sys.exit()


            a1flagM = a1flagM * ~ma.masked_invalid(a1metric).mask
            #----------------
            if imetric==0:
                a1xmetric = a1metric
            elif imetric==1:
                a1ymetric = a1metric

        #-------------------------------
        o1precpmw.extend(list(a1precpmw))
        o1precrad.extend(list(a1precrad))
        o1xmetric.extend(list(a1xmetric))
        o1ymetric.extend(list(a1ymetric))
        o1lat    .extend(list(a1lat))
        o1lon    .extend(list(a1lon))

    #--- Screen by metrics -----
    a1maskM1 = ma.masked_invalid(o1xmetric).mask
    a1maskM2 = ma.masked_invalid(o1ymetric).mask
    a1maskM  = a1maskM1 + a1maskM2
    #a1flagM  = ~( a1flagM1 * a1flagM2 )
    a1flagM  = ~a1maskM

    a1precpmw = np.array(o1precpmw)[a1flagM]
    a1precrad = np.array(o1precrad)[a1flagM]
    a1xmetric = np.array(o1xmetric)[a1flagM]
    a1ymetric = np.array(o1ymetric)[a1flagM]
    a1lat     = np.array(o1lat    )[a1flagM]
    a1lon     = np.array(o1lon    )[a1flagM]

    #--- Classify by surface --- 
    d1flagS = {} 
    d1flagS['sea'] = ma.masked_equal(a1surftype,1).mask
    d1flagS['veg'] = ma.masked_inside(a1surftype,3,7).mask
    d1flagS['snow']= ma.masked_inside(a1surftype,8,11).mask
    d1flagS['coast']=ma.masked_equal(a1surftype,13).mask
    d1flagS['land'] =d1flagS['veg'] + d1flagS['snow']
    d1flagS['all'] = np.array([True])

    for surf in lsurf:
        a1flagS = d1flagS[surf]

        #--- Classify by prec ------
        for preclev in dprecrange.keys():
            iprec,eprec = dprecrange[preclev]
            a1flagP = ma.masked_inside(a1precrad, iprec, eprec).mask

            a1flag  = a1flagS * a1flagP       
            if a1flag.sum()==0: continue
    
            a1xmetricTmp = a1xmetric[a1flag]
            a1ymetricTmp = a1ymetric[a1flag]
    
            lxbnd = dxbnd[xmetric]
            lybnd = dxbnd[ymetric]

            #-- 2-D Histogram --
            H,xedges,yedges = np.histogram2d(a1xmetricTmp, a1ymetricTmp, bins = [lxbnd,lybnd])
            H = H.T
            dout['hist2d',surf,preclev] = H

            #-- 1-D Histogram --
            bins = np.sort(lxbnd)
            a1xhist,a1bnd,_ = stats.binned_statistic(x=a1xmetricTmp, values=None, statistic='count', bins=bins)

            bins = np.sort(lybnd)
            a1yhist,a1bnd,_ = stats.binned_statistic(x=a1ymetricTmp, values=None, statistic='count', bins=bins)
            dout['xhist',surf,preclev]= a1xhist
            dout['yhist',surf,preclev]= a1yhist

            #--- Mean of y binned by x ---
            bins = np.sort(lxbnd)
            a1yave,a1bnd,_ = stats.binned_statistic(x=a1xmetricTmp, values=a1ymetricTmp, statistic='mean', bins=bins)

            dout['yave',surf,preclev] = a1yave

            #--- Mean of x binned by y ---
            bins = np.sort(lybnd)
            a1xave,a1bnd,_ = stats.binned_statistic(x=a1ymetricTmp, values=a1xmetricTmp, statistic='mean', bins=bins)
            dout['xave',surf,preclev] = a1xave

            #--- Median of y binned by x ---
            bins = np.sort(lxbnd)
            a1ymed,a1bnd,_ = stats.binned_statistic(x=a1xmetricTmp, values=a1ymetricTmp, statistic='median', bins=bins)

            dout['ymed',surf,preclev] = a1ymed

            #--- Median of x binned by y ---
            bins = np.sort(lybnd)
            a1xmed,a1bnd,_ = stats.binned_statistic(x=a1ymetricTmp, values=a1xmetricTmp, statistic='median', bins=bins)
            dout['xmed',surf,preclev] = a1xmed



    #--- Other parameters ----
    dout['dprecrange'] = dprecrange
    dout['xbnd'] = lxbnd
    dout['ybnd'] = lybnd

    #--- Save -----
    outDir  = tankbaseDir + '/utsumi/PMM/validprof/metrix.vs.metrix/pickle'
    print outDir
    util.mk_dir(outDir)

    histPath= outDir + '/hist2d.%s.%s.%s.pickle'%(xmetric,ymetric,rettype)

    with open(histPath, 'wb') as f:
        pickle.dump(dout, f)
    print histPath



##************************************
# Figure
#************************************
for irettype,rettype in enumerate(lrettype):
    #*************************************
    # Read
    #*************************************
    srcDir  = tankbaseDir + '/utsumi/PMM/validprof/metrix.vs.metrix/pickle'
    histPath= srcDir + '/hist2d.%s.%s.%s.pickle'%(xmetric,ymetric,rettype)
    with open(histPath, 'rb') as f:
        dvar = pickle.load(f)

    lxbnd = dvar['xbnd']
    lybnd = dvar['ybnd']
    dprecrange = dvar['dprecrange']
    ##************************************
    ## Plot (xmetric vs ymetric)
    ##************************************
    dlabel = {'dprec':'Surface precipitation error (mm/h)'
             ,'ndprec':'Norm. surf. precip. error'
             ,'rmse':'Profile RMSE (g/m3)'
             ,'cc'  :'Profile Correlation Coef.'
             ,'dwatNorm' :'Mean condensed water\n error (Normalized)'
             ,'dconvfrac':'Convective fraction error'
             ,'dstop-prof':'Storm top height error\n(profile-based) (km)'
             ,'dstop':'Storm top height error\n(top-weighted) (km)'
             ,'dtqv' :'Total water vapor error (kg/m2)'
             ,'dvfracConv' :'convective fraction error'
             }

    #************************************
    # Plot 2D histogram + lines
    #************************************
    for surf in lsurf:
        H  = dvar['hist2d',surf,-1]
        X,Y= np.meshgrid(lxbnd, lybnd)

        fig = plt.figure(figsize=(6,6))
        ax = fig.add_axes([0.25,0.20,0.45,0.45])

        if irettype==0: 
            vmax= H.max()
        im = ax.pcolormesh(X,Y,H, norm=matplotlib.colors.LogNorm(), cmap='jet', vmin=10, vmax=vmax)

   
        if ymetric =='cc':
            #--- Lines (y-mean) ---
            preclev = -1
            linestyle= '-'
            linewidth= 2
            mycolor = 'k'

            a1x = (lxbnd[:-1] + lxbnd[1:])*0.5
            #a1y = dvar['ymed',surf,preclev]
            a1y = dvar['yave',surf,preclev]

            ax.plot(a1x, a1y, linestyle=linestyle, linewidth=linewidth, color=mycolor)

            ##--- Lines (x-mean) ---
            #preclev = -1
            #linestyle= '-'
            #linewidth= 2
            #mycolor = 'k'
            #a1x = dvar['xave',surf,preclev]
            #a1y = (lybnd[:-1] + lybnd[1:])*0.5

            #ax.plot(a1x, a1y, linestyle=linestyle, linewidth=linewidth, color=mycolor)




        ##--- Lines (y-mean) ---
        #for preclev in dprecrange.keys():
        #    #if preclev ==-1: continue

        #    iprec,eprec = dprecrange[preclev]
        #    slabel = '%d-%dmm/h'%(iprec,eprec)
        #    linestyle= ['-', '--', '-'][preclev]
        #    linewidth= [1, 2, 2][preclev]
        #    #mycolor = ['darkblue','orange','crimson'][preclev]
        #    mycolor = ['k','k','k'][preclev]
    
        #    a1x = (lxbnd[:-1] + lxbnd[1:])*0.5
        #    a1y = dvar['yave',surf,preclev]

        #    ax.plot(a1x, a1y, linestyle=linestyle, linewidth=linewidth, color=mycolor, label=slabel)

        ##--- Lines (x-mean) ---
        #for preclev in dprecrange.keys():
        #    if preclev ==-1: continue

        #    iprec,eprec = dprecrange[preclev]
        #    slabel = '%d-%dmm/h'%(iprec,eprec)
        #    linestyle= ['-', '--', '-'][preclev]
        #    linewidth= [1, 2, 2][preclev]
        #    #mycolor = ['darkblue','orange','crimson'][preclev]
        #    mycolor = ['k','k','k'][preclev]
    
        #    a1x = dvar['xave',surf,preclev]
        #    a1y = (lybnd[:-1] + lybnd[1:])*0.5

        #    ax.plot(a1x, a1y, linestyle=linestyle, linewidth=linewidth, color=mycolor, label=slabel)
        #------------

        #------------

        #iprec,eprec = dprecrange[preclev]
 
        
        xlabel = dlabel[xmetric]
        ylabel = dlabel[ymetric]
        
        ymin,ymax = dminmax[ymetric]
        xmin,xmax = dminmax[xmetric]
 
        ax.set_ylim([ymin,ymax])
        ax.set_xlim([xmin,xmax])
 
        ax.set_ylabel(ylabel, fontsize=20)
        ax.set_xlabel(xlabel, fontsize=20)

        ax.xaxis.set_tick_params(labelsize=16)
        ax.yaxis.set_tick_params(labelsize=16)
        ax.axvline(x=0, linestyle=':', color='gray', linewidth=4)
        ax.axhline(y=0, linestyle=':', color='gray', linewidth=4)


        ##-------------------------------
        ##-- PDF : top
        ##-------------------------------  
        #a1x = (lxbnd[:-1] + lxbnd[1:])*0.5

        #axpdfx = fig.add_axes([0.25,0.67, 0.45, 0.05])
        #xpdf = ma.masked_less(H,0).sum(axis=0)
        #xpdf = xpdf / float(xpdf.max())

        #axpdfx.plot(a1x, xpdf, '-', linewidth=2 , color='k')
        #axpdfx.set_xlim([xmin,xmax])
        #axpdfx.tick_params(labelbottom=False, labelleft=False)
        #axpdfx.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        #axpdfx.ticklabel_format(style="sci",  axis="y",scilimits=(0,0))
        #-------------------------------
        #-- PDF : Right
        #-------------------------------  
        a1y = (lybnd[:-1] + lybnd[1:])*0.5
        
        axpdfy = fig.add_axes([0.72,0.20, 0.05, 0.45])
        ypdf = ma.masked_less(H,0).sum(axis=1)
        ypdf = ypdf / float(ypdf.max())

        axpdfy.plot(ypdf, a1y, '-', linewidth=2 , color='k')
        axpdfy.set_ylim([ymin,ymax])
        axpdfy.tick_params(labelbottom=False, labelleft=False)
        axpdfy.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        axpdfy.ticklabel_format(style="sci",  axis="x",scilimits=(0,0))

        #-------------------------------
        # Title
        #-------------------------------
        rettypename = {'epc':'EPC','gprof-shift':'GPROF'}[rettype]
        #stitle = '%s %s'%(string.upper(rettype), surf)
        stitle = '%s'%(rettypename)
        ax.set_title(stitle, fontsize=20, pad=10)



        #--- Separate legend -----
        figleg = plt.figure(figsize=(6,2)) 
        axleg  = figleg.add_axes([0.1,0.1,0.8,0.8])
        axleg.legend(*ax.get_legend_handles_labels(), loc='center',fontsize=20)
        axleg.xaxis.set_visible(False)
        axleg.yaxis.set_visible(False)
        #-------------------------

        if xmetric=='rmse':
            xmin, xmax= 0, 0.7
        elif xmetric=='dwatNorm':
            xmin, xmax= -1.2, 5.5
        else:
            xmin, xmax= None,None
        ax.set_xlim([xmin,xmax])

        #-- colorbar ---
        cax = fig.add_axes([0.79,0.20,0.02, 0.44])
        cbar=plt.colorbar(im, orientation='vertical', cax=cax)
        cbar.ax.tick_params(labelsize=16)


        figDir = '/home/utsumi/temp/ret'
        figPath = figDir + '/plot.scatter.%s.vs.%s.%s.%s.png'%(xmetric,ymetric, rettype,surf)
        fig.savefig(figPath)

        legPath = figDir + '/legend.plot.%s.vs.%s.%s.%s.png'%(xmetric,ymetric, rettype,surf)
        figleg.savefig(legPath)


        plt.clf()
        plt.close()
        print figPath


# %%
