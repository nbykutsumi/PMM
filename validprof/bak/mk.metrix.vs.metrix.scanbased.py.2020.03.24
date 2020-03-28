import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.ticker import ScalarFormatter
import numpy as np
import os, sys
from numpy import ma
import myfunc.util as util
from numpy import *
import calendar, glob
from datetime import datetime, timedelta
import random
import socket
import scipy.stats as stats
import pickle
import string

season = 'ALL'
#season = 'JJADJF'
#season = 6
#calcflag = True
calcflag = False
lrettype = ['epc','gprof']
#lrettype = ['epc']
#lrettype = ['gprof']
DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
#nsample = 18 # sampling orbits per day
nsample = 4 # sampling orbits per day

lskipdates = [[2014,6,4],[2014,6,6],[2014,10,22],[2014,10,23],[2014,10,24],[2014,12,9],[2014,12,10],[2014,11,25]]


#xymetric = ['dvfracConv','ndprec']
#xymetric = ['dwatNorm','ndprec']
#xymetric = ['cc','ndprec']
#xymetric = ['dstop-prof','ndprec']
#xymetric = ['dstop','ndprec']

#xymetric = ['ndprec','cc']
#xymetric = ['ndprec','dwatNorm']
#xymetric = ['ndprec','dstop-prof']
xymetric = ['ndprec','dstop']
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
nz    = 25


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
    a1ref = ma.masked_less(a2ref,0).filled(0).sum(axis=1)
    a1dat = ma.masked_less(a2dat,0).filled(0).sum(axis=1)
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


def ret_lmon(season):
    if   season=='ALL':
        lYM = util.ret_lYM([2014,6],[2015,5])

    elif season=='JJADJF':
        lYM = util.ret_lYM([2014,6],[2014,8]) + util.ret_lYM([2014,12],[2015,2])

    elif season=='JJA':
        lYM = util.ret_lYM([2014,6],[2014,8])
    elif season=='SON':
        lYM = util.ret_lYM([2014,9],[2014,11])
    elif season=='DJF':
        lYM = util.ret_lYM([2014,12],[2015,2])
    elif season=='MAM':
        lYM = util.ret_lYM([2015,3],[2015,5])
    elif type(season)==int:
        if season <6:
            lYM = [[2015,season]]
        else:
            lYM = [[2014,season]]
    elif season =='JJ':
        lYM = util.ret_lYM([2014,6],[2014,7])
    elif season =='DJ':
        lYM = util.ret_lYM([2014,12],[2015,1])

    else:
        print 'check season',season
        sys.exit()
    return lYM


def get_mycm(cmap=plt.cm.jet):
    cm = cmap
    cm_list = cm(np.arange(cm.N))[:-20]
    return ListedColormap(cm_list)

#*****************************************
lYM = ret_lmon(season)
xmetric, ymetric = xymetric
for rettype in lrettype:

    #-- Initialize ------
    da1sum = {}
    da1sum2= {}
    da1num = {}
    da2hist= {}
    lxbnd = dxbnd[xmetric]
    lybnd = dxbnd[ymetric]
    for surf in lsurf:
        for preclev in dprecrange.keys():
            da1sum [surf,preclev]= np.zeros(len(lxbnd)-1, float32)
            da1sum2[surf,preclev]= np.zeros(len(lxbnd)-1, float32)
            da1num [surf,preclev]= np.zeros(len(lxbnd)-1, int32)
            da2hist[surf,preclev]= np.zeros([len(lybnd)-1, len(lxbnd)-1], float32)
    #--------------------

    for Year,Mon in lYM:
        if calcflag ==False: continue
    
    
        eDay = calendar.monthrange(Year,Mon)[1]
        iDTime = datetime(Year,Mon,1)
        eDTime = datetime(Year,Mon,eDay)
        dDTime = timedelta(days=1)
        lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)
        #lDTime = lDTime[:10]  # test 
        for DTime in lDTime:
            Year,Mon,Day = DTime.timetuple()[:3]
            if [Year,Mon,Day] in lskipdates: continue
    
            if rettype =='epc':
                pairDir = tankbaseDir + '/utsumi/PMM/validprof/pair/epc.%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)
            elif rettype=='gprof':
                pairDir = tankbaseDir + '/utsumi/PMM/validprof/pair/gprof/%04d/%02d/%02d'%(Year,Mon,Day)
            else:
                print 'check rettype',rettype
                sys.exit()

            ssearch= pairDir + '/Latitude.*.npy'
            llatPath= sort(glob.glob(ssearch))
            #-- random sampling -----
            if len(llatPath)>=nsample:
                llatPath = random.sample(llatPath, nsample)
            else:
                pass
            #------------------------
            
            for latPath in llatPath:
                #print latPath
                oid = int(latPath.split('.')[-2])
                print rettype, xymetric , DTime
                #if oid !=1918: continue  # test

                a1lat = np.load(pairDir + '/Latitude.%06d.npy'%(oid))
                a1lon = np.load(pairDir + '/Longitude.%06d.npy'%(oid))
    
                a1precrad = np.load(pairDir+ '/precrad.%06d.npy'%(oid))
                a1precpmw = np.load(pairDir+ '/precpmw.%06d.npy'%(oid))
                a1surftype= np.load(pairDir + '/surfaceTypeIndex.%06d.npy'%(oid))
                a1elev    = np.load(pairDir + '/surfaceElevationrad.%06d.npy'%(oid))

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
                        a2rad = np.load(pairDir+ '/profrad.%06d.npy'%(oid))[:,4:nz]  # bottom to top
                        a2pmw = np.load(pairDir+ '/profpmw.%06d.npy'%(oid))[:,4:nz]  # bottom to top
                        a1metric = calc_rmse(a2rad, a2pmw)

                    elif metric=='cc':
                        a2rad= np.load(pairDir+ '/profrad.%06d.npy'%(oid))[:,4:nz]  # bottom to top
                        a2pmw= np.load(pairDir+ '/profpmw.%06d.npy'%(oid))[:,4:nz]  # bottom to top
                        a1metric = calc_cc(a2rad, a2pmw, axis=1)


                    elif metric=='dwatNorm':
                        a2pmw = np.load(pairDir+ '/profrad.%06d.npy'%(oid))[:,4:nz]  # bottom to top
                        a2rad = np.load(pairDir+ '/profpmw.%06d.npy'%(oid))[:,4:nz]  # bottom to top
                        a1metric = calc_dwatNorm(a2pmw, a2rad)
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
                        a1pmw = np.load(pairDir + '/stop-profpmw.%06d.npy'%(oid))
                        a1rad = np.load(pairDir + '/stop-profrad.%06d.npy'%(oid))

                        a1metric = (ma.masked_less_equal(a1pmw,0) - ma.masked_less_equal(a1rad,0)) * 0.001  # [km]

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
                #*****************************************
                if rettype == 'gprof':
                    a1qflag   = np.load(pairDir + '/qualityFlag.%06d.npy'%(oid))

                    a1flagS = ma.masked_equal(np.load(pairDir + '/surfaceTypeIndex.%06d.npy'%(oid)) ,2)
                    a1flagS = ma.masked_inside(a1flagS,8,11)
                    a1flagS = ma.masked_equal(a1flagS,14)
                    a1flagS = ~(a1flagS.mask)

    
                #--- Screen missing surface precipitation ---
                a1flagP1 = ma.masked_greater_equal(a1precpmw,0).mask
                a1flagP2 = ma.masked_greater_equal(a1precrad,0).mask
                a1flagP  = a1flagP1 * a1flagP2 

                #--- Screen high mountains ------------------
                a1flagE  = ma.masked_less(a1elev, 1000).mask

                #--- Check quality flag --
                if rettype == 'gprof':
                    a1flagQ = ma.masked_equal(a1qflag,0).mask   # Good quality (flag=0)


                    #print a1flagM.shape, a1flagS.shape, a1flagP.shape, a1flagQ.shape, a1flagE.shape
                    a1flag  = a1flagM * a1flagS * a1flagP *a1flagQ * a1flagE
                else:
                    a1flag  = a1flagM * a1flagP * a1flagE
            
                ##---------------------
                if a1flag.sum()==0:
                    continue

                a1precpmw = a1precpmw[a1flag]
                a1precrad = a1precrad[a1flag]
                a1xmetric = a1xmetric[a1flag]
                a1ymetric = a1ymetric[a1flag]

                a1surftype= a1surftype[a1flag]
  
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
      
                        #-- bin by x-bnd --
                        bins = np.sort(lxbnd)
                        a1sum,a1bnd,_ = stats.binned_statistic(a1xmetricTmp, a1ymetricTmp, statistic='sum', bins=bins)
                        a1num,a1bnd,_ = stats.binned_statistic(a1xmetricTmp, a1ymetricTmp, statistic='count', bins=bins)
                        a1sum2    = a1sum**2
       
                        a1sum  = ma.masked_invalid(a1sum).filled(0.0)
                        a1sum2 = ma.masked_invalid(a1sum2).filled(0.0)
                        a1num  = ma.masked_invalid(a1num).filled(0.0)
                       
                        da1sum [surf,preclev] = da1sum [surf,preclev] + a1sum
                        da1sum2[surf,preclev] = da1sum2[surf,preclev] + a1sum2
                        da1num [surf,preclev] = da1num [surf,preclev] + a1num

                        H,xedges,yedges = np.histogram2d(a1xmetricTmp, a1ymetricTmp, bins = [lxbnd,lybnd])
                        H = H.T
                        da2hist[surf,preclev] = da2hist[surf,preclev] + H

      
    #--- Save -----
    xmetric,ymetrci = xymetric
    outDir  = tankbaseDir + '/utsumi/PMM/validprof/metrix.vs.metrix/pickle'
    print outDir
    util.mk_dir(outDir)
    sumPath = outDir + '/sum.%s.%s.%s.%s.binned.bfile'%(xmetric,ymetric,rettype,season)
    numPath = outDir + '/num.%s.%s.%s.%s.binned.bfile'%(xmetric,ymetric,rettype,season) 
    sum2Path = outDir + '/sum2.%s.%s.%s.%s.binned.bfile'%(xmetric,ymetric,rettype,season)
    precrangePath = outDir + '/precrange.%s.%s.%s.%s.bfile'%(xmetric,ymetric,rettype,season)
    xbndPath= outDir + '/xbnd.%s.%s.%s.%s.npy'%(xmetric,ymetric,rettype,season) 
    ybndPath= outDir + '/ybnd.%s.%s.%s.%s.npy'%(xmetric,ymetric,rettype,season) 

    histPath= outDir + '/hist.%s.%s.%s.%s.binned.bfile'%(xmetric,ymetric,rettype,season)



    if calcflag == True:
        with open(sumPath, 'wb') as f:
            pickle.dump(da1sum, f)
        with open(numPath, 'wb') as f:
            pickle.dump(da1num, f)
        with open(sum2Path, 'wb') as f:
            pickle.dump(da1sum2, f) 
        with open(precrangePath, 'wb') as f:
            pickle.dump(dprecrange, f)
        with open(histPath, 'wb') as f:
            pickle.dump(da2hist, f)


        lxbnd.tofile(xbndPath)
        lybnd.tofile(ybndPath)
        print sumPath

    ##************************************
    # Read
    #************************************
    with open(sumPath, 'rb') as f:
        da1sum = pickle.load(f)
    with open(numPath, 'rb') as f:
        da1num = pickle.load(f)
    with open(sum2Path,'rb') as f:
        da1sum2 = pickle.load(f) 
    with open(precrangePath, 'rb') as f:
        dprecrange = pickle.load(f)
    with open(histPath, 'rb') as f:
        da2hist = pickle.load(f)

    #lxbnd = np.fromfile(xbndPath)
    #lybnd = np.fromfile(ybndPath)

    lxbnd = dxbnd[xmetric]
    lybnd = dxbnd[ymetric]

    ##************************************
    ## Plot (xmetric vs ymetric)
    ##************************************
    dlabel = {'dprec':'Surface precipitation error (mm/h)'
             ,'ndprec':'Norm. surf. precip. error'
             ,'rmse':'Profile RMSE (g/m3)'
             ,'cc'  :'Profile Correlation Coef.'
             ,'dwatNorm' :'Total condensed water \n difference (Normed)'
             ,'dconvfrac':'Convective fraction error'
             ,'dstop-prof':'Storm top height error (km)\n(profile-based)'
             ,'dstop':'Storm top height error (km)\n(top-weighted)'
             ,'dtqv' :'Total water vapor error (kg/m2)'
             ,'dvfracConv' :'convective fraction error'
             }

    #for surf in lsurf:
    #    fig = plt.figure(figsize=(6,6))
    #    ax = fig.add_axes([0.25,0.20,0.7,0.7])
    #
    #    a1x = (lxbnd[:-1] + lxbnd[1:])*0.5
    #    for preclev in dprecrange.keys():
    #        if preclev ==-1: continue

    #        iprec,eprec = dprecrange[preclev]
    #        slabel = '%d-%dmm/h'%(iprec,eprec)
    #        linestyle= ['-', '--', '-'][preclev]
    #        linewidth= [1, 2, 2][preclev]
    #        #mycolor = ['darkblue','orange','crimson'][preclev]
    #        mycolor = ['k','k','k'][preclev]
    #
    #        #a1y = ma.masked_invalid(da1sum[surf,preclev] / da1num[surf,preclev])
    #        a1y = ma.masked_where(da1num[surf,preclev] <100, da1sum[surf,preclev]) / da1num[surf,preclev]

    #        ax.plot(a1x, a1y, linestyle=linestyle, linewidth=linewidth, color=mycolor, label=slabel)
    #
    #    stitle = '%s %s'%(string.upper(rettype), surf)
    #    plt.title(stitle, fontsize=20)
    #
    #    xlabel = dlabel[xmetric]
    #    ylabel = dlabel[ymetric]
    #
    #
    #    ax.set_ylabel(ylabel, fontsize=20)
    #    ax.set_xlabel(xlabel, fontsize=20)
    #    ax.legend(fontsize=20)
    #    ax.xaxis.set_tick_params(labelsize=16)
    #    ax.yaxis.set_tick_params(labelsize=16)
    #    ax.axvline(x=0, linestyle=':', color='gray', linewidth=3)
    #    ax.axhline(y=0, linestyle=':', color='gray', linewidth=3)
    #
    #    if xmetric=='rmse':
    #        xmin, xmax= 0, 0.7
    #    elif xmetric=='dwatNorm':
    #        xmin, xmax= -1.2, 5.5
    #    else:
    #        xmin, xmax= None,None
    #    ax.set_xlim([xmin,xmax])
    #
    #    figDir = '/home/utsumi/temp/ret'
    #    figPath = figDir + '/plot.%s.vs.%s.%s.%s.%s.png'%(xmetric,ymetric, rettype,surf,season)
    #    plt.savefig(figPath)
    #    plt.clf()
    #    plt.close()
    #    print figPath
    
    #************************************
    # Plot (xmetric vs num)
    #************************************
    '''
    for surf in lsurf:
        fig = plt.figure(figsize=(6,6))
        ax = fig.add_axes([0.25,0.20,0.7,0.7])
    
        a1x = (lxbnd[:-1] + lxbnd[1:])*0.5
        for preclev in dprecrange.keys():

            if preclev ==-1: continue

            iprec,eprec = dprecrange[preclev]
            slabel = '%d-%dmm/h'%(iprec,eprec)
            linestyle= ['-', '--', '-'][preclev]
            linewidth= [1, 2, 2][preclev]
            #mycolor = ['darkblue','orange','crimson'][preclev]
            mycolor = ['k','k','k'][preclev]
    
            a1y = ma.masked_invalid(da1num[surf,preclev]) 
    
    
            ax.plot(a1x, a1y, linestyle=linestyle, linewidth=linewidth, color=mycolor, label=slabel)
    
        stitle = '%s %s'%(string.upper(rettype),surf)
        plt.title(stitle, fontsize=20)
    
        xlabel = dlabel[xmetric]
        ylabel = 'count'
    
    
        ax.set_ylabel(ylabel, fontsize=20)
        ax.set_xlabel(xlabel, fontsize=20)
        ax.legend(fontsize=20)
        ax.xaxis.set_tick_params(labelsize=16)
        ax.yaxis.set_tick_params(labelsize=16)
        ax.axvline(x=0, linestyle=':', color='gray', linewidth=2)
        ax.axhline(y=0, linestyle=':', color='gray', linewidth=2)
    
        if xmetric=='rmse':
            xmin, xmax= 0, 0.7
        elif xmetric=='dwatNorm':
            xmin, xmax= -1.2, 5.5
        else:
            xmin, xmax= None,None
        ax.set_xlim([xmin,xmax])
    
        figDir = '/home/utsumi/temp/ret'
        figPath = figDir + '/plot.num-%s.vs.%s.%s.%s.%s.png'%(xmetric,ymetric, rettype,surf,season)
        plt.savefig(figPath)
        plt.clf()
        plt.close()
        print figPath
    '''
    #************************************
    # Plot 2D histogram + lines
    #************************************
    a1x = (lxbnd[:-1] + lxbnd[1:])*0.5
    for surf in lsurf:
        H  = da2hist[surf,-1]
        X,Y= np.meshgrid(lxbnd, lybnd)

        fig = plt.figure(figsize=(6,6))
        ax = fig.add_axes([0.25,0.20,0.45,0.45])
       
        vmax= H.max()
        im = ax.pcolormesh(X,Y,H, norm=matplotlib.colors.LogNorm(), cmap='jet', vmin=10, vmax=vmax)

        #--- Lines ---
        for preclev in dprecrange.keys():
            if preclev ==-1: continue

            iprec,eprec = dprecrange[preclev]
            slabel = '%d-%dmm/h'%(iprec,eprec)
            linestyle= ['-', '--', '-'][preclev]
            linewidth= [1, 2, 2][preclev]
            #mycolor = ['darkblue','orange','crimson'][preclev]
            mycolor = ['k','k','k'][preclev]
    
            a1y = ma.masked_where(da1num[surf,preclev] <100, da1sum[surf,preclev]) / da1num[surf,preclev]

            ax.plot(a1x, a1y, linestyle=linestyle, linewidth=linewidth, color=mycolor, label=slabel)
        #------------

        iprec,eprec = dprecrange[preclev]
 
        
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
        #-- PDF : hight
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
        stitle = '%s %s'%(string.upper(rettype), surf)
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
        figPath = figDir + '/plot.scatter.%s.vs.%s.%s.%s.%s.png'%(xmetric,ymetric, rettype,surf,season)
        fig.savefig(figPath)

        legPath = figDir + '/legend.plot.%s.vs.%s.%s.%s.%s.png'%(xmetric,ymetric, rettype,surf,season)
        figleg.savefig(legPath)


        plt.clf()
        plt.close()
        print figPath


    '''
    #************************************
    # Plot scatterplot (each prec-range)
    #************************************
    for surf in lsurf:
        for preclev in dprecrange.keys():

            H  = da2hist[surf,preclev]
            X,Y= np.meshgrid(lxbnd, lybnd)

            fig = plt.figure(figsize=(6,6))
            ax = fig.add_axes([0.25,0.20,0.55,0.55])
       
            vmax= H.max()
            im = ax.pcolormesh(X,Y,H, norm=matplotlib.colors.LogNorm(), cmap='jet', vmin=1, vmax=vmax)


            a1x = (lxbnd[:-1] + lxbnd[1:])*0.5
            iprec,eprec = dprecrange[preclev]
            slabel = '%d-%dmm/h'%(iprec,eprec)
            #mycolor = ['darkblue','orange','crimson'][preclev]
            a1y = ma.masked_where(da1num[surf,preclev] <100, da1sum[surf,preclev]) / da1num[surf,preclev]
            ax.plot(a1x, a1y, linestyle='-', linewidth=2, color='k')


            iprec,eprec = dprecrange[preclev]
 
            stitle = '%s %s %d-%dmm/h'%(string.upper(rettype),surf, iprec, eprec)

            plt.title(stitle, fontsize=20)
        
            xlabel = dlabel[xmetric]
            ylabel = dlabel[ymetric]
        
        
            ax.set_ylabel(ylabel, fontsize=20)
            ax.set_xlabel(xlabel, fontsize=20)

            ax.xaxis.set_tick_params(labelsize=16)
            ax.yaxis.set_tick_params(labelsize=16)
            ax.axvline(x=0, linestyle=':', color='gray', linewidth=2)
            ax.axhline(y=0, linestyle=':', color='gray', linewidth=2)
        
            if xmetric=='rmse':
                xmin, xmax= 0, 0.7
            elif xmetric=='dwatNorm':
                xmin, xmax= -1.2, 5.5
            else:
                xmin, xmax= None,None
            ax.set_xlim([xmin,xmax])

            #-- colorbar ---
            cax = fig.add_axes([0.82,0.15,0.02, 0.6])
            cbar=plt.colorbar(im, orientation='vertical', cax=cax)
            cbar.ax.tick_params(labelsize=16)
            #---------------

        
            figDir = '/home/utsumi/temp/ret'
            figPath = figDir + '/scatter.%s.vs.%s.%s.%s.%s.%d-%dmmh.png'%(xmetric,ymetric, rettype,surf,season,iprec,eprec)
            plt.savefig(figPath)
            plt.clf()
            plt.close()
            print figPath

    '''
