import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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

iYM = [2014,6]
eYM = [2014,6]
lYM = util.ret_lYM(iYM,eYM)
#calcflag = True
calcflag = False
lrettype = ['gprof']
DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
#metric= 'rmse'
metric= 'unsimilarity'

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



nsample = 18 # sampling orbits per day
#lprec = arange(0.25, 20+0.01, 0.5)
#lbiaslev = [0,1,2,3]
#dbiasrange = {0:[0,1],1:[1,2],2:[2,5],3:[5,999]}
dprecrange = {0:[1,3], 1:[5,7], 2:[10,12]}
aprecrange = dprecrange.values()
#lbias = np.arange(0.,20+0.1, 1.0)
lbias = np.arange(-20,20+0.1,1.0)
nz    = 25
#*****************************************
def taylor_index(a2ref, a2dat):
    a2mask1 = ma.masked_less(a2ref, 0).mask
    a2mask2 = ma.masked_less(a2dat, 0).mask
    a2mask  = a2mask1 + a2mask2
    a2ref = ma.masked_where(a2mask, a2ref)
    a2dat = ma.masked_where(a2mask, a2dat)

    a1num = (~a2mask).sum(axis=1)
    a1stdref = a2ref.std(axis=1)
    a1stddat = a2dat.std(axis=1)
    a1mref   = a2ref.mean(axis=1).reshape(-1,1)
    a1mdat   = a2dat.mean(axis=1).reshape(-1,1)

    a1cov    = ((a2ref - a1mref)*(a2dat - a1mdat)).sum(axis=1)/a1num
    a1corr = a1cov / (a1stdref * a1stddat)
    corrmax= 1.0

    S = 4*(1.0+a1corr)**4 /((a1stdref/a1stddat + a1stddat/a1stdref)**2) / (1.0+corrmax)**4
    S = ma.masked_invalid(S)
    return S


def calc_rmse(a2ref, a2dat):
    a2ref = ma.masked_less(a2ref,0)
    a2dat = ma.masked_less(a2dat,0)
    return np.sqrt( ((a2dat-a2ref)**2).mean(axis=1) )    


#*****************************************
for rettype in lrettype:

    #-- Initialize ------
    da1sum = {}
    da1sum2= {}
    da1num = {}
    for preclev in dprecrange.keys():
        da1sum [preclev]= np.zeros(len(lbias)-1, float32)
        da1sum2[preclev]= np.zeros(len(lbias)-1, float32)
        da1num [preclev]= np.zeros(len(lbias)-1, int32)
    #--------------------

    for Year,Mon in lYM:
        if calcflag ==False: continue
    
    
        eDay = calendar.monthrange(Year,Mon)[1]
        iDTime = datetime(Year,Mon,1)
        eDTime = datetime(Year,Mon,eDay)
        dDTime = timedelta(days=1)
        lDTime = util.ret_lDTime(iDTime, eDTime, dDTime)
        #lDTime = lDTime[:2]  # test 
        for DTime in lDTime:
            Year,Mon,Day = DTime.timetuple()[:3]
    
            pairDir = tankbaseDir + '/utsumi/PMM/validprof/pair/epc.%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)
            ssearch= pairDir + '/Latitude.*.npy'
            llatPath= sort(glob.glob(ssearch))
            #-- random sampling -----
            if len(llatPath)>=nsample:
                llatPath = random.sample(llatPath, nsample)
            else:
                pass
            #------------------------
            
            for latPath in llatPath:
                print latPath
                oid = int(latPath.split('.')[-2])
                a1lat = np.load(pairDir + '/Latitude.%06d.npy'%(oid))
                a1lon = np.load(pairDir + '/Longitude.%06d.npy'%(oid))
    
                a1precradTmp= np.load(pairDir+ '/precrad.%06d.npy'%(oid))
                a1precpmwTmp= np.load(pairDir+ '/precpmw.%06d.npy'%(oid))
                a2profradTmp= np.load(pairDir+ '/profrad.%06d.npy'%(oid))[:,2:nz] 
                a2profpmwTmp= np.load(pairDir+ '/profpmw.%06d.npy'%(oid))[:,2:nz]
                a1surftype= np.load(pairDir + '/surfaceTypeIndex.%06d.npy'%(oid))
                if rettype == 'gprof':
                    a1qflag   = np.load(pairDir + '/qualityFlag.%06d.npy'%(oid))
    
                #--- Check surface type --
                a1flagS1 = ma.masked_equal(a1surftype,1).mask   # Ocean
                a1flagS2 = ma.masked_inside(a1surftype,3,7).mask # Vegetation
                a1flagS3 = ma.masked_inside(a1surftype,12,13).mask # Standing water and rivers & Water/Coast boundary
                a1flagS  = a1flagS1 + a1flagS2 + a1flagS3
    
                #--- Check quality flag --
                if rettype == 'gprof':
                    a1flagQ = ma.masked_equal(a1qflag,0).mask   # Good quality (flag=0)
                    a1flag  = a1flagQ * a1flagS 
                else:
                    a1flag  = a1flagS   
             
                a1precpmwTmp = a1precpmwTmp[a1flag]
                a1precradTmp = a1precradTmp[a1flag]
                a2profpmwTmp = a2profpmwTmp[a1flag,:]
                a2profradTmp = a2profradTmp[a1flag,:] 
    
                #--- Mask missing data --
                a1precpmwTmp = ma.masked_less(a1precpmwTmp,0)
                a1precradTmp = ma.masked_less(a1precradTmp,0)
    
                #a2profpmwTmp = ma.masked_outside(a2profpmwTmp,0,30)
                a2profpmwTmp = ma.masked_less(a2profpmwTmp,0)
                a2profradTmp = ma.masked_less(a2profradTmp,0)
    
    
                if a1precpmwTmp.shape[0]==0:
                    continue
    
    
                ##--- Bias -------
                a1prbiasTmp = a1precpmwTmp - a1precradTmp
    
                #--- Classify by prec --
                for preclev in dprecrange.keys():
                    iprec,eprec = dprecrange[preclev]
                    a1flagP = ma.masked_inside(a1precradTmp, iprec, eprec).mask
   
                    if a1flagP.sum()==0: continue
 
                    a1prbiasTmp2 = a1prbiasTmp[a1flagP]
                    a2profpmwTmp2= a2profpmwTmp[a1flagP]
                    a2profradTmp2= a2profradTmp[a1flagP]
    
                    if metric=='rmse':
                        a1metricTmp2 = calc_rmse(a2profradTmp2, a2profpmwTmp2)
                    elif metric=='unsimilarity':
                        a1metricTmp2 = 1 - taylor_index(a2profradTmp2, a2profpmwTmp2)
                    else:
                        print 'check metric',metric
                        sys.exit()
   
                    #-- Bin by bias --
                    bins = lbias
                    a1sum,a1bnd,_ = stats.binned_statistic(a1prbiasTmp2, a1metricTmp2, statistic='sum', bins=bins)
                    a1num,a1bnd,_ = stats.binned_statistic(a1prbiasTmp2, a1metricTmp2, statistic='count', bins=bins)
                    a1sum2    = a1sum**2
   
                    a1sum  = ma.masked_invalid(a1sum).filled(0.0)
                    a1sum2 = ma.masked_invalid(a1sum2).filled(0.0)
                    a1num  = ma.masked_invalid(a1num).filled(0.0)
                   
                    da1sum [preclev] = da1sum [preclev] + a1sum
                    da1sum2[preclev] = da1sum2[preclev] + a1sum2
                    da1num [preclev] = da1num [preclev] + a1num
       
    #--- Save -----
    outDir  = tankbaseDir + '/utsumi/PMM/validprof/metrix.vs.prec/pickle.tmp'
    print outDir
    util.mk_dir(outDir)
    sumPath = outDir + '/sum.%s.binned.bfile'%(metric)
    numPath = outDir + '/num.%s.binned.bfile'%(metric) 
    sum2Path = outDir + '/sum2.%s.binned.bfile'%(metric)
    precrangePath = outDir + '/precrange.%s.bfile'%(metric)
    biasbndPath= outDir + '/bnd.bias.%s.npy'%(metric) 

    if calcflag == True:
        with open(sumPath, 'wb') as f:
            pickle.dump(da1sum, f)
        with open(numPath, 'wb') as f:
            pickle.dump(da1num, f)
        with open(sum2Path, 'wb') as f:
            pickle.dump(da1sum2, f) 
        with open(precrangePath, 'wb') as f:
            pickle.dump(dprecrange, f)

        lbias.tofile(biasbndPath)
        print sumPath

    #************************************
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

    lbias = np.fromfile(biasbndPath)

    #************************************
    # Plot
    #************************************
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_axes([0.15,0.15,0.7,0.7])

    a1x = (lbias[:-1] + lbias[1:])*0.5
    for preclev in dprecrange.keys():
        iprec,eprec = dprecrange[preclev]
        slabel = '%d-%dmm/h'%(iprec,eprec)
        linestyle= '-'
        mycolor = ['darkblue','orange','crimson'][preclev]

        a1y = ma.masked_invalid(da1sum[preclev] / da1num[preclev]) 
        ax.plot(a1x, a1y, linestyle=linestyle, color=mycolor, label=slabel)

    ylabel = 'Profile ' + {'rmse':'RMSE (g/m3)'
                          ,'unsimilarity':'Unsimilarity'}[metric]
    xlabel = 'Surface precip. absolute bias (mm/h)'
    ax.set_ylabel(ylabel, fontsize=15)
    ax.set_xlabel(xlabel, fontsize=15)
    ax.legend(fontsize=15)

    figDir = '/home/utsumi/temp/ret'
    figPath = figDir + '/plot.%s.vs.bias.%s.png'%(metric, rettype)
    plt.savefig(figPath)
    print figPath
    print 'sum'
    for i in range(a1x.shape[0]):
        print a1x[i],da1sum[2][i], da1num[2][i]

