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
import string

season = 'JJADJF'
calcflag = True
#calcflag = False
#lrettype = ['epc']
lrettype = ['epc','gprof']
#lrettype = ['gprof']
DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
#metric= 'rmse'
metric= 'cc'

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
dprecrange = {0:[1,3], 1:[5,7], 2:[10,12]}
#dprecrange = {2:[10,12]}
aprecrange = dprecrange.values()
#lbias = np.arange(0.,20+0.1, 1.0)
lbias = np.arange(-20,20+0.1,1.0)
#lbias = np.arange(-10,10+0.1,1.0)
nz    = 25
#*****************************************
def calc_rmse(a2ref, a2dat):
    a2ref = ma.masked_less(a2ref,0)
    a2dat = ma.masked_less(a2dat,0)
    return np.sqrt( ((a2dat-a2ref)**2).mean(axis=1) )    

def calc_cc(x,y,axis):
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


def ret_lmon(season):
    if season=='JJADJF':
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


#*****************************************
lYM = ret_lmon(season)
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

                    a1flagS = ma.masked_equal(np.load(pairDir + '/surfaceTypeIndex.%06d.npy'%(oid)) ,2)
                    a1flagS = ma.masked_inside(a1flagS,8,11)
                    a1flagS = ma.masked_equal(a1flagS,14)
                    a1flagS = ~(a1flagS.mask)


    
                ##--- Check surface type --
                #a1flagS1 = ma.masked_equal(a1surftype,1).mask   # Ocean
                #a1flagS2 = ma.masked_inside(a1surftype,3,7).mask # Vegetation
                #a1flagS3 = ma.masked_inside(a1surftype,12,13).mask # Standing water and rivers & Water/Coast boundary
                #a1flagS  = a1flagS1 + a1flagS2 + a1flagS3
    
                #--- Screen missing surface precipitation ---
                a1flagP1 = ma.masked_greater_equal(a1precpmwTmp,0).mask
                a1flagP2 = ma.masked_greater_equal(a1precradTmp,0).mask
                a1flagP  = a1flagP1 * a1flagP2 

                #--- Check quality flag --
                if rettype == 'gprof':
                    a1flagQ = ma.masked_equal(a1qflag,0).mask   # Good quality (flag=0)
                    a1flag  = a1flagS * a1flagP *a1flagQ
                else:
                    a1flag  = a1flagP 
            

                a1precpmwTmp = a1precpmwTmp[a1flag]
                a1precradTmp = a1precradTmp[a1flag]
                a2profpmwTmp = a2profpmwTmp[a1flag,:]
                a2profradTmp = a2profradTmp[a1flag,:] 
    
                #--- Mask missing data --
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
                    elif metric=='cc':
                        a1metricTmp2 = calc_cc(a2profradTmp2, a2profpmwTmp2, axis=1)

                    else:
                        print 'check metric',metric
                        sys.exit()
   
                    #-- Bin by bias --
                    bins = np.sort(lbias)
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
    outDir  = tankbaseDir + '/utsumi/PMM/validprof/metrix.vs.prec/pickle'
    print outDir
    util.mk_dir(outDir)
    sumPath = outDir + '/sum.%s.%s.%s.binned.bfile'%(metric,rettype,season)
    numPath = outDir + '/num.%s.%s.%s.binned.bfile'%(metric,rettype,season) 
    sum2Path = outDir + '/sum2.%s.%s.%s.binned.bfile'%(metric,rettype,season)
    precrangePath = outDir + '/precrange.%s.%s.%s.bfile'%(metric,rettype,season)
    biasbndPath= outDir + '/bnd.bias.%s.%s.%s.npy'%(metric,rettype,season) 

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

    lbias = np.fromfile(biasbndPath)

    #************************************
    # Plot
    #************************************
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_axes([0.2,0.15,0.7,0.7])

    a1x = (lbias[:-1] + lbias[1:])*0.5
    for preclev in dprecrange.keys():
        iprec,eprec = dprecrange[preclev]
        slabel = '%d-%dmm/h'%(iprec,eprec)
        linestyle= ['-', '--', '-'][preclev]
        linewidth= [1, 2, 2][preclev]
        #mycolor = ['darkblue','orange','crimson'][preclev]
        mycolor = ['k','k','k'][preclev]

        a1y = ma.masked_invalid(da1sum[preclev] / da1num[preclev]) 
        ax.plot(a1x, a1y, linestyle=linestyle, linewidth=linewidth, color=mycolor, label=slabel)

    stitle = '%s'%(string.upper(rettype))
    plt.title(stitle, fontsize=20)

    ylabel = 'Profile ' + {'rmse':'RMSE (g/m3)'
                          ,'cc'  :'Correlation Coef.'}[metric]
    xlabel = 'Surface precip. bias (mm/h)'
    ax.set_ylabel(ylabel, fontsize=20)
    ax.set_xlabel(xlabel, fontsize=20)
    ax.legend(fontsize=20)
    ax.xaxis.set_tick_params(labelsize=16)
    ax.yaxis.set_tick_params(labelsize=16)
    ax.axvline(x=0, linestyle=':', color='gray', linewidth=2)

    figDir = '/home/utsumi/temp/ret'
    figPath = figDir + '/plot.%s.vs.bias.%s.%s.png'%(metric, rettype,season)
    plt.savefig(figPath)
    print figPath

