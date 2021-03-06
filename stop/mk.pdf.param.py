import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import pylab as pl
import matplotlib.gridspec as gridspec
from glob import glob

import numpy.ma as ma
import sys,os
from datetime import datetime, timedelta
import myfunc.util as util
from collections import deque
from scipy.optimize import curve_fit
import scipy.stats as stats

import calendar
import socket
import pickle
#get_ipython().magic(u'matplotlib inline')


hostname = socket.gethostname()
if hostname == 'shui':
    stopbaseDir= '/tank/utsumi/PMM/stop'
    tankDir    = '/tank'
    figDir = '/home/utsumi/temp/stop'
elif hostname == 'well':
    stopbaseDir= '/home/utsumi/mnt/lab_tank/utsumi/PMM/stop'
    tankDir    = '/home/utsumi/mnt/lab_tank'
    #figDir = '/home/utsumi/mnt/lab_home_rainbow/utsumi/public_html/tempfig/stop'
    figDir = '/home/utsumi/temp/stop'

else:
    stopbaseDir= '/mnt/j/PMM/stop'
    figDir = '/mnt/c/ubuntu/fig'


Year  = 2017
season = 0
argvs = sys.argv
if len(argvs)>1:
    print argvs
    lsurf = argvs[1:]
else:
    #lisurf = range(1,14+1)
    lsurf = ['ocean','seaice','vege','snow','swater','coast','siedge']


lsurfTmp = []
for x in lsurf:
    if x.isdigit(): x = int(x)
    lsurfTmp.append(x)
lsurf = lsurfTmp
#******************************************
def ret_lisurf(surf):
    if surf in ['ocean','seaice','vege','snow','swater','coast','siedge']:
        if   surf=='ocean':  lisurf=[1]
        elif surf=='seaice': lisurf=[2]
        elif surf=='vege':   lisurf=[3,4,5,6,7]
        elif surf=='snow':   lisurf=[8,9,10,11]
        elif surf=='swater': lisurf=[12]
        elif surf=='coast':  lisurf=[13]
        elif surf=='siedge': lisurf=[14]
    elif surf in range(1,14+1):
        lisurf = [surf]
    else:
        print 'check surf',surf
        sys.exit()

    return lisurf



def read_var_collect(varName=None, lDTime=None, ldydx=None, surf=None, samplerate=None):

    lisurf = ret_lisurf(surf)
    a2var = deque([])
    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]

        for isurf in lisurf:

            a2varTmp = None
            for idydx,(dy,dx) in enumerate(ldydx):
                #srcDir = '/work/hk01/utsumi/PMM/stop/data/Tc/%04d/%02d/%02d'%(Year,Mon,Day)
                srcDir = stopbaseDir + '/data/%s/%04d/%02d/%02d'%(varName,Year,Mon,Day)
                srcPath=srcDir + '/%s.%ddy.%ddx.%02dsurf.npy'%(varName,dy,dx,isurf)
                if not os.path.exists(srcPath):
                    print 'No file',srcPath
                    continue
                avar = np.load(srcPath)
    
                if avar.shape[0]==0:
                    continue
    
                if a2varTmp is None:
                    a2varTmp = avar
                else:
                    a2varTmp = np.c_[a2varTmp, avar]
    
            #**********************
            # Resample
            #**********************
            if samplerate is not None:
                np.random.seed(0)  # Do not change !!
                aidx = np.random.choice(range(a2varTmp.shape[0]), int(a2varTmp.shape[0]*samplerate), replace=False)
    
                a2varTmp = a2varTmp[aidx] 
            #**********************
    
    
            if a2varTmp is None:
                continue
            else:
                a2varTmp = np.array(a2varTmp)
            #**********************
            a2var.extend(a2varTmp)
    return np.array(a2var)




def split2batchs(a,bsize):
    n=len(a)/bsize + 1
    return [a[i*bsize:(i+1)*bsize] for i in range(n)] 


def lognorm(x,m,s):
    return 1.0 / (np.sqrt(2*np.pi)*s*x) \
            * np.exp( -np.square(np.log(x) - m) / (2*s*s)    ) 

def fnorm(x,m,s):
    return 1.0 / np.sqrt(2*np.pi*s*s) \
            * np.exp( -np.square(x - m) / (2*s*s)    ) 


def twolognorm(x, amp1, m1, s1, amp2, m2, s2):
    y1 = amp1*1.0 / (np.sqrt(2*np.pi)*s1*x) \
            * np.exp( -np.square(np.log(x) - m1) / (2*s1*s1)    ) 
    y2 = amp2*1.0 / (np.sqrt(2*np.pi)*s1*x) \
            * np.exp( -np.square(np.log(x) - m1) / (2*s1*s1)    ) 

    return (y1 + y2)/(amp1+amp2)

def normlognorm(x, m1, s1, m2, s2):
    y1 = 1.0 / np.sqrt(2*np.pi*s*s) \
            * np.exp( -np.square(x - m) / (2*s*s)    ) 

    y2 = 1.0 / (np.sqrt(2*np.pi)*s1*x) \
            * np.exp( -np.square(np.log(x) - m1) / (2*s1*s1)    ) 

    return (y1 + y2)/2.0



def loggauss(x, *params):

    num_func = int(len(params)/3)

    y_list = []
    for i in range(num_func):
        y = np.zeros_like(x)
        param_range = list(range(3*i,3*(i+1),1))
        amp = params[int(param_range[0])]
        ctr = params[int(param_range[1])]
        wid = params[int(param_range[2])]
        y = y + amp/x * np.exp( -((np.log(x) - ctr)/wid)**2)

        y_list.append(y)

    y_sum = np.zeros_like(x)
    for i in y_list:
        y_sum = y_sum + i

    y_sum = y_sum + params[-1]

    return y_sum


def gauss(x, *params):

    num_func = int(len(params)/3)

    y_list = []
    for i in range(num_func):
        y = np.zeros_like(x)
        param_range = list(range(3*i,3*(i+1),1))
        amp = params[int(param_range[0])]
        ctr = params[int(param_range[1])]
        wid = params[int(param_range[2])]
        y = y + amp * np.exp( -((x - ctr)/wid)**2)
        y_list.append(y)

    y_sum = np.zeros_like(x)
    for i in y_list:
        y_sum = y_sum + i

    y_sum = y_sum + params[-1]

    return y_sum

#******************************************
for surf in lsurf:
    if surf==1:
        bsize = 10
        samplerate = 0.5
    else:
        bsize = 80
        samplerate = None

    wbnd = 200  # m
    bnds = np.arange(0,32000+1,wbnd)
    nbnds= len(bnds)
    mbnds= 0.5*(bnds[:-1] + bnds[1:])

    lDTime_train = np.load('/home/utsumi/bin/PMM/stop/ldtime-%04d-train.npy'%(Year), allow_pickle=True)


    llDTime = split2batchs(lDTime_train, bsize)

    ahist = np.zeros(nbnds-1).astype('int32')

    for klDTime,lDTime in enumerate(llDTime):
        a1stop = read_var_collect(varName='stop', lDTime=lDTime, ldydx=[[0,0]], surf=surf, samplerate=None)

        ahistTmp,_ = np.histogram(a1stop, bins=bnds)

      
        ahist = ahist + ahistTmp 

        #----- sampling data ----
        if klDTime == 0:
            asample = a1stop
        #------------------------

    #*********************
    # PDF
    #*********************
    apdf = ahist / float(ahist.sum()*wbnd)

    #*********************
    # Fitting
    #*********************
    f = lognorm
    (m,s),_ = curve_fit(f, mbnds, apdf)
    alognorm = lognorm(mbnds, m, s)
    print m, s
    #*********************

    # [Amp, center(m), width] + [Amp, center(m), width] + [background]
    dipara    = {}
    dipara[1] = [0.00012, 2000, 1000] + [0.00020, 2600, 1000] + [0.0005, 6000, 2000] + [0]
    dipara[2] = [0.00012, 2000, 1000] + [0.00020, 2600, 1000] + [0]
    dipara[3] = [0.0017, 3000,500] + [0.0002, 5000, 500] + [0.00005, 7000, 1000] + [0]
    dipara[4] = [0.00015, 3000, 1000] + [0.00015, 5000, 1000] + [0.00005, 7000, 2000] + [0]
    dipara[5] = [0.00012, 3000, 1000] + [0.00020, 5000, 1000] + [0.00003, 9000, 2000] + [0]
    dipara[6] = [0.00020, 5000, 1000] + [0.00020, 5000, 1000] + [0]
    dipara[7] = [0.00010, 3000, 500]  + [0.00025, 7000,500] + [0.00002, 9000, 2000] + [0]
    dipara[8] = [0.00020, 2000, 500]  + [0.00005, 3000,1000] + [0]
    dipara[9] = [0.00020, 2000, 500]  + [0.00005, 3000,1000] + [0]
    dipara[10]= [0.00020, 2000, 500]  + [0.00005, 3000,1000] + [0.00002,7000, 1000] + [0]
    dipara[11]= [0.00020, 2000, 500]  + [0.00005, 3000,1000] + [0.00002,8000, 2000] + [0]
    dipara[12]= [0.00020, 3000, 500]  + [0.00020, 5000,1000] + [0.00010,8000, 1000] + [0]
    dipara[13]= [0.00020, 2500, 500]  + [0.00010, 5000,500] + [0.00005,7000, 2000] + [0]
    dipara[14]= [0.00050, 2000, 500]  + [0.00030, 3000,500] + [0.00005, 4000, 2000] + [0]

    dipara['ocean'] = [0.00012, 2000, 1000] + [0.00020, 2600, 1000] + [0.0005, 6000, 2000] + [0]
    dipara['seaice']= [0.00012, 2000, 1000] + [0.00020, 2600, 1000] + [0]
    dipara['vege']  = [0.00012, 3000, 1000] + [0.00020, 5000, 1000] + [0.00003, 9000, 2000] + [0]
    dipara['snow']  = [0.00020, 2000, 500]  + [0.00005, 3000,1000] + [0.00005, 6000, 2000] + [0]

    dipara['swater']= [0.00020, 3000, 500]  + [0.00020, 5000,1000] + [0.00010,8000, 1000] + [0]
    dipara['coast'] = [0.00020, 2500, 500]  + [0.00010, 5000,500] + [0.00005,7000, 2000] + [0]
    dipara['siedge']= [0.00050, 2000, 500]  + [0.00030, 3000,500] + [0.00005, 4000, 2000] + [0]





    func = gauss
    p0 = dipara[surf]
    para, _ = curve_fit(func, mbnds, apdf, p0=p0)
    print 
    print para
    afit  = func(mbnds, *para)
    pdfmax = afit.max()
    para   = np.concatenate([para, np.array([pdfmax])])
    print para
    print type(para)
    #*********************
    paramDir = tankDir + '/utsumi/PMM/stop/pdfparam'
    paramPath= paramDir + '/gauss.param.surf-%s.npy'%(surf)
    np.save(paramPath, np.array(para))
    print paramPath

    #*********************
    figDir = '/home/utsumi/temp/stop'
    figPath = figDir + '/hist.surf-%s.png'%(surf)
    fig = plt.figure(figsize=(4,4))
    ax  = fig.add_axes([0.2,0.2,0.6,0.6])
    x   = 0.5*(bnds[:-1] + bnds[1:]) *0.001

    #ax.plot(x, ahist)
    ax.plot(x, apdf, '-', color='k')
    ax.plot(x, alognorm, '-',color='r')
    ax.plot(x, afit, '-',color='b')

    ax.set_xlim([0,16])
    plt.savefig(figPath)
    print figPath

    #******************
    # PDF




 
