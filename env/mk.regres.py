import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
from datetime import datetime, timedelta
import myfunc.util as util
import os, sys, glob, socket
import numpy as np
import calendar
import pickle
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
calcflag = True
#calcflag = False
#season = 7
iDTime = datetime(2017,7,1)
eDTime = datetime(2017,7,31)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

dprver = 'V06'
dprverfull='V06A'
myhost = socket.gethostname()
if myhost =='shui':
    pairbaseDir = '/tank/utsumi/env/pair'
    pickleDir   = '/tank/utsumi/env/pickle'
    #figDir = '/home.rainbow/utsumi/public_html/tempfig'
    figDir = '/home/utsumi/temp/env'
elif myhost == 'well':
    pairbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/env/pair'
    pickleDir = '/home/utsumi/mnt/lab_tank/utsumi/env/pickle'
    #figDir = '/home/utsumi/mnt/lab_home_rainbow/public_html/tempfig'
    figDir = '/home/utsumi/temp/env'
else:
    print 'check myhost'
    sys.exit()

#targetvar = 'cvfrac'
targetvar = 'heightStormTop'
lts = ['cold','cool','warm','hot']
dts = {'cold':[-999,10+273.15], 'cool':[10+273.15,20+273.15], 'warm':[20+273.15,26+273.15], 'hot':[26+273.15,999]}
lsurf = ['land','sea']

lvar = ['dtv_low','dtv_mid','dept_low','dept_mid','r_15','r_45','r_75','skt']
#lvar = ['r_15','r_45','r_75','skt']

dt  = -1

lsurftype = ['sea','land']
#lsurftype = ['sea']
#lsurftype = ['land']
miss_out = -9999.

lskipoid = [
 19104  # 2017/7/7 UTC17-23, invalid era5 t 
,19105
,19106
,19107
,19108
,19109
]
#***********************8
def season_lDTime(season):
    Year = 2017
    lMon = util.ret_lmon(season)
    lDTime = []
    for Mon in lMon:
        eDay = calendar.monthrange(Year,Mon)[1]
        iDTime = datetime(Year,Mon,1)
        eDTime = datetime(Year,Mon,eDay)
        lDTimeTmp = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
        lDTime = lDTime + lDTimeTmp
    return lDTime

def ret_var(var):
    varName,lev = var.split('_')
    if lev in ['low','mid']:
        if lev=='low':
            lev0 = 1.5
            lev1 = 4.5
        elif lev=='mid':
            lev0 = 4.5
            lev1 = 7.5

        varNameIn = varName[1:]
        srcPath0 = pairbaseDir + '/%s/%04d/%02d/%02d/%s.%03dh.%04.1fkm.%06d.npy'%(varNameIn,Year,Mon,Day,varNameIn,dt,lev0,oid)
        srcPath1 = pairbaseDir + '/%s/%04d/%02d/%02d/%s.%03dh.%04.1fkm.%06d.npy'%(varNameIn,Year,Mon,Day,varNameIn,dt,lev1,oid)

        a1var0 = np.load(srcPath0)
        a1var1 = np.load(srcPath1)
        a1var = (ma.masked_less(a1var1,0) - ma.masked_less(a1var0,0))/(lev1-lev0)


    else:
        lev = float(lev)/10
        srcPath = pairbaseDir + '/%s/%04d/%02d/%02d/%s.%03dh.%04.1fkm.%06d.npy'%(varName,Year,Mon,Day,varName,dt,lev,oid)
        a1var = np.load(srcPath)
    return a1var

#*************************
# Make PCA coefficients
#*************************
if calcflag == True:
    #lDTime = season_lDTime(season)
    X  = None
    for var in lvar:
    
        a1var = np.array([])
    
        for DTime in lDTime:
            Year,Mon,Day = DTime.timetuple()[:3]
            convDir = pairbaseDir + '/cvfrac/%04d/%02d/%02d'%(Year,Mon,Day)
            lconvPath = np.sort(glob.glob(convDir + '/*.npy'))
        
            for convPath in lconvPath:
                oid = int(convPath.split('.')[-2])
            
                if int(oid) in lskipoid:
                    print 'skip',oid
                    continue
        
                #** Read variable ****
                if len(var.split('_'))==1:
                    varName  = var
                    lev      = 0
                    a1varTmp = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.%03dh.%04.1fkm.%06d.npy'%(varName,Year,Mon,Day,varName,dt,lev,oid))
    
                elif var.split('_')[1] in ['low','mid']:
                    a1varTmp = ret_var(var)
    
                else:
                    varName, lev = var.split('_')
                    lev      = float(lev)/10. 
                    a1varTmp = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.%03dh.%04.1fkm.%06d.npy'%(varName,Year,Mon,Day,varName,dt,lev,oid))
    
                a1var = np.concatenate([a1var, a1varTmp])
    
        #** Replace invalid data *****
        varName = var.split('_')[0] 
        if varName in ['dtv','dept']:
            a1var = ma.masked_outside(a1var, -30, 30).filled(miss_out)
       
        elif varName in ['r']:
            a1var = ma.masked_inside(a1var, -100, 0).filled(0)
            a1var = ma.masked_inside(a1var, 100, 200).filled(100)
             
    
        #*****************************
        if X is None:
            X = a1var
        else:
            X = np.vstack([X, a1var])
    
    X = X.T
    #** Mask ****
    a1flag = ~ma.masked_equal(X,miss_out).mask.any(axis=1)
    X  = X[a1flag,:]
    #** Normalize **
    xave = X.mean(axis=0)
    xstd = X.std(axis=0)
    
    X  = (X-xave) / xstd
    
    pca = PCA(n_components= len(lvar))
    pca.fit(X)
    egvec = pca.components_  # eigen vector for ith PC = egvec[i]
    egval = pca.explained_variance_ratio_
    
    print egvec.shape 

#*****************************************
# Save PCA coefficients
#-----------------------------------------
evecPath = pickleDir + '/pca.era.evec.%dvars.npy'%(len(lvar))
evalPath = pickleDir + '/pca.era.eval.%dvars.npy'%(len(lvar))
avePath = pickleDir + '/pca.era.ave.%dvars.npy'%(len(lvar))
stdPath = pickleDir + '/pca.era.std.%dvars.npy'%(len(lvar))

if calcflag==True:
    np.save(evecPath, egvec) 
    np.save(evalPath, egval)
    np.save(avePath , xave)
    np.save(stdPath , xstd)

#*****************************************
# Load PCA coefficients
#-----------------------------------------
egvec = np.load(evecPath) 
egval = np.load(evalPath)
xave  = np.load(avePath )
xstd  = np.load(stdPath )
 
#*****************************************
# Regression with princibal components 
#*****************************************
# Read environmental variables
X = None
for var in lvar:
    a1var = np.array([])

    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]
        convDir = pairbaseDir + '/cvfrac/%04d/%02d/%02d'%(Year,Mon,Day)
        lconvPath = np.sort(glob.glob(convDir + '/*.npy'))
    
        for convPath in lconvPath:
            oid = int(convPath.split('.')[-2])
        
            if int(oid) in lskipoid:
                print 'skip',oid
                continue
    
            #** Read variable ****
            if len(var.split('_'))==1:
                varName  = var
                lev      = 0
                a1varTmp = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.%03dh.%04.1fkm.%06d.npy'%(varName,Year,Mon,Day,varName,dt,lev,oid))

            elif var.split('_')[1] in ['low','mid']:
                a1varTmp = ret_var(var)

            else:
                varName, lev = var.split('_')
                lev      = float(lev)/10. 
                a1varTmp = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.%03dh.%04.1fkm.%06d.npy'%(varName,Year,Mon,Day,varName,dt,lev,oid))

            a1var = np.concatenate([a1var, a1varTmp])

    #** Replace invalid data *****
    varName = var.split('_')[0] 
    if varName in ['dtv','dept']:
        a1var = ma.masked_outside(a1var, -30, 30).filled(miss_out)
   
    elif varName in ['r']:
        a1var = ma.masked_inside(a1var, -100, 0).filled(0)
        a1var = ma.masked_inside(a1var, 100, 200).filled(100)
         

    #*****************************
    if X is None:
        X = a1var
    else:
        X = np.vstack([X, a1var])
X = X.T
#-------------------------------
# Read other variables
#-------------------------------
d1var = {}
lev = 0
for varName in ['skt','landSurfaceType'] + [targetvar]:
    a1var = np.array([])
    print lDTime
    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]
        targetDir = pairbaseDir + '/%s/%04d/%02d/%02d'%(targetvar,Year,Mon,Day)
        ltargetPath = np.sort(glob.glob(targetDir + '/*.npy'))
   
        print ltargetPath 
        for targetPath in ltargetPath:
            oid = int(targetPath.split('.')[-2])
        
            if int(oid) in lskipoid:
                print 'skip',oid
                continue

            if varName=='skt':     
                srcPath = pairbaseDir + '/%s/%04d/%02d/%02d/%s.%03dh.%04.1fkm.%06d.npy'%(varName,Year,Mon,Day,varName,dt,lev,oid)
            else:
                srcPath = pairbaseDir + '/%s/%04d/%02d/%02d/%s.%04.1fkm.%06d.npy'%(varName,Year,Mon,Day,varName,lev,oid)
            a1tmp = np.load(srcPath)

            a1var = np.concatenate([a1var, a1tmp])

    d1var[varName] = a1var

#-------------------------------
# Screen invalid data
#-------------------------------
a1flag = ~ma.masked_equal(X,miss_out).mask.any(axis=1)
for varName in ['skt','landSurfaceType'] + [targetvar]:
    a1flag = a1flag * ma.masked_not_equal(d1var[varName], miss_out).mask


X = X[a1flag,:]
for varName in ['skt','landSurfaceType'] + [targetvar]:
    d1var[varName] = d1var[varName][a1flag]


#-------------------------------
# Regression
#-------------------------------
X = (X-xave) / xstd

for surf in lsurf:
    if surf =='land':
        a1flags = ma.masked_inside(d1var['landSurfaceType'],0,99).mask
    elif surf=='sea':
        a1flags = ma.masked_inside(d1var['landSurfaceType'],100,199).mask
    
    else:
        print 'check surf',surf
        sys.exit()

    for ts in lts:
        ts0,ts1 = dts[ts]

        a1flagTmp = ma.masked_inside(d1var['skt'], ts0, ts1).mask
        a1flag  = a1flags * a1flagTmp

        if a1flag.sum() <10:
            print 'No enough events',surf,ts
            continue
        #------------------
        xtmp = X[a1flag,:]
        
        xtmp = np.dot(xtmp, egvec.T) # PCA transformation  (nsample, ncomponents)

        atarget= d1var[targetvar][a1flag]

    
        reg = LinearRegression()
        reg = reg.fit(xtmp, atarget.reshape(-1,1))   

        #print 'coef',reg.coef_
        #print 'intercept',reg.intercept_
        apred = reg.predict(xtmp).flatten()

        print np.corrcoef(atarget, apred)
        #--- Plot ---
        fig = plt.figure(figsize=(4,4))
        ax  = fig.add_axes([0.1,0.1,0.8,0.8])

        #ax.scatter(atarget, apred)
        hb =ax.hexbin(atarget, apred, gridsize=30, norm=matplotlib.colors.LogNorm())
        fig.colorbar(hb)
        stitle = '%s %s %s'%(targetvar, surf, ts)
        plt.title(stitle)
        figPath = figDir + '/scatter.pred-era.%s.%s.%s.png'%(targetvar,surf,ts)
        plt.savefig(figPath)
        plt.clf()
        print figPath
