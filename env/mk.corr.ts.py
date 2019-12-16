import matplotlib
matplotlib.use('Agg')
from numpy import *
from datetime import datetime, timedelta
import myfunc.util as util
import os, sys, glob, socket
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
#from collections import deque
import calendar
import pickle

calcflag = True
#calcflag = False
lseason = [7]

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


lprtype= ['conv']   # 'conv', 'strat'
#lheight= ['high','low','all']
lheight= ['all']
lts = ['cold','cool','warm','hot']
dts = {'cold':[-999,10], 'cool':[10,20], 'warm':[20,26], 'hot':[26,999]}
lsurf = ['land','sea']

#lvar = ['pre.dtvdz_low','pre.dtvdz_mid']
#lvar = ['dtv_low','dtv_mid']
#lvar = ['w_1.5','w_4.5','w_7.5']
#lvar = ['r_15','r_45','r_75']
lvar = ['r_75']
#lvar = ['dept_mid']
#lvar = ['dept_low']
#lvar = ['mvimd']
#lvar = ['tcwv']
#lvar = ['skt','r_15','r_45','dtv_low','dtv_mid','dept_low','dept_mid']
#lvar = ['cape']
#lvar = ['pre.deptdz_low','pre.deptdz_mid']

#lvar = ['w_1.5','w_4.5','w_7.5']
#lvar = lvar+['r_1.5','r_4.5','r_7.5']
#lvar = lvar+['mvimd','tcwv']
#lvar = lvar+['cape','pre.cape']
#lvar = ['pre.deptdz_low','pre.deptdz_mid']
#lvar = ['pre.dtvdz_low','pre.dtvdz_mid']
#lvar = ['deptdz_low','deptdz_mid']
#lvar = lvar+['dtvdz_low','dtvdz_mid']
#lvar = ['cape']

ldtenv  = [-1]

lsurftype = ['sea','land']
#lsurftype = ['sea']
#lsurftype = ['land']
miss_out = -9999.

iDTime = datetime(2017,7,1)
eDTime = datetime(2017,7,30)
lDTime = util.ret_lDTime(iDTime,eDTime, timedelta(days=1))

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





#***********************8
for season in lseason:
    for var in lvar:
        for prtype in lprtype:
            for dt in ldtenv:
                if calcflag is not True: continue
                #****************************
                # Initialize
                #----------------------------
                dvar = {}
                dfrac= {}
                for height in lheight:
                    for ts in lts:
                        for surf in lsurf:
                            dvar [height,ts,surf] = np.array([])
                            dfrac[height,ts,surf] = np.array([]) 
                #----------------------------
                lDTime = season_lDTime(season)
                #lDTime = lDTime[:5] 
                for DTime in lDTime:
                    print DTime
                    Year,Mon,Day = DTime.timetuple()[:3]
                    convDir = pairbaseDir + '/cvfrac/%04d/%02d/%02d'%(Year,Mon,Day)
                    stratDir= pairbaseDir + '/svfrac/%04d/%02d/%02d'%(Year,Mon,Day)
                    varDir  = pairbaseDir + '/%s/%04d/%02d/%02d'%(var,Year,Mon,Day)
                    surfDir = pairbaseDir + '/%s/%04d/%02d/%02d'%('landSurfaceType',Year,Mon,Day)
                    stopDir = pairbaseDir + '/%s/%04d/%02d/%02d'%('heightStormTop',Year,Mon,Day)
                    tsDir   = pairbaseDir + '/%s/%04d/%02d/%02d'%('skt',Year,Mon,Day)               
                    lconvPath = np.sort(glob.glob(convDir + '/*.npy'))
                
                    for convPath in lconvPath:
                        oid = int(convPath.split('.')[-2])
        
                        if int(oid) in lskipoid:
                            print 'skip',oid
                            continue
         
                        #** Read variable ****
                        if len(var.split('_'))==1:
                            a1varTmp = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.%03dh.00.0km.%06d.npy'%(var,Year,Mon,Day,var,dt,oid))

                            #try:
                            #    a1aveTmp = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/ave.%03dh.00.0km.%06d.npy'%(var,Year,Mon,Day,0,oid))
                            #    a1stdTmp = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/std.%03dh.00.0km.%06d.npy'%(var,Year,Mon,Day,0,oid))
                            #except IOError:
                            #    continue
                        elif var.split('_')[1] in ['low','mid']:
                            a1varTmp = ret_var(var)

                        else:
                            varName, lev = var.split('_')
                            lev      = float(lev)/10. 
                            a1varTmp = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.%03dh.%04.1fkm.%06d.npy'%(varName,Year,Mon,Day,varName,dt,lev,oid))

                            #try:
                            #    a1aveTmp = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/ave.%03dh.%04.1fkm.%06d.npy'%(varName,Year,Mon,Day,0,lev,oid))
                            #    a1stdTmp = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/std.%03dh.%04.1fkm.%06d.npy'%(varName,Year,Mon,Day,0,lev,oid))
                            #except IOError:
                            #    continue

                        #--------------------- 
                        #a1aveTmp = ma.masked_equal(a1aveTmp,-9999)
                        #a1stdTmp = ma.masked_equal(a1stdTmp,-9999)
                        #a1varTmp = a1varTmp - a1aveTmp
 

                        a1surfTmp= np.load(surfDir + '/landSurfaceType.00.0km.%06d.npy'%(oid))
                        a1stopTmp= np.load(stopDir + '/heightStormTop.00.0km.%06d.npy'%(oid))
                        a1tsTmp  = np.load(tsDir   + '/skt.00.0km.%06d.npy'%(oid))


        
                        if prtype =='conv':
                            a1fracTmp = np.load(convDir  + '/cvfrac.00.0km.%06d.npy'%(oid))
                        elif prtype == 'strat':
                            a1fracTmp = np.load(stratDir + '/svfrac.00.0km.%06d.npy'%(oid))
                        else:
                            print 'check prtype',prtype
                            sys.exit() 
                        #****************
                        # Flags
                        #----------------
                
                        d1flagh = {}
                        d1flagh['high'] = ma.masked_greater(a1stopTmp,6000).mask
                        d1flagh['mid' ] = ma.masked_inside(a1stopTmp,4000,6000).mask
                        d1flagh['low' ] = ma.masked_less(a1stopTmp,4000).mask
                        d1flagh['all' ] = np.array([True])
                
                        d1flagt = {}
                        for ts in lts:
                            ts0,ts1 = np.array(dts[ts]) + 273.15
                            d1flagt[ts]    = ma.masked_inside(a1tsTmp, ts0, ts1).mask

                
                        d1flags = {}
                        d1flags['land']= ma.masked_inside(a1surfTmp,0,99).mask
                        d1flags['sea'] = ma.masked_inside(a1surfTmp,100,199).mask
        
                        a1flagvar = ma.masked_not_equal(a1varTmp,-9999.).mask
        
                        lkey = [[height,ts,surf]
                                for height in lheight
                                for ts     in lts
                                for surf   in lsurf
                                ]
                        for (height,ts,surf) in lkey:
                            a1flagh = d1flagh[height]
                            a1flagt = d1flagt[ts]
                            a1flags = d1flags[surf]
        
                            a1flag = a1flagvar * a1flagh * a1flagt * a1flags
                            #*******************
                            # Stack
                            #-------------------
                            dvar[height,ts,surf] = np.concatenate([dvar[height,ts,surf], a1varTmp[a1flag]])
        
                            dfrac[height,ts,surf] = np.concatenate([dfrac[height,ts,surf], a1fracTmp[a1flag]])
    
            #**********************
            # Save
            #----------------------
            varPath  = pickleDir + '/stack.var.%s.dt%03d.%s.bfile'%(var,dt,season)
            fracPath = pickleDir + '/stack.frac%s.%s.%s.dt%03d.bfile'%(prtype,var,dt,season)
    
            if calcflag == True:
                with open(varPath, 'wb') as f:
                    pickle.dump(dvar, f)
            
                with open(fracPath, 'wb') as f:
                    pickle.dump(dfrac, f)

            #*********************
            # Load
            #---------------------
            with open(varPath, 'rb') as f:
                dvar  = pickle.load(f)
            
            with open(fracPath, 'rb') as f:
                dfrac = pickle.load(f)
     
            #**********************
            #----------------------
            lkey = [[height,ts,surf]
                    for height in lheight
                    for ts     in lts
                    for surf   in lsurf
                    ]
            for (height,ts,surf) in lkey:
                print height, ts, surf
                a1var  = dvar [height,ts,surf]
                a1frac = dfrac[height,ts,surf]
                if a1var.shape[0]==0: continue
                #-- Histograms --------
                xmin  = np.percentile(a1var,10)
                xmax  = np.percentile(a1var,90)
                xbins  = np.linspace(xmin,xmax,20)
                ybins  = np.linspace(0,1,20)
                
                H,xedges,yedges = np.histogram2d(a1var, a1frac, bins = [xbins,ybins])
                H = H.T
                H = ma.masked_equal(H,0)
                X,Y = np.meshgrid(xbins,ybins)
                
                #-- Binned average -----
                a1fracmean,_,_ = stats.binned_statistic(a1var, a1frac, bins = xbins)
                a1xcenter = 0.5*(xbins[:-1]+xbins[1:])
                
                #*** Prep for figure **
                if   var in ['dept_low','dept_mid']:
                
                    lev = var.split('_')[1]
                    varlongname = r'$d\theta$' +'e/dz %s'%(lev)
                
                    sunit = 'K/km'
                
                elif   var in ['dtv_low','dtv_mid','dtv_low','dtv_mid']:
                
                    lev = var.split('_')[1]
                    varlongname = 'Tv/dz %s'%(lev)
                
                    sunit = 'K/km'
                
                elif var in ['skt']:
                    varlongname='Skin temperature'
                    sunit = 'K'
                    
                elif var in ['cape','pre.cape']:
                    varlongname=var
                    sunit = 'J/kg'
                elif var in ['tcwv']:
                    varlongname='Total column water vapor'
                    sunit = 'kg/m2'
                elif var in ['mvimd']:
                    varlongname='Total moisture divergence'
                    sunit = 'kg/m2/s'
                elif var.split('_')[0] == 'w':
                    lev = float(var.split('_')[1])
                    varlongname='Upward velocity (%3.1fkm)'%(lev)
                    sunit = 'm/s'
                elif var.split('_')[0] == 'r':
                    lev = float(var.split('_')[1])
                    varlongname='Relative Humid. (%3.1fkm)'%(lev)
                    sunit = '%'
                
                
                else:
                    print 'check var in varlongname and sunit',var
                    sys.exit() 
                
                corr = np.corrcoef(a1var, a1frac)[0][1]
                #*** Figure ******
                fig = plt.figure(figsize=(4,4))
                ax  = fig.add_axes([0.17,0.15,0.65,0.65])
                
                im  = ax.pcolormesh(X,Y,H, cmap='jet', norm=matplotlib.colors.LogNorm())
                
                #-- binned mean convective ratio --
                ax.plot(a1xcenter, a1fracmean, '-',color='k')
                
                #----------------------------------
                ax.set_xlabel(sunit, fontsize=20)
                ax.set_ylabel('convective ratio', fontsize=20)
                
                #stitle = '%s %s'%(varlongname,surftype)
                #stitle = varlongname + ' ' + surf
                stitle = '%s %s %s %s'%(var, height, ts, surf)
                stitle = stitle + '\n' + 'corr=%.2f'%(corr)
                plt.title(stitle, fontsize=19)
                
                cax = fig.add_axes([0.84,0.17,0.02, 0.6])
                cbar=plt.colorbar(im, orientation='vertical', cax=cax)
                cbar.ax.tick_params(labelsize=16)
                
                
                
                figPath = figDir + '/scatter.frac.%s.%s.dt%03d.%s.%s.%s.%s.png'%(prtype,var,dt,season,height,ts,surf)
                plt.savefig(figPath)
                print figPath
                plt.clf() 
                
                print 'corr=',corr
                #-- write corrcoef ---
                corrDir = figDir + '/txt'
                util.mk_dir(corrDir)
                corrPath= corrDir + '/cc.%s.%s.dt%03d.%s.%s.%s.%s.txt'%(prtype,var,dt,season,height,ts,surf)
                f=open(corrPath,'w');f.write('%.2f'%(corr)); f.close()
    
