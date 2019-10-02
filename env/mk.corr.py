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

calcflag = True:
lseason = [7]

dprver = 'V06'
dprverfull='V06A'
myhost = socket.gethostname()
if myhost =='shui':
    pairbaseDir = '/tank/utsumi/env/pair'
    pickleDir   = '/tank/utsumi/env/pickle'
    figDir = '/home.rainbow/utsumi/public_html/tempfig'
elif myhost == 'well':
    pairbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/env/pair'
    pickleDir = '/home/utsumi/mnt/lab_tank/utsumi/env/pickle'
    figDir = '/home/utsumi/mnt/lab_home_rainbow/public_html/tempfig'
else:
    print 'check myhost'
    sys.exit()


lprtype= ['conv']   # 'conv', 'strat'
lheight= ['high','low']
lregion= ['t','s','m']
lsurf  = ['land','sea']


#lvar = ['pre.dtvdz_low','pre.dtvdz_mid']
#lvar = ['dtvdz_low','dtvdz_mid']
#lvar = ['w_1.5','w_4.5','w_7.5']
#lvar = ['r_1.5','r_4.5','r_7.5']
#lvar = ['r_1.5']
#lvar = ['mvimd']
#lvar = ['tcwv']
#lvar = ['cape','pre.cape']
#lvar = ['pre.deptdz_low','pre.deptdz_mid']
#lvar = ['deptdz_low','deptdz_mid']

#lvar = ['w_1.5','w_4.5','w_7.5']
#lvar = lvar+['r_1.5','r_4.5','r_7.5']
#lvar = lvar+['mvimd','tcwv']
#lvar = lvar+['cape','pre.cape']
#lvar = ['pre.deptdz_low','pre.deptdz_mid']
#lvar = ['pre.dtvdz_low','pre.dtvdz_mid']
lvar = ['deptdz_low','deptdz_mid']
lvar = lvar+['dtvdz_low','dtvdz_mid']



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
    if var in ['deptdz_low','deptdz_mid','pre.deptdz_low','pre.deptdz_mid', 'dtvdz_low','dtvdz_mid','pre.dtvdz_low','pre.dtvdz_mid']:

        if var.split('_')[-1] == 'low':
            lev0, lev1 = 1.5, 4.5
        elif var.split('_')[-1] == 'mid':
            lev0, lev1 = 4.5, 7.5
        else:
            print 'check lev in var',var
            sys.exit() 
        
        if var.split('.')[0]=='pre':
            shead = 'pre.' 
        else:
            shead = ''

        if var in ['deptdz_low','deptdz_mid','pre.deptdz_low','pre.deptdz_mid']:
            varName = 'ept'

        elif var in ['dtvdz_low','dtvdz_mid','pre.dtvdz_low','pre.dtvdz_mid']:
            varName = 'tv'
 
        srcDir = pairbaseDir + '/%s/%04d/%02d/%02d'%(varName, Year,Mon,Day)
        a1var0 = np.load(srcDir + '/%s%s.%04.1fkm.%s.npy'%(shead,varName, lev0,oid))
        a1var1 = np.load(srcDir + '/%s%s.%04.1fkm.%s.npy'%(shead,varName,lev1,oid))
        a1var = (ma.masked_less(a1var1,0) - ma.masked_less(a1var0,0))/(lev1-lev0)


    elif var in ['cape','tcwv','mvimd','pre.cape','skt']:
        if var.split('.')[0]=='pre':
            vardirname = var.split('.')[1]
        else:
            vardirname = var
        a1var = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.00.0km.%s.npy'%(vardirname,Year,Mon,Day,var,oid))

    elif var.split('_')[0] in ['w','r']:
        lev = float(var.split('_')[1])
        varHead = var.split('_')[0]        
        a1var = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.%04.1fkm.%s.npy'%(varHead,Year,Mon,Day,varHead,lev,oid))

 
    else:
        print 'check var',var
        sys.exit()
    return a1var


#***********************8
for season in lseason:
    for var in lvar:
        for prtype in lprtype:
            for dt in ldtenv:
                #****************************
                # Initialize
                #----------------------------
                d1var = {}
                d1frac= {}
                for height in lheight:
                    for region in lregion:
                        for surf in lsurf:
                            dvar [height,region,surf] = np.array([])
                            dfrac[height,region,surf] = np.array([]) 
                #----------------------------
                lDTime = season_lDTime(season)
    
                for DTime in lDTime:
                    print DTime
                    Year,Mon,Day = DTime.timetuple()[:3]
                    convDir = pairbaseDir + '/convfrac/%04d/%02d/%02d'%(Year,Mon,Day)
                    stratDir= pairbaseDir + '/stratfrac/%04d/%02d/%02d'%(Year,Mon,Day)
                    varDir  = pairbaseDir + '/%s/%04d/%02d/%02d'%(var,Year,Mon,Day)
                    surfDir = pairbaseDir + '/%s/%04d/%02d/%02d'%('landSurfaceType',Year,Mon,Day)
                    stopDir = pairbaseDir + '/%s/%04d/%02d/%02d'%('heightStormTop',Year,Mon,Day)
                    latDir  = pairbaseDir + '/%s/%04d/%02d/%02d'%('Latitude',Year,Mon,Day)
                    lonDir  = pairbaseDir + '/%s/%04d/%02d/%02d'%('Longitude',Year,Mon,Day)
                
                    lconvPath = np.sort(glob.glob(convDir + '/*.npy'))
                
                    for convPath in lconvPath:
                        oid = convPath.split('.')[-2]
        
                        if int(oid) in lskipoid:
                            print 'skip',oid
                            continue
         
                        a1varTmp = ret_var(var)
                        a1surfTmp= np.load(surfDir + '/landSurfaceType.00.0km.%s.npy'%(oid))
                        a1stopTmp= np.load(stopDir + '/heightStormTop.00.0km.%s.npy'%(oid))
                        a1latTmp = np.load(stopDir + '/Latitude.00.0km.%s.npy'%(oid))
                        a1lonTmp = np.load(stopDir + '/Longitude.00.0km.%s.npy'%(oid))
        
                        if prtype =='conv':
                            a1fracTmp = np.load(convDir  + '/convfrac.00.0km.%s.npy'%(oid))
                        elif prtype == 'strat':
                            a1fracTmp = np.load(stratDir + '/stratfrac.00.0km.%s.npy'%(oid))
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
                
                        d1flagr = {}
                        d1flagr['t']    = ma.masked_less(np.abs(a1latTmp), 15).mask
                        d1flagr['s']    = ma.masked_inside(np.abs(a1latTmp), 20, 30).mask
                        d1flagr['m']    = ma.masked_greater(np.abs(a1latTmp), 35).mask
                
                        d1flags = {}
                        d1flags['land']= ma.masked_inside(a1surfTmp,0,99).mask
                        d1flags['sea'] = ma.masked_inside(a1surfTmp,100,199).mask
        
                        aflagvar = ma.masked_not_equal(a1var,-9999.).mask
        
                        lkey = [[height,region,surf]
                                for height in lheight
                                for region in lregion
                                for surf   in lsurf
                                ]
                        for (height,region,surf) in lkey:
                            a1flagh = d1flagh[height]
                            a1flagr = d1flagr[region]
                            a1flags = d1flags[surf]
        
                            a1flag = a1flagvar * a1flagh * a1flagr * a1flags * a1flage
                            #*******************
                            # Stack
                            #-------------------
                            d1var[height,region,surf] = np.concatenate([d1var[height,region,surf], a1varTmp[a1flag]])
        
                            d1frac[height,region,surf] = np.concatenate([d1frac[height,region,surf], a1convTmp[a1flag]])
    
            #**********************
            # Save
            #----------------------
            varPath  = pickleDir + '/stack.var.%s.dt%03d.%s.bfile'%(var,dt,season)
            fracPath = pickleDir + '/stack.frac%s.%s.%s.dt%03d.bfile'%(prtype,var,dt,season)
    
            if calcflag == True:
                with open(varPath, 'wb') as f:
                    pickle.dump(dvar, f)
            
                with open(fracPath, 'wb') as f:
                    pickle.dump(dconv, f)

            #*********************
            # Load
            #---------------------
            with open(varPath, 'rb') as f:
                dvar = pickle.load(f)
            
            with open(fracPath, 'rb') as f:
                dconv = pickle.load(f)
     
            #**********************
            #----------------------
            lkey = [[height,region,surf]
                    for height in lheight
                    for region in lregion
                    for surf   in lsurf
                    ]
            for (height,region,surf) in lkey:
                a1var  = d1var [height,region,surf]
                a1frac = d1frac[height,region,surf]
                #-- Histograms --------
                xmin  = np.percentile(a1var,10)
                xmax  = np.percentile(a1var,90)
                xbins  = np.linspace(xmin,xmax,20)
                ybins  = np.linspace(0,1,20)
                
                H,xedges,yedges = np.histogram2d(a1var, a1conv, bins = [xbins,ybins])
                H = H.T
                H = ma.masked_equal(H,0)
                X,Y = np.meshgrid(xbins,ybins)
                
                #-- Binned average -----
                a1convmean,_,_ = stats.binned_statistic(a1var, a1conv, bins = xbins)
                a1xcenter = 0.5*(xbins[:-1]+xbins[1:])
                
                #*** Prep for figure **
                if   var in ['deptdz_low','deptdz_mid','pre.deptdz_low','pre.deptdz_mid']:
                
                    lev = var.split('_')[1]
                    varlongname = r'$d\theta$' +'e/dz %s'%(lev)
                    if var.split('.')[0]=='pre':
                        varlongname = varlongname + ' pre'
                
                    sunit = 'K/km'
                
                elif   var in ['dtvdz_low','dtvdz_mid','pre.dtvdz_low','pre.dtvdz_mid']:
                
                    lev = var.split('_')[1]
                    varlongname = 'Tv/dz %s'%(lev)
                    if var.split('.')[0]=='pre':
                        varlongname = varlongname + ' pre'
                
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
                
                corr = np.corrcoef(a1var, a1conv)[0][1]
                #*** Figure ******
                fig = plt.figure(figsize=(4,4))
                ax  = fig.add_axes([0.17,0.15,0.65,0.65])
                
                im  = ax.pcolormesh(X,Y,H, cmap='jet', norm=matplotlib.colors.LogNorm())
                
                #-- binned mean convective ratio --
                ax.plot(a1xcenter, a1convmean, '-',color='k')
                
                #----------------------------------
                ax.set_xlabel(sunit, fontsize=20)
                ax.set_ylabel('convective ratio', fontsize=20)
                
                #stitle = '%s %s'%(varlongname,surftype)
                stitle = varlongname + ' ' + surftype
                stitle = stitle + '\n' + 'corr=%.2f'%(corr)
                plt.title(stitle, fontsize=19)
                
                cax = fig.add_axes([0.84,0.17,0.02, 0.6])
                cbar=plt.colorbar(im, orientation='vertical', cax=cax)
                cbar.ax.tick_params(labelsize=16)
                
                
                
                figPath = figDir + '/scatter.frac.%s.%s.dt%03d.%s.%s.%s.%s.txt'%(prtype,var,dt,season,height,region,surftype)
                plt.savefig(figPath)
                print figPath
                plt.clf() 
                
                print 'corr=',corr
                #-- write corrcoef ---
                corrDir = figDir + '/txt'
                util.mk_dir(corrDir)
                corrPath= corrDir + '/cc.%s.%s.dt%03d.%s.%s.%s.%s.txt'%(prtype,var,dt,season,height,region,surftype)
                f=open(corrPath,'w');f.write('%.2f'%(corr)); f.close()
    
