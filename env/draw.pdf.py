import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import numpy as np
from datetime import datetime, timedelta
import os, sys, socket
import glob, calendar
import myfunc.util as util
import pickle

calcflag = True
#calcflag = False
lseason = [7]
dprver = 'V06'
dprverfull='V06A'
lvar = ['cape']
#lvar =['cape','r2m','q2m']
#lvar = ['tcwv','mvimd']
lprtype= ['conv','strat']
lheight= ['high','low'] 
lregion= ['t','s','m']
lsurf  = ['land','sea']
ldtenv = [0]

myhost = socket.gethostname()
if myhost =='shui':
    pairbaseDir = '/tank/utsumi/env/pair'
    pickleDir   = '/tank/utsumi/env/pickle'
elif myhost == 'well':
    pairbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/env/pair'
    pickleDir = '/home/utsumi/mnt/lab_tank/utsumi/env/pickle'
else:
    print 'check myhost'
    sys.exit()

miss_out = -9999.
#****************************
# Functions
#----------------------------
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

#****************************
# Calc
#----------------------------
for season in lseason:
    for var in lvar:
        for dt in ldtenv:
            dsum = {}
            dnum = {}
            for prtype in lprtype:
                for height in lheight: 
                    for region in lregion:
                        for surf in lsurf:
                            for dt in ldtenv:
                                dsum[var,prtype,height,region,surf,dt] = 0.0
                                dnum[var,prtype,height,region,surf,dt] = 0
    
            #****************************
            # Initialize
            #----------------------------
            lDTime = season_lDTime(season)
            for DTime in lDTime:
                if calcflag != True: continue
            
                Year,Mon,Day = DTime.timetuple()[:3]
                latDir = pairbaseDir + '/Latitude/%04d/%02d/%02d'%(Year,Mon,Day)
                llatPath = np.sort(glob.glob(latDir + '/Latitude.00.0km.??????.npy'))
                
                for latPath in llatPath:
                    oid = int(latPath.split('.')[-2])
                    print DTime,oid
                
                    a1conv = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.00.0km.%06d.npy'%('convfrac',Year,Mon,Day,'convfrac',oid))
            
                    a1strat= np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.00.0km.%06d.npy'%('stratfrac',Year,Mon,Day,'stratfrac',oid))
                   
                    a1lat  = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.00.0km.%06d.npy'%('Latitude',Year,Mon,Day,'Latitude',oid))
                    a1lon  = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.00.0km.%06d.npy'%('Longitude',Year,Mon,Day,'Longitude',oid))
                    a1stop = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.00.0km.%06d.npy'%('heightStormTop',Year,Mon,Day,'heightStormTop',oid))
                    a1surf = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.00.0km.%06d.npy'%('landSurfaceType',Year,Mon,Day,'landSurfaceType',oid))
            
                    #****************
                    # Flags
                    #----------------
                    d1flagp = {}
                    d1flagp['conv']  = ma.masked_greater(a1conv,0.8).mask
                    d1flagp['strat'] = ma.masked_greater(a1strat,0.8).mask
            
                    d1flagh = {}
                    d1flagh['high'] = ma.masked_greater(a1stop,6000).mask
                    d1flagh['mid' ] = ma.masked_inside(a1stop,4000,6000).mask
                    d1flagh['low' ] = ma.masked_less(a1stop,4000).mask
            
                    d1flagr = {}
                    d1flagr['t']    = ma.masked_less(np.abs(a1lat), 15).mask
            
                    d1flagr['s']    = ma.masked_inside(np.abs(a1lat), 20, 30).mask
                    d1flagr['m']    = ma.masked_greater(np.abs(a1lat), 35).mask
            
                    d1flags = {}
                    d1flags['land'] = ma.masked_inside(a1surf,0,99).mask
                    d1flags['sea'] = ma.masked_inside(a1surf,100,199).mask
                    d1flags['coast'] = ma.masked_inside(a1surf,200,299).mask
            
                    a1env = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.%03dh.00.0km.%06d.npy'%(var,Year,Mon,Day,var,dt,oid))
                    a1env = ma.masked_invalid(a1env)
                    a1env = ma.masked_equal(a1env,-9999)
                    a1flage = ~(a1env.mask)
                    for prtype in lprtype:
                        for height in lheight:
                            for region in lregion:
                                for surf in lsurf:
                                    a1flagp = d1flagp[prtype]
                                    a1flagh = d1flagh[height]
                                    a1flagr = d1flagr[region]
                                    a1flags = d1flags[surf]
            
                                    a1flag = a1flagp * a1flagh * a1flagr * a1flags * a1flage
                
                                    #* Make sum and num **
                                    n = a1flag.sum()
                                    if n==0:
                                        s=0.0
                                    else:
                                        s=a1env[a1flag].sum()
            
                                    dsum[var,prtype,height,region,surf,dt] = dsum[var,prtype,height,region,surf,dt] + s
                                    dnum[var,prtype,height,region,surf,dt] = dnum[var,prtype,height,region,surf,dt] + n
            
            ##*********************
            ## Save
            ##--------------------- 
            #sumPath = pickleDir + '/time.intersect.sum.bfile'
            #numPath = pickleDir + '/time.intersect.num.bfile'
            #if calcflag == True:
            #    with open(sumPath, 'wb') as f:
            #        pickle.dump(dsum, f)
            #    
            #    with open(numPath, 'wb') as f:
            #        pickle.dump(dnum, f)
            #
            ##*********************
            ## Load
            ##--------------------- 
            #with open(sumPath, 'rb') as f:
            #    dsum = pickle.load(f)
            #
            #with open(numPath, 'rb') as f:
            #    dnum = pickle.load(f)
            
            #*********************
            # Figure
            #--------------------- 
            
            for var in lvar:
                for height in lheight:
                    for surf in lsurf:
            
                        #-- vmin, vmax -----
                        if var == 'cape':
                            vmin,vmax = 0, 1100
                        elif var in ['r2m']:
                            vmin,vmax = 0.6,1.0
                        elif var in ['q2m']:
                            vmin,vmax = 0.004,0.020
                        else:
                            vmin,vmax = None,None
            
                        a1x = ldtenv
                        d1y = {}
                        for prtype in lprtype:
                            for region in lregion:
                                a1s = np.array([dsum[var,prtype,height,region,surf,dt] for dt in ldtenv])
                                a1n = np.array([dnum[var,prtype,height,region,surf,dt] for dt in ldtenv])
                                d1y[prtype,region] = ma.masked_invalid(ma.masked_where(a1n==0, a1s) / a1n)
                        
                        
                        fig = plt.figure(figsize = (6,3))
                        ax  = fig.add_axes([0.1,0.1,0.8,0.8])
                        
                        ax.plot(a1x, d1y['conv','t'], '-',  color='r', label='conv/trop')
                        ax.plot(a1x, d1y['strat' ,'t'], '--', color='r', label='strat/trop')
                        ax.plot(a1x, d1y['conv','s'], '-',  color='orange', label='conv/subtrop')
                        ax.plot(a1x, d1y['strat' ,'s'], '--', color='orange', label='strat/subtrop')
                        ax.plot(a1x, d1y['conv','m'], '-',  color='blue', label='conv/mid lat')
                        ax.plot(a1x, d1y['strat' ,'m'], '--', color='blue', label='strat/mid lat')
                       
                        plt.ylim([vmin,vmax]) 
                        #- x=0 ---
                        plt.axvline(x=0,linestyle='--',color='gray', linewidth=2)
                        
                        a1xticks = [-48,-36,-24,-12,-6,-3,0,3,6,12,24,36,48]
                        #a1xticks = ldtenv
                        ax.set_xticks(a1xticks)
                        ax.set_xticklabels(a1xticks)
                        
                        ax.grid(which = "major", axis = "x", color = "gray", alpha = 0.8, linestyle = "--", linewidth = 1)
                        
                        plt.title('%s %s-cloud.(%s)'%(var, height, surf))
                        plt.legend()    
                        figPath = '/media/disk2/home_temp/env/time.itst.%s.%s-stop.%s.png'%(var,height,surf)
                        plt.savefig(figPath)
                        print figPath
                        plt.clf()
            

