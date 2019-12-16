import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import numpy as np
from datetime import datetime, timedelta
import os, sys, socket
import glob
import myfunc.util as util
import pickle
import calendar

calcflag = True
#calcflag = False

ldtenv = [-48,-36,-24,-12,-6,-3,-1,0,1,3,6,12,24,36,48]
#ldtenv = [-6,-1,0,6]
#ldtenv = np.arange(-12,12+1,2)
iYM   = [2017,7]
eYM   = [2017,7]
lYM   = util.ret_lYM(iYM,eYM)

dprver = 'V06'
dprverfull='V06A'
#lvar =['cape','r2m','q2m']
lvar =['skt','tcwv','mvimd','cape','r2m','q2m']
#lvar =['cape']
#lvar = ['r_75']
#lvar = ['dtvdz_low']
#lvar = ['deptdz_low']

lprtype= ['conv','strat']
lheight= ['high','low'] 
lregion= ['t','s','m']
lsurf  = ['land','sea']

#lprtype= ['conv']
#lheight= ['high'] 
#lregion= ['l']
#lsurf  = ['sea']



myhost = socket.gethostname()
if myhost =='shui':
    pairbaseDir = '/tank/utsumi/env/pair'
    pickleDir   = '/tank/utsumi/env/pickle'
elif myhost == 'well':
    pairbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/env/pair'
    meanbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/env/eramean'
    pickleDir = '/home/utsumi/mnt/lab_tank/utsumi/env/pickle'
else:
    print 'check myhost'
    sys.exit()

miss_out = -9999.

#****************************
# Functions
#----------------------------
def ret_var(var):
    varName,lev = var.split('_')
    if lev in ['low','mid']: 
        if lev=='low':
            lev0 = 1.5
            lev1 = 4.5
        elif lev=='mid':
            lev0 = 4.5
            lev1 = 7.5

        varNameIn = varName[1:-2]
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


#****************************
# Initialize
#----------------------------
dsum = {}
dnum = {}
dclm = {}
dstd = {}

for var in lvar:
    for prtype in lprtype:
        for height in lheight: 
            for region in lregion:
                for surf in lsurf:
                    for dt in ldtenv:
                        dsum[var,prtype,height,region,surf,dt] = 0.0
                        dnum[var,prtype,height,region,surf,dt] = 0
                        dclm[var,prtype,height,region,surf,dt] = 0
                        dstd[var,prtype,height,region,surf,dt] = 0

#****************************

#****************************
# Calculation
#----------------------------
for (Year,Mon) in lYM:
    eDay   = calendar.monthrange(Year,Mon)[1]
    iDTime = datetime(Year,Mon,1)
    #eDTime = datetime(Year,Mon,eDay)
    eDTime = datetime(Year,Mon,eDay) - timedelta(days=1)
    #eDTime = datetime(Year,Mon,5)
    lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
    
    for DTime in lDTime:
        if calcflag != True: continue
    
        Year,Mon,Day = DTime.timetuple()[:3]
        latDir = pairbaseDir + '/Latitude/%04d/%02d/%02d'%(Year,Mon,Day)
        llatPath = np.sort(glob.glob(latDir + '/Latitude.00.0km.??????.npy'))
        
        for latPath in llatPath:
            oid = int(latPath.split('.')[-2])
            print DTime,oid
        
            #a1conv = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.00.0km.%06d.npy'%('convfrac',Year,Mon,Day,'convfrac',oid))
    
            #a1strat= np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.00.0km.%06d.npy'%('stratfrac',Year,Mon,Day,'stratfrac',oid))

            a1conv = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.00.0km.%06d.npy'%('cvfrac',Year,Mon,Day,'cvfrac',oid))
    
            a1strat= np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.00.0km.%06d.npy'%('svfrac',Year,Mon,Day,'svfrac',oid))

           
            a1lat  = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.00.0km.%06d.npy'%('Latitude',Year,Mon,Day,'Latitude',oid))
            a1lon  = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.00.0km.%06d.npy'%('Longitude',Year,Mon,Day,'Longitude',oid))
            a1stop = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.00.0km.%06d.npy'%('heightStormTop',Year,Mon,Day,'heightStormTop',oid))
            a1surf = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.00.0km.%06d.npy'%('landSurfaceType',Year,Mon,Day,'landSurfaceType',oid))
    
            #****************
            # Flags
            #----------------
            d1flagp = {}
            #d1flagp['conv']  = ma.masked_greater(a1conv,0.8).mask
            #d1flagp['strat'] = ma.masked_greater(a1strat,0.8).mask
            d1flagp['conv']  = ma.masked_greater(a1conv,0.5).mask
            d1flagp['strat'] = ma.masked_greater(a1strat,0.5).mask

    
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
    
            for var in lvar:
                for dt in ldtenv:
                    if len(var.split('_'))==2:
                        a1env = ret_var(var)
    
                    else:
                        a1env = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/%s.%03dh.00.0km.%06d.npy'%(var,Year,Mon,Day,var,dt,oid))
                        a1ave = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/ave.%03dh.00.0km.%06d.npy'%(var,Year,Mon,Day,0,oid))
                        a1std = np.load(pairbaseDir + '/%s/%04d/%02d/%02d/std.%03dh.00.0km.%06d.npy'%(var,Year,Mon,Day,0,oid))

                    a1ave = ma.masked_equal(a1ave,-9999)
                    a1std = ma.masked_equal(a1std,-9999)

                    a1env = ma.masked_where(a1std==0, a1env - a1ave)/a1std
                    #a1env = ma.masked_where(a1ave==0, a1env - a1ave)
                    a1env = ma.masked_invalid(a1env)
                    a1env = ma.masked_equal(a1env,-9999)
                    
                    #-------------------
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
                                        s   = 0.0
                                        clm = 0.0
                                        std = 0.0 

                                    else:
                                        s   = a1env[a1flag].sum()
                                        clm = a1ave[a1flag].sum()  
                                        std = a1std[a1flag].sum()  
 
                                    dsum[var,prtype,height,region,surf,dt] = dsum[var,prtype,height,region,surf,dt] + s
                                    dnum[var,prtype,height,region,surf,dt] = dnum[var,prtype,height,region,surf,dt] + n


                                    dclm[var,prtype,height,region,surf,dt] = dclm[var,prtype,height,region,surf,dt] + clm
                                    dstd[var,prtype,height,region,surf,dt] = dstd[var,prtype,height,region,surf,dt] + std

                #print var,prtype,height,region,surf,dt,dsum[var,prtype,height,region,surf,dt],dnum[var,prtype,height,region,surf,dt]
                                
#*********************
# Save
#--------------------- 
sumPath = pickleDir + '/time.intersect.sum.bfile'
numPath = pickleDir + '/time.intersect.num.bfile'
if calcflag == True:
    with open(sumPath, 'wb') as f:
        pickle.dump(dsum, f)
    
    with open(numPath, 'wb') as f:
        pickle.dump(dnum, f)

#*********************
# Load
#--------------------- 
with open(sumPath, 'rb') as f:
    dsum = pickle.load(f)

with open(numPath, 'rb') as f:
    dnum = pickle.load(f)

#*********************
# Figure
#--------------------- 

for var in lvar:
    for height in lheight:
        for surf in lsurf:

            #-- vmin, vmax -----
            vmin, vmax = None,None
            #if var == 'cape':
            #    vmin,vmax = 0, 1100
            #elif var in ['r2m']:
            #    vmin,vmax = 0.6,1.0
            #elif var in ['q2m']:
            #    vmin,vmax = 0.004,0.020
            #else:
            #    vmin,vmax = None,None

            #-- change unit -----
            coef = 1.0
            #if var == 'mvimd':
            #    coef = 60*60. #kg m**-2 s**-1 --> kg m**-2 h**-1
            #else:
            #    coef = 1.0




            a1x = ldtenv
            d1y = {}
            d1n = {}
            d1clm = {}
            d1std = {}
            for prtype in lprtype:
                for region in lregion:
                    a1s = np.array([dsum[var,prtype,height,region,surf,dt] for dt in ldtenv])
                    a1n = np.array([dnum[var,prtype,height,region,surf,dt] for dt in ldtenv])

                    a1clm=np.array([dclm[var,prtype,height,region,surf,dt] for dt in ldtenv])
                    a1std=np.array([dstd[var,prtype,height,region,surf,dt] for dt in ldtenv])

                    d1y[prtype,region] = ma.masked_invalid(ma.masked_where(a1n==0, a1s) / a1n) * coef
           
                    d1n[prtype,region] = a1n 
                    d1clm[prtype,region] = ma.masked_where(a1n==0, a1clm)/a1n
                    d1std[prtype,region] = ma.masked_where(a1n==0, a1std)/a1n

            #-- Figure variables ---
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
            
            #a1xticks = [-48,-36,-24,-12,-6,-3,0,3,6,12,24,36,48]
            a1xticks = ldtenv
            ax.set_xticks(a1xticks)
            ax.set_xticklabels(a1xticks)
            
            ax.grid(which = "major", axis = "x", color = "gray", alpha = 0.8, linestyle = "--", linewidth = 1)
            
            plt.title('%s %s-cloud.(%s) '%(var, height, surf))
            plt.legend()    
            figPath = '/media/disk2/home_temp/env/time.itst.%s.%s-stop.%s.png'%(var,height,surf)
            plt.savefig(figPath)
            print figPath
            plt.clf()
           
            ##******************** 
            ##-- Figure Clim ---
            #fig = plt.figure(figsize = (6,3))
            #ax  = fig.add_axes([0.1,0.1,0.8,0.8])
            #
            #ax.plot(a1x, d1clm['conv','t'], '-',  color='r', label='conv/trop')
            #ax.plot(a1x, d1clm['strat' ,'t'], '--', color='r', label='strat/trop')
            #ax.plot(a1x, d1clm['conv','s'], '-',  color='orange', label='conv/subtrop')
            #ax.plot(a1x, d1clm['strat' ,'s'], '--', color='orange', label='strat/subtrop')
            #ax.plot(a1x, d1clm['conv','m'], '-',  color='blue', label='conv/mid lat')
            #ax.plot(a1x, d1clm['strat' ,'m'], '--', color='blue', label='strat/mid lat')
           
            #plt.ylim([vmin,vmax]) 
            ##plt.ylim([vmin,100000]) 
            ##- x=0 ---
            #plt.axvline(x=0,linestyle='--',color='gray', linewidth=2)
            #
            ##a1xticks = [-48,-36,-24,-12,-6,-3,0,3,6,12,24,36,48]
            #a1xticks = ldtenv
            #ax.set_xticks(a1xticks)
            #ax.set_xticklabels(a1xticks)
            #
            #ax.grid(which = "major", axis = "x", color = "gray", alpha = 0.8, linestyle = "--", linewidth = 1)
            #
            #plt.title('Clim %s %s-cloud.(%s) '%(var, height, surf))
            #plt.legend()    
            #figPath = '/media/disk2/home_temp/env/time.itst.clm-%s.%s-stop.%s.png'%(var,height,surf)
            #plt.savefig(figPath)
            #print figPath
            #plt.clf()
            #

            ##******************** 
            ##-- Figure Std ---
            #fig = plt.figure(figsize = (6,3))
            #ax  = fig.add_axes([0.1,0.1,0.8,0.8])
            #
            #ax.plot(a1x, d1std['conv','t'], '-',  color='r', label='conv/trop')
            #ax.plot(a1x, d1std['strat' ,'t'], '--', color='r', label='strat/trop')
            #ax.plot(a1x, d1std['conv','s'], '-',  color='orange', label='conv/subtrop')
            #ax.plot(a1x, d1std['strat' ,'s'], '--', color='orange', label='strat/subtrop')
            #ax.plot(a1x, d1std['conv','m'], '-',  color='blue', label='conv/mid lat')
            #ax.plot(a1x, d1std['strat' ,'m'], '--', color='blue', label='strat/mid lat')
           
            #plt.ylim([vmin,vmax]) 
            ##plt.ylim([vmin,100000]) 
            ##- x=0 ---
            #plt.axvline(x=0,linestyle='--',color='gray', linewidth=2)
            #
            ##a1xticks = [-48,-36,-24,-12,-6,-3,0,3,6,12,24,36,48]
            #a1xticks = ldtenv
            #ax.set_xticks(a1xticks)
            #ax.set_xticklabels(a1xticks)
            #
            #ax.grid(which = "major", axis = "x", color = "gray", alpha = 0.8, linestyle = "--", linewidth = 1)
            #
            #plt.title('std %s %s-cloud.(%s) '%(var, height, surf))
            #plt.legend()    
            #figPath = '/media/disk2/home_temp/env/time.itst.std-%s.%s-stop.%s.png'%(var,height,surf)
            #plt.savefig(figPath)
            #print figPath
            #plt.clf()
            




