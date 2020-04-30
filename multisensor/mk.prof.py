# %%
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import ListedColormap
#%matplotlib inline

from numpy import *
import myfunc.util as util
from datetime import datetime, timedelta
import glob, os, sys
import numpy as np
import calendar
import random

tankDir = '/home/utsumi/mnt/lab_tank'
figDir  = '/home/utsumi/temp/multi'
util.mk_dir(figDir)

miss   = -9999.
thpr = 0.5  # mm/h
expr = 'prf'

lsensor = ['AMSR2','SSMIS','ATMS','MHS']
#lsensor = ['AMSR2']
ltrange  = [(0,7), (20,20+7), (26,26+7)]
lstype = ['sea','veg']
lregion= ['tro','mid']
lptype = ['conv','stra']

#lsensor  = ['AMSR2']
#ltrange  = [(26,26+7)]
#lstype = ['sea']
#lregion= ['tro']
#lptype = ['stra']


#***************************************
lkey= [(sensor,trange,stype,region,ptype)
    for sensor  in lsensor
    for trange in ltrange
    for stype   in lstype
    for region  in lregion
    for ptype   in lptype
    ]

for (sensor,trange,stype,region,ptype) in lkey:
    if region=='tro':
        if trange !=(26,26+7):
            print 'Skip',region,'T2m=',trange
            continue
    elif region=='mid':
        if trange ==(26,26+33):
            print 'Skip',region,'T2m=',trange
            continue

    print sensor
    #*********************************
    # Make mean profiles
    #---------------------------------
    srcdir = tankDir + '/utsumi/PMM/multi/pair-prof-rs/%s.%s'%(expr,sensor)
    a1lat    = np.load(srcdir + '/lat.npy')
    a1lon    = np.load(srcdir + '/lon.npy')
    a1precrad= np.load(srcdir + '/precrad.npy')
    a1precepc= np.load(srcdir + '/precepc.npy')
    a1precgpr= np.load(srcdir + '/precgpr.npy')
    a2profrad= np.load(srcdir + '/profrad.npy') # Bottom to top
    a2profepc= np.load(srcdir + '/profepc.npy') # Bottom to top
    a2profgpr= np.load(srcdir + '/profgpr.npy') # Bottom to top
    a1stype  = np.load(srcdir + '/stype.npy')
    a1t2m    = np.load(srcdir + '/t2m.npy')
    a1inc    = np.load(srcdir + '/inc.npy')
    a1vfracconv = np.load(srcdir + '/vfracconv.npy')

    #-- Region ---------------
    a1flagtro  = ma.masked_inside(a1lat,-15,15).mask
    a1flagmid  = ma.masked_inside(np.abs(a1lat),35,60).mask

    d1flagreg = {}
    d1flagreg['tro'] = a1flagtro
    d1flagreg['mid'] = a1flagmid

    #-- Precip type ----------
    a1flagconv = ma.masked_greater_equal(a1vfracconv,0.6).mask
    a1flagstra = ma.masked_less(a1vfracconv,0.4).mask

    d1flagptype= {}
    d1flagptype['conv']=a1flagconv
    d1flagptype['stra']=a1flagstra

    #-- Surface type --
    a1flagsea  = ma.masked_equal(a1stype,1).mask
    a1flagveg  = ma.masked_inside(a1stype,3,7).mask
    a1flagsnow = ma.masked_inside(a1stype,8,11).mask
    a1flagcoast= ma.masked_equal(a1stype,13).mask
    a1flagland = a1flagveg + a1flagsnow

    d1flagstype= {}
    d1flagstype['sea'] = a1flagsea
    d1flagstype['veg'] = a1flagveg
    d1flagstype['snow'] = a1flagsnow
    d1flagstype['coast'] = a1flagcoast
    d1flagstype['land'] = a1flagland

    d1flagt = {}
    t0,t1 = trange
    d1flagt[trange] = ma.masked_inside(a1t2m, 273.15+t0, 273.15+t1).mask

    a1flag = d1flagreg[region] * d1flagptype[ptype] * d1flagstype[stype] * d1flagt[trange]

    print a1flag.shape[0],d1flagreg[region].sum() , d1flagptype[ptype].sum() , d1flagstype[stype].sum() , d1flagt[trange].sum()


    #********************************************** 
    # Figure
    #********************************************** 
    stamp ='%s-s-%s.p-%s.r-%s.t-%d-%d'%(sensor, stype,ptype,region,t0,t1)

    a2radtmp = ma.masked_invalid(ma.masked_less(a2profrad[a1flag],0))
    a2epctmp = ma.masked_invalid(ma.masked_less(a2profepc[a1flag],0))
    a2gprtmp = ma.masked_invalid(ma.masked_less(a2profgpr[a1flag],0))


    a1profrad = a2radtmp.mean(axis=0)
    a1profepc = a2epctmp.mean(axis=0)
    a1profgpr = a2gprtmp.mean(axis=0)

    a1numrad = a2radtmp.count(axis=0)
    a1numepc = a2epctmp.count(axis=0)
    a1numgpr = a2gprtmp.count(axis=0)


    #-- mask with height --
    a1hgt = np.arange(0.5,10+0.01, 0.5)
    a1profrad[a1hgt<=1] = np.nan
    a1profepc[a1hgt<=1] = np.nan
    a1profgpr[a1hgt<=1] = np.nan

    lvar = ['ave','num']
    for var in lvar:
        if var=='ave':
            a1varrad = a1profrad
            a1varepc = a1profepc
            a1vargpr = a1profgpr
        elif var=='num':
            a1varrad = a1numrad
            a1varepc = a1numepc
            a1vargpr = a1numgpr

        if a1varrad.shape[0]==0:
            print 'No profile',key
            continue
            #sys.exit()

        fig = plt.figure(figsize=(2.5,3.2))
        #fig = plt.figure(figsize=(3,4))
        ax  = fig.add_axes([0.2,0.15,0.65,0.7])

        a1y = a1hgt

        ax.plot( a1varrad, a1y, '-',  c='k', linewidth=2, label='CMB', clip_on=False) 
        ax.plot( a1varepc, a1y, '-',  c='k', linewidth=1, label='EPC', clip_on=False)
        ax.plot( a1vargpr, a1y, '--', c='k', linewidth=1.3, label='GPROF', clip_on=False)

        xmax = {'ave':0.65, 'num':None}[var]
        ax.set_ylim([0.9,10.3])
        ax.set_xlim([0,xmax])
        plt.xlabel('(g/m3)',fontsize=12)
        plt.ylabel('(km)',fontsize=12)
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)


        stitle = '%s %s %s %s %s'%(sensor,region, stype, ptype, var)
        stitle = stitle +'\n' +'T2m=%d-%d (deg.C)'%(t0,t1)
        plt.title(stitle, fontsize=11)
        #plt.legend()   
        plt.show()
        ##---------------------------------
        figPath = figDir + '/prof.%s.%s.png'%(stamp,var)
        plt.savefig(figPath)
        print figPath



    # %%
