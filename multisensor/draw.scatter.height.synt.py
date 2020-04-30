import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as cm
import numpy as np
from numpy import *
import h5py
from bisect import bisect_left
import epcfunc
import sys, os, glob, socket
import myfunc.util as util
import calendar
import pickle
import JPLDB, EPCDB

calcflag  = True
#calcflag  = False
coefflag  = 'nocoef'  # 'nocoef', 'wcoef'
#coefflag  = 'wcoef'  # 'nocoef', 'wcoef'
DB_MAXREC = 10000
DB_MINREC = 1000
nsample   = 1000
lsensor = ['GMI','AMSR2','SSMIS','ATMS','MHS']
#lsensor = ['SSMIS','ATMS','MHS']
#lsensor = ['AMSR2','ATMS','MHS']
#lsensor = ['GMI']

#** Constants ******

#lsurftype = ['all','ocean','vegetation','coast','snow']
lsurftype = ['all','ocean','vegetation','coast']
#lsurftype = ['ocean']
lthwat   = [0.033, 0.05, 0.2]
dsurflabel={ 'all': 'All surface' 
            ,'ocean':'Class1 (Ocean)'
            ,'vegetation':'Classes 3-7 (Vegetation)'
            ,'snow':'Classes 8-11 (Snow)'
            ,'coast':'Class 13 (Coast)'
            }

dlsatid = {'GMI':[1],'AMSR2':[30], 'SSMIS':[16,17,18,19], 'ATMS':[100,101], 'MHS':[201,202,318,319]}

dsatname = {999:'ALLSATE',0:'TRMM',1:'GPM',16:'F16',17:'F17',18:'F18',19:'F19',30:'GCOMW',100:'NPP',101:'NOAA20',201:'METOP-A',202:'METOP-B',318:'NOAA18',319:'NOAA19', 400:'SAPHIR'}
#**************************************
# Functions
#--------------------------------------
def read_orbitlist(Year,Mon):
    listPath = listDir + '/overpass.GPM.%04d.%02d.csv'%(Year,Mon)
    f=open(listPath); lines=f.readlines(); f.close()
    lout=[]
    for line in lines:
        line = map(int,line.strip().split(','))
        lout.append(line)
    return lout

#--------------------------------------

dvnummax = {}
for (sensor,thwat) in [(sensor,thwat)
                for sensor  in lsensor
                for thwat   in lthwat]:

    if sensor=='GMI':
        dbtype='my'
        db = EPCDB.EPCDB()
    else:
        dbtype='jpl'
        db = JPLDB.JPLDB(sensor)

    myhost = socket.gethostname()
    if myhost =='shui':
        if dbtype=='jpl':
            countDir= '/media/disk2/share/PMM/JPLDB/list'
            tankDir = '/tank'
            figDir   = '/home/utsumi/temp/multi'

        elif dbtype=='my':
            countDir= '/work/hk01/utsumi/PMM/EPCDB/list'
            tankDir = '/tank'
            figDir   = '/home/utsumi/temp/multi'

    elif myhost == 'well':
        if dbtype=='jpl':
            countDir= '/media/disk2/share/PMM/JPLDB/list'
            tankDir = '/home/utsumi/mnt/lab_tank'
            figDir   = '/home/utsumi/temp/multi'
        elif dbtype=='my':
            countDir= '/media/disk2/share/PMM/EPCDB/list'
            tankDir = '/home/utsumi/mnt/lab_tank'
            figDir   = '/home/utsumi/temp/multi'

    else:
        print 'check hostname',myhost
        sys.exit()



    expr = 'prf.%s.smp%d'%(sensor,nsample)
    exprhead = expr.split('.')[0]

    retbaseDir = tankDir + '/utsumi/PMM/retsynt/%s'%(expr)
    if dbtype=='jpl':
        dbDir = '/media/disk2/share/PMM/JPLDB/EPC_DB/%s_EPC_DATABASE'%(sensor)
    else:
        dbDir = tankDir + '/utsumi/PMM/EPCDB/samp.%05d.%s.V05A.S1.ABp103-117.01-12'%(DB_MAXREC, sensor)

    #***************************
    # Read paired data
    #---------------------------
    pairdir = tankDir + '/utsumi/PMM/multi/pair-prof-rs/%s.%s'%(exprhead,sensor)
    a1lat     = np.load(pairdir + '/lat.npy')
    a1lon     = np.load(pairdir + '/lon.npy')
    a1lvfracconv = np.load(pairdir + '/vfracconv.npy')
    a1stype   = np.load(pairdir + '/stype.npy')
    #a1t2m     = np.load(pairdir + '/t2m.npy')
    #a1inc     = np.load(pairdir + '/inc.npy')
    #a1precrad = np.load(pairdir + '/precrad.npy')
    #a1precepc = np.load(pairdir + '/precepc.npy')
    #a1precgpr = np.load(pairdir + '/precgpr.npy')
    a2profrad = np.load(pairdir + '/profrad.npy')   # Bottom to top, 500m res
    a2profepc = np.load(pairdir + '/profepc.npy')   # Bottom to top, 500m res
    a2profgpr = np.load(pairdir + '/profgpr.npy')   # Bottom to top, 500m res

    #***************************
    # Storm height
    #---------------------------
    ny,nz = a2profrad.shape
    a2iz   = np.array( range(nz) * ny).reshape(ny,nz)
    a1hrad = ma.masked_where(a2profrad<thwat, a2iz).argmax(axis=1)*0.5 
    a1hepc = ma.masked_where(a2profepc<thwat, a2iz).argmax(axis=1)*0.5 
    a1hgpr = ma.masked_where(a2profgpr<thwat, a2iz).argmax(axis=1)*0.5 


    ##** Read histogram of EPC DB records **
    #if coefflag == 'wcoef':
    #    d1coef = {}
    #    for satid in [999] + lsatid:
    #        if ((dbtype=='my')and(sensor=='GMI')):
    #            countPath = countDir +'/count.epc.csv'
    #        elif dbtype=='jpl':
    #            countPath = countDir +'/count.jpldb.%s.%d.csv'%(sensor, satid)
    #        else:
    #            print 'check dbtype,sensor',dbtype,sensor
    #            sys.exit()

    #        f=open(countPath,'r'); lines=f.readlines(); f.close()
    #        ncol = len(lines[0].strip().split(',')) - 1
    #        d1coef[satid] = np.zeros(29*29*29, float32)
    #        
    #        for line in lines[1:]:
    #            line =  map(int,line.strip().split(','))
    #            epcid= line[0]
    #            nall = line[1]
    #            if nall < nsample:
    #                coef = 1.0
    #            else:
    #                coef = float(nall)/nsample
    #            d1coef[satid][epcid] = coef


 
    #elif coefflag =='nocoef':
    #    d1coef[satid] = np.ones(29*29*29, float32)
    #else:
    #    print 'check coefflag', coefflag
    #    sys.exit()
    
    ##---------------------------
    a1masklat = ma.masked_outside(a1lat, -65,65).mask

    for surftype in lsurftype:
        #-- Screen by surface types -------------------
        if   surftype=='all':
            a1masksurf = np.array([False])
        elif surftype=='ocean':
            a1masksurf = ma.masked_not_equal(a1stype,1).mask
        elif surftype=='vegetation':
            a1masksurf = ma.masked_outside(a1stype,3,7).mask
        elif surftype=='snow':
            a1masksurf = ma.masked_outside(a1stype,8,11).mask
        elif surftype=='coast':
            a1masksurf = ma.masked_not_equal(a1stype,13).mask
        else:
            print '\n'+'check surftype',surftype
            sys.exit() 
    
        a1mask  = a1masklat + a1masksurf
 
        if a1mask.sum() == 0:
            a1radTmp = a1hrad
            a1epcTmp = a1hepc
            a1gprTmp = a1hgpr

        else:
            a1radTmp = ma.masked_where(a1mask, a1hrad).compressed()
            a1epcTmp = ma.masked_where(a1mask, a1hepc).compressed()
            a1gprTmp = ma.masked_where(a1mask, a1hgpr).compressed()

        if a1radTmp.shape[0]==0:
            print 'skip',sensor, rettype, surftype
            continue


        #-- Histograms log-scale -------
        dvnummax = {}
        for rettype in ['EPC','GPROF']:
            a1obs = a1radTmp
            if rettype=='EPC':
                a1ret = a1epcTmp
            elif rettype=='GPROF':
                a1ret = a1gprTmp

            a1bin = np.arange(0,10+0.1,0.5)
            axmin,axmax = a1bin.min(), a1bin.max()
            H,xedges,yedges = np.histogram2d(a1obs, a1ret, bins = [a1bin, a1bin])
            H = H.T
            #H = H * d1coef[satid][idx_db]
            X,Y = np.meshgrid(xedges, yedges)
            if rettype=='EPC':
                vnummax = H.max()
            #*******************************
            # Figure
            #-------------------------------
            fig = plt.figure(figsize=[6,6])
            ax  = fig.add_axes([0.2,0.13,0.61,0.61])
            im  = ax.pcolormesh(X,Y,H, norm=matplotlib.colors.LogNorm(), cmap='jet', vmin=1, vmax=vnummax)
    
            if H.max() ==0:
                print 'No precipitating case',satid, surftype
                print 'Skip'
                continue 
            
            #-- plot 1:1 line
            ax.plot(array([axmin,axmax]),array([axmin,axmax]),'-',color='k',linewidth=0.5)
            #-- axis labels ---
            lticks = a1bin[::2].astype('int32')
            lticklabels = a1bin[::2].astype('int32')

            ax.set_xticks(lticks)
            ax.set_xticklabels(lticklabels, fontsize=16)
            ax.set_yticks(lticks)
            ax.set_yticklabels(lticklabels, fontsize=16)

            ax.set_xlabel('DPR/NScmb height [km]', fontsize=20)

            ax.set_ylabel('%s height [km]'%(rettype), fontsize=20)
            ax.set_ylim([axmin,axmax])
            ax.set_xlim([axmin,axmax])

            stitle =  '%s (%s) %.3fg/m3\n%s'%(rettype, sensor, thwat, dsurflabel[surftype])
            plt.title(stitle, fontsize=22)

            cax = fig.add_axes([0.82,0.15,0.02, 0.6])
            cbar=plt.colorbar(im, orientation='vertical', cax=cax)
            cbar.ax.tick_params(labelsize=16)

            figPath= figDir + '/scatter.height.%s.%s.%s.%s.th-%.3f.png'%(expr, coefflag, rettype,surftype,thwat)
            util.mk_dir(figDir)

            plt.savefig(figPath)
            print figPath
