import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#import matplotlib.cm as cm
import numpy as np
from numpy import *
import h5py
from bisect import bisect_left
import epcfunc
import sys, os, glob, socket
import myfunc.util as util
import calendar
import pickle

#calcflag  = True
calcflag  = False
#coefflag  = 'nocoef'  # 'nocoef', 'wcoef'
coefflag  = 'wcoef'  # 'nocoef', 'wcoef'
DB_MAXREC = 10000
DB_MINREC = 1000
nsample   = 1000
dbtype = 'my'
#lrettype = ['NS','MS','NScmb','MScmb']
#lrettype = ['NS','MS','NScmb','MScmb','GPROF']
#lrettype = ['GPROF']  
#lrettype = ['MS','NS','NScmb','MScmb']
lrettype = ['NScmb','GPROF']
#prmin = 0.1
#prmin = 0.01
prmin = 0.0
expr = 'org.smp%d'%(nsample)
lidx_db = range(29*29*29)[1:]
#lidx_db = range(29*29*4)
#lidx_db = range(29*29*4,29*29*8)
#lidx_db = range(29*29*8,29*29*12)
#lidx_db = range(29*29*12,29*29*16)
#lidx_db = range(29*29*16,29*29*20)
#lidx_db = range(29*29*20,29*29*24)
#lidx_db = range(29*29*24,29*29*29)
#** Constants ******
myhost = socket.gethostname()
if myhost =='shui':
    if dbtype=='JPL':
        dbDir   = '/work/hk01/utsumi/JPLDB/EPC_DB/GMI_EPC_DATABASE_TEST29'
        countDir= '/work/hk01/utsumi/PMM/EPCDB/list'

    elif dbtype=='my':
        dbDir   = '/work/hk01/utsumi/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
        countDir= '/work/hk01/utsumi/PMM/EPCDB/list'
        retbaseDir = '/tank/utsumi/PMM/retsynt/%s'%(expr)

elif myhost == 'well':
    if dbtype=='JPL':
        dbDir   = '/work/hk01/utsumi/JPLDB/EPC_DB/GMI_EPC_DATABASE_TEST29'
        countDir= '/work/hk01/utsumi/PMM/EPCDB/list'
        figDir   = '/home/utsumi/temp/ret'
    elif dbtype=='my':
        dbDir   = '/media/disk2/share/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
        countDir= '/media/disk2/share/PMM/EPCDB/list'
        retbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retsynt/%s'%(expr)
        figDir   = '/home/utsumi/temp/ret'

else:
    print 'check hostname',myhost
    sys.exit()


lsurftype = ['ocean','vegetation','coast','snow']
#lsurftype = ['ocean']

dsurflabel={ 'ocean':'Class1 (Ocean)'
            ,'vegetation':'Classes 3-7 (Vegetation)'
            ,'snow':'Classes 8-11 (Snow)'
            ,'coast':'Class 13 (Coast)'
            }

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

#** Read histogram of EPC DB records **
if coefflag == 'wcoef':
    countPath = countDir +'/count.epc.csv'
    f=open(countPath,'r'); lines=f.readlines(); f.close()
    ncol = len(lines[0].strip().split(',')) - 1
    a1nall = np.zeros(29*29*29, int32)
    a1coef = np.zeros(29*29*29, float32)
    
    for line in lines[1:]:
        line =  map(int,line.strip().split(','))
        epcid= line[0]
        nall = line[1]
        a1nall[epcid] = nall
        if nall < nsample:
            coef = 1.0
        else:
            coef = float(nall)/nsample
        a1coef[epcid] = coef

elif coefflag =='nocoef':
    a1coef = np.ones(29*29*29, float32)
else:
    print 'check coefflag', coefflag
    sys.exit()

#--------------------------------------
dvnummax = {}
for rettype in lrettype:

    #***************************
    # Initialize histogram
    #---------------------------
    logvmin, logvmax = -1.2, 2
    bins   = np.arange(-2,logvmax+0.001,0.025)
    nbins  = len(bins)
    d2freq = {}
    for surftype in lsurftype:
        d2freq[surftype] = np.zeros([nbins-1,nbins-1],float64)  # x:obs  y:ret

    #---------------------------
    for idx_db in lidx_db:
        if calcflag ==False:
            continue

        if rettype != 'GPROF':
            print 'idx_db=',idx_db
            obsPath = retbaseDir + '/%05d/nsurf%s.obs.%05d.npy'%(idx_db,rettype,idx_db)
            if not os.path.exists(obsPath):
                print 'No file'
                print obsPath
                continue

            a1irec= np.load(retbaseDir + '/%05d/%s.obs.%05d.npy'%(idx_db,'irec',idx_db)).astype(int32)
            a1obs = np.load(retbaseDir + '/%05d/nsurf%s.obs.%05d.npy'%(idx_db,rettype,idx_db))
            a1ret = np.load(retbaseDir + '/%05d/nsurf%s.est.%05d.npy'%(idx_db,rettype,idx_db))

    
        elif rettype =='GPROF':
            obsPath = retbaseDir + '/%05d/nsurf%s.obs.%05d.npy'%(idx_db,'NScmb',idx_db)
            if not os.path.exists(obsPath):
                print 'No file'
                print obsPath
                continue
        
            a1irec= np.load(retbaseDir + '/%05d/%s.obs.%05d.npy'%(idx_db,'irec',idx_db)).astype(int32)
            a1obs = np.load(retbaseDir + '/%05d/nsurf%s.obs.%05d.npy'%(idx_db,'NScmb',idx_db))

            a1ret = np.load(dbDir + '/surfacePrecipitation/%05d.npy'%(idx_db))[a1irec]
            #print 'No GPROF DB'
            #print dbDir + '/surfacePrecipitation/%05d.npy'%(idx_db)
            #    continue

        else:
            print 'check rettype',rettype
            sys.exit()

        #***************************
        # Use irec to read EPC database
        #---------------------------
        a1surftype = np.load(dbDir + '/surfaceTypeIndex/%05d.npy'%(idx_db))[a1irec]

        #---------------------------
        for surftype in lsurftype:
            #-- Screen by surface types -------------------
            if surftype=='ocean':
                a1masksurf = ma.masked_not_equal(a1surftype,1).mask
            elif surftype=='vegetation':
                a1masksurf = ma.masked_outside(a1surftype,3,7).mask
            elif surftype=='snow':
                a1masksurf = ma.masked_outside(a1surftype,8,11).mask
            elif surftype=='coast':
                a1masksurf = ma.masked_not_equal(a1surftype,13).mask
            else:
                print '\n'+'check surftype',surftype
                sys.exit() 

            #-- Screeen no-precip cases for both datasets--
            a1mask1 = ma.masked_less_equal(a1ret, 0).mask
            a1mask2 = ma.masked_less_equal(a1obs, 0).mask
            a1mask  = a1mask1 * a1mask2
            a1mask  = a1mask + a1masksurf

            a1retTmp = ma.masked_where(a1mask, a1ret).compressed()
            a1obsTmp = ma.masked_where(a1mask, a1obs).compressed()
   
            if len(a1retTmp) ==0: continue

            #-- Histograms log-scale -------
            #-- shift 0.01mm/h for log-scale --
            logvmin, logvmax = -1.2, 2
            a1retTmp  = np.log10(a1retTmp + 0.01)
            a1obsTmp  = np.log10(a1obsTmp + 0.01)
            H,xedges,yedges = np.histogram2d(a1obsTmp, a1retTmp, bins = bins)
            H = H.T
            H = H * a1coef[idx_db]
            d2freq[surftype] = d2freq[surftype] + H

    #*******************************
    # Save
    #-------------------------------
    idx_db0 = min(lidx_db)
    idx_db1 = max(lidx_db)

    pickleDir = '/'.join(retbaseDir.split('/')[:-1] ) + '/pickle'
    histoPath  = pickleDir + '/histo.%s.%s.%05d-%05d.%s.bfile'%(expr, coefflag, idx_db0, idx_db1, rettype)
    binsPath   = pickleDir + '/bins.%s.%s.%05d-%05d.%s.npy'%(expr, coefflag, idx_db0, idx_db1, rettype)
    if calcflag == True:
        with open(histoPath, 'wb') as f:
            pickle.dump(d2freq, f)

        np.save(binsPath, bins)
    #*******************************
    # Load 
    #-------------------------------
    with open(histoPath, 'r') as f:
        d2freq = pickle.load(f)
    bins = np.load(binsPath)
    #*******************************
    # Figure
    #-------------------------------
    for surftype in lsurftype:
        H   = d2freq[surftype]
        X,Y = np.meshgrid(bins,bins) 
        #-- Figure density plot ----
        if rettype == lrettype[0]:
            dvnummax[surftype] = np.percentile(H,95)
         
    
        fig = plt.figure(figsize=[6,6])
        ax  = fig.add_axes([0.15,0.13,0.68,0.68])
        im  = ax.pcolormesh(X,Y,H, norm=matplotlib.colors.LogNorm(), cmap='jet', vmin=1, vmax=dvnummax[surftype])
    
        #-- plot 1:1 line
        ax.plot(array([logvmin,logvmax]),array([logvmin,logvmax]),'-',color='k',linewidth=0.5)
        #-- axis labels ---
        lticks = [-1,0,1,2]
        lticklabels = [0.1, 1, 10, 100]
    
        ax.set_xticks(lticks)
        ax.set_xticklabels(lticklabels, fontsize=16)
        ax.set_yticks(lticks)
        ax.set_yticklabels(lticklabels, fontsize=16)
   
        if rettype=='GPROF': 
            ax.set_xlabel('DPR/NScmb [mm/hour]', fontsize=22)
        else:
            ax.set_xlabel('DPR/%s [mm/hour]'%(rettype), fontsize=22)

        ax.set_ylabel('%s [mm/hour]'%(rettype), fontsize=22)
        ax.set_ylim([logvmin,logvmax])
        ax.set_xlim([logvmin,logvmax])
    
        plt.title('%s (%s)\n%s'%(rettype, expr, dsurflabel[surftype]), fontsize=24)
    
        cax = fig.add_axes([0.84,0.15,0.02, 0.6])
        cbar=plt.colorbar(im, orientation='vertical', cax=cax)
        cbar.ax.tick_params(labelsize=16)
    
        figPath= figDir + '/scatter.%s.%s.%05d-%05d.%s.%s.png'%(expr,coefflag, idx_db0, idx_db1, rettype,surftype)
        util.mk_dir(figDir)
    
        plt.savefig(figPath)
        print figPath
