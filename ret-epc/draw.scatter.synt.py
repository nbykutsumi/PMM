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
#coefflag  = 'nocoef'  # 'nocoef', 'wcoef'
coefflag  = 'wcoef'  # 'nocoef', 'wcoef'
DB_MAXREC = 10000
DB_MINREC = 1000
nsample   = 1000
#lsensor = ['GMI','AMSR2','SSMIS','ATMS','MHS']
lsensor = ['SSMIS','ATMS','MHS']
#lsensor = ['GMI']
lrettype = ['GPROF']  
#lrettype = ['MS','NS','NScmb','MScmb']
#lrettype = ['NScmb']
#prmin = 0.1
#prmin = 0.01
prmin = 0.0
lidx_db = range(29*29*29)[1:]
#lidx_db = range(810,812)
#lidx_db = range(6000,7000)
#lidx_db = range(29*29*4,29*29*8)
#lidx_db = range(29*29*8,29*29*12)
#lidx_db = range(29*29*12,29*29*16)
#lidx_db = range(29*29*16,29*29*20)
#lidx_db = range(29*29*20,29*29*24)
#lidx_db = range(29*29*24,29*29*29)
#** Constants ******

lsurftype = ['all','ocean','vegetation','coast','snow']
#lsurftype = ['ocean']

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
for (sensor,rettype) in [(sensor,rettype)
                for sensor  in lsensor
                for rettype in lrettype]:

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
    retbaseDir = tankDir + '/utsumi/PMM/retsynt/%s'%(expr)
    if dbtype=='jpl':
        dbDir = '/media/disk2/share/PMM/JPLDB/EPC_DB/%s_EPC_DATABASE'%(sensor)
    else:
        dbDir = tankDir + '/utsumi/PMM/EPCDB/samp.%05d.%s.V05A.S1.ABp103-117.01-12'%(DB_MAXREC, sensor)

    #***************************
    # Initialize histogram
    #---------------------------
    logvmin, logvmax = -1.2, 2
    bins   = np.arange(-2,logvmax+0.001,0.025)
    nbins  = len(bins)
    d2freq = {}

    lsatid = dlsatid[sensor]

    #** Read histogram of EPC DB records **
    if coefflag == 'wcoef':
        d1coef = {}
        for satid in [999] + lsatid:
            if ((dbtype=='my')and(sensor=='GMI')):
                countPath = countDir +'/count.epc.csv'
            elif dbtype=='jpl':
                countPath = countDir +'/count.jpldb.%s.%d.csv'%(sensor, satid)
            else:
                print 'check dbtype,sensor',dbtype,sensor
                sys.exit()

            f=open(countPath,'r'); lines=f.readlines(); f.close()
            ncol = len(lines[0].strip().split(',')) - 1
            d1coef[satid] = np.zeros(29*29*29, float32)
            
            for line in lines[1:]:
                line =  map(int,line.strip().split(','))
                epcid= line[0]
                nall = line[1]
                if nall < nsample:
                    coef = 1.0
                else:
                    coef = float(nall)/nsample
                d1coef[satid][epcid] = coef


 
    elif coefflag =='nocoef':
        d1coef[satid] = np.ones(29*29*29, float32)
    else:
        print 'check coefflag', coefflag
        sys.exit()
    
    #--------------------------------------
    for satid in [999] + lsatid:
        for surftype in lsurftype:
            d2freq[satid,surftype] = np.zeros([nbins-1,nbins-1],float64)  # x:obs  y:ret

    #---------------------------
    for idx_db in lidx_db:
        if calcflag ==False:
            continue

        if dbtype=='jpl':
            dbPath = dbDir + '/db_%05d.bin'%(idx_db)
            
            if not os.path.exists(dbPath):
                continue

            db.set_file(dbPath)

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

            if sensor in ['SSMIS','ATMS','MHS']:
                a1satid = np.load(retbaseDir + '/%05d/satid.obs.%05d.npy'%(idx_db,idx_db))


    
        elif rettype =='GPROF':
            if dbtype =='my':
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

            elif dbtype == 'jpl':
                irecPath = retbaseDir + '/%05d/%s.obs.%05d.npy'%(idx_db,'irec',idx_db)
                if not os.path.exists(irecPath):
                    print 'No file'
                    print irecPath
                    continue

                            
                a1irec= np.load(retbaseDir + '/%05d/%s.obs.%05d.npy'%(idx_db,'irec',idx_db)).astype(int32)
                a1obs = db.get_var('precip_NS_cmb')[a1irec] 
                a1ret = db.get_var('precip_GPROF')[a1irec]

                if sensor in ['SSMIS','ATMS','MHS']:
                    a1satid = np.load(retbaseDir + '/%05d/satid.obs.%05d.npy'%(idx_db,idx_db))



        else:
            print 'check rettype',rettype
            sys.exit()

        #***************************
        # Use irec to read EPC database
        #---------------------------
        if   dbtype=='jpl':
            a1surftype = db.get_var('sfc_class')[a1irec]

        elif dbtype=='my':
            a1surftype = np.load(dbDir + '/surfaceTypeIndex/%05d.npy'%(idx_db))[a1irec]



        #---------------------------
        for satid in [999] + lsatid:
            if sensor in ['GMI','AMSR2']:
                if satid ==999: continue
                else: a1masksat = np.array([False])

            else:
                if satid ==999:
                    a1masksat = np.array([False])
                else:
                    a1masksat = ma.masked_not_equal(a1satid, satid).mask

            for surftype in lsurftype:
                #-- Screen by surface types -------------------
                if   surftype=='all':
                    a1masksurf = np.array([False])
                elif surftype=='ocean':
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
                a1mask  = a1mask + a1masksurf + a1masksat
 

                if a1mask.sum() == 0:
                    a1retTmp = a1ret
                    a1obsTmp = a1obs

                else:
                    a1retTmp = ma.masked_where(a1mask, a1ret).compressed()
                    a1obsTmp = ma.masked_where(a1mask, a1obs).compressed()

                if a1retTmp.shape[0]==0:
                    #print 'a1retTmp.shape[0]=0: Skip',satid,surftype
                    continue


                #-- Histograms log-scale -------
                #-- shift 0.01mm/h for log-scale --
                logvmin, logvmax = -1.2, 2
                a1retTmp  = np.log10(a1retTmp + 0.01)
                a1obsTmp  = np.log10(a1obsTmp + 0.01)
                H,xedges,yedges = np.histogram2d(a1obsTmp, a1retTmp, bins = bins)
                H = H.T
                H = H * d1coef[satid][idx_db]
                d2freq[satid,surftype] = d2freq[satid,surftype] + H



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
    for satid in [999] + lsatid:
        for surftype in lsurftype:
            H   = d2freq[satid,surftype]
            X,Y = np.meshgrid(bins,bins) 

            print satid, surftype, H.max()    
            #-- Figure density plot ----
            if rettype == lrettype[0]:
                vpercentile = np.percentile(H,95)
    
                if vpercentile !=0:
                    dvnummax[satid, surftype] = vpercentile
                else:
                    dvnummax[satid, surftype] = H.max()
                
            
            fig = plt.figure(figsize=[6,6])
            ax  = fig.add_axes([0.15,0.13,0.66,0.66])
            im  = ax.pcolormesh(X,Y,H, norm=matplotlib.colors.LogNorm(), cmap='jet', vmin=1, vmax=dvnummax[satid, surftype])
   
            if H.max() ==0:
                print 'No precipitating case',satid, surftype
                print 'Skip'
                continue 
        
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
      
            satname = dsatname[satid] 
            stitle =  '%s (%s)\n%s %s'%(rettype, expr, satname, dsurflabel[surftype])
            plt.title(stitle, fontsize=24)
        
            cax = fig.add_axes([0.82,0.15,0.02, 0.6])
            cbar=plt.colorbar(im, orientation='vertical', cax=cax)
            cbar.ax.tick_params(labelsize=16)
       
            figPath= figDir + '/scatter.%s.%s.%s.%05d-%05d.%s.%s.png'%(expr, satname, coefflag, idx_db0, idx_db1, rettype,surftype)
            util.mk_dir(figDir)
        
            plt.savefig(figPath)
            print figPath
