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

calcflag  = True
#calcflag  = False
#coefflag  = 'nocoef'  # 'nocoef', 'wcoef'
coefflag  = 'wcoef'  # 'nocoef', 'wcoef'
DB_MAXREC = 10000
DB_MINREC = 1000
nsample   = 1000
dbtype = 'my'
#prmin = 0.1
#prmin = 0.01
prmin = 0.0
nz   = 50
expr = 'org.smp%d'%(nsample)
lidx_db = range(29*29*29)[1:]
#lidx_db = range(29*29*29)[1:2000]
#lidx_db = [1681]
#lidx_db = [10820]

#** Constants ******
stamp0  = 'region'
nkeys   = 4
lthpr   = [0.1, 10]
#lseason = ['ALL','DJF','JJA']
lseason = np.arange(1,12+1)
#lseason_fig= ['ALL','DJF','JJA']
#lseason = ['ALL']
#lseason = [6]
lstop   = ['All','High','Mid','Low']
#lstop   = ['Low']

dsumtoprange= { 'All':[0,30000], 'Low':[0,4000],'Mid':[4000,6000],'High':[6000,30000]}
lregion = ['AMZ','CUS','EUS','TIB','NETP','SETP','NTP','STP','WTP','ETI','WMP','WMA','TAF','NEA']
#lregion = ['CUS','TIB']
#lregion = ['GLB']
dBBox = {
         'GLB':  [[-65,-180],[65,180]]
        ,'NHM':  [[0,-180],[65,180]]
        ,'AMZ':  [[-5,-65],[5,-55]]
        ,'CUS':  [[35,-105],[45,-95]]
        ,'EUS':  [[30,-90],[40,-80]]
        ,'TIB':  [[30,80],[40,90]]
        ,'NETP':  [[0,-120],[10,-110]]
        ,'SETP':  [[10,-120],[0,-110]]
        ,'NTP' : [[0,-35],[10,-25]]
        ,'STP' : [[-10,-35],[0,-25]]
        ,'WTP':  [[0,140],[10,150]]
        ,'ETI':  [[-5,85],[5,95]]
        ,'WMP':  [[22,135],[32,145]]
        ,'WMA':  [[25,-70],[35,-60]]
        ,'TAF':  [[0, 15],[10,25]]
        ,'NEA':  [[40,120],[50,130]]
        }



myhost = socket.gethostname()
if myhost =='shui':
    dbDir   = '/work/hk01/utsumi/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
    countDir= '/work/hk01/utsumi/PMM/EPCDB/list'
    retbaseDir = '/tank/utsumi/PMM/retsynt/%s'%(expr)
    tankbaseDir= '/tank'
    figDir   = '/home/utsumi/temp/ret'


elif myhost == 'well':
    dbDir   = '/media/disk2/share/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
    countDir= '/media/disk2/share/PMM/EPCDB/list'
    retbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retsynt/%s'%(expr)
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    figDir   = '/home/utsumi/temp/ret'

else:
    print 'check hostname',myhost
    sys.exit()


lsurftype = ['ocean','vegetation','coast']
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


#***************************
#* Elevation **
a2orog = np.load(tankbaseDir + '/utsumi/PMM/validprof/const/orog.meter.sp.dec.1800x3600.npy')
dorog = {}
for region in lregion:
    [[lat0,lon0],[lat1,lon1]] = dBBox[region]
    iy0 = int(floor((lat0-(-90))/0.1))
    iy1 = int(floor((lat1-(-90))/0.1))
    ix0 = int(floor((lon0-(-180))/0.1))
    ix1 = int(floor((lon1-(-180))/0.1))

    dorog[region] = ma.masked_less(a2orog[iy0:iy1+1,ix0:ix1+1], -1000).mean()
#***************************
# Initialize
#---------------------------
dsobs = {}
dsobs2= {}
dnobs = {}

dsest = {}
dsest2= {}
dnest = {}

dsumtop = {}
dsumtop2= {}
dnumtop = {}

dsstopobs= {}
dsstopest= {}

dsprecobs= {}
dsprecest= {}

#---------------------------
for idx_db in lidx_db:
    if calcflag is not True:
        continue

    coef = a1coef[idx_db]
    #irecobsPath = retbaseDir + '/%05d/irec.obs.%05d.npy'%(idx_db,idx_db)
    #if not os.path.exists(irecobsPath):
    #    print 'No file'
    #    print irecobsPath
    #    continue

    irectopPath = retbaseDir + '/%05d/irec.top.%05d.npy'%(idx_db,idx_db)
    if not os.path.exists(irectopPath):
        print 'No file', idx_db
        #print irectopPath
        continue
    print 'idx_db=',idx_db

    a1irec= np.load(retbaseDir + '/%05d/irec.obs.%05d.npy'%(idx_db,idx_db)).astype(int32)

    #a1irectop = np.load(retbaseDir + '/%05d/irec.top.%05d.npy'%(idx_db,idx_db)).astype(int32)
    #a1idxdbtop = np.load(retbaseDir + '/%05d/idx_db.top.%05d.npy'%(idx_db,idx_db)).astype(int32)

    a1precobs = np.load(retbaseDir + '/%05d/nsurfNScmb.obs.%05d.npy'%(idx_db,idx_db))
    a1latobs  = np.load(retbaseDir + '/%05d/Latitude.obs.%05d.npy'%(idx_db,idx_db))   
    a1lonobs  = np.load(retbaseDir + '/%05d/Longitude.obs.%05d.npy'%(idx_db,idx_db))   
    #a1stopobs = np.load(retbaseDir + '/%05d/Ku_NS_heightStormTop.obs.%05d.npy'%(idx_db,idx_db))   

    a2prwatobs = np.load(retbaseDir + '/%05d/precip_water_prof_NS.obs.%05d.npy'%(idx_db,idx_db))[:,-nz:][:,::-1]  # bottom to top
    a2prwatest = np.load(retbaseDir + '/%05d/precip_water_prof_NS.est.%05d.npy'%(idx_db,idx_db))[:,-nz:][:,::-1]  # bottom to top
    #a2prwattop = np.load(retbaseDir + '/%05d/precip_water_prof_NS.top.%05d.npy'%(idx_db,idx_db))[:,-nz:][:,::-1]  # bottom to top

    a1surftypeobs= np.load(retbaseDir + '/%05d/surfaceTypeIndex.obs.%05d.npy'%(idx_db,idx_db))   
    a1stopobs    = np.load(retbaseDir + '/%05d/Ku_NS_heightStormTop.obs.%05d.npy'%(idx_db,idx_db))   
    a1prtypeobs  = np.load(retbaseDir + '/%05d/Ku_NS_typePrecip.obs.%05d.npy'%(idx_db,idx_db))   

    a1mdhmsobs   = np.load(retbaseDir + '/%05d/mdhms.obs.%05d.npy'%(idx_db,idx_db))
    a1monobs     = a1mdhmsobs[:,0] 

    #** Invalid values **
    a2prwatobs = ma.masked_invalid(a2prwatobs).filled(-9999.)
    a2prwatest = ma.masked_invalid(a2prwatest).filled(-9999.)
    #a2prwattop = ma.masked_invalid(a2prwattop).filled(-9999.)

    #** Make bit and replace missing by zero **
    a2bitobs    = ma.masked_greater_equal(a2prwatobs,0).mask.astype(int32) * coef
    a2bitest    = ma.masked_greater_equal(a2prwatest,0).mask.astype(int32) * coef

    #a2bittop    = (ma.masked_greater_equal(a2prwattop,0)*0).filled(1).astype(int32) * coef

    a2prwatobs  = ma.masked_less(a2prwatobs,0).filled(0.) * coef
    a2prwatest  = ma.masked_less(a2prwatest,0).filled(0.) * coef
    #a2prwattop  = ma.masked_less(a2prwattop,0).filled(0.) * coef

    lkey = [(region,season,thpr,stop) for region in lregion
                            for season in lseason
                            for thpr   in lthpr
                            for stop   in lstop
                            ]

    if len(lkey[0]) !=nkeys:
        print 'check nkeys',nkeys
        sys.exit()


    for key in lkey:
        (region,season,thpr,stop) = key

        #* Region **
        [[lat0,lon0],[lat1,lon1]] = dBBox[region]
        a1flaglat = ma.masked_inside(a1latobs,lat0,lat1).mask
        a1flaglon = ma.masked_inside(a1lonobs,lon0,lon1).mask
        a1flaglatlon = a1flaglat * a1flaglon

        #* Season **
        if season =='ALL':
            a1flagmon = True
        else:
            lMon = util.ret_lmon(season)
            a1flagmon = False
            for Mon in lMon:
                a1flagmon = a1flagmon + ma.masked_equal(a1monobs, Mon).mask

        #* Precipitation **
        a1flagpr = ma.masked_greater_equal(a1precobs, thpr).mask
        
        #* Storm top height **
        stop0, stop1 = dsumtoprange[stop]
        a1flagstop = ma.masked_inside(a1stopobs, stop0, stop1).mask
        #---------------
        a1flag = a1flaglatlon * a1flagmon * a1flagpr * a1flagstop


        if a1flag is np.bool_(False): continue
        if a1flag.sum()==0: continue
        
        a1prwatobsTmp  = a2prwatobs[a1flag].sum(axis=0)
        a1prwatestTmp  = a2prwatest[a1flag].sum(axis=0)
        #a1prwattopTmp  = a2prwattop[a1flag].sum(axis=0)

        a1bitobsTmp    = a2bitobs[a1flag].sum(axis=0)
        a1bitestTmp    = a2bitest[a1flag].sum(axis=0)
        #a1bittopTmp    = a2bittop[a1flag].sum(axis=0)

        if key not in dsobs.keys():
            dcount[key]   = 1

            dsobs [key]   = a1prwatobsTmp 
            dsobs2[key]   = a1prwatobsTmp**2
            dnobs [key]   = a1bitobsTmp

            dsest [key]   = a1prwatestTmp 
            dsest2[key]   = a1prwatestTmp**2
            dnest [key]   = a1bitestTmp

            #dsumtop [key]   = a1prwattopTmp
            #dsumtop2[key]   = a1prwattopTmp**2
            #dnumtop [key]   = a1bittopTmp
    
        else:
            dcount[key]  += 1

            dsobs [key]  += a1prwatobsTmp 
            dsobs2[key]  += a1prwatobsTmp**2
            dnobs [key]  += a1bitobsTmp

            dsest [key]  += a1prwatestTmp 
            dsest2[key]  += a1prwatestTmp**2
            dnest [key]  += a1bitestTmp

            #dsumtop [key]  += a1prwattopTmp
            #dsumtop2[key]  += a1prwattopTmp**2
            #dnumtop [key]  += a1bittopTmp

##*******************************
## Save
##-------------------------------
stamp = '%dkey.%s'%(nkeys, stamp0)
idx_db0 = lidx_db[0]
idx_db1 = lidx_db[-1]
pickleDir = '/'.join(retbaseDir.split('/')[:-1] ) + '/pickle/%s'%(expr)
util.mk_dir(pickleDir)

countPath  = pickleDir + '/count.%s.%05d-%05d.bfile' %(stamp, idx_db0, idx_db1)

sobsPath  = pickleDir + '/sum.obs.%s.%05d-%05d.bfile' %(stamp, idx_db0, idx_db1)
sobs2Path = pickleDir + '/sum2.obs.%s.%05d-%05d.bfile'%(stamp, idx_db0, idx_db1)
nobsPath  = pickleDir + '/num.obs.%s.%05d-%05d.bfile' %(stamp, idx_db0, idx_db1)

sestPath  = pickleDir + '/sum.est.%s.%05d-%05d.bfile' %(stamp, idx_db0, idx_db1)
sest2Path = pickleDir + '/sum2.est.%s.%05d-%05d.bfile'%(stamp, idx_db0, idx_db1)
nestPath  = pickleDir + '/num.est.%s.%05d-%05d.bfile' %(stamp, idx_db0, idx_db1)

#stopPath  = pickleDir + '/sum.top.%s.%05d-%05d.bfile' %(stamp, idx_db0, idx_db1)
#stop2Path = pickleDir + '/sum2.top.%s.%05d-%05d.bfile'%(stamp, idx_db0, idx_db1)
#ntopPath  = pickleDir + '/num.top.%s.%05d-%05d.bfile' %(stamp, idx_db0, idx_db1)

if calcflag == True:
    with open(countPath, 'wb') as f:
        pickle.dump(dcount, f)

    with open(sobsPath, 'wb') as f:
        pickle.dump(dsobs, f)

    with open(sobs2Path, 'wb') as f:
        pickle.dump(dsobs2, f)

    with open(nobsPath, 'wb') as f:
        pickle.dump(dnobs, f)

    with open(sestPath, 'wb') as f:
        pickle.dump(dsest, f)

    with open(sest2Path, 'wb') as f:
        pickle.dump(dsest2, f)

    with open(nestPath, 'wb') as f:
        pickle.dump(dnest, f)

    #with open(stopPath, 'wb') as f:
    #    pickle.dump(dsumtop, f)

    #with open(stop2Path, 'wb') as f:
    #    pickle.dump(dsumtop2, f)

    #with open(ntopPath, 'wb') as f:
    #    pickle.dump(dnumtop, f)


##*******************************
## Load 
##-------------------------------
print '---------------Load _--------------'

with open(countPath, 'r') as f:
    dcount = pickle.load(f)

with open(sobsPath, 'r') as f:
    dsobs = pickle.load(f)

with open(sobs2Path, 'r') as f:
    dsobs2 = pickle.load(f)

with open(nobsPath, 'r') as f:
    dnobs = pickle.load(f)

with open(sestPath, 'r') as f:
    dsest = pickle.load(f)

with open(sest2Path, 'r') as f:
    dsest2 = pickle.load(f)

with open(nestPath, 'r') as f:
    dnest = pickle.load(f)

#with open(stopPath, 'r') as f:
#    dsumtop = pickle.load(f)
#
#with open(stop2Path, 'r') as f:
#    dsumtop2 = pickle.load(f)
#
#with open(ntopPath, 'r') as f:
#    dnumtop = pickle.load(f)


#for season in lseason_fig:
#    lmon = 
#    lkey = [(region,season,thpr,stop) for region in lregion
#                            for season in lseason
#                            for thpr   in lthpr
#                            for stop   in lstop
#                            ]
#
#sys.exit()
print '------------- figure ---------------'
lkey = dsobs.keys()
for key in lkey:    
    print key
    region, season, thpr, stop = key

    #if season !=6: continue
    fig = plt.figure(figsize=(2.5,6))
    ax  = fig.add_axes([0.3,0.15,0.65,0.7])

    a1y = 0.25 + np.arange(nz) * 0.25 # [km]
    a1obs = ma.masked_invalid(dsobs[key] / dnobs[key])
    a1est = ma.masked_invalid(dsest[key] / dnest[key])
    #a1top = ma.masked_invalid(dsumtop[key] / dnumtop[key])

    #-- Mask lowest 1.5km ----
    a1y   = a1y[6:]
    a1obs = a1obs[6:]
    a1est = a1est[6:]

    a1nobs= dnobs[key][6:]
    a1nest= dnest[key][6:]
    #-----------------------
    for i,y in enumerate(a1y):
        print i,y,'%.3f, %.3f, %.3f, %.3f'%( a1obs[i], a1est[i], a1nobs[i], a1nest[i])

    ax.plot( a1obs, a1y, '-', c='k', linewidth=2, label='CMB')
    ax.plot( a1est, a1y, '-', c='k', linewidth=1, label='PMW')
    #ax.plot( a1top, a1y, '--', c='k', linewidth=1, label='PMW 1st')

    ax.axhline(y=dorog[region]*0.001, color='k', linestyle='-.')

    stitle = '%s %s %s'%(region, stop, season)
    plt.title(stitle)
    plt.legend()

    plt.ylim([0,12])
    plt.xlabel('precipitation water (g/m3)')
    plt.ylabel('hight (above sea level) (km)')
    figPath = figDir + '/prof.synt.pr%.1f.%s.%s.sh-%s.%05d-%05d.png'%(thpr, region, season, stop, lidx_db[0], lidx_db[-1])
    plt.savefig(figPath)
    print figPath



