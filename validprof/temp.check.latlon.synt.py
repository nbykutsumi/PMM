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
#lidx_db = range(29*29*29)[1:]
lidx_db = range(29*29*29)[1:5000]
#lidx_db = [1681]
#lidx_db = [10820]

#** Constants ******
#lregion = ['GLB']
#lregion = ['NHM']
lregion = ['AMZ','CUS','EUS','TIB','NETP','SETP','NTP','STP','WTP','ETI','WMP','WMA','TAF','NEA']
dBBox = {
         'GLB':  [[-65,-180],[65,180]]
        ,'NHM':  [[0,-180],[65,180]]
        ,'AMZ':  [[-5,-65],[5,-55]]
        ,'CUS':  [[35,-105],[45,-95]]
        ,'EUS':  [[30,-90],[40,-80]]
        ,'TIB':  [[30,-90],[40,-80]]
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
# Initialize
#---------------------------
dnobs = {}

#---------------------------
for idx_db in lidx_db:
    print 'idx_db=',idx_db

    #irecobsPath = retbaseDir + '/%05d/irec.obs.%05d.npy'%(idx_db,idx_db)
    #if not os.path.exists(irecobsPath):
    #    print 'No file'
    #    print irecobsPath
    #    continue

    irectopPath = retbaseDir + '/%05d/irec.top.%05d.npy'%(idx_db,idx_db)
    if not os.path.exists(irectopPath):
        print 'No file'
        print irectopPath
        continue

    a1latobs  = np.load(retbaseDir + '/%05d/Latitude.obs.%05d.npy'%(idx_db,idx_db))   
    a1lonobs  = np.load(retbaseDir + '/%05d/Longitude.obs.%05d.npy'%(idx_db,idx_db))   

    for region in lregion:
        [[lat0,lon0],[lat1,lon1]] = dBBox[region]

        a1flaglat = ma.masked_inside(a1latobs,lat0,lat1).mask
        a1flaglon = ma.masked_inside(a1lonobs,lon0,lon1).mask
        a1flaglatlon = a1flaglat * a1flaglon

        a1flag = a1flaglatlon
        if a1flag is np.bool_(False): continue

        a1tmp = a1latobs[a1flag]
        key = (region)
        if key not in dnobs.keys():
            dnobs [key]   = len(a1tmp)
       
        else:
            dnobs [key]  += len(a1tmp)
print dnobs
for region in lregion:
    print region,dnobs[region]


    ##---------------------------
    #for surftype in lsurftype:
    #    #-- Screen by surface types -------------------
    #    if surftype=='ocean':
    #        a1masksurf = ma.masked_not_equal(a1surftype,1).mask
    #    elif surftype=='vegetation':
    #        a1masksurf = ma.masked_outside(a1surftype,3,7).mask
    #    elif surftype=='snow':
    #        a1masksurf = ma.masked_outside(a1surftype,8,11).mask
    #    elif surftype=='coast':
    #        a1masksurf = ma.masked_not_equal(a1surftype,13).mask
    #    else:
    #        print '\n'+'check surftype',surftype
    #        sys.exit() 

    #    #-- Screeen no-precip cases for both datasets--
    #    a1mask1 = ma.masked_less_equal(a1ret, 0).mask
    #    a1mask2 = ma.masked_less_equal(a1obs, 0).mask
    #    a1mask  = a1mask1 * a1mask2
    #    a1mask  = a1mask + a1masksurf

    #    a1retTmp = ma.masked_where(a1mask, a1ret).compressed()
    #    a1obsTmp = ma.masked_where(a1mask, a1obs).compressed()

    #    if len(a1retTmp) ==0: continue


##*******************************
## Save
##-------------------------------
#idx_db0 = min(lidx_db)
#idx_db1 = max(lidx_db)
#
#pickleDir = '/'.join(retbaseDir.split('/')[:-1] ) + '/pickle'
#histoPath  = pickleDir + '/histo.%s.%s.%05d-%05d.%s.bfile'%(expr, coefflag, idx_db0, idx_db1, rettype)
#binsPath   = pickleDir + '/bins.%s.%s.%05d-%05d.%s.npy'%(expr, coefflag, idx_db0, idx_db1, rettype)
#if calcflag == True:
#    with open(histoPath, 'wb') as f:
#        pickle.dump(d2freq, f)
#
#    np.save(binsPath, bins)
##*******************************
## Load 
##-------------------------------
#with open(histoPath, 'r') as f:
#    d2freq = pickle.load(f)
#bins = np.load(binsPath)
##*******************************
## Figure
##-------------------------------
#for surftype in lsurftype:
#    H   = d2freq[surftype]
#    X,Y = np.meshgrid(bins,bins) 
#    #-- Figure density plot ----
#    dvnummax[surftype] = np.percentile(H.max(),70)
#
#    fig = plt.figure(figsize=[6,6])
#    ax  = fig.add_axes([0.15,0.13,0.68,0.68])
#    im  = ax.pcolormesh(X,Y,H, norm=matplotlib.colors.LogNorm(), cmap='jet', vmin=1, vmax=dvnummax[surftype])
#
#    #-- plot 1:1 line
#    ax.plot(array([logvmin,logvmax]),array([logvmin,logvmax]),'-',color='k',linewidth=0.5)
#    #-- axis labels ---
#    lticks = [-1,0,1,2]
#    lticklabels = [0.1, 1, 10, 100]
#
#    ax.set_xticks(lticks)
#    ax.set_xticklabels(lticklabels, fontsize=16)
#    ax.set_yticks(lticks)
#    ax.set_yticklabels(lticklabels, fontsize=16)
#
#    if rettype=='GPROF': 
#        ax.set_xlabel('DPR/NScmb [mm/hour]', fontsize=22)
#    else:
#        ax.set_xlabel('DPR/%s [mm/hour]'%(rettype), fontsize=22)
#
#    ax.set_ylabel('%s [mm/hour]'%(rettype), fontsize=22)
#    ax.set_ylim([logvmin,logvmax])
#    ax.set_xlim([logvmin,logvmax])
#
#    plt.title('%s (%s)\n%s'%(rettype, expr, dsurflabel[surftype]), fontsize=24)
#
#    cax = fig.add_axes([0.84,0.15,0.02, 0.6])
#    cbar=plt.colorbar(im, orientation='vertical', cax=cax)
#    cbar.ax.tick_params(labelsize=16)
#
#    figPath= figDir + '/scatter.%s.%s.%05d-%05d.%s.%s.png'%(expr,coefflag, idx_db0, idx_db1, rettype,surftype)
#    util.mk_dir(figDir)
#
#    plt.savefig(figPath)
#    print figPath
