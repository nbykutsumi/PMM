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


myhost = socket.gethostname()
if myhost == 'shui':
    listDir    = '/work/hk01/utsumi/PMM/US/obtlist'
    srcbaseDir = '/tank/utsumi/PMM/retepc'
    gprofbaseDir = '/work/hk01/PMM/NASA/GPM.GMI/2A/V05'

    mrmsDir  = '/work/hk01/PMM/MRMS/match-GMI-orbit'
    figDir   = '/home/utsumi/temp/ret'

elif myhost == 'well':
    listDir    = '/home/utsumi/mnt/lab_tank/utsumi/PMM/US/obtlist'
    srcbaseDir = '/media/disk2/share/PMM/retepc'
    gprofbaseDir = '/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.GMI/2A/V05'
    mrmsDir  = '/home/utsumi/mnt/lab_work/hk01/PMM/MRMS/match-GMI-orbit'
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

iYM = [2014,7]
eYM = [2014,8]
lYM = util.ret_lYM(iYM,eYM)
lseason = ['JJA']
nrec = 20000
#lrettype = ['NS','MS','NScmb','MScmb']
lrettype = ['NS','MS','NScmb','MScmb','GPROF']
#lrettype = ['GPROF'] + 
#lrettype = ['NS']
#prmin = 0.1
prmin = 0.01
expr  = 'glb.wprof.org'

def read_orbitlist(Year,Mon):
    listPath = listDir + '/overpass.GPM.%04d.%02d.csv'%(Year,Mon)
    f=open(listPath); lines=f.readlines(); f.close()
    lout=[]
    for line in lines:
        line = map(int,line.strip().split(','))
        lout.append(line)
    return lout


for season in lseason:
    lMon = util.ret_lmon(season)
    dvnummax = {}
    for rettype in lrettype:
        for surftype in lsurftype:
    
            a1ret = array([])
            a1mrms= array([])
            for Year,Mon in lYM:
                if Mon not in lMon:
                    continue

                print surftype,Year,Mon

                lorbit = read_orbitlist(Year,Mon)
                for orbit in lorbit:
                    Day,oid,iy,ey = orbit[2:]
                    if Day>23: continue   # test
                    #print Year,Mon,Day,oid 
    
        
                    #-- Read surface type ----
                    ssearch = gprofbaseDir + '/%04d/%02d/%02d/2A.GPM.GMI.GPROF2017v1.20140701-S030429-E043702.001921.V05A.HDF5'
                    ssearch = gprofbaseDir + '/%04d/%02d/%02d/2A.GPM.GMI.GPROF2017v1.*.%06d.V05A.HDF5'%(Year,Mon,Day,oid)
                    gprofPath= glob.glob(ssearch)[0]

                    with h5py.File(gprofPath,'r') as h:
                        a2surftype = h['S1/surfaceTypeIndex'][iy:ey+1]

                    #-- Read retrieved precipitation ----
                    if rettype =='GPROF':
                        with h5py.File(gprofPath,'r') as h:
                            a2ret = h['S1/surfacePrecipitation'][iy:ey+1]

                    else:
                        retDir = srcbaseDir + '/%s/%04d/%02d/%02d'%(expr,Year,Mon,Day)
                        a2ret  = np.load(retDir + '/nsurf%s.%06d.y-9999--9999.nrec%d.npy'%(rettype, oid, nrec))[iy:ey+1]
                        a2ret = ma.masked_less(a2ret, prmin).filled(0.0) 


                    #-- Read MRSM over the orbit ------
                    #a2mrms = np.load(mrmsDir + '/GMI.MRMS.130W_55W_20N_55N.20141021.003663.00931-01189.npy')
                    ssearch = mrmsDir + '/GMI.MRMS.130W_55W_20N_55N.????????.%06d.?????-?????.npy'%(oid)
                    mrmsPath= glob.glob(ssearch)[0]
                    a2mrms  = np.load(mrmsPath)
                    a2mrms  = ma.masked_less(a2mrms,0).filled(0) 
                    
                    #-- Screen by surface types -------------------
                    if surftype=='ocean':
                        a2masksurf = ma.masked_not_equal(a2surftype,1).mask
                    elif surftype=='vegetation':
                        a2masksurf = ma.masked_outside(a2surftype,3,7).mask
                    elif surftype=='snow':
                        a2masksurf = ma.masked_outside(a2surftype,8,11).mask
                    elif surftype=='coast':
                        a2masksurf = ma.masked_not_equal(a2surftype,13).mask
                    else:
                        print '\n'+'check surftype',surftype
                        sys.exit() 
                    #-- Screeen no-precip cases for both datasets--
                    a2mask1 = ma.masked_less_equal(a2ret, 0).mask
                    a2mask2 = ma.masked_less_equal(a2mrms,0).mask
                    a2mask  = a2mask1 * a2mask2
                    a2mask  = a2mask + a2masksurf

                    a1retTmp  = ma.masked_where(a2mask, a2ret).compressed()
                    a1mrmsTmp = ma.masked_where(a2mask, a2mrms).compressed()
    
                    a1ret  = np.concatenate([a1ret, a1retTmp])
                    a1mrms = np.concatenate([a1mrms,a1mrmsTmp])
    
    
    
            #-- Histograms log-scale -------
            #-- shift 0.01mm/h for log-scale --
            logvmin, logvmax = -1.2, 2
            a1ret  = np.log10(a1ret + 0.01)
            a1mrms = np.log10(a1mrms+ 0.01)
            bins   = np.arange(-2,logvmax+0.001,0.025)
            H,xedges,yedges = np.histogram2d(a1mrms, a1ret, bins = bins)
            H = H.T
            H = ma.masked_equal(H,0)
            X,Y = np.meshgrid(bins,bins) 
            #-- Figure density plot ----
            if rettype==lrettype[0]:
                dvnummax[surftype] = np.percentile(H.max(),70)

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
    
            ax.set_xlabel('MRMS [mm/hour]', fontsize=22)
            ax.set_ylabel('%s [mm/hour]'%(rettype), fontsize=22)
            ax.set_ylim([logvmin,logvmax])
            ax.set_xlim([logvmin,logvmax])
    
            plt.title('%s\n%s %s'%(rettype, dsurflabel[surftype],season), fontsize=24)
    
            cax = fig.add_axes([0.84,0.15,0.02, 0.6])
            cbar=plt.colorbar(im, orientation='vertical', cax=cax)
            cbar.ax.tick_params(labelsize=16)
    
            figPath= figDir + '/scatter.log.%s.%s.%s.png'%(rettype,surftype,season)
            util.mk_dir(figPath)

            plt.savefig(figPath)
            print figPath
