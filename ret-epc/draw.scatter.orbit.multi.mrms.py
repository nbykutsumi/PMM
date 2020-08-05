
import matplotlib
matplotlib.use('Agg')
#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
import glob
import os, sys
import h5py
import myfunc.util as util
from datetime import datetime, timedelta

#----------------------------
gmi       = ["GPM","GMI","1C"]
amsr2     = ["GCOMW1","AMSR2"]
ssmis_f16 = ["F16","SSMIS"]
ssmis_f17 = ["F17","SSMIS"]      # 37V is missing after April 2016
ssmis_f18 = ["F18","SSMIS"]
atms_npp  = ["NPP","ATMS"]
atms_noaa20= ["NOAA20","ATMS"]   # MRMS is not available

mhs_metopa= ["METOPA","MHS"]
mhs_metopb= ["METOPB","MHS"]
mhs_noaa18= ["NOAA18","MHS"]
mhs_noaa19= ["NOAA19","MHS"]

iDTime = datetime(2018,1,24,0)
eDTime = datetime(2018,1,31,23)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
#lspec = [amsr2, ssmis_f16, ssmis_f17, ssmis_f18, atms_npp, atms_noaa20, mhs_metopa, mhs_metopb, mhs_noaa18, mhs_noaa19]
#lspec = [amsr2, ssmis_f16, ssmis_f18, atms_npp, mhs_metopa, mhs_metopb, mhs_noaa18, mhs_noaa19, gmi]
#lspec = [mhs_noaa18]
lspec = [amsr2]
lrettype=['epc']
mrmsbasedir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/MRMS/level2-pixel-match'
epcbasedir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/us.v01'

epcbasedir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/us.v01'
gprbasedir = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA'

tankdir = '/home/utsumi/mnt/lab_tank'
thpr = 0.1  # mm/h

lsurftype = ['ocean','vegetation','coast','snow']
#lsurftype = ['ocean']

dsurflabel={ 'all':'All surface'
            ,'ocean':'Class1 (Ocean)'
            ,'vegetation':'Classes 3-7 (Vegetation)'
            ,'snow':'Classes 8-11 (Snow)'
            ,'coast':'Class 13 (Coast)'
            }

dretname = {'epc':'EPC', 'gpr':'GPROF'}

#*********************************
# a1xpick
#---- AMSR2 ----
nxlw  = 243
a1xpick_amsr2 = (np.arange(nxlw)*2+2).astype('int16')
a1xpick_amsr2[-1] = nxlw*2-1   # [2,4,6,...482,484,485]

#---- SSMIS ----
nxlw = 90
a1xpick_ssmis = (np.arange(nxlw)*2+1).astype('int16') # [1,3,5,...,177,179]
#--------------------------


dvnummax = {}
for rettype,spec in [[rettype,spec] for spec in lspec for rettype in lrettype]:
    sate,sensor = spec

    #***********************************
    # Make MRMS list
    #-----------------------------------
    iY,iM = iDTime.timetuple()[:2]
    eY,eM = eDTime.timetuple()[:2]
    lYM   = util.ret_lYM([iY,iM],[eY,eM])

    lmrmspath = []
    lgoodpath = []
    liescan = []
    for (Year,Mon) in lYM:
        listPath = tankdir + '/utsumi/PMM/US/obtlist/overpass.%s.%s.%04d.%02d.csv'%(sensor,sate,Year,Mon)
        f=open(listPath,'r'); lines=f.readlines(); f.close()
        for line in lines:
            _,_,Day,oid,iscan,escan = map(int,line.strip().split(','))

            #if (sensor=='AMSR2')&(oid < 30388): continue # test

            if datetime(Year,Mon,Day) < iDTime: continue
            if eDTime < datetime(Year,Mon,Day): continue

            #--- MRMS -----
            ssearch = '/home/utsumi/mnt/lab_tank/utsumi/PMM/MRMS/level2-pixel-match/%s.%s.%s/%04d/%02d/%02d/mrms.%06d.*'%(rettype, sensor, sate,Year,Mon,Day,oid)
            lmrmspath_tmp= glob.glob(ssearch)

            ssearch = '/home/utsumi/mnt/lab_tank/utsumi/PMM/MRMS/level2-pixel-match/%s.%s.%s/%04d/%02d/%02d/goodfrac.%06d.*'%(rettype,sensor, sate,Year,Mon,Day,oid)
            lgoodpath_tmp= glob.glob(ssearch)

            if len(lmrmspath_tmp)==0: continue
            if len(lgoodpath_tmp)==0: continue

            mrmspath = lmrmspath_tmp[0]
            goodpath = lgoodpath_tmp[0] 

            lmrmspath.append(mrmspath)
            lgoodpath.append(goodpath)
            liescan.append([Year,Mon,Day,oid,iscan,escan])
            print mrmspath

    #**********************************
    # Start loop
    #----------------------------------
    d1obs  = {}
    d1pmw  = {}
    for surftype in lsurftype:
        d1obs[surftype]  = []
        d1pmw[surftype]  = []

    for mrmspath, goodpath, iescan in zip(lmrmspath, lgoodpath,liescan):

        Year,Mon,Day,oid,iscan,escan = iescan

        #**************************
        # Read PMW
        #--------------------------
        if rettype=='epc':
            ssearch = epcbasedir + '/%s.%s/%04d/%02d/%02d/nsurfNScmb.%06d.*.npy'%(sensor,sate,Year,Mon,Day,oid)
  
            lpmwpath = glob.glob(ssearch)
            print ssearch
            if len(lpmwpath)==0: continue
            pmwpath= lpmwpath[0]
            a2pmw = np.load(pmwpath)
            print a2pmw.shape

            #if (sensor=='AMSR2')&(oid < 30388): continue # test

        #**************************
        # Read surface type (and precip) from GPROF
        #--------------------------
        ssearch = gprbasedir + '/%s.%s/2A-CLIM/V05/%04d/%02d/%02d/2A-CLIM.*.%06d.????.HDF5'%(sate,sensor,Year,Mon,Day,oid)
        lgprpath = glob.glob(ssearch)

        if len(lgprpath)==0:
            print 'No GPROF file for',Year,Mon,Day,oid
            continue
        gprpath= lgprpath[0]

        with h5py.File(gprpath,'r') as h:
            a2surftype = h['S1/surfaceTypeIndex'][iscan:escan+1]

            if rettype=='gpr':
                a2pmw = h['S1//surfacePrecipitation'][iscan:escan+1]


        #-- if rettype=='epc' ----
        if rettype=='epc':
            if sensor in ['GMI','ATMS','MHS']:
                pass
            elif sensor=='AMSR2':
                a2surftype = a2surftype[:, a1xpick_amsr2]

            elif sensor=='SSMIS':
                a2surftype = a2surftype[:, a1xpick_ssmis]

            else:
                print 'check sensor',sensor
                sys.exit()

        elif rettype=='gpr':
            pass
        else:
            print 'rettype',rettype 

        print a2surftype.shape, a2pmw.shape
        #**************************
        # Read MRMS
        #--------------------------
        a2obs  = np.load(mrmspath)
        a2good = np.load(goodpath)

        a2obs = ma.masked_less(a2obs, thpr).filled(0)
        a2pmw = ma.masked_less(a2pmw, thpr).filled(0)

        a2mask1 = ma.masked_less(a2obs, 0).mask
        a2mask2 = ma.masked_less(a2pmw, 0).mask
        a2mask3 = ma.masked_less(a2good, 0.5).mask

        a2mask4 = ma.masked_equal(a2pmw,0).mask * ma.masked_equal(a2obs,0).mask


        for surftype in lsurftype:
            if surftype=='all':
                a2masksurf = np.array([False])
            elif surftype=='ocean':
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
            

            a2mask  = a2mask1 + a2mask2 + a2mask3 + a2mask4 + a2masksurf


            d1obs[surftype].extend( ma.masked_where(a2mask, a2obs).compressed())
            d1pmw[surftype].extend( ma.masked_where(a2mask, a2pmw ).compressed())


    for surftype in lsurftype:
        a1obs = np.array(d1obs[surftype])
        a1pmw = np.array(d1pmw[surftype])

        #***************************
        # Histogram
        #---------------------------
        cc   = np.corrcoef(a1obs, a1pmw)[0,1]
        rmse = np.sqrt(((a1pmw - a1obs)**2).mean())
        rbias = (a1pmw - a1obs).mean() / a1obs.mean()

        logvmin, logvmax = -1, 2
        #a1pmw = np.log10(a1pmw+ 0.01)
        #a1obs = np.log10(a1obs+ 0.01)

        a1pmw = np.log10(a1pmw)
        a1obs = np.log10(a1obs)

        bins  = np.arange(-2,logvmax+0.001,0.025)
        H,xedges,yedges = np.histogram2d(a1obs, a1pmw, bins = bins)
        H = H.T

        X,Y = np.meshgrid(bins,bins)

        if rettype==lrettype[0]:
            dvnummax[surftype] = np.percentile(H,99)

        fig = plt.figure(figsize=[7,6])
        ax  = fig.add_axes([0.15,0.13,0.68,0.68])
        im  = ax.pcolormesh(X,Y,H, norm=matplotlib.colors.LogNorm(), cmap='jet', vmin=1, vmax=dvnummax[surftype])


        ##-- Text ----------
        #t=ax.text(0.02, 0.90, 'RMSE:%.2f'%(rmse), size=32, fontweight='bold', path_effects=[pe.withStroke(linewidth=5, foreground='w')], transform=ax.transAxes)
        #t=ax.text(0.02, 0.80, 'CC  :%.2f'%(cc), size=32,fontweight='bold', path_effects=[pe.withStroke(linewidth=5, foreground='w')], transform=ax.transAxes)

        #-- plot 1:1 line
        ax.plot(array([logvmin,logvmax]),array([logvmin,logvmax]),'-',color='k',linewidth=0.5)
        #-- axis labels ---
        lticks = [-1,0,1,2]
        lticklabels = [0.1, 1, 10, 100]

        #lticks = np.arange(-1,0.5+0.01, 0.2)
        #lticklabels = ['%.1f'%(10**k) for k in np.arange(-1,0.5+0.01,0.2)]

        ax.set_xticks(lticks)
        ax.set_xticklabels(lticklabels, fontsize=16)
        ax.set_yticks(lticks)
        ax.set_yticklabels(lticklabels, fontsize=16)


        ax.set_xlabel('MRMS [mm/hour]', fontsize=22)
        ax.set_ylabel('%s [mm/hour]'%(dretname[rettype]), fontsize=22)
        ax.set_ylim([logvmin,logvmax])
        ax.set_xlim([logvmin,logvmax])

        plt.title('%s %s'%(dretname[rettype], dsurflabel[surftype]), fontsize=24)

        cax = fig.add_axes([0.84,0.15,0.02, 0.6])
        cbar=plt.colorbar(im, orientation='vertical', cax=cax)
        cbar.ax.tick_params(labelsize=16)

        figDir = '/home/utsumi/temp/multi/scatter'
        figPath= figDir + '/scatter.MRMS.%s.%s.%s.png'%(sensor,sate,surftype)
        util.mk_dir(figDir)

        plt.savefig(figPath)
        print figPath


        #a1obs  = np.array(a1obs ) 
        #a1pmw  = np.array(a1pmw ) 
        #plt.scatter(a1obs, a1pmw)
        #plt.ylim([0,10])
        #plt.xlim([0,10])
        #plt.show()
        #print len(a1obs)
        #print len(a1pmw)

        #sys.exit()


# %%
