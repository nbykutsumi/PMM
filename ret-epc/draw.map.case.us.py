import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
from numpy import *
import h5py
from bisect import bisect_left
import epcfunc
import sys, os, glob, socket, glob
import myfunc.util as util
from datetime import datetime, timedelta

iDTime = datetime(2018,1,1)
eDTime = datetime(2018,1,31)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))
expr   = 'us.v01'
thpr = 0.5  # Minimum threshold for EPC
[[lllat,lllon],[urlat,urlon]] = [[20,-130],[55,-55]]  # MRMS

myhost = socket.gethostname()

gmi       = ["GPM","GMI","1C","1C","V05"]
amsr2     = ["GCOMW1","AMSR2","1C","1C","V05"]
ssmis_f16 = ["F16","SSMIS","1C","1C","V05"]
ssmis_f17 = ["F17","SSMIS","1C","1C","V05"]
ssmis_f18 = ["F18","SSMIS","1C","1C","V05"]
atms_npp  = ["NPP","ATMS","1C","1C","V05"]
atms_noaa20= ["NOAA20","ATMS","1C","1C","V05"]

mhs_metopa= ["METOPA","MHS","1C","1C","V05"]
mhs_metopb= ["METOPB","MHS","1C","1C","V05"]
mhs_noaa18= ["NOAA18","MHS","1C","1C","V05"]
mhs_noaa19= ["NOAA19","MHS","1C","1C","V05"]

#lspec = [amsr2, ssmis_f16, ssmis_f17, ssmis_f18, atms_npp, atms_noaa20, mhs_metopa, mhs_metopb, mhs_noaa18, mhs_noaa19]

lspec = [amsr2]

#*****************************************
# Start sate,sensor loop
#-----------------------------------------
for spec in lspec:
    sate      = spec[0]
    sensor    = spec[1]
    prdName   = spec[2]
    prj       = spec[3]
    ver       = spec[4]

    if myhost == 'shui':
        epcbaseDir   = '/tank/utsumi/PMM/retepc/%s/%s.%s'%(expr,sensor,sate)
        gprbaseDir   = '/work/hk02/PMM/NASA/%s.%s/2A-CLIM/%s'%(sate,sensor,ver)
        tbbaseDir    = '/work/hk02/PMM/NASA/%s.%s/1C/%s'%(sate,sensor,ver)
        workbaseDir  = '/work'
        mrmsDir      = '/work/hk02/PMM/MRMS/match-GMI-orbit'
        figbaseDir   = '/home/utsumi/temp/multi'

    elif myhost == 'well':
        epcbaseDir   = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/%s/%s.%s'%(expr,sensor,sate)
        gprbaseDir   = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/%s.%s/2A-CLIM/%s'%(sate,sensor,ver)
        tbbaseDir    = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/%s.%s/1C/%s'%(sate,sensor,ver)
        workbaseDir  = '/home/utsumi/mnt/lab_work'
        mrmsDir      = '/home/utsumi/mnt/lab_work/hk02/PMM/MRMS/match-GMI-orbit'
        figbaseDir   = '/home/utsumi/temp/multi'

    else:
        print 'check hostname',myhost
        sys.exit()

    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]
        epcsearch = epcbaseDir + '/%04d/%02d/%02d/nsurfNScmb.??????.*.npy'%(Year,Mon,Day)
        lepcpath  = glob.glob(epcsearch)
        print lepcpath

        for epcpath in lepcpath:
            epcdir = os.path.dirname(epcpath) 
            oid = int(os.path.basename(epcpath).split('.')[1])
            iey = os.path.basename(epcpath).split('.')[2]
            iy  = int(iey.split('-')[0][1:])
            ey  = int(iey.split('-')[1])

            stamp = '.'.join(os.path.basename(epcpath).split('.')[-3:-1])
            latpath   = epcdir + '/lat.%06d.%s.npy'%(oid,stamp)
            lonpath   = epcdir + '/lon.%06d.%s.npy'%(oid,stamp)

            a2precepc = np.load(epcpath)
            a2latepc  = np.load(latpath)
            a2lonepc  = np.load(lonpath)

            #-- GPROF --
            gprDir   = gprbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
            ssearch  = gprDir + '/2A-CLIM.%s.%s.GPROF*.%06d.????.HDF5'%(sate,sensor,oid)
            print ssearch
            gprPath  = glob.glob(ssearch)[0]
            print 'Read GPROF'
            with h5py.File(gprPath) as h:
                a2precgpr = h['S1/surfacePrecipitation'][iy:ey+1]
                a2latgpr  = h['S1/Latitude'][iy:ey+1]
                a2longpr  = h['S1/Longitude'][iy:ey+1]

            ##-- MRMS --
            #ssearch  = mrmsDir + '/GMI.MRMS.130W_55W_20N_55N.%04d%02d%02d.%06d.*.npy'%(Year,Mon,Day,oid)
            #print ssearch
            #mrmsPath = glob.glob(ssearch)[0]
            #print 'Read MRMS'
            #iymr,eymr = map(int, mrmsPath.split('.')[-2].split('-'))
            #a2mrms  = np.load(mrmsPath)
            #a2latmr = a2latgpOrg[iymr:eymr+1]
            #a2lonmr = a2longpOrg[iymr:eymr+1]


            #********************************
            #-- Draw figure ---
            print 'Draw figure'
            fig   = plt.figure(figsize=(6,6))
            vmin,vmax = thpr,10

            #-- My retrieval --
            for i in range(2):
                if i==1:
                    ax = fig.add_axes([0.08,0.05,0.8,0.35])
                    a2dat = ma.masked_less_equal(a2precgpr,thpr)
                    a2lat = a2latgpr
                    a2lon = a2longpr
                    stype = 'GPROF' + ' (>%.2fmm/h)'%(thpr)

                elif i==0:
                    ax = fig.add_axes([0.08,0.55,0.8,0.35])
                    a2dat = ma.masked_less_equal(a2precepc,thpr)
                    a2lat = a2latepc
                    a2lon = a2lonepc
                    stype = 'EPC' + ' (>%.2fmm/h)'%(thpr)

                ssize = 2
                mycm  = 'gist_ncar_r' 
                M     = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)
                im    = M.scatter(a2lon, a2lat, c=a2dat, cmap=mycm, s=ssize, vmin=vmin, vmax=vmax)

                M.bluemarble()
                M.drawcoastlines(color='0.8')

                #---swath boundary lines ---
                a1latTmp= a2lat[:,0]
                a1lonTmp= a2lon[:,0]
                M.plot(a1lonTmp, a1latTmp, '-', color='w',linewidth=0.5)

                a1latTmp= a2lat[:,-1]
                a1lonTmp= a2lon[:,-1]
                M.plot(a1lonTmp, a1latTmp, '-', color='w',linewidth=0.5)

                a1latTmp= a2lat[0,:]
                a1lonTmp= a2lon[0,:]
                M.plot(a1lonTmp, a1latTmp, '-', color='w',linewidth=0.5)

                a1latTmp= a2lat[-1,:]
                a1lonTmp= a2lon[-1,:]
                M.plot(a1lonTmp, a1latTmp, '-', color='w',linewidth=0.5)
                #---------------------------

                dgrid      = 10
                parallels  = arange(-90,90, dgrid)
                meridians  = arange(-180,180,dgrid)
                M.drawparallels(parallels, labels=[1,0,0,0], fontsize=12, linewidth=0.5, fmt='%d', color='0.8')
                M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=12, linewidth=0.5, fmt='%d', color='0.8')

                plt.title(stype, fontsize=14)



            #-- Colorbar (Shared) ----
            cax = fig.add_axes([0.88, 0.2, 0.02, 0.6])
            cb  = plt.colorbar(im, orientation='vertical', cax=cax)
            cb.ax.tick_params(labelsize=14)
            #-- Suptitle -------------
            ssuptitle = '%04d/%02d/%02d #%06d (%s)'%(Year,Mon,Day,oid,expr)
            plt.suptitle(ssuptitle, fontsize=12)

            ##------------
            figDir = figbaseDir + '/%s.%s/%04d/%02d/%02d'%(sensor,sate,Year,Mon,Day)
            util.mk_dir(figDir)
            outPath  = figDir + '/prcp.map.US.%s.%06d.y%d-%d.png'%(expr,oid,iy,ey)
            plt.savefig(outPath)
            print outPath