# %%
import matplotlib
matplotlib.use('Agg')
#%matplotlib inline
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from datetime import datetime, timedelta
import myfunc.util as util
import glob
import h5py
from numpy import *

iDTime  = datetime(2018,1,1)
eDTime  = datetime(2018,1,1)
dDTime  = timedelta(days=1)
lDTime  = util.ret_lDTime(iDTime, eDTime, dDTime)

gpr_amsr2      = ["GCOMW1","AMSR2"]
gpr_ssmis_f16  = ["F16","SSMIS"]
gpr_atms_noaa20= ["NOAA20","ATMS"]
amsr2     = ["GCOMW1","AMSR2"]
ssmis_f16 = ["F16","SSMIS"]
ssmis_f17 = ["F17","SSMIS"]
ssmis_f18 = ["F18","SSMIS"]
atms_npp  = ["NPP","ATMS"]
atms_noaa20= ["NOAA20","ATMS"]
mhs_metopa= ["METOPA","MHS"]
mhs_metopb= ["METOPB","MHS"]

lspec = [atms_noaa20]

#[[lllat,lllon],[urlat,urlon]] = [[-90,-180],[90,180]]
[[lllat,lllon],[urlat,urlon]] = [[20,-130],[55,-55]]  # MRMS
for (sate,sensor) in lspec:
    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]
        srcbasedir   = "/home/utsumi/mnt/lab_work/hk02/PMM/NASA"
        srcdir = srcbasedir + '/%s.%s/2A-CLIM/V05/%04d/%02d/%02d'%(sate,sensor,Year,Mon,Day)
        search= srcdir + '/2A-CLIM.%s.%s.*.HDF5'%(sate,sensor)
        lsrcpath = glob.glob(search)
        lsrcpath = sort(lsrcpath)
        print lsrcpath
        print search
        for srcpath in lsrcpath:
            with h5py.File(srcpath,'r') as h:
                a2prec = h['S1/surfacePrecipitation'][:]
                a2lat  = h['S1/Latitude'][:]
                a2lon  = h['S1/Longitude'][:]

            oid = int(srcpath.split('.')[-3])

            #--- test ---
            if oid ==623:
                iscan,escan = 1470, 1569
                a2prec = a2prec[iscan:escan+1]
                a2lat  = a2lat [iscan:escan+1]
                a2lon  = a2lon [iscan:escan+1]


            #-------------------
            a2dat = ma.masked_less_equal(a2prec,0.1)

            fig = plt.figure(figsize=(10,5)) 
            ax  = fig.add_axes([0.1,0.1,0.8,0.8])
            M   = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)

            #- test ---
            #M.plot(-173.4, 8, 'x', color='k')
            #----------

            im = M.scatter(a2lon, a2lat, c=a2dat, cmap='jet', s=1, vmin=0, vmax=10)

            dgrid = 2
            parallels  = arange(-90,90, dgrid)
            meridians  = arange(-180,180,dgrid)
            M.drawparallels(parallels, labels=[1,0,0,0], fontsize=12, linewidth=0.5, fmt='%d', color='0.8')
            M.drawmeridians(meridians, labels=[0,0,0,1], fontsize=12, linewidth=0.5, fmt='%d', color='0.8', rotation=30)

            M.drawcoastlines()

            stitle = '%s.%s %04d/%02d/%02d #%06d'%(sate,sensor,Year,Mon,Day,oid)
            plt.title(stitle, fontsize=14)
            plt.colorbar(im)

            figdir = '/home/utsumi/temp/multi'
            figpath= figdir + '/prec.%s.%s.%04d.%02d.%02d.id.%06d.png'%(sate,sensor,Year,Mon,Day,oid)
            plt.savefig(figpath)
            plt.show()
            print figpath
    # %%
