from numpy import *
import numpy as np
import myfunc.IO.GPM as GPM
from gv_fsub import *

srcDir = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.2A-CLIM/FLORIDA-SFL-N/200905'

dtime   = np.load(srcDir + '/dtime.npy')
sateLat = np.load(srcDir + '/sateLat.npy')
sateLon = np.load(srcDir + '/sateLon.npy')
gLat    = np.load(srcDir + '/gLat.npy')
gLon    = np.load(srcDir + '/gLon.npy')
gName   = np.load(srcDir + '/gName.npy')
gvprcp  = np.load(srcDir + '/gvprcp.npy')
eSurf   = np.load(srcDir + '/eSurf.npy')

dt = 0
for i in range(len(dtime)):
    if not ((gvprcp[i][dt+15]==0)&(eSurf[i]==0)):
        print 'dt=',dt,dtime[i],gName[i],sateLat[i],sateLon[i],gvprcp[i][dt+15],eSurf[i]

'''
dt= 0 2009-05-30 22:17:00 0134 25.3298 -80.5325 137.13 1.90017
dt= 2 2009-05-30 22:17:00 0134 25.3298 -80.5325 182.84 1.90017

'''
#-- load satellite ---
sensor   = 'TRMM.TMI'
prdName  = '2A-CLIM'
prdVer   = 'V05'
minorVer = 'A'
gpm = GPM.L2A_GPROF_HDF5(sensor=sensor, prdName=prdName, version=prdVer, minorversion=minorVer)

graDir  = '/home/utsumi/mnt/wellshare/data/GPM/TRMM.TMI/2A-CLIM/V05/2009/05'
graPath = graDir + '/2A-CLIM.TRMM.TMI.GPROF2017v2.20090530-S214139-E231401.065742.V05A.HDF5'

satedat = gpm.load_var_granule(graPath, 'S1/surfacePrecipitation')
satelat = gpm.load_var_granule(graPath, 'S1/Latitude')
satelon = gpm.load_var_granule(graPath, 'S1/Longitude')
satetime= gpm.load_dtime_granule(graPath)

glat = 25.3298
glon = -80.5325
thdist= 2.5 # [km]
a1xsate, a1ysate = gv_fsub.gauge_match_pyxy(satelon.T, satelat.T, glon, glat, thdist)

a1xsate = ma.masked_less(a1xsate,0).compressed()
a1ysate = ma.masked_less(a1ysate,0).compressed()

print a1xsate, a1ysate
slat = satelat[a1ysate,a1xsate]
slon = satelon[a1ysate,a1xsate]
sdat = satedat[a1ysate,a1xsate]
stime= satetime[a1ysate]
print slat, slon, sdat, stime

