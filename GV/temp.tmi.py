from numpy import *
import numpy as np
import myfunc.IO.GPM as GPM
from gv_fsub import *
import pickle

srcDir = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.2A-CLIM/FLORIDA-SFL-N/200905'

dtime   = np.load(srcDir + '/p_dtime.npy')
sateLat = np.load(srcDir + '/p_sateLat.npy')
sateLon = np.load(srcDir + '/p_sateLon.npy')
gvprcp  = np.load(srcDir + '/p_gvprcp.npy')
eSurf   = np.load(srcDir + '/p_eSurf.npy')
ngv     = np.load(srcDir + '/p_ngv.npy')
gbin    = np.load(srcDir + '/p_groundBin.npy')
eSurf0  = np.load(srcDir + '/eSurf.npy')

with open(srcDir + '/p_gLat.pickle', 'r') as f:
    gLat = pickle.load(f)
with open(srcDir + '/p_gLon.pickle', 'r') as f:
    gLon = pickle.load(f)
with open(srcDir + '/p_gName.pickle', 'r') as f:
    gName= pickle.load(f)


print gvprcp.shape
print dtime.shape
print sateLat.shape
print eSurf.shape
print eSurf0.shape
print sum(ngv)
dt = 0
for i in range(len(dtime)):
    if not ((gvprcp[i][dt+15]==0)&(eSurf[i]==0)):
        print 'dt=',dt,dtime[i],ngv[i],gName[i],gbin[i],sateLat[i],sateLon[i],gvprcp[i][dt+15],eSurf[i]

'''
dt= 0 2009-05-26 01:18:00 ['0063', '0144'] 26.6946 -80.8204 0.555 3.19759
0063 2009 05 26 146 01 19 00    1.11
0144 0
'''

'''
#-- load satellite ---
sensor   = 'TRMM.TMI'
prdName  = '2A-CLIM'
prdVer   = 'V05'
minorVer = 'A'
gpm = GPM.L2A_GPROF_HDF5(sensor=sensor, prdName=prdName, version=prdVer, minorversion=minorVer)

graDir  = '/home/utsumi/mnt/wellshare/data/GPM/TRMM.TMI/2A-CLIM/V05/2009/04'
graPath = graDir + '/2A-CLIM.TRMM.TMI.GPROF2017v2.20090414-S213805-E231028.065025.V05A.HDF5'

satedat = gpm.load_var_granule(graPath, 'S1/surfacePrecipitation')
satelat = gpm.load_var_granule(graPath, 'S1/Latitude')
satelon = gpm.load_var_granule(graPath, 'S1/Longitude')
satetime= gpm.load_dtime_granule(graPath)

glat = 27.6432
glon = -81.1205
thdist= 7.5 # [km]
a1xsate, a1ysate = gv_fsub.gauge_match_pyxy(satelon.T, satelat.T, glon, glat, thdist)

a1xsate = ma.masked_less(a1xsate,0).compressed()
a1ysate = ma.masked_less(a1ysate,0).compressed()

print a1xsate, a1ysate
slat = satelat[a1ysate,a1xsate]
slon = satelon[a1ysate,a1xsate]
sdat = satedat[a1ysate,a1xsate]
stime= satetime[a1ysate]
print slat, slon, sdat, stime

'''
