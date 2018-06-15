import matplotlib
matplotlib.use('Agg')
from numpy import *
import numpy as np
import myfunc.IO.GPM as GPM
from gv_fsub import *
import myfunc.util as util
import GPMGV
from collections import deque
import sys, os
import matplotlib.pyplot as plt

prdName = 'L2A25'


gvPath0    = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.L2A25/FLORIDA-SFL-N/201404/gvprcp.npy'
esurfPath0 = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.L2A25/FLORIDA-SFL-N/201404/eSurf.npy'

gvPath1    = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.L2A25/N.Carolina-IPHEx_NASA/201404/gvprcp.npy'
esurfPath1 = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.L2A25/N.Carolina-IPHEx_NASA/201404/eSurf.npy'

gvprcp0 = np.load(gvPath0)
eSurf0  = np.load(esurfPath0)

gvprcp1 = np.load(gvPath1)
eSurf1  = np.load(esurfPath1)


domain = 'FLORIDA-SFL-N' 
figDir = '/work/a01/utsumi/GPMGV/fig'
figPath= figDir + '/temp.simple.plot.%s.%s.png'%(prdName,domain)
#plt.plot(gvprcp0,eSurf0,'o', color='k')
plt.loglog(gvprcp0,eSurf0,'o', color='k')
plt.title('%s %s'%(prdName, domain))

vmax  = 100
vmin  = 0.1
plt.ylim([vmin,vmax])
plt.xlim([vmin,vmax])

#-- 1:1 line -
plt.plot([vmin,vmax],[vmin,vmax],'-',color='k')

plt.savefig(figPath)
plt.clf()
print figPath



domain = 'N.Carolina-IPHEx_NASA' 
figDir = '/work/a01/utsumi/GPMGV/fig'
figPath= figDir + '/temp.simple.plot.%s.%s.png'%(prdName,domain)
#plt.plot(gvprcp1,eSurf1,'o', color='k')
plt.loglog(gvprcp1,eSurf1,'o', color='k')
plt.title('%s %s'%(prdName, domain))

vmax  = 100
vmin  = 0.1
plt.ylim([vmin,vmax])
plt.xlim([vmin,vmax])

#-- 1:1 line -
plt.plot([vmin,vmax],[vmin,vmax],'-',color='k')

plt.savefig(figPath)
plt.clf()
print figPath




# lat & lon
latPath0    = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.L2A25/FLORIDA-SFL-N/201404/sateLat.npy'
lonPath0 = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.L2A25/FLORIDA-SFL-N/201404/sateLon.npy'

latPath1    = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.L2A25/N.Carolina-IPHEx_NASA/201404/sateLat.npy'
lonPath1 = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.L2A25/N.Carolina-IPHEx_NASA/201404/sateLon.npy'

lat0 = np.load(latPath0)
lon0  = np.load(lonPath0)
lat1 = np.load(latPath1)
lon1  = np.load(lonPath1)



domain = 'FLORIDA-SFL-N' 
figDir = '/work/a01/utsumi/GPMGV/fig'
figPath= figDir + '/temp.simple.latlon.%s.%s.png'%(prdName,domain)
plt.plot(lon0,lat0,'o', color='k')
plt.title('%s %s'%(prdName, domain))

#plt.ylim([vmin,vmax])
#plt.xlim([vmin,vmax])

#-- 1:1 line -
plt.plot([vmin,vmax],[vmin,vmax],'-',color='k')

plt.savefig(figPath)
plt.clf()
print figPath



domain = 'N.Carolina-IPHEx_NASA' 
figDir = '/work/a01/utsumi/GPMGV/fig'
figPath= figDir + '/temp.simple.latlon.%s.%s.png'%(prdName,domain)
plt.plot(lon1,lat1,'o', color='k')
plt.title('%s %s'%(prdName, domain))

#plt.ylim([vmin,vmax])
#plt.xlim([vmin,vmax])

#-- 1:1 line -
plt.plot([vmin,vmax],[vmin,vmax],'-',color='k')

plt.savefig(figPath)
plt.clf()
print figPath











