import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py
import myfunc.util as util
import numpy as np
from numpy import *

srcDir = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/GPM.DPRGMI/2B/V06/2014/10/14'
srcPath= srcDir + '/2B.GPM.DPRGMI.CORRA2018.20141014-S050829-E064102.003556.V06A.HDF5'

#---------------------------------------------
def ave_9grids_3d(a3in, a1y, a1x, miss):
    '''
    returns 2-d array with the size of (nl,nz)
    a3in: (ny,nx,nz)
    nl = len(a1y)=len(a1x)
    output: (nl, nz)

    copied from mk.epcdb.dprvar.py, 2020/01/28
    '''

    if ma.is_masked(a3in):
        a3in = a3in.filled(miss)  # 2019/12/02
    #-- Average 9 grids (over Linearlized Z)--
    nydpr,nxdpr,nzdpr= a3in.shape
    ldydx = [[dy,dx] for dy in [-1,0,1] for dx in [-1,0,1]]

    a1dprmask   = False

    a3datTmp    = empty([9,len(a1y),nzdpr], float32)

    for itmp, [dy,dx] in enumerate(ldydx):
        a1yTmp  = ma.masked_outside(a1y + dy, 0,nydpr-1)
        a1xTmp  = ma.masked_outside(a1x + dx, 0,nxdpr-1)
        a1dprmask= a1dprmask + a1yTmp.mask + a1xTmp.mask

        a2datTmp= a3in[a1yTmp.filled(0),a1xTmp.filled(0),:]

        a3datTmp[itmp,:] = a2datTmp


    a2datTmp = ma.masked_equal(a3datTmp,miss).mean(axis=0)
    a2datTmp[a1dprmask,:] = miss
    return a2datTmp


#---------------------------------------------

# SE.US case, oid=003556, 2014/10/14
iy, ey = 987, 1047
#iy,ey = 917,1117
#iy,ey = 917,1117
xpos    = 100  # x-position for cross section

#--- Xdpr Ydpr --
Year,Mon,Day = 2014,10,14
oid = 3556
srcDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.GMI.V05A/S1.ABp083-137.Ku.V06A.IDX/%04d/%02d/%02d'%(Year,Mon,Day)
xPath= srcDir + '/Xpy.1.%06d.npy'%(oid)
yPath= srcDir + '/Ypy.1.%06d.npy'%(oid)

##-- extract center and domain --
a1dprx = np.load(xPath)[iy:ey+1,xpos-83]
a1dpry = np.load(yPath)[iy:ey+1,xpos-83]


with h5py.File(srcPath,'r') as h:
    a3dat = h['NS/precipTotWaterCont'][:]

print a3dat.shape
a2dat = a3dat[a1dpry, a1dprx,::-1]  # top to bottom --> bottom to top
a2dat = ma.masked_less(a2dat,0)

a1alt = np.arange(88)*0.25
for i in range(a2dat.shape[0]):
    plt.plot(a2dat[i], a1alt, '-')

#-- Single average line ---
a1dat = ma.masked_less(a2dat,0).mean(axis=0)
plt.plot(a1dat, a1alt, '-', color='k', linewidth=3)
#--------------------------

plt.ylim([0,10])
plt.xlim([0,6])
plt.title('raw %06d'%(oid))
print a2dat 
print a2dat.max()

figPath = '/home/utsumi/temp/ret/temp.plot.prof-raw.%06d.png'%(oid)
plt.savefig(figPath)
plt.clf()
print figPath

#****************************************
# 3x3 average
#----------------------------------------
print ''
print ''
a2ave = ave_9grids_3d(a3dat, a1dpry, a1dprx, miss=-9999.9)
a2ave = a2ave[:,::-1]  # top to bottom -> bottom to top
print a2ave
a1alt = np.arange(88)*0.25
for i in range(a2ave.shape[0]):
    plt.plot(a2ave[i], a1alt, '-')

#-- Single average line ---
a1ave = ma.masked_less(a2ave,0).mean(axis=0)
plt.plot(a1ave, a1alt, '-', color='k', linewidth=3)
#--------------------------



plt.ylim([0,10])
plt.xlim([0,6])
plt.title('ave %06d'%(oid))
figPath = '/home/utsumi/temp/ret/temp.plot.prof-ave.%06d.png'%(oid)
plt.savefig(figPath)
plt.clf()
print figPath


