import numpy as np
import pylab as pl
#import matplotlib.gridspec as gridspec
from glob import glob
import numpy.ma as ma
import sys,os, glob
from datetime import datetime, timedelta
import myfunc.util as util
from collections import deque
import socket
import shutil

ldy   = [-1,0,1]
ldx   = [-3,-2,-1,0,1,2,3]
#ldx   = [-2,-1,0,1,2]
ldydx = [[dy,dx] for dy in ldy for dx in ldx]
ntc1  = 9
ntc2  = 4
ncomb = (ntc1+ntc2)* len(ldydx)
imid  = int((len(ldy)*len(ldx)-1)/2)
lisurf = range(1,14)

nsample = 10000
expr = 'stop-ml.smp%d'%(nsample)
stopexpr = 'wtpdf01'
stopact  = 'LTQZ'

iDTime = datetime(2014,8,1)
eDTime = datetime(2014,8,31)
#iDTime = datetime(2014,12,1)
#eDTime = datetime(2014,12,31)

lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

myhost = socket.gethostname()
if myhost =='shui':
    pmmbaseDir ='/tank/utsumi/PMM'
    hdfbaseDir ='/work/hk01/PMM/NASA'
    stopbaseDir= '/tank/utsumi/PMM/stop'

    figDir = '/home.rainbow/utsumi/public_html/tempfig/stop'
elif myhost =='well':
    pmmbaseDir ='/home/utsumi/mnt/lab_tank/utsumi/PMM'
    hdfbaseDir ='/home/utsumi/mnt/lab_work/hk01/PMM/NASA'
    stopbaseDir='/home/utsumi/mnt/lab_tank/utsumi/PMM/stop'

else:
    print 'check myhost',myhost


tcmin, tcmax = 50, 350
dvarminmax = {'T':[200,350],'Q':[0,120],'Z':[-100,8800]}

#get_ipython().magic(u'matplotlib inline')
#******************************************
# Functions
#******************************************
def shift_array(ain=None, dy=None,dx=None,miss=-9999):
    ny,nx,nz = ain.shape
    aout = np.ones([ny,nx,nz]).astype(ain.dtype)*miss
    if   dy<=0: iy0=0; ey0=ny-abs(dy); iy1=abs(dy); ey1=ny
    elif dy> 0: iy0=abs(dy); ey0=ny; iy1=0; ey1=ny-abs(dy)
    if   dx<=0: ix0=0; ex0=nx-abs(dx); ix1=abs(dx); ex1=nx
    elif dx> 0: ix0=abs(dx); ex0=nx; ix1=0; ex1=nx-abs(dx)

    aout[iy0:ey0,ix0:ex0] = ain[iy1:ey1,ix1:ex1]
    return aout


def my_unit(x,Min,Max):
    return (x-Min)/(Max-Min)

#******************************************
outbaseDir= stopbaseDir + '/orbit/%s-%s-ssn%s'%(stopexpr,stopact, 0)
#******************************************
for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]
    #** Read HDF GMI Tb ****
    gmiDir = hdfbaseDir + '/GPM.GMI/1C/V05/%04d/%02d/%02d'%(Year,Mon,Day)

    ssearch = gmiDir + '/1C.GPM.GMI.*.??????.????.HDF5'

    lgmiPath= glob.glob(ssearch)

    for gmiPath in lgmiPath:
        oid = int(gmiPath.split('.')[-3])

        inDir = outbaseDir + '/%04d/%02d'%(Year,Mon)
        inPath= inDir + '/stop.%06d.npy'%(oid)
        if not os.path.exists(inPath): continue
        print ''
        print ''
        print inPath

        outDir = outbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
        outPath= outDir + '/stop.%06d.npy'%(oid)
        util.mk_dir(outDir)
        shutil.move(inPath, outPath)
        print outPath

        

