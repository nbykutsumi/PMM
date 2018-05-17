from numpy import *
from gv_fsub import *
import numpy as np

srcDir = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.L2A25/FLORIDA-STJ/201401'

dtimepath   = srcDir + '/dtime.npy'
glatpath    = srcDir + '/gLat.npy'
glonpath    = srcDir + '/gLon.npy'
gnamepath   = srcDir + '/gName.npy'
profpath    = srcDir + '/prof.npy'
gvpath      = srcDir + '/gvprcp.npy'
satelatpath = srcDir + '/sateLat.npy'
satelonpath = srcDir + '/sateLon.npy'


adtime    = np.load(dtimepath)
aglat     = np.load(glatpath)
aglon     = np.load(glonpath)
agname    = np.load(gnamepath)
aprof     = np.load(profpath)
agvprcp   = np.load(gvpath)
asatelat  = np.load(satelatpath)
asatelon  = np.load(satelonpath)


print adtime.shape
print aglat.shape
print agname.shape
print agvprcp.shape
print asatelat.shape

for i in range(len(adtime)):
    print ''
    print '-'*50
    print adtime[i],agname[i],aglat[i], aglon[i], asatelat[i], asatelon[i]
    print agvprcp[i]
    print aprof[i]
