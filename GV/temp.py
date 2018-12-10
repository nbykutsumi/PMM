from numpy import *
import numpy as np
srcDir = '/media/disk2/share/GPMGV/DOM.L2A25/cbin.11/TEXAS-HAR/200109'
elevPath = srcDir + '/elev.npy'
latPath  = srcDir + '/sateLat.npy'
lonPath  = srcDir + '/sateLon.npy'

a1elev = np.load(elevPath)
a1lat  = np.load(latPath)
a1lon  = np.load(lonPath)


for i,elev in enumerate(a1elev):
    lat = a1lat[i]
    lon = a1lon[i]

    print lat,lon,elev
