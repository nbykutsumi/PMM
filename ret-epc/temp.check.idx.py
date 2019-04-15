import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import numpy as np
import h5py
import epcfunc

clat    = 14
clon    = 2  # -180 - +180

dlatlon = 6
#dscan   = 55
dscan   = 25

srcPath0 = '/home/utsumi/temp/out/idx-jpl-full-002421.npy'
srcPath1 = 'GPM_EPC_002421_20140802_0726.NS_MS.nc'

a2idx0 = np.load(srcPath0)[2034:2084+1]
#*****************
#- Read JPL data ----
with h5py.File(srcPath1) as h:
    a2esurfjpl = h['MS/precip'][:]
    a2dpr      = h['DPR/precip_NS'][:]
    a2latjpl   = h['latitude'][:]
    a2lonjpl   = h['longitude'][:]
    a2idxjpl   = h['db_index'][:]

a2idxjpl, a2latjpl, a2lonjpl = epcfunc.extract_domain_2D(a2idxjpl, a2latjpl, a2lonjpl, clat, clon, dlatlon, dscan)

print a2idx0.shape, a2idxjpl.shape

plt.scatter(a2idxjpl, a2idx0)
plt.savefig('/home/utsumi/temp/out/temp.scatter.idx.png')
