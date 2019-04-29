from numpy import *
import h5py

srcPath = '/work/hk01/PMM/NASA/GPM.Ku/2A/V06/2017/01/01/2A.GPM.Ku.V8-20180723.20170101-S003624-E020857.016155.V06A.HDF5'

with h5py.File(srcPath) as h:
    a3dat = h['/NS/SLV/precipRate'][:,15:49-15+1,:]
    a2surf= h['/NS/SLV/precipRateESurface'][:,15:49-15+1]
print a3dat.shape
a2flag = ma.masked_greater(a3dat,0).mask.any(axis=2)
print a2flag.sum()
print ''

a1flag = a2flag.flatten()
a2dat  = a3dat.reshape(-1,176)[a1flag]
a1surf = a2surf.flatten()[a1flag]
for i in range(100):
    print a1surf[i],'  ', a2dat[i,:]

'''
a3dat  = (ma.masked_less(a3dat,0)*100).astype(int16)
a3dat  = a3dat.filled(-9999)
a2flag = ma.masked_greater(a3dat,0).mask.any(axis=2)
print a3dat.shape
a2flag = ma.masked_greater(a3dat,0).mask.any(axis=2)
print a2flag.sum()
'''
