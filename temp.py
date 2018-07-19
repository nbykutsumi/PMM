from numpy import *
import h5py

srcDir  = '/data4/common/GPM/GPM.DPR/L2.DPR/05/2014/05'
srcPath = srcDir + '/GPMCOR_DPR_1405310502_0634_001440_L2S_DD2_05A.h5'

varName = 'NS/DSD/phase'
h5 = h5py.File(srcPath, 'r')
h5Var = h5[varName]
aOut  = h5Var[:]
a = ma.masked_equal(aOut, 255).compressed()
print a
print a.min(), a.max()
