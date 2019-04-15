import os, sys
import glob
import subprocess
from datetime import datetime, timedelta
import myfunc.util as util
import socket

#- compile --------
scmd = 'python compile_ret.py jplret-test.c'
lcmd = scmd.split()
subprocess.call(lcmd)
#------------------

progDir  = '/home/utsumi/bin/PMM/ret-epc'
prog    = progDir + '/jplret-test'

#srcPath = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/2017/01/01/1C.GPM.GMI.XCAL2016-C.20170101-S020858-E034131.016156.V05A.HDF5'
#srcPath = '/home/utsumi/temp/temp.ext.016156.HDF5'
srcPath = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/2014/08/02/1C.GPM.GMI.XCAL2016-C.20140802-S062222-E075455.002421.V05A.HDF5'
modelPath='/home/utsumi/bin/JPLCODE/EPC_ret_20190221/1C.GPM.GMI.XCAL2016-C.20140802-S062222-E075455.002421.V05A.HDF5.MERRA2.nc'

outDir  = '/home/utsumi/temp/EPC'
util.mk_dir(outDir)
clat, clon = 14, 2
#clat,clon = -8.861334, -15.027719
#lcmd = map(str, [prog, srcPath, clat, clon, outDir])
lcmd = map(str, [prog, srcPath, clat, clon, outDir, modelPath])
print lcmd 
#p = subprocess.Popen(lcmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#p = subprocess.Popen(lcmd)
#p.wait()
#p = subprocess.call(lcmd)
print ' '.join(lcmd)
#subprocess.call(lcmd)
