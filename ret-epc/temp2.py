from numpy import *
import subprocess
import netCDF4

'''
db_idx = 2601
#dbPath = '/work/hk01/utsumi/JPLDB/EPC_DB/GMI_EPC_DATABASE/db_02601.bin'
dbPath = '/work/hk01/utsumi/JPLDB/EPC_DB/GMI_EPC_DATABASE/db_%05d.bin'%(db_idx)

prog = '/home/utsumi/bin/JPLCODE/EPC_ret_20190221/EPC_DB_conv2netcdf'
soutPath = '/home/utsumi/temp/out/out.nc'
scmd = '%s %s %s'%(prog, dbPath, soutPath)
lcmd = scmd.split()
print lcmd
subprocess.call(lcmd)
'''

srcDir = '/home/utsumi/bin/PMM/ret-epc'
srcPath= srcDir + '/test_netcdf.nc'
nc = netCDF4.Dataset(srcPath,'r')
print nc

a3prof = nc.groups['DPR'].variables['precip_prof_MS'][:]
print a3prof.min(), a3prof.max()
print a3prof.shape

