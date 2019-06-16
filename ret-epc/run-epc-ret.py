import sys, os
import myfunc.util as util
import glob
import subprocess

dargv = {}
#** Constants ******
sensor  = 'GMI'
coefDir = '/work/hk01/utsumi/JPLDB/EPC_COEF/%s'%(sensor)

#dbDir   = '/work/hk01/utsumi/JPLDB/EPC_DB/GMI_EPC_DATABASE'
dbDir   = '/work/hk01/utsumi/JPLDB/EPC_DB/GMI_EPC_DATABASE_TEST29'

dargv['sensor']=sensor
dargv['coefDir']=coefDir
dargv['dbDir']=dbDir
#------------
#oid     = 3556
#clat    = 34    # SE.US case. oid = 003556
#clon    = -86   # 2014/10/14  05:42:03 UTC
#
#dlatlon = 3  # used to search the domain center
#dscan   = 90
##dscan   = 55
##dscan   = 5
Year, Mon, Day= 2014, 10, 14
dargv['oid'] = 3556
dargv['clat'] = 34
dargv['clon'] = -86
dargv['dlatlon'] = 3
dargv['iscan'] = 1010
dargv['escan'] = 1022

dargv['dscan'] = 5
#------------
dargv['NEM'] = 12
dargv['NTBREG'] = 13
dargv['NEM_USE'] = 3
dargv['NPCHIST'] = 29
dargv['NLEV_DPR'] = 88    # extract this number of layers
dargv['NLEV_PRECIP'] = 22
dargv['thwtmin'] = 0.01
dargv['miss'] = -9999.
#------------
dargv['DB_MAXREC'] = 1000
dargv['DB_MINREC'] = 5000
dargv['NDB_EXPAND'] = 20
dargv['DB_RAINFRAC'] = 0.01 # minimum fraction of precipitating events (>=1mm/h) in the DB required for retrieval
dargv['MAX_T2M_DIFF'] = 20

dargv['outDir'] = '/home/utsumi/temp/out'
#------------
oid = dargv['oid']
# 2014/10/14 SE.US GMI=003556 **
dargv['srcPath'] = glob.glob('/work/hk01/PMM/NASA/GPM.GMI/1C/V05/%04d/%02d/%02d/1C.GPM.GMI.*.%06d.????.HDF5'%(Year,Mon,Day,oid))[0]
dargv['srcPath'] = glob.glob('/work/hk01/PMM/NASA/GPM.GMI/1C/V05/%04d/%02d/%02d/1C.GPM.GMI.XCAL2016-C.*.%06d.????.HDF5'%(Year,Mon,Day,oid))[0]
dargv['s2xPath'] = glob.glob('/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d/Xpy.1.%06d.npy'%(Year,Mon,Day,oid))[0]
dargv['s2yPath'] = glob.glob('/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d/Ypy.1.%06d.npy'%(Year,Mon,Day,oid))[0]
dargv['tsPath']  = glob.glob('/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.MERRA2.t2m/%04d/%02d/%02d/t2m.%06d.npy'%(Year,Mon,Day,oid))[0]
dargv['elevPath']= glob.glob('/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.gtopo/%04d/%02d/%02d/gtopo.%06d.npy'%(Year,Mon,Day,oid))[0]

#dargv['srcPath'] = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/2014/10/14/1C.GPM.GMI.XCAL2016-C.20141014-S050829-E064102.003556.V05A.HDF5'
#dargv['srcPath'] = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/2014/10/14/1C.GPM.GMI.XCAL2016-C.20141014-S050829-E064102.003556.V05A.HDF5'
#dargv['s2xPath'] = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/2014/10/14/Xpy.1.003556.npy'
#dargv['s2yPath'] = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/2014/10/14/Ypy.1.003556.npy'
#dargv['tsPath']   = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.MERRA2.t2m/2014/10/14/t2m.003556.npy'
#dargv['elevPath']= '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.gtopo/2014/10/14/gtopo.003556.npy'


#-------------
sargv = ['%s=%s'%(key, dargv[key]) for key in dargv.keys()]
sargv = ' '.join(sargv)


prog = 'ret-epc-29bins.py'
lcmd = ['python', prog, sargv]
print lcmd
subprocess.call(lcmd)
