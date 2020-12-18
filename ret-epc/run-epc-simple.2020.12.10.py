import subprocess
import numpy as np
##*****************************
## Make GMI S2 match-up location data
##*****************************
#progtmp = './mk.match.idx.gmiS2.gmi.fullswath.py'
#
#dtmp = {}
#dtmp['Year']    = 2014
#dtmp['Mon']     = 6
#dtmp['Day']     = 1
#dtmp['pmwPath'] = './1C.GPM.GMI.XCAL2016-C.20140601-S010534-E023807.001453.V05A.HDF5'
#dtmp['obaseDir']= '.'
#
#dydxDir = '.'
#dtmp['dyPath0'  ]= dydxDir + '/dy.000.npy'
#dtmp['dxPath0'  ]= dydxDir + '/dx.000.npy'
#dtmp['dyPath180']= dydxDir + '/dy.180.npy'
#dtmp['dxPath180']= dydxDir + '/dx.180.npy'
#
#sargv = ['%s=%s'%(key, dtmp[key]) for key in list(dtmp.keys())]
#sargv = ' '.join(sargv)
#lcmd = ['python', progtmp, sargv]
#print(lcmd)
#returncode = subprocess.call(lcmd)

#*****************************
# Retrieval
#*****************************
dargv = {}

dargv['dbtype']     = 'my'
dargv['sate']       = 'GPM'
dargv['sensor']     = 'GMI'

sensor = dargv['sensor']
dargv['nscan'] = {'GMI':2, 'AMSR2':5, 'SSMIS':4, 'ATMS':4, 'MHS':1}[sensor]   # number of scan classes (S1, S2, ..)
dargv['mainscan'] = {'GMI':1, 'AMSR2':1, 'SSMIS':1, 'ATMS':1, 'MHS':1, 'SAPHIR':1}[sensor]  # Precipitation variables are estimated at this scan class footprints
dargv['NEM'] =  {'GMI':12, 'AMSR2':13, 'SSMIS':10, 'ATMS':6, 'MHS':4, 'SAPHIR':4}[sensor]
dargv['NTBREG'] = {'GMI':13, 'AMSR2':10, 'SSMIS':11, 'ATMS':9, 'MHS':5, 'SAPHIR':6}[sensor]
dargv['NEM_USE'] = 3
dargv['NPCHIST'] = 29
dargv['NLEV_DPR'] = 50    # extract this number of layers
dargv['NLEV_PRECIP'] = 50
dargv['thwtmin'] = 0.1
dargv['miss'] = -9999.
dargv['miss_int32'] = np.int32(-9999)
dargv['miss_int16'] = np.int16(-9999)
#------------
dargv['DB_MAXREC'] = 10000
dargv['DB_MINREC'] = 1000
dargv['DB_USE_MINREC'] = 2
dargv['NDB_EXPAND'] = 30
dargv['DB_RAINFRAC'] = 0.0001 # minimum fraction of precipitating events (>=1mm/h) in the DB required for retrieval

dargv['MAX_INC_DIFF'] = 20  # degrees
dargv['MAX_T2M_DIFF'] = 10
dargv['MAX_TQV_DIFF'] = 10
dargv['MAX_STOP_DIFF'] = 2000  # [meter]
dargv['MAX_TB_RMSD'] = -9999 # max TB RMS difference between observed and DB pixel
dargv['MAX_RMA_0'] = -9999.  # Maximum Ratio of missing amount (0>=mm/h) acceptable for rain / no-rain classification # -9999. --> No screening

dargv['MIN_RNR'] = 0.1  # [mm/h] for Rain/No-rain based on first guess precipitation
dargv['STD_STORMTOP'] = 2100 # [m] stop

dargv['flag_crosstrack'] = {'GMI':0, 'AMSR2':0, 'SSMIS':0, 'ATMS':1, 'MHS':1, 'SAPHIR':1}[sensor]
dargv['flag_top_var'] = 1  # 0: No top-ranked vars. 1: Retrieve top-ranked vars.
dargv['flag_rel_surf'] = 1  # 0: Not relative to surface 1: Profiles that are relative to surface
dargv['flag_out_conv'] = 1  # 0: Not retrieve convective

#------------
dargv['coefDir']    = '/media/disk2/share/PMM/EPCDB/EPC_COEF/%s'%(sensor)
dargv['dbDir']      = '/media/disk2/share/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(dargv['DB_MAXREC'])
dargv['relprofDir'] = ''

#dargv['oid'] = 1453
dargv['oid'] = 0000
dargv['clat'] = -9999   # center latitude. Set to -9999 for global
dargv['clon'] = -9999 # center longitude. Set to -9999 for global
dargv['dlatlon'] = -9999 # parameter for center location searching (degree). Set to -9999 for global
dargv['iscan'] = -9999
dargv['escan'] = -9999
dargv['dscan'] = -9999  # dscan*2 = number of retrieved scans. Set to -9999 for global


#dargv['srcPath'] = '../data/1C.GPM.GMI.XCAL2016-C.20140601-S010534-E023807.001453.V05A.HDF5'  # L1C TB file
#dargv['s2xPath'] = '../data/Xpy.1.001453.npy'
#dargv['s2yPath'] = '../data/Ypy.1.001453.npy'
dargv['srcPath'] = '../data2/tb.HDF5'  # TB file
dargv['s2xPath'] = '../data2/Xpy.npy'
dargv['s2yPath'] = '../data2/Ypy.npy'

dargv['t2mPath'] = ''
dargv['tqvPath'] = ''
dargv['elevPath']= ''
dargv['stopPath']= ''
dargv['rnrPath'] = ''
dargv['outDir']  = '../out'   # output directory
#-------------
sargv = ['%s=%s'%(key, dargv[key]) for key in list(dargv.keys())]
sargv = ' '.join(sargv)

prog = './ret-epc-multi.2020.12.10.py'
lcmd = ['python', prog, sargv]
print(lcmd)
returncode = subprocess.call(lcmd)
if returncode==1:
    print('Error for',lcmd)


