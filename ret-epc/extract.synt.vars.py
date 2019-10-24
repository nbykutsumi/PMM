import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#import matplotlib.cm as cm
import numpy as np
from numpy import *
import h5py
from bisect import bisect_left
import epcfunc
import sys, os, glob, socket
import myfunc.util as util
import calendar
import pickle
import EPCDB

calcflag  = True
DB_MAXREC = 10000
DB_MINREC = 1000
nsample   = 1000
prmin = 0.0
expr = 'org.smp%d'%(nsample)
lidx_db = range(29*29*29)[1:]

#lvar = [['surfaceTypeIndex','obs']]
#lvar = [['t2m','obs']]
lvar = [['mdhms','obs']]
#lvar = [['Ku_NS_heightStormTop','obs']]
#lvar = [['Latitude','obs'],['Longitude','obs']]
#lvar = [['Latitude','obs'],['Longitude','obs']]
#lvar = [['Ku_NS_heightStormTop','obs'],['Ku_NS_typePrecip','obs'],['t2m','obs']]
#lvar = [['Ku_NS_typePrecip','obs']]
#lvar = [['precip_water_prof_NS','obs'],['Latitude','obs'],['Longitude','obs'],['surfaceTypeIndex','obs']]
ddbvarName = {'surfaceTypeIndex':'surfaceTypeIndex'
         ,'t2m': 't2m'
         ,'Ku_NS_heightStormTop':'Ku_NS_heightStormTop'
         ,'precip_water_prof_NS':'DPRGMI_NS_precipTotWaterCont'
          }

#** Constants ******
myhost = socket.gethostname()
if myhost =='shui':
    dbDir   = '/work/hk01/utsumi/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
    countDir= '/work/hk01/utsumi/PMM/EPCDB/list'
    retbaseDir = '/tank/utsumi/PMM/retsynt/%s'%(expr)

elif myhost == 'well':
    dbDir   = '/media/disk2/share/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
    countDir= '/media/disk2/share/PMM/EPCDB/list'
    retbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retsynt/%s'%(expr)
    figDir   = '/home/utsumi/temp/ret'

else:
    print 'check hostname',myhost
    sys.exit()

db = EPCDB.EPCDB()
#**************************************
#--------------------------------------
#---------------------------
for idx_db in lidx_db:
    print 'idx_db=',idx_db
    irecobsPath = retbaseDir + '/%05d/irec.obs.%05d.npy'%(idx_db,idx_db)
    if not os.path.exists(irecobsPath):
        print 'No file'
        print irecobsPath
        continue

    a1irecobs= np.load(retbaseDir + '/%05d/%s.obs.%05d.npy'%(idx_db,'irec',idx_db)).astype(int32)
    try:
        a1irectop= np.load(retbaseDir + '/%05d/irec.top.%05d.npy'%(idx_db,idx_db)).astype(int32)
        a1idxdbtop= np.load(retbaseDir + '/%05d/idx_db.top.%05d.npy'%(idx_db,idx_db)).astype(int32)
    except IOError:
        print 'No top info'
        print retbaseDir + '/%05d/irec.top.%05d.npy'%(idx_db,idx_db)
        print 'Skip'
    #***************************
    # Use irec to read EPC database
    #---------------------------
    for varName,vartype in lvar:
        if varName in ddbvarName.keys():
            dbvarName = ddbvarName[varName]
        else:
            dbvarName = varName

        if vartype == 'obs':
            avar = np.load(dbDir + '/%s/%05d.npy'%(dbvarName,idx_db))[a1irecobs]

        elif vartype=='top':
            set_idxdbtop = np.sort(list(set(a1idxdbtop)))

            avar = None
            for idx_db_top in set_idxdbtop:
                a1flag = ma.masked_equal(a1idxdbtop, idx_db_top).mask
                a1irectopTmp = a1irectop[a1flag]
                avarTmp  = np.load(dbDir + '/%s/%05d.npy'%(dbvarName,idx_db_top))[a1irectopTmp] 

                if avar is None:
                    nlen = a1irectop.shape[0]
                    dattype= avarTmp.dtype
                    if avarTmp.ndim==1:
                        avar = np.empty(nlen).astype(dattype)
                    else:
                        nz   = avarTmp.shape[1]
                        avar = np.empty([nlen,nz]).astype(dattype)

                avar[a1flag] = avarTmp 


        outPath = retbaseDir + '/%05d/%s.%s.%05d.npy'%(idx_db,varName,vartype,idx_db)
        np.save(outPath, avar)
        print outPath
    #---------------------------

#*******************************
# Save
#-------------------------------

