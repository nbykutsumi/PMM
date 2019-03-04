from numpy import *
import h5py
import sys, os
import numpy as np
import JPLDB
from bisect import bisect_left
import epcfunc
#from math import acos, cos, sin

sensor  = 'GMI'
clat    = 30.00   # SE.US, oid=016166
clon    = 269.0 - 360 # -180 - +180, SE.US, oid=016166
#clat    = 31.562       
#clon    = 272.903 -360  # -180 - +180
#clat    = -47.08
#clon    = 286.121 -360  # -180 - +180
#clat    = -9999.
#clon    = -9999.

dlatlon = 3  # used to search the domain center
dscan   = 30

NEM     = 12
NTBREG  = 13
NEM_USE = 3
NPCHIST = 25

db      = JPLDB.JPLDB()

thwtmin = 0.01 
miss    = -9999.


DB_MAXREC   = 20000
DB_MINREC   = 100
NDB_EXPAND  = 20
DB_RAINFRAC = 0.01  # minimum fraction of precipitating events (>=1mm/h) in the DB required for retrieval


#srcPath = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/2017/01/05/1C.GPM.GMI.XCAL2016-C.20170105-S045326-E062600.016220.V05A.HDF5'
#srcPath = '/home/utsumi/temp/1B.GPM.GMI.TB2016.20171130-S205705-E222939.021348.V05A.HDF5'
#srcPath = '/home/utsumi/temp/1B.GPM.GMI.TB2016.20171206-S141617-E154850.021437.V05A.HDF5'
#srcPath = '/home/utsumi/temp/1C.GPM.GMI.XCAL2016-C.20171206-S141617-E154850.021437.V05A.HDF5'
srcPath = '/home/utsumi/temp/1C.GPM.GMI.XCAL2016-C.20170101-S173441-E190714.016166.V05A.HDF5'
#srcPath = '/home/utsumi/temp/1C.GPM.GMI.XCAL2016-C.20170101-S190715-E203948.016167.V05A.HDF5'

s2xPath= '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/2017/01/01/Xpy.1.016166.npy'
s2yPath= '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/2017/01/01/Ypy.1.016166.npy'

coefDir = '/work/hk01/utsumi/JPLDB/EPC_COEF/%s'%(sensor)
dbDir   = '/work/hk01/utsumi/JPLDB/EPC_DB/GMI_EPC_DATABASE'

#-- functions -----
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except OSError:
    pass


def read_table(srcPath, type=float):
    f=open(srcPath,'r'); lines=f.readlines(); f.close()
    lout = []
    for line in lines:
        line = line.strip().split()
        line = map(type, line)
        lout.append(line)
    return array(lout)


def read_nrain(idx_db):
    '''
     /* first six= Ntotal then Nrain exceeding 1 5 10 20 50 mm/hr when T2m < 278K */
     /* second six= same for when T2m > 278K */
    '''
    srcPath = dbDir + '/db_%05d.bin.nrain.txt'%(idx_db)
    a1nrain = read_table(srcPath,type=int32)[0]
    a1nrain_cold = a1nrain[:6]
    a1nrain_warm = a1nrain[6:]
    return a1nrain_warm, a1nrain_cold 

#def dist_simple(lat1, lon1, lat2, lon2):
#    RADEARTH = 6371
#    RTD = 57.29578
#    DTR = 0.017453
#
#    dist = RADEARTH*acos(cos(DTR*lon1-DTR*lon2)*cos(DTR*lat1)*cos(DTR*lat2) + sin(DTR*lat1)*sin(DTR*lat2))
#    return dist



#-- Read PC coefficient file --
coefPath = coefDir + '/coef_pc.txt'
#coefPath = '/home/utsumi/bin/ENSPR/coef_pc.txt'
a2coef   = read_table(coefPath)
a2coef   = a2coef[:,1:]

#-- Read EPC range files --
rangePath = coefDir + '/PC_MIN_MAX_25_no_overlap.txt'
#rangePath = '/home/utsumi/bin/ENSPR/PC_MIN_MAX_10_no_overlap.txt'
a2tmp     = read_table(rangePath)
a2tmp     = a2tmp[:,3:]

nytmp, nxtmp = a2tmp.shape
a2pc_edge     = zeros([nytmp,nxtmp/2+1]) # 12*(10+1)
a2pc_edge_min = zeros([nytmp,nxtmp/2])   # 12*10
a2pc_edge_max = zeros([nytmp,nxtmp/2])   # 12*10
for ytmp in range(nytmp):
    for xtmp in range(0,nxtmp,2):
        a2pc_edge_min[ytmp,xtmp/2] = a2tmp[ytmp,xtmp]
        a2pc_edge_max[ytmp,xtmp/2] = a2tmp[ytmp,xtmp+1]

a2pc_edge[:,:-1] = a2pc_edge_min
a2pc_edge[:,-1]  = a2pc_edge_max[:,-1]

#-- Read PC ave and std file --
pcavePath  = coefDir + '/ave_pc.txt'
#pcavePath  = '/home/utsumi/bin/ENSPR/ave_pc.txt'
a2pc_avestd= read_table(pcavePath)
a1pc_ave   = a2pc_avestd[:,1]
a1pc_std   = a2pc_avestd[:,2]


#-- Read granule data data --
with h5py.File(srcPath, 'r') as h5:
    a3tb1    = h5['/S1/Tc'][:]
    a3tb2org = h5['/S2/Tc'][:]
    a2lat = h5['/S1/Latitude'][:]
    a2lon = h5['/S1/Longitude'][:]

    #a2lat2= h5['/S2/Latitude'][:]
    #a2lon2= h5['/S2/Longitude'][:]

#-- Matchup and Joint S1 and S2 Tb --
a1x2  = np.load(s2xPath).flatten()
a1y2  = np.load(s2yPath).flatten()

a1mask= ma.masked_less(a1x2,0).mask
a1x2  = ma.masked_less(a1x2,0).filled(0)
a1y2  = ma.masked_less(a1y2,0).filled(0)

nytmp, nxtmp, ztmp = a3tb2org.shape
a2tb2 = a3tb2org[a1y2, a1x2]
a2tb2[a1mask] = miss
a3tb2 = a2tb2.reshape(nytmp,nxtmp,-1)
a3tb = concatenate([a3tb1, a3tb2],axis=2)

#****************************************************
# Extract target domain
#----------------------------------------------------
nyTmp, nxTmp = a2lat.shape
a1lat = a2lat[:,nxTmp/2]
a1lon = a2lon[:,nxTmp/2]

idx_latmax = np.argmax(a1lat)
a1lat0 = a1lat[:idx_latmax+1]
a1lat1 = a1lat[idx_latmax+1:]
a1lon0 = a1lon[:idx_latmax+1]
a1lon1 = a1lon[idx_latmax+1:]

if (-180<=clat)and(clat <=180):
    #-- search first half: ascending --
    found_domain = 0
    idx_c  = bisect_left(a1lat0, clat)
    latTmp = a1lat0[idx_c]
    lonTmp = a1lon0[idx_c]
    if (clat-dlatlon<=latTmp)&(latTmp <=clat+dlatlon)&(clon-dlatlon<=lonTmp)&(lonTmp<=clon+dlatlon):
        found_domain = 1
    else:
        #-- search second half: descending --
        idx_c  = bisect_left(a1lat1[::-1], clat)
        idx_c  = len(a1lat) - idx_c -1
        latTmp = a1lat[idx_c]
        lonTmp = a1lon[idx_c]

        if (clat-dlatlon<=latTmp)&(latTmp <=clat+dlatlon)&(clon-dlatlon<=lonTmp)&(lonTmp<=clon+dlatlon):
            found_domain =1
    
    if found_domain==1:
        idx_first = idx_c - dscan
        idx_last  = idx_c + dscan    
        a3tb    = a3tb [idx_first:idx_last+1,:] 
        a2lat   = a2lat[idx_first:idx_last+1,:]        
        a2lon   = a2lon[idx_first:idx_last+1,:]

    else:
        print 'No matching scans in the target domain are found.'
        print 'Exit'
        sys.exit()

    print 'Extract target domain'
    print 'Extracted array size=', a3tb.shape

else:
    pass

#-- test ---------------
#a3tb  = a3tb[:10,:,:]
#iy    = 5
#ix    = 215
#emtemp= 0
##-----------------------
#****************************************************
# Make mask data
#----------------------------------------------------
#-- missing tb ---
a2mask = np.any(ma.masked_outside(a3tb, 50, 350).mask, axis=2)

#****************************************************
# Convert Tb to EPC
#----------------------------------------------------
print 'calc epc'
a3epc = epcfunc.mk_epc_12pc(a3tb, a2coef)
print 'calc epc done'

#****************************************************
# Find EPC bin numbers 
#----------------------------------------------------
print 'calc idx'
a2idx_db = epcfunc.mk_epc_id_25bins(a3epc, a2pc_edge)
print 'calc idx done'

print a2mask.shape
print a2idx_db.shape
a2idx_db = ma.masked_where(a2mask, a2idx_db).filled(miss)
#-----------------
lidxset  = list(set(a2idx_db.flatten()))
lidxset  = sort(lidxset)
#****************************************************
# DBidx-y-x mapping
#----------------------------------------------------
nyout,nxout = a2idx_db.shape

#-- Initialize output array ---
a2esurf_NS_cmb = ones([nyout,nxout],float32)*miss

#-- Start retrieve --
X,Y = meshgrid(range(nxout),range(nyout))

for i,idx_db in enumerate(lidxset):
    if idx_db==-9999: continue

    a2bool = ma.masked_equal(a2idx_db, idx_db).mask
    a1x    = X[a2bool]
    a1y    = Y[a2bool]
    print i,'/',len(lidxset), 'idx_db=%d pixels=%d'%(idx_db, len(a1x))

    #******************************
    #- check # of DB entries (expand to neighborhood, if necessary)
    #------------------------------
    lidx_db_expand_tmp = [idx_db] + [idx_db + i*sign for i in range(NDB_EXPAND) for sign in [-1,1]]

    lidx_db_expand = []    
    nevent_all   = 0
    nevent_cold  = 0
    nevent_warm  = 0
    nevent_cold1 = 0
    nevent_warm1 = 0

    for idx_db_expand in lidx_db_expand_tmp:
        #-- If idx_db == -9999 --
        if ((idx_db_expand <0)or(pow(NPCHIST, NEM_USE)-1<idx_db_expand)):
            print 'No matching database'
            print 'idx_db=',idx_db_expand
            continue
    
        #-- Read nrain file --
        try:
            a1nrain_warm, a1nrain_cold = read_nrain(idx_db_expand)
        except IOError:
            print 'No DB file for idx_db=',idx_db_expand
            print 'SKIP'
            continue
    
        nevent_all  = nevent_all  + a1nrain_cold[0] + a1nrain_warm[0]
        nevent_cold = nevent_cold + a1nrain_cold[0]
        nevent_warm = nevent_warm + a1nrain_warm[0]
        nevent_cold1= nevent_cold1+ a1nrain_cold[1]
        nevent_warm1= nevent_warm1+ a1nrain_warm[1]

        lidx_db_expand.append(idx_db_expand)
        if nevent_all >DB_MINREC: break

    frac0 = nevent_cold1/float(nevent_cold)
    frac1 = nevent_warm1/float(nevent_warm)

    #******************************
    #- Check rain factio in DB
    #------------------------------
    if ((frac0<DB_RAINFRAC)&(frac1<DB_RAINFRAC)):
        print 'Insufficient rain>1 frac0=%.3f frac1=%.3f Skip idx_db=%d'%(frac0, frac1, idx_db_expand)
        print '%.3f %.3f DB_RAINFRAC=%.3f'%(frac0,frac1, DB_RAINFRAC)
        continue

    print 'ndbrec =',
    #******************************
    #- search DB (expand to neighborhood, if necessary)
    #------------------------------
    for idx_db_expand in lidx_db_expand:
        #-- If idx_db == -9999 --
    
        #-- Read database file --
        #print 'set file'   
        dbPath = dbDir + '/db_%05d.bin'%(idx_db_expand)
        db.set_file(dbPath)
        #print 'set file done' 
        #print 'read DB'
        a2epcdb = db.get_var('pc_emis', nrec=DB_MAXREC)[:,:NEM]  # (nrec, 12)
        #print 'read DB done'
        a1esurf_NS_cmb = ma.masked_less(db.get_var('precip_NS_cmb', nrec=DB_MAXREC), 0).filled(0.0)
        #print a2epcdb
    
        #-- Extract subset --
        if (idx_db_expand==idx_db):
            a2epcdb = a2epcdb[:DB_MAXREC,:]
            a1esurf_NS_cmb = a1esurf_NS_cmb[:DB_MAXREC]
        else:
            a2epcdb = concatenate([a2epcdb,a2epcdb[:DB_MAXREC,:]],axis=0)
            a1esurf_NS_cmb = concatenate([a1esurf_NS_cmb, a1esurf_NS_cmb[:DB_MAXREC]])


    #-- Start loop over y,x with the same idx_db --
    '''
    # a1epc   : (11)
    # a3epcdb : (nrec,11) 
    '''
    print 'a2epcdb.shape, len(a1x)=',a2epcdb.shape, len(a1x)
    for (y,x) in zip(a1y,a1x):
        #print 'idx_db y x=',idx_db,y,x
        a1epc = a3epc[y,x,:]
        a1rmsd= np.sqrt(np.square((a2epcdb - a1epc)/a1pc_std).sum(axis=1)/NEM)
        idxtop= np.argmin(a1rmsd)
        rmsd_min= a1rmsd[idxtop]

        #-- Weight --
        a1wt = np.exp(-0.5*np.square(a1rmsd/rmsd_min))
        a1wt[idxtop] = 1.0
        a1boolwt = ma.masked_greater_equal(a1wt, thwtmin).mask
        a1wt = a1wt[a1boolwt]
        wtsum= a1wt.sum()

        #-- Weighting average --
        esurf_NS_cmb = (a1esurf_NS_cmb[a1boolwt] * a1wt).sum() / wtsum
        a2esurf_NS_cmb[y,x] = esurf_NS_cmb
        #print '1esurf=',a1esurf_NS_cmb
        #print 'esurf=',esurf_NS_cmb

    #sys.exit()

#--- save (temporary)--
outDir = '/home/utsumi/temp/out'
mk_dir(outDir)
esurfPath = outDir + '/esurf.%06d.npy'%(oid)
latPath   = outDir + '/lat.%06d.npy'%(oid)
lonPath   = outDir + '/lon.%06d.npy'%(oid)

np.save(esurfPath, a2esurf_NS_cmb)
np.save(latPath, a2lat)
np.save(lonPath, a2lon)
print esurfPath
print lidxset

