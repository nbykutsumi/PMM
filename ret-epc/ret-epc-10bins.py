from numpy import *
import h5py
import sys, os
import numpy as np
import JPLDB
from bisect import bisect_left

clat    = 31.562       
clon    = 272.903 -360  # -180 - +180
dlatlon = 3  # used to search the domain center
dscan   = 30

NEM     = 11
NTBREG  = 9
NEM_USE = 4
db      = JPLDB.JPLDB()

thwtmin = 0.01 
miss    = -9999.

DB_MAXREC   = 20000
DB_RAINFRAC = 0.1  # minimum fraction of precipitating events in the DB required for retrieval


#srcPath = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/2017/01/05/1C.GPM.GMI.XCAL2016-C.20170105-S045326-E062600.016220.V05A.HDF5'
#srcPath = '/home/utsumi/temp/1B.GPM.GMI.TB2016.20171130-S205705-E222939.021348.V05A.HDF5'
srcPath = '/home/utsumi/temp/1B.GPM.GMI.TB2016.20171206-S141617-E154850.021437.V05A.HDF5'

coefDir = '/home/utsumi/bin/ENSPR'
dbDir   = '/work/a01/utsumi/data/JPLDB/GMI_dbase_PC_prof2_V5_NPC4_BIN10_OVERLAP1'

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


#-- Read PC coefficient file --
coefPath = coefDir + '/coef_pc.txt'
a2coef   = read_table(coefPath)
a2coef   = a2coef[:,1:]

#-- Read EPC range files --
rangePath = coefDir + '/PC_MIN_MAX_10_no_overlap.txt'
a2tmp     = read_table(rangePath)
a2tmp     = a2tmp[:,3:]

nytmp, nxtmp = a2tmp.shape
a2pc_range     = zeros([nytmp,nxtmp/2+1]) # 11*(10+1)
a2pc_range_min = zeros([nytmp,nxtmp/2])   # 11*10
a2pc_range_max = zeros([nytmp,nxtmp/2])   # 11*10
for ytmp in range(nytmp):
    for xtmp in range(0,nxtmp,2):
        a2pc_range_min[ytmp,xtmp/2] = a2tmp[ytmp,xtmp]
        a2pc_range_max[ytmp,xtmp/2] = a2tmp[ytmp,xtmp+1]

a2pc_range[:,:-1] = a2pc_range_min
a2pc_range[:,-1]  = a2pc_range_max[:,-1]

#-- Read PC ave and std file --
pcavePath  = coefDir + '/ave_pc.txt'
a2pc_avestd= read_table(pcavePath)
a1pc_ave   = a2pc_avestd[:,1]
a1pc_std   = a2pc_avestd[:,2]

#-- Read granule data data --
h5   = h5py.File(srcPath, 'r')
#a3tb_full = h5['/S1/Tc'][:]
a3tb  = h5['/S1/Tb'][:]
a2lat = h5['/S1/Latitude'][:]
a2lon = h5['/S1/Longitude'][:]

h5.close()

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

if (-180<=clat)or(clat <=180):
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
else:
    pass

print 'Extract target domain'
print 'Extracted array size=', a3tb.shape
#-- test ---------------
#a3tb  = a3tb[:10,:,:]
#iy    = 5
#ix    = 215
#emtemp= 0
##-----------------------
#****************************************************
# Convert Tb to EPC
#----------------------------------------------------
#'''
ny,nx,ntmp = a3tb.shape
a3epc = zeros([ny,nx,NEM])

# +a0
kt  = 0
for iem in range(NEM):
    a3epc[:,:,iem] = a2coef[kt,iem]

for i in range(NTBREG):   # NTBRERG=9
    print 'iTBREG=',i,'out of',NTBREG
    kt =kt+1

    # + bi*Tbi
    for iem in range(NEM): # NEM = 11
        b = a2coef[kt, iem]
        a3epc[:,:,iem] = a3epc[:,:,iem] + b*a3tb[:,:,i]

        ##-- test --
        #if iem==emtemp:
        #    print 'i,b, Tb',i,b,a3tb[iy,ix,i]
        ##----------

    # + cij * Tbi * Tbj
    for j in range(i,NTBREG):
        kt = kt+1
        #print i,j
        a2tmp = a3tb[:,:,i] * a3tb[:,:,j]

        for iem in range(NEM):
            c = a2coef[kt, iem]
            a3epc[:,:,iem] = a3epc[:,:,iem] + c*a3tb[:,:,i]*a3tb[:,:,j]

            ##-- test --
            #if iem==emtemp:
            #    print 'i,j,c, Tbi*Tbj',i,j,c,a3tb[iy,ix,i]*a3tb[iy,ix,j]
            ##----------
    
# + d * (V-H)/(V+H)
kt = kt+1
for iem in range(NEM):
    d = a2coef[kt,iem]
    a3epc[:,:,iem] = a3epc[:,:,iem] \
         +d*(a3tb[:,:,0]-a3tb[:,:,1])/(a3tb[:,:,0]+a3tb[:,:,1])

kt = kt+1
for iem in range(NEM):
    d = a2coef[kt,iem]
    a3epc[:,:,iem] = a3epc[:,:,iem] \
         +d*(a3tb[:,:,2]-a3tb[:,:,3])/(a3tb[:,:,2]+a3tb[:,:,3])

kt = kt+1
for iem in range(NEM):
    d = a2coef[kt,iem]
    a3epc[:,:,iem] = a3epc[:,:,iem] \
         +d*(a3tb[:,:,5]-a3tb[:,:,6])/(a3tb[:,:,5]+a3tb[:,:,6])

#****************************************************
# Find EPC bin numbers 
#----------------------------------------------------
a2idx_db  = zeros([ny,nx],int32)
for iem in range(NEM_USE):
    a1bin = a2pc_range[iem]
    a2idxTmp = np.digitize(a3epc[:,:,iem], a1bin, right=False) - 1
    print a2idx_db.shape, a2idxTmp.shape

    a2idxTmp = ma.masked_outside(a2idxTmp,0,9)
    a2idx_db = a2idx_db + a2idxTmp*pow(10, iem)

a2idx_db = a2idx_db.filled(-9999)

#-- for test ----
epcPath  = '/home/utsumi/temp/temp.epc.npy'
np.save(epcPath, a3epc)

idxPath  = '/home/utsumi/temp/temp.idx.npy'
np.save(idxPath, a2idx_db)
#'''

epcPath  = '/home/utsumi/temp/temp.epc.npy'
a3epc    = np.load(epcPath)
idxPath  = '/home/utsumi/temp/temp.idx.npy'
a2idx_db = np.load(idxPath)
lidxset  = list(set(a2idx_db.flatten()))
lidxset  = sort(lidxset)
#****************************************************
# DBidx-y-x mapping
#----------------------------------------------------
nyout,nxout = a2idx_db.shape

#-- Initialize output array ---
a2esurf_MS_cmb = ones([nyout,nxout],float32)*miss

#-- Start retrieve --
X,Y = meshgrid(range(nxout),range(nyout))

for i,idx_db in enumerate(lidxset):
    a2bool = ma.masked_equal(a2idx_db, idx_db).mask
    a1x    = X[a2bool]
    a1y    = Y[a2bool]
    print 'all i idx_db len=',len(lidxset), i, idx_db, len(a1x)

    #-- If idx_db == -9999 --
    if idx_db==miss:
        print 'No matching database'
        print 'idx_db=',idx_db
        continue

    #-- Read nrain file --
    try:
        a1nrain_warm, a1nrain_cold = read_nrain(idx_db)
    except IOError:
        print 'No DB file for idx_db=',idx_db
        print 'SKIP'
        continue

    nnorain_warm = a1nrain_warm[0]
    nnorain_cold = a1nrain_cold[0]
    if nnorain_warm==0:
        prain_warm = 0.
    else:
        prain_warm = float(a1nrain_warm[1:].sum())/a1nrain_warm[0]
    if nnorain_cold==0:
        prain_cold = 0.
    else:
        prain_cold = float(a1nrain_cold[1:].sum())/a1nrain_cold[0]

    #-- Set maxrec to be used for reatrievals --
    maxrec = DB_MAXREC
    if   (a1nrain_warm[4]!=0)or(a1nrain_cold[4]!=0):
        maxrec = int(DB_MAXREC)
    elif (a1nrain_warm[3]!=0)or(a1nrain_cold[3]!=0):
        maxrec = int(0.1*DB_MAXREC)
    elif (a1nrain_warm[2]!=0)or(a1nrain_cold[2]!=0):
        if (prain_warm<DB_RAINFRAC)&(prain_cold<DB_RAINFRAC):
            a2esurf_MS_cmb[a1y,a1x]=0.0
            print 'too small number of precip event in DB:%d'%(idx_db)
            continue
        else:
            maxrec = int(0.01*DB_MAXREC)
    else:
            maxrec = int(0.01*DB_MAXREC)

    print 'maxrec=',maxrec,DB_MAXREC

    #-- Read database file --
    print 'set file'   
    dbPath = dbDir + '/db_%05d.bin'%(idx_db)
    db.set_file(dbPath)
    print 'set file done' 
    print 'read DB'
    a2epcdb = db.get_var('pc_emis', nrec=DB_MAXREC)  # (nrec, 11)
    print 'read DB done'
    a1esurf_MS_cmb = ma.masked_less(db.get_var('precip_MS_cmb', nrec=DB_MAXREC), 0).filled(0.0)
    #print a2epcdb

    #-- Extract subset --
    a2epcdb = a2epcdb[:maxrec,:]
    a1esurf_MS_cmb = a1esurf_MS_cmb[:maxrec]

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
        esurf_MS_cmb = (a1esurf_MS_cmb[a1boolwt] * a1wt).sum() / wtsum
        a2esurf_MS_cmb[y,x] = esurf_MS_cmb
        #print '1esurf=',a1esurf_MS_cmb
        #print 'esurf=',esurf_MS_cmb

    #sys.exit()

#--- save (temporary)--
outDir = '/home/utsumi/temp/out'
mk_dir(outDir)
esurfPath = outDir + '/esurf.npy'
latPath   = outDir + '/lat.npy'
lonPath   = outDir + '/lon.npy'

np.save(esurfPath, a2esurf_MS_cmb)
np.save(latPath, a2lat)
np.save(lonPath, a2lon)


print a3epc.shape
#sys.exit()
##-- test ---
#print '------------------------'
#a1tb   = a3tb[iy,ix,:]
#a1coef = a2coef[:,emtemp]
#
#epc =0
#kt = 0
#a  = a1coef[0]
#
#epc = epc+a
#
#for i in range(9):
#    kt = kt +1
#    b  = a1coef[kt]
#    tb = a1tb[i]
#    print 'i,b,Tb', i,b,tb
#
#    epc=epc + b*tb
#    for j in range(i,9):
#        kt = kt+1
#        c  = a1coef[kt]
#        tb0= a1tb[i]
#        tb1= a1tb[j]
#        print 'i,j,c, Tbi*Tbj',i,j,c,tb0*tb1
#        epc=epc+ c*tb0*tb1
#
#kt = kt+1
#d  = a1coef[kt]
#tb0= a1tb[0]
#tb1= a1tb[1] 
#epc = epc + d*(tb0-tb1)/(tb0+tb1)
#print 'i,d, (V-H)/(V+H)',i,d,(tb0-tb1)/(tb0+tb1)
#
#kt = kt+1
#d  = a1coef[kt]
#tb0= a1tb[2]
#tb1= a1tb[3] 
#epc = epc + d*(tb0-tb1)/(tb0+tb1)
#print 'i,d, (V-H)/(V+H)',i,d,(tb0-tb1)/(tb0+tb1)
#
#kt = kt+1
#d  = a1coef[kt]
#tb0= a1tb[5]
#tb1= a1tb[6] 
#epc = epc + d*(tb0-tb1)/(tb0+tb1)
#print 'i,d, (V-H)/(V+H)',i,d,(tb0-tb1)/(tb0+tb1)
#
#print ''
#print 'epc0=',a3epc[iy,ix,emtemp]
#print 'epc1=',epc
#print 'iy,ix=',iy,ix
#print a3epc[iy,ix,:]
#
##print a2coef
##--- test end ----     

