import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import h5py
import sys, os
import numpy as np
import JPLDB
from bisect import bisect_left
import epcfunc
#from math import acos, cos, sin
import netCDF4

sensor  = 'GMI'
#clat    = 30.00   # SE.US, oid=016166
#clon    = 269.0 - 360 # -180 - +180, SE.US, oid=016166

#clat    = 14    # Africa. oid = 002421
#clon    = 2     # 2014/8/2

clat    = 32    # QJRMS case. oid = 012149
clon    = -94   # 2016/4/18



#clat    = 31.562       
#clon    = 272.903 -360  # -180 - +180
#clat    = -47.08
#clon    = 286.121 -360  # -180 - +180
#clat    = -9999.
#clon    = -9999.

dlatlon = 3  # used to search the domain center
#dscan   = 55
dscan   = 30
#dscan   = 0

NEM     = 12
NTBREG  = 13
NEM_USE = 3
NPCHIST = 25
NLEV_PRECIP = 22

db      = JPLDB.JPLDB()

thwtmin = 0.01
miss    = -9999.


DB_MAXREC   = 20000
#DB_MAXREC   = 2000
DB_MINREC   = 5000
NDB_EXPAND  = 10
DB_RAINFRAC = 0.01  # minimum fraction of precipitating events (>=1mm/h) in the DB required for retrieval
MAX_T2M_DIFF= 20

#srcPath = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/2017/01/05/1C.GPM.GMI.XCAL2016-C.20170105-S045326-E062600.016220.V05A.HDF5'
#srcPath = '/home/utsumi/temp/1B.GPM.GMI.TB2016.20171130-S205705-E222939.021348.V05A.HDF5'
#srcPath = '/home/utsumi/temp/1B.GPM.GMI.TB2016.20171206-S141617-E154850.021437.V05A.HDF5'
#srcPath = '/home/utsumi/temp/1C.GPM.GMI.XCAL2016-C.20171206-S141617-E154850.021437.V05A.HDF5'
#srcPath = '/home/utsumi/temp/1C.GPM.GMI.XCAL2016-C.20170101-S173441-E190714.016166.V05A.HDF5'
#srcPath = '/home/utsumi/temp/1C.GPM.GMI.XCAL2016-C.20170101-S190715-E203948.016167.V05A.HDF5'
#elevPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.gtopo/2017/01/01/gtopo.002421.npy'

# 2016/4/18 QJRMS GMI=012149 **
srcPath = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/2016/04/18/1C.GPM.GMI.XCAL2016-C.20160418-S115529-E132803.012149.V05A.HDF5'
s2xPath= '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/2016/04/18/Xpy.1.012149.npy'
s2yPath= '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/2016/04/18/Ypy.1.012149.npy'
tsPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.MERRA2.t2m/2016/04/18/t2m.012149.npy'
elevPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.gtopo/2016/04/18/gtopo.012149.npy'

## 2014/8/2 Africa GMI=002421 **
#srcPath = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/2014/08/02/1C.GPM.GMI.XCAL2016-C.20140802-S062222-E075455.002421.V05A.HDF5'
#s2xPath= '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/2014/08/02/Xpy.1.002421.npy'
#s2yPath= '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/2014/08/02/Ypy.1.002421.npy'
#tsPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.MERRA2.t2m/2014/08/02/t2m.002421.npy'
##tsPath = '/home/utsumi/bin/JPLCODE/EPC_ret_20190221/1C.GPM.GMI.XCAL2016-C.20140802-S062222-E075455.002421.V05A.HDF5.MERRA2.nc'
#elevPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.gtopo/2014/08/02/gtopo.002421.npy'

#** Constants ******
coefDir = '/work/hk01/utsumi/JPLDB/EPC_COEF/%s'%(sensor)
dbDir   = '/work/hk01/utsumi/JPLDB/EPC_DB/GMI_EPC_DATABASE'
#*******************

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

def ret_domain_cy(a2lat, a2lon, clat, dlatlon):
    nyTmp, nxTmp = a2lat.shape
    a1lat = a2lat[:,nxTmp/2]
    a1lon = a2lon[:,nxTmp/2]
    
    idx_latmax = np.argmax(a1lat)
    a1lat0 = a1lat[:idx_latmax+1]
    a1lat1 = a1lat[idx_latmax+1:]
    a1lon0 = a1lon[:idx_latmax+1]
    a1lon1 = a1lon[idx_latmax+1:]
    
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
        return idx_c

    else:
        print 'No matching scans in the target domain are found.'
        print 'Exit'
        sys.exit()
 


#-- Orbit id -------------
oid = int(srcPath.split('/')[-1].split('.')[-3])

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

# expand the lowest and highest ranges
a2pc_edge[:,0]   = a2pc_edge[:,0] - 1.e6
a2pc_edge[:,-1]  = a2pc_edge[:,-1]+ 1.e6
a2pc_edge_min[:,0] = a2pc_edge_min[:,0] - 1.e6
a2pc_edge_max[:,-1]= a2pc_edge_max[:,-1]+ 1.e6


#-- Read PC ave and std file --
pcavePath  = coefDir + '/ave_pc.txt'
#pcavePath  = '/home/utsumi/bin/ENSPR/ave_pc.txt'
a2pc_avestd= read_table(pcavePath)
a1pc_ave   = a2pc_avestd[:,1]
a1pc_std   = a2pc_avestd[:,2]


#-- Read granule data --
with h5py.File(srcPath, 'r') as h5:
    a3tb1    = h5['/S1/Tc'][:]
    a3tb2org = h5['/S2/Tc'][:]
    a2lat = h5['/S1/Latitude'][:]
    a2lon = h5['/S1/Longitude'][:]


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


##--- test: replace tb data with JPL's --------
#print '*'*40
#print 'CAUTION! Tb data is replaced.'
#print '*'*40
#tmpPath = '/home/utsumi/temp/out/tb-jpl-full-002421.npy'
#a3tb    = np.load(tmpPath)

#-- Read MERRA2 data ---------
a2ts = np.load(tsPath)
#nc = netCDF4.Dataset(tsPath)
#a2ts = nc.variables['t2m'][:]
#-- Read elevation data ---------
a2elev = np.load(elevPath)

#****************************************************
# Extract target domain
#----------------------------------------------------
if (-180<=clat)and(clat <=180):
    idx_c = ret_domain_cy(a2lat, a2lon, clat, dlatlon)

    #idx_c = 2059 #  test
    
    a3tb  = a3tb  [idx_c-dscan: idx_c+dscan+1]
    a2lat = a2lat [idx_c-dscan: idx_c+dscan+1]
    a2lon = a2lon [idx_c-dscan: idx_c+dscan+1]
    a2ts  = a2ts  [idx_c-dscan: idx_c+dscan+1]
    a2elev= a2elev[idx_c-dscan: idx_c+dscan+1]

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

##-- test: replace idx_db with JPL's ---
#print '*'*40
#print 'CAUTION!!'
#print 'IDX is replaced for test'
#print '*'*40
#a2idx_db = np.load('/home/utsumi/temp/out/idx-jpl-full-002421.npy')
#a2idx_db = a2idx_db[idx_c-dscan: idx_c+dscan+1]
##--------------------------------------

a2idx_db = ma.masked_where(a2mask, a2idx_db).filled(miss)

##-- test ---------
#print 'lat,lon=',a2lat[50,125], a2lon[50,125]
#print 'epc-vect',a3epc[50,125]
#print 'idx_db=    ',a2idx_db[50,125]
#
#print 'idx_db_man=',25**2*5 + 25**1*24 + 0
#plt.imshow(ma.masked_less(a2idx_db,0))
#plt.colorbar()
#plt.savefig('/home/utsumi/temp/out/temp.idx.png')
#
#np.save('/home/utsumi/temp/out/temp.idx.npy', a2idx_db)
#
#sys.exit()
#-----------------
lidxset  = list(set(a2idx_db.flatten()))
lidxset  = sort(lidxset)
#****************************************************
# DBidx-y-x mapping
#----------------------------------------------------
nyout,nxout = a2idx_db.shape

#-- Initialize output array ---
a2esurf  = ones([nyout,nxout],float32)*miss
a3prprof = ones([nyout,nxout,NLEV_PRECIP],float32)*miss

#-- Start retrieve --
X,Y = meshgrid(range(nxout),range(nyout))
for i,idx_db in enumerate(lidxset):
    if idx_db==-9999: continue


    ##***** test **************
    #if idx_db !=3077: continue
    ##***** test **************


    a2bool = ma.masked_equal(a2idx_db, idx_db).mask
    a1x    = X[a2bool]
    a1y    = Y[a2bool]
    print i,'/',len(lidxset), 'idx_db=%d pixels=%d'%(idx_db, len(a1x))

    #******************************
    #- check Num of DB entries (expand to neighborhood, if necessary)
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
    #- Check rain ratio in DB
    #------------------------------
    if ((frac0<DB_RAINFRAC)&(frac1<DB_RAINFRAC)):
        print 'Insufficient rain>1 frac0=%.3f frac1=%.3f Skip idx_db=%d'%(frac0, frac1, idx_db_expand)
        print '%.3f %.3f DB_RAINFRAC=%.3f'%(frac0,frac1, DB_RAINFRAC)
        continue

    print 'ndbrec =',
    #******************************
    #- Read DB (expand to neighborhood, if necessary)
    #------------------------------
    for iidx_db, idx_db_expand in enumerate(lidx_db_expand):
        #-- If idx_db == -9999 --
    
        #-- Read database file --
        #print 'set file'   
        dbPath = dbDir + '/db_%05d.bin'%(idx_db_expand)
        db.set_file(dbPath)
        print 'set file done' 
        print 'read DB'
        a2epcdbTmp = db.get_var('pc_emis', nrec=DB_MAXREC)[:,:NEM]  # (nrec, 12)
        a1esurfMScmbTmp = db.get_var('precip_MS_cmb', nrec=DB_MAXREC)
        a1esurfNScmbTmp = db.get_var('precip_NS_cmb', nrec=DB_MAXREC)
        a1esurfMSTmp    = db.get_var('precip_esfc_MS', nrec=DB_MAXREC)
        a1esurfNSTmp    = db.get_var('precip_esfc_NS', nrec=DB_MAXREC)

        a2prprofMSTmp   = ma.masked_less(db.get_var('precip_prof_MS',     nrec=DB_MAXREC), 0).filled(0.0)
        a2prprofNSTmp   = ma.masked_less(db.get_var('precip_prof_MS',     nrec=DB_MAXREC), 0).filled(0.0)


        a1tsdbTmp  = db.get_var('t2m', nrec=DB_MAXREC) 
        a1revdbTmp = db.get_var('rev', nrec=DB_MAXREC) 
        a1elevdbTmp= db.get_var('elev', nrec=DB_MAXREC) 
        print 'read DB done'
        print 'imported DB record length=%d'%(len(a1esurfTmp))
        a1irecTmp  = arange(len(a1esurfTmp)).astype(int32) 
 
        #-- Stack data --
        if (iidx_db==0):
            a2epcdb = a2epcdbTmp
            a1esurfMScmb = a1esurfMScmbTmp
            #a1esurfNScmb = a1esurfNScmbTmp
            a1esurfMS    = a1esurfMSTmp
            #a1esurfNS    = a1esurfNSTmp

            a2prprof= a2prprofTmp

            a1tsdb  = a1tsdbTmp
            a1revdb = a1revdbTmp
            a1elevdb= a1elevdbTmp
            a1irec  = a1irecTmp

        else:
            a2epcdb = concatenate([a2epcdb,  a2epcdbTmp],axis=0)
            a1esurfMScmb = concatenate([a1esurfMScmb,  a1esurfMScmbTmp], axis=0)
            a1esurfNScmb = concatenate([a1esurfNScmb,  a1esurfNScmbTmp], axis=0)
            a1esurfMS    = concatenate([a1esurfMS,     a1esurfMSTmp], axis=0)
            a1esurfNS    = concatenate([a1esurfNS,     a1esurfNSTmp], axis=0)

            a2prprof= concatenate([a2prprof, a2prprofTmp], axis=0)

            a1tsdb  = concatenate([a1tsdb,   a1tsdbTmp], axis=0) 
            a1revdb = concatenate([a1revdb,  a1revdbTmp], axis=0) 
            a1elevdb= concatenate([a1elevdb, a1elevdbTmp], axis=0) 
            a1irec  = concatenate([a1irec,   a1irecTmp], axis=0) 

    #******************************
    #-- Start loop over y,x with the same idx_db --
    #******************************
    '''
    # a1epc   : (11)
    # a3epcdb : (nrec,11) 
    '''
    print 'lidx_db_expand=',lidx_db_expand
    print 'a2epcdb.shape, len(a1x)=',a2epcdb.shape, len(a1x)

    for (y,x) in zip(a1y,a1x):  # in idx_db loop
        ####***** test **************
        #if x !=98: continue
        ####***** test **************

        #-- Obs EPC --
        #print 'idx_db y x=',idx_db,y,x
        a1epc = a3epc[y,x,:]

        #********************
        # Constrain candidates from DB
        #********************
        #-- Discard entries from same granule (revolution) --
        a1revflag = ma.masked_not_equal(a1revdb, oid).mask

        #-- Only valid precip entries --
        a1prflag  = ma.masked_greater_equal(a1esurf,0).mask

        #-- Ts --
        ts    = a2ts[y,x]

        a1tsflag = ma.masked_inside( a1tsdb-ts, -MAX_T2M_DIFF, MAX_T2M_DIFF).mask
        #if ts < 0:
        #    a1tsflag = ma.masked_less(a1tsdb, 0).mask
        #else:
        #    a1tsflag = ma.masked_greater_equal(a1tsdb, 0).mask

        ##-- Elevation --
        #elev = a2elev[y,x]
        #if elev < 500:
        #    a1elevflag = ma.masked_less(a1elevdb, 500).mask
        #elif (500 <=elev)and(elev < 1000):
        #    a1elevflag = ma.masked_inside(a1elevdb, 500, 1000).mask
        #elif 1000 <=elev:
        #    a1elevflag = ma.masked_greater_equal(a1elevdb, 1000).mask

        #else:
        #    print 'check elev',elev
        #    sys.exit()
        
        ##-- test -----------
        #irec = 19187
        #print 'prflag=',a1prflag[irec]
        #print 'tsflag=',a1tsflag[irec]
        #print 'revflag=',a1revflag[irec]
        #sys.exit()

        #-- Screen DB candidates --
        a1flag    = a1prflag * a1tsflag * a1revflag
        #a1flag    = a1tsflag * a1revflag

        a2epcdbSC = a2epcdb[a1flag]
        a1esurfSC = a1esurf[a1flag]
        a2prprofSC= a2prprof[a1flag]
        a1irecSC  = a1irec[a1flag]

        #print 'screened a2epcdb.shape',a2epcdbSC.shape

        #-- RMSE --
        a1rmsd= np.sqrt(np.square((a2epcdbSC - a1epc)/a1pc_std).sum(axis=1)/NEM)
        idxtop= np.argmin(a1rmsd)
        rmsd_min= a1rmsd[idxtop]

           
        #-- Weight --
        a1wt = np.exp(-0.5*np.square(a1rmsd/rmsd_min))
        a1wt[idxtop] = 1.0


        ###-- test -----
        #a1tmpidx= arange(len(a1wt)).astype(int32)
        #a2tmp   = zip(a1wt, a1tmpidx)
        #a2tmp.sort(key=lambda x: x[0], reverse=True)
        #a1iorg = array(zip(*a2tmp)[1])
        #for itmp,iorg in enumerate(a1iorg):
        #    irec = a1irecSC[iorg]
        #    #lirec = [17732,17938,18891]
        #    #if irec in lirec:
        #    #    print i, a1irecSC[i], 'esurf=%.3f'%a1esurfSC[i], 'wt=%.3f'%(a1wt[i])
        #    print itmp, iorg, a1irecSC[iorg], 'esurf=%.3f'%a1esurfSC[iorg], 'wt=%.3f'%(a1wt[iorg])
        #    if a1wt[iorg]<thwtmin: break
        ##---------------

        a1boolwt = ma.masked_greater_equal(a1wt, thwtmin).mask
        a1wt = a1wt[a1boolwt]
        wtsum= a1wt.sum()

        #-- Weighting average --
        esurf = (a1esurfSC[a1boolwt] * a1wt).sum() / wtsum
        a2esurf[y,x] = esurf

        prprof= (a2prprofSC[a1boolwt] * a1wt.reshape(-1,1)).sum(axis=0) / wtsum
        a3prprof[y,x,:] = prprof


        ##***** test **************
        #print 'Tb(my)=',['%.2f'%x for x in a3tb[y,x,:]]
        #print 'a1wt.shape',a1wt.shape
        #print 'wtsum=',wtsum
        #print 'idx_db=',idx_db
        #print 'precip=',esurf
        #if x ==100: sys.exit()
        ##***** test **************

#--- save (temporary)--
outDir = '/home/utsumi/temp/out'
mk_dir(outDir)
esurfPath = outDir + '/esurf.%06d.y%04d-%04d.nrec%d.npy'%(oid, idx_c-dscan, idx_c+dscan, DB_MAXREC)
prprofPath= outDir + '/prprof.%06d.y%04d-%04d.nrec%d.npy'%(oid, idx_c-dscan, idx_c+dscan,DB_MAXREC)
latPath   = outDir + '/lat.%06d.y%04d-%04d.nrec%d.npy'%(oid, idx_c-dscan, idx_c+dscan,DB_MAXREC)
lonPath   = outDir + '/lon.%06d.y%04d-%04d.nrec%d.npy'%(oid, idx_c-dscan, idx_c+dscan,DB_MAXREC)

np.save(esurfPath, a2esurf)
np.save(prprofPath, a3prprof)
np.save(latPath, a2lat)
np.save(lonPath, a2lon)
print esurfPath
#print lidxset

