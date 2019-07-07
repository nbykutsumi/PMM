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




#*******************

#**************************************************************
# main
#--------------------------------------------------------------
db    = JPLDB.JPLDB()
argvs = sys.argv
if len(argvs)==1:
    print '***********************************'
    print ''
    print 'No standard Input'
    print 'Default files and parameters are used'
    print ''
    print '***********************************'
    #** Constants ******
    sensor  = 'GMI'
    coefDir = '/work/hk01/utsumi/JPLDB/EPC_COEF/%s'%(sensor)
    
    #dbDir   = '/work/hk01/utsumi/JPLDB/EPC_DB/GMI_EPC_DATABASE'
    dbDir   = '/work/hk01/utsumi/JPLDB/EPC_DB/GMI_EPC_DATABASE_TEST29'
    
    #-- Single Run ---------
    #clat    = 30.00   # SE.US, oid=016166
    #clon    = 269.0 - 360 # -180 - +180, SE.US, oid=016166
    
    #clat    = 14    # Africa. oid = 002421
    #clon    = 2     # 2014/8/2
    
    oid     = 3556
    clat    = 34    # SE.US case. oid = 003556
    clon    = -86   # 2014/10/14  05:42:03 UTC
    
    #oid     = 12149
    #clat    = 32    # QJRMS case. oid = 012149
    #clon    = -94   # 2016/4/18
    
    
    #clat    = 31.562       
    #clon    = 272.903 -360  # -180 - +180
    #clat    = -47.08
    #clon    = 286.121 -360  # -180 - +180
    #clat    = -9999.
    #clon    = -9999.
    
    
    dlatlon = 3  # used to search the domain center
    iscan   = -9999
    escan   = -9999
    dscan   = 90   # set iscan and dscan =-9999 for entire orbit
    #dscan   = 55
    #dscan   = 5
    
    NEM     = 12
    NTBREG  = 13
    NEM_USE = 3
    #NPCHIST = 25
    NPCHIST = 29
    #NLEV_DPR    = 3  # extract this number of layers
    NLEV_DPR    = 88  # extract this number of layers
    #NLEV_DPR    = 44  # extract this number of layers
    NLEV_PRECIP = 22
    
    
    thwtmin = 0.01
    miss    = -9999.
    
    
    #DB_MAXREC   = 20000
    DB_MAXREC   = 10000
    #DB_MAXREC   = 5000
    DB_MINREC   = 5000
    #NDB_EXPAND  = 10
    NDB_EXPAND  = 20
    #NDB_EXPAND = 0 # test
    DB_RAINFRAC = 0.01  # minimum fraction of precipitating events (>=1mm/h) in the DB required for retrieval
    MAX_T2M_DIFF= 20
        
    #srcPath = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/2017/01/05/1C.GPM.GMI.XCAL2016-C.20170105-S045326-E062600.016220.V05A.HDF5'
    #srcPath = '/home/utsumi/temp/1B.GPM.GMI.TB2016.20171130-S205705-E222939.021348.V05A.HDF5'
    #srcPath = '/home/utsumi/temp/1B.GPM.GMI.TB2016.20171206-S141617-E154850.021437.V05A.HDF5'
    #srcPath = '/home/utsumi/temp/1C.GPM.GMI.XCAL2016-C.20171206-S141617-E154850.021437.V05A.HDF5'
    #srcPath = '/home/utsumi/temp/1C.GPM.GMI.XCAL2016-C.20170101-S173441-E190714.016166.V05A.HDF5'
    #srcPath = '/home/utsumi/temp/1C.GPM.GMI.XCAL2016-C.20170101-S190715-E203948.016167.V05A.HDF5'
    #elevPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.gtopo/2017/01/01/gtopo.002421.npy'
    
    ## 2014/8/2 Africa GMI=002421 **
    #srcPath = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/2014/08/02/1C.GPM.GMI.XCAL2016-C.20140802-S062222-E075455.002421.V05A.HDF5'
    #s2xPath= '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/2014/08/02/Xpy.1.002421.npy'
    #s2yPath= '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/2014/08/02/Ypy.1.002421.npy'
    #tsPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.MERRA2.t2m/2014/08/02/t2m.002421.npy'
    ##tsPath = '/home/utsumi/bin/JPLCODE/EPC_ret_20190221/1C.GPM.GMI.XCAL2016-C.20140802-S062222-E075455.002421.V05A.HDF5.MERRA2.nc'
    #elevPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.gtopo/2014/08/02/gtopo.002421.npy'
        
     
    # 2014/10/14 SE.US GMI=003556 **
    srcPath = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/2014/10/14/1C.GPM.GMI.XCAL2016-C.20141014-S050829-E064102.003556.V05A.HDF5'
    s2xPath= '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/2014/10/14/Xpy.1.003556.npy'
    s2yPath= '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/2014/10/14/Ypy.1.003556.npy'
    tsPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.MERRA2.t2m/2014/10/14/t2m.003556.npy'
    elevPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.gtopo/2014/10/14/gtopo.003556.npy'
    
    
    ## 2016/4/18 QJRMS GMI=012149 **
    #srcPath = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/2016/04/18/1C.GPM.GMI.XCAL2016-C.20160418-S115529-E132803.012149.V05A.HDF5'
    #s2xPath= '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/2016/04/18/Xpy.1.012149.npy'
    #s2yPath= '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/2016/04/18/Ypy.1.012149.npy'
    #tsPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.MERRA2.t2m/2016/04/18/t2m.012149.npy'
    #elevPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.gtopo/2016/04/18/gtopo.012149.npy'
    #

    outDir = '/home/utsumi/temp/out'

#**************************************************************
# Standard input
#**************************************************************
elif len(argvs)>2:
    print 'Too many standart input'
    print 'Python [prog] [parameter-string]'
    sys.exit()
else:
    largvs = argvs[1].split()
    dargv = {}
    for argvs in largvs:
        key,param = argvs.split('=')
        dargv[key] = param

    sensor = dargv['sensor']
    coefDir= dargv['coefDir']
    dbDir  = dargv['dbDir']
    oid    = int(dargv['oid'])
    clat   = float(dargv['clat'])
    clon   = float(dargv['clon'])
    dlatlon= float(dargv['dlatlon'])
    iscan  = int(dargv['iscan'])
    escan  = int(dargv['escan'])
    dscan  = int(dargv['dscan'])
    NEM    = int(dargv['NEM'])
    NTBREG = int(dargv['NTBREG'])
    NEM_USE= int(dargv['NEM_USE']) 
    NPCHIST= int(dargv['NPCHIST'])
    NLEV_DPR = int(dargv['NLEV_DPR'])
    NLEV_PRECIP =int(dargv['NLEV_PRECIP'])
    thwtmin = float(dargv['thwtmin'])
    miss    = float(dargv['miss'])
    

    DB_MAXREC = int(dargv['DB_MAXREC']) 
    DB_MINREC = int(dargv['DB_MINREC']) 
    NDB_EXPAND= int(dargv['NDB_EXPAND'])
    DB_RAINFRAC = float(dargv['DB_RAINFRAC']) # minimum fraction of precipitating events (>=1mm/h) in the DB required for retrieval
    MAX_T2M_DIFF= int(dargv['MAX_T2M_DIFF'])
    
    srcPath   = dargv['srcPath'] 
    s2xPath   = dargv['s2xPath'] 
    s2yPath   = dargv['s2yPath'] 
    tsPath    = dargv['tsPath']   
    elevPath  = dargv['elevPath']
    outDir    = dargv['outDir']

#**************************************************************
# Read parameters
#--------------------------------------------------------------
#-- Read PC coefficient file --
coefPath = coefDir + '/coef_pc.txt'
a2coef   = read_table(coefPath)
a2coef   = a2coef[:,1:]
#-- Read EPC range files --
rangePath = coefDir + '/PC_MIN_MAX_29.txt'
a2pc_edge = read_table(rangePath)

# expand the lowest and highest ranges
a2pc_edge[:,0]   = a2pc_edge[:,0] - 1.e6
a2pc_edge[:,-1]  = a2pc_edge[:,-1]+ 1.e6


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
if (iscan<0)and( (clat<-180)or(180<clat)):
    pass
else:
    if (iscan<0)and(-180<=clat)and(clat <=180):
        idx_c = ret_domain_cy(a2lat, a2lon, clat, dlatlon)
        iscan = idx_c-dscan
        escan = idx_c+dscan

    a3tb  = a3tb  [iscan: escan+1]
    a2lat = a2lat [iscan: escan+1]
    a2lon = a2lon [iscan: escan+1]
    a2ts  = a2ts  [iscan: escan+1]
    a2elev= a2elev[iscan: escan+1]

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
#a2idx_db = epcfunc.mk_epc_id_25bins(a3epc, a2pc_edge)
a2idx_db = epcfunc.mk_epc_id_nbins(a3epc, a2pc_edge, NPCHIST)
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
a2nsurfMS    = ones([nyout,nxout],float32)*miss
a2nsurfNS    = ones([nyout,nxout],float32)*miss
a2nsurfMScmb = ones([nyout,nxout],float32)*miss
a2nsurfNScmb = ones([nyout,nxout],float32)*miss

a3prprofNS   = ones([nyout,nxout,NLEV_PRECIP],float32)*miss
a3prprofNScmb= ones([nyout,nxout,NLEV_PRECIP],float32)*miss

a2top_idxdbMS   = ones([nyout,nxout],int32)*miss
a2top_idxdbNS   = ones([nyout,nxout],int32)*miss

a2top_irecMS   = ones([nyout,nxout],int32)*miss
a2top_irecNS   = ones([nyout,nxout],int32)*miss

a3top_zmMS    = ones([nyout,nxout,NLEV_DPR],int32)*miss
a3top_zmNS    = ones([nyout,nxout,NLEV_DPR],int32)*miss

a3top_prprofNS    = ones([nyout,nxout,NLEV_PRECIP],float32)*miss
a3top_prprofNScmb = ones([nyout,nxout,NLEV_PRECIP],float32)*miss

a3top_tbMS      = ones([nyout,nxout,NTBREG],float32)*miss
a3top_tbNS      = ones([nyout,nxout,NTBREG],float32)*miss


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
    lidx_db_expand_tmp = [idx_db] + [idx_db + i*sign for i in range(1,NDB_EXPAND) for sign in [-1,1]]

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
        a1nsurfMScmbTmp = db.get_var('precip_MS_cmb', nrec=DB_MAXREC)
        a1nsurfNScmbTmp = db.get_var('precip_NS_cmb', nrec=DB_MAXREC)
        a1nsurfMSTmp    = db.get_var('precip_esfc_MS', nrec=DB_MAXREC)
        a1nsurfNSTmp    = db.get_var('precip_esfc_NS', nrec=DB_MAXREC)

        #a2prprofNSTmp   = ma.masked_less(db.get_var('precip_prof_MS',     nrec=DB_MAXREC), 0).filled(0.0)
        a2prprofNSTmp   = db.get_var('precip_prof_NS',     nrec=DB_MAXREC)  # test
        a2prprofNScmbTmp= ma.masked_less(db.get_var('precip_prof_NS_cmb', nrec=DB_MAXREC), 0).filled(0.0)

        a1t2mdbTmp = db.get_var('t2m', nrec=DB_MAXREC) 
        a1tsdbTmp  = db.get_var('ts',  nrec=DB_MAXREC) 
        a1tqvdbTmp = db.get_var('tqv', nrec=DB_MAXREC) 
        a1revdbTmp = db.get_var('rev', nrec=DB_MAXREC) 
        a1elevdbTmp= db.get_var('elev', nrec=DB_MAXREC) 

        a1idxdbTmp = np.ones(a2epcdbTmp.shape[0]).astype(int32)*idx_db_expand
        a1irecTmp  = np.arange(a2epcdbTmp.shape[0]).astype(int32)

        print 'read DB done'
        print 'MS: imported DB record length=%d'%(len(a1nsurfMSTmp))
        print 'NS: imported DB record length=%d'%(len(a1nsurfMSTmp))
 
        #-- Stack data --
        if (iidx_db==0):
            a2epcdb = a2epcdbTmp
            a1nsurfMScmb = a1nsurfMScmbTmp
            a1nsurfNScmb = a1nsurfNScmbTmp
            a1nsurfMS    = a1nsurfMSTmp
            a1nsurfNS    = a1nsurfNSTmp

            a2prprofNS    = a2prprofNSTmp
            a2prprofNScmb = a2prprofNScmbTmp

            a1tsdb  = a1tsdbTmp
            a1revdb = a1revdbTmp
            a1elevdb= a1elevdbTmp

            a1idxdb = a1idxdbTmp
            a1irec  = a1irecTmp

        else:
            a2epcdb = concatenate([a2epcdb,  a2epcdbTmp],axis=0)
            a1nsurfMScmb = concatenate([a1nsurfMScmb,  a1nsurfMScmbTmp], axis=0)
            a1nsurfNScmb = concatenate([a1nsurfNScmb,  a1nsurfNScmbTmp], axis=0)
            a1nsurfMS    = concatenate([a1nsurfMS,     a1nsurfMSTmp], axis=0)
            a1nsurfNS    = concatenate([a1nsurfNS,     a1nsurfNSTmp], axis=0)

            a2prprofNS    = concatenate([a2prprofNS, a2prprofNSTmp], axis=0)
            a2prprofNScmb = concatenate([a2prprofNScmb, a2prprofNScmbTmp], axis=0)

            a1tsdb  = concatenate([a1tsdb,   a1tsdbTmp], axis=0) 
            a1revdb = concatenate([a1revdb,  a1revdbTmp], axis=0) 
            a1elevdb= concatenate([a1elevdb, a1elevdbTmp], axis=0) 

            a1idxdb = concatenate([a1idxdb,  a1idxdbTmp], axis=0)
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

        #-- Only valid precip entries (two scans: NS & MS)--
        ### Make only 2 types(NS & MS) based on DPR  ###
        ### Share for DPR and combined               ###

        a1prflagNS = ma.masked_greater_equal(a1nsurfNS,0).mask
        a1prflagMS = ma.masked_greater_equal(a1nsurfMS,0).mask

        #-- Ts --
        ts    = a2ts[y,x]

        a1tsflag = ma.masked_inside( a1tsdb-ts, -MAX_T2M_DIFF, MAX_T2M_DIFF).mask

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
        
        #-- Screen DB candidates --
        a1flagNS   = a1prflagNS * a1tsflag * a1revflag
        a1flagMS   = a1prflagMS * a1tsflag * a1revflag

        a2epcdbMSSC = a2epcdb[a1flagMS]
        a2epcdbNSSC = a2epcdb[a1flagNS]

        a1nsurfMScmbSC = a1nsurfMScmb[a1flagMS]
        a1nsurfNScmbSC = a1nsurfNScmb[a1flagNS]
        a1nsurfMSSC    = a1nsurfMS   [a1flagMS]
        a1nsurfNSSC    = a1nsurfNS   [a1flagNS]

        a2prprofNSSC   = a2prprofNS   [a1flagNS]
        a2prprofNScmbSC= a2prprofNScmb[a1flagNS]

        a1idxdbMSSC    = a1idxdb[a1flagMS]
        a1idxdbNSSC    = a1idxdb[a1flagNS]

        a1irecMSSC    = a1irec[a1flagMS]
        a1irecNSSC    = a1irec[a1flagNS]



        #print 'screened a2epcdb.shape',a2epcdbSC.shape

        #-- RMSE --
        a1rmsdMS = np.sqrt(np.square((a2epcdbMSSC - a1epc)/a1pc_std).sum(axis=1)/NEM)
        a1rmsdNS = np.sqrt(np.square((a2epcdbNSSC - a1epc)/a1pc_std).sum(axis=1)/NEM)

        idxtopMS = np.argmin(a1rmsdMS)
        idxtopNS = np.argmin(a1rmsdNS)
        rmsd_minMS = a1rmsdMS[idxtopMS]
        rmsd_minNS = a1rmsdNS[idxtopNS]

        #-- Top ranked entris ---
        a2top_idxdbMS[y,x] = a1idxdbMSSC[idxtopMS]
        a2top_idxdbNS[y,x] = a1idxdbNSSC[idxtopNS]
        a2top_irecMS[y,x] = a1irecMSSC[idxtopMS]
        a2top_irecNS[y,x] = a1irecNSSC[idxtopNS]
 
        topidxdbMS = a1idxdbMSSC[idxtopMS]
        topidxdbNS = a1idxdbNSSC[idxtopNS]
        topirecMS  = a1irecMSSC[idxtopMS]
        topirecNS  = a1irecNSSC[idxtopNS]


        # Read top-db file (for MS) --
        dbtopPath = dbDir + '/db_%05d.bin'%(topidxdbMS)
        db.set_file(dbtopPath)
        topirec = topirecMS

        a3top_zmMS[y,x,:] = db.get_var('z_ka', nrec=1, origin=topirec).flatten()[-NLEV_DPR:] 
        a3top_tbMS[y,x,:] = db.get_var('tb', nrec=1, origin=topirec).flatten() 

        # Read top-db file (for NS) --
        dbtopPath = dbDir + '/db_%05d.bin'%(topidxdbNS)
        db.set_file(dbtopPath)
        topirec = topirecNS

        a3top_zmNS[y,x,:] = db.get_var('z_ku', nrec=1, origin=topirec).flatten()[-NLEV_DPR:] 
        a3top_tbNS[y,x,:] = db.get_var('tb', nrec=1, origin=topirec).flatten()

        a3top_prprofNS[y,x,:] = db.get_var('precip_prof_NS', nrec=1, origin=topirecNS)
        a3top_prprofNScmb[y,x,:] = db.get_var('precip_prof_NS_cmb', nrec=1, origin=topirecNS)


        #-- Weight --
        a1wtMS = np.exp(-0.5*np.square(a1rmsdMS/rmsd_minMS))
        a1wtNS = np.exp(-0.5*np.square(a1rmsdNS/rmsd_minNS))
        a1wtMS[idxtopMS] = 1.0
        a1wtNS[idxtopNS] = 1.0

        a1boolwtMS = ma.masked_greater_equal(a1wtMS, thwtmin).mask
        a1boolwtNS = ma.masked_greater_equal(a1wtNS, thwtmin).mask

        a1wtMS = a1wtMS[a1boolwtMS]
        a1wtNS = a1wtNS[a1boolwtNS]

        wtsumMS= a1wtMS.sum()
        wtsumNS= a1wtNS.sum()

        #-- Weighting average --
        nsurfMS    = (a1nsurfMSSC[a1boolwtMS] * a1wtMS).sum() / wtsumMS
        nsurfNS    = (a1nsurfNSSC[a1boolwtNS] * a1wtNS).sum() / wtsumNS
        nsurfNScmb = (a1nsurfNScmbSC[a1boolwtNS] * a1wtNS).sum() / wtsumNS
        nsurfMScmb = (a1nsurfMScmbSC[a1boolwtMS] * a1wtMS).sum() / wtsumMS


        a2nsurfMS[y,x] = nsurfMS
        a2nsurfNS[y,x] = nsurfNS
        a2nsurfMScmb[y,x] = nsurfMScmb
        a2nsurfNScmb[y,x] = nsurfNScmb

        prprofNS   = (a2prprofNSSC[a1boolwtNS] * a1wtNS.reshape(-1,1)).sum(axis=0) / wtsumNS
        prprofNScmb= (a2prprofNScmbSC[a1boolwtNS] * a1wtNS.reshape(-1,1)).sum(axis=0) / wtsumNS

        a3prprofNS[y,x,:]    = prprofNS
        a3prprofNScmb[y,x,:] = prprofNScmb


        #print '%.1f'%nsurfMS, a2prprofNSTmp.min(),a2prprofNSTmp.max()


#--- save (temporary)--
mk_dir(outDir)

stamp = '%06d.y%04d-%04d.nrec%d'%(oid, iscan, escan,DB_MAXREC) 

np.save(outDir + '/nsurfMS.%s.npy'%(stamp), a2nsurfMS)
np.save(outDir + '/nsurfNS.%s.npy'%(stamp), a2nsurfNS)
np.save(outDir + '/nsurfMScmb.%s.npy'%(stamp), a2nsurfMScmb)
np.save(outDir + '/nsurfNScmb.%s.npy'%(stamp), a2nsurfNScmb)

np.save(outDir + '/prprofNS.%s.npy'%(stamp), a3prprofNS)
np.save(outDir + '/prprofNScmb.%s.npy'%(stamp), a3prprofNScmb)
np.save(outDir + '/lat.%s.npy'%(stamp), a2lat)
np.save(outDir + '/lon.%s.npy'%(stamp), a2lon)

np.save(outDir + '/top-idxdbMS.%s.npy'%(stamp), a2top_idxdbMS)
np.save(outDir + '/top-idxdbNS.%s.npy'%(stamp), a2top_idxdbNS)

np.save(outDir + '/top-irecMS.%s.npy'%(stamp), a2top_irecMS)
np.save(outDir + '/top-irecNS.%s.npy'%(stamp), a2top_irecNS)

np.save(outDir + '/top-zmMS.%s.npy'%(stamp), a3top_zmMS)
np.save(outDir + '/top-zmNS.%s.npy'%(stamp), a3top_zmNS)

np.save(outDir + '/top-prprofNS.%s.npy'%(stamp), a3top_prprofNS)
np.save(outDir + '/top-prprofNScmb.%s.npy'%(stamp), a3top_prprofNScmb)

np.save(outDir + '/top-tbMS.%s.npy'%(stamp), a3top_tbMS)
np.save(outDir + '/top-tbNS.%s.npy'%(stamp), a3top_tbNS)



#nsurfMSPath = outDir + '/nsurfMS.%s.npy'%(stamp)
#nsurfNSPath = outDir + '/nsurfNS.%s.npy'%(stamp)
#nsurfMScmbPath = outDir + '/nsurfMScmb.%s.npy'%(stamp)
#nsurfNScmbPath = outDir + '/nsurfNScmb.%s.npy'%(stamp)
#
#prprofNSPath   = outDir + '/prprofNS.%s.npy'%(stamp)
#prprofNScmbPath= outDir + '/prprofNScmb.%s.npy'%(stamp)
#
#latPath   = outDir + '/lat.%s.npy'%(stamp)
#lonPath   = outDir + '/lon.%s.npy'%(stamp)
#
#
#topidxdbNSPath = outDir + '/top-idxdbNS.%s.npy'%(stamp)
#topidxdbMSPath = outDir + '/top-idxdbMS.%s.npy'%(stamp)
#
#topirecNSPath = outDir + '/top-irecNS.%s.npy'%(stamp)
#topirecMSPath = outDir + '/top-irecMS.%s.npy'%(stamp)
#
#topzmNSPath = outDir + '/top-zmNS.%s.npy'%(stamp)
#topzmMSPath = outDir + '/top-zmMS.%s.npy'%(stamp)
#
#topprprofNSPath = outDir + '/top-prprofNS.%s.npy'%(stamp)
#topprprofNScmbPath = outDir + '/top-prprofNScmb.%s.npy'%(stamp)
#
#
#
#np.save(nsurfMSPath, a2nsurfMS)
#np.save(nsurfNSPath, a2nsurfNS)
#np.save(nsurfMScmbPath, a2nsurfMScmb)
#np.save(nsurfNScmbPath, a2nsurfMScmb)
#
#np.save(prprofNSPath, a3prprofNS)
#np.save(prprofNScmbPath, a3prprofNScmb)
#np.save(latPath, a2lat)
#np.save(lonPath, a2lon)
#
#np.save(topidxdbMSPath, a2top_idxdbMS)
#np.save(topidxdbNSPath, a2top_idxdbNS)
#
#np.save(topirecMSPath, a2top_irecMS)
#np.save(topirecNSPath, a2top_irecNS)
#
#np.save(topzmMSPath, a3top_zmMS)
#np.save(topzmNSPath, a3top_zmNS)
#
#np.save(topprprofNSPath, a3top_prprofNS)
#np.save(topprprofNScmbPath, a3top_prprofNScmb)

print outDir
print stamp
