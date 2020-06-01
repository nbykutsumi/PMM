# %%
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
from numpy import *
import h5py
import sys, os, socket
import numpy as np
import JPLDB
import EPCDB
from bisect import bisect_left
import epcfunc
import random

#lidx_db = [1681]
#lidx_db = [23484]  # N ~60,000
#lidx_db = [9680]   # N ~10,000
#lidx_db = [5498]   # N ~10,000
#lidx_db = [14155]  # N ~3500
#lidx_db = [4721]
lidx_db = range(29*29*29)[1:]
#lidx_db = range(29*29*29)[4731:]
#lidx_db = range(3000,3500)
#lidx_db = range(23516,29*29*29)[1:]
#lidx_db = range(23543,29*29*29)[1:]
#lidx_db = range(1000)
#lidx_db = [5810]
#lidx_db = [3005]
nsample = 1000

#dbtype  = 'my'
dbtype  = 'jpl'

DB_MAXREC    = 10000
DB_MINREC    = 1000
#sensor  = 'GMI'
#sensor  = 'AMSR2'
#sensor  = 'SSMIS'
#sensor  = 'ATMS'
sensor  = 'MHS'
#expr = 'org.%s.smp%d'%(sensor,nsample)
expr = 'prf.%s.smp%d'%(sensor,nsample)
#expr = 'no-same-rev.smp%d'%(nsample)

#** Constants ******
myhost = socket.gethostname()
if myhost =='shui':
    dbDir   = '/tank/utsumi/PMM/JPLDB/EPC_DB/%s_EPC_DATABASE'%(sensor)
    coefDir = '/tank/utsumi/PMM/JPLDB/EPC_COEF/%s'%(sensor)
    outbaseDir = '/tank/utsumi/PMM/retsynt/%s'%(expr)

elif myhost == 'well':
    #dbDir   = '/home/utsumi/mnt/lab_tank/utsumi/PMM/JPLDB/EPC_DB/%s_EPC_DATABASE'%(sensor)
    #coefDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/JPLDB/EPC_COEF/%s'%(sensor)
    #countDir= '/home/utsumi/mnt/lab_tank/utsumi/PMM/EPCDB/list'
    dbDir   = '/media/disk2/share/PMM/JPLDB/EPC_DB/%s_EPC_DATABASE'%(sensor)
    coefDir = '/media/disk2/share/PMM/JPLDB/EPC_COEF/%s'%(sensor)
    relprofDir = '/media/disk2/share/PMM/JPLDB/EPC_DB/%s_rs_precip_water_prof_NS'%(sensor)

    outbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retsynt/%s'%(expr)

else:
    print 'check dtype',dtype
    sys.exit()


#** Parameters ************    

tqvflag = 0
elevflag= 0

if sensor =='GMI':
    NEM     = 12
    NTBREG  = 13
elif sensor == 'AMSR2':
    NEM     = 13
    NTBREG  = 10
elif sensor == 'SSMIS':
    NEM     = 10
    NTBREG  = 7
elif sensor == 'ATMS':
    NEM     = 6
    NTBREG  = 9
elif sensor == 'MHS':
    NEM     = 4
    NTBREG  = 5


if sensor in ['MHS','ATMS']:
    flag_crosstrack = 1
else:
    flag_crosstrack = 0


NEM_USE = 3
#NPCHIST = 25
NPCHIST = 29
if dbtype=='jpl':
    NLEV_PRECIP = 50
elif dbtype=='my':
    NLEV_PRECIP = 50

#N2MAX = 200
N2MAX = -9999
#thwtmin = 0.5 # test
thwtmin = 0.2

miss    = -9999.

DB_USE_MINREC= 200
NDB_EXPAND  = 20
#NDB_EXPAND = 0 # test
DB_RAINFRAC = 0.01  # minimum fraction of precipitating events (>=1mm/h) in the DB required for retrieval
MAX_INC_DIFF = 20  # degree
MAX_T2M_DIFF= 10  # K
MAX_TQV_DIFF= 10  # kg/m2
MAX_RMA_0   = 0.05
flag_top_var= 0
flag_prof   = 1
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
    #srcPath = dbDir + '/db_%05d.bin.nrain.txt'%(idx_db)
    nrainDir= dbDir 
    srcPath = nrainDir + '/db_%05d.bin.nrain.txt'%(idx_db)
    a1nrain = read_table(srcPath,type=int32)[0]
    a1nrain_cold = a1nrain[:6]
    a1nrain_warm = a1nrain[6:]
    return a1nrain_warm, a1nrain_cold 


def read_rnr_table(idx_db):
    '''
    Rain/No Rain screening
    --- Contents of files ----
    calc 0:Skip all pixels  1:Screening  2:Calc all pixels
    D :  Linear discriminator threshold.
         If discriminator for the pixelsi smaller than this threshold, assign 'no-rain' to the pixel.
    SR:  Skip Ratio (expected fraction of skipped pixels)
    WER: Wet Event Ratio
    RMA: Ratio of Missing Amount  (>=0, 0.1, 1, 5, 10mm/h)
    m_org-s_no: First 12: For 12 EPCs. 13th: T2m
    '''

    rnrDir  = dbDir + '/rnr.minrec%d'%(DB_MINREC)
    srcPath = rnrDir + '/rnr.%05d.txt'%(idx_db)
    f=open(srcPath,'r'); lines=f.readlines(); f.close()
    #print srcPath
    rnrflag = int(lines[0].split('\t')[1])
    thD     = float(lines[1].split('\t')[1])
    SR      = float(lines[2].split('\t')[1])
    WER     = float(lines[3].split('\t')[1])
    lRMA    = map(float, lines[4].split('\t')[1:])
    ave_org   = np.array(map(float, lines[5].split('\t')[1:]))  # for normalization
    ave_rain  = np.array(map(float, lines[6].split('\t')[1:]))  # Averages of normalized variables
    ave_no    = np.array(map(float, lines[7].split('\t')[1:]))  # Averages of normalized variables
    std_org   = np.array(map(float, lines[8].split('\t')[1:]))  # for normalization
    std_rain  = np.array(map(float, lines[9].split('\t')[1:]))  # Averages of normalized variables
    std_no    = np.array(map(float, lines[10].split('\t')[1:]))  # Averages of normalized variables


    return rnrflag, thD, SR, WER, lRMA, ave_org, ave_rain, ave_no, std_org, std_rain, std_no



#*******************

#**************************************************************
# main
#--------------------------------------------------------------
if dbtype =='jpl':
    db    = JPLDB.JPLDB(sensor)
elif dbtype=='my':
    db    = EPCDB.EPCDB()
else:
    print 'check dbtype',dbtype

#**************************************************************
# Read parameters
#--------------------------------------------------------------
##-- Read PC coefficient file --
#coefPath = coefDir + '/coef_pc.txt'
#a2coef   = read_table(coefPath)
#a2coef   = a2coef[:,1:]
##-- Read EPC range files --
#rangePath = coefDir + '/PC_MIN_MAX_29.txt'
#a2pc_edge = read_table(rangePath)
#
## expand the lowest and highest ranges
#a2pc_edge[:,0]   = a2pc_edge[:,0] - 1.e6
#a2pc_edge[:,-1]  = a2pc_edge[:,-1]+ 1.e6


#-- Read PC ave and std file --
pcavePath  = coefDir + '/ave_pc.txt'
#pcavePath  = '/home/utsumi/bin/ENSPR/ave_pc.txt'
a2pc_avestd= read_table(pcavePath)
a1pc_ave   = a2pc_avestd[:,1]
a1pc_std   = a2pc_avestd[:,2]


#****************************************************
# Start idx_db loop
#****************************************************
for idx_db in lidx_db:
    print '************************************'
    print 'idx_db for primary loop =',idx_db

    try:
        a1nrain_warm, a1nrain_cold = read_nrain(idx_db)
    except IOError:
        print 'No primary DB file for idx_db=',idx_db
        print 'SKIP'
        continue

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
            #print 'No matching database'
            #print 'idx_db=',idx_db_expand
            continue
    
        #-- Read nrain file --
        try:
            a1nrain_warm, a1nrain_cold = read_nrain(idx_db_expand)
        except IOError:
            #print 'No DB file for idx_db=',idx_db_expand
            #print 'SKIP'
            continue
    
        nevent_all  = nevent_all  + a1nrain_cold[0] + a1nrain_warm[0]
        nevent_cold = nevent_cold + a1nrain_cold[0]
        nevent_warm = nevent_warm + a1nrain_warm[0]
        nevent_cold1= nevent_cold1+ a1nrain_cold[1]
        nevent_warm1= nevent_warm1+ a1nrain_warm[1]

        lidx_db_expand.append(idx_db_expand)
        if nevent_all >DB_MINREC: break

    if nevent_cold==0:
        frac0 = 0.0
    else:
        frac0 = nevent_cold1/float(nevent_cold)
    if nevent_warm==0:
        frac1 = 0.0
    else:
        frac1 = nevent_warm1/float(nevent_warm)


    ##******************************
    ##- Check rain ratio in DB
    ##------------------------------
    #if ((frac0<DB_RAINFRAC)&(frac1<DB_RAINFRAC)):
    #    print 'Insufficient rain>1 frac0=%.3f frac1=%.3f Skip idx_db=%d'%(frac0, frac1, idx_db_expand)
    #    print '%.3f %.3f DB_RAINFRAC=%.3f'%(frac0,frac1, DB_RAINFRAC)
    #    continue

    #******************************
    #- Read DB (expand to neighborhood, if necessary)
    #------------------------------
    for iidx_db, idx_db_expand in enumerate(lidx_db_expand):
        #-- Read database file --
        #print 'set file'   
        dbPath = dbDir + '/db_%05d.bin'%(idx_db_expand)
        if   dbtype == 'jpl':
            db.set_file(dbPath)
        elif dbtype == 'my':
            db.set_idx_db(dbDir, idx_db_expand)

        #print 'read DB'
        a2epcdbTmp = db.get_var('pc_emis')[:,:NEM]  # (nrec, 12)
        a1nsurfMScmbTmp = db.get_var('precip_MS_cmb')
        a1nsurfNScmbTmp = db.get_var('precip_NS_cmb')
        a1nsurfMSTmp    = db.get_var('precip_nsfc_MS')
        a1nsurfNSTmp    = db.get_var('precip_nsfc_NS')

        if flag_prof == 1:
            #a2prwatprofNSTmp = ma.masked_invalid(db.get_var('precip_water_prof_NS')[:,-NLEV_PRECIP:]).filled(-9999.)  # 2019/11/12

            a2prwatprofNSTmp = ma.masked_invalid(np.load(relprofDir + '/%05d.npy'%(idx_db_expand))[:,-NLEV_PRECIP:])
            a2prwatprofNSTmp = (ma.masked_less(a2prwatprofNSTmp,0)/1000.).filled(-9999.)

        elif flag_prof==0:
            nrecTmp = a1nsurfNScmbTmp.shape[0] 
            a2prwatprofNSTmp = np.ones(nrecTmp).reshape(-1,NLEV_PRECIP)*-9999.
        else:
            print 'check flag_prof', flag_prof
            sys.exit()

           

        #a1tsdbTmp  = db.get_var('ts') 
        a1t2mdbTmp = db.get_var('t2m') 
        a1revdbTmp = db.get_var('rev') 
        if tqvflag ==1:
            a1tqvdbTmp = db.get_var('tqv') 
        if elevflag ==1:
            a1elevdbTmp= db.get_var('elev') 
        if flag_crosstrack==1:
            a1incdbTmp   = db.get_var('inc_S1')


        a1idxdbTmp = np.ones(a2epcdbTmp.shape[0]).astype(int32)*idx_db_expand
        a1irecTmp  = np.arange(a2epcdbTmp.shape[0]).astype(int32)

        #-- Stack data --
        if (iidx_db==0):
            a2epcdb = a2epcdbTmp
            a1nsurfMScmb = a1nsurfMScmbTmp
            a1nsurfNScmb = a1nsurfNScmbTmp
            a1nsurfMS    = a1nsurfMSTmp
            a1nsurfNS    = a1nsurfNSTmp

            #a2prprofNS    = a2prprofNSTmp
            #a2prprofNScmb = a2prprofNScmbTmp
            a2prwatprofNS = a2prwatprofNSTmp

            #a1tsdb  = a1tsdbTmp
            a1t2mdb = a1t2mdbTmp
            a1revdb = a1revdbTmp
            if tqvflag ==1:
                a1tqvdb = a1tqvdbTmp
            if elevflag ==1:
                a1elevdb= a1elevdbTmp
            if flag_crosstrack==1:
                a1incdb = a1incdbTmp


            a1idxdb = a1idxdbTmp
            a1irecsubdb = a1irecTmp

        else:
            a2epcdb = concatenate([a2epcdb,  a2epcdbTmp],axis=0)
            a1nsurfMScmb = concatenate([a1nsurfMScmb,  a1nsurfMScmbTmp], axis=0)
            a1nsurfNScmb = concatenate([a1nsurfNScmb,  a1nsurfNScmbTmp], axis=0)
            a1nsurfMS    = concatenate([a1nsurfMS,     a1nsurfMSTmp], axis=0)
            a1nsurfNS    = concatenate([a1nsurfNS,     a1nsurfNSTmp], axis=0)

            #a2prprofNS    = concatenate([a2prprofNS, a2prprofNSTmp], axis=0)
            #a2prprofNScmb = concatenate([a2prprofNScmb, a2prprofNScmbTmp], axis=0)
            a2prwatprofNS = concatenate([a2prwatprofNS, a2prwatprofNSTmp])

            #a1tsdb  = concatenate([a1tsdb,   a1tsdbTmp], axis=0) 
            a1t2mdb  = concatenate([a1t2mdb,  a1t2mdbTmp], axis=0) 
            a1revdb  = concatenate([a1revdb,  a1revdbTmp], axis=0) 
            if tqvflag ==1:
                a1tqvdb  = concatenate([a1tqvdb,  a1tqvdbTmp], axis=0) 
            if elevflag ==1:
                a1elevdb= concatenate([a1elevdb, a1elevdbTmp], axis=0) 
            if flag_crosstrack==1:
                a1incdb = concatenate([a1incdb, a1incdbTmp],   axis=0)



            a1idxdb = concatenate([a1idxdb,  a1idxdbTmp], axis=0)
            a1irecsubdb = concatenate([a1irecsubdb,   a1irecTmp], axis=0)


    #***********************************************
    # Read Primary idx database file --
    # assumed to be observation
    #***********************************************
    
    #print 'set file'   
    dbPath = dbDir + '/db_%05d.bin'%(idx_db)
    if   dbtype == 'jpl':
        db.set_file(dbPath)
    elif dbtype == 'my':
        db.set_idx_db(dbDir, idx_db)
    else:
        print 'check dtype', dbtype
        sys.exit()

    #print 'read DB'
    a2epcObs = db.get_var('pc_emis')[:,:NEM]  # (nrec, 12)
    a1nsurfMScmbObs = db.get_var('precip_MS_cmb')
    a1nsurfNScmbObs = db.get_var('precip_NS_cmb')
    a1nsurfMSObs    = db.get_var('precip_nsfc_MS')
    a1nsurfNSObs    = db.get_var('precip_nsfc_NS')
    #a2prwatprofNSObs = db.get_var('precip_water_prof_NS')[:,-NLEV_PRECIP:]
    if flag_prof==1:
        #a2prwatprofNSObs = ma.masked_invalid(db.get_var('precip_water_prof_NS')[:,-NLEV_PRECIP:]).filled(-9999.) # 2019/11/12

        a2prwatprofNSObs = ma.masked_invalid(np.load(relprofDir + '/%05d.npy'%(idx_db))[:,-NLEV_PRECIP:])
        a2prwatprofNSObs = (ma.masked_less(a2prwatprofNSObs,0)/1000.).filled(-9999.)

        print relprofDir + '/%05d.npy'%(idx_db)

    elif flag_prof==0:
        nrecTmp = a1nsurfNScmbObs.shape[0] 
        a2prwatprofNSObs = np.ones(nrecTmp).reshape(-1,NLEV_PRECIP)*-9999.

    
    #a1tsdbTmp  = db.get_var('ts') 
    a1t2mObs = db.get_var('t2m') 
    a1revObs = db.get_var('rev') 
    if tqvflag ==1:
        a1tqvObs = db.get_var('tqv') 
    if elevflag ==1:
        a1elevObs= db.get_var('elev') 
    if flag_crosstrack==1:
        a1incObs = db.get_var('inc_S1') 
    if sensor not in ['GMI','AMSR2']:
        a1satidObs=db.get_var('satid')

    #****************************************************
    # Mask invalid db data
    #****************************************************
    a2epcdb      = ma.masked_invalid(a2epcdb     )
    a1nsurfMScmb = ma.masked_invalid(a1nsurfMScmb)
    a1nsurfNScmb = ma.masked_invalid(a1nsurfNScmb)
    a1nsurfMS    = ma.masked_invalid(a1nsurfMS   )
    a1nsurfNS    = ma.masked_invalid(a1nsurfNS   )

    a1nsurfMScmb = ma.masked_less(a1nsurfMScmb,0) 
    a1nsurfNScmb = ma.masked_less(a1nsurfNScmb,0) 
    a1nsurfMS    = ma.masked_less(a1nsurfMS   ,0) 
    a1nsurfNS    = ma.masked_less(a1nsurfNS   ,0)     

    a1epcmaskdb = a2epcdb.mask.any(axis=1)
    a1nsurfmask1 = a1nsurfNS.mask
    a1nsurfmask2 = a1nsurfMS.mask
    a1nsurfmask3 = a1nsurfNScmb.mask
    a1nsurfmask4 = a1nsurfMScmb.mask
    a1nsurfmaskdb  = a1nsurfmask1 + a1nsurfmask2 + a1nsurfmask3 + a1nsurfmask4

    a1maskdb = a1epcmaskdb + a1nsurfmaskdb
    
    a2epcdb[a1maskdb] = miss
    a2epcdb = ma.masked_equal(a2epcdb, miss)



    #****************************************************
    # Mask invalid obs data
    #****************************************************
    a2epcObs     = ma.masked_invalid(a2epcObs     )
    a1nsurfMScmbObs = ma.masked_invalid(a1nsurfMScmbObs)
    a1nsurfNScmbObs = ma.masked_invalid(a1nsurfNScmbObs)
    a1nsurfMSObs    = ma.masked_invalid(a1nsurfMSObs   )
    a1nsurfNSObs    = ma.masked_invalid(a1nsurfNSObs   )

    a1nsurfmask1 = ma.masked_less(a1nsurfMScmbObs,0).mask 
    a1nsurfmask2 = ma.masked_less(a1nsurfNScmbObs,0).mask 
    a1nsurfmask3 = ma.masked_less(a1nsurfMSObs   ,0).mask 
    a1nsurfmask4 = ma.masked_less(a1nsurfNSObs   ,0).mask     

    a1epcmaskObs = a2epcObs.mask.any(axis=1)
    a1nsurfmaskObs  = a1nsurfmask1 + a1nsurfmask2 + a1nsurfmask3 + a1nsurfmask4

    a1maskObs = a1epcmaskObs + a1nsurfmaskObs

    #**********************************************
    # Loop over obaservation (OBS DB)
    #**********************************************
    random.seed(0)
    lirecOrg = np.arange(a2epcObs.shape[0]).astype(int32)
    lirecOrg = ma.masked_where(a1maskObs, lirecOrg).compressed()
    nsampleTmp = min(nsample, len(lirecOrg))
    lirec = np.sort(random.sample(lirecOrg, nsampleTmp)) # No duplication

    nrecTmp = len(lirec)

    if nrecTmp ==0: continue


    #****************************************************
    # Initialize output variables
    #****************************************************
    a1nsurfNSest     = np.ones(len(lirec), float32)*miss
    a1nsurfMSest     = np.ones(len(lirec), float32)*miss
    a1nsurfNScmbest  = np.ones(len(lirec), float32)*miss
    a1nsurfMScmbest  = np.ones(len(lirec), float32)*miss
    a2prwatprofNSest = np.ones([len(lirec),NLEV_PRECIP], float32)*miss

    # Comment out 2020/03/07
    #a1nsurfNSobs     = ma.masked_invalid(a1nsurfNS).filled(miss)[lirec]
    #a1nsurfMSobs     = ma.masked_invalid(a1nsurfMS).filled(miss)[lirec]
    #a1nsurfNScmbobs  = ma.masked_invalid(a1nsurfNScmb).filled(miss)[lirec]
    #a1nsurfMScmbobs  = ma.masked_invalid(a1nsurfMScmb).filled(miss)[lirec]

    # 2020/03/07
    a1nsurfNSObs     = ma.masked_invalid(a1nsurfNSObs).filled(miss)[lirec]
    a1nsurfMSObs     = ma.masked_invalid(a1nsurfMSObs).filled(miss)[lirec]
    a1nsurfNScmbObs  = ma.masked_invalid(a1nsurfNScmbObs).filled(miss)[lirec]
    a1nsurfMScmbObs  = ma.masked_invalid(a1nsurfMScmbObs).filled(miss)[lirec]
    a2prwatprofNSObs = ma.masked_invalid(a2prwatprofNSObs).filled(miss)[lirec]
    #

    a1irec           = np.array(lirec).astype(int32)
    a1topirec        = np.ones(len(lirec), int32)*miss
    a1topidxdb       = np.ones(len(lirec), int32)*miss

    if sensor not in ['GMI','AMSR2']:
        a1satidout   = np.ones(len(lirec), float32)*miss

    #****************************************************
    for iirec,irec in enumerate(lirec):
        a1epcObs = a2epcObs[irec]
    
        #**********************************
        # Mask same oid (rev) entries
        #---------------------------------- 
        rev    = a1revObs[irec]
        a1maskrev = ma.masked_equal(a1revdb, rev).mask
        a1maskirec= ma.masked_equal(a1irecsubdb, irec).mask

        a1masktrue= a1maskrev * a1maskirec
        #a1masktrue= a1maskrev  # for no-same-rev
        #**********************************
        # Mask with t2m
        #---------------------------------- 
        t2m = a1t2mObs[irec]
        a1maskt2m = ma.masked_outside(a1t2mdb-t2m, -MAX_T2M_DIFF, MAX_T2M_DIFF).mask 


        #**********************************
        # Mask with incident angle difference
        #---------------------------------- 
        if flag_crosstrack ==1:
            incObs = a1incObs[irec]
            a1maskinc = ma.masked_greater(np.abs(a1incdb) - abs(incObs), MAX_INC_DIFF).mask

        else:
            a1maskinc = np.array([False])
        #**********************************
        # Total Mask
        #---------------------------------- 
        a1mask = a1masktrue + a1maskt2m + a1maskinc
        #a1mask = False  # test
        #---------------------------------- 


        a1rmsd = np.sqrt(np.square((a2epcdb - a1epcObs)/a1pc_std).sum(axis=1)/NEM)
        a1rmsd = ma.masked_where(a1mask, a1rmsd) 
        rmsd_min= np.min(a1rmsd)
        if rmsd_min==0:
            rmsd_min = ma.masked_less_equal(a1rmsd,0).min()

        irectop= np.argmin(a1rmsd)

        a1wt   = np.exp(-0.5*np.square(a1rmsd/rmsd_min))
        a1wt   = ma.masked_less(a1wt, thwtmin)
        a1wt[irectop] = 0.607

        #**********************************
        # Top-ranked entry info
        #---------------------------------- 
        a1topirec[iirec] = a1irecsubdb[irectop]
        a1topidxdb[iirec]= a1idxdb[irectop]

        #**********************************
        # Constrain maximum number of entries
        #---------------------------------- 
        if N2MAX >0:
            if a1wt.count >N2MAX: 
                a1wtsort = np.sort(a1wt.compressed())[::-1]
                if a1wtsort.shape[0] > N2MAX:
                    thwtTmp = a1wtsort[N2MAX]
                    a1wt    = ma.masked_less(a1wt, thwtTmp)
        #----------------------------------

        wtsum  = a1wt.sum()
        #-- Weighting average --
        nsurfMS    = (a1nsurfMS * a1wt).sum() / wtsum
        nsurfNS    = (a1nsurfNS * a1wt).sum() / wtsum
        nsurfNScmb = (a1nsurfNScmb * a1wt).sum() / wtsum
        nsurfMScmb = (a1nsurfMScmb * a1wt).sum() / wtsum

        #nsurfNStop = a1nsurfNS[irectop]

        a1prwatprofNS = (ma.masked_less(a2prwatprofNS,0) * a1wt.reshape(-1,1)).sum(axis=0) / wtsum

        #--- Output array ------
        a1nsurfNSest    [iirec] = nsurfNS 
        a1nsurfMSest    [iirec] = nsurfMS 
        a1nsurfNScmbest [iirec] = nsurfNScmb
        a1nsurfMScmbest [iirec] = nsurfMScmb
        a2prwatprofNSest[iirec,:] = a1prwatprofNS 

        if sensor not in ['GMI','AMSR2']:
            satid    = a1satidObs[irec]
            a1satidout  [iirec] = satid

        #print iirec,'/',nrecTmp,'%.2f, %.2f, %.2f'%(a1nsurfNS[irec], nsurfNS, nsurfNStop)

    a1nsurfNSest     = ma.masked_invalid(a1nsurfNSest     ).filled(miss)
    a1nsurfMSest     = ma.masked_invalid(a1nsurfMSest     ).filled(miss)
    a1nsurfNScmbest  = ma.masked_invalid(a1nsurfNScmbest  ).filled(miss)
    a1nsurfMScmbest  = ma.masked_invalid(a1nsurfMScmbest  ).filled(miss)
    a2prwatprofNSest = ma.masked_invalid(a2prwatprofNSest ).filled(miss)

    #**********************************
    # Convert bottom=surface to bottom=elipsoid
    #----------------------------------
    dbPath = dbDir + '/db_%05d.bin'%(idx_db)
    if   dbtype == 'jpl':
        db.set_file(dbPath)
    elif dbtype == 'my':
        db.set_idx_db(dbDir, idx_db)
    else:
        print 'check dtype', dbtype
        sys.exit()


    vres     = 250 # m
    tmpbottom = -1000 # m
    tmpnz     = 2*NLEV_PRECIP + 4  # 104

    a1elev = db.get_var('elev')[a1irec]
    a1surfbin = ((a1elev - tmpbottom) / vres).astype('int16')
    a1surfbin = tmpnz - 1 - a1surfbin # top = 0, bottom = tmpnz -1

    ny = a2prwatprofNSest.shape[0]

    #-- Est ----
    a2tmp = np.ones([ny,tmpnz], float32) * miss
    a1y   = range(ny)
    for iz in range(NLEV_PRECIP):
        a1k = a1surfbin - NLEV_PRECIP +1 + iz
        a2tmp[a1y,a1k] = a2prwatprofNSest[:,iz]
    a2prwatprofNSest = a2tmp[:,NLEV_PRECIP:2*NLEV_PRECIP] 

    #-- Obs ----
    a2tmp = np.ones([ny,tmpnz], float32) * miss
    a1y   = range(ny)
    for iz in range(NLEV_PRECIP):
        a1k = a1surfbin - NLEV_PRECIP +1 + iz
        a2tmp[a1y,a1k] = a2prwatprofNSObs[:,iz]
    a2prwatprofNSObs = a2tmp[:,NLEV_PRECIP:2*NLEV_PRECIP] 

    #for i in range(a2tmp.shape[0]):
    #    elev= a1elev[i]
    #    print ''
    #    print i,elev, a2prwatprofNSObs[i]
    #    print i,elev, a2pr


    #**********************************
    # Save output
    #----------------------------------
    outDir = outbaseDir + '/%05d'%(idx_db)
    mk_dir(outDir)

    np.save(outDir + '/nsurfNS.est.%05d.npy'%(idx_db), a1nsurfNSest) 
    np.save(outDir + '/nsurfNScmb.est.%05d.npy'%(idx_db), a1nsurfNScmbest) 
    np.save(outDir + '/precip_water_prof_NS.est.%05d.npy'%(idx_db), a2prwatprofNSest) 

    np.save(outDir + '/nsurfNS.obs.%05d.npy'%(idx_db), a1nsurfNSObs) 
    np.save(outDir + '/nsurfNScmb.obs.%05d.npy'%(idx_db), a1nsurfNScmbObs) 
    np.save(outDir + '/precip_water_prof_NS.obs.%05d.npy'%(idx_db), a2prwatprofNSObs) 

    np.save(outDir + '/irec.obs.%05d.npy'%(idx_db), a1irec) 
    np.save(outDir + '/irec.top.%05d.npy'%(idx_db), a1topirec)
    np.save(outDir + '/idx_db.top.%05d.npy'%(idx_db), a1topidxdb)

    if sensor not in ['GMI','AMSR2']:
        np.save(outDir + '/satid.obs.%05d.npy'%(idx_db), a1satidout.astype(int16))

    print outDir + '/nsurfNS.est.%05d.npy'%(idx_db)





# %%
