# %%
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
from numpy import *
import h5py
import sys, os
import numpy as np
import JPLDB
import EPCDB
from bisect import bisect_left
import epcfunc

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
        #line = map(type, line)
        lout.append(line)
    return array(lout).astype(type)


def read_nrain(idx_db):
    '''
     /* first six= Ntotal then Nrain exceeding 1 5 10 20 50 mm/hr when T2m < 278K */
     /* second six= same for when T2m > 278K */
    '''
    if dbtype=='my':
        nrainDir= dbDir + '/nrain'
    elif dbtype=='JPL':
        nrainDir= dbDir 
    else:
        print('stop in read_nrain: dbtype=', dbtype)
        sys.exit()

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
    lRMA    = list(map(float, lines[4].split('\t')[1:]))
    ave_org   = np.array(list(map(float, lines[5].split('\t')[1:])))  # for normalization
    ave_rain  = np.array(list(map(float, lines[6].split('\t')[1:])))  # Averages of normalized variables
    ave_no    = np.array(list(map(float, lines[7].split('\t')[1:])))  # Averages of normalized variables
    std_org   = np.array(list(map(float, lines[8].split('\t')[1:])))  # for normalization
    std_rain  = np.array(list(map(float, lines[9].split('\t')[1:])))  # Averages of normalized variables
    std_no    = np.array(list(map(float, lines[10].split('\t')[1:])))  # Averages of normalized variables


    return rnrflag, thD, SR, WER, lRMA, ave_org, ave_rain, ave_no, std_org, std_rain, std_no


def ret_domain_cy(a2lat, a2lon, clat, clon, dlatlon):
    nyTmp, nxTmp = a2lat.shape
    a1lat = a2lat[:,int(nxTmp/2)]
    a1lon = a2lon[:,int(nxTmp/2)]
    
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
        print('No matching scans in the target domain are found.')
        print('Exit')
        sys.exit()


def pickyx_low2high_amsr2(nyhi):
    nxhi  = 486
    axpick = ((np.arange(nxhi) -1)/2).astype('int16')
    axpick[0] = 0    # [0,0,0,1,1,2,2,3,3,...]

    a2xpick, a2ypick = np.meshgrid(axpick, np.arange(nyhi))
    return a2ypick, a2xpick

def pickyx_low2high_ssmis(nyhi):
    nxhi = 180
    axpick = (np.arange(nxhi) /2).astype('int16')  # [0,0,1,1,2,2,3,3,...,89,89]

    a2xpick, a2ypick = np.meshgrid(axpick, np.arange(nyhi))
    return a2ypick, a2xpick

def pickyx_high2low_amsr2(nyhi):
    nxlw  = 243
    axpick = (np.arange(nxlw)*2+2).astype('int16')
    axpick[-1] = nxlw*2-1   # [2,4,6,...482,484,485]

    a2xpick, a2ypick = np.meshgrid(axpick, np.arange(nyhi))
    return a2ypick, a2xpick

def pickyx_high2low_ssmis(nyhi):
    nxlw = 90
    axpick = (np.arange(nxlw)*2+1).astype('int16') # [1,3,5,...,177,179]

    a2xpick, a2ypick = np.meshgrid(axpick, np.arange(nyhi))
    return a2ypick, a2xpick

#**************************************************************
# main
#--------------------------------------------------------------
dtb_compare= {'GMI'    :[1,1,1,1,1,1,1,1,1,1,1,1,1],
              'AMSR2'  :[1,1,1,1,1,1,1,1,1,1],
              'SSMIS'  :[1,1,1,1,1,1,1,0,0,0,0],
              'ATMS'   :[1,1,1,1,1,1,1,1,1],
              'MHS'    :[1,1,1,1,1],
              'SAPHIR' :[1,1,1,1,1,1],
            } 

argvs = sys.argv
if len(argvs)==1:
    print('***********************************')
    print('')
    print('No standard Input')
    print('Default files and parameters are used')
    print('')
    print('***********************************')
    #** Constants ******
    dbtype= 'JPL'
    #dbtype= 'my'

    #sate    = 'NOAA20'
    #sensor  = 'ATMS'
    sate    = 'GCOMW1'
    sensor  = 'AMSR2'

    coefDir = '/media/disk2/share/PMM/JPLDB/JPLDB/EPC_COEF/%s'%(sensor)
    miss_int32= np.int32(-9999)
    miss_int16= np.int16(-9999)

    if dbtype=='JPL':  
        #dbDir   = '/work/hk01/utsumi/JPLDB/EPC_DB/GMI_EPC_DATABASE_TEST29'
        ##nrainDIr= '/work/hk01/utsumi/JPLDB/EPC_DB/GMI_EPC_DATABASE_TEST29'
        #coefDir = '/tank/utsumi/PMM/JPLDB/EPC_COEF/%s'%(sensor)

        dbDir   = '/media/disk2/share/PMM/JPLDB/EPC_DB/%s_EPC_DATABASE'%(sensor)
        coefDir = '/media/disk2/share/PMM/JPLDB/EPC_COEF/%s'%(sensor)
        relprofDir = '/media/disk2/share/PMM/JPLDB/EPC_DB/%s_rs_precip_water_prof_NS'%(sensor) 


    elif dbtype=='my':
        dbDir   = '/work/hk01/utsumi/PMM/EPCDB/samp.20000.GMI.V05A.S1.ABp103-117.01-12' 
        #nrainDir= '/work/hk01/utsumi/PMM/EPCDB/samp.20000.GMI.V05A.S1.ABp103-117.01-12/nrain' 
        coefDir = '/media/disk2/share/PMM/JPLDB/EPC_COEF/%s'%(sensor)
        relprofDir = '' 

    #-- Single Run ---------
    flag_crosstrack = {'GMI':0, 'AMSR2':0, 'SSMIS':0, 'ATMS':1, 'MHS':1, 'SAPHIR':1}[sensor] 
    flag_top_var = 0
    flag_out_conv= 1
    flag_rel_surf= 1
    flag_crosstrack=1

    nscan = {'GMI':2, 'AMSR2':5, 'SSMIS':4, 'ATMS':4, 'MHS':1}[sensor]

    mainscan = {'GMI':1, 'AMSR2':1, 'SSMIS':1, 'ATMS':1, 'MHS':1, 'SAPHIR':1}[sensor]
    NEM     = {'GMI':12, 'AMSR2':13, 'SSMIS':10, 'ATMS':6, 'MHS':4, 'SAPHIR':4}[sensor]
    NTBREG  = {'GMI':13, 'AMSR2':10, 'SSMIS':11, 'ATMS':9, 'MHS':5, 'SAPHIR':6}[sensor]

    NEM_USE = 3
    NPCHIST = 29
    if dbtype=='JPL':
        NLEV_DPR    = 88  # extract this number of layers
        NLEV_PRECIP = 22
    elif dbtype=='my':
        NLEV_DPR    = 50  # extract this number of layers
        NLEV_PRECIP = 50


    MAX_INC_DIFF= 20  # degree
    MAX_T2M_DIFF= 10  # K
    MAX_TQV_DIFF= 10  # kg/m2
    MAX_TB_RMSD = -9999.
    MAX_RMA_0   = -9999.  # Maximum Ratio of missing amount (0>=mm/h) acceptable for rain / no-rain classification # -9999. --> No screening
    MIN_RNR     = -9999. # 2020/02/07 If the rain/no-rain index exceeds this number, do retrieval
    DB_MAXREC   = 10000
    DB_MINREC    = 2000
    DB_USE_MINREC= 2
    NDB_EXPAND  = 30
    DB_RAINFRAC = 0.01  # minimum fraction of precipitating events (>=1mm/h) in the DB required for retrieval
    thwtmin = 0.1
    miss    = -9999.


    clat, clon = (35,150) #GCOMW1.AMSR2 2018/1/1 #029920
    #clat, clon = (0,20) #F16.SSMIS 2018/1/1 #073295
    #clat, clon = (5,105) #METOPA.MHS 2018/1/1 #058125
    #clat, clon = (10,-170) #NOAA20.ATMS 2018/1/1 #000620
    dlatlon = 3  # used to search the domain center
    iscan   = -9999
    escan   = -9999
    dscan   = 90   # set clat and iscan =-9999 for entire orbit
    #dscan   = 55
    #dscan   = 5


    # GCOMW1.AMSR2 2018/1/1 #029920
    if (sate=='GCOMW1')&(sensor=='AMSR2'):
        srcPath = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/GCOMW1.AMSR2/1C/V05/2018/01/01/1C.GCOMW1.AMSR2.XCAL2016-V.20180101-S023755-E041647.029920.V05A.HDF5'

    # F16.SSMIS 2018/1/1 #073295
    elif (sate=='F16')&(sensor=='SSMIS'):
        srcPath = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/F16.SSMIS/1C/V05/2018/01/01/1C.F16.SSMIS.XCAL2016-V.20180101-S012148-E030342.073295.V05B.HDF5'

    elif (sate=='NOAA20')&(sensor=='ATMS'):
        srcPath = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/NOAA20.ATMS/1C/V05/2018/01/01/1C.NOAA20.ATMS.XCAL2018-V.20180101-S001651-E015819.000620.V05A.HDF5'

    elif (sate=='METOPA')&(sensor=='MHS'):
        srcPath = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/METOPA.MHS/1C/V05/2018/01/01/1C.METOPA.MHS.XCAL2016-V.20180101-S011542-E025703.058125.V05A.HDF5'

    print(srcPath)
    s2xPath= ''
    s2yPath= ''
    #t2mPath = '/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.GCOMW1.AMSR2.V05/S1.ABp000-242.MERRA2.t2m/2018/01/01/t2m.029920.npy'
    t2mPath = ''
    tqvPath = ''
    elevPath = ''
    rnrPath = ''
    outDir = '/home/utsumi/temp/out/%s/%s.%s'%(dbtype,sensor,sate)
    oid    = int(srcPath.split('.')[-3])

    ## 2016/4/18 QJRMS GMI=012149 **
    #srcPath = '/work/hk01/PMM/NASA/GPM.GMI/1C/V05/2016/04/18/1C.GPM.GMI.XCAL2016-C.20160418-S115529-E132803.012149.V05A.HDF5'
    #s2xPath= '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/2016/04/18/Xpy.1.012149.npy'
    #s2yPath= '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.GMI.S2.IDX/2016/04/18/Ypy.1.012149.npy'
    #t2mPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.MERRA2.t2m/2016/04/18/t2m.012149.npy'
    #elevPath = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp000-220.gtopo/2016/04/18/gtopo.012149.npy'
    #


#**************************************************************
# Standard input
#**************************************************************
elif len(argvs)>2:
    print('Too many standard input')
    print('Python [prog] [parameter-string]')
    sys.exit()

elif (len(argvs)==2)&(argvs[1]=='test'):
    dbtype     = 'my'
    sate       = 'GPM'
    sensor     = 'GMI'
    nscan   ={'GMI':2, 'AMSR2':5, 'SSMIS':4, 'ATMS':4, 'MHS':1}[sensor]   # number of scan classes (S1, S2, ..)

    mainscan={'GMI':1, 'AMSR2':1, 'SSMIS':1, 'ATMS':1, 'MHS':1, 'SAPHIR':1}[sensor]
    NEM    = {'GMI':12, 'AMSR2':13, 'SSMIS':10, 'ATMS':6, 'MHS':4, 'SAPHIR':4}[sensor]
    NTBREG = {'GMI':13, 'AMSR2':10, 'SSMIS':11, 'ATMS':9, 'MHS':5, 'SAPHIR':6}[sensor]
    NEM_USE= 3
    NPCHIST= 29
    NLEV_DPR = 50 # extract this number of layers
    NLEV_PRECIP = 50
    thwtmin = 0.1
    miss    = -9999.
    miss_int32= np.int32(-9999)
    miss_int16= np.int16(-9999)

    DB_MAXREC = 10000
    DB_MINREC = 1000
    DB_USE_MINREC = 2
    NDB_EXPAND= 30
    DB_RAINFRAC = 0.0001 # minimum fraction of precipitating events (>=1mm/h) in the DB required for retrieval
    MAX_INC_DIFF= 20
    MAX_T2M_DIFF= 10
    MAX_TQV_DIFF= 10
    MAX_TB_RMSD = -9999
    MAX_RMA_0   = -9999
    MIN_RNR     = 0.1  # 2020/02/07 If the rain/no-rain index exceeds this number, do retrieval
    flag_crosstrack= {'GMI':0, 'AMSR2':0, 'SSMIS':0, 'ATMS':1, 'MHS':1, 'SAPHIR':1}[sensor]
    flag_top_var   = 0
    flag_rel_surf  = 1
    flag_out_conv  = 1


    coefDir    = '/media/disk2/share/PMM/EPCDB/EPC_COEF/%s'%(sensor)
    dbDir      = '/media/disk2/share/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
    relprofDir = ''


    oid    = 1453
    clat   = -9999
    clon   = -9999
    dlatlon= -9999
    iscan  = 2000
    escan  = 2050
    dscan  = -9999

    srcPath   = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/GPM.GMI/1C/V05/2014/06/01/1C.GPM.GMI.XCAL2016-C.20140601-S010534-E023807.001453.V05A.HDF5'
    s2xPath   = '/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.GMI.GPM.V05/S1.ABp000-220.GMI.S2.IDX/2014/06/01/Xpy.1.001453.npy'
    s2yPath   = '/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.GMI.GPM.V05/S1.ABp000-220.GMI.S2.IDX/2014/06/01/Ypy.1.001453.npy'
    #tsPath    = dargv['tsPath']   
    t2mPath   = '' 
    tqvPath   = ''
    elevPath  = ''
    #cnvfrcinPath= dargv['cnvfrcinPath']
    rnrPath   = ''  # 2020/02/07 rain/no-rain index file
    outDir    = '/home/utsumi/temp/ret/test'




else:
    largvs = argvs[1].split()
    dargv = {}
    for argvs in largvs:
        key,param = argvs.split('=')
        dargv[key] = param

    dbtype     = dargv['dbtype']
    sate       = dargv['sate']
    sensor     = dargv['sensor']
    coefDir    = dargv['coefDir']
    dbDir      = dargv['dbDir']
    relprofDir = dargv['relprofDir']
    oid    = int(dargv['oid'])
    clat   = float(dargv['clat'])
    clon   = float(dargv['clon'])
    dlatlon= float(dargv['dlatlon'])
    iscan  = int(dargv['iscan'])
    escan  = int(dargv['escan'])
    dscan  = int(dargv['dscan'])
    nscan  = int(dargv['nscan'])
    mainscan=int(dargv['mainscan'])
    NEM    = int(dargv['NEM'])
    NTBREG = int(dargv['NTBREG'])
    NEM_USE= int(dargv['NEM_USE']) 
    NPCHIST= int(dargv['NPCHIST'])
    NLEV_DPR = int(dargv['NLEV_DPR'])
    NLEV_PRECIP =int(dargv['NLEV_PRECIP'])
    thwtmin = float(dargv['thwtmin'])
    miss    = float(dargv['miss'])
    miss_int32= float(dargv['miss_int32'])
    miss_int16= float(dargv['miss_int16'])

    DB_MAXREC = int(dargv['DB_MAXREC']) 
    DB_MINREC = int(dargv['DB_MINREC']) 
    DB_USE_MINREC = int(dargv['DB_USE_MINREC']) 
    NDB_EXPAND= int(dargv['NDB_EXPAND'])
    DB_RAINFRAC = float(dargv['DB_RAINFRAC']) # minimum fraction of precipitating events (>=1mm/h) in the DB required for retrieval
    MAX_INC_DIFF= float(dargv['MAX_INC_DIFF'])
    MAX_T2M_DIFF= float(dargv['MAX_T2M_DIFF'])
    MAX_TQV_DIFF= float(dargv['MAX_TQV_DIFF'])
    MAX_TB_RMSD = float(dargv['MAX_TB_RMSD'])
    MAX_RMA_0   = float(dargv['MAX_RMA_0'])
    MIN_RNR     = float(dargv['MIN_RNR'])  # 2020/02/07 If the rain/no-rain index exceeds this number, do retrieval
    flag_crosstrack= int(dargv['flag_crosstrack'])
    flag_top_var   = int(dargv['flag_top_var'])
    flag_rel_surf  = int(dargv['flag_rel_surf'])   
    flag_out_conv  = int(dargv['flag_out_conv'])


    srcPath   = dargv['srcPath'] 
    s2xPath   = dargv['s2xPath'] 
    s2yPath   = dargv['s2yPath'] 
    #tsPath    = dargv['tsPath']   
    t2mPath   = dargv['t2mPath']   
    tqvPath   = dargv['tqvPath']   
    elevPath  = dargv['elevPath']
    #cnvfrcinPath= dargv['cnvfrcinPath']
    rnrPath   = dargv['rnrPath']  # 2020/02/07 rain/no-rain index file
    outDir    = dargv['outDir']


#**************************************************************
# Load database module
#**************************************************************
if dbtype =='JPL':
    db    = JPLDB.JPLDB(sensor)
elif dbtype=='my':
    db    = EPCDB.EPCDB()
else:
    print('check dbtype',dbtype)


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
a2pc_edge = a2pc_edge[:,1:]    #   Added. 2020/11/11
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
#dnscan = {'GMI':2, 'AMSR2':5, 'SSMIS':4, 'ATMS':4, 'MHS':1}
#dmainscan = {'GMI':1, 'AMSR2':1, 'SSMIS':1, 'ATMS':1, 'MHS':1}

ltb_compare= np.array(dtb_compare[sensor]).astype(np.bool)

d3tb = {}
with h5py.File(srcPath, 'r') as h5:
    for scan in range(1,nscan+1):
        d3tb[scan] = h5['/S%d/Tc'%(scan)][:]

    a2lat = h5['/S%d/Latitude'%(mainscan)][:]
    a2lon = h5['/S%d/Longitude'%(mainscan)][:]
    a2inc = h5['/S%d/incidenceAngle'%(mainscan)][:,:,0]
    nyobt = a2lat.shape[0]




#-- Matchup and Joint S1 and S2 Tb --
if sensor =='GMI':
    a3tb1    = d3tb[1]
    a3tb2org = d3tb[2]

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

elif sensor=='AMSR2':
    if mainscan==1:
        a2ypick, a2xpick = pickyx_high2low_amsr2(nyobt)
        a3tb = concatenate([d3tb[1], d3tb[2], d3tb[3], d3tb[4], d3tb[5][a2ypick, a2xpick]], axis=2) 

    elif mainscan==5:
        a2ypick, a2xpick = pickyx_low2high_amsr2(nyobt)
        a3tb = concatenate([d3tb[1][a2ypick, a2xpick],
                            d3tb[2][a2ypick, a2xpick],
                            d3tb[3][a2ypick, a2xpick],
                            d3tb[4][a2ypick, a2xpick],
                            d3tb[5]], axis=2) 
    else:
        print('check mainscan', sensor, mainscan)

elif sensor=='SSMIS':
    if mainscan==1:
        a2ypick, a2xpick = pickyx_high2low_ssmis(nyobt)
        a3tb = concatenate([d3tb[1],
                            d3tb[2],
                            d3tb[4][a2ypick, a2xpick],
                            d3tb[3][a2ypick, a2xpick],
                            ], axis=2)

    elif mainscan==4:
        a2ypick, a2xpick = pickyx_low2high_ssmis(nyobt)
        a3tb = concatenate([d3tb[1][a2ypick, a2xpick],
                            d3tb[2][a2ypick, a2xpick],
                            d3tb[4],
                            d3tb[3],
                            ], axis=2) 
    else:
        print('check mainscan', sensor, mainscan)

elif sensor in ['ATMS','MHS']:
    a3tb = concatenate([d3tb[i] for i in range(1,nscan+1)],axis=2)

else:
    print('check sensor',sensor)
    sys.exit()

#-- TB_COMPARE -------------------
a3tb = a3tb[:,:,ltb_compare]

#-- Read MERRA2 data ---------
if t2mPath !='':
    a2t2m = np.load(t2mPath)
if tqvPath !='':
    a2tqv = np.load(tqvPath)

#-- Read elevation data ---------
if elevPath !='':
    a2elev = np.load(elevPath)

##-- Read convective fraction ----
#if cnvfrcinPath !='':
#    a2cnvfrc_in = np.load(cnvfrcinPath)


#****************************************************
# Extract target domain
#----------------------------------------------------
if (iscan<0)and( (clat<-180)or(180<clat)):
    pass
else:
    #if (iscan<0)and(-180<=clat)and(clat <=180):
    if (-180<=clat)and(clat <=180):
        idx_c = ret_domain_cy(a2lat, a2lon, clat, clon, dlatlon)
        iscan = idx_c-dscan
        escan = idx_c+dscan

    a3tb  = a3tb  [iscan: escan+1]
    a2lat = a2lat [iscan: escan+1]
    a2lon = a2lon [iscan: escan+1]
    a2inc = a2inc [iscan: escan+1]

    if t2mPath !='':
        a2t2m  = a2t2m[iscan: escan+1]

    if tqvPath !='':
        a2tqv  = a2tqv[iscan: escan+1]
   
    if elevPath !='': 
        a2elev= a2elev[iscan: escan+1]

    print(a3tb.shape)

    #if cnvfrcinPath !='':
    #    a2cnvfrc_in = a2cnvfrc_in[iscan: escan+1]
#****************************************************
# Make mask data
#----------------------------------------------------

#-- missing tb ---
try:
    a2mask = np.any(ma.masked_outside(a3tb, 50, 350).mask, axis=2)
except np.AxisError:
    a2mask = np.zeros(a3tb.shape[:2], dtype=bool)

#-- Check Rain/No-rain (RNR) file ---
if rnrPath !='':
    a2rnr = np.load(rnrPath)
    a2mask = ma.masked_where(a2rnr < MIN_RNR, a2mask).mask

#****************************************************
# Convert Tb to EPC
#----------------------------------------------------
print('calc epc')
#a3epc = epcfunc.mk_epc_12pc(a3tb, a2coef)
print(coefPath)
print(a3tb.shape, a2coef.shape)
a3epc = epcfunc.mk_epc(a3tb, a2coef)
print('calc epc done')

#****************************************************
# Find EPC bin numbers 
#----------------------------------------------------
print('calc idx')
#a2idx_db = epcfunc.mk_epc_id_25bins(a3epc, a2pc_edge)
a2idx_db = epcfunc.mk_epc_id_nbins(a3epc, a2pc_edge, NPCHIST)
print('calc idx done')

##--------------------------------------

a2idx_db = ma.masked_where(a2mask, a2idx_db).filled(miss)


#****************************************************
lidxset  = list(set(a2idx_db.flatten()))
lidxset  = sort(lidxset)
#****************************************************
# Initialize output variables
#----------------------------------------------------
nyout,nxout = a2idx_db.shape

#-- Initialize output array ---
a2nsurfMS    = ones([nyout,nxout],float32)*miss
a2nsurfNS    = ones([nyout,nxout],float32)*miss
a2nsurfMScmb = ones([nyout,nxout],float32)*miss
a2nsurfNScmb = ones([nyout,nxout],float32)*miss
a2cnvnsurfNScmb = ones([nyout,nxout],float32)*miss

#a3prprofNS   = ones([nyout,nxout,NLEV_PRECIP],float32)*miss
#a3prprofNScmb= ones([nyout,nxout,NLEV_PRECIP],float32)*miss
a3prwatprofNS= ones([nyout,nxout,NLEV_PRECIP],float32)*miss

a2top_idxdbMS   = ones([nyout,nxout],int32)*miss_int32
a2top_idxdbNS   = ones([nyout,nxout],int32)*miss_int32

a2top_irecMS   = ones([nyout,nxout],int32)*miss_int32
a2top_irecNS   = ones([nyout,nxout],int32)*miss_int32

if flag_top_var == 1:
    #a2top_nsurfMS = ones([nyout,nxout],float32)*miss
    #a2top_nsurfNS = ones([nyout,nxout],float32)*miss
    #a2top_nsurfMScmb = ones([nyout,nxout],float32)*miss
    #a2top_nsurfNScmb = ones([nyout,nxout],float32)*miss

    a3top_zmMS    = ones([nyout,nxout,NLEV_DPR],int32)*miss
    a3top_zmNS    = ones([nyout,nxout,NLEV_DPR],int32)*miss
    
    #a3top_prprofNS    = ones([nyout,nxout,NLEV_PRECIP],float32)*miss
    #a3top_prprofNScmb = ones([nyout,nxout,NLEV_PRECIP],float32)*miss
    a3top_prwatprofNS = ones([nyout,nxout,NLEV_PRECIP],float32)*miss
    
    a3top_tbMS      = ones([nyout,nxout,NTBREG],float32)*miss
    a3top_tbNS      = ones([nyout,nxout,NTBREG],float32)*miss

elif flag_top_var ==0:
    print('Top-ranked variables will not be retrieved')
else:
    print('check flag_top_var',flag_top_var)
    sys.exit()



##----- test ----------------
#lat, lon = 8, -173.4
#a2x, a2y = np.meshgrid(range(nxout), range(nyout))
#a2masklon = ma.masked_outside(a2lon, lon-0.1, lon+0.1).mask
#a2masklat = ma.masked_outside(a2lat, lat-0.1, lat+0.1).mask
#a2masklatlon = a2masklat + a2masklon
#lx = ma.masked_where(a2masklatlon, a2x).compressed()
#ly = ma.masked_where(a2masklatlon, a2y).compressed()
#print zip(ly,lx)
#sys.exit()
#print a2idx_db[79,32]
#sys.exit()
#-- Start retrieval --
X,Y = meshgrid(list(range(nxout)),list(range(nyout)))
print('lidxset=')
print(lidxset)
#sys.exit()   # test
for i,idx_db in enumerate(lidxset):
    print('')
    print('************************************')
    print('idx_db for primary loop =',idx_db)
    if idx_db==-9999: continue


    ##***** test **************
    #if idx_db !=3077: continue
    ##***** test **************

    #****************************************************
    a2bool = ma.masked_equal(a2idx_db, idx_db).mask
    a1x    = X[a2bool]
    a1y    = Y[a2bool]
    print(i,'/',len(lidxset), 'idx_db=%d pixels=%d'%(idx_db, len(a1x)))


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
            print('No matching database')
            print('idx_db=',idx_db_expand)
            continue
    
        #-- Read nrain file --
        try:
            a1nrain_warm, a1nrain_cold = read_nrain(idx_db_expand)
        except IOError:
            print('No DB file for idx_db=',idx_db_expand)
            print('SKIP')
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

    #******************************
    #- Check rain ratio in DB
    #------------------------------

    if ((frac0<DB_RAINFRAC)&(frac1<DB_RAINFRAC)):
        print('Insufficient rain>1 frac0=%.3f frac1=%.3f Skip idx_db=%d'%(frac0, frac1, idx_db_expand))
        print('%.3f %.3f DB_RAINFRAC=%.3f'%(frac0,frac1, DB_RAINFRAC))
        continue


    #****************************************************
    # Rain/No Rain (RNR) Screening : Conservative
    #----------------------------------------------------
    '''
    Based on the first (closest) DB in lidx_db_expand
    '''
    if rnrPath =='':
        if MAX_RMA_0 > 0:
            idx_db_tmp = lidx_db_expand[0]
     
            rnrflag, thD, SR, WER, lRMA, ave_org, ave_rain, ave_no, std_org, std_rain, std_no = read_rnr_table(idx_db_tmp)
    
    
        if MAX_RMA_0 <= 0:
            ''' No screening. Calc all pixels. '''
            pass
    
        elif rnrflag == 0:
            print('RNR screening: SKIP all pixels for idx_db = ',idx_db)
            print('Determined based on neigborhood idx_db=',idx_db_tmp)
            continue
    
        elif rnrflag==2:
            print('RNR screening: CALC all pixels for idx_db=',idx_db)
            print('Determined based on neigborhood idx_db=',idx_db_tmp)
            pass
    
        elif (rnrflag==1) and (lRMA[0] > MAX_RMA_0):
            print('RNR screening: CALC all pixels for idx_db=',idx_db)
            print('Determined based on neigborhood idx_db=',idx_db_tmp)
            print('RMA(>=0mm/h) is too large for screening: RMA(>=0mm/h)=',lRMA[0])
            pass 
    
    
        elif rnrflag==1:  # Do screening
            a2epcTmp = a3epc[a1y,a1x,:]
            a1t2mTmp = a2t2m[a1y,a1x]
            a2epcTmp = np.concatenate([a2epcTmp, a1t2mTmp.reshape(-1,1)],axis=1)
            a2epcTmp = (a2epcTmp - ave_org)/std_org
    
            a1discriminant = ( (ave_no - ave_rain)/(std_no + std_rain)*(ave_no - a2epcTmp) ).sum(axis=1)
            a1flag_rain = ma.masked_greater_equal(a1discriminant, thD).mask
    
            #-- If all pixels are no-rain --
            if a1flag_rain is np.bool_(False):
                print('RNR screening: SKIP all pixels for idx_db = ',idx_db)
                print('Determined based on neigborhood idx_db=',idx_db_tmp)
                print('Based on screening')
                continue
    
            #-- Keep only raining pixels --
            print('Screen: keep', len(a1flag_rain),'/',len(a1y))
            #if len(a1flag_rain) != len(a1y):
            #    sys.exit()
            a1y = a1y[a1flag_rain]
            a1x = a1x[a1flag_rain]

    #******************************
    #- Read DB (expand to neighborhood, if necessary)
    #------------------------------
    for iidx_db, idx_db_expand in enumerate(lidx_db_expand):
        #-- If idx_db == -9999 --
    
        #-- Read database file --
        #print 'set file'   
        dbPath = dbDir + '/db_%05d.bin'%(idx_db_expand)
        if   dbtype == 'JPL':
            db.set_file(dbPath)
        elif dbtype == 'my':
            db.set_idx_db(dbDir, idx_db_expand)


        #print 'set file done' 
        #print 'read DB'
        a2epcdbTmp = db.get_var('pc_emis')[:,:NEM]  # (nrec, 12)
        a1nsurfMScmbTmp = db.get_var('precip_MS_cmb')
        a1nsurfNScmbTmp = db.get_var('precip_NS_cmb')
        a1nsurfMSTmp    = db.get_var('precip_nsfc_MS')
        a1nsurfNSTmp    = db.get_var('precip_nsfc_NS')
        if MAX_TB_RMSD>0:
            a2tbdbTmp  = db.get_var('tb')


        #a2prprofNSTmp   = ma.masked_less(db.get_var('precip_prof_MS'), 0).filled(0.0)[:,-NLEV_PRECIP:]
        #a2prprofNSTmp   = db.get_var('precip_prof_NS')[:,-NLEV_PRECIP:]  # test
        #a2prprofNScmbTmp= ma.masked_less(db.get_var('precip_prof_NS_cmb'), 0).filled(0.0)[:,-NLEV_PRECIP:]

        if flag_rel_surf ==1:
            if dbtype=='my':
                a2prwatprofNSTmp = ma.masked_invalid(db.get_var('precip_water_prof_NS_relsurf')[:,-NLEV_PRECIP:]).filled(-9999.)

            elif dbtype=='JPL':
                a2prwatprofNSTmp = np.load(relprofDir + '/%05d.npy'%(idx_db_expand))[:,-NLEV_PRECIP:].astype('float32')   # Scaled by 1000

        else:
            a2prwatprofNSTmp = ma.masked_invalid(db.get_var('precip_water_prof_NS')[:,-NLEV_PRECIP:]).filled(-9999.)

        a1revdbTmp = db.get_var('rev') 

        if flag_crosstrack==1:
            a1incdbTmp = db.get_var('inc_S1')
        if t2mPath !='':
            a1t2mdbTmp = db.get_var('t2m') 
        if tqvPath !='':
            a1tqvdbTmp = db.get_var('tqv') 
        if elevPath !='':
            a1elevdbTmp= db.get_var('elev') 
        
        #if (cnvfrcinPath !='')or(flag_out_conv==1):
        if flag_out_conv==1:
            a1cnvfrcdbTmp = db.get_var('vfrac_conv_NS_cmb')

        a1idxdbTmp = np.ones(a2epcdbTmp.shape[0]).astype(int32)*idx_db_expand
        a1irecTmp  = np.arange(a2epcdbTmp.shape[0]).astype(int32)

        #print 'read DB done'
        #print 'MS: imported DB record length=%d'%(len(a1nsurfMSTmp))
        #print 'NS: imported DB record length=%d'%(len(a1nsurfMSTmp))
 
        #-- Stack data --
        if (iidx_db==0):
            a2epcdb = a2epcdbTmp
            a1nsurfMScmb = a1nsurfMScmbTmp
            a1nsurfNScmb = a1nsurfNScmbTmp
            a1nsurfMS    = a1nsurfMSTmp
            a1nsurfNS    = a1nsurfNSTmp
            if MAX_TB_RMSD>0:
                a2tbdb  = a2tbdbTmp

            #a2prprofNS    = a2prprofNSTmp
            #a2prprofNScmb = a2prprofNScmbTmp
            a2prwatprofNS = a2prwatprofNSTmp

            #a1tsdb  = a1tsdbTmp
            a1revdb = a1revdbTmp

            if flag_crosstrack==1:
                a1incdb = a1incdbTmp
            if t2mPath !='':
                a1t2mdb = a1t2mdbTmp
            if tqvPath !='':
                a1tqvdb = a1tqvdbTmp
            if elevPath !='':
                a1elevdb= a1elevdbTmp
            #if (cnvfrcinPath !='')or(flag_out_conv==1):
            if flag_out_conv==1:
                a1cnvfrcdb= a1cnvfrcdbTmp

            a1idxdb = a1idxdbTmp
            a1irec  = a1irecTmp

        else:
            a2epcdb = concatenate([a2epcdb,  a2epcdbTmp],axis=0)
            a1nsurfMScmb = concatenate([a1nsurfMScmb,  a1nsurfMScmbTmp], axis=0)
            a1nsurfNScmb = concatenate([a1nsurfNScmb,  a1nsurfNScmbTmp], axis=0)
            a1nsurfMS    = concatenate([a1nsurfMS,     a1nsurfMSTmp], axis=0)
            a1nsurfNS    = concatenate([a1nsurfNS,     a1nsurfNSTmp], axis=0)
            if MAX_TB_RMSD>0:
                a2tbdb  = concatenate([a2tbdb,   a2tbdbTmp],axis=0)

            #a2prprofNS    = concatenate([a2prprofNS, a2prprofNSTmp], axis=0)
            #a2prprofNScmb = concatenate([a2prprofNScmb, a2prprofNScmbTmp], axis=0)
            a2prwatprofNS = concatenate([a2prwatprofNS, a2prwatprofNSTmp])

            #a1tsdb  = concatenate([a1tsdb,   a1tsdbTmp], axis=0) 
            a1revdb = concatenate([a1revdb,  a1revdbTmp], axis=0) 

            if flag_crosstrack==1:
                a1incdb  = concatenate([a1incdb,  a1incdbTmp], axis=0) 
            if t2mPath !='':
                a1t2mdb  = concatenate([a1t2mdb,  a1t2mdbTmp], axis=0) 
            if tqvPath !='':
                a1tqvdb  = concatenate([a1tqvdb,  a1tqvdbTmp], axis=0) 
            if elevPath !='':
                a1elevdb= concatenate([a1elevdb, a1elevdbTmp], axis=0) 
            #if (cnvfrcinPath !='')or(flag_out_conv==1):
            if flag_out_conv==1:
                a1cnvfrcdb= concatenate([a1cnvfrcdb, a1cnvfrcdbTmp], axis=0)

            a1idxdb = concatenate([a1idxdb,  a1idxdbTmp], axis=0)
            a1irec  = concatenate([a1irec,   a1irecTmp], axis=0)
             

    #******************************
    #-- Start loop over y,x with the same idx_db --
    #******************************
    '''
    # a1epc   : (11)
    # a3epcdb : (nrec,11) 
    '''

    #-- Discard entries from same granule (revolution) --
    a1revflag = ma.masked_not_equal(a1revdb, oid).mask

    #-- Only valid precip entries (two scans: NS & MS)--
    ### Make only 2 types(NS & MS) based on DPR  ###
    ### Share for DPR and combined               ###

    a1prflagNS1 = ma.masked_greater_equal(a1nsurfNS,0).mask
    a1prflagMS1 = ma.masked_greater_equal(a1nsurfMS,0).mask
    a1prflagNS2 = ma.masked_greater_equal(a1nsurfNScmb,0).mask
    a1prflagMS2 = ma.masked_greater_equal(a1nsurfMScmb,0).mask

    a1prflagNS = a1prflagNS1 + a1prflagNS2
    a1prflagMS = a1prflagMS1 + a1prflagMS2

    #a1prflagNS = ma.masked_greater_equal(a1nsurfNScmb,0).mask
    #a1prflagMS = ma.masked_greater_equal(a1nsurfMScmb,0).mask

    for (y,x) in zip(a1y,a1x):  # in idx_db loop
        ####***** test **************
        #if x !=98: continue
        ####***** test **************

        #-- Obs EPC --
        #print 'idx_db y x=',idx_db,y,x
        a1epc = a3epc[y,x,:]

        #-- Obs TB  (for 'MAX_TB_RMSD')--
        if MAX_TB_RMSD>0:
            a1tb  = a3tb[y,x,:]

        #********************
        # Constrain candidates from DB
        #********************
        ##-- Discard entries from same granule (revolution) --
        #a1revflag = ma.masked_not_equal(a1revdb, oid).mask

        ##-- Only valid precip entries (two scans: NS & MS)--
        #### Make only 2 types(NS & MS) based on DPR  ###
        #### Share for DPR and combined               ###

        #a1prflagNS1 = ma.masked_greater_equal(a1nsurfNS,0).mask
        #a1prflagMS1 = ma.masked_greater_equal(a1nsurfMS,0).mask
        #a1prflagNS2 = ma.masked_greater_equal(a1nsurfNScmb,0).mask
        #a1prflagMS2 = ma.masked_greater_equal(a1nsurfMScmb,0).mask

        #a1prflagNS = a1prflagNS1 + a1prflagNS2
        #a1prflagMS = a1prflagMS1 + a1prflagMS2

        ##a1prflagNS = ma.masked_greater_equal(a1nsurfNScmb,0).mask
        ##a1prflagMS = ma.masked_greater_equal(a1nsurfMScmb,0).mask

        ##-- Ts --
        #ts    = a2ts[y,x]

        #a1tsflag = ma.masked_inside( a1tsdb-ts, -MAX_T2M_DIFF, MAX_T2M_DIFF).mask

        #-- Incidence angle --
        if flag_crosstrack==1:
            inc  = a2inc[y,x]
            a1incflag = ma.masked_inside(np.abs(a1incdb) - abs(inc), -MAX_INC_DIFF, MAX_INC_DIFF).mask
        else:
            a1incflag = True

        ##-- T2m --
        if t2mPath !='':
            t2m  = a2t2m[y,x]
            a1t2mflag = ma.masked_inside( a1t2mdb-t2m, -MAX_T2M_DIFF, MAX_T2M_DIFF).mask
        else:
            a1t2mflag = True

        ##-- tqv --
        if tqvPath !='':
            tqv  = a2tqv[y,x]
            a1tqvflag = ma.masked_inside( a1tqvdb-tqv, -MAX_TQV_DIFF, MAX_TQV_DIFF).mask
        else:
            a1tqvflag = True

        ##-- Elevation --
        if elevPath !='':
            elev = a2elev[y,x]
            if elev < 500:
                a1elevflag = ma.masked_less(a1elevdb, 500).mask
            elif (500 <=elev)and(elev < 1000):
                a1elevflag = ma.masked_inside(a1elevdb, 500, 1000).mask
            elif 1000 <=elev:
                a1elevflag = ma.masked_greater_equal(a1elevdb, 1000).mask
    
            else:
                print('check elev',elev)
                sys.exit()
        else:
            a1elevflag = True
        
        #-- Screen DB candidates --
        a1flagNS   = a1prflagNS * a1incflag * a1t2mflag * a1tqvflag * a1elevflag * a1revflag
        a1flagMS   = a1prflagMS * a1incflag * a1t2mflag * a1tqvflag * a1elevflag * a1revflag



        #if ((len(a1flagNS)<DB_USE_MINREC) or (len(a1flagMS)<DB_USE_MINREC)):
        if (a1flagNS.sum()<DB_USE_MINREC):
            print('the Number of records are too small')
            print('Skip') 
            continue
        if (a1flagMS.sum()<DB_USE_MINREC):
            print('the Number of records are too small')
            print('Skip') 
            continue
        #------------------------------

        #if idx_db==3337:  # test
        #    print 'after'
        #    print len(a1flagNS),a1flagNS.sum(), DB_USE_MINREC
        #    sys.exit()

        a2epcdbMSSC = a2epcdb[a1flagMS]
        a2epcdbNSSC = a2epcdb[a1flagNS]


        if MAX_TB_RMSD>0:
            a2tbdbMSSC  = a2tbdb[a1flagMS]
            a2tbdbNSSC  = a2tbdb[a1flagNS]

        a1nsurfMScmbSC = a1nsurfMScmb[a1flagMS]
        a1nsurfNScmbSC = a1nsurfNScmb[a1flagNS]
        a1nsurfMSSC    = a1nsurfMS   [a1flagMS]
        a1nsurfNSSC    = a1nsurfNS   [a1flagNS]

        if flag_out_conv==1:
            a1cnvfrcdbSC     = a1cnvfrcdb[a1flagNS]

        #a2prprofNSSC   = a2prprofNS   [a1flagNS]
        #a2prprofNScmbSC= a2prprofNScmb[a1flagNS]
        a2prwatprofNSSC= a2prwatprofNS[a1flagNS]

        a1idxdbMSSC    = a1idxdb[a1flagMS]
        a1idxdbNSSC    = a1idxdb[a1flagNS]

        a1irecMSSC    = a1irec[a1flagMS]
        a1irecNSSC    = a1irec[a1flagNS]

        #-- RMSD (EPC) --
        a1rmsdMS = np.sqrt(np.square((a2epcdbMSSC - a1epc)/a1pc_std).sum(axis=1)/NEM)
        a1rmsdNS = np.sqrt(np.square((a2epcdbNSSC - a1epc)/a1pc_std).sum(axis=1)/NEM)

        #-- Mask by RMSD (TB) ---
        if MAX_TB_RMSD>0:
            a1rmsdTBMS= np.sqrt(np.square(a2tbdbMSSC - a1tb).mean(axis=1))
            a1rmsdTBNS= np.sqrt(np.square(a2tbdbNSSC - a1tb).mean(axis=1))
           
            a1rmsdMS = ma.masked_where( a1rmsdTBMS > MAX_TB_RMSD, a1rmsdMS)
            a1rmsdNS = ma.masked_where( a1rmsdTBNS > MAX_TB_RMSD, a1rmsdNS)

        #-- ------------------- 

        idxtopMS = np.argmin(a1rmsdMS)
        idxtopNS = np.argmin(a1rmsdNS)
        rmsd_minMS = a1rmsdMS[idxtopMS]
        rmsd_minNS = a1rmsdNS[idxtopNS]


        ##-- test use 2nd, not 1st top ----------
        #idxtopMS = np.argmin(ma.masked_less_equal(a1rmsdMS,rmsd_minMS))
        #idxtopNS = np.argmin(ma.masked_less_equal(a1rmsdNS,rmsd_minNS))
        #rmsd_minMS = a1rmsdMS[idxtopMS]
        #rmsd_minNS = a1rmsdNS[idxtopNS]

        ##-- test use 3rd, not 1st top ----------
        #idxtopMS = np.argmin(ma.masked_less_equal(a1rmsdMS,rmsd_minMS))
        #idxtopNS = np.argmin(ma.masked_less_equal(a1rmsdNS,rmsd_minNS))
        #rmsd_minMS = a1rmsdMS[idxtopMS]
        #rmsd_minNS = a1rmsdNS[idxtopNS]

        ##-- test use 4th not 1st top ----------
        #idxtopMS = np.argmin(ma.masked_less_equal(a1rmsdMS,rmsd_minMS))
        #idxtopNS = np.argmin(ma.masked_less_equal(a1rmsdNS,rmsd_minNS))
        #rmsd_minMS = a1rmsdMS[idxtopMS]
        #rmsd_minNS = a1rmsdNS[idxtopNS]


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
        if dbtype=='JPL':
            dbtopPath = dbDir + '/db_%05d.bin'%(topidxdbMS)
            db.set_file(dbtopPath)
        elif dbtype=='my':
            db.set_idx_db(dbDir, topidxdbMS)
        else:
            print('check dbtype', dbtype)
            sys.exit()

        topirec = topirecMS

        # Top ranked db entries -------------------
        if flag_top_var ==1:
            a3top_zmMS[y,x,:] = db.get_var('z_ka', nrec=1, origin=topirec).flatten()[-NLEV_DPR:] 
            a3top_tbMS[y,x,:] = db.get_var('tb', nrec=1, origin=topirec).flatten() 
    
            # Read top-db file (for NS) --
            if dbtype=='JPL':
                dbtopPath = dbDir + '/db_%05d.bin'%(topidxdbNS)
                db.set_file(dbtopPath)
            elif dbtype=='my':
                db.set_idx_db(dbDir, topidxdbNS)
            else:
                print('check dbtype', dbtype)
                sys.exit()
    
            topirec = topirecNS
    
            a3top_zmNS[y,x,:] = db.get_var('z_ku', nrec=1, origin=topirec).flatten()[-NLEV_DPR:] 
            a3top_tbNS[y,x,:] = db.get_var('tb', nrec=1, origin=topirec).flatten()
    
            #a3top_prprofNS[y,x,:] = db.get_var('precip_prof_NS', nrec=1, origin=topirecNS).flatten()[-NLEV_PRECIP:]
            #a3top_prprofNScmb[y,x,:] = db.get_var('precip_prof_NS_cmb', nrec=1, origin=topirecNS).flatten()[-NLEV_PRECIP:]
            a3top_prwatprofNS[y,x,:] = db.get_var('precip_water_prof_NS', nrec=1, origin=topirecNS).flatten()[-NLEV_PRECIP:]
    
            #a2top_nsurfMS[y,x] = db.get_var('precip_nsfc_MS', nrec=1, origin=topirecNS)
            #a2top_nsurfNS[y,x] = db.get_var('precip_nsfc_NS', nrec=1, origin=topirecNS)
    
            #a2top_nsurfMScmb[y,x] = db.get_var('precip_MS_cmb', nrec=1, origin=topirecNS)
            #a2top_nsurfNScmb[y,x] = db.get_var('precip_NS_cmb', nrec=1, origin=topirecNS)

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

        if flag_out_conv==1:
            a1cnvnsurfNScmbSC = a1nsurfNScmbSC * ma.masked_outside(a1cnvfrcdbSC,0,1)
            cnvnsurfNScmb = (a1cnvnsurfNScmbSC[a1boolwtNS] * a1wtNS).sum() / wtsumNS
            a2cnvnsurfNScmb[y,x] = cnvnsurfNScmb

        #prprofNS   = (a2prprofNSSC[a1boolwtNS] * a1wtNS.reshape(-1,1)).sum(axis=0) / wtsumNS
        #prprofNScmb= (a2prprofNScmbSC[a1boolwtNS] * a1wtNS.reshape(-1,1)).sum(axis=0) / wtsumNS
        #a3prprofNS[y,x,:]    = prprofNS
        #a3prprofNScmb[y,x,:] = prprofNScmb

        prwatprofNS   = (ma.masked_less(a2prwatprofNSSC[a1boolwtNS],0) * a1wtNS.reshape(-1,1)).sum(axis=0) / wtsumNS
        a3prwatprofNS[y,x,:]    = prwatprofNS.filled(-9999.)


        #if ((y==3)&(x==100)):
        #    print a1wtNS
        #    sys.exit()


#--- Replace nan in cnvnsurfNScmb ----
a2flagvalid= ma.masked_not_equal(a2nsurfNScmb,miss).mask
a2flagmiss = ma.masked_equal(a2nsurfNScmb,miss).mask
a2flagnan  = ma.masked_invalid(a2cnvnsurfNScmb).mask

a2flagvalid = a2flagvalid * a2flagnan
a2flagmiss  = a2flagmiss  * a2flagnan
a2cnvnsurfNScmb[a2flagvalid] = 0.0
a2cnvnsurfNScmb[a2flagmiss]  = miss

#--- Change data type of profiles ---
int16max = 32767 
if dbtype=='my':
    a3flagmax  = ma.masked_greater(a3prwatprofNS, (int16max-1)/1000.).mask
    a3prwatprofNS = (ma.masked_equal(a3prwatprofNS, miss)*1000.).filled(miss_int16).astype('int16')  # Scale by 1000
    a3prwatprofNS[a3flagmax] = np.int16(int16max)

elif dbtype=='JPL':
    a3prwatprofNS = ma.masked_invalid(ma.masked_greater(a3prwatprofNS, int16max)).filled(int16max).astype('int16')  # Already scaled by 1000


#--- save -----
mk_dir(outDir)

stamp = '%06d.y%04d-%04d.nrec%d'%(oid, iscan, escan,DB_MAXREC) 

np.save(outDir + '/nsurfMS.%s.npy'%(stamp), a2nsurfMS)
np.save(outDir + '/nsurfNS.%s.npy'%(stamp), a2nsurfNS)
np.save(outDir + '/nsurfMScmb.%s.npy'%(stamp), a2nsurfMScmb)
np.save(outDir + '/nsurfNScmb.%s.npy'%(stamp), a2nsurfNScmb)

if flag_out_conv==1:
    np.save(outDir + '/nsurfConvNScmb.%s.npy'%(stamp), a2cnvnsurfNScmb)

if flag_rel_surf==0:
    np.save(outDir + '/prwatprofNS.%s.npy'%(stamp), a3prwatprofNS)  # scaled by 1000
elif flag_rel_surf==1:
    np.save(outDir + '/prwatprofNS-rs.%s.npy'%(stamp), a3prwatprofNS)  # scaled by 1000

np.save(outDir + '/lat.%s.npy'%(stamp), a2lat)
np.save(outDir + '/lon.%s.npy'%(stamp), a2lon)

np.save(outDir + '/top-idxdbMS.%s.npy'%(stamp), a2top_idxdbMS)
np.save(outDir + '/top-idxdbNS.%s.npy'%(stamp), a2top_idxdbNS)

np.save(outDir + '/top-irecMS.%s.npy'%(stamp), a2top_irecMS)
np.save(outDir + '/top-irecNS.%s.npy'%(stamp), a2top_irecNS)
    
if flag_top_var ==1:
    np.save(outDir + '/top-zmMS.%s.npy'%(stamp), a3top_zmMS)
    np.save(outDir + '/top-zmNS.%s.npy'%(stamp), a3top_zmNS)
    
    #np.save(outDir + '/top-prprofNS.%s.npy'%(stamp), a3top_prprofNS)
    #np.save(outDir + '/top-prprofNScmb.%s.npy'%(stamp), a3top_prprofNScmb)
    np.save(outDir + '/top-prwatprofNS.%s.npy'%(stamp), a3top_prwatprofNS)
    
    np.save(outDir + '/top-tbMS.%s.npy'%(stamp), a3top_tbMS)
    np.save(outDir + '/top-tbNS.%s.npy'%(stamp), a3top_tbNS)
    
    
    #np.save(outDir + '/top-nsurfMS.%s.npy'%(stamp), a2top_nsurfMS)
    #np.save(outDir + '/top-nsurfNS.%s.npy'%(stamp), a2top_nsurfNS)
    #
    #np.save(outDir + '/top-nsurfMScmb.%s.npy'%(stamp), a2top_nsurfMScmb)
    #np.save(outDir + '/top-nsurfNScmb.%s.npy'%(stamp), a2top_nsurfNScmb)



print(outDir)
print(stamp)

# %%
