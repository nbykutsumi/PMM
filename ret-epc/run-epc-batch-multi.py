# %%
import sys, os, shutil, socket
import glob
import subprocess
from datetime import datetime, timedelta
import numpy as np
import numpy.ma as ma
import h5py
import shutil

#iDTime = datetime(2014,6,1)
#eDTime = datetime(2014,11,30)
iDTime = datetime(2018,1,1)   # MRMS: 2018/1/23 -. No NOAA20
eDTime = datetime(2018,12,31)
#eDTime = datetime(2018,10,31)


dDTime = timedelta(days=1)
makeS2IDX= True
maket2m  = True
#------------
def ret_lDTime(iDTime,eDTime,dDTime):
  total_steps = int( (eDTime - iDTime).total_seconds() / dDTime.total_seconds() + 1 )
  return [iDTime + dDTime*i for i in range(total_steps)]

def ret_lYM(iYM, eYM):
  """
  iYM = [iYear, iMon], eYM = [eYear, eMon]
  """
  iYear, iMon = iYM
  eYear, eMon = eYM
  lYM = []
  for Year in range(iYear, eYear+1):
    if iYear == eYear:
      lMon = range(iMon,eMon+1)
    elif Year == iYear:
      lMon = range(iMon,12+1)
    elif Year == eYear:
      lMon = range(1,eMon+1)
    else:
      lMon = range(1,12+1)

    for Mon in lMon:
      lYM.append([Year,Mon])
  return lYM
#------------
lDTime = ret_lDTime(iDTime,eDTime,dDTime)
#useorblist = False
#useorblist = True
useorblist = 'regional'

batchsize = 6
#batchsize = 1

DB_MAXREC = 10000
DB_MINREC = 1000
dbtype    = 'JPL'
#dbtype    = 'my'

#--- Storm top parameter ----
#stoptype = 'ret'
#stoptype = 'cor'
stoptype = 'no'
#stoptype = 'obs'
if stoptype=='ret':
    stopstamp = 'best01-HTQZ-ssn0'
elif stoptype=='cor':
    stopstamp = 'best01cr-HTQZ-ssn0'
#----------------------------
gmi       = ["GPM","GMI","1C","1C","V05"]
amsr2     = ["GCOMW1","AMSR2","1C","1C","V05"]
ssmis_f16 = ["F16","SSMIS","1C","1C","V05"]
ssmis_f17 = ["F17","SSMIS","1C","1C","V05"]   # 37V is missing since Apr 2016.
ssmis_f18 = ["F18","SSMIS","1C","1C","V05"]
atms_npp  = ["NPP","ATMS","1C","1C","V05"]
atms_noaa20= ["NOAA20","ATMS","1C","1C","V05"]   # MRMS does not have NOAA20 for eny year.

mhs_metopa= ["METOPA","MHS","1C","1C","V05"]
mhs_metopb= ["METOPB","MHS","1C","1C","V05"]
mhs_noaa18= ["NOAA18","MHS","1C","1C","V05"]  # Not available at arthurhou.pps after 2018/10/21
mhs_noaa19= ["NOAA19","MHS","1C","1C","V05"]

#lspec = [amsr2, ssmis_f16, ssmis_f18, atms_npp, mhs_metopa, mhs_metopb, mhs_noaa18, mhs_noaa19]
#lspec = [amsr2, ssmis_f16, ssmis_f18, atms_npp, mhs_metopa, mhs_metopb, mhs_noaa18, mhs_noaa19]
#lspec = [mhs_noaa18, mhs_noaa19]
#lspec = [mhs_noaa18]
lspec = [mhs_noaa19]



dnscan = {'GMI':2, 'AMSR2':5, 'SSMIS':4, 'ATMS':4, 'MHS':1}
#** Constants ******
#expr = 'glb.stop01.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
#expr = 'glb.stop-wgt-obs-01.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
#expr = 'glb.stop-rng-obs-01.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
#expr = 'glb.stop-wgt-%s-01.minrec%d.maxrec%d'%(stoptype,DB_MINREC,DB_MAXREC)
#expr = 'glb.stop-rng-%s-01.minrec%d.maxrec%d'%(stoptype,DB_MINREC,DB_MAXREC)
#expr = 'glb.relsurf01.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
#expr = 'glb.relsurf02.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)
#expr = 'test'
expr = 'us.v01'

if stoptype !='no':
    type_stop = expr.split('.')[1].split('-')[1]
else:
    type_stop = ''

prog = 'ret-epc-multi.py'
#prog = 'ret-myepc-29bins.py'
#prog = 'ret-myepc-stop.py'

#*******************
# Function
#*******************
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except OSError:
    pass

#*******************
# Copy program
#*******************
copyDir = './progtemp'
mk_dir(copyDir)
progtime = datetime.now().strftime('%Y-%m-%d-%H:%M-%S-%f')
progcopy = copyDir + '/%s.%s'%(prog,progtime)
shutil.copy(prog, progcopy)
print progcopy

#*****************************************
# Start sate,sensor loop
#-----------------------------------------
for spec in lspec:
    print 'spec=',spec
    print ''
    sate      = spec[0]
    sensor    = spec[1]
    prdName   = spec[2]
    prj       = spec[3]
    ver       = spec[4]

    myhost = socket.gethostname()
    if myhost =="well":
        #tbbaseDir  = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/GPM.GMI/1C/V05'
        #tbbaseDir  = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/%s.%s/1C/%s'%(sate,sensor,ver)
        tbbaseDir  = '/media/disk2/data/PMM/NASA/%s.%s/1C/%s'%(sate,sensor,ver)
        matchbaseDir= '/media/disk2/share/PMM/MATCH.%s.%s.%s'%(sensor,sate,ver)
        rnrbaseDir  = ''
        outbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/%s/%s.%s'%(expr,sensor,sate)
        tankDir = '/home/utsumi/mnt/lab_tank'
        if dbtype=='JPL':
            coefDir   = '/media/disk2/share/PMM/JPLDB/EPC_COEF/%s'%(sensor)
            dbDir     = '/media/disk2/share/PMM/JPLDB/EPC_DB/%s_EPC_DATABASE'%(sensor)
            relprofDir= '/media/disk2/share/PMM/JPLDB/EPC_DB/%s_rs_precip_water_prof_NS'%(sensor)

        elif dbtype=='my':
            coefDir   = '/media/disk2/share/PMM/EPCDB/EPC_COEF/%s'%(sensor)
            dbDir     = '/media/disk2/share/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
            relprofDir= ''

    else:
        tbbaseDir  = '/home/utsumi/mnt/lab_work/hk02/PMM/NASA/GPM.GMI/1C/V05'
        matchbaseDir= '/media/disk2/share/PMM/MATCH.GMI.V05A'
        rnrbaseDir  = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/glb.v03.minrec1000.maxrec10000'
        coefDir = '/media/disk2/share/PMM/EPCDB/EPC_COEF/%s'%(sensor)
        dbDir   = '/media/disk2/share/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
        outbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/%s/%s.%s'%(expr,sensor,sate)
        tankDir = '/home/utsumi/mnt/lab_tank'


    #*******************
    # Make level 1C Tb path list
    #*******************
    if useorblist is True:
        orblistPath= tankDir + '/utsumi/PMM/retepc/list-orbit/list.GMI.well.201406-201505.1000obts.txt'
        #orblistPath= tankDir + '/utsumi/PMM/retepc/list-orbit/list.GMI.well.201406-201505.2obts.txt'

        f=open(orblistPath,'r'); lines=f.readlines(); f.close()
        ltbPathAll = []
        liescanAll = []
        for tbPath in lines:
            tbPath = tbPath.strip()
            Year,Mon,Day= map(int, os.path.dirname(tbPath).split('/')[-3:])
            DTimeTmp = datetime(Year,Mon,Day)
            if DTimeTmp < iDTime: continue
            if DTimeTmp > eDTime: continue
            ltbPathAll.append(tbPath)
            liescanAll.append([-9999,-9999]) 

    elif useorblist is 'regional':
        iY,iM = iDTime.timetuple()[:2]
        eY,eM = eDTime.timetuple()[:2]
        lYM   = ret_lYM([iY,iM],[eY,eM])
        ltbPathAll = []
        liescanAll = []
        for (Year,Mon) in lYM:
            listPath=tankDir + '/utsumi/PMM/US/obtlist/overpass.%s.%s.%04d.%02d.csv'%(sensor,sate,Year,Mon)
            f=open(listPath,'r'); lines=f.readlines(); f.close()
            for line in lines:
                _,_,Day,oid,iscan,escan = map(int,line.strip().split(','))

                if (sensor=='GMI')&(oid<=24473): continue
                #if (sensor=='AMSR2')&(oid <= 31235): continue # test
                #if (sate=='F16')&(oid <= 74768): continue # test
                #if (sate=='METOPB')&(oid <= 28746): continue # test
                if (sate=='NOAA18')&(oid <= 66711): continue # test
                if (sate=='NOAA18')&(datetime(2018,10,22)<=datetime(Year,Mon,Day)): continue
                if (sate=='NOAA19')&(oid <= 47550): continue # test

                if (iDTime<=datetime(Year,Mon,Day))&(datetime(Year,Mon,Day)<=eDTime):
                    ltbPathTmp = sorted(glob.glob(tbbaseDir + '/%04d/%02d/%02d/*.%06d.????.HDF5'%(Year,Mon,Day,oid)))
                    ltbPathAll = ltbPathAll + ltbPathTmp
                    liescanAll.append([iscan,escan]) 

    else:
        ltbPathAll = []
        liescanAll = []
        for DTime in lDTime:
            Year,Mon,Day = DTime.timetuple()[:3]
            ltbPathTmp = sorted(glob.glob(tbbaseDir + '/%04d/%02d/%02d/1C.GPM.GMI.XCAL2016-C.*.??????.????.HDF5'%(Year,Mon,Day)))
            ltbPathAll = ltbPathAll + ltbPathTmp
            liescanAll.append([-9999,-9999]) 
    #*******************
    # Make batch list 
    #*******************
    lltbPath = [ltbPathAll[i*batchsize:(i+1)*batchsize] for i in range(int(len(ltbPathAll)/batchsize) +1)]
    lliescan = [liescanAll[i*batchsize:(i+1)*batchsize] for i in range(int(len(liescanAll)/batchsize) +1)]


    ##-- test --
    #lltbPath = [['/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.GMI/1C/V05/2014/07/01/1C.GPM.GMI.XCAL2016-C.20140701-S121951-E135224.001927.V05A.HDF5']]
    ##----------
    icount = -1
    for ibatch,(ltbPath,liescan) in enumerate(zip(lltbPath, lliescan)):
        #** Loop in each batch ***
        nscan= dnscan[sensor]
        lscan= range(1,nscan+1)
        d3tb = {i:[] for i in lscan}
        a3tb1 = []
        a3tb2 = []
        a2lat = []
        a2lon = []
        a2inc = [] 
        a2x   = []
        a2y   = []
        a2t2m = []
        a2tqv = []
        a2elev= []
        a2stop= []
        a2rnr = []
        loid  = []
        lny   = [] 
        lymd  = []

        yoffset = 0
        for tbPath in ltbPath:
            icount = icount+1
            iscan,escan = liescanAll[icount]
            if iscan<0:
                iscan,escan = None, None
            else:
                escan = escan+1

            oid = int(tbPath.split('.')[-3])
            Year,Mon,Day = map(int, os.path.dirname(tbPath).split('/')[-3:])
            print 'oid=',oid
            #if oid <=2780: continue  # test

            #*****************************
            # Make GMI S2 position index data
            #*****************************
            if (sensor=='GMI')&(makeS2IDX is True):
                progtmp = './mk.match.idx.gmiS2.gmi.fullswath.py'

                dtmp = {}
                dtmp['Year']    = Year
                dtmp['Mon']     = Mon
                dtmp['Day']     = Day
                dtmp['pmwPath'] = tbPath
                dtmp['obaseDir']= matchbaseDir + '/S1.ABp000-220.GMI.S2.IDX'

                dydxDir = matchbaseDir + '/S1.ABp000-220.GMI.S2.dydx'
                dtmp['dyPath0'  ]= dydxDir + '/dy.000.npy'
                dtmp['dxPath0'  ]= dydxDir + '/dx.000.npy'
                dtmp['dyPath180']= dydxDir + '/dy.180.npy'
                dtmp['dxPath180']= dydxDir + '/dx.180.npy' 

                sargv = ['%s=%s'%(key, dtmp[key]) for key in dtmp.keys()]
                sargv = ' '.join(sargv)
                lcmd = ['python', progtmp, sargv]
                print lcmd
                returncode = subprocess.call(lcmd)
                if returncode==1:
                    print 'Error for',lcmd
                    sys.exit()



            #*****************************
            # Make T2m data
            #*****************************
            if maket2m is True:
                progtmp = './extract.merra2.orbit.multi.py'    
            
                dtmp= {}
                dtmp['sate']   = sate
                dtmp['sensor'] = sensor
                dtmp['rabaseDir'] = '/home/utsumi/mnt/lab_tank/utsumi/data/MERRA2'
                dtmp['pmwPath']   = tbPath
                dtmp['varName'] = 't2m'
                dtmp['Year']    = Year
                dtmp['Mon']     = Mon
                dtmp['Day']     = Day
                sargv = ['%s=%s'%(key, dtmp[key]) for key in dtmp.keys()]
                sargv = ' '.join(sargv)
                lcmd = ['python', progtmp, sargv]
                print lcmd
                returncode = subprocess.call(lcmd)
                if returncode==1:
                    print 'Error for',lcmd
                    sys.exit()

            #*****************************

            lymd.append([Year,Mon,Day])
            loid.append(oid)
            #------------
            srcPath = tbPath
            if sensor=='GMI':
                s2xPath = glob.glob(matchbaseDir + '/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d/Xpy.1.%06d.npy'%(Year,Mon,Day,oid))[0]
                s2yPath = glob.glob(matchbaseDir + '/S1.ABp000-220.GMI.S2.IDX/%04d/%02d/%02d/Ypy.1.%06d.npy'%(Year,Mon,Day,oid))[0]
            else:
                s2xPath = ''
                s2yPath = ''
            #t2mPath = glob.glob(matchbaseDir + '/S1.ABp000-220.MERRA2.t2m/%04d/%02d/%02d/t2m.%06d.npy'%(Year,Mon,Day,oid))[0]
            t2mPath = ''
            #tqvPathTmp = glob.glob(matchbaseDir + '/S1.ABp000-220.MERRA2.tqv/%04d/%02d/%02d/tqv.%06d.npy'%(Year,Mon,Day,oid))[0]
            tqvPath = ''
            #cnvfrcinPath= ''
            #elevPathTmp = glob.glob(matchbaseDir + '/S1.ABp000-220.gtopo/%04d/%02d/%02d/gtopo.%06d.npy'%(Year,Mon,Day,oid))[0]
            elevPath = ''
            try:
                if stoptype=='obs':
                    stopPath = glob.glob(matchbaseDir + '/S1.ABp083-137.Ku.V06A.9ave.heightStormTop/%04d/%02d/%02d/heightStormTop.%06d.npy'%(Year,Mon,Day,oid))[0]
                elif stoptype in ['ret','cor']:
                    stopPath = glob.glob(tankDir + '/utsumi/PMM/stop/orbit/%s/%04d/%02d/%02d/stop.%06d.npy'%(stopstamp,Year,Mon,Day,oid))[0]
                elif stoptype == 'no':
                    stopPath = ''

            except IndexError:
                print 'skip oid=',oid
                lymd.pop(-1)
                loid.pop(-1)
                continue
            #rnrPath=glob.glob(rnrbaseDir + '/%04d/%02d/%02d/nsurfNScmb.%06d.y-9999--9999.nrec10000.npy'%(Year,Mon,Day,oid))[0]
            rnrPath = ''

            #-- Append HDF file contents-----
            with h5py.File(srcPath,'r') as h: 
                for scan in lscan:
                    d3tb[scan].append( h['/S%s/Tc'%(scan)][iscan:escan])
                a2lat.append( h['/S1/Latitude'][iscan:escan])
                a2lon.append( h['/S1/Longitude'][iscan:escan])
                a2inc.append( h['/S1/incidenceAngle'][iscan:escan])
            #-- Append npy files -----
            if s2xPath !='':
                a2xTmp = np.load(s2xPath).astype('int32')[iscan:escan]
                a2yTmp = np.load(s2yPath).astype('int32')[iscan:escan]
                a2yTmp = (ma.masked_outside(a2yTmp,iscan,escan-1) -iscan + yoffset ).filled(-9999)  # add yoffset

                a2x.append(a2xTmp)
                a2y.append(a2yTmp)

            if t2mPath !='':
                a2t2m.append(np.load(t2mPath)[iscan:escan])
            if tqvPath !='':
                a2tqv.append(np.load(tqvPath)[iscan:escan])
            if elevPath !='':
                a2elev.append(np.load(elevPath)[iscan:escan])
            if stopPath !='':
                a2stop.append(np.load(stopPath)[iscan:escan])
            if rnrPath !='':
                a2rnr.append(np.load(rnrPath)[iscan:escan])

            #-- Renew yoffset ---------------
            yoffset = yoffset + d3tb[1][-1].shape[0]

            #-- Append NREC -----------------
            lny.append(d3tb[1][-1].shape[0])
            #--------------------------------

        if len(loid)==0: continue

        #-- Concatenate data ----
        for scan in lscan:
            d3tb[scan] = np.concatenate(d3tb[scan], axis=0)
        a2lat = np.concatenate(a2lat, axis=0)
        a2lon = np.concatenate(a2lon, axis=0) 
        a2inc = np.concatenate(a2inc, axis=0) 

        if s2xPath !='':
            a2x   = np.concatenate(a2x, axis=0)
            a2y   = np.concatenate(a2y, axis=0)
        if t2mPath !='':
            a2t2m = np.concatenate(a2t2m,axis=0)
        if tqvPath !='':
            a2tqv  = np.concatenate(a2tqv, axis=0) 
        if elevPath !='':
            a2elev = np.concatenate(a2elev,axis=0)
        if stopPath !='':
            a2stop = np.concatenate(a2stop,axis=0)
            nytmp,nxtmp = a2stop.shape
            print a2stop.shape
            if (sensor=='GMI')and(nxtmp != 221):
                a2stopTmp = (np.ones([nytmp,221])*(-9999.)).astype('int32')
                #a2stopTmp[:,103:117+1] = a2stop
                a2stopTmp[:,83:137+1] = a2stop
                a2stop = a2stopTmp
        if rnrPath !='':
            a2rnr=np.concatenate(a2rnr, axis=0)

        #----------------------- ----
        # Make temporary directory
        #----------------------------
        oid0 = ltbPath[0].split('.')[-3]
        oid1 = ltbPath[-1].split('.')[-3]
        #outDirTmp  = outbaseDir + '/%04d/%02d/%02d/temp'%(Year,Mon,Day)
        outDirTmp  = outbaseDir + '/tempfiles/tmp.%s.%s-%s'%(progtime,oid0,oid1)
        mk_dir(outDirTmp)

        dargv = {}
        #-- Batch file names ----
        oid = loid[0]
        srcPathTmp = outDirTmp + '/tb.%06d.HDF5'%(oid)

        if s2xPath !='':
            s2xPathTmp = outDirTmp + '/Xpy.%06d.npy'%(oid)
            s2yPathTmp = outDirTmp + '/Ypy.%06d.npy'%(oid)
        else:
            s2xPathTmp = ''
            s2yPathTmp = ''
        if t2mPath !='':
            t2mPathTmp = outDirTmp + '/t2m.%06d.npy'%(oid)
        else:
            t2mPathTmp = ''
        if tqvPath !='':
            tqvPathTmp = outDirTmp + '/tqv.%06d.npy'%(oid)
        else:
            tqvPathTmp = ''
        if elevPath !='':
            elevPathTmp= outDirTmp + '/elev.%06d.npy'%(oid)
        else:
            elevPathTmp = ''
        #if cnvfrcinPath !='':
        #    cnvfrcinPathTmp = outDirTmp + '/cnvfrc.%06d.npy'
        #else:
        #    cnvfrcinPathTmp = ''
        if stopPath !='':
            stopPathTmp= outDirTmp + '/stop.%06d.npy'%(oid)
        else:
            stopPathTmp = ''
        if rnrPath !='':
            rnrPathTmp= outDirTmp + '/rnr.%06d.npy'%(oid)
        else:
            rnrPathTmp = ''


        #-- Save HDF file -------
        with h5py.File(srcPathTmp, 'w') as h:
            for scan in lscan:
                h.create_group('S%d'%(scan))
                h.create_dataset('S%d/Tc'%(scan), data=d3tb[scan])

            h.create_dataset('S1/Latitude', data=a2lat)
            h.create_dataset('S1/Longitude',data=a2lon)
            h.create_dataset('S1/incidenceAngle',data=a2inc)
            h.flush()

        #-- Save npy files ------
        if s2xPath !='':
            np.save(s2xPathTmp, a2x)
            np.save(s2yPathTmp, a2y)
        if t2mPath !='':
            np.save(t2mPathTmp, a2t2m)
        if tqvPath !='':
            np.save(tqvPathTmp, a2tqv)
        if elevPath !='':
            np.save(elevPathTmp, a2elev)
        if stopPath !='':
            np.save(stopPathTmp, a2stop)
        if rnrPath !='':
            np.save(rnrPathTmp, a2rnr)


        #** Save program ******
        progDir = outbaseDir + '/prog'
        mk_dir(progDir)
        stime    = datetime.now().strftime('%Y.%m.%d_%H:%M:%S')

        iprog    = os.path.basename(__file__)
        oprog    = progDir + '/%s.%s'%(iprog, stime)
        shutil.copyfile(iprog, oprog)
        print oprog  

        iprog    = os.path.basename(prog)
        oprog    = progDir + '/%s.%s'%(iprog, stime)

        if ibatch==0:
            shutil.copyfile(iprog, oprog)
            print oprog  

        #***** Set parameter dictionary *************
        #------------
        dargv['dbtype']     = dbtype
        dargv['sate']       = sate
        dargv['sensor']     = sensor
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
        dargv['DB_MAXREC'] = DB_MAXREC
        dargv['DB_MINREC'] = DB_MINREC
        dargv['DB_USE_MINREC'] = 2
        dargv['NDB_EXPAND'] = 30
        dargv['DB_RAINFRAC'] = 0.0001 # minimum fraction of precipitating events (>=1mm/h) in the DB required for retrieval
        dargv['MAX_INC_DIFF'] = 20  # degrees
        dargv['MAX_T2M_DIFF'] = 10
        dargv['MAX_TQV_DIFF'] = 10
        dargv['MAX_STOP_DIFF'] = 2000  # [meter]
        dargv['MAX_TB_RMSD'] = -9999 # max TB RMS difference between observed and DB pixel 
        dargv['MAX_RMA_0'] = -9999.  # Maximum Ratio of missing amount (0>=mm/h) acceptable for rain / no-rain classification # -9999. --> No screening

        #dargv['MIN_RNR'] = 0.2  # [mm/h] for Rain/No-rain based on first guess precipitation
        dargv['MIN_RNR'] = 0.1  # [mm/h] for Rain/No-rain based on first guess precipitation
        dargv['STD_STORMTOP'] = 2100 # [m] stop

        dargv['flag_crosstrack'] = {'GMI':0, 'AMSR2':0, 'SSMIS':0, 'ATMS':1, 'MHS':1, 'SAPHIR':1}[sensor]        
        dargv['flag_top_var'] = 0  # 0: No top-ranked vars. 1: Retrieve top-ranked vars.
        dargv['flag_rel_surf'] = 1  # 0: Not relative to surface 1: Profiles that are relative to surface
        dargv['flag_out_conv'] = 1  # 0: Not retrieve convective

        dargv['type_stop'] = type_stop  # stop


        #------------
        dargv['coefDir']    = coefDir
        dargv['dbDir']      = dbDir
        dargv['relprofDir'] = relprofDir
        dargv['outDir']     = outDirTmp

        dargv['oid'] = oid    # have to be list?
        dargv['clat'] = -9999
        dargv['clon'] = -9999
        dargv['dlatlon'] = -9999
        dargv['iscan'] = -9999
        dargv['escan'] = -9999
        dargv['dscan'] = -9999

        dargv['srcPath'] = srcPathTmp
        dargv['s2xPath'] = s2xPathTmp
        dargv['s2yPath'] = s2yPathTmp
        dargv['t2mPath'] = t2mPathTmp
        dargv['tqvPath'] = tqvPathTmp
        dargv['elevPath']= elevPathTmp
        #dargv['cnvfrcinPath']=cnvfrcinPathTmp
        dargv['stopPath']= stopPathTmp
        dargv['rnrPath'] = rnrPathTmp 
        #-------------
        sargv = ['%s=%s'%(key, dargv[key]) for key in dargv.keys()]
        sargv = ' '.join(sargv)


        lcmd = ['python', progcopy, sargv]
        print lcmd
        returncode = subprocess.call(lcmd)
        if returncode==1:
            print 'Error for',lcmd
            sys.exit()
        print returncode
        #***************************
        # Split output files 
        #***************************
        DB_MAXREC = dargv['DB_MAXREC']
        ssearch = outDirTmp + '/*.nrec%d.npy'%(DB_MAXREC)
        ltmpPath = glob.glob(ssearch)

        print ltmpPath
        print len(ltmpPath)
        print ''
        print liescan
        for tmpPath in ltmpPath:
            varName = os.path.basename(tmpPath).split('.')[0]

            print varName
            atmp = np.load(tmpPath)
            y0 = 0
            for i in range(len(loid)):
                iscan,escan = liescan[i]
                Year,Mon,Day = lymd[i]
                oid = loid[i]
                ny  = lny[i]
                aout= atmp[y0:y0+ny]
                y0  = y0+ny
                stamp = '%06d.y%04d-%04d.nrec%d'%(oid, iscan, escan,DB_MAXREC)

                outDir  = outbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
                mk_dir(outDir)
                outPath = outDir + '/%s.%s.npy'%(varName, stamp)
                np.save(outPath, aout)
                print outPath


        #***************************
        # Clean temporaty files
        #***************************
        ssearch = outDirTmp + '/*'
        ltmpPath = glob.glob(ssearch)
        for tmpPath in ltmpPath:
            print 'REMOVE    ',os.path.exists(tmpPath), tmpPath

            os.remove(tmpPath)
        os.rmdir(outDirTmp)
        print ''

        #sys.exit()  # test 
    #os.remove(progcopy)
    #print 'remove progcopy',progcopy



# %%
