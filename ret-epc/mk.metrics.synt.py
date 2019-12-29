import numpy as np
from numpy import *
import h5py
from bisect import bisect_left
import epcfunc
import sys, os, glob, socket
import myfunc.util as util
import calendar
import pickle
import JPLDB

calcflag  = True
#calcflag  = False
#coefflag  = 'nocoef'  # 'nocoef', 'wcoef'
coefflag  = 'wcoef'  # 'nocoef', 'wcoef'
DB_MAXREC = 10000
DB_MINREC = 1000
nsample   = 1000
dbtype = 'my'
#dbtype = 'jpl'
sensor = 'GMI'
#sensor = 'AMSR2'
#sensor = 'SSMIS'
#sensor = 'ATMS'
#sensor = 'MHS'
#lrettype = ['NS','MS','NScmb','MScmb']
#lrettype = ['NS','MS','NScmb','MScmb','GPROF']
#lrettype = ['GPROF']  
#lrettype = ['MS','NS','NScmb','MScmb']
lrettype = ['NScmb']
thpr = 0.2
expr = 'org.%s.smp%d'%(sensor,nsample)
lidx_db = range(29*29*29)[1:]
#lidx_db = range(5829,5850)
#lidx_db = range(6000,6100)
#lidx_db = range(29*29*4,29*29*8)
#lidx_db = range(29*29*8,29*29*12)
#lidx_db = range(29*29*12,29*29*16)
#lidx_db = range(29*29*16,29*29*20)
#lidx_db = range(29*29*20,29*29*24)
#lidx_db = range(29*29*24,29*29*29)
#** Constants ******
myhost = socket.gethostname()
if myhost =='shui':
    if dbtype=='jpl':
        dbDir   = '/tank/utsumi/PMM/JPLDB/EPC_DB/%s_EPC_DATABASE_TEST29'%(sensor)
        countDir= '/tank/utsumi/PMM/JPLDB/list'
        retbaseDir = '/tank/utsumi/PMM/retsynt/%s'%(expr)

    elif dbtype=='my':
        dbDir   = '/work/hk01/utsumi/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
        countDir= '/work/hk01/utsumi/PMM/EPCDB/list'
        retbaseDir = '/tank/utsumi/PMM/retsynt/%s'%(expr)

elif myhost == 'well':
    if dbtype=='jpl':
        dbDir   = '/home/utsumi/mnt/lab_tank/utsumi/PMM/JPLDB/EPC_DB/%s_EPC_DATABASE_TEST29'%(sensor)
        countDir= '/home/utsumi/mnt/lab_tank/utsumi/PMM/JPLDB/list'
        retbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retsynt/%s'%(expr)
    elif dbtype=='my':
        dbDir   = '/media/disk2/share/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
        countDir= '/media/disk2/share/PMM/EPCDB/list'
        retbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retsynt/%s'%(expr)


else:
    print 'check hostname',myhost
    sys.exit()


lsurftype = ['all','ocean','vegetation','coast','snow']
#lsurftype = ['ocean']

dsurflabel={ 'all'  :'All surface'
            ,'ocean':'Class1 (Ocean)'
            ,'vegetation':'Classes 3-7 (Vegetation)'
            ,'snow':'Classes 8-11 (Snow)'
            ,'coast':'Class 13 (Coast)'
            }

dlsatid = {'GMI':[1],'AMSR2':[30], 'SSMIS':[16,17,18,19], 'ATMS':[100,101], 'MHS':[201,202,318,319]}

dsatname = {999:'ALL',0:'TRMM',1:'GPM',16:'F16',17:'F17',18:'F18',19:'F19',30:'GCOMW',100:'NPP',101:'NOAA20',201:'METOP-A',202:'METOP-B',318:'NOAA18',319:'NOAA19', 400:'SAPHIR'}
#**************************************
# Functions
#--------------------------------------
def read_orbitlist(Year,Mon):
    listPath = listDir + '/overpass.GPM.%04d.%02d.csv'%(Year,Mon)
    f=open(listPath); lines=f.readlines(); f.close()
    lout=[]
    for line in lines:
        line = map(int,line.strip().split(','))
        lout.append(line)
    return lout

#--------------------------------------
if dbtype=='jpl':
    db = JPLDB.JPLDB(sensor)


lsatid = dlsatid[sensor]
for rettype in lrettype:

    if calcflag is not True:
        continue

    #** Read count file of DB records **
    if coefflag == 'wcoef':
        d1coef = {}
        for satid in [999] + lsatid:
            if ((dbtype=='my')and(sensor=='GMI')):
                countPath = countDir +'/count.epc.csv'
            elif dbtype=='jpl':
                countPath = countDir +'/count.jpldb.%s.%d.csv'%(sensor, satid)
            else:
                print 'check dbtype,sensor',dbtype,sensor
                sys.exit()


            f=open(countPath,'r'); lines=f.readlines(); f.close()
            ncol = len(lines[0].strip().split(',')) - 1
            d1coef[satid] = np.zeros(29*29*29, float32)
            
            for line in lines[1:]:
                line =  map(int,line.strip().split(','))
                epcid= line[0]
                nall = line[1]
                if nall < nsample:
                    coef = 1.0
                else:
                    coef = float(nall)/nsample
                d1coef[satid][epcid] = coef
        
    elif coefflag =='nocoef':
        d1coef[satid] = np.ones(29*29*29, float32)
    else:
        print 'check coefflag', coefflag
        sys.exit()

    #***************************
    # Initialize
    #---------------------------
    h,m,f,c = {}, {}, {}, {}
    un, uy, ux, uyy, uxx, uxy = {}, {}, {}, {}, {}, {}
    cn, cy, cx, cyy, cxx, cxy = {}, {}, {}, {}, {}, {}

    for satid in [999] + lsatid:
        for surftype in lsurftype:
            key = satid,surftype
            h[key], m[key], f[key], c[key] = 0, 0, 0, 0
            un[key], uy[key], ux[key], uyy[key], uxx[key], uxy[key] = 0, 0, 0, 0, 0, 0
            cn[key], cy[key], cx[key], cyy[key], cxx[key], cxy[key] = 0, 0, 0, 0, 0, 0
            
    #---------------------------
    for idx_db in lidx_db:
        if dbtype=='jpl':
            dbPath = dbDir + '/db_%05d.bin'%(idx_db)
            
            if not os.path.exists(dbPath):
                continue

            db.set_file(dbPath)

        if rettype != 'GPROF':
            print 'idx_db=',idx_db
            obsPath = retbaseDir + '/%05d/nsurf%s.obs.%05d.npy'%(idx_db,rettype,idx_db)
            if not os.path.exists(obsPath):
                print 'No file'
                print obsPath
                continue

            a1irec= np.load(retbaseDir + '/%05d/%s.obs.%05d.npy'%(idx_db,'irec',idx_db)).astype(int32)
            a1obs = np.load(retbaseDir + '/%05d/nsurf%s.obs.%05d.npy'%(idx_db,rettype,idx_db))
            a1ret = np.load(retbaseDir + '/%05d/nsurf%s.est.%05d.npy'%(idx_db,rettype,idx_db))

            if sensor in ['SSMIS','ATMS','MHS']:
                a1satid = np.load(retbaseDir + '/%05d/satid.obs.%05d.npy'%(idx_db,idx_db))
    
        elif rettype =='GPROF':
            if dbtype =='my':
                obsPath = retbaseDir + '/%05d/nsurf%s.obs.%05d.npy'%(idx_db,'NScmb',idx_db)
                if not os.path.exists(obsPath):
                    print 'No file'
                    print obsPath
                    continue
            
                a1irec= np.load(retbaseDir + '/%05d/%s.obs.%05d.npy'%(idx_db,'irec',idx_db)).astype(int32)
                a1obs = np.load(retbaseDir + '/%05d/nsurf%s.obs.%05d.npy'%(idx_db,'NScmb',idx_db))
    
                a1ret = np.load(dbDir + '/surfacePrecipitation/%05d.npy'%(idx_db))[a1irec]

            elif dbtype == 'jpl':
                irecPath = retbaseDir + '/%05d/%s.obs.%05d.npy'%(idx_db,'irec',idx_db)
                if not os.path.exists(irecPath):
                    print 'No file'
                    print irecPath
                    continue

                            
                a1irec= np.load(retbaseDir + '/%05d/%s.obs.%05d.npy'%(idx_db,'irec',idx_db)).astype(int32)
                a1obs = db.get_var('precip_NS_cmb')[a1irec] 
                a1ret = db.get_var('precip_GPROF')[a1irec]

        else:
            print 'check rettype',rettype
            sys.exit()



        #***************************
        # Use irec to read database
        #---------------------------
        if   dbtype=='jpl':
            a1surftype = db.get_var('sfc_class')[a1irec]

        elif dbtype=='my':
            a1surftype = np.load(dbDir + '/surfaceTypeIndex/%05d.npy'%(idx_db))[a1irec]


        #---------------------------

        for satid in [999] + lsatid:
            if satid == 999:
                a1masksat = np.array([False])

            else:
                if sensor in ['GMI','AMSR2']:
                    a1masksat = np.array([False])
                else:
                    a1masksat = ma.masked_not_equal(a1satid, satid).mask

            for surftype in lsurftype:
                #-- Screen by surface types -------------------
                if   surftype=='all':
                    a1masksurf = np.array([False])
                elif surftype=='ocean':
                    a1masksurf = ma.masked_not_equal(a1surftype,1).mask
                elif surftype=='vegetation':
                    a1masksurf = ma.masked_outside(a1surftype,3,7).mask
                elif surftype=='snow':
                    a1masksurf = ma.masked_outside(a1surftype,8,11).mask
                elif surftype=='coast':
                    a1masksurf = ma.masked_not_equal(a1surftype,13).mask
                else:
                    print '\n'+'check surftype',surftype
                    sys.exit() 
    
                #-- Screeen no-precip cases for both datasets--
                a1mask1 = ma.masked_less(a1ret, 0).mask
                a1mask2 = ma.masked_less(a1obs, 0).mask
                a1mask  = a1mask1 * a1mask2
                a1mask  = a1mask + a1masksurf + a1masksat
  

                if a1mask.sum() == 0:
                    a1retTmp = a1ret
                    a1obsTmp = a1obs

                else:
                    a1retTmp = ma.masked_where(a1mask, a1ret).compressed()
                    a1obsTmp = ma.masked_where(a1mask, a1obs).compressed()

                if a1retTmp.shape[0]==0: continue

                #************************************
                # Metrics
                #------------------------------------
                key = satid,surftype
                coef = d1coef[satid][idx_db]

                h[key] = h[key] + ((a1retTmp >= thpr)&(a1obsTmp >= thpr)).sum() *coef
                m[key] = m[key] + ((a1retTmp <  thpr)&(a1obsTmp >= thpr)).sum() *coef
                f[key] = f[key] + ((a1retTmp >= thpr)&(a1obsTmp < thpr)).sum()  *coef
                c[key] = c[key] + ((a1retTmp < thpr )&(a1obsTmp < thpr)).sum()  *coef

                a1maskprec1 = ma.masked_less(a1retTmp, thpr).mask
                a1maskprec2 = ma.masked_less(a1obsTmp, thpr).mask
                a1maskprec  = a1maskprec1 + a1maskprec2

                a1retCond   = ma.masked_where(a1maskprec, a1retTmp).compressed()
                a1obsCond   = ma.masked_where(a1maskprec, a1obsTmp).compressed()

                un [key]= un  [key]+ a1retTmp.shape[0] * coef
                uy [key]= uy  [key]+ a1retTmp.sum() * coef
                ux [key]= ux  [key]+ a1obsTmp.sum() * coef
                uyy[key]= uyy [key]+ (a1retTmp**2).sum() * coef
                uxx[key]= uxx [key]+ (a1obsTmp**2).sum() * coef
                uxy[key]= uxy [key]+ (a1retTmp * a1obsTmp).sum() * coef

                if a1retCond.shape[0] ==0: continue

                cn [key]= cn  [key]+ a1retCond.shape[0] * coef
                cy [key]= cy [key] + a1retCond.sum() * coef
                cx [key]= cx [key] + a1obsCond.sum() * coef
                cyy[key]= cyy[key] + (a1retCond**2).sum() * coef
                cxx[key]= cxx[key] + (a1obsCond**2).sum() * coef
                cxy[key]= cxy[key] + (a1retCond * a1obsCond).sum() * coef


    #*******************************
    lout = []
    for satid in [999] + lsatid:
        for surftype in lsurftype:
            key = satid,surftype 
            lout.append( [satid, surftype, 'h', h[key]]) 
            lout.append( [satid, surftype, 'm', m[key]]) 
            lout.append( [satid, surftype, 'f', f[key]]) 
            lout.append( [satid, surftype, 'c', c[key]]) 
            lout.append( [satid, surftype, 'un',  un[key]]) 
            lout.append( [satid, surftype, 'uy',  uy[key]]) 
            lout.append( [satid, surftype, 'ux',  ux[key]]) 
            lout.append( [satid, surftype, 'uyy', uyy[key]]) 
            lout.append( [satid, surftype, 'uxx', uxx[key]]) 
            lout.append( [satid, surftype, 'uxy', uxy[key]]) 
            lout.append( [satid, surftype, 'cn',  cn[key]]) 
            lout.append( [satid, surftype, 'cy',  cy[key]]) 
            lout.append( [satid, surftype, 'cx',  cx[key]])
            lout.append( [satid, surftype, 'cyy', cyy[key]]) 
            lout.append( [satid, surftype, 'cxx', cxx[key]]) 
            lout.append( [satid, surftype, 'cxy', cxy[key]]) 

    sout = util.list2csv(lout)

    #*******************************
    # Save
    #-------------------------------
    idx_db0 = min(lidx_db)
    idx_db1 = max(lidx_db)

    outDir  = '/home/utsumi/temp/ret'
    csvPath = outDir + '/multi.metrics.%s.%s.%05d-%05d.%s.csv'%(expr, coefflag, idx_db0, idx_db1, rettype)

    f=open(csvPath,'w'); f.write(sout); f.close()
    print csvPath

#*******************************
# Load  and calc metrics
#-------------------------------
for rettype in lrettype:
    idx_db0 = min(lidx_db)
    idx_db1 = max(lidx_db)

    outDir  = '/home/utsumi/temp/ret'
    csvPath = outDir + '/multi.metrics.%s.%s.%05d-%05d.%s.csv'%(expr, coefflag, idx_db0, idx_db1, rettype)

    f=open(csvPath,'r'); lines=f.readlines(); f.close()
    d = {}
    for line in lines:
        satid,surftype,varname,dat = line.strip().split(',')        
        satid = int(satid)
        dat   = float(dat)
 
        d[satid,surftype,varname] = dat

    lout = []

    for satid in [999] + lsatid:
        for surftype in lsurftype:
            h  = d[satid,surftype,'h'] 
            m  = d[satid,surftype,'m'] 
            f  = d[satid,surftype,'f'] 
            c  = d[satid,surftype,'c'] 

            un = d[satid,surftype,'un'] 
            uy = d[satid,surftype,'uy'] 
            ux = d[satid,surftype,'ux'] 
            uyy= d[satid,surftype,'uyy'] 
            uxx= d[satid,surftype,'uxx'] 
            uxy= d[satid,surftype,'uxy'] 

            cn = d[satid,surftype,'cn'] 
            cy = d[satid,surftype,'cy'] 
            cx = d[satid,surftype,'cx'] 
            cyy= d[satid,surftype,'cyy'] 
            cxx= d[satid,surftype,'cxx'] 
            cxy= d[satid,surftype,'cxy'] 

            urmse = np.sqrt( (uyy - 2*uxy + uxx)/un )
            crmse = np.sqrt( (cyy - 2*cxy + cxx)/cn )

            uxmean= ux/un
            cxmean= cx/cn
            uymean= uy/un
            cymean= cy/cn

            unrmse= urmse / uxmean
            cnrmse= crmse / cxmean

            ubias = (uy -ux)/un / uxmean
            cbias = (cy -cx)/cn / cxmean

            ucov  = uxy - uymean*ux - uxmean*uy + un*uxmean*uymean
            ccov  = cxy - cymean*cx - cxmean*cy + cn*cxmean*cymean

            uxvar = uxx - 2*uxmean*ux + un * uxmean**2
            cxvar = cxx - 2*cxmean*cx + cn * cxmean**2

            uyvar = uyy - 2*uymean*uy + un * uymean**2
            cyvar = cyy - 2*cymean*cy + cn * cymean**2

            ucc   = ucov / (uxvar*uyvar)**0.5
            ccc   = ccov / (cxvar*cyvar)**0.5

            pod   = float(h) / (h + m)
            far   = float(f) / (h + f)

            satName = dsatname[satid]
            ltmp = [sensor,satName,surftype,un,cn,unrmse,cnrmse,ubias,cbias,ucc,ccc,pod,far]
            lout.append(ltmp)

    lout = [['sensor','satellite','surface','un','cn','unrmse','cnrmse','unbias','cnbias','ucc','ccc','pod','far']] + lout

    sout = util.list2csv(lout)

    outPath = outDir + '/multi.metrics.summary.%s.%s.%05d-%05d.%s.csv'%(expr, coefflag, idx_db0, idx_db1, rettype)
    f=open(outPath, 'w'); f.write(sout); f.close()
    print outPath
  
