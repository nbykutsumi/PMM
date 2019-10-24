from numpy import *
import numpy as np
import myfunc.util as util
import os, sys, glob, socket


lidx_db = np.arange(29*29*29)[1:].astype(int32)
#lidx_db = np.arange(29*29*29)[1:][:2000].astype(int32) # test
#lidx_db = [2477]
nsample = 1000

dbtype  = 'my'

DB_MAXREC    = 10000
DB_MINREC    = 1000
sensor  = 'GMI'
expr = 'org.smp%d'%(nsample)
#** Constants ******
myhost = socket.gethostname()
if myhost =='shui':
    if dbtype=='JPL':
        dbDir   = '/work/hk01/utsumi/JPLDB/EPC_DB/GMI_EPC_DATABASE_TEST29'
        countDir= '/work/hk01/utsumi/PMM/EPCDB/list'

    elif dbtype=='my':
        dbDir   = '/work/hk01/utsumi/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
        countDir= '/work/hk01/utsumi/PMM/EPCDB/list'
        retbaseDir = '/tank/utsumi/PMM/retsynt/%s'%(expr)
        outDir  = '/home/utsumi/temp/ret'

elif myhost == 'well':
    if dbtype=='JPL':
        dbDir   = '/work/hk01/utsumi/JPLDB/EPC_DB/GMI_EPC_DATABASE_TEST29'
        countDir= '/work/hk01/utsumi/PMM/EPCDB/list'
    elif dbtype=='my':
        dbDir   = '/media/disk2/share/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
        countDir= '/media/disk2/share/PMM/EPCDB/list'
        retbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retsynt/%s'%(expr)
        outDir  = '/home/utsumi/temp/ret'

#** Read histogram of EPC DB records **
countPath = countDir +'/count.epc.csv'
f=open(countPath,'r'); lines=f.readlines(); f.close()
ncol = len(lines[0].strip().split(',')) - 1
a2hist = np.zeros([29*29*29, ncol], int32)
for line in lines[1:]:
    line =  map(int,line.strip().split(','))
    epcid= line[0]
    a2hist[epcid] = line[1:]

#**************************************
sout = ''
for idx_db in lidx_db:
    srcDir = retbaseDir + '/%05d'%(idx_db)
    if not os.path.exists(srcDir): continue
    airec   = np.load(srcDir + '/irec.obs.%05d.npy'%(idx_db)).astype('int32')

    nsest    = np.load(srcDir + '/nsurfNS.est.%05d.npy'%(idx_db))
    msest    = np.load(srcDir + '/nsurfMS.est.%05d.npy'%(idx_db))
    nscmbest = np.load(srcDir + '/nsurfNScmb.est.%05d.npy'%(idx_db))
    mscmbest = np.load(srcDir + '/nsurfMScmb.est.%05d.npy'%(idx_db))

    nsobs    = np.load(srcDir + '/nsurfNS.obs.%05d.npy'%(idx_db))
    msobs    = np.load(srcDir + '/nsurfMS.obs.%05d.npy'%(idx_db))
    nscmbobs = np.load(srcDir + '/nsurfNScmb.obs.%05d.npy'%(idx_db))
    mscmbobs = np.load(srcDir + '/nsurfMScmb.obs.%05d.npy'%(idx_db))
    
    #-----------------------------------

    nsest    = ma.masked_invalid( ma.masked_less( nsest, 0))
    msest    = ma.masked_invalid( ma.masked_less( msest, 0))
    nscmbest = ma.masked_invalid( ma.masked_less( nscmbest, 0))
    mscmbest = ma.masked_invalid( ma.masked_less( mscmbest, 0))

    nsobs    = ma.masked_invalid( ma.masked_less( nsobs, 0))
    msobs    = ma.masked_invalid( ma.masked_less( msobs, 0))
    nscmbobs = ma.masked_invalid( ma.masked_less( nscmbobs, 0))
    mscmbobs = ma.masked_invalid( ma.masked_less( mscmbobs, 0))

    ccns     = np.ma.corrcoef(nsest,nsobs ,allow_masked=True)[0,1]
    ccms     = np.ma.corrcoef(msest,msobs ,allow_masked=True)[0,1]
    ccnscmb  = np.ma.corrcoef(nscmbest,nscmbobs ,allow_masked=True)[0,1]
    ccmscmb  = np.ma.corrcoef(mscmbest,mscmbobs ,allow_masked=True)[0,1]

    n0   = ma.masked_less(nscmbobs,0).count()
    n01  = ma.masked_less(nscmbobs,0.1).count()
    n1   = ma.masked_less(nscmbobs,1).count()
    n5   = ma.masked_less(nscmbobs,5).count()
    n10  = ma.masked_less(nscmbobs,10).count()
    n30  = ma.masked_less(nscmbobs,30).count()

    bratns    = ma.masked_invalid(nsest.mean() / nsobs.mean() )
    bratms    = ma.masked_invalid(msest.mean() / msobs.mean() )
    bratnscmb = ma.masked_invalid(nscmbest.mean() / nscmbobs.mean() )
    bratmscmb = ma.masked_invalid(mscmbest.mean() / mscmbobs.mean() )

    meanobs   = nscmbobs.mean()


    #-- variables in the sampled DB ----
    t2m      = np.load(dbDir + '/t2m/%05d.npy'%(idx_db))[airec]
    elev     = np.load(dbDir + '/gtopo/%05d.npy'%(idx_db))[airec]
    surftype = np.load(dbDir + '/surfaceTypeIndex/%05d.npy'%(idx_db))[airec]


    a1maskinvalid = ma.masked_invalid(nscmbobs).mask
    t2m      = ma.masked_where(a1maskinvalid, t2m).compressed()
    elev     = ma.masked_where(a1maskinvalid, elev).compressed()
    surftype = ma.masked_where(a1maskinvalid, surftype).compressed()
 
    #-- Count environmental variables --
    nt0 = ma.masked_greater_equal(t2m, 273.15+0).count()
    nt5 = ma.masked_outside(t2m,273.15+0, 273.15+5).count()
    nt10= ma.masked_outside(t2m,273.15+5, 273.15+10).count()
    nt15= ma.masked_outside(t2m,273.15+10,273.15+15).count()
    nt20= ma.masked_outside(t2m,273.15+15,273.15+20).count()
    nt25= ma.masked_outside(t2m,273.15+20,273.15+25).count()
    nt30= ma.masked_outside(t2m,273.15+25,273.15+30).count()
    nt99= ma.masked_outside(t2m,273.15+30,273.15+99).count()

    nelev100= ma.masked_greater(elev, 100).count()
    nelev500 = ma.masked_outside(elev,100,500).count()
    nelev1000= ma.masked_outside(elev,500,1000).count()
    nelev2000= ma.masked_outside(elev,1000,2000).count()
    nelev9000= ma.masked_less(elev,2000).count()

    nsea     = ma.masked_not_equal(surftype,1).count()
    nice     = ma.masked_not_equal(surftype,2).count()  
    nvegh    = ma.masked_outside(surftype,3,5).count() 
    nvegl    = ma.masked_outside(surftype,6,7).count() 
    nsnow    = ma.masked_outside(surftype,8,11).count()
    nswater  = ma.masked_not_equal(surftype,12).count()
    ncoast   = ma.masked_not_equal(surftype,13).count()
    nseaiceedge = ma.masked_not_equal(surftype,14).count()



    #-- Linear regression ---
    if ccnscmb is np.ma.masked:
        a,b = np.ma.masked, np.ma.masked
    else:
        a,b = np.polyfit( nscmbobs, nscmbest, 1)    



    #ltmp = [idx_db, n0, n01, n1, n5, n10, n30, meanobs, ccns, ccms, ccnscmb, ccmscmb, bratns, bratms, bratnscmb, bratmscmb, a, b]
    ltmp = [idx_db, n0, n01, n1, n5, n10, n30, nt0, nt5, nt10, nt15, nt20, nt25, nt30, nt99, nelev100, nelev500, nelev1000, nelev2000, nelev9000, nsea, nice, nvegh, nvegl, nsnow, nswater, ncoast, nseaiceedge, meanobs, ccns, ccms, ccnscmb, ccmscmb, bratns, bratms, bratnscmb, bratmscmb, a, b]
    sout = sout + ','.join( map(str, ltmp) ) + '\n'
    print idx_db,a,b

##**********************-
## Save
##-----------------------
slabel = ','.join(['idx_db', 'n0', 'n0.1', 'n1', 'n5', 'n10', 'n30',  '<=0degC','0-5degC','5-10degC','10-15degC','15-20degC','20-25degC','25-30degC','>30degC', '<=100', '100-500m', '500-1000m', '1000-2000m', '2000-9000m', 'sea', 'ice', 'veg-high', 'veg-low', 'snow', 'standing-water','Coast', 'Sea-ice-edge', 'meanobs', 'cc_ns', 'cc_ms', 'cc_nscmb', 'cc_mscmb', 'brat_ns', 'brat_ms', 'brat_nscmb', 'brat_mscmb','a','b'])

sout = slabel + '\n' + sout
csvPath = outDir + '/synt.%s.csv'%(expr)
f=open(csvPath, 'w'); f.write(sout); f.close()
print csvPath

