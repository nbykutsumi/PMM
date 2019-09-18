#import matplotlib
#matplotlib.use('Agg')
#impor://www.atmarkit.co.jp/ait/articles/1906/27/news016.htmlt matplotlib.pyplot as plt
from numpy import *
from scipy.stats import gaussian_kde
import socket
import random
import os, sys
import numpy as np
import epcfunc
import EPCDB

sensor = 'GMI'
myhost = socket.gethostname()

if myhost =="shui":
    dbDir   = '/tank/utsumi/PMM/EPCDB/samp.20000.GMI.V05A.S1.ABp103-117.01-12'
    outDir  = '/home/utsumi/temp/ret'
elif myhost =="well":
    coefDir = '/media/disk2/share/PMM/EPCDB/EPC_COEF/%s'%(sensor)
    dbDir   = '/media/disk2/share/PMM/EPCDB/samp.20000.GMI.V05A.S1.ABp103-117.01-12'
    outDir  = '/home/utsumi/temp/ret'


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
    nrainDir= dbDir + '/nrain'
    srcPath = nrainDir + '/db_%05d.bin.nrain.txt'%(idx_db)
    a1nrain = read_table(srcPath,type=int32)[0]
    a1nrain_cold = a1nrain[:6]
    a1nrain_warm = a1nrain[6:]
    return a1nrain_warm, a1nrain_cold


#* Chose database *****
dbtype= 'my'
if dbtype =='JPL':
    db    = JPLDB.JPLDB()
elif dbtype=='my':
    db    = EPCDB.EPCDB()
else:
    print 'check dbtype',dbtype
#**********************
DB_MAXREC = 20000
NEM = 12


#** Make idx_db list ***
ndb = 29*29*29  # 24388
lidx_db = []
for idx_db in range(ndb):
    nrainDir= dbDir + '/nrain'
    srcPath = nrainDir + '/db_%05d.bin.nrain.txt'%(idx_db)
    if not os.path.exists(srcPath):
        continue
    _,a1nrain = read_nrain(idx_db)
    nall = a1nrain[0]
    if nall >10000:
        lidx_db.append(idx_db)
#***********************
lidx_db = np.sort(lidx_db)
#lidx_db = np.sort(random.sample(lidx_db, 10))


#lidx_db = [3312,3948,4750,8859,9503,9865,13668,13752,14425,19618,20666,22617]
#lidx_db = [22617]
#lidx_db = [20666]
#lidx_db = [3312]
#lidx_db = [4750]
#lidx_db = [9865,4750]
icount = 0
for idx_db in lidx_db:
    print 'idx_db=',idx_db

    dbPath = dbDir + '/db_%05d.bin'%(idx_db)
    if   dbtype == 'JPL':
        db.set_file(dbPath)
    elif dbtype == 'my':
        db.set_idx_db(dbDir, idx_db)

    a2epcdb = db.get_var('pc_emis', nrec=DB_MAXREC)[:,:NEM]  # (nrec, 12)
    a1nsurfNScmb = db.get_var('precip_NS_cmb', nrec=DB_MAXREC)
    a1t2m        = db.get_var('t2m', nrec=DB_MAXREC)

    #-- Screening ----------
    a1flagepc= ma.masked_not_equal(a2epcdb,-9999.).all(axis=1).mask
    a1flagT  = ma.masked_greater(a1t2m,0).mask   # screen T2m=0
    a1flag_good = a1flagepc * a1flagT

    a1flag_rain = ma.masked_greater(a1nsurfNScmb, 0).mask
    a1flag_no   = ma.masked_equal(a1nsurfNScmb,0).mask

    a1flag_rain = a1flag_rain * a1flag_good
    a1flag_no   = a1flag_no   * a1flag_good

    #-------------
    if a1flag_rain.sum() + a1flag_no.sum() ==0:
        print 'No valid data in idx_db',idx_db
        continue

    icount = icount + 1
    #-------------
    p_rain     = float(a1flag_rain.sum()) / ( a1flag_rain.sum() + a1flag_no.sum() )
    p_no       = float(a1flag_no.sum()) / ( a1flag_rain.sum() + a1flag_no.sum() )

    avmin = []
    avmax = []
    a2bnd = []  # just for info
    a2p_x      = []
    a2p_x_rain = []
    a2p_x_no   = []
    #for iepc in [6]:  # test
    for iepc in range(NEM+1):
        if iepc==NEM:
            a1epc_all  = a1t2m[a1flag_good]
            a1epc_rain = a1t2m[a1flag_rain]
            a1epc_no   = a1t2m[a1flag_no  ]
        else:
            a1epc_all  = a2epcdb[a1flag_good,iepc]
            a1epc_rain = a2epcdb[a1flag_rain,iepc]
            a1epc_no   = a2epcdb[a1flag_no  ,iepc]



        #vmin, vmax = a1epc_all.min(), a1epc_all.max()
        vmin,vmax = np.percentile(a1epc_all,0.05), np.percentile(a1epc_all, 99.95)
        nbin = 50
        abnd = np.linspace(vmin,vmax,nbin+1)
        x    = 0.5*(abnd[:-1]+abnd[1:])
        w    = abnd[1]-abnd[0]
        #-- KDE ----------------
        kde_all  = gaussian_kde(a1epc_all)
        kde_rain = gaussian_kde(a1epc_rain)
        kde_no   = gaussian_kde(a1epc_no)

        #-- Avoid zero-probability ----
        ap_x      = kde_all(x)
        ap_x_rain = kde_rain(x)
        ap_x_no   = kde_no(x)

        ap_x      = ap_x + 1e-20
        ap_x_rain = ap_x_rain + 1e-20
        ap_x_no   = ap_x_no   + 1e-20

        ap_x      = ap_x / ap_x.sum()
        ap_x_rain = ap_x_rain / ap_x_rain.sum()
        ap_x_no   = ap_x_no   / ap_x_no.sum()
        #------------------------------

        avmin.append(vmin)
        avmax.append(vmax)       
        a2bnd.append(abnd)
    
        a2p_x     .append(ap_x)
        a2p_x_rain.append(ap_x_rain)
        a2p_x_no  .append(ap_x_no)

        ##-- test ---------        
        ix  = 0
        p_rain_x = ap_x_rain[ix] * p_rain / ap_x[ix]
        p_no_x   = ap_x_no[ix]   * p_no   / ap_x[ix]
        #print iepc,p_rain_x, p_no_x, p_rain_x + p_no_x

    avmin      = np.array(avmin)  # NEM+1
    avmax      = np.array(avmax)  # NEM+1
    a2bnd      = np.array(a2bnd)
    a2p_x      = np.array(a2p_x)       # (NEM+1, NBIN)
    a2p_x_rain = np.array(a2p_x_rain)  # (NEM+1, NBIN)
    a2p_x_no   = np.array(a2p_x_no)    # (NEM+1, NBIN)

    #print a2p_x_rain.shape

    ##-- test --
    a2x = np.concatenate([a2epcdb, array(a1t2m).reshape(-1,1)],axis=1)  # (NREC,NEM+1)
    a1w = (avmax-avmin)/nbin               # NEM+1
    a2k = ((a2x-avmin)/a1w).astype(int32)  # (NREC, NEM+1)
    a2k = ma.masked_less(a2k,0).filled(0)
    a2k = ma.masked_greater(a2k,nbin-1).filled(nbin-1)

    a2p_xTmp      = np.array([a2p_x[iepc, a2k[:,iepc]] for iepc in range(NEM+1)]).T  # (NREC, NEM+1)
    a2p_x_rainTmp = np.array([a2p_x_rain[iepc, a2k[:,iepc]] for iepc in range(NEM+1)]).T  # (NREC, NEM+1)
    a2p_x_noTmp   = np.array([a2p_x_no[iepc, a2k[:,iepc]] for iepc in range(NEM+1)]).T  # (NREC, NEM+1)



    logp_rain_x = np.log(a2p_x_rainTmp).sum(axis=1) + np.log(p_rain)
    logp_no_x   = np.log(a2p_x_noTmp  ).sum(axis=1) + np.log(p_no  )

    thdif = log(99.99)
    #thdif = log(1)

    lthpr =[0,0.1,1,5,10]
    dpod = {}
    dfar = {}
    dmr  = {}
    drma = {}
    sline = '%d'%(idx_db)
    for thpr in lthpr: 
        a1flag_prec0 = ma.masked_greater_equal(a1nsurfNScmb,thpr).mask
        a1flag_prec1 = ~(ma.masked_invalid(a1nsurfNScmb).mask)
        a1flag_prec  = a1flag_prec0 * a1flag_prec1

        a1flag_est_no = ma.masked_greater(logp_no_x - logp_rain_x, thdif).mask
        a1flag_est_rain = ~a1flag_est_no

        a1flag_est_no = a1flag_est_no * a1flag_prec
        a1flag_est_rain = a1flag_est_rain * a1flag_prec
        
        OE = (a1flag_rain * a1flag_est_rain).sum()
        Oe = (a1flag_rain * a1flag_est_no).sum()
        oE = (a1flag_no   * a1flag_est_rain).sum()
        oe = (a1flag_no   * a1flag_est_no).sum()
    
        pod = OE/float(OE+Oe)
        far = oE/float(OE+oE)
        mr  = Oe/float(OE+Oe)
   
        rma = ma.masked_where( ~(a1flag_est_no *a1flag_prec), a1nsurfNScmb ).sum() / ma.masked_where(~a1flag_prec, a1nsurfNScmb).sum()
 
        nall= OE+Oe+oe+oE
        #print '----------------------'
        #print 'thrp=',thpr
        #print 'nall=',nall
        #print 'OE  =',OE
        #print 'Oe  =',Oe
        #print 'oe  =',oe
        #print 'oE  =',oE
        #print ''
        #print 'POD =',pod
        #print 'FAR =',far 
        #print 'MR  =',mr
        #print 'RMA =',rma
        #print 'random, OE/(OE+oE)', p_rain, OE/float(OE+oE)
        #print ''
        dpod[thpr] = pod
        dfar[thpr] = far
        dmr [thpr] = mr
        drma[thpr] = rma

        lout  = [nall,OE,Oe,oe,oE]
        sline = sline + ',' + ','.join(map(str,lout))
    
    lout = [d[thpr] for d in [dpod,dfar,dmr,drma] for thpr in lthpr]
    sline = sline + ',' + ','.join(map(str,lout))


    if (icount==1):
        llabel = ['db_idx'] + ['%s %.1f'%(var,thpr) for thpr in lthpr for var in ['nall','OE','Oe','oe','oE']]
        llabel = llabel + ['%s %.1f'%(var,thpr) for var in ['POD','FAR','MR','RMA'] for thpr in lthpr]
        slabel = ','.join(map(str,llabel))
        sout   = slabel + '\n'
        sout   = sout + sline + '\n'
    else:
        sout = sout + sline + '\n'


#-- Save -------------
outPath = outDir + '/RNR.bayes.csv'
f=open(outPath,'w'); f.write(sout); f.close()
print outPath



    ##-- test --
    #a2x = np.concatenate([a2epcdb, array(a1t2m).reshape(-1,1)],axis=1)[irec].reshape(1,-1)  # (NREC,NEM+1)
    #iepc0 = 6
    #iepc1 = 8
    #avmin = avmin[iepc0:iepc1]
    #avmax = avmax[iepc0:iepc1]
    #a2p_x = a2p_x[iepc0:iepc1,:]
    #a2p_x_rain = a2p_x_rain[iepc0:iepc1,:]
    #a2p_x_no   = a2p_x_no[iepc0:iepc1,:]
    #NEM = iepc1-iepc0-1
    #a2x = a2x[:,iepc0:iepc1]


    #a1w = (avmax-avmin)/nbin               # NEM+1
    #a2k = ((a2x-avmin)/a1w).astype(int32)  # (NREC, NEM+1)
    #a2k = ma.masked_less(a2k,0).filled(0)
    #a2k = ma.masked_greater(a2k,nbin-1).filled(nbin-1)

    #a2p_xTmp      = np.array([a2p_x[iepc, a2k[:,iepc]] for iepc in range(NEM+1)]).T  # (NREC, NEM+1)
    #a2p_x_rainTmp = np.array([a2p_x_rain[iepc, a2k[:,iepc]] for iepc in range(NEM+1)]).T  # (NREC, NEM+1)
    #a2p_x_noTmp   = np.array([a2p_x_no[iepc, a2k[:,iepc]] for iepc in range(NEM+1)]).T  # (NREC, NEM+1)

    #logp_rain_x = np.log(a2p_x_rainTmp).sum(axis=1) + np.log(p_rain) - np.log(a2p_xTmp).sum(axis=1)
    #logp_no_x   = np.log(a2p_x_noTmp  ).sum(axis=1) + np.log(p_no  ) - np.log(a2p_xTmp).sum(axis=1)



    #p_rain_x = np.exp(logp_rain_x)
    #p_no_x   = np.exp(logp_no_x)
    #print p_rain_x, p_no_x, p_rain_x +p_no_x



    ##-- test --
    #a2x = np.concatenate([a2epcdb, array(a1t2m).reshape(-1,1)],axis=1)[irec].reshape(1,-1)  # (NREC,NEM+1)
    #avmin = avmin[6:7]
    #avmax = avmax[6:7]
    #a2p_x = a2p_x[6:7,:]
    #a2p_x_rain = a2p_x_rain[6:7,:]
    #a2p_x_no   = a2p_x_no[6:7,:]
    #NEM = 0
    #a2x = a2x[:,6:7]


    #a1w = (avmax-avmin)/nbin               # NEM+1
    #a2k = ((a2x-avmin)/a1w).astype(int32)  # (NREC, NEM+1)
    #a2k = ma.masked_less(a2k,0).filled(0)
    #a2k = ma.masked_greater(a2k,nbin-1).filled(nbin-1)

    #a2p_xTmp      = np.array([a2p_x[iepc, a2k[:,iepc]] for iepc in range(NEM+1)]).T  # (NREC, NEM+1)
    #a2p_x_rainTmp = np.array([a2p_x_rain[iepc, a2k[:,iepc]] for iepc in range(NEM+1)]).T  # (NREC, NEM+1)
    #a2p_x_noTmp   = np.array([a2p_x_no[iepc, a2k[:,iepc]] for iepc in range(NEM+1)]).T  # (NREC, NEM+1)

    #logp_rain_x = np.log(a2p_x_rainTmp) + np.log(p_rain) - np.log(a2p_xTmp)
    #logp_no_x   = np.log(a2p_x_noTmp  ) + np.log(p_no  ) - np.log(a2p_xTmp)

    #p_rain_x = np.exp(logp_rain_x)
    #p_no_x   = np.exp(logp_no_x)
    #print p_rain_x, p_no_x, p_rain_x +p_no_x


    ##-- test --
    #irec = 0    # rain
    ##irec = 3    # no-rain
    ##irec = 100  # rain
    ##irec = 200  # T2m=0
    ##irec = 300  # rain
    ##print a2epcdb[irec], a1t2m[irec]
    #a1x = np.concatenate([a2epcdb[irec], array([a1t2m[irec]])])  # NREC+1
    #a1w = (avmax-avmin)/nbin   # NREC+1
    #a1k = ((a1x-avmin)/a1w).astype(int32)  # NREC+1
    #a1iepc = np.arange(NEM+1)
    #ap_x = a2p_x[a1iepc,a1k]  # NREC+1
    #ap_x_rain = a2p_x_rain[a1iepc,a1k]  # NREC+1
    #ap_x_no   = a2p_x_no  [a1iepc,a1k]  # NREC+1

    ##ap_x      = ap_x + 1e-30

    #logp_rain_x = np.log(ap_x_rain).sum() + np.log(p_rain) - np.log(ap_x).sum()
    #logp_no_x   = np.log(ap_x_no  ).sum() + np.log(p_no  ) - np.log(ap_x).sum()

    #print logp_rain_x, np.log(ap_x_rain).sum() ,np.log(p_rain) , np.log(ap_x).sum()

    #p_rain_x = np.exp(logp_rain_x)
    #p_no_x   = np.exp(logp_no_x)
    #print p_rain_x, p_no_x, p_rain_x +p_no_x
