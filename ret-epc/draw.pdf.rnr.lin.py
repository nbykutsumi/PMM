import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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


###** Make idx_db list ***
#ndb = 29*29*29  # 24388
#lidx_db = []
#for idx_db in range(ndb):
#    nrainDir= dbDir + '/nrain'
#    srcPath = nrainDir + '/db_%05d.bin.nrain.txt'%(idx_db)
#    if not os.path.exists(srcPath):
#        continue
#    _,a1nrain = read_nrain(idx_db)
#    nall = a1nrain[0]
#    if nall >10000:
#        lidx_db.append(idx_db)
#lidx_db = np.sort(lidx_db)
##***********************
#lidx_db = np.sort(random.sample(lidx_db, 10))


#lidx_db = [3312,3948,4750,8859,9503,9865,13668,13752,14425,19618,20666,22617]
#lidx_db = [22617]
#lidx_db = [20666]
#lidx_db = [3312]
#lidx_db = [4750]
#lidx_db = [9865,4750]
#lidx_db = [11719]  # Precip ratio>0.9
lidx_db = [4120]
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
    dd_all = {}
    dd_rain= {}
    dd_no  = {}
    a1d_all = None
    a1d_rain= None
    a1d_no  = None
    for iepc in range(NEM+1):
        if iepc==NEM:
            a1epc_all  = a1t2m[a1flag_good]
            a1epc_rain = a1t2m[a1flag_rain]
            a1epc_no   = a1t2m[a1flag_no  ]
        else:
            a1epc_all  = a2epcdb[a1flag_good,iepc]
            a1epc_rain = a2epcdb[a1flag_rain,iepc]
            a1epc_no   = a2epcdb[a1flag_no  ,iepc]

        #-- Normalize ----
        m_org   = a1epc_all.mean()
        std_org = a1epc_all.std()
        a1epc_all = (a1epc_all - m_org)/std_org
        a1epc_rain= (a1epc_rain- m_org)/std_org
        a1epc_no  = (a1epc_no  - m_org)/std_org

        m_rain  = a1epc_rain.mean()
        m_no    = a1epc_no.mean()
        std_rain= a1epc_rain.std()
        std_no  = a1epc_no.std()

        a1d_allTmp  = (m_no - m_rain)/(std_no + std_rain)*(m_no - a1epc_all)
        a1d_rainTmp = (m_no - m_rain)/(std_no + std_rain)*(m_no - a1epc_rain)
        a1d_noTmp   = (m_no - m_rain)/(std_no + std_rain)*(m_no - a1epc_no)

        dd_all[iepc] = a1d_allTmp
        dd_rain[iepc]= a1d_rainTmp
        dd_no[iepc]  = a1d_noTmp

        if a1d_all is None:
            a1d_all  = a1d_allTmp
            a1d_rain = a1d_rainTmp
            a1d_no   = a1d_noTmp
        else: 
            a1d_all  = a1d_all  + a1d_allTmp
            a1d_rain = a1d_rain + a1d_rainTmp
            a1d_no   = a1d_no   + a1d_noTmp

    dd_all[NEM+1]  = a1d_all
    dd_rain[NEM+1] = a1d_rain
    dd_no[NEM+1]   = a1d_no
    #-------------
    #-- Figure of PDF ----------
    fig = plt.figure(figsize=(17,10))
    h = 0.23
    w = 0.16
    for iepc in range(NEM+2):

        #-- Histogram ----------
        a1d_all = dd_all[iepc]
        a1d_rain= dd_rain[iepc]
        a1d_no  = dd_no[iepc]

        #vmin, vmax = a1epc_all.min(), a1epc_all.max()
        vmin,vmax = np.percentile(a1d_all,0.05), np.percentile(a1d_all, 99.95)
        abnd = np.linspace(vmin,vmax,50)
        afreq_all,  abnds_all  = np.histogram(a1d_all,  bins=abnd)
        afreq_rain, abnds_rain = np.histogram(a1d_rain, bins=abnd)
        afreq_no,   abnds_no   = np.histogram(a1d_no,   bins=abnd)
        #afreq_rain, abnds_rain = np.histogram(a1epc_rain, bins=100)
        #afreq_no,   abnds_no   = np.histogram(a1epc_no,bins=100)

        y_all  = afreq_all.astype(float32) / afreq_all.sum()
        y_rain = afreq_rain.astype(float32) / afreq_rain.sum()
        y_no   = afreq_no.astype(float32) / afreq_no.sum()
        #y_rain = afreq_rain
        #y_no   = afreq_no

        x_all  = 0.5*(abnds_all[:-1] +abnds_all[1:])
        x_rain = 0.5*(abnds_rain[:-1]+abnds_rain[1:])
        x_no   = 0.5*(abnds_no[:-1]  +abnds_no[1:])

        ##-- KDE ----------------
        #kde_rain = gaussian_kde(a1d_rain)
        #kde_no   = gaussian_kde(a1d_no)
        #k_rain   = kde_rain(x_rain)
        #k_no     = kde_no(x_no)



        #-- Draw ---------------
        if int(iepc/5)   ==0:
            y0 = 0.05 + (h+0.05)*2
        elif int(iepc/5) ==1:
            y0 = 0.05 + (h+0.05)*1
        elif int(iepc/5) ==2:
            y0 = 0.05

        x0 = 0.05 +  (iepc%5)*(w+0.03)
        ax = fig.add_axes([x0,y0,w,h])

        ax.plot(x_all,  y_all,  '-',color='k',linewidth=2)
        ax.plot(x_rain, y_rain, '--',color='r')
        ax.plot(x_no,   y_no, '--',color='k')

        #ax2= ax.twinx()
        #ax2.plot(x_rain, k_rain, '-',color='r')
        #ax2.plot(x_no,   k_no, '-',color='k')
        if iepc==NEM:
            plt.title('D(T2m)')
        elif iepc==NEM+1:
            plt.title('D(All)')

        else:
            plt.title('D(PC%d)'%(iepc+1))
    plt.suptitle('DB=%05d Rain=%d NoRain=%d'%(idx_db, len(a1epc_rain), len(a1epc_no)))
    #-- Save figure ---
    figPath = outDir + '/pdf.NRN.lin.%05d.png'%(idx_db)
    plt.savefig(figPath)
    print figPath
                                




    '''
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
    '''
