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


##** Make idx_db list ***
#ndb = 29*29*29  # 24388
#lidx_db = []
#for idx_db in range(ndb):
#    nrainDir= dbDir + '/nrain'
#    srcPath = nrainDir + '/db_%05d.bin.nrain.txt'%(idx_db)
#    if not os.path.exists(srcPath):
#        continue
#    _,a1nrain = read_nrain(idx_db)
#    nall = a1nrain[0]
#    if nall >19000:
#        lidx_db.append(idx_db)
##***********************
#lidx_db = np.sort(random.sample(lidx_db, 10))

lidx_db = [3312,3948,4750,8859,9503,9865,13668,13752,14425,19618,20666,22617]
#lidx_db = [3948]

for idx_db in lidx_db:
    print idx_db

    dbPath = dbDir + '/db_%05d.bin'%(idx_db)
    if   dbtype == 'JPL':
        db.set_file(dbPath)
    elif dbtype == 'my':
        db.set_idx_db(dbDir, idx_db)

    #dbPath = dbDir + '/db_%05d.bin'%(idx_db_expand)
    #db.set_idx_db(dbDir, idx_db_expand)
    a2epcdb = db.get_var('pc_emis', nrec=DB_MAXREC)[:,:NEM]  # (nrec, 12)
    #a1nsurfMScmb = db.get_var('precip_MS_cmb', nrec=DB_MAXREC)
    a1nsurfNScmb = db.get_var('precip_NS_cmb', nrec=DB_MAXREC)
    #a1nsurfMS    = db.get_var('precip_nsfc_MS', nrec=DB_MAXREC)
    #a1nsurfNS    = db.get_var('precip_nsfc_NS', nrec=DB_MAXREC)
    a1t2m        = db.get_var('t2m', nrec=DB_MAXREC)

    a1flag_rain = ma.masked_greater(a1nsurfNScmb, 0).mask
    a1flag_no   = ma.masked_equal(a1nsurfNScmb,0).mask

    #-- Figure of PDF ----------
    fig = plt.figure(figsize=(17,10))
    h = 0.23
    w = 0.16
    for iepc in range(NEM+1):
        if iepc==NEM:
            a1flagepc  = ma.masked_greater(a1t2m,0).mask   # screen T2m=0
            a1epc_all  = a1t2m[a1flagepc]
            a1epc_rain = a1t2m[a1flag_rain*a1flagepc]
            a1epc_no   = a1t2m[a1flag_no  *a1flagepc]
        else:
            a1flagepc  = ma.masked_not_equal(a2epcdb,-9999.).all(axis=1).mask  # screen epc=-9999.
            a1epc_all  = a2epcdb[a1flagepc,iepc]
            a1epc_rain = a2epcdb[a1flag_rain*a1flagepc,iepc]
            a1epc_no   = a2epcdb[a1flag_no*a1flagepc  ,iepc]



        #-- Histogram ----------
        #vmin, vmax = a1epc_all.min(), a1epc_all.max()
        vmin,vmax = np.percentile(a1epc_all,0.05), np.percentile(a1epc_all, 99.95)
        abnd = np.linspace(vmin,vmax,50)
        afreq_rain, abnds_rain = np.histogram(a1epc_rain, bins=abnd)
        afreq_no,   abnds_no   = np.histogram(a1epc_no,bins=abnd)
        #afreq_rain, abnds_rain = np.histogram(a1epc_rain, bins=100)
        #afreq_no,   abnds_no   = np.histogram(a1epc_no,bins=100)

        y_rain = afreq_rain.astype(float32) / afreq_rain.sum()
        y_no   = afreq_no.astype(float32) / afreq_no.sum()
        #y_rain = afreq_rain
        #y_no   = afreq_no

        x_rain = 0.5*(abnds_rain[:-1]+abnds_rain[1:])        
        x_no   = 0.5*(abnds_no[:-1]  +abnds_no[1:]) 

        #-- KDE ----------------
        kde_rain = gaussian_kde(a1epc_rain)
        kde_no   = gaussian_kde(a1epc_no)
        k_rain   = kde_rain(x_rain)
        k_no     = kde_no(x_no)

        #-- Draw ---------------
        if int(iepc/5)   ==0:
            y0 = 0.05 + (h+0.05)*2
        elif int(iepc/5) ==1:
            y0 = 0.05 + (h+0.05)*1
        elif int(iepc/5) ==2:
            y0 = 0.05

        x0 = 0.05 +  (iepc%5)*(w+0.03)
        ax = fig.add_axes([x0,y0,w,h])

        ax.plot(x_rain, y_rain, '--',color='r')
        ax.plot(x_no,   y_no, '--',color='k')

        ax2= ax.twinx()
        ax2.plot(x_rain, k_rain, '-',color='r')
        ax2.plot(x_no,   k_no, '-',color='k')
        if iepc==NEM:
            plt.title('T2m')
        else:
            plt.title('PC%d'%(iepc+1))
    plt.suptitle('DB=%05d Rain=%d NoRain=%d'%(idx_db, len(a1epc_rain), len(a1epc_no)))
    #-- Save figure ---
    figPath = outDir + '/pdf.NRN.%05d.png'%(idx_db)
    plt.savefig(figPath)
    print figPath
