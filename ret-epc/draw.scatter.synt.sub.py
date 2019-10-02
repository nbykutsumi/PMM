import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import numpy as np
import myfunc.util as util
import os, sys, glob, socket

lidx_db = [1652]

nsample = 1000
DB_MAXREC = 10000
DB_MINREC = 1000

dbtype  = 'my'
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


for idx_db in lidx_db:
    retDir = retbaseDir + '/%05d'%(idx_db)
    a1est = np.load( retDir + '/nsurfNScmb.est.%05d.npy'%(idx_db))
    a1obs = np.load( retDir + '/nsurfNScmb.obs.%05d.npy'%(idx_db))


    vmax= 50
    fig = plt.figure(figsize=(4,4))
    ax  = fig.add_axes([0.15,0.15,0.8,0.8])
    ax.scatter(a1obs, a1est, s=1, color='k')
    ax.plot([0,vmax], [0,vmax], '-', color='k')
    ax.set_xlabel('Obs (NScmb)')
    ax.set_ylabel('Est (NScmb)') 
    ax.set_ylim([0,vmax])
    ax.set_xlim([0,vmax])
    
    stitle = 'DB:%05d %s'%(idx_db, expr)
    plt.title(stitle)

    figPath = outDir + '/scatter.%s.%05d.png'%(expr,idx_db)
    plt.savefig(figPath)
    print figPath
