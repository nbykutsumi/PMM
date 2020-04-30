# %%
%matplotlib inline
import numpy as np
from numpy import *
import glob, h5py
import myfunc.util as util
from datetime import datetime, timedelta
import EPCDB
import matplotlib.pyplot as plt

DB_MAXREC = 10000
dbDir   = '/media/disk2/share/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
db = EPCDB.EPCDB()

idx_db = 5810
db.set_idx_db(dbDir, idx_db)

srcdir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retsynt/prf.smp1000/%05d'%(idx_db)
a2profest = np.load(srcdir + '/precip_water_prof_NS_rs.est.%05d.npy'%(idx_db))
a1irec = np.load(srcdir + '/irec.obs.%05d.npy'%(idx_db))

a2profobs = db.get_var('precip_water_prof_NS_relsurf')[a1irec,-50:]
a1precobs = db.get_var('precip_NS_cmb')[a1irec]
a1elev    = db.get_var('elev')[a1irec]

a2profest = ma.masked_less(a2profest,0)
a2profobs = ma.masked_less(a2profobs,0)


#for i in range(len(a1elev)):
#    elev = a1elev[i]
#    prec = a1precobs[i]
#    if (elev > 500)&(prec>1):
#        break

a1profest = a2profest[i]
a1profobs = a2profobs[i]

fig = plt.figure()
plt.plot(a1profobs, color='k')
plt.plot(a1profest, color='r')
plt.show()

# %%
