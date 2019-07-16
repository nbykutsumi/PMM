from numpy import *
import EPCDB

db = EPCDB.EPCDB()


baseDir = '/work/hk01/utsumi/PMM/EPCDB/samp.5000.GMI.V05A.S1.ABp103-117'
vname = 'pc_emis'
#vname = 'elev'
idx_db = 1652
#idx_db = 1649
db.set_idx_db(baseDir, idx_db=idx_db)

dat = db.get_var(vname)
print dat
print dat.shape 
#dictvars = db.dictvars
#for key, val in dictvars.iteritems():
#    print key, val
