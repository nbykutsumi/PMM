from numpy import *
import JPLDB

#vname = 'z_ku'
vname = 'elev'
db_idx = 2601
#dbPath = '/work/hk01/utsumi/JPLDB/EPC_DB/GMI_EPC_DATABASE/db_02601.bin'
dbPath = '/work/hk01/utsumi/JPLDB/EPC_DB/GMI_EPC_DATABASE_TEST29/db_%05d.bin'%(db_idx)

db = JPLDB.JPLDB()
db.set_file(dbPath)
dat = db.get_var(vname)
print dat.shape
print dat
