import JPLDB

jpldb = JPLDB.JPLDB()

idx_db = 624
srcPath ='/work/hk01/utsumi/JPLDB/EPC_DB/GMI_EPC_DATABASE/db_%05d.bin'%(idx_db)

varName = 'precip_MS_cmb'
jpldb.set_file(srcPath)
a1prcp  = jpldb.get_var(varName, nrec=20000)
a1jns   = jpldb.get_var('j_NS',  nrec=20000)

for i in range(len(a1jns)):
    jns = a1jns[i]
    prcp= a1prcp[i]
    if (jns<12)or(36<jns):
        print i,jns, prcp


