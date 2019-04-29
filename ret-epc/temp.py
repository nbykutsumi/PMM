from numpy import *
import JPLDB
import h5py

db = JPLDB.JPLDB()
#DB_MAXREC = 20000
#for idx_db in range(1,10000,50):
lvar = ['precip_prof_MS','precip_prof_NS','precip_prof_MS_cmb','precip_prof_NS_cmb']

for var in lvar:
    for idx_db in [2601]:
        dbDir   = '/work/hk01/utsumi/JPLDB/EPC_DB/GMI_EPC_DATABASE'
        srcPath = dbDir + '/db_%05d.bin'%(idx_db)
        db.set_file(srcPath)
        #a1nsurfMS  = db.get_var('precip_MS_cmb',  nrec=DB_MAXREC)
        #a2prprofMS = db.get_var('precip_prof_MS', nrec=DB_MAXREC)
    
        a1nsurfMS  = db.get_var('precip_MS_cmb')
        a2prprofMS = db.get_var(var)
        #a2prprofMS = db.get_var('precip_prof_MS')
        print ''
        print var
        print idx_db, a1nsurfMS.max(), a2prprofMS.min(), a2prprofMS.max()
    

#srcDir = '/home/utsumi/temp/EPC'
#srcPath= srcDir + '/GPM_EPC_002421_20140802_0726.NS_MS.nc'
#with h5py.File(srcPath) as h:
#    a3profjpl   = h['MS/precip_prof'][:]
#
#print a3profjpl.min(), a3profjpl.max()
#
#srcPath = '/home/utsumi/bin/PMM/ret-epc/GPM_EPC_012149_20160418_1228.NS_MS.nc'
#with h5py.File(srcPath) as h:
#    a3profjpl   = h['MS/precip_prof'][:]
#
#print a3profjpl.min(), a3profjpl.max()
#
#
