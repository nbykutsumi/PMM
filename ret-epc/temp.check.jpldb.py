from numpy import *
import JPLDB
import socket
import os, sys

#sensor = 'AMSR2'
#sensor = 'SSMIS'
sensor = 'MHS'
#sensor = 'ATMS'
dbtype = 'jpl'
#** Constants ******
myhost = socket.gethostname()
if myhost =='shui':
    if dbtype=='jpl':
        dbDir   = '/tank/utsumi/PMM/JPLDB/EPC_DB/%s_EPC_DATABASE_TEST29'%(sensor)
        coefDir = '/tank/utsumi/PMM/JPLDB/EPC_COEF/%s'%(sensor)
        countDir= '/tank/utsumi/PMM/EPCDB/list'

    elif dbtype=='my':
        dbDir   = '/tank/utsumi/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
        coefDir = '/tank/utsumi/PMM/EPCDB/EPC_COEF/%s'%(sensor)
        countDir= '/tank/utsumi/PMM/EPCDB/list'

elif myhost == 'well':
    if dbtype=='jpl':
        dbDir   = '/home/utsumi/mnt/lab_tank/utsumi/PMM/JPLDB/EPC_DB/%s_EPC_DATABASE_TEST29'%(sensor)
        coefDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/JPLDB/EPC_COEF/%s'%(sensor)
        countDir= '/home/utsumi/mnt/lab_tank/utsumi/PMM/EPCDB/list'


    elif dbtype=='my':
        dbDir   = '/media/disk2/share/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
        coefDir = '/media/disk2/share/PMM/EPCDB/EPC_COEF/%s'%(sensor)
        countDir= '/media/disk2/share/PMM/EPCDB/list'
else:
    print 'check dtype',dtype
    sys.exit()


db = JPLDB.JPLDB(sensor)
#for idx_db in range(784,802+1):
for idx_db in range(20000,24000):

    dbPath = dbDir + '/db_%05d.bin'%(idx_db)
    if not os.path.exists(dbPath): continue

    db.set_file(dbPath)

    sat1 = db.get_var('satid')
    sat2 = db.get_var('satid2')


    print set(sat1), set(sat2)

    #d = {}
    #for ivar,var in enumerate(db.vars):

    #    #if ivar<27: continue
    #    #if ivar > 40: continue


    #    d[var] = db.get_var(var)

    #    print ''
    #    print '-----------------------------'
    #    print ivar,var
    #    print d[var][0]
    #    print d[var].shape

    #sys.exit()
