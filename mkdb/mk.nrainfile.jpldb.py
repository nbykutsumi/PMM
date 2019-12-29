from numpy import *
import JPLDB
import sys, os, socket
import numpy
import myfunc.util as util


myhost = socket.gethostname()
#lsensor = ['AMSR2','SSMIS','MHS','ATMS']
#lsensor = ['SSMIS','MHS','ATMS']
lsensor = ['AMSR2']

'''
 /* first six= Ntotal then Nrain exceeding 1 5 10 20 50 mm/hr when T2m < 278K */
 /* second six= same for when T2m > 278K */

Warm/Cold separation is not considered for now.
'''


#lepcid = range(29*29*29)
lepcid = range(23543,29*29*29)
#lepcid = range(1000)
for sensor in lsensor:
    if myhost == 'shui':
        baseDir = '/work/hk01/utsumi/PMM/JPLDB/EPC_DB/%s_EPC_DATABASE_TEST29'%(sensor)
    elif myhost=='well':
        baseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/JPLDB/EPC_DB/%s_EPC_DATABASE_TEST29'%(sensor)


    db = JPLDB.JPLDB(sensor)


    for epcid in lepcid:
        dbPath = baseDir + '/db_%05d.bin'%(epcid)
        if not os.path.exists(dbPath):
            continue

        db.set_file(dbPath)

        #print epcid, os.path.exists(dbPath)

        avar = db.get_var('precip_NS_cmb')

        nwarmtot = len(avar)
        nwarm1   = ma.count( ma.masked_less_equal(avar,1) )
        nwarm5   = ma.count( ma.masked_less_equal(avar,5) )
        nwarm10  = ma.count( ma.masked_less_equal(avar,10) )
        nwarm20  = ma.count( ma.masked_less_equal(avar,20) )
        nwarm50  = ma.count( ma.masked_less_equal(avar,50) )

        lout = [nwarmtot,nwarm1,nwarm5,nwarm10,nwarm20,nwarm50,0,0,0,0,0,0]
        sout = '\t'.join(map(str, lout))

        print sensor,epcid,lout
        #-- Save to file --
        outDir  = baseDir + '/nrain'
        util.mk_dir(outDir)
        outPath = outDir + '/db_%05d.bin.nrain.txt'%(epcid)
        f=open(outPath,'w'); f.write(sout); f.close()

