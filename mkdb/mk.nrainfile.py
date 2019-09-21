from numpy import *
import myfunc.util as util
import numpy as np
import os,sys
import socket

iMon = 1
eMon = 12
NREC = 10000
varName = 'DPRGMI_NS_surfPrecipTotRate'  # surface precipitation

sensor = 'GMI'
myhost = socket.gethostname()
if myhost =="shui":
    baseDir = '/work/hk01/utsumi/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.%02d-%02d'%(NREC,iMon,eMon)
elif myhost =="well":
    baseDir = '/media/disk2/share/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.%02d-%02d'%(NREC,iMon,eMon)



#lepcid  = range(0,25*25*25) # 25*25*25=15625
lepcid  = range(0,29*29*29) # 29*29*29 = 24389

'''
 /* first six= Ntotal then Nrain exceeding 1 5 10 20 50 mm/hr when T2m < 278K */
 /* second six= same for when T2m > 278K */

Warm/Cold separation is not considered for now.
'''


#for epcid in lepcid[:200]:
for epcid in lepcid:
    srcDir  = baseDir + '/%s'%(varName)
    srcPath = srcDir + '/%05d.npy'%(epcid)
    if not os.path.exists(srcPath):
        continue

    avar     = np.load(srcPath)
    nwarmtot = len(avar)
    nwarm1   = ma.count( ma.masked_less_equal(avar,1) )
    nwarm5   = ma.count( ma.masked_less_equal(avar,5) )
    nwarm10  = ma.count( ma.masked_less_equal(avar,10) )
    nwarm20  = ma.count( ma.masked_less_equal(avar,20) )
    nwarm50  = ma.count( ma.masked_less_equal(avar,50) )

    lout = [nwarmtot,nwarm1,nwarm5,nwarm10,nwarm20,nwarm50,0,0,0,0,0,0]
    sout = '\t'.join(map(str, lout))

    print epcid,lout
    #-- Save to file --
    outDir  = baseDir + '/nrain'
    util.mk_dir(outDir)
    outPath = outDir + '/db_%05d.bin.nrain.txt'%(epcid)
    f=open(outPath,'w'); f.write(sout); f.close()
    print outPath
