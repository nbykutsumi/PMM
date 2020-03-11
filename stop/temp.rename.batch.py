import os, shutil

srcDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/EPCDB/wetcase.samp.5000.GMI.V05A.S1.ABp103-117/batch.epc'
for i in range(10):
    iPath = srcDir + '/epc.%02d.npy'%(i)
    oPath = srcDir + '/%02d.npy'%(i)
    os.rename(iPath, oPath) 
    print oPath
