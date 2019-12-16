import glob
import os


#baseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/validprof/pair/gprof'
baseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/validprof/pair/epc.glb.v03.minrec1000.maxrec10000'

lvar = ['stoprad','top-stoppmw']
for var in lvar:
    for Year in [2014,2015]:
        for Mon in range(1,12+1):
            ssearch = baseDir + '/%04d/%02d/??/%s.??????.npy'%(Year,Mon,var)
            
            lsrcPath = glob.glob(ssearch)
            if len(lsrcPath)==0: continue
            

            for srcPath in lsrcPath:
                print srcPath
                #os.remove(srcPath)
