import os, shutil,glob

#baseDir = '/tank/utsumi/PMM/validprof/pair/gprof'
baseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/validprof/pair/gprof'

ssearch = baseDir + '/????/??/??/stoprad.??????.npy'
lsrcPath = glob.glob(ssearch)
print ssearch
for srcPath in lsrcPath:
    oid = int(srcPath.split('.')[-2])
    outDir = os.path.dirname(srcPath)
    outPath= outDir + '/heightStormToprad.%06d.npy'%(oid)
    print ''
    print srcPath
    print outPath, oid
    os.rename(srcPath, outPath)
