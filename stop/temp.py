from numpy import *
import numpy as np
import glob
atmp = np.load('/tank/utsumi/PMM/stop/data/t2m/2017/02/22/t2m.0dy.0dx.01surf.npy')
atqv = np.load('/tank/utsumi/PMM/stop/data/tqv/2017/02/22/tqv.0dy.0dx.01surf.npy')
astop= np.load('/tank/utsumi/PMM/stop/data/stop/2017/02/22/stop.0dy.0dx.01surf.npy')


lsrcPath = glob.glob('/tank/utsumi/PMM/stop/data/tqv/2017/*/*/tqv.*.npy')
lsrcPath = np.sort(lsrcPath)
for srcPath in lsrcPath:
    Year,Mon,Day = map(int, srcPath.split('/')[7:7+3])
    isurf = srcPath.split('.')[-2]

    tcPath = '/tank/utsumi/PMM/stop/data/Tc/%04d/%02d/%02d/Tc1.0dy.0dx.%s.npy'%(Year,Mon,Day,isurf)
    a = np.load(srcPath)
    tc= np.load(tcPath)
    if a.shape[0]==0: continue
    print Year,Mon,Day,a.shape, tc.shape
    if a.shape[0] != tc.shape[0]:
        print srcPath
        sys.exit()
