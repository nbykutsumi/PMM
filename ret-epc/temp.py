# %%
import numpy as np
import glob
import os
srcdir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retepc/us.v01/AMSR2.GCOMW1/2018/01/20'

search = srcdir + '/*.030200.*.npy'
lsrcpath = glob.glob(search)
print sorted(lsrcpath)

for srcpath in lsrcpath:
    var = os.path.basename(srcpath).split('.')[0]
    a=np.load(srcpath)
    print var,a.min(),a.max()

# %%
