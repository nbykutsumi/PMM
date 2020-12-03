# %%
import numpy as np
import os
import GSMaP_GMI_V7


srcPath = '/home/utsumi/temp/multi/gsmap_l2/GPMGW1_AM2_1801010146_0235_000000_1BS_MMR_04C_retrv.dat'

#srcPath = '/home/utsumi/temp/multi/gsmap_l2/GPMCOR_GMI_1801010109_0242_000000_1BS_MMR_04C_retrv.dat.gz'
compressed = False
#compressed = True
lvname = ['lon','lat','itoil','irflg']
nrec = 1
#origin = 2207 -1
origin = 2865 -1
#origin = 2864 -1


gmi = GSMaP_GMI_V7.level2()
lout = gmi.get_var(srcPath=srcPath, lvname=['lon','lat','dtime'], nrec=nrec, origin=origin, compressed=True)
print lout

# %%
