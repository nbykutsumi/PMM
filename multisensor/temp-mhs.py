# %%
import numpy as np
import os
import GSMaP_MHS_V7


#srcPath = '/home/utsumi/temp/multi/gsmap_l2/GPMCOR_GMI_1802010057_0229_000000_1BS_MMR_04C_retrv.dat'
#srcPath = '/home/utsumi/temp/multi/gsmap_l2/GPMCOR_GMI_1801010109_0242_000000_1BS_MMR_04C_retrv.dat.gz'
srcPath = '/home/utsumi/temp/multi/gsmap_l2/GPMMTA_MHS_1801010023_0207_000000_1BS_MMR_04C_retrv.dat'

#srcPath = '/home/utsumi/temp/multi/gsmap_l2/GPMCOR_GMI_1801010109_0242_000000_1BS_MMR_04C_retrv.dat.gz'
compressed = False
#compressed = True
#lvname = ['year','day','utc','startdate','enddate','lat','lon','sfc','si','lz','qc','snowm','rr']
lvname = ['utc']
nrec = 1
#origin = 2207 -1
origin = 0


gsmap = GSMaP_MHS_V7.level2()
#lout = gsmap.get_var(srcPath=srcPath, lvname=lvname, nrec=nrec, origin=origin, compressed=False)
lout = gsmap.get_var(srcPath=srcPath, lvname=lvname, compressed=False)
#print lout

# %%
