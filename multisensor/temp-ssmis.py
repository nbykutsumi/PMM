# %%
import numpy as np
import os
import GSMaP_SSMIS_V7


#srcPath = '/home/utsumi/temp/multi/gsmap_l2/GPMCOR_GMI_1802010057_0229_000000_1BS_MMR_04C_retrv.dat'
srcPath = './GPMF16_MIS_1801010009_0204_000000_1CS_MMR_04C_retrv.dat'

#srcPath = '/home/utsumi/temp/multi/gsmap_l2/GPMCOR_GMI_1801010109_0242_000000_1BS_MMR_04C_retrv.dat.gz'
compressed = False
#compressed = True
lvname = ['lon','lat','dtime','itoil','irflg','rainfg','snowprb']
nrec = 1

origin = 3568


gsmap = GSMaP_SSMIS_V7.level2()
#lout = gsmap.get_var(srcPath=srcPath, lvname=lvname, nrec=nrec, origin=origin, compressed=True)
#lout = gsmap.get_var(srcPath=srcPath, lvname=lvname,compressed=False)
lout = gsmap.get_var(srcPath=srcPath, lvname=lvname,compressed=False, nrec=1, origin=origin)
lon= lout[0]
lat= lout[1]
dtime=lout[2]
itoil=lout[3]
irflg=lout[4]
rainfg=lout[5]
snowprb=lout[6]
y,x = 0,127
print(lon[y,x],lat[y,x],dtime[y],itoil[y,x],rainfg[y,x])
#y,x = 0,1
#print(lon[y,x],lat[y,x],dtime[y],itoil[y,x],rainfg[y,x])
#y,x = 0,2
#print(lon[y,x],lat[y,x],dtime[y],itoil[y,x],rainfg[y,x])

# %%
