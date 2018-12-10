from numpy import *
import myfunc.util as util
import numpy as np
import os

domain   = 'FLORIDA-SFL-N'
#iYM = [2010,4]
iYM = [2010,4]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]

sout = 'rat,esurf,gv,gvave,dtime,lat,lon'+',' + ','.join(map(str,range(-15,30+1))) + '\n'
for YM in lYM:
    Year,Mon = YM
    srcDir   = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.L2A25/5.0km/%s/%04d%2d'%(domain,Year,Mon)
    gvPath   = srcDir + '/p_gvprcp.npy'
    profPath = srcDir + '/p_prof.npy'
    esurfPath= srcDir + '/p_eSurf.npy'
    latPath  = srcDir + '/p_sateLat.npy'
    lonPath  = srcDir + '/p_sateLon.npy'
    timePath = srcDir + '/p_dtime.npy'

    if not os.path.exists(gvPath): continue

    a2gv     = abs(np.load(gvPath))
    a2prof   = np.load(profPath)
    a2prof   = ma.masked_less(a2prof,0)*0.01
    a2prof   = a2prof.data
    a1esurf  = np.load(esurfPath)
    a1lat    = np.load(latPath)
    a1lon    = np.load(lonPath)
    a1dtime  = np.load(timePath)

    a1rat    = ma.masked_invalid(a2gv[:,15]/a1esurf)
    a1gvave   = ma.masked_less(a2gv[:,15-5:15-5+30],0).mean(axis=1)

    a1idx = arange(a2gv.shape[0])
    a1maskP   = ma.masked_less(a1esurf,2).mask
    a1maskRat = ma.masked_outside(a1rat, 0.8,1.2).mask

    a1mask   = a1maskP + a1maskRat
    a1idxTmp = ma.masked_where(a1mask, a1idx).compressed()


    a2gvTmp  =  a2gv[a1idxTmp,:]
    a2profTmp=  a2prof[a1idxTmp,:]
    a1esurfTmp= a1esurf[a1idxTmp]
    a1gvaveTmp= a1gvave[a1idxTmp]
    a1ratTmp  = a1rat[a1idxTmp]
    a1latTmp  = a1lat[a1idxTmp]
    a1lonTmp  = a1lon[a1idxTmp]
    a1dtimeTmp= a1dtime[a1idxTmp]
     
    for i in range(len(a1ratTmp)):
        l1=[a1ratTmp[i],a1esurfTmp[i],a2gvTmp[i,15],a1gvaveTmp[i],a1dtimeTmp[i],a1latTmp[i],a1lonTmp[i]]

        l2= map(str,a2gvTmp[i,:])
        l = l1 + l2
        stmp = ','.join(map(str,l)) + '\n'
        sout = sout + stmp



outDir  = '/home/utsumi/temp'
outPath = outDir + '/temp.csv'
f=open(outPath,'w'); f.write(sout); f.close()
print outPath
