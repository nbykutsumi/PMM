from numpy import *
import os, sys
import Image
import myfunc.util as util


liYM = [[2012,4],[2011,4],[2010,4],[2009,4],[2008,4],[2007,4],[2006,4],[2005,4]]
eYM = [2014,10]
thdist = 5.0
minNum = 3
prdName = 'L2A25'
basepr = 'gv'
dattype= 'cc'
prtype = 'all'
iy  = 0  # top
ey  = -1
ix  = 0
ex  = -1


for rhtype in ['hum','all','dry']:

    da2dat = {}
    for i,iYM in enumerate(liYM):
        figDir  = '/work/a01/utsumi/GPMGV/fig/%04d.%02d-%04d.%02d'%(iYM[0],iYM[1],eYM[0],eYM[1])
        figPath = figDir + '/plot.dt-lev.RH.runmean.%s.%.1fkm.minNum.%d.base.%s.%s.%s.%s.png'%(prdName, thdist, minNum, basepr, dattype, rhtype, prtype)

        iimg   = Image.open(figPath)
        a2array= asarray(iimg)[iy:ey,ix:ex]
        da2dat[i] = a2array

    a2line0 = concatenate([da2dat[0],da2dat[1],da2dat[2],da2dat[3]],axis=1)
    a2line1 = concatenate([da2dat[4],da2dat[5],da2dat[6],da2dat[7]],axis=1)

    a2oarray= concatenate([a2line0, a2line1], axis=0)   

    oimg   = Image.fromarray(a2oarray) 
    outDir = '/work/a01/utsumi/GPMGV/fig'
    oPath  = outDir + '/join.tmp.RH.%s.png'%(rhtype)
   
    oimg.save(oPath)
    print oPath 
    
