import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import os, sys
from numpy import ma
import myfunc.util as util
from numpy import *

iYM = [2017,1]
eYM = [2017,12]
lYM = util.ret_lYM(iYM,eYM)
lbiaslev = [0,1,2,3]
lcmap = ['darkblue','darkseagreen','orange','crimson']

fig = plt.figure(figsize=(6,6))
ax  = fig.add_axes([0.15,0.15,0.75,0.75])

for ibiaslev, biaslev in enumerate(lbiaslev):
    for sgn in [1,-1]:
        a1sum = []
        a1sum2= []
        a1num = []

        for Year,Mon in lYM:
            sym   = '%04d.%02d'%(Year,Mon)
            srcDir = '/tank/utsumi/validprof/pr.vs.metrics'
            a1prbin = np.load(srcDir + '/taylor.by.bias.prbin.%s.npy'%(sym))
            abiasrange = np.load(srcDir + '/taylor.by.bias.biasrange.%s.npy'%(sym))
            if len(lbiaslev) != len(abiasrange):
                print 'check lbiaslevl'
                sys.exit()
        
            a1sumTmp = np.load( srcDir + '/taylor.sum.bias=%dx%d.%s.npy'%(sgn,biaslev, sym))
            a1sum2Tmp= np.load( srcDir + '/taylor.sum2.bias=%dx%d.%s.npy'%(sgn,biaslev, sym))
            a1numTmp = np.load( srcDir + '/taylor.num.bias=%dx%d.%s.npy'%(sgn,biaslev, sym))
    
            if a1sum ==[]:
                a1sum = a1sumTmp
                a1sum2= a1sum2Tmp
                a1num = a1numTmp
            else:
                a1sum  = a1sum  + a1sumTmp
                a1sum2 = a1sum2 + a1sum2Tmp
                a1num  = a1num  + a1numTmp
    
    
        a1ave = ma.masked_invalid(a1sum / a1num)
        if ibiaslev==len(lbiaslev)-1:
            if sgn==1:
                slabel = 'bias=%d mm/day<'%(abiasrange[ibiaslev,0])
            elif sgn==-1:
                slabel = None
        else:
            if sgn==1:
                slabel = 'bias=%d-%d mm/day'%(abiasrange[ibiaslev,0],abiasrange[ibiaslev,1])
            if sgn==-1:
                slabel = None

        cmap = lcmap[biaslev]
        if sgn==-1:
            linestyle='--'
        else:
            linestyle='-'    

        ax.plot(a1prbin, a1ave, linestyle=linestyle, color=cmap, label=slabel)

ax.set_ylabel('Similarity (S-Index)',fontsize=20)
ax.set_xlabel('Mean precipitation intensity (mm/day)',fontsize=20)

ax.set_xlim([0,20])
plt.legend()
figDir = '/tank/utsumi/hometemp/validprof'
figPath = figDir + '/plot.byBias.scan.prec.vs.taylor.png'
plt.savefig(figPath)
print figPath
