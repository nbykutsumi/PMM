import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import os, sys
from numpy import ma
import myfunc.util as util

iYM = [2017,1]
eYM = [2017,12]
lYM = util.ret_lYM(iYM,eYM)
#lBBox = [((-50,-30),(-40,-20)),((-5,-30),(5,-20)) ]
lBBox = [((-50,-170),(-40,-160)),((-5,-170),(5,-160)) ]
lvar = ['profS','precbias','precrad','precpmw']
#lvar = ['profS']
lat0 = -60
lon0 = -180
dlatlon=1.0
ddat = {}
for var in lvar:
    for BBox in lBBox:
        ldat = []
        for (Year,Mon) in lYM:
            [[lllat,lllon],[urlat,urlon]] = BBox
            srcDir = '/tank/utsumi/validprof/maperror'
            if var in ['precbias']:
                a2varrad = ma.masked_less_equal(np.load(srcDir + '/ave.%s.%04d%02d.npy'%('precpmw', Year,Mon)),-9999.)
                a2varpmw = ma.masked_less_equal(np.load(srcDir + '/ave.%s.%04d%02d.npy'%('precrad', Year,Mon)),-9999.)
                a2var = a2varpmw - a2varrad
                
            else:
                a2var = np.load(srcDir + '/ave.%s.%04d%02d.npy'%(var, Year,Mon))
            iy = int((lllat-lat0)/dlatlon)
            ey = int((urlat-lat0)/dlatlon)
            ix = int((lllon-lon0)/dlatlon)
            ex = int((urlon-lon0)/dlatlon)
            dat = ma.masked_less_equal(a2var[iy:ey+1, ix:ex+1] ,-9999.).mean()
            ldat.append(dat)
        ddat[var,BBox] = ldat


#--- Figure ----
fig = plt.figure(figsize=(8,5))
w = 0.7
h = 0.2
x0 = 0.2
for i,var in enumerate(lvar):
    y0 = 0.05 + i*h
    ax = fig.add_axes([x0,y0,w,h])

    for BBox in lBBox:
        dat = ddat[var,BBox]
        ax.plot(dat)
    ax.text(0.1,0.8, var, transform=ax.transAxes)

figDir = '/tank/utsumi/hometemp/validprof'
figPath= figDir + '/plot.ts.metrices.png'
plt.savefig(figPath)
print figPath



