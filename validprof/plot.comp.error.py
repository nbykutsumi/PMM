from numpy import *
import myfunc.util as util
from datetime import datetime, timedelta
import os, sys
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

iDTime = datetime(2017,1,1)
eDTime = datetime(2017,1,10)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

baseDir = '/tank/utsumi/validprof/pair.gprof'

def taylor_index(a2ref, a2dat, miss):
    a2mask1 = ma.masked_less_equal(a2ref, miss).mask
    a2mask2 = ma.masked_less_equal(a2dat, miss).mask
    a2mask  = a2mask1 + a2mask2
    a2ref = ma.masked_where(a2mask, a2ref)
    a2dat = ma.masked_where(a2mask, a2dat)
 
    a1num = (~a2mask).sum(axis=1)
    a1stdref = a2ref.std(axis=1)
    a1stddat = a2dat.std(axis=1)
    a1mref   = a2ref.mean(axis=1).reshape(-1,1)
    a1mdat   = a2dat.mean(axis=1).reshape(-1,1)

    a1cov    = ((a2ref - a1mref)*(a2dat - a1mdat)).sum(axis=1)/a1num
    a1corr = a1cov / (a1stdref * a1stddat)
    corrmax= 1.0
 
    S = 4*(1.0+a1corr)**4 /((a1stdref/a1stddat + a1stddat/a1stdref)**2) / (1.0+corrmax)**4
    S = ma.masked_invalid(S).filled(miss)
    return S

a1profS    = []
a1profrmse = []
a1precbias = []
a1precbrat = []
a1tmp      = []
for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]
    srcDir = baseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
    ssearch = srcDir + '/precpmw.*.npy'
    lprecpmwPath = sort(glob.glob(ssearch))
    for precpmwPath in lprecpmwPath:
        oid = int(precpmwPath.split('.')[-2])
        print oid, precpmwPath
        a1precpmw = np.load(srcDir + '/precpmw.%06d.npy'%(oid))
        a1precrad = np.load(srcDir + '/precrad.%06d.npy'%(oid)) 
        a2profpmw = np.load(srcDir + '/profpmw.%06d.npy'%(oid))
        a2profrad = np.load(srcDir + '/profrad.%06d.npy'%(oid))  # miss=-9999.

        a2profrad = ma.masked_less(a2profrad,0)


        #-- Set condition on surface prec by radar --
        #a1precrad = ma.masked_less(a1precrad, 1)

        #-- profile error ----
        a1profrmseTmp = np.sqrt(np.square(a2profpmw - a2profrad).mean(axis=1))
        a1profSTmp    = 1-taylor_index(a2profrad, a2profpmw, miss=-9999.)
        #-- surface error ----
        a1precbiasTmp = np.abs(a1precpmw - a1precrad)
        a1precbratTmp = (a1precpmw - a1precrad)/ a1precrad
        #-- append -----------
        a1profS.append(a1profSTmp)
        a1profrmse.append(a1profrmseTmp)
        a1precbias.append(a1precbiasTmp)
        a1precbrat.append(a1precbratTmp)
        a1tmp.append(a1precrad)


a1profS    = np.concatenate(a1profS)
a1profrmse = np.concatenate(a1profrmse)
a1precbias = np.concatenate(a1precbias)
a1precbrat = np.concatenate(a1precbrat)
a1tmp      = np.concatenate(a1tmp)

#*** Figure ****************
figsize = (6,6)
labelsize=15
bratlim = (-1,5)


fig = plt.figure(figsize=figsize)
ax  = fig.add_axes([0.1,0.1,0.8,0.8])
ax.plot(a1precbrat, a1profrmse, '.')
plt.xlabel('Surface precip error(bias ratio)', fontsize=labelsize)
plt.ylabel('Profile error(RMSE)', fontsize=labelsize)

plt.ylim(0,5)
plt.xlim(bratlim)

figPath = '/home/utsumi/temp/validprof/plot.brat.rmse.png'
plt.savefig(figPath)
plt.clf()
print figPath

#--------------
fig = plt.figure(figsize=figsize)
ax  = fig.add_axes([0.1,0.1,0.8,0.8])
ax.plot(a1precbrat, a1profS, '.')
plt.xlabel('Surface precip error(bias ratio)', fontsize=labelsize)
plt.ylabel('Profile similarity(S-Index)', fontsize=labelsize)
plt.ylim(0,1)
plt.xlim(bratlim)

figPath = '/home/utsumi/temp/validprof/plot.brat.sindex.png'
plt.savefig(figPath)
plt.clf()
print figPath

#-- 2D-hist ---
fig = plt.figure(figsize=figsize)
ax  = fig.add_axes([0.1,0.1,0.8,0.8])
H, xbnd, ybnd = np.histogram2d(a1precbrat, a1profS, bins=[np.linspace(-1,10,22+1), np.linspace(0, 1, 20+1)])
H   = H/H.sum()
#H   = np.log10(ma.masked_equal(H,0))
X,Y = np.meshgrid(xbnd,ybnd)
#im=plt.pcolormesh(X,Y,H.T)
im=plt.pcolormesh(X,Y,H.T, vmin=0, vmax=0.01)

plt.xlabel('Surface precip error(bias ratio)', fontsize=labelsize)
plt.ylabel('Profile error(1- SIndex)', fontsize=labelsize)
plt.ylim(0,1)
plt.xlim(bratlim)
plt.colorbar(im)

figPath = '/home/utsumi/temp/validprof/hist2d.brat.sindex.png'
plt.savefig(figPath)
plt.clf()
print figPath



#--------------
fig = plt.figure(figsize=figsize)
ax  = fig.add_axes([0.1,0.1,0.8,0.8])
ax.plot(a1tmp, a1profS, '.')
plt.xlabel('Surface precip [mm/h]', fontsize=labelsize)
plt.ylabel('Profile error(1 - SIndex)', fontsize=labelsize)
plt.ylim(0,1)
figPath = '/home/utsumi/temp/validprof/plot.surfaceprec.sindex.png'
plt.savefig(figPath)
plt.clf()
print figPath

#--------------
fig = plt.figure(figsize=figsize)
ax  = fig.add_axes([0.1,0.1,0.8,0.8])
ax.plot(a1tmp, a1profrmse, '.')
plt.xlabel('Surface precip [mm/h]', fontsize=labelsize)
plt.ylabel('Profile error(RMSE)', fontsize=labelsize)
plt.ylim(0,5)
figPath = '/home/utsumi/temp/validprof/plot.surfaceprec.rmse.png'
plt.savefig(figPath)
plt.clf()
print figPath




