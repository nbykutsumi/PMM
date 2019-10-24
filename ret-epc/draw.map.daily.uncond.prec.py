import matplotlib
matplotlib.use('Agg')
import numpy as np
import os, sys
import myfunc.util as util
from numpy import *
import calendar, glob
from datetime import datetime, timedelta
import socket
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

#*******************************
myhost = socket.gethostname()
if myhost =='shui':
    tankbaseDir = '/tank'
    workbaseDir = '/work'
    figDir      = '/home.rainbow/utsumi/public_html/tempfig/validprof'

elif myhost == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    workbaseDir = '/home/utsumi/mnt/lab_work'
    #figDir      = '/home/utsumi/mnt/lab_home_rainbow/public_html/tempfig/ret'
    figDir      = '/home/utsumi/temp/ret'

else:
    print 'check myhost'
    sys.exit()
#*******************************
iDTime = datetime(2014,6,1)
eDTime = datetime(2014,6,15)
lDTime = util.ret_lDTime(iDTime,eDTime,timedelta(days=1))

[[ilat,ilon],[elat,elon]] = [[30,-110],[50,-80]]
#[[ilat,ilon],[elat,elon]] = [[-5,-80],[10,-45]]

ny,nx = 120,360
DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

dprof = {}
#lrettype = ['dpr','epc','gprof']
lrettype = ['dpr','epcNScmb','gprof']
#lrettype = ['epcNScmb']
#*******************************
def draw_map(a2dat):
    fig = plt.figure(figsize=(8,4))
    ax  = fig.add_axes([0.08,0.1,0.8,0.75])
    a1lat = np.arange(-59.5,59.5+0.01,1.0)
    a1lon = np.arange(-179.5,179.5+0.01,1.0)
    X,Y = np.meshgrid(a1lon,a1lat)
    M = Basemap(resolution='l', llcrnrlat=ilat, llcrnrlon=ilon, urcrnrlat=elat, urcrnrlon=elon, ax=ax)
    im  = M.pcolormesh(X, Y, a2dat, vmin=vmin, vmax=vmax, cmap=mycm)

    plt.title(stitle, fontsize=15)
    M.drawcoastlines(linewidth=1)

    M.drawmeridians(np.arange(-180,180+1,30), labels=[0,0,0,1], fontsize=10, linewidth=0.5, fmt='%d',rotation=50)
    M.drawparallels(np.arange(-60,60+1,30), labels=[1,0,0,0], fontsize=10, linewidth=0.5, fmt='%d')

    cax = fig.add_axes([0.89,0.21,0.02,0.55])
    cbar= plt.colorbar(im, orientation='vertical', cax=cax)
    cax.tick_params(labelsize=13)
    plt.savefig(figPath)
    print figPath
    plt.clf()  

#-- Elevation --------
a2orog = np.load(tankbaseDir + '/utsumi/validprof/const/orog.meter.sp.one.180x360.npy')
a2orog = a2orog[30:30+120,:]

#*******************************
for DTime in lDTime:
    Year,Mon,Day = DTime.timetuple()[:3]
    ddat = {} 
    dsum = {} # test
    dnum = {} # test
    for rettype in lrettype:
        if rettype in ['epcNScmb','epcNS']:
            srcDir = tankbaseDir + '/utsumi/validprof/map-daily-uncond/%s.%s'%(rettype, expr)
        else:
            srcDir = tankbaseDir + '/utsumi/validprof/map-daily-uncond/%s'%(rettype)
    
        sumPath= srcDir  + '/prec.sum.%04d%02d%02d.sp.one.npy' %(Year,Mon,Day)
        numPath= srcDir  + '/prec.num.%04d%02d%02d.sp.one.npy' %(Year,Mon,Day)
    
        a2sum = np.load(sumPath)
        a2num = np.load(numPath)
    
        ddat[rettype] = (ma.masked_where(a2num==0, a2sum) / a2num)
        ddat[rettype] = ddat[rettype] * 24  # mm/h --> mm/day

        #dsum[rettype] = ma.masked_less(a2sum,0)
        #dnum[rettype] = ma.masked_equal(a2num,0)

    #*** Precipitation rate ********
    vmin,vmax=0,50
    #mycm  = 'gist_rainbow'
    mycm  = 'rainbow'
    # DPR 
    a2dat  = ddat['dpr']
    stitle= 'precipitation rate (mm/day)  %04d/%02d/%02d\n(COMB)'%(Year,Mon,Day)
    figPath= figDir + '/mmap.uncond.prec.dpr.%04d%02d%02d.png'%(Year,Mon,Day)
    draw_map(a2dat)

    # EPC
    a2dat  = ddat['epcNScmb']
    stitle= 'precipitation rate (mm/day)  %04d/%02d/%02d\n(EPC(NScmb))'%(Year,Mon,Day)
    figPath= figDir + '/mmap.uncond.prec.%s.epcNScmb.%04d%02d%02d.png'%(expr,Year,Mon,Day)
    draw_map(a2dat)

    # GPROF
    a2dat  = ddat['gprof']
    stitle= 'precipitation rate (mm/day)  %04d/%02d/%02d\n(GPROF)'%(Year,Mon,Day)
    figPath= figDir + '/mmap.uncond.prec.gprof.%04d%02d%02d.png'%(Year,Mon,Day)
    draw_map(a2dat)

    # EPC-DPR
    vmin,vmax = -5,5
    a2dat  = ddat['epcNScmb'] - ddat['dpr']
    stitle= 'precipitation difference (mm/day) %04d/%02d/%02d\n(EPC(NScmb) - COMB)'%(Year,Mon,Day)
    figPath= figDir + '/mmap.uncond.difprec.epcNScmb.%s.%04d%02d%02d.png'%(expr,Year,Mon,Day)
    draw_map(a2dat)

    # GPROF-DPR
    vmin,vmax = -5,5
    a2dat  = ddat['gprof'] - ddat['dpr']
    stitle= 'precipitation difference (mm/day) %04d/%02d/%02d\n(GPROF - COMB)'%(Year,Mon,Day)
    figPath= figDir + '/mmap.uncond.difprec.gprof.%04d%02d%02d.png'%(Year,Mon,Day)
    draw_map(a2dat)


