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
    #figDir      = '/home.rainbow/utsumi/public_html/tempfig/validprof'
    figDir      = '/home/utsumi/temp/ret'

elif myhost == 'well':
    tankbaseDir = '/home/utsumi/mnt/lab_tank'
    workbaseDir = '/home/utsumi/mnt/lab_work'
    #figDir      = '/home/utsumi/mnt/lab_home_rainbow/public_html/tempfig/ret'
    figDir      = '/home/utsumi/temp/ret'

else:
    print 'check myhost'
    sys.exit()
#*******************************
#lseason=['JJA']
#lseason=['DJF']
lseason=['JJADJF']
ny,nx = 120,360
DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

dprof = {}
#lrettype = ['dpr','epc','gprof']
lrettype = ['dpr','epcNScmb','gprof']
#*******************************
def draw_map(a2dat):
    fig = plt.figure(figsize=(8,4))
    ax  = fig.add_axes([0.08,0.1,0.8,0.8])
    a1lat = np.arange(-59.5,59.5+0.01,1.0)
    a1lon = np.arange(-179.5,179.5+0.01,1.0)
    X,Y = np.meshgrid(a1lon,a1lat)
    M = Basemap(resolution='l', llcrnrlat=-60, llcrnrlon=-180, urcrnrlat=60, urcrnrlon=180, ax=ax)
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
a2orog = np.load(tankbaseDir + '/utsumi/PMM/validprof/const/orog.meter.sp.one.180x360.npy')
a2orog = a2orog[30:30+120,:]

#*******************************
for season in lseason:
    if season =='JJADJF':
        lYM = util.ret_lYM([2014,6],[2014,8]) + util.ret_lYM([2014,12],[2015,2])
    elif season=='JJA':
        lYM = util.ret_lYM([2014,6],[2014,8])
    elif season=='SON':
        lYM = util.ret_lYM([2014,9],[2014,11])
    elif season=='DJF':
        lYM = util.ret_lYM([2014,12],[2015,2])
    elif season=='MAM':
        lYM = util.ret_lYM([2015,3],[2015,5])
    elif type(season)==int:
        if season <6:
            lYM = [[2015,season]]
        else:
            lYM = [[2014,season]]
    if season=='JJ':
        lYM = util.ret_lYM([2014,6],[2014,7])
    if season=='DJ':
        lYM = util.ret_lYM([2014,12],[2015,1])



    ddat = {} 
    for rettype in lrettype:

        a2sum = zeros([ny,nx], float32)
        a2num = zeros([ny,nx], int32)
    
        for Year,Mon in lYM:
            if rettype in ['epcNScmb','epcNS']:
                srcDir = tankbaseDir + '/utsumi/PMM/validprof/map-daily-uncond/%s.%s'%(rettype, expr)
            else:
                srcDir = tankbaseDir + '/utsumi/PMM/validprof/map-daily-uncond/%s'%(rettype)
    
            iDay = 1 
            #eDay = 15
            eDay = calendar.monthrange(Year,Mon)[1]
            for Day in range(iDay,eDay+1):

                sumPath= srcDir  + '/prec.sum.%04d%02d%02d.sp.one.npy' %(Year,Mon,Day)
                numPath= srcDir  + '/prec.num.%04d%02d%02d.sp.one.npy' %(Year,Mon,Day)
            
                a2sumTmp = np.load(sumPath)
                a2numTmp = np.load(numPath)
        
                a2sum = a2sum + a2sumTmp
                a2num = a2num + a2numTmp
    
        ddat[rettype] = (ma.masked_where(a2num==0, a2sum) / a2num).filled(0.0)
        ddat[rettype] = ddat[rettype] * 24  # mm/h --> mm/day
    

    #*** Precipitation rate ********
    vmin,vmax=0,10
    #mycm  = 'gist_rainbow'
    mycm  = 'rainbow'
    # DPR 
    a2dat  = ddat['dpr']
    stitle= 'precipitation rate (mm/day)  %s\n(COMB)'%(season)
    figPath= figDir + '/mmap.uncond.prec.dpr.%s.png'%(season)
    draw_map(a2dat)

    # EPC
    a2dat  = ddat['epcNScmb']
    stitle= 'precipitation rate (mm/day)  %s\n(EPC(NScmb))'%(season)
    figPath= figDir + '/mmap.uncond.prec.epcNScmb.%s.%s.png'%(expr,season)
    draw_map(a2dat)

    # GPROF
    a2dat  = ddat['gprof']
    stitle= 'precipitation rate (mm/day)  %s\n(GPROF)'%(season)
    figPath= figDir + '/mmap.uncond.prec.gprof.%s.png'%(season)
    draw_map(a2dat)

    # EPC-DPR
    mycm  = 'Spectral_r'
    vmin,vmax = -4,4
    a2dat  = ddat['epcNScmb'] - ddat['dpr']
    stitle= 'precipitation difference (mm/day)  %s\n(EPC(NScmb) - COMB)'%(season)
    figPath= figDir + '/mmap.uncond.difprec.epcNScmb.%s.%s.png'%(expr,season)
    draw_map(a2dat)

    # GPROF-DPR
    mycm  = 'Spectral_r'
    vmin,vmax = -4,4
    a2dat  = ddat['gprof'] - ddat['dpr']
    stitle= 'precipitation difference (mm/day)  %s\n(GPROF - COMB)'%(season)
    figPath= figDir + '/mmap.uncond.difprec.gprof.%s.png'%(season)
    draw_map(a2dat)


