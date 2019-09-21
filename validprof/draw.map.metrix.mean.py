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
    figDir      = '/home/utsumi/mnt/lab_home_rainbow/public_html/tempfig/validprof'

else:
    print 'check myhost'
    sys.exit()
#*******************************
lseason=['JJA']
ny,nx = 120,360
dprof = {}
lrettype = ['rad-epc','rad-gprof','epc','gprof']
#*******************************
def calc_cc(x,y,axis):
    ny,nx,nz = x.shape
    xm = x.mean(axis=axis).reshape(ny,nx,1)
    ym = y.mean(axis=axis).reshape(ny,nx,1)
    A  = ((x-xm)*(y-ym)).sum(axis=axis)
    B  = ((x-xm)**2).sum(axis=axis)
    C  = ((y-ym)**2).sum(axis=axis)
    return A/( np.sqrt(B*C) )

def calc_rmse(x,y,axis):
    ny,nx,nz = x.shape
    return np.sqrt(((x-y)**2).sum(axis=axis)/nz)

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
a2orog = np.load(tankbaseDir + '/utsumi/validprof/const/orog.meter.sp.one.180x360.npy')
a2orog = a2orog[30:30+120,:]

#*******************************
for season in lseason:
    if season=='JJA':
        lYM = util.ret_lYM([2014,6],[2014,8])
        #lYM = util.ret_lYM([2014,6],[2014,6])
    elif season=='SON':
        lYM = util.ret_lYM([2014,9],[2014,11])
    elif season=='DJF':
        lYM = util.ret_lYM([2014,12],[2015,2])
    elif season=='MAM':
        lYM = util.ret_lYM([2015,3],[2015,5])
 
    for rettype in lrettype:
        if rettype =='rad':
            #dirtype = 'gprof'
            dirtype = 'epc'
            filetype= 'rad'
            profname= 'prof'
            nz = 36
        elif rettype =='rad-epc':
            dirtype = 'epc'
            filetype= 'rad'
            profname= 'prof'
            nz = 36
        elif rettype =='rad-gprof':
            dirtype = 'gprof'
            filetype= 'rad'
            profname= 'prof'
            nz = 36
        elif rettype=='gprof':
            dirtype = 'gprof'
            filetype= 'pmw'
            profname= 'prof'
            nz = 36
        elif rettype=='epc':
            dirtype = 'epc'
            filetype= 'pmw'
            profname= 'prof'
            nz = 25
        elif rettype=='epc-top':
            dirtype = 'epc'
            filetype= 'pmw'
            profname= 'top-prof'
            nz = 25
   

        a3sum = zeros([ny,nx,25], float32)
        a2num = zeros([ny,nx], int32)
    
        for Year,Mon in lYM:
            outDir = tankbaseDir + '/utsumi/validprof/mapprof/%s'%(rettype)
            util.mk_dir(outDir)
        
            sumPath= outDir  + '/prof.sum.%04d%02d.sp.one.npy' %(Year,Mon)
            numPath= outDir  + '/prof.num.%04d%02d.sp.one.npy' %(Year,Mon)
            numprofPath= outDir  + '/prof.numprof.%04d%02d.sp.one.npy' %(Year,Mon)
            sum2Path= outDir + '/prof.sum2.%04d%02d.sp.one.npy'%(Year,Mon)
        
            a3sumTmp = np.load(sumPath)[:,:,:25]
            a2numTmp = np.load(numPath)
    
            a3sum = a3sum + a3sumTmp
            a2num = a2num + a2numTmp
    
    
        dprof[rettype] = a3sum / a2num.reshape(ny,nx,1)
    

    #*** CC ********
    # DPR vs EPC
    a2cc  = calc_cc(dprof['rad-epc'], dprof['epc'], axis=2)
    vmin,vmax=0,1 
    mycm  = 'jet'
    stitle= 'Correlation of profles. %s\n(COMB & EPC)'%(season)
    figPath= figDir + '/mmap.cc.prof.epc.%s.png'%(season)
    draw_map(a2cc)

    # DPR vs GPROF
    a2cc = calc_cc(dprof['rad-gprof'], dprof['gprof'], axis=2)
    stitle= 'Correlation of profiles %s\n(COMB & GPROF)'%(season)
    figPath= figDir + '/mmap.cc.prof.gprof.%s.png'%(season)
    draw_map(a2cc)


    ##*** RMSE ******
    #vmin,vmax=0,0.1
    #mycm  = 'jet'
    ## DPR vs EPC
    #a2rmse = calc_rmse(dprof['rad-epc'], dprof['epc'], axis=2)
    #stitle= 'RMSE of prof. (g/m3) %s\n(COMB & EPC)'%(season)
    #figPath= figDir + '/mmap.rmse.prof.epc.%s.png'%(season)
    #draw_map(a2rmse)

    ## DPR vs GPROF
    #a2rmse = calc_rmse(dprof['rad-gprof'], dprof['gprof'], axis=2)
    #stitle= 'RMSE of prof. (g/m3) %s\n(COMB & GPROF)'%(season)
    #figPath= figDir + '/mmap.rmse.prof.gprof.%s.png'%(season)
    #draw_map(a2rmse)

    #*** Peak height ******
    vmin,vmax=0,8
    mycm  = 'jet'
    # DPR
    a2ph = dprof['rad-epc'].argmax(axis=2) * 0.5
    a2mask= ma.masked_invalid( dprof['rad-epc'].max(axis=2) ).mask
    a2ph = ma.masked_where(a2mask, a2ph)
    stitle= 'Peak height (Above sea) (km) %s\n(COMB)'%(season)
    figPath= figDir + '/mmap.peakh-asl.prof.dpr.%s.png'%(season)
    draw_map(a2ph)

    a2ph = a2ph-a2orog*0.001
    stitle= 'Peak height (Above surface) (km) %s\n(COMB)'%(season)
    figPath= figDir + '/mmap.peakh-agl.prof.dpr.%s.png'%(season)
    draw_map(a2ph)

    # EPC
    a2ph = dprof['epc'].argmax(axis=2) * 0.5
    a2mask= ma.masked_invalid( dprof['epc'].max(axis=2) ).mask
    a2ph = ma.masked_where(a2mask, a2ph)
    stitle= 'Peak height (Above sea) (km) %s\n(EPC)'%(season)
    figPath= figDir + '/mmap.peakh-asl.prof.epc.%s.png'%(season)
    draw_map(a2ph)

    a2ph = a2ph-a2orog*0.001
    stitle= 'Peak height (Above surface) (km) %s\n(EPC)'%(season)
    figPath= figDir + '/mmap.peakh-agl.prof.epc.%s.png'%(season)
    draw_map(a2ph)

    # GPROF
    a2ph = dprof['gprof'].argmax(axis=2) * 0.5
    a2mask= ma.masked_invalid( dprof['epc'].max(axis=2) ).mask
    a2ph = ma.masked_where(a2mask, a2ph)
    stitle= 'Peak height (Above sea) (km) %s\n(GPROF)'%(season)
    figPath= figDir + '/mmap.peakh-asl.prof.gprof.%s.png'%(season)
    draw_map(a2ph)

    a2ph = a2ph-a2orog*0.001
    stitle= 'Peak height (Above surface) (km) %s\n(GPROF)'%(season)
    figPath= figDir + '/mmap.peakh-agl.prof.gprof.%s.png'%(season)
    draw_map(a2ph)

    #*** Peak water content ***
    #vmin,vmax = 0,0.2
    #mycm = 'jet'
    ## DPR
    #a2pw = dprof['rad-epc'].max(axis=2)
    #stitle= 'Peak Wat. Cont (g/m3) %s\n(COMB)'%(season)
    #figPath= figDir + '/mmap.peakw.prof.dpr.%s.png'%(season)
    #draw_map(a2pw)

    ## DPR
    #a2pw = dprof['epc'].max(axis=2)
    #stitle= 'Peak Wat. Cont (g/m3) %s\n(EPC)'%(season)
    #figPath= figDir + '/mmap.peakw.prof.epc.%s.png'%(season)
    #draw_map(a2pw)

    ## GPROF
    #a2pw = dprof['gprof'].max(axis=2)
    #stitle= 'Peak Wat. Cont (g/m3) %s\n(GPROF)'%(season)
    #figPath= figDir + '/mmap.peakw.prof.gprof.%s.png'%(season)
    #draw_map(a2pw)


    ##*** Tot. precip. water ***
    #vmin,vmax = 0,4
    #mycm = 'jet'
    ## DPR
    #a2ptot = dprof['rad-epc'][:,:,2:24+1].sum(axis=2)
    #stitle= 'Tot. Prec. Water (1-12km) (g/m2) %s\n(COMB)'%(season)
    #figPath= figDir + '/mmap.totwat.prof.dpr.%s.png'%(season)
    #draw_map(a2ptot)

    ## EPC
    #a2ptot = dprof['epc'][:,:,2:24+1].sum(axis=2)
    #stitle= 'Tot. Prec. Water (1-12km) (g/m2) %s\n(EPC)'%(season)
    #figPath= figDir + '/mmap.totwat.prof.epc.%s.png'%(season)
    #draw_map(a2ptot)

    ## GPROF
    #a2ptot = dprof['gprof'][:,:,2:24+1].sum(axis=2)
    #stitle= 'Tot. Prec. Water (1-12km) (g/m2) %s\n(GPROF)'%(season)
    #figPath= figDir + '/mmap.totwat.prof.gprof.%s.png'%(season)
    #draw_map(a2ptot)





 

    x = 180-20
    y = 60-10
    a0 = dprof['rad-epc'][y,x]
    a1 = dprof['epc'][y,x]
    a2 = dprof['gprof'][y,x]
    plt.plot(a0, np.arange(25)*0.5, '-',color='k',linewidth=3)
    plt.plot(a1, np.arange(25)*0.5, '-',color='k',linewidth=1)
    plt.plot(a2, np.arange(25)*0.5, '--',color='k',linewidth=1)
    plt.savefig(figDir + '/temp.png')
    plt.clf()
    print a0
    print a1
