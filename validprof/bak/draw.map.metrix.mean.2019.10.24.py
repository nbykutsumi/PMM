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
import string
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
    figDir  = '/home/utsumi/temp/ret'

else:
    print 'check myhost'
    sys.exit()
#*******************************
#lseason=['JJA','DJF']
#lseason=['JJA']
#lseason=['DJF']
lseason=['JJADJF','JJA','DJF']
DB_MAXREC = 10000
DB_MINREC = 1000
expr = 'glb.v03.minrec%d.maxrec%d'%(DB_MINREC,DB_MAXREC)

ny,nx = 120,360
dvar = {}
#lrettype = ['epc','gprof']
#lrettype = ['gprof']
lrettype = ['epc']
#lvar = ['profrad','profpmw','top-profpmw']
lvar = ['profrad','profpmw']
#lvar = ['stoprad','top-stoppmw']
#lvar = ['precrad','precpmw']
#lvar = ['convfreqrad','stratfreqrad','convfreqpmw','stratfreqpmw']
#lvar = ['stoprad','top-stoppmw','convfreqrad','stratfreqrad','convfreqpmw','stratfreqpmw']
#lvar = ['zeroDegAltituderad']
#lthpr= [0.1, 10]
lthpr = [0.5]
#lthpr = [0.5,10]
#lthpr = [10]
#lthpr = [0.1,0.5]
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
def draw_map(a2dat, textcbartop=None):
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

    if textcbartop is not None:
        cax.text(0.1,1.05,textcbartop,fontsize=13)

    plt.savefig(figPath)
    print figPath
    plt.clf()  

#-- Elevation --------
a2orog = np.load(tankbaseDir + '/utsumi/PMM/validprof/const/orog.meter.sp.one.180x360.npy')
a2orog = a2orog[30:30+120,:]

#*******************************
for thpr in lthpr:
    for season in lseason:
        if   season=='JJADJF':
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
        elif season =='JJ':
            lYM = util.ret_lYM([2014,6],[2014,7])
        elif season =='DJ':
            lYM = util.ret_lYM([2014,12],[2015,1])

        else:
            print 'check season',season
            sys.exit()
    
        for rettype in lrettype:
            for var in lvar:
                if (rettype=='gprof')and(var=='top-profpmw'):
                    continue
                if (rettype=='gprof')and(var in ['stoprad','top-stoppmw']):
                    continue
    
    
                if var in ['profpmw','profrad','top-profpmw']:
                    nz = 25
                else:
                    nz = 1 
    
                #** Initialize ******
                if nz ==1:
                    a2sum = zeros([ny,nx],float32)
                    a2num = zeros([ny,nx], int32)
                else:
                    a3sum = zeros([ny,nx,nz-4], float32)  # 2-12.5km (500m layers)
                    a2num = zeros([ny,nx], int32)
            
                for Year,Mon in lYM:
                    if rettype =='epc': 
                        outDir = tankbaseDir + '/utsumi/PMM/validprof/map.pair/epc.%s'%(expr)
        
                    elif rettype =='gprof':
                        outDir = tankbaseDir + '/utsumi/PMM/validprof/map.pair/gprof'
        
                    else:
                        print 'check retype',rettype
        
                    util.mk_dir(outDir)
              
                    stamp  = 'pr%.1f.%04d%02d'%(thpr,Year,Mon)
                    if var in ['profpmw','profrad','top-profpmw']: 
                        sumPath= outDir  + '/%s.sum.%s.sp.one.npy' %(var, stamp)
                        numPath= outDir  + '/%s.num.%s.sp.one.npy' %(var, stamp)
                        sum2Path= outDir + '/%s.sum2.%s.sp.one.npy'%(var, stamp)
                
                        a3sumTmp = np.load(sumPath)[:,:,4:nz]  # 2-12.5km
                        a2numTmp = np.load(numPath)
    
                        a3sum = a3sum + a3sumTmp
                        a2num = a2num + a2numTmp
                    else:
                        sumPath= outDir  + '/%s.sum.%s.sp.one.npy' %(var,stamp)
                        numPath= outDir  + '/%s.num.%s.sp.one.npy' %(var,stamp)
                        sum2Path= outDir + '/%s.sum2.%s.sp.one.npy'%(var,stamp)
                
                        a2sumTmp = np.load(sumPath)
                        a2numTmp = np.load(numPath)
            
                        a2sum = a2sum + a2sumTmp
                        a2num = a2num + a2numTmp
                        print Mon,a2sum.max(), a2num.max() 
            
                if var in ['profpmw','profrad','top-profpmw']: 
                    dvar[var] = a3sum / a2num.reshape(ny,nx,1)
                else:
                    dvar[var] = a2sum / a2num.reshape(ny,nx)
    
                dvar[var] = ma.masked_invalid(dvar[var])
   
            ##*** Surface precipitation ***
            if 'zeroDegAltituderad' in lvar:
                # DPR
                a2var = dvar['zeroDegAltituderad'] * 0.001   # m-->km
                vmin,vmax= 0, 8
                mycm  = 'rainbow'
                stitle= 'Freezing level [km] %s\n(CMB)'%(season)
                figPath= figDir + '/mmap.freezh.%s.pr%.1f.%s.png'%(rettype,thpr,season)
                draw_map(a2var)
        

    
            ##*** Surface precipitation ***
            if 'precpmw' in lvar:
                # DPR
                a2var = dvar['precrad']
                vmin,vmax= 0, 4
                mycm  = 'rainbow'
                stitle= 'Surface precip (conditional) [mm/h] %s\n(CMB)'%(season)
                figPath= figDir + '/mmap.precrad.%s.pr%.1f.%s.png'%(rettype,thpr,season)
                draw_map(a2var)
        
                # EPC
                a2var = dvar['precpmw']
                vmin,vmax= 0, 4
                mycm  = 'rainbow'
                stitle= 'Surface precip (conditional) [mm/h] %s\n(%s)'%(season,string.upper(rettype))
                figPath= figDir + '/mmap.precpmw.%s.pr%.1f.%s.png'%(rettype,thpr,season)
                draw_map(a2var)
        
                # Difference
                a2var = dvar['precpmw'] - dvar['precrad']
                vmin,vmax= -1,1
                mycm  = 'Spectral_r'
                stitle= 'Surface precip difference [mm/h] %s\n(%s - CMB)'%(season,string.upper(rettype))
                figPath= figDir + '/mmap.prec.dif.%s.pr%.1f.%s.png'%(rettype,thpr,season)
                draw_map(a2var)
        
        
            #*** Convective / Stratiform ***
            if (rettype =='epc')and('convfreqrad' in lvar):
                vmin,vmax = 0, 1
                # DPR convective
                a2var = dvar['convfreqrad']
                vmin,vmax= None,None
                mycm  = 'rainbow'
                stitle= 'Convective pixel fraction (>%.1f mm/h) %s\n(DPR-Ku)'%(thpr, season)
                figPath= figDir + '/mmap.convpixelrad.pr%.1f.%s.png'%(thpr, season)
                draw_map(a2var)

                # DPR stratiform
                a2var = dvar['stratfreqrad']
                vmin,vmax= 0,1
                mycm  = 'rainbow'
                stitle= 'Stratiform pixel fraction (>%.1f mm/h) %s\n(DPR-Ku)'%(thpr, season)
                figPath= figDir + '/mmap.stratpixelrad.pr%.1f.%s.png'%(thpr, season)
                draw_map(a2var)

                # DPR convective
                a2var = dvar['convfreqpmw']
                vmin,vmax= 0,1
                mycm  = 'rainbow'
                stitle= 'Convective pixel fraction (>%.1f mm/h) %s\n(EPC top-ranked)'%(thpr, season)
                figPath= figDir + '/mmap.convpixelpmw.pr%.1f.%s.png'%(thpr, season)
                draw_map(a2var)

                # DPR stratiform
                a2var = dvar['stratfreqpmw']
                vmin,vmax= 0,1
                mycm  = 'rainbow'
                stitle= 'Stratiform pixel fraction (>%.1f mm/h) %s\n(EPC top-ranked)'%(thpr, season)
                figPath= figDir + '/mmap.stratpixelpmw.pr%.1f.%s.png'%(thpr, season)
                draw_map(a2var)

                # Difference convective
                a2var = ma.masked_less(dvar['convfreqpmw'],0) - ma.masked_less(dvar['convfreqrad'],0)
                a2var = ma.masked_invalid(a2var)
                vmin,vmax= -0.2,0.2
                mycm  = 'Spectral_r'
                stitle= 'Difference of conv pixel fraction (>%.1f mm/h) %s\n(EPC - CMB)'%(thpr, season)
                figPath= figDir + '/mmap.convpixel.diff.pr%.1f.%s.png'%(thpr, season)
                draw_map(a2var)

                # Difference stratiform
                a2var = ma.masked_less(dvar['stratfreqpmw'],0) - ma.masked_less(dvar['stratfreqrad'],0)
                a2var = ma.masked_invalid(a2var)
                vmin,vmax= -0.2,0.2
                mycm  = 'Spectral_r'
                stitle= 'Difference of strat pixel fraction (>%.1f mm/h) %s\n(EPC - CMB)'%(thpr, season)
                figPath= figDir + '/mmap.stratpixel.diff.pr%.1f.%s.png'%(thpr, season)
                draw_map(a2var)

           

 
                 
    
            #*** Storm top height ***
            if (rettype =='epc')and('stoprad' in lvar):
                # DPR
                a2var = dvar['stoprad']*0.001
                vmin,vmax= 0, 10
                mycm  = 'jet'
                stitle= 'Precip top height [km] (>%.1f mm/h) %s\n(DPR-Ku)'%(thpr, season)
                figPath= figDir + '/mmap.stoprad.pr%.1f.%s.png'%(thpr, season)
                draw_map(a2var)
        
                # EPC Top-rank
                a2var = dvar['top-stoppmw']*0.001
                vmin,vmax= 0,10
                mycm  = 'jet'
                stitle= 'Precip top height [km] (>%.1f mm/h) %s\n(EPC (Top-ranked))'%(thpr, season)
                figPath= figDir + '/mmap.stoppmw-top.pr%.1f.%s.png'%(thpr,season)
                draw_map(a2var)
        
                # Storm top height difference
                a2var = (dvar['top-stoppmw'] - dvar['stoprad'])*0.001
                vmin,vmax= -2,2
                mycm  = 'Spectral_r'
                stitle= 'Precip top height difference [km] (>%.1f mm/h) %s\n(EPC(Top ranked)-DPR/Ku)'%(thpr, season)
                figPath= figDir + '/mmap.stop.dif.pr%.1f.%s.png'%(thpr,season)
                draw_map(a2var)
    
    
            if 'profpmw' in lvar:
                #*** CC ********
                vmin,vmax=0,1 
                mycm  = 'jet'
                # DPR vs PMW
                a2cc  = calc_cc(dvar['profrad'], dvar['profpmw'], axis=2)
                stitle= 'Correlation of profles. %s\n(COMB & %s) %s'%(season, string.upper(rettype), expr)
                figPath= figDir + '/mmap.cc.prof.%s.pr%.1f.%s.png'%(rettype,thpr,season)
                draw_map(a2cc)
        
                # DPR vs PME top-ranked
                if (rettype=='epc')and('top-profpmw' in dvar.keys()):
                    a2cc  = calc_cc(dvar['profrad'], dvar['top-profpmw'], axis=2)
                    stitle= 'Correlation of profles. %s\n(COMB & %s top-ranked) %s'%(season, string.upper(rettype), expr)
                    figPath= figDir + '/mmap.cc.prof-top.%s.pr%.1f.%s.png'%(rettype,thpr,season)
                    draw_map(a2cc)
            
        
                
                #*** RMSE ******
                vmin,vmax=0,0.1
                mycm  = 'rainbow'
                a2rmse = calc_rmse(dvar['profrad'], dvar['profpmw'], axis=2)
                stitle= 'RMSE of prof. (g/m3) %s\n(COMB & %s)'%(season, string.upper(rettype))
                figPath= figDir + '/mmap.rmse.prof.%s.pr%.1f.%s.png'%(rettype,thpr, season)
                draw_map(a2rmse)
        
                # Top-ranked
                if (rettype=='epc')and('top-profpmw' in dvar.keys()):
                    a2rmse = calc_rmse(dvar['profrad'], dvar['top-profpmw'], axis=2)
                    stitle= 'RMSE of prof. (g/m3) %s\n(COMB & %s top-ranked)'%(season, string.upper(rettype))
                    figPath= figDir + '/mmap.rmse.prof-top.%s.pr%.1f.%s.png'%(rettype,thpr,season)
                    draw_map(a2rmse)
         
                
                
                #*** Peak height ******
                vmin,vmax=0,8
                mycm  = 'rainbow'
                # DPR
                a2ph = dvar['profrad'].argmax(axis=2) * 0.5
                a2mask= ma.masked_invalid( dvar['profrad'].max(axis=2) ).mask
                a2ph = ma.masked_where(a2mask, a2ph)
                stitle= 'Peak height (Above sea) (km) %s\n(COMB)'%(season)
                figPath= figDir + '/mmap.peakh-asl.profrad.%s.pr%.1f.%s.png'%(rettype,thpr,season)
                draw_map(a2ph)
                
        
                # PMW
                a2ph = dvar['profpmw'].argmax(axis=2) * 0.5
                a2mask= ma.masked_invalid( dvar['profpmw'].max(axis=2) ).mask
                a2ph = ma.masked_where(a2mask, a2ph)
                stitle= 'Peak height (Above sea) (km) %s\n(%s)'%(season, string.upper(rettype))
                figPath= figDir + '/mmap.peakh-asl.profpmw.%s.pr%.1f.%s.png'%(rettype,thpr,season)
                draw_map(a2ph)
                
        
                # Top PMW
                if (rettype =='epc')and('top-profpmw' in dvar.keys()):
                    a2ph = dvar['top-profpmw'].argmax(axis=2) * 0.5
                    a2mask= ma.masked_invalid( dvar['profpmw'].max(axis=2) ).mask
                    a2ph = ma.masked_where(a2mask, a2ph)
                    stitle= 'Peak height (Above sea) (km) %s\n(%s top-ranked)'%(season, string.upper(rettype))
                    figPath= figDir + '/mmap.peakh-asl.profpmw-top.%s.pr%.1f.%s.png'%(rettype,thpr,season)
                    draw_map(a2ph)
             

                # PMW
                vmin,vmax=-4,4
                mycm  = 'Spectral_r'
                a2ph = dvar['profpmw'].argmax(axis=2) * 0.5 - dvar['profrad'].argmax(axis=2) * 0.5
                a2mask1= ma.masked_invalid( dvar['profpmw'].max(axis=2) ).mask
                a2mask2= ma.masked_invalid( dvar['profrad'].max(axis=2) ).mask
                a2mask = a2mask1 + a2mask2
                a2ph = ma.masked_where(a2mask, a2ph)
                stitle= 'Peak height difference (km) %s\n(%s - CMB)'%(season,string.upper(rettype))
                figPath= figDir + '/mmap.peakh-asl.dif.%s.pr%.1f.%s.png'%(rettype,thpr,season)
                draw_map(a2ph)
                       
 
                
                ##*** Peak water content ***
                #vmin,vmax = 0,0.2
                #mycm = 'rainbow'
                ## DPR
                #a2pw = dvar['profrad'].max(axis=2)
                #stitle= 'Peak Wat. Cont (g/m3) %s\n(COMB)'%(season)
                #figPath= figDir + '/mmap.peakw.profrad.%s.pr%.1f.%s.png'%(rettype, thpr, season)
                #draw_map(a2pw)
                #
                ## PMW
                #a2pw = dvar['profpmw'].max(axis=2)
                #stitle= 'Peak Wat. Cont (g/m3) %s\n(%s) %s'%(season, string.upper(rettype), expr)
                #figPath= figDir + '/mmap.peakw.profpmw.%s.pr%.1f.%s.png'%(rettype, thpr, season)
                #draw_map(a2pw)
        
                ## PMW
                #if (rettype == 'epc')and('top-profpmw' in dvar.keys()):
                #    a2pw = dvar['top-profpmw'].max(axis=2)
                #    stitle= 'Peak Wat. Cont (g/m3) %s\n(%s top-ranked) %s'%(season, string.upper(rettype), expr)
                #    figPath= figDir + '/mmap.peakw.profpmw-top.%s.pr%.1f.%s.png'%(rettype, thpr, season)
                #    draw_map(a2pw)
                
                
                
                
                #*** Mean. precip. water ***
                #vmin,vmax = None,None
                vmin,vmax = 0,10
                mycm = 'rainbow'
                textcbartop = ur'($\times10^{-2}$)'
                # DPR
                a2ptot = dvar['profrad'].mean(axis=2)*100
                stitle= 'Mean prec. size hydrometeor (2-12km) (g/m3) %s\nCMB'%(season)
                figPath= figDir + '/mmap.totwat.rad.%s.pr%.1f.%s.png'%(rettype,thpr,season)
                draw_map(a2ptot,textcbartop=textcbartop)
               
                # PMW
                a2ptot = dvar['profpmw'].mean(axis=2)*100
                stitle= 'Mean prec. size hydrometeor (2-12km) (g/m3) %s\n%s'%(season,string.upper(rettype))
                figPath= figDir + '/mmap.totwat.pmw.%s.pr%.1f.%s.png'%(rettype,thpr,season)
                draw_map(a2ptot)
                draw_map(a2ptot,textcbartop=textcbartop)
                
       
                # Difference  
                a2var = (dvar['profpmw'].mean(axis=2) - dvar['profrad'].mean(axis=2))*100
                vmin,vmax= -5,5
                mycm  = 'Spectral_r'
                stitle= 'Mean precip. size hydrometeor difference (g/m3) %s\n(%s - CMB)'%(season, string.upper(rettype))
                figPath= figDir + '/mmap.totwat.dif.%s.pr%.1f.%s.png'%(rettype,thpr,season)
                draw_map(a2var,textcbartop=textcbartop)
 
