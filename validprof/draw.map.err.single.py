import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import os, sys
from numpy import ma
import myfunc.util as util
from numpy import *

Year = 2017
lseason = ['ALL']
#lseason = ['DJF','MAM','JJA','SON']
#lseason = range(1,12+1)
#lseason = [1]
[[lllat,lllon],[urlat,urlon]] = [[-60,-180],[60,180]]
a2lon, a2lat = np.meshgrid(np.arange(lllon,urlon+0.01,1), np.arange(lllat,urlat+0.01,1))

ny,nx = 120,360
srcDir = '/tank/utsumi/validprof/maperror'
#lvar = ['precbias','precbrat','profS','profrmse']
#lvar = ['stoprad','stoppmw']
lvar = ['stoprad','stoppmw','peakhrad','peakhpmw']
#a2precbias = np.load(srcDir + '/ave.precbias.%04d%02d.npy'%(Year,Mon))
#a2precbrat = np.load(srcDir + '/ave.precbrat.%04d%02d.npy'%(Year,Mon))
#a2profS    = np.load(srcDir + '/ave.profS.%04d%02d.npy'%(Year,Mon))
#a2profrmse = np.load(srcDir + '/ave.profrmse.%04d%02d.npy'%(Year,Mon))


for season in lseason:
    fig = plt.figure(figsize=(9,3.5))
    for ivar, var in enumerate(lvar):
        lMon = util.ret_lmon(season)

        if var == 'precbias':
            a2sumrad = np.zeros([ny,nx],float32)
            a2sumpmw = np.zeros([ny,nx],float32)
            a2numrad = np.zeros([ny,nx],int32)
            a2numpmw = np.zeros([ny,nx],int32)

        else:
            a2sum = np.zeros([ny,nx],float32)
            a2num = np.zeros([ny,nx],int32)

        for Mon in lMon:
            if var=='precbias':
                a2varrad = np.load(srcDir + '/ave.%s.%04d%02d.npy'%('precrad', Year,Mon))
                a2varpmw = np.load(srcDir + '/ave.%s.%04d%02d.npy'%('precpmw', Year,Mon))
                a2numradTmp = np.load(srcDir + '/num.%s.%04d%02d.npy'%('precrad', Year,Mon))
                a2numpmwTmp = np.load(srcDir + '/num.%s.%04d%02d.npy'%('precpmw', Year,Mon))

                a2sumrad = a2sumrad + a2varrad * a2numradTmp
                a2sumpmw = a2sumpmw + a2varpmw * a2numpmwTmp
                a2numrad = a2numrad + a2numradTmp
                a2numpmw = a2numpmw + a2numpmwTmp

         
            else:     
                a2var = np.load(srcDir + '/ave.%s.%04d%02d.npy'%(var, Year,Mon))
                a2numTmp = np.load(srcDir + '/num.%s.%04d%02d.npy'%(var, Year,Mon))
                a2sum = a2sum + a2var * a2numTmp
                a2num = a2num + a2numTmp

       
        if var =='precbias':
            a2precrad = ma.masked_where(a2numrad==0, a2sumrad)/a2numrad
            a2precpmw = ma.masked_where(a2numpmw==0, a2sumpmw)/a2numpmw
            a2fig = a2precpmw - a2precrad
        else:
            a2fig = ma.masked_where(a2num==0, a2sum)/a2num
       
 
        a2fig= ma.masked_equal(a2fig, -9999.)
        ax  = fig.add_axes([0.1,0.15,0.8,0.75])
        print 'min,max=',a2fig.min(), a2fig.max()
    
        if var=='precbias':
            vmin,vmax = [-1,1]
            mycm = 'PuOr_r'
        elif var=='profS':
            vmin,vmax = [0,0.7]
            mycm = 'jet'
        elif var=='precrad':
            vmin,vmax = [0,2]
            mycm = 'rainbow'
        elif var=='precpmw':
            vmin,vmax = [0,2]
            mycm = 'rainbow'
        elif var in ['stoprad','stoppmw']:
            vmin,vmax = [0,12]
            mycm = 'rainbow'
        elif var in ['peakhrad','peakhpmw']:
            vmin,vmax = [0,12]
            mycm = 'rainbow'
    
        else:
            print 'check var=',var
            sys.exit()
    
        cmap = plt.get_cmap(mycm)
        cmap.set_bad(color='gray', alpha=1.)
        M   = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)
        im  = M.pcolormesh(a2lon, a2lat, a2fig, vmin=vmin, vmax=vmax, cmap=cmap)
        stamp = '%s'%(season)
        if var=='profS':
            stitle = 'taylorS' + ' '+stamp
        else:
            stitle = '%s'%(var)+ ' '+stamp
       
        plt.title(stitle, fontsize=15) 
        M.drawcoastlines(linewidth=1)
        
        M.drawmeridians(np.arange(0,360,10), labels=[0,0,0,1], fontsize=10, linewidth=0.5, fmt='%d',rotation=50)
        M.drawparallels(np.arange(-60,60,10), labels=[1,0,0,0], fontsize=10, linewidth=0.5, fmt='%d')

        cax = fig.add_axes([0.91, 0.2, 0.02, 0.64]) 
        cbar= plt.colorbar(im, orientation='vertical', cax=cax)
        cax.tick_params(labelsize=13) 
 
        figPath = '/home/utsumi/temp/validprof/map.%s.%s.png'%(var,season)
        plt.savefig(figPath)
        print figPath
        plt.clf() 
