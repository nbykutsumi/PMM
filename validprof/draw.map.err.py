import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import os, sys
from numpy import ma

Year,Mon = 2017,1
[[lllat,lllon],[urlat,urlon]] = [[-60,-180],[60,180]]
a2lon, a2lat = np.meshgrid(np.arange(lllon,urlon+0.01,1), np.arange(lllat,urlat+0.01,1))


srcDir = '/tank/utsumi/validprof/maperror'
#lvar = ['precbias','precbrat','profS','profrmse']
lvar = ['precbias','profS','precrad','precpmw']
#a2precbias = np.load(srcDir + '/ave.precbias.%04d%02d.npy'%(Year,Mon))
#a2precbrat = np.load(srcDir + '/ave.precbrat.%04d%02d.npy'%(Year,Mon))
#a2profS    = np.load(srcDir + '/ave.profS.%04d%02d.npy'%(Year,Mon))
#a2profrmse = np.load(srcDir + '/ave.profrmse.%04d%02d.npy'%(Year,Mon))


fig = plt.figure(figsize=(18,10))
for ivar, var in enumerate(lvar):
    a2var = np.load(srcDir + '/ave.%s.%04d%02d.npy'%(var, Year,Mon))
    a2fig= ma.masked_equal(a2var, -9999.)
    #a2fig= a2precbias
    print a2fig.min(), a2fig.max()

    if var=='precbias':
        ax  = fig.add_axes([0.05,0.1,0.35,0.35])
        vmin,vmax = [-3,3]
        mycm = 'PuOr_r'
    elif var=='profS':
        ax  = fig.add_axes([0.5,0.1,0.35,0.35])
        vmin,vmax = [0,1]
        mycm = 'jet'
    elif var=='precrad':
        ax  = fig.add_axes([0.05,0.5,0.35,0.35])
        vmin,vmax = [0,2]
        mycm = 'rainbow'
    elif var=='precpmw':
        ax  = fig.add_axes([0.5,0.5,0.35,0.35])
        vmin,vmax = [0,2]
        mycm = 'rainbow'



#    elif var=='precbrat':
#        ax  = fig.add_axes([0.1,0.5,0.3,0.3])
#        vmin,vmax = [-2,2]
#        mycm = 'PuOr_r'
#    elif var=='profrmse':
#        ax  = fig.add_axes([0.5,0.5,0.3,0.3])
#        vmin,vmax = [0,10]
#        mycm = 'jet'
    else:
        print 'check var=',var
        sys.exit()

    cmap = plt.get_cmap(mycm)
    cmap.set_bad(color='gray', alpha=1.)
    M   = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)
    im  = M.pcolormesh(a2lon, a2lat, a2fig, vmin=vmin, vmax=vmax, cmap=cmap)
    if var=='profS':
        plt.title('1 - taylorS')
    else:
        plt.title('%s'%(var))
    
    M.drawcoastlines(linewidth=1)
    
    M.drawmeridians(np.arange(0,360,10), labels=[0,0,0,1], fontsize=8, linewidth=0.5, fmt='%d',rotation=50)
    M.drawparallels(np.arange(-60,60,10), labels=[1,0,0,0], fontsize=8, linewidth=0.5, fmt='%d')
    plt.colorbar(im, orientation='horizontal')
    
    season  = '%02d'%(Mon)
figPath = '/home/utsumi/temp/validprof/map.err.%s.png'%(season)
plt.savefig(figPath)
print figPath
    
