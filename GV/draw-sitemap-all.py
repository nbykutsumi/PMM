import matplotlib
matplotlib.use('Agg')
from numpy import *
import myfunc.util as util
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
#import GPMGV

iYM = [2005,4]
#iYM = [2014,4]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]

#gv = GPMGV.GPMGV()
#gv.load_sitelist_reclassified()
#dgName = gv.ret_ddomYM2gName()
#ldomain = gv.domains
ldomain = ['FLORIDA-KSC','FLORIDA-SFL-N','TEXAS-HAR','N.Carolina-IPHEx_Duke','N.Carolina-IPHEx_NASA','KWAJALEIN-KWA']

dBBox = {
 'FLORIDA-KSC':[[28.42,-80.75],[28.70,-80.53]]
,'FLORIDA-SFL-N':[[25.21,-81.88],[28.31,-80.12]]
,'TEXAS-HAR':[[29.46,-95.89],[30.43,-94.99]]
,'N.Carolina-IPHEx_Duke':[[35.37,-83.26],[35.77,-82.92]]
,'N.Carolina-IPHEx_NASA':[[35.44,-83.59],[35.69,-82.57]]
,'KWAJALEIN-KWA':[[9.39,167.46],[9.40,167.48]]
}
#-- Read site location ------
listDir = '/work/a01/utsumi/data/GPMGV/sitelist'
listPath= listDir + '/sitelist-usedsite.csv'
f=open(listPath,'r'); lines=f.readlines(); f.close()
dllat = {}
dllon = {}
for line in lines:
    line = line.strip().split(',')
    domain = line[0]
    lat,lon= map(float, line[1:2+1])
    try:
        dllat[domain].append(lat)
        dllon[domain].append(lon)
    except KeyError:
        dllat[domain] = [lat]
        dllon[domain] = [lon]

#-- Draw hemisphere ---------
fig = plt.figure(figsize=(11,4))
ax1 = fig.add_axes([0.1,0.1,0.8,0.8])
[[lllat,lllon],[urlat,urlon]] = [[8,-180-18],[37,-76]]
M   = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax1)

for domain in ldomain:
    #if domain =='KWAJALEIN-KWA': continue
    [[lat0,lon0],[lat1,lon1]] = dBBox[domain]
    M.plot([lon0,lon1],[lat0,lat0],color='r')
    M.plot([lon0,lon1],[lat1,lat1],color='r')
    M.plot([lon0,lon0],[lat0,lat1],color='r')
    M.plot([lon1,lon1],[lat0,lat1],color='r')

    if domain == 'KWAJALEIN-KWA':
        M.plot(lon0-360,lat1,'.',color='r')

#- draw mainland domain --
[[lat0,lon0],[lat1,lon1]] = [[24,-96],[36,-80]]
M.plot([lon0,lon1],[lat0,lat0],color='k')
M.plot([lon0,lon1],[lat1,lat1],color='k')
M.plot([lon0,lon0],[lat0,lat1],color='k')
M.plot([lon1,lon1],[lat0,lat1],color='k')

#- draw kwajaleing domain -
[[lat0,lon0],[lat1,lon1]] = dBBox['KWAJALEIN-KWA']
lat0 = lat0 - 0.4
lat1 = lat1 + 0.4
lon0 = lon0 - 360 -0.4
lon1 = lon1 - 360 +0.4

M.plot([lon0,lon1],[lat0,lat0],color='k')
M.plot([lon0,lon1],[lat1,lat1],color='k')
M.plot([lon0,lon0],[lat0,lat1],color='k')
M.plot([lon1,lon1],[lat0,lat1],color='k')


#-------------------------
gridlatlon = 5  
parallels = arange(0,40+1,10)
meridians = arange(-200,-70,10)
M.drawparallels(parallels,labels=[1,0,0,0], dashes=[6,900],fontsize=15)
M.drawmeridians(meridians, labels=[0,0,0,1], dashes=[6,900],fontsize=15, rotation=60)

M.drawcoastlines()

figDir = '/work/a01/utsumi/GPMGV/fig'
figPath= figDir + '/sitemap-hemisphere.png'
plt.savefig(figPath)
print figPath




#-- Draw mainland ---------
#fig = plt.figure(figsize=(8,5))
fig = plt.figure(figsize=(8,5))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
[[lllat,lllon],[urlat,urlon]] = [[24,-96],[36,-80]]
M   = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)

for domain in ldomain:
    [[lat0,lon0],[lat1,lon1]] = dBBox[domain]
    if domain =='KWAJALEIN-KWA':
        continue
    else:
        if domain != 'FLORIDA-KSC':
            lat0 = lat0-0.2
    
        lon0 = lon0-0.2
    
        if domain != 'FLORIDA-SFL-N':
            lat1 = lat1+0.2
        lon1 = lon1+0.2

    M.plot([lon0,lon1],[lat0,lat0],color='r')
    M.plot([lon0,lon1],[lat1,lat1],color='r')
    M.plot([lon0,lon0],[lat0,lat1],color='r')
    M.plot([lon1,lon1],[lat0,lat1],color='r')

    #-- Plot gauges -
    Lat = dllat[domain]
    Lon = dllon[domain]
    M.plot(Lon,Lat,'.',color='r')

gridlatlon = 5  
parallels = arange(25,40+1,5)
meridians = arange(-95,-70,5)
M.drawparallels(parallels,labels=[1,0,0,0], dashes=[6,900],fontsize=15)
M.drawmeridians(meridians, labels=[0,0,0,1], dashes=[6,900],fontsize=15)

M.drawcoastlines()

figDir = '/work/a01/utsumi/GPMGV/fig'
figPath= figDir + '/sitemap-mainland.png'
plt.savefig(figPath)
print figPath



#-- Draw Kwajalein ---------
domain ='KWAJALEIN-KWA'
fig = plt.figure(figsize=(3,3))
ax = fig.add_axes([0.2,0.2,0.7,0.7])
[[lllat,lllon],[urlat,urlon]] = [[9.38,167.45],[9.41,167.48]]
M   = Basemap(resolution='f', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)

#-- Plot gauges -
Lat = dllat[domain]
Lon = dllon[domain]
M.plot(Lon,Lat,'.',color='r')
for i,lat in enumerate(Lat):
    print Lat[i],Lon[i]

gridlatlon = 0.02
parallels = [9.38,9.40]
meridians = [167.46, 167.48]
M.drawparallels(parallels,labels=[1,0,0,0], dashes=[6,900], fontsize=15)
M.drawmeridians(meridians, labels=[0,0,0,1], dashes=[6,900], fontsize=15)

M.drawcoastlines()

figDir = '/work/a01/utsumi/GPMGV/fig'
figPath= figDir + '/sitemap-Kwajalein.png'
plt.savefig(figPath)
print figPath






