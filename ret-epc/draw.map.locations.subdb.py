import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import *
import numpy as np
import myfunc.util as util
import os, sys, glob, socket
from mpl_toolkits.basemap import Basemap, cm
import random
#lidx_db = [1681] # bad
lidx_db = [10820] # good
#lidx_db = [10849] # good
#lidx_db = [17954] # good
#lidx_db = [10878] # good, but bimodal
#lidx_db = [13198] # good
#lidx_db = [20129] # good, but not good for extreme
#lidx_db = [23547] # bad

nsample = 1000
DB_MAXREC = 10000
DB_MINREC = 1000

dbtype  = 'my'
sensor  = 'GMI'
expr = 'org.smp%d'%(nsample)
#********
myhost = socket.gethostname()
if myhost =='shui':
    dbDir   = '/work/hk01/utsumi/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
    countDir= '/work/hk01/utsumi/PMM/EPCDB/list'
    retbaseDir = '/tank/utsumi/PMM/retsynt/%s'%(expr)
    figDir  = '/home/utsumi/temp/ret'

elif myhost == 'well':
    dbDir   = '/media/disk2/share/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
    countDir= '/media/disk2/share/PMM/EPCDB/list'
    retbaseDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/retsynt/%s'%(expr)
    figDir  = '/home/utsumi/temp/ret'


for idx_db in lidx_db:
    retDir = retbaseDir + '/%05d'%(idx_db)
    a1est = np.load( retDir + '/nsurfNScmb.est.%05d.npy'%(idx_db))
    a1obs = np.load( retDir + '/nsurfNScmb.obs.%05d.npy'%(idx_db))
    a1latobs = np.load( retDir + '/Latitude.obs.%05d.npy'%(idx_db))
    a1lonobs = np.load( retDir + '/Longitude.obs.%05d.npy'%(idx_db))
    a2mdhmsobs =np.load( retDir + '/mdhms.obs.%05d.npy'%(idx_db))
    a1monobs = a2mdhmsobs[:,0]



    fig = plt.figure()
    plt.plot(a1obs, a1est, '.',color='k')
    plt.ylim([0,50])
    plt.xlim([0,50])
    plt.savefig(figDir + '/temp.png')
    plt.clf()
    #**** Sampling ***********
    a1idx  = np.arange(len(a1obs)).astype(int32)
    a1bias = ma.masked_less(a1est,0) - ma.masked_less(a1obs,0)

    #a1mask = ma.masked_inside(a1bias, -20,20).mask
    #a1idxTmp= ma.masked_where(a1mask, a1idx).compressed()

    a1idxTmp= random.sample(a1idx, 30) 


    a1latTmp = a1latobs[a1idxTmp]
    a1lonTmp = a1lonobs[a1idxTmp] 
    a1monTmp = a1monobs[a1idxTmp] 
    a1biasTmp= a1bias[a1idxTmp]
    
    fig = plt.figure(figsize=(16,8),dpi=100)
    #M = Basemap(projection='robin', resolution = 'l', area_thresh = 1000.0,
    #              lat_0=0, lon_0=0)

    M = Basemap(llcrnrlat=-67, llcrnrlon=-180, urcrnrlat=67, urcrnrlon=180, resolution='l')
    M.drawcoastlines()
    M.drawcountries()
    M.fillcontinents(color = 'gray')
    M.drawmapboundary()
    M.drawmeridians(np.arange(0, 360, 30))
    M.drawparallels(np.arange(-90, 90, 30))

    def get_marker_color(time):
        return ('ro')
    
    
    for i in range(a1latTmp.shape[0]):
        lat= a1latTmp[i]
        lon= a1lonTmp[i]
        #if ( lon < 0 ): lon+= 360
        mm= a1monTmp[i]
        bias    = a1biasTmp[i]
        absbias = abs(a1biasTmp[i])
      
        if ( mm > 11 or mm < 3 ):
          mycm = 'b'
        if ( mm > 2 and mm < 6 ):
          mycm = 'y'
        if ( mm > 5 and mm < 9 ):
          mycm = 'r'
        if ( mm > 8 and mm < 12 ):
          mycm = 'c'
        if   ( absbias <25):   msize = 10
        elif ( absbias <30  ): msize = 14
        elif ( absbias <40  ): msize = 16
        else:             msize = 20

        if bias>0 : fcolor='None'; mystyle='o'
        else:       fcolor='None'; mystyle='v'
        #x,y = M(lon, lat)
        #marker_string = get_marker_color(sat)
        #marker_string = 'ro'
        #M.plot(x, y, marker_string, markersize=msize)
        x,y = M(lon, lat)
        #msize = 50
        #mycm  = 'r'
        #fcolor = 'r'
        
        M.plot(x, y, marker=mystyle, markeredgecolor=mycm, markerfacecolor=fcolor, markeredgewidth=3, markersize=msize)
        print lon,lat,bias 


    x,y = M(80, -55)
    plt.text(x, y, "Marker size proportional to weight",
            verticalalignment='bottom', horizontalalignment='right',
            color='black', fontsize=12, fontname='Arial')
     
    x,y = M(150, -65)
    plt.text(x, y, 'Jun-Jul-Aug',
            verticalalignment='bottom', horizontalalignment='right',
            color='red', fontsize=18, fontname='Arial', weight='bold')
    x,y = M(-100, -65)
    plt.text(x, y, 'Sep-Oct-Nov',
            verticalalignment='bottom', horizontalalignment='right',
            color='cyan', fontsize=18, fontname='Arial', weight='bold')
    x,y = M(-20, -65)
    plt.text(x, y, 'Dec-Jan-Feb',
            verticalalignment='bottom', horizontalalignment='right',
            color='blue', fontsize=18, fontname='Arial', weight='bold')
    x,y = M(60, -65)
    plt.text(x, y, 'Mar-Apr-May',
            verticalalignment='bottom', horizontalalignment='right',
            color='yellow', fontsize=18, fontname='Arial', weight='bold')
    
    #stitle = 'DB#:%05d (Large bias entries)'%(idx_db)
    stitle = 'DB#:%05d'%(idx_db)
    plt.title(stitle, fontname='Arial', fontsize=18)
     
    figPath = figDir + '/map.loc.subdb.%05d.png'%(idx_db)
    plt.savefig(figPath, transparent='True', bbox_inches='tight', pad_inches=0.05)
    print figPath
    #plt.savefig('db_locations.png', transparent='True', bbox_inches='tight', pad_inches=0.05)
    #plt.show()


 


'''
plats= 1.0*np.array(nc['latitude'])
plons= 1.0*np.array(nc['longitude'])
tbs= 1.0*np.array(nc['Tb'])
rains= 1.0*np.array(nc['NS/precip'])
lats= 1.0*np.array(nc['NS/lat_1'])
lons= 1.0*np.array(nc['NS/lon_1'])
wts= 1.0*np.array(nc['NS/wt_1'])
mms= 1.0*np.array(nc['NS/mm_1'])
scs= 1.0*np.array(nc['NS/sfc_class_1'])

clat = getattr(nc,'center_lat')
clon = getattr(nc,'center_lon')
cdate = getattr(nc,'center_date')
satname = getattr(nc,'satname')
sensor = getattr(nc,'sensor')
rev= getattr(nc,'orbit_rev')

cdate2 = cdate.split()[0]
ctime2 = cdate.split()[1].split(':')
cdate = cdate2 + ' ' + ctime2[0] + ctime2[1] + ' UTC' 


nx = lats.shape[0]
ny = lats.shape[1]
ndb = lats.shape[2]
print (lats.shape)
print (cdate2, ctime2, rev)

nc.close()

lats[lats < -90]= np.nan
lons[lons < -180]= np.nan
wts[wts < 0]= np.nan
mms[mms < 0]= np.nan
scs[scs < 0]= np.nan

ipix= jpix= -1
rmax= -1.0
for i in range(0,nx):
  for j in range(0,ny):
  #for j in range(90,135):
    if ( rains[i,j] > rmax ):
      rmax= rains[i,j]
      ipix= i
      jpix= j

if ( ipix < 0 ):
  sys.exit()

#ipix = 135

ipix = 65
#ipix = 32
jpix = 121

#ipix = 36
#jpix = 61

rmax= round(rains[ipix,jpix],2)

plat= plats[ipix,jpix]
plon= plons[ipix,jpix]

print("Max= ",rmax," located at ",plat,plon,ipix,jpix)
#print(scs[ipix,jpix,:])

cls= np.zeros(15)
for j in range(0,15): cls[j]= 0
for i in range(0,ndb):
  sc= scs[ipix,jpix,i]
  for j in range(1,15):
    if ( sc == j ): cls[j]+= 1
for j in range(0,15): 
  print("cls= ",j,int(cls[j]))

print("TB= ",tbs[ipix,jpix,:])

fig = plt.figure(figsize=(16,8),dpi=100)
map = Basemap(projection='robin', resolution = 'l', area_thresh = 1000.0,
              lat_0=0, lon_0=0)
map.drawcoastlines()
map.drawcountries()
map.fillcontinents(color = 'gray')
map.drawmapboundary()
map.drawmeridians(np.arange(0, 360, 30))
map.drawparallels(np.arange(-90, 90, 30))
 

def get_marker_color(time):
  return ('ro')

#x,y = map(plon, plat)
#msize = 8
#marker_string = 'go'
#map.plot(x, y, marker_string, markersize=msize)

for i in range(0,ndb):
  lat= lats[ipix,jpix,i]
  lon= lons[ipix,jpix,i]
  #if ( lon < 0 ): lon+= 360
  mm= mms[ipix,jpix,i]
  wt= wts[ipix,jpix,i]

  if ( mm > 11 or mm < 3 ):
    marker_string= 'bo'
  if ( mm > 2 and mm < 6 ):
    marker_string= 'yo'
  if ( mm > 5 and mm < 9 ):
    marker_string= 'ro'
  if ( mm > 8 and mm < 12 ):
    marker_string= 'co'
  msize = 8
  if ( wt > 0.2 ): msize = 10
  if ( wt > 0.4 ): msize = 14
  if ( wt > 0.6 ): msize = 16
  if ( wt > 0.9 ): msize = 20
  x,y = map(lon, lat)
  #marker_string = get_marker_color(sat)
  #marker_string = 'ro'
  map.plot(x, y, marker_string, markersize=msize)

x,y = map(plon, plat)
map.plot(x, y, 'kP', markersize=16)

x,y = map(80, -55)
plt.text(x, y, "Marker size proportional to weight",
        verticalalignment='bottom', horizontalalignment='right',
        color='black', fontsize=12, fontname='Arial')
 
x,y = map(150, -87)
plt.text(x, y, 'Jun-Jul-Aug',
        verticalalignment='bottom', horizontalalignment='right',
        color='red', fontsize=16, fontname='Arial', weight='bold')
x,y = map(-100, -87)
plt.text(x, y, 'Sep-Oct-Nov',
        verticalalignment='bottom', horizontalalignment='right',
        color='cyan', fontsize=16, fontname='Arial', weight='bold')
x,y = map(-20, -87)
plt.text(x, y, 'Dec-Jan-Feb',
        verticalalignment='bottom', horizontalalignment='right',
        color='blue', fontsize=16, fontname='Arial', weight='bold')
x,y = map(60, -87)
plt.text(x, y, 'Mar-Apr-May',
        verticalalignment='bottom', horizontalalignment='right',
        color='yellow', fontsize=16, fontname='Arial', weight='bold')

title_string = "Locations of top DB candidates  " + satname + " " + sensor + "  Rev=" + rev + '  ' + cdate + '  (' + str(plat) + ', ' + str(plon) + ')' + "  R=" + str(rmax)
#title_string = "Similar locations to GMI on " + cdate + ' ' + ctime + '  (' + str(plat) + ', ' + str(plon) + ')' + '  Ts=' + str(pts) + '  Vap=' + str(ptqv) + '  Class=' + str(pcls)
#title_string += "%s through %s" % (timestrings[-1], timestrings[0])
plt.title(title_string, fontname='Arial')
 
outfile = 'db_locations_' + str(ipix) + '_' + str(jpix) + '.png'
plt.savefig(outfile, transparent='True', bbox_inches='tight', pad_inches=0.05)
#plt.savefig('db_locations.png', transparent='True', bbox_inches='tight', pad_inches=0.05)
#plt.show()
'''
