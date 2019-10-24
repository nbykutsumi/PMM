#! /opt/local/bin/python

from mpl_toolkits.basemap import Basemap, cm
# requires netcdf4-python (netcdf4-python.googlecode.com)
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib.pyplot as plt
#import h5py
import sys
import os
import matplotlib.dates as mdate
import datetime 

ncfile = sys.argv[1]
nc = NetCDFFile(ncfile, 'r')

#ipix = 200
#jpix = 110
#ipix = 32
#jpix = 61

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
