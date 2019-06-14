#! /opt/local/bin/python

from mpl_toolkits.basemap import Basemap, cm
# requires netcdf4-python (netcdf4-python.googlecode.com)
from netCDF4 import Dataset as NetCDFFile
from math import sin, cos, sqrt, atan2, radians
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import os

#  Map of all four retrievals from EPC

outfile = "ts_compare.png"
if ( len(sys.argv) == 3 ):
  outfile = sys.argv[2]


ncfile = sys.argv[1]
f = NetCDFFile(ncfile)
for attr in f.ncattrs():
  print (attr, '=', getattr(f, attr))

lats = 1.0*np.array(f['latitude'])
lons = 1.0*np.array(f['longitude'])
tb= 1.0*np.array(f['Tb'])


ts1= 1.0*np.array(f['NS/ts'])
t2m1= 1.0*np.array(f['NS/t2m'])
tqv1= 1.0*np.array(f['NS/tqv'])
frz1= 0.001*np.array(f['NS/h273'])

ts2= 1.0*np.array(f['MS/ts'])
t2m2= 1.0*np.array(f['MS/t2m'])
tqv2= 1.0*np.array(f['MS/tqv'])
frz2= 0.001*np.array(f['MS/h273'])

ts3= 1.0*np.array(f['MERRA2/ts'])
t2m3= 1.0*np.array(f['MERRA2/t2m'])
tqv3= 1.0*np.array(f['MERRA2/tqv'])
frz3= 0.001*np.array(f['MERRA2/h273'])


clat = getattr(f,'center_lat')
clon = getattr(f,'center_lon')
cdate = getattr(f,'center_date')
satname = getattr(f,'satname')
sensor = getattr(f,'sensor')
tbnames = getattr(f,'tb_order')

f.close()

#clat = 21
#clon = -95

nx = lats.shape[0]
ny = lons.shape[1]
print (nx,ny)

cdate2 = cdate.split()[0]
ctime2 = cdate.split()[1].split(':')
cdate = satname + ' ' + sensor + ' ' + cdate2 + ' ' + ctime2[0] + ctime2[1] + ' UTC' 


#---------------------------------------------

nx = lats.shape[0]
ny = lats.shape[1]
print ("dims=",nx,ny)

msize = 5

if ( sensor == 'ATMS' ):
  tb1 = tb[:,:,1]
  tb1name = tbnames[1]
  tb2 = tb[:,:,2]
  tb2name = tbnames[2]
  tb3 = tb[:,:,3]
  tb3name = tbnames[3]
elif ( sensor == 'SSMIS' ):
  tb1 = tb[:,:,1]
  tb1name = tbnames[1]
  tb2 = tb[:,:,4]
  tb2name = tbnames[4]
  tb3 = tb[:,:,6]
  tb3name = tbnames[6]
elif ( sensor == 'AMSR2' ):
  tb1 = tb[:,:,1]
  tb1name = tbnames[1]
  tb2 = tb[:,:,7]
  tb2name = tbnames[7]
  tb3 = tb[:,:,9]
  tb3name = tbnames[9]
elif ( sensor == 'GMI' ):
  tb1 = tb[:,:,1]
  tb1name = tbnames[1]
  tb2 = tb[:,:,8]
  tb2name = tbnames[8]
  tb3 = tb[:,:,10]
  tb3name = tbnames[10]
elif ( sensor == 'MHS' ):
  tb1 = tb[:,:,0]
  tb1name = tbnames[0]
  tb2 = tb[:,:,2]
  tb2name = tbnames[2]
  tb3 = tb[:,:,4]
  tb3name = tbnames[4]
elif ( sensor == 'SAPHIR' ):
  tb1 = tb[:,:,3]
  tb1name = tbnames[3]
  tb2 = tb[:,:,4]
  tb2name = tbnames[4]
  tb3 = tb[:,:,5]
  tb3name = tbnames[5]
else:
  print ('Not recognized ', sensor)
  sys.exit(1)

tb1[tb1 < 0]= np.nan
tb2[tb2 < 0]= np.nan
tb3[tb3 < 0]= np.nan

ts1[ts1< 0]= np.nan
t2m1[t2m1< 0]= np.nan
tqv1[tqv1< 0]= np.nan
frz1[frz1< 0]= np.nan

ts2[ts2< 0]= np.nan
t2m2[t2m2< 0]= np.nan
tqv2[tqv2< 0]= np.nan
frz2[frz2< 0]= np.nan

ts3[ts3< 0]= np.nan
t2m3[t2m3< 0]= np.nan
tqv3[tqv3< 0]= np.nan
frz3[frz3< 0]= np.nan


edge_lat1 = lats[:,0]
edge_lon1 = lons[:,0]
edge_lat2 = lats[:,ny-1]
edge_lon2 = lons[:,ny-1]
edge_lat3 = lats[0,:]
edge_lon3 = lons[0,:]
edge_lat4 = lats[nx-1,:]
edge_lon4 = lons[nx-1,:]


# create figure and axes instances
fig = plt.figure(figsize=(21,15),dpi=80)

ax = fig.add_axes([0.1,0.1,0.8,0.8])
# create polar stereographic Basemap instance.
# 2400-km per side is OK for 5-min GMI
#m = Basemap(width=3600000,height=3600000,
m = Basemap(width=1400000,height=1400000,
            resolution='l',projection='stere',\
            lat_ts=clat,lat_0=clat,lon_0=clon)
#m = Basemap(width=6000000,height=4000000,
#            resolution='l',projection='stere',\
#            lat_ts=clat,lat_0=clat,lon_0=clon)
#m = Basemap(projection='cyl',llcrnrlat=27,urcrnrlat=35,\
#            llcrnrlon=-101,urcrnrlon=-93,resolution='i')


x, y = m(lons, lats)     # EPC map proj coordinates.

x1, y1 = m(edge_lon1, edge_lat1) 
x2, y2 = m(edge_lon2, edge_lat2) 
x3, y3 = m(edge_lon3, edge_lat3) 
x4, y4 = m(edge_lon4, edge_lat4) 

dlat = 2
dlon = 3


tsfc1 = 253
tsfc2 = 293
tvap1 = 0
tvap2 = 40


ax=plt.subplot(3, 4, 1)

m.drawcoastlines()
m.drawstates()
m.drawcountries()
parallels = np.arange(-90,90,dlat)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
meridians = np.arange(0,360,dlon)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)

svar = 'Ts (MERRA2)'
clevs = np.arange(tsfc1,tsfc2,1)
clevs2 = np.arange(tsfc1,tsfc2,10)
#cs = m.contourf(x,y,ts3,clevs,extend='both',cmap='jet')
cs=plt.scatter(x, y, c=ts3, s=msize, alpha=1.0, vmin=tsfc1, vmax=tsfc2, cmap='jet')
plt.plot(x1,y1,color='black',linewidth=1.5)
plt.plot(x2,y2,color='black',linewidth=1.5)
plt.plot(x3,y3,color='black',linewidth=1.5)
plt.plot(x4,y4,color='black',linewidth=1.5)
cbar = m.colorbar(cs,location='bottom',pad="8%")
cbar.set_ticks(clevs2)
plt.title(svar, fontsize=14, fontweight='bold')

print ("frame 1")


ax=plt.subplot(3, 4, 2)

m.drawcoastlines()
m.drawstates()
m.drawcountries()
parallels = np.arange(-90,90,dlat)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
meridians = np.arange(0,360,dlon)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)

svar = 'T2m (MERRA2)'
clevs = np.arange(tsfc1,tsfc2,1)
clevs2 = np.arange(tsfc1,tsfc2,10)
#cs = m.contourf(x,y,t2m3,clevs,extend='both',cmap='jet')
cs=plt.scatter(x, y, c=t2m3, s=msize, alpha=1.0, vmin=tsfc1, vmax=tsfc2, cmap='jet')
plt.plot(x1,y1,color='black',linewidth=1.5)
plt.plot(x2,y2,color='black',linewidth=1.5)
plt.plot(x3,y3,color='black',linewidth=1.5)
plt.plot(x4,y4,color='black',linewidth=1.5)
cbar = m.colorbar(cs,location='bottom',pad="8%")
cbar.set_ticks(clevs2)
plt.title(svar, fontsize=14, fontweight='bold')

print ("frame 2")


ax=plt.subplot(3, 4, 3)

m.drawcoastlines()
m.drawstates()
m.drawcountries()
parallels = np.arange(-90,90,dlat)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
meridians = np.arange(0,360,dlon)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)

svar = 'TotVap (MERRA2)'
clevs = np.arange(tvap1,tvap2,1)
clevs2 = np.arange(tvap1,tvap2,10)
#cs = m.contourf(x,y,tqv3,clevs,extend='both',cmap='jet')
cs=plt.scatter(x, y, c=tqv3, s=msize, alpha=1.0, vmin=tvap1, vmax=tvap2, cmap='jet')
plt.plot(x1,y1,color='black',linewidth=1.5)
plt.plot(x2,y2,color='black',linewidth=1.5)
plt.plot(x3,y3,color='black',linewidth=1.5)
plt.plot(x4,y4,color='black',linewidth=1.5)
cbar = m.colorbar(cs,location='bottom',pad="8%")
cbar.set_ticks(clevs2)
plt.title(svar, fontsize=14, fontweight='bold')

print ("frame 3")


ax=plt.subplot(3, 4, 4)

m.drawcoastlines()
m.drawstates()
m.drawcountries()
parallels = np.arange(-90,90,dlat)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
meridians = np.arange(0,360,dlon)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)

svar = sensor + ' ' + tb1name + ' (K)'
clevs = np.arange(160,290,1)
clevs2 = np.arange(160,290,20)
#cs = m.contourf(x,y,tb1,clevs,extend='both',cmap='jet')
cs=plt.scatter(x, y, c=tb1, s=msize, alpha=1.0, vmin=160, vmax=290, cmap='jet')
plt.plot(x1,y1,color='black',linewidth=1.5)
plt.plot(x2,y2,color='black',linewidth=1.5)
plt.plot(x3,y3,color='black',linewidth=1.5)
plt.plot(x4,y4,color='black',linewidth=1.5)

cbar = m.colorbar(cs,location='bottom',pad="8%")
cbar.set_ticks(clevs2)
plt.title(svar, fontsize=14, fontweight='bold')


print ("frame 4")




ax=plt.subplot(3, 4, 5)

m.drawcoastlines()
m.drawstates()
m.drawcountries()
parallels = np.arange(-90,90,dlat)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
meridians = np.arange(0,360,dlon)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)

svar = 'Ts (EPC-NS)'
clevs = np.arange(tsfc1,tsfc2,1)
clevs2 = np.arange(tsfc1,tsfc2,10)
#cs = m.contourf(x,y,ts1,clevs,extend='both',cmap='jet')
cs=plt.scatter(x, y, c=ts1, s=msize, alpha=1.0, vmin=tsfc1, vmax=tsfc2, cmap='jet')
plt.plot(x1,y1,color='black',linewidth=1.5)
plt.plot(x2,y2,color='black',linewidth=1.5)
plt.plot(x3,y3,color='black',linewidth=1.5)
plt.plot(x4,y4,color='black',linewidth=1.5)
cbar = m.colorbar(cs,location='bottom',pad="8%")
cbar.set_ticks(clevs2)
plt.title(svar, fontsize=14, fontweight='bold')

print ("frame 5")


ax=plt.subplot(3, 4, 6)

m.drawcoastlines()
m.drawstates()
m.drawcountries()
parallels = np.arange(-90,90,dlat)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
meridians = np.arange(0,360,dlon)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)

svar = 'T2m (EPC-NS)'
clevs = np.arange(tsfc1,tsfc2,1)
clevs2 = np.arange(tsfc1,tsfc2,10)
#cs = m.contourf(x,y,t2m1,clevs,extend='both',cmap='jet')
cs=plt.scatter(x, y, c=t2m1, s=msize, alpha=1.0, vmin=tsfc1, vmax=tsfc2, cmap='jet')
plt.plot(x1,y1,color='black',linewidth=1.5)
plt.plot(x2,y2,color='black',linewidth=1.5)
plt.plot(x3,y3,color='black',linewidth=1.5)
plt.plot(x4,y4,color='black',linewidth=1.5)
cbar = m.colorbar(cs,location='bottom',pad="8%")
cbar.set_ticks(clevs2)
plt.title(svar, fontsize=14, fontweight='bold')

print ("frame 6")


ax=plt.subplot(3, 4, 7)

m.drawcoastlines()
m.drawstates()
m.drawcountries()
parallels = np.arange(-90,90,dlat)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
meridians = np.arange(0,360,dlon)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)

svar = 'TotVap (EPC-NS)'
clevs = np.arange(tvap1,tvap2,1)
clevs2 = np.arange(tvap1,tvap2,10)
#cs = m.contourf(x,y,tqv1,clevs,extend='both',cmap='jet')
cs=plt.scatter(x, y, c=tqv1, s=msize, alpha=1.0, vmin=tvap1, vmax=tvap2, cmap='jet')
plt.plot(x1,y1,color='black',linewidth=1.5)
plt.plot(x2,y2,color='black',linewidth=1.5)
plt.plot(x3,y3,color='black',linewidth=1.5)
plt.plot(x4,y4,color='black',linewidth=1.5)
cbar = m.colorbar(cs,location='bottom',pad="8%")
cbar.set_ticks(clevs2)
plt.title(svar, fontsize=14, fontweight='bold')

print ("frame 7")


ax=plt.subplot(3, 4, 8)

m.drawcoastlines()
m.drawstates()
m.drawcountries()
parallels = np.arange(-90,90,dlat)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
meridians = np.arange(0,360,dlon)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)

svar = sensor + ' ' + tb2name + ' (K)'
clevs = np.arange(160,290,1)
clevs2 = np.arange(160,290,20)
#cs = m.contourf(x,y,tb2,clevs,extend='both',cmap='jet')
cs=plt.scatter(x, y, c=tb2, s=msize, alpha=1.0, vmin=160, vmax=290, cmap='jet')
plt.plot(x1,y1,color='black',linewidth=1.5)
plt.plot(x2,y2,color='black',linewidth=1.5)
plt.plot(x3,y3,color='black',linewidth=1.5)
plt.plot(x4,y4,color='black',linewidth=1.5)

cbar = m.colorbar(cs,location='bottom',pad="8%")
cbar.set_ticks(clevs2)
plt.title(svar, fontsize=14, fontweight='bold')

print ("frame 8")


ax=plt.subplot(3, 4, 9)

m.drawcoastlines()
m.drawstates()
m.drawcountries()
parallels = np.arange(-90,90,dlat)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
meridians = np.arange(0,360,dlon)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)

svar = 'Ts (EPC-MS)'
clevs = np.arange(tsfc1,tsfc2,1)
clevs2 = np.arange(tsfc1,tsfc2,10)
#cs = m.contourf(x,y,ts2,clevs,extend='both',cmap='jet')
cs=plt.scatter(x, y, c=ts2, s=msize, alpha=1.0, vmin=tsfc1, vmax=tsfc2, cmap='jet')
plt.plot(x1,y1,color='black',linewidth=1.5)
plt.plot(x2,y2,color='black',linewidth=1.5)
plt.plot(x3,y3,color='black',linewidth=1.5)
plt.plot(x4,y4,color='black',linewidth=1.5)
cbar = m.colorbar(cs,location='bottom',pad="8%")
cbar.set_ticks(clevs2)
plt.title(svar, fontsize=14, fontweight='bold')

print ("frame 9")


ax=plt.subplot(3, 4, 10)

m.drawcoastlines()
m.drawstates()
m.drawcountries()
parallels = np.arange(-90,90,dlat)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
meridians = np.arange(0,360,dlon)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)

svar = 'T2m (EPC-MS)'
clevs = np.arange(tsfc1,tsfc2,1)
clevs2 = np.arange(tsfc1,tsfc2,10)
#cs = m.contourf(x,y,t2m2,clevs,extend='both',cmap='jet')
cs=plt.scatter(x, y, c=t2m2, s=msize, alpha=1.0, vmin=tsfc1, vmax=tsfc2, cmap='jet')
plt.plot(x1,y1,color='black',linewidth=1.5)
plt.plot(x2,y2,color='black',linewidth=1.5)
plt.plot(x3,y3,color='black',linewidth=1.5)
plt.plot(x4,y4,color='black',linewidth=1.5)
cbar = m.colorbar(cs,location='bottom',pad="8%")
cbar.set_ticks(clevs2)
plt.title(svar, fontsize=14, fontweight='bold')

print ("frame 10")


ax=plt.subplot(3, 4, 11)

m.drawcoastlines()
m.drawstates()
m.drawcountries()
parallels = np.arange(-90,90,dlat)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
meridians = np.arange(0,360,dlon)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)

svar = 'TotVap (EPC-MS)'
clevs = np.arange(tvap1,tvap2,1)
clevs2 = np.arange(tvap1,tvap2,10)
#cs = m.contourf(x,y,tqv2,clevs,extend='both',cmap='jet')
cs=plt.scatter(x, y, c=tqv2, s=msize, alpha=1.0, vmin=tvap1, vmax=tvap2, cmap='jet')
plt.plot(x1,y1,color='black',linewidth=1.5)
plt.plot(x2,y2,color='black',linewidth=1.5)
plt.plot(x3,y3,color='black',linewidth=1.5)
plt.plot(x4,y4,color='black',linewidth=1.5)
cbar = m.colorbar(cs,location='bottom',pad="8%")
cbar.set_ticks(clevs2)
plt.title(svar, fontsize=14, fontweight='bold')

print ("frame 11")


ax=plt.subplot(3, 4, 12)

m.drawcoastlines()
m.drawstates()
m.drawcountries()
parallels = np.arange(-90,90,dlat)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
meridians = np.arange(0,360,dlon)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)

svar = sensor + ' ' + tb3name + ' (K)'
clevs = np.arange(160,290,1)
clevs2 = np.arange(160,290,20)
#cs = m.contourf(x,y,tb3,clevs,extend='both',cmap='jet')
cs=plt.scatter(x, y, c=tb3, s=msize, alpha=1.0, vmin=160, vmax=290, cmap='jet')
plt.plot(x1,y1,color='black',linewidth=1.5)
plt.plot(x2,y2,color='black',linewidth=1.5)
plt.plot(x3,y3,color='black',linewidth=1.5)
plt.plot(x4,y4,color='black',linewidth=1.5)

cbar = m.colorbar(cs,location='bottom',pad="8%")
cbar.set_ticks(clevs2)
plt.title(svar, fontsize=14, fontweight='bold')


print ("frame 12")















































plt.suptitle(cdate, y=0.93, fontsize=16, fontweight='bold')

#plt.show()
plt.savefig(outfile, dpi=80, transparent='True', bbox_inches='tight', pad_inches=0.1)
sys.exit()




