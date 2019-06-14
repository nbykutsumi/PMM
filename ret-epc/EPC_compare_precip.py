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

outfile = "precip_compare.png"
if ( len(sys.argv) == 4 ):
  outfile = sys.argv[3]

beam = int(sys.argv[2])
#beam = 60

ncfile = sys.argv[1]
f = NetCDFFile(ncfile)
for attr in f.ncattrs():
  print (attr, '=', getattr(f, attr))

lats = 1.0*np.array(f['latitude'])
lons = 1.0*np.array(f['longitude'])
emis= 1.0*np.array(f['emis'])
tb= 1.0*np.array(f['Tb'])

precip1= 1.0*np.array(f['NS/precip'])
precip2= 1.0*np.array(f['NS/precip2'])
precip3= 1.0*np.array(f['MS/precip'])
precip4= 1.0*np.array(f['MS/precip2'])
q_precip1= np.array(f['NS/quality'])
q_precip3= np.array(f['MS/quality'])
precipg= 1.0*np.array(f['GPROF/precip'])
#sfc= np.array(f['GPROF/sfc_class'])

clat = getattr(f,'center_lat')
clon = getattr(f,'center_lon')
cdate = getattr(f,'center_date')
satname = getattr(f,'satname')
sensor = getattr(f,'sensor')
tbnames = getattr(f,'tb_order')
dprfile = getattr(f,'DPR')
cmbfile = getattr(f,'CMB')
gproffile = getattr(f,'GPROF')

print ('TB= ', tbnames)
print ('DPR= ',dprfile)
print ('CMB= ',cmbfile)

f.close()

#clat = -40
#clon = -63
#clat = 14
#clon = 3

nx = lats.shape[0]
ny = lons.shape[1]
print (nx,ny)

imid= int((nx/2)-1)
jmid= int(ny/2)
x1 = lats[imid,jmid]
x2 = lats[imid,jmid]
orb = 'Descending'
if ( x1 < x2 ): orb = 'Ascending' 
print (x1, x2, orb)


cdate2 = cdate.split()[0]
ctime2 = cdate.split()[1].split(':')
cdate = satname + ' ' + sensor + ' ' + cdate2 + ' ' + ctime2[0] + ctime2[1] + ' UTC  ' + orb 

f = h5py.File(cmbfile,'r')
lats_cmb_NS = 1.0*np.array(f['NS/Latitude'])
lons_cmb_NS = 1.0*np.array(f['NS/Longitude'])
precip_cmb_NS = 1.0*np.array(f['NS/surfPrecipTotRate'])
lats_cmb_MS = 1.0*np.array(f['MS/Latitude'])
lons_cmb_MS = 1.0*np.array(f['MS/Longitude'])
precip_cmb_MS = 1.0*np.array(f['MS/surfPrecipTotRate'])
f.close()

f = h5py.File(dprfile,'r')
lats_dpr_NS = 1.0*np.array(f['NS/Latitude'])
lons_dpr_NS = 1.0*np.array(f['NS/Longitude'])
precip_dpr_NS = 1.0*np.array(f['NS/SLV/precipRateESurface'])
lats_dpr_MS = 1.0*np.array(f['MS/Latitude'])
lons_dpr_MS = 1.0*np.array(f['MS/Longitude'])
precip_dpr_MS = 1.0*np.array(f['MS/SLV/precipRateESurface'])
f.close()

f = h5py.File(gproffile,'r')
lats_gprof = 1.0*np.array(f['S1/Latitude'])
lons_gprof = 1.0*np.array(f['S1/Longitude'])
precip_gprof = 1.0*np.array(f['S1/surfacePrecipitation'])
#status_gprof = np.array(f['S1/pixelStatus'])
status_gprof = np.array(f['S1/qualityFlag'])
f.close()

gnx = lats_gprof.shape[0]
gny = lats_gprof.shape[1]
qlats_gprof, qlons_gprof = [], []
for i in range(0,gnx):
  for j in range(0,gny):
    if ( status_gprof[i,j] != 0 ):
      qlats_gprof.append(lats_gprof[i,j])
      qlons_gprof.append(lons_gprof[i,j])

#---------------------------------------------


nx = tb.shape[0]
ny = tb.shape[1]
nc = tb.shape[2]
print ("TB dims=",nx,ny,nc)

for i in range(0,nx):
  for j in range(0,ny):
    for k in range(0,nc):
      if ( tb[i,j,k] < -20 and tb[i,j,k] > -400 ):
        tb[i,j,k]*= -1.0
tb[tb < 0]= np.nan
 
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
  tb1 = tb[:,:,3]
  tb1name = tbnames[3]
  tb2 = tb[:,:,7]
  tb2name = tbnames[7]
  tb3 = tb[:,:,9]
  tb3name = tbnames[9]
elif ( sensor == 'GMI' ):
  tb1 = tb[:,:,3]
  tb1name = tbnames[3]
  tb2 = tb[:,:,8]
  tb2name = tbnames[8]
  tb3 = tb[:,:,10]
  tb3name = tbnames[10]
elif ( sensor == 'MHS' ):
  tb1 = tb[:,:,0]
  tb1name = tbnames[0]
  tb2 = tb[:,:,1]
  tb2name = tbnames[1]
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

precip1[precip1 < 0.3]= np.nan
precip2[precip2 < 0.3]= np.nan
precip3[precip3 < 0.3]= np.nan
precip4[precip4 < 0.3]= np.nan
precipg[precipg < 0.3]= np.nan
precip_cmb_NS[precip_cmb_NS < 0.3]= np.nan
precip_cmb_MS[precip_cmb_MS < 0.3]= np.nan
precip_dpr_NS[precip_dpr_NS < 0.3]= np.nan
precip_dpr_MS[precip_dpr_MS < 0.3]= np.nan

precip_gprof[precip_gprof < 0.3]= np.nan

qlats_NS, qlons_NS = [], []
qlats_MS, qlons_MS = [], []
for i in range(0,nx):
  for j in range(0,ny):
    if ( q_precip1[i,j] == -2 ):
      qlats_NS.append(lats[i,j])
      qlons_NS.append(lons[i,j])
      print (lats[i,j], lons[i,j])
    if ( q_precip3[i,j] == -2 ):
      qlats_MS.append(lats[i,j])
      qlons_MS.append(lons[i,j])
      print (lats[i,j], lons[i,j])

beam_lat1 = lats[:,beam]
beam_lon1 = lons[:,beam]

edge_lat1 = lats[:,0]
edge_lon1 = lons[:,0]
edge_lat2 = lats[:,ny-1]
edge_lon2 = lons[:,ny-1]
edge_lat3 = lats[0,:]
edge_lon3 = lons[0,:]
edge_lat4 = lats[nx-1,:]
edge_lon4 = lons[nx-1,:]

nx2 = lats_cmb_NS.shape[0]
ny2 = lats_cmb_NS.shape[1]
print ("DPR NS dims=",nx2,ny2)
edge2_lat1 = lats_cmb_NS[:,0]
edge2_lon1 = lons_cmb_NS[:,0]
edge2_lat2 = lats_cmb_NS[:,ny2-1]
edge2_lon2 = lons_cmb_NS[:,ny2-1]

nx3 = lats_cmb_MS.shape[0]
ny3 = lats_cmb_MS.shape[1]
print ("DPR MS dims=",nx3,ny3)
edge3_lat1 = lats_cmb_MS[:,0]
edge3_lon1 = lons_cmb_MS[:,0]
edge3_lat2 = lats_cmb_MS[:,ny3-1]
edge3_lon2 = lons_cmb_MS[:,ny3-1]


# create figure and axes instances
fig = plt.figure(figsize=(20,16),dpi=80)

ax = fig.add_axes([0.1,0.1,0.8,0.8])
# create polar stereographic Basemap instance.
# 2400-km per side is OK for 5-min GMI
#m = Basemap(width=1000000,height=1000000,
m = Basemap(width=1600000,height=1600000,
            resolution='l',projection='stere',\
            lat_ts=clat,lat_0=clat,lon_0=clon)
#m = Basemap(width=6000000,height=4000000,
#            resolution='l',projection='stere',\
#            lat_ts=clat,lat_0=clat,lon_0=clon)
#m = Basemap(projection='cyl',llcrnrlat=27,urcrnrlat=35,\
#            llcrnrlon=-101,urcrnrlon=-93,resolution='i')


x, y = m(lons, lats)     # EPC map proj coordinates.
qx_NS, qy_NS = m(qlons_NS, qlats_NS)     # EPC map proj coordinates where DB was not located
qx_MS, qy_MS = m(qlons_MS, qlats_MS) 

gx, gy = m(lons_gprof, lats_gprof)     # EPC map proj coordinates.
qx_gprof, qy_gprof = m(qlons_gprof, qlats_gprof) 

xn, yn = m(lons_cmb_NS, lats_cmb_NS) # NS map proj coordinates.
xm, ym = m(lons_cmb_MS, lats_cmb_MS) # MS map proj coordinates.

bx1, by1 = m(beam_lon1, beam_lat1) 

x1, y1 = m(edge_lon1, edge_lat1) 
x2, y2 = m(edge_lon2, edge_lat2) 
x3, y3 = m(edge_lon3, edge_lat3) 
x4, y4 = m(edge_lon4, edge_lat4) 

dx1, dy1 = m(edge2_lon1, edge2_lat1) 
dx2, dy2 = m(edge2_lon2, edge2_lat2) 

dx3, dy3 = m(edge3_lon1, edge3_lat1) 
dx4, dy4 = m(edge3_lon2, edge3_lat2) 

dlat = 2
dlon = 3

msize = 4
msize2 = 1
precipmax = 12

ax=plt.subplot(3, 4, 1)

m.drawcoastlines()
m.drawstates()
m.drawcountries()
parallels = np.arange(-90,90,dlat)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=12)
meridians = np.arange(0,360,dlon)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=12)

svar = 'EPC Precip (CMB-NS)'
clevs = np.arange(0,20.0,0.2)
clevs2 = np.arange(0,20.0,5.0)
#clevs = np.arange(0,15,1)
#clevs2 = np.arange(0,15,1)
#levels = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]
#cs = m.contourf(x,y,precip1,clevs,extend='both',cmap='jet')
cs=plt.scatter(x, y, c=precip1, s=msize, alpha=1.0, vmin=0, vmax=precipmax, cmap='jet')
plt.scatter(qx_NS,qy_NS, color='magenta', s=msize, marker="s")
plt.plot(x1,y1,color='black',linewidth=1.0)
plt.plot(x2,y2,color='black',linewidth=1.0)
plt.plot(x3,y3,color='black',linewidth=1.0)
plt.plot(x4,y4,color='black',linewidth=1.0)
plt.plot(dx1,dy1,color='black',linewidth=1.0)
plt.plot(dx2,dy2,color='black',linewidth=1.0)

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

svar = 'EPC Precip (DPR-NS)'
clevs = np.arange(0,20.0,0.2)
clevs2 = np.arange(0,20.0,5.0)
#cs = m.contourf(x,y,precip2,clevs,extend='both',cmap='jet')
cs=plt.scatter(x, y, c=precip2, s=msize, alpha=1.0, vmin=0, vmax=precipmax, cmap='jet')
plt.scatter(qx_NS,qy_NS, color='magenta', s=msize, marker="s")
plt.plot(x1,y1,color='black',linewidth=1.0)
plt.plot(x2,y2,color='black',linewidth=1.0)
plt.plot(x3,y3,color='black',linewidth=1.0)
plt.plot(x4,y4,color='black',linewidth=1.0)
plt.plot(dx1,dy1,color='black',linewidth=1.0)
plt.plot(dx2,dy2,color='black',linewidth=1.0)

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

svar = 'EPC Precip (CMB-MS)'
clevs = np.arange(0,20.0,0.2)
clevs2 = np.arange(0,20.0,5.0)
#cs = m.contourf(x,y,precip3,clevs,extend='both',cmap='jet')
cs=plt.scatter(x, y, c=precip3, s=msize, alpha=1.0, vmin=0, vmax=precipmax, cmap='jet')
plt.scatter(qx_MS,qy_MS, color='magenta', s=msize, marker="s")
plt.plot(x1,y1,color='black',linewidth=1.0)
plt.plot(x2,y2,color='black',linewidth=1.0)
plt.plot(x3,y3,color='black',linewidth=1.0)
plt.plot(x4,y4,color='black',linewidth=1.0)
plt.plot(dx3,dy3,color='black',linewidth=1.0)
plt.plot(dx4,dy4,color='black',linewidth=1.0)

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

svar = 'EPC Precip (DPR-MS)'
clevs = np.arange(0,20.0,0.2)
clevs2 = np.arange(0,20.0,5.0)
#cs = m.contourf(x,y,precip4,clevs,extend='both',cmap='jet')
cs=plt.scatter(x, y, c=precip4, s=msize, alpha=1.0, vmin=0, vmax=precipmax, cmap='jet')
plt.scatter(qx_MS,qy_MS, color='magenta', s=msize, marker="s")
plt.plot(x1,y1,color='black',linewidth=1.0)
plt.plot(x2,y2,color='black',linewidth=1.0)
plt.plot(x3,y3,color='black',linewidth=1.0)
plt.plot(x4,y4,color='black',linewidth=1.0)
plt.plot(dx3,dy3,color='black',linewidth=1.0)
plt.plot(dx4,dy4,color='black',linewidth=1.0)

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

svar = 'Precip (CMB-NS)'
clevs = np.arange(0,20.0,0.2)
clevs2 = np.arange(0,20.0,5.0)
#cs = m.contourf(xn,yn,precip_cmb_NS,clevs,extend='both',cmap='jet')
cs=plt.scatter(xn, yn, c=precip_cmb_NS, s=msize2, alpha=1.0, vmin=0, vmax=precipmax, cmap='jet')
plt.plot(x1,y1,color='black',linewidth=1.0)
plt.plot(x2,y2,color='black',linewidth=1.0)
plt.plot(x3,y3,color='black',linewidth=1.0)
plt.plot(x4,y4,color='black',linewidth=1.0)
plt.plot(dx1,dy1,color='black',linewidth=1.0)
plt.plot(dx2,dy2,color='black',linewidth=1.0)

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

svar = 'Precip (DPR-NS)'
clevs = np.arange(0,20.0,0.2)
clevs2 = np.arange(0,20.0,5.0)
#cs = m.contourf(xn,yn,precip_dpr_NS,clevs,extend='both',cmap='jet')
cs=plt.scatter(xn, yn, c=precip_dpr_NS, s=msize2, alpha=1.0, vmin=0, vmax=precipmax, cmap='jet')
plt.plot(x1,y1,color='black',linewidth=1.0)
plt.plot(x2,y2,color='black',linewidth=1.0)
plt.plot(x3,y3,color='black',linewidth=1.0)
plt.plot(x4,y4,color='black',linewidth=1.0)
plt.plot(dx1,dy1,color='black',linewidth=1.0)
plt.plot(dx2,dy2,color='black',linewidth=1.0)

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

svar = 'Precip (CMB-MS)'
clevs = np.arange(0,20.0,0.2)
clevs2 = np.arange(0,20.0,5.0)
#cs = m.contourf(xm,ym,precip_cmb_MS,clevs,extend='both',cmap='jet')
cs=plt.scatter(xm, ym, c=precip_cmb_MS, s=msize2, alpha=1.0, vmin=0, vmax=precipmax, cmap='jet')
plt.plot(x1,y1,color='black',linewidth=1.0)
plt.plot(x2,y2,color='black',linewidth=1.0)
plt.plot(x3,y3,color='black',linewidth=1.0)
plt.plot(x4,y4,color='black',linewidth=1.0)
plt.plot(dx3,dy3,color='black',linewidth=1.0)
plt.plot(dx4,dy4,color='black',linewidth=1.0)

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

svar = 'Precip (DPR-MS)'
clevs = np.arange(0,20.0,0.2)
clevs2 = np.arange(0,20.0,5.0)
#cs = m.contourf(xm,ym,precip_dpr_MS,clevs,extend='both',cmap='jet')
cs=plt.scatter(xm, ym, c=precip_dpr_MS, s=msize2, alpha=1.0, vmin=0, vmax=precipmax, cmap='jet')
plt.plot(x1,y1,color='black',linewidth=1.0)
plt.plot(x2,y2,color='black',linewidth=1.0)
plt.plot(x3,y3,color='black',linewidth=1.0)
plt.plot(x4,y4,color='black',linewidth=1.0)
plt.plot(dx3,dy3,color='black',linewidth=1.0)
plt.plot(dx4,dy4,color='black',linewidth=1.0)

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

msize = 4

svar = 'Precip (GPROF-V5)'
clevs = np.arange(0,20.0,0.2)
clevs2 = np.arange(0,20.0,5.0)
#cs = m.contourf(x,y,precipg,clevs,extend='both',cmap='jet')
#cs=plt.scatter(x, y, c=precipg, s=msize, alpha=1.0, vmin=0, vmax=precipmax, cmap='jet')
cs=plt.scatter(gx, gy, c=precip_gprof, s=msize, alpha=1.0, vmin=0, vmax=precipmax, cmap='jet')
plt.scatter(qx_gprof,qy_gprof, color='magenta', s=msize, marker="s")
#plt.plot(x1,y1,color='black',linewidth=1.0)
#plt.plot(x2,y2,color='black',linewidth=1.0)
#plt.plot(x3,y3,color='black',linewidth=1.0)
#plt.plot(x4,y4,color='black',linewidth=1.0)
plt.plot(dx1,dy1,color='black',linewidth=1.0)
plt.plot(dx2,dy2,color='black',linewidth=1.0)
plt.plot(dx3,dy3,color='black',linewidth=1.0)
plt.plot(dx4,dy4,color='black',linewidth=1.0)

cbar = m.colorbar(cs,location='bottom',pad="8%")
cbar.set_ticks(clevs2)
plt.title(svar, fontsize=14, fontweight='bold')

print ("frame 9")

msize = 4

ax=plt.subplot(3, 4, 10)

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
plt.plot(bx1,by1,color='magenta',linewidth=1.5)
plt.plot(x1,y1,color='black',linewidth=1.0)
plt.plot(x2,y2,color='black',linewidth=1.0)
plt.plot(x3,y3,color='black',linewidth=1.0)
plt.plot(x4,y4,color='black',linewidth=1.0)
plt.plot(dx1,dy1,color='black',linewidth=1.0)
plt.plot(dx2,dy2,color='black',linewidth=1.0)
plt.plot(dx3,dy3,color='black',linewidth=1.0)
plt.plot(dx4,dy4,color='black',linewidth=1.0)

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

svar = sensor + ' ' + tb2name + ' (K)'
clevs = np.arange(160,290,1)
clevs2 = np.arange(160,290,20)
#cs = m.contourf(x,y,tb2,clevs,extend='both',cmap='jet')
cs=plt.scatter(x, y, c=tb2, s=msize, alpha=1.0, vmin=160, vmax=290, cmap='jet')
plt.plot(bx1,by1,color='magenta',linewidth=1.5)
plt.plot(x1,y1,color='black',linewidth=1.0)
plt.plot(x2,y2,color='black',linewidth=1.0)
plt.plot(x3,y3,color='black',linewidth=1.0)
plt.plot(x4,y4,color='black',linewidth=1.0)
plt.plot(dx1,dy1,color='black',linewidth=1.0)
plt.plot(dx2,dy2,color='black',linewidth=1.0)
plt.plot(dx3,dy3,color='black',linewidth=1.0)
plt.plot(dx4,dy4,color='black',linewidth=1.0)

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
plt.plot(bx1,by1,color='magenta',linewidth=1.5)
plt.plot(x1,y1,color='black',linewidth=1.0)
plt.plot(x2,y2,color='black',linewidth=1.0)
plt.plot(x3,y3,color='black',linewidth=1.0)
plt.plot(x4,y4,color='black',linewidth=1.0)
plt.plot(dx1,dy1,color='black',linewidth=1.0)
plt.plot(dx2,dy2,color='black',linewidth=1.0)
plt.plot(dx3,dy3,color='black',linewidth=1.0)
plt.plot(dx4,dy4,color='black',linewidth=1.0)

cbar = m.colorbar(cs,location='bottom',pad="8%")
cbar.set_ticks(clevs2)
plt.title(svar, fontsize=14, fontweight='bold')

print ("frame 12")







plt.suptitle(cdate, y=0.93, fontsize=16, fontweight='bold')
#plt.subplots_adjust(top=0.85)



#plt.show()
plt.savefig(outfile, dpi=80, transparent='True', bbox_inches='tight', pad_inches=0.1)
sys.exit()


