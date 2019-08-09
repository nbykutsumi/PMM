import matplotlib
matplotlib.use('Agg')
from numpy import *
import numpy as np
#import metpy
import netCDF4
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

#g = 9.80665
erabaseDir = '/tank/utsumi/era5'


def interp_vertical(a3gph, a2orog, a3dat,a2surfdat, lh): 
    a3depth = a3gph - a2orog  # meter from suface
    nlayer,ny,nx = a3depth.shape
    X,Y = np.meshgrid(range(nx),range(ny))
    print 'a3dat[0] before lower filling -------------'
    print a3dat[0]
    a3out = empty([len(lh),ny,nx],float32)
    for ihintp, hintp in enumerate(lh):
        a3datTmp = a3dat.copy()
        #*** Fill under suface layer with surface variables ****
        for ilayer in range(nlayer):
            a2flag = ma.masked_less(a3depth[ilayer], 0).mask
            a3depth[ilayer][a2flag] = 0
            a3datTmp[ilayer][a2flag] = a2surfdat[a2flag]
        #*** Find lower and upper layer ***
        a2ilow = ma.masked_less(a3depth, hintp).argmin(axis=0) -1
        a2iup  = a2ilow + 1
        
        #*** Interpolation ***
        print 'a3depth[-1]'
        print a3depth[-1]
        a2mask = ma.masked_less(a3depth[-1], hintp).mask
        print 'a2mask'
        print a2mask
        a2iup[a2mask] = a2ilow[a2mask]   # temporally replace
       
        #*** lower later ***** 
        a2depthlow = a3depth[a2ilow,Y,X]-1
        a2datlow   = a3datTmp[a2ilow,Y,X]-1
        
        a2flag = ma.masked_equal(a2ilow,-1).mask
        a2depthlow[a2flag] = 0
        a2datlow[a2flag] = a2surfdat[a2flag]

        #*** upper later *****
        a2depthup  = a3depth[a2iup, Y,X]+1
        a2datup    = a3datTmp[a2iup, Y,X]+1
    
        #*** Interpolation *** 
        a2datintp = ((hintp-a2depthlow)*a2datup + (a2depthup-hintp)*a2datlow)/(a2depthup-a2depthlow)
        a3out[ihintp] = ma.masked_where(a2mask, a2datintp).filled(miss_out)
  
        print 'orog-------------------------'
        print a2orog
        print 'a3gph[0]------------------'
        print a3gph[0]
        print 'a2surfdat--------------------'
        print a2surfdat 
        print 'a3datTmp[0]---------------------'
        print a3datTmp[0]
        print ''
        print 'hintp=',hintp,'m'
        print ''
        print 'a2iup'
        print a2iup
        print 'a2ilow'
        print a2ilow
        print 'depthup'
        print a2depthup
        
        print 'depthlow'
        print a2depthlow
        print 'interp'
        print a2datintp
        print 'a3depth[0]'
        print a3depth[0]
    return a3datTmp[0]

def read_var_3d(var,Year,Mon,Day):
    srcDir = erabaseDir + '/%s/%04d%02d'%(var,Year,Mon)
    srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(var,Year,Mon,Day)
    with netCDF4.Dataset(srcPath) as np:
        a3var = np.variables[var][:]
    return a3var

   
def read_var_surf(var,Year,Mon,Day):
    srcDir = erabaseDir + '/%s/%04d%02d'%(var,Year,Mon)
    srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(var,Year,Mon,Day)
    ncvar = {'2t':'t2m', '2d':'d2m','sp':'sp'}[var]
 
    with netCDF4.Dataset(srcPath) as np:
        a2var = np.variables[ncvar][:]
    return a2var

def read_zmeter(Year,Mon,Day):
    ''' geopotential height [m] '''
    g = 9.80665
    var = 'z'
    srcDir = erabaseDir + '/%s/%04d%02d'%(var,Year,Mon)
    srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(var,Year,Mon,Day)
    with netCDF4.Dataset(srcPath) as np:
        a2var = np.variables[var][:]
    return a2var/g

def read_orogmeter():
    return np.load(erabaseDir + '/orog/orog.meter.na.npy')


##a3z = array([[8000,8500],[6000,6500],[3000,3500],[1000,1500]])
#a3z = array([[[1000,1500]],[[3000,3500]],[[6000,6500]],[[8000,8500]]]) 
#a2orog = array([[3200,3200]]) 
#a2surfdat = array([[0,0]])
var = 't'
surfvar = '2t'
Year,Mon,Day =2017,7,3

a2orog  = read_orogmeter()
a3zmeter= read_zmeter(Year,Mon,Day)[0][::-1,:,:]
a3dat   = read_var_3d(var,Year,Mon,Day)[0][::-1,:,:]
a2surfdat = read_var_surf(surfvar,Year,Mon,Day)[0]
##-- test for pressure --
#a2surfdat = a2surfdat /100.
#lp = [975,900,825,750,600,450,300,200]
#for i,p in enumerate(lp):
#    a3dat[i] = p
##-----------------------
#*** Point **************
ilat = 90
ilon = 0.0
dlat = -0.25
dlon = 0.25
lp = [975,900,825,750,600,450,300,200]


lat,lon =45, 158
y = int((lat - ilat)/dlat)
x = int((lon - ilon)/dlon)
nlayer = a3dat.shape[0]

print 'surf',a2orog[y,x], a2surfdat[y,x]
for ilayer in range(nlayer):
    print lp[ilayer],a3zmeter[ilayer][y,x],a3dat[ilayer][y,x]
#************************

'''
miss_out=-9999.
lh = [100] # [m]
print 'out'
a2dat0 = interp_vertical(a3zmeter, a2orog, a3dat, a2surfdat,lh)

lllat = -90
urlat = 90
lllon = 0.0
urlon = 360-0.25
X,Y = np.meshgrid(arange(0,359.75+0.01,0.25), arange(90,-90-0.01,-0.25))
a2fig = a2surfdat - a2dat0
print a2fig.min(), a2fig.max()
#a2fig = ma.masked_less(a2fig,-2)

fig = plt.figure(figsize=(8,5))
ax  = fig.add_axes([0.1,0.1,0.8,0.8])
M = Basemap(resolution='l', llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=ax)

im = M.pcolormesh(X,Y, a2fig, cmap='seismic',vmin=-10,vmax=10)
plt.colorbar(im)

figPath = '/home/utsumi/temp/temp.png'
plt.savefig(figPath)
plt.clf()
print figPath


'''
