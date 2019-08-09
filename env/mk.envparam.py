from numpy import *
import numpy as np
import netCDF4

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
        #*** Fill under suface layer with surface variables ****
        a3datTmp = a3dat.copy()
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
    return a3out

def read_var_3d(var,Year,Mon,Day):
    srcDir = erabaseDir + '/%s/%04d%02d'%(var,Year,Mon)
    srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(var,Year,Mon,Day)
    with netCDF4.Dataset(srcPath) as np:
        a3var = np.variables[var][:]
    return a3var

   
def read_var_surf(var,Year,Mon,Day):
    srcDir = erabaseDir + '/%s/%04d%02d'%(var,Year,Mon)
    srcPath= srcDir + '/%s.%04d.%02d.%02d.nc'%(var,Year,Mon,Day)
    ncvar = {'2t':'t2m', '2d':'d2m', 'sp':'sp'}[var]
 
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

#***********************************************
var = 't'
surfvar = '2t'
Year,Mon,Day =2017,7,3
a2orog  = read_orogmeter()
a4zmeter= read_zmeter(Year,Mon,Day)[:,::-1,:,:]
a3zmeter= a4zmeter[0]
#a3dat   = read_var_3d(var,Year,Mon,Day)[0][::-1,:,:]
#a2surfdat = read_var_surf(surfvar,Year,Mon,Day)[0]

if var=='tv':
    a4t = read_var_3d('t',Year,Mon,Day)[0][::-1,:,:]
    a4q = read_var_3d('q',Year,Mon,Day)[0][::-1,:,:]
    a4mr= 
    a2surft = read_var_surf('2t',Year,Mon,Day)[0]


##---------

miss_out=-9999.
lh = [100] # [m]
print 'out'
a3intp = interp_vertical(a3zmeter, a2orog, a3dat, a2surfdat,lh)
print 'a3intp[0]'
print a3intp[0]
