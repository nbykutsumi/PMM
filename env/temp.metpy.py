from numpy import *
import numpy as np
#import metpy
import metpy.calc as mpcalc
from metpy.units import units
import netCDF4

##--- vertual temperature  ----
#at = array([273,273])   *units.K  # [K]
#aq = array([0.01,0.01]) *units.g/units.g # [kg/kg] specific humidity
#aw = mpcalc.mixing_ratio_from_specific_humidity(aq)  # mixing ratio
#tv = mpcalc.virtual_temperature(at,  aw)  # virtual temperature
#
##--- equivalent potential temperature ---
#at = at  # [K] temperature
#aq = aq  # [kg/kg] specific humidity
#ap = array([1000,1000]) *units.hPa  # [hPa]
#adew=mpcalc.dewpoint_from_specific_humidity(aq, at, ap)  # degC
#print adew
#atheta_e = mpcalc.equivalent_potential_temperature(ap, at, adew)
#print atheta_e



#**** ERA5 *********************
srcDir = '/home/utsumi/temp/env'
a3zmeter = np.load(srcDir + '/zmeter.npy')[::-1,:,:] * units.m
a3t = np.load(srcDir + '/t.npy')[::-1,:,:] * units.K
a3q = np.load(srcDir + '/q.npy')[::-1,:,:] * units.g/units.g

a2tsurf = np.load(srcDir + '/2t.npy') * units.K
a2dsurf = np.load(srcDir + '/2d.npy') * units.K
a2sp= np.load(srcDir + '/sp.npy') / 100. *units.hPa  # [hPa]

a1p = np.array([975,900,825,750,600,450,300,200]).reshape(-1,1,1) * units.hPa
#*******************************
nplev,ny,nx = a3zmeter.shape
a4zmeter = a3zmeter.reshape(1,nplev,ny,nx)
a4t      = a3t.reshape(1,nplev,ny,nx)
a4q      = a3q.reshape(1,nplev,ny,nx)

a3tsurf  = a2tsurf.reshape(1,ny,nx)
a3dsurf  = a2dsurf.reshape(1,ny,nx)
a3sp     = a2sp.reshape(1,ny,nx)

#*******************************

print 't'
print a3t[0,300,580]
print 'q'
print a3q[0,300,580]
print 'sp'
print a2sp[300,580]


#--- equivalent potential temperature ---
at = a3t  # [K] temperature
aq = a3q  # [kg/kg] specific humidity
ap = a1p # [hPa]
adew=mpcalc.dewpoint_from_specific_humidity(aq, at, ap)  # degC
print 'dew'
print adew[0,300,580]
atheta_e = mpcalc.equivalent_potential_temperature(ap, at, adew)
print 'theta_e'
print atheta_e[0,300,580]

#--- equivalent potential temperature (surf) ---
print ''
print 'surface ----------------'
at = a2tsurf    # [K] temperature
adew = a2dsurf  #
ap = a2sp   # [hPa]
print 'dew'
print adew[300,580]
atheta_e = mpcalc.equivalent_potential_temperature(ap, at, adew)
print 'theta_e'
print atheta_e[300,580]



