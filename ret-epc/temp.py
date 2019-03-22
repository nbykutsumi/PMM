import h5py
with h5py.File('GPM_EPC_002421_20140802_0726.NS_MS.nc') as h:
    a2lat = h['latitude'][:]
    a2lon = h['longitude'][:]
    a2idx = h['db_index'][:]
    a3Tb  = h['Tb'][:]
    a3epc = h['pc_emis'][:]
    #print h.items()

ny,nx = a2lat.shape
for iy in range(ny):
    for ix in range(nx):
        lat = a2lat[iy,ix]
        lon = a2lon[iy,ix]
        if (14.997182-0.001<=lat)&(lat <=14.997182+0.001):
            if (2.0716197-0.001<=lon)&(lon <=2.0716197+0.001):
                y = iy
                x = ix
                print y,x
                break 

print 'Tb=    ', a3Tb[y,x,:]
print ''
print 'EPC=   ', a3epc[y,x,:]
print ''
print 'db_idx=', a2idx[y,x]
