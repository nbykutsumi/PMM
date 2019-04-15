import h5py
srcPath = 'GPM_EPC_002421_20140802_0726.NS_MS.nc'
with h5py.File(srcPath) as h:
    litem = h.items()
    for (sitem,memo) in litem:
        print '%s \t\t %s'%(sitem,memo)
        try:
            litem2 = h[sitem].items()
        except:
            litem2 = []

        if len(litem2)>0:
            for (sitem2,memo2) in litem2:
                print '\t%s \t %s'%(sitem2,memo2)
