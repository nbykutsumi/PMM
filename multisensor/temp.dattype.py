import numpy as np

expr = 'prf'
lsensor = ['AMSR2','SSMIS','ATMS','MHS']
#lsensor = ['AMSR2']
lvar = ['lat','lon','precrad','precepc','precgpr','profrad','profepc','profgpr','stype','t2m','inc','vfracconv']

for sensor in lsensor:
    for var in lvar:
        tankDir= '/home/utsumi/mnt/lab_tank'
        srcdir = tankDir + '/utsumi/PMM/multi/pair-prof-rs/%s.%s'%(expr,sensor)
        #srcdir = tankDir + '/utsumi/PMM/multi/pair-prof-rs/test-%s.%s'%(expr,sensor)
        srcpath= srcdir + '/%s.npy'%(var)
        adat = np.load(srcpath, allow_pickle=True)

        if var=='stype':
            outtype='int16'
        else:
            outtype='float32'
        adat = adat.astype(outtype) 
        print sensor,var, adat.dtype 
        #np.save(srcpath, adat)

