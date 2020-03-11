# %%
import numpy as np
from numpy import *
import myfunc.util as util
import JPLDB
import os, sys

lidx_db = range(29*29*29)
#lidx_db = lidx_db[:1000]
#lidx_db = [5810]
#lsensor = ['AMSR2','ATMS','SSMIS','MHS']
lsensor = ['ATMS','SSMIS','MHS']

#******************
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except OSError:
    pass
#******************
for sensor in lsensor:
    db = JPLDB.JPLDB(sensor)
    srcDir = '/media/disk2/share/PMM/JPLDB/EPC_DB/%s_EPC_DATABASE'%(sensor)

    for idx_db in lidx_db:
        srcPth = srcDir + '/db_%05d.bin'%(idx_db)
        if not os.path.exists(srcPth): continue
        db.set_file(srcPth)

        a1elev = db.get_var('elev')
        a2prof= db.get_var('precip_water_prof_NS')
        a1sfcbin= db.get_var('bin_sfc_ku')

        #-- Find lowest clutter free bin ---
        ny,nx = a2prof.shape
        a2idxtmp= np.array(range(nx)*ny).reshape(ny,nx).astype(int32)
        a1clutbin = ma.masked_where(a2prof<0, a2idxtmp).max(axis=1).astype(int32)  # clutter free bin (top=0, bottom=87)

        #-- Fill cltter range -----
        a1y   = np.arange(ny).astype(int32)
        a1valid = a2prof[a1y, a1clutbin]
        setclutbin = list(set(a1clutbin))

        for clutbin in setclutbin:
            a1yTmp = ma.masked_where(a1clutbin !=clutbin, a1y).compressed()
            a2prof[a1yTmp, clutbin:] = a1valid[a1yTmp].reshape(-1,1)

        #-- Make double-size array --
        a2prof = np.concatenate([np.zeros([ny,nx], int16), a2prof], axis=1)
        a1sfcbin = a1sfcbin + nx

        #-- Relative to acctual surface ---
        a2oprof = np.zeros([ny,nx], int16)

        for i in range(nx):
            a1k = a1sfcbin - nx + 1 + i  # i=0-87
            a2oprof[:,i] = a2prof[a1y,a1k]

        #-- Save ---------------------
        outDir = os.path.dirname(srcDir) + '/%s_rs_precip_water_prof_NS'%(sensor)
        outPath = outDir + '/%05d.npy'%(idx_db)
        mk_dir(outDir)
        np.save(outPath, a2oprof)
        print outPath
        #sys.exit()


        #a2proftmp1 = ma.masked_less(a2prof[:,-22:],0)
        #a2proftmp2 = ma.masked_less(a2oprof[:,-22:],0)

        #for i in range(ny):
        #    print i,'prof0',88*2-1-a1sfcbin[i], a2proftmp1[i]
        #    print i,'prof1',88*2-1-a1sfcbin[i], a2proftmp2[i]
        #    print ''
        #sys.exit()

# %%
