# %%
%matplotlib inline
from numpy import *
import numpy as np
import matplotlib.pyplot as plt
import glob
import pandas as pd
import os, sys

def ret_lpath(ltmpPath, varname):
    loutPath = []
    for tmpPath in ltmpPath:
        outDir, ifname = os.path.split(tmpPath)
        ofname = varname + '.' + '.'.join(ifname.split('.')[1:])
        outPath= outDir + '/' + ofname
        loutPath.append(outPath)
    return loutPath

lrettype = ['epc','gprof']
#lrettype = ['gprof']
lvar = ['pmw','rad']
#lvar = ['pmw']
lkey = [[rettype,var] for rettype in lrettype for var in lvar]
thpr = 1 # mm/h

dfqtile=None
for [rettype,var] in lkey:
    print rettype,var
    if rettype=='epc':
        srcDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/validprof/pair/epc.glb.relsurf01.minrec1000.maxrec10000/*/*/*'
    elif rettype=='gprof':
        srcDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/validprof/pair/gprof-shift/*/*/*'
    else:
        print 'check rettype', rettype
        sys.exit()

    searchPth= srcDir + '/prof%s.*.npy'%(var)
    lprofPth = glob.glob(searchPth)
    np.random.seed(0)
    np.random.shuffle(lprofPth)
    lprofPth = sort(lprofPth[:10])
    lprecPth = ret_lpath(lprofPth, 'prec%s'%(var))
    print lprofPth
    aprof=np.concatenate([np.load(profPth) for profPth in lprofPth], axis=0)[:,:25]
    aprec=np.concatenate([np.load(precPth) for precPth in lprecPth], axis=0)

    #-- DataFrame ----
    df = pd.DataFrame(aprof, columns=np.arange(25)*0.5)
    df = pd.concat([pd.DataFrame(aprec,columns=['prec']), df],axis=1)
    df = df[df['prec']>thpr]
    df[df < 0] = np.nan
    df.columns = pd.MultiIndex.from_product([[rettype],[var],df.columns])

    #-- Quantiles ----
    dftmp = df.loc[:,(rettype,var,[2.0, 4.5, 12.0])] .quantile([0.5,0.9,0.99,0.999,0.9999])
    if dfqtile is None:
        dfqtile = dftmp
    else:
        dfqtile = pd.concat([dfqtile,dftmp], axis=1)
print dfqtile
# %%


# %%
