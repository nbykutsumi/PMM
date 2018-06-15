from numpy import *
import numpy as np
import myfunc.IO.GPM as GPM
from gv_fsub import *
import pickle
import myfunc.util as util
import GPMGV
import os, sys


iYM = [2014,8]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
lYM = lYM[::-1]
print lYM
thdist  = 7.5
prdName = '2A-CLIM'
gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()
dgName  = gv.ret_ddomYM2gName()
ldomain = gv.domains

minnum = 2
for domain in ldomain:
    for YM in lYM:
        Year,Mon = YM
        srcDir = '/home/utsumi/mnt/wellshare/GPMGV/MATCH.2A-CLIM/%.1fkm/%s/%04d%02d'%(thdist, domain, Year, Mon)
   
        if not os.path.exists(srcDir): continue
        ngv  = np.load(srcDir + '/p_ngv.npy')
        ngv  = ma.masked_less(ngv,minnum).compressed()
        mngv = ngv.mean()
        n    = len(ngv)
        print domain, YM, n, mngv



 
