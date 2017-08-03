from numpy import *
from PMM.pmm_fsub import *
from myfunc.IO import GPyM
from datetime import datetime, timedelta
import myfunc.util as util
import os

#prj     = 'GPM.KuPR'
#prdLv   = 'L2'
#prdVer  = '03'
#var     = 'NS/SLV/precipRateESurface'

#prj     = 'GPM.KaPR'
#prdLv   = 'L2'
#prdVer  = '03'
##var     = 'HS/SLV/precipRateESurface'
#var     = 'MS/SLV/precipRateESurface'

prj     = 'GPM.GMI'
prdLv   = 'L2'
prdVer  = '03'
var     = 'S1/surfacePrecipitation'

#prj     = "TRMM.PR"
#prdLv   = "L2A25"
#prdVer  = "07"
#var     = "e_SurfRain"

#prj     = "TRMM.TMI"
#prdLv   = "L2A12"
#prdVer  = "07"
#var     = "surfacePrecipitation"

#BBox    = [[-0.1, 113.875],[52.1, 180.125]]
BBox    = [[-89., 0],[89.9, 359.9]]
ny,nx   = 261, 265
miss_in = -9999.9
miss_out= -9999.

gpm     = GPyM.GPM(prj, prdLv, prdVer)

iDTime  = datetime(2014,4,1,0)
eDTime  = datetime(2014,4,1,2)
gpmobt  = gpm(var, iDTime, eDTime, BBox=BBox)

dlat    = gpmobt.lat[1:] - gpmobt.lat[:-1]
dlon    = gpmobt.lon[1:] - gpmobt.lon[:-1]

