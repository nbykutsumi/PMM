from numpy import *
from datetime import datetime, timedelta
import TagPMM
import ioPMM
from cf2.io import GPM

prj     = 'GPM.KuPR'
prdLv   = 'L2'
prdVer  = '03'
var     = 'NS/SLV/precipRateESurface'

BBox    = [[20.,118.],[48.,150.]]  # RadarAMeDAS
res     = 0.1

iDTime  = datetime(2014,12,1,0)
eDTime  = datetime(2014,12,3,0)
gpm     = GPM.GPM(prj, prdLv, prdVer)
gpmdom  = gpm(var, iDTime, eDTime, BBox, res)
