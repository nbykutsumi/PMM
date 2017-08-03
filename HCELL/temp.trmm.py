from   numpy     import *
from   myfunc.IO import GPyM
from   datetime  import datetime, timedelta
import myfunc.util as util
import os

iDTime = datetime(2001,1,1,0)
eDTime = datetime(2001,1,1,1)

BBox   = [[-36,0],[36,120]]

#prj     = 'GPM.KuPR'
#prdLv   = 'L2'
#prdVer  = '03'
#var     = 'NS/SLV/precipRateESurface'

#prj     = 'GPM.KaPR'
#prdLv   = 'L2'
#prdVer  = '03'
##var     = 'HS/SLV/precipRateESurface'
#var     = 'MS/SLV/precipRateESurface'

#prj     = 'GPM.GMI'
#prdLv   = 'L2'
#prdVer  = '03'
#var     = 'S1/surfacePrecipitation'

#prj     = "TRMM.PR"
#prdLv   = "L2A25"
#prdVer  = "07"
#var     = "e_SurfRain"

#prj     = "TRMM.TMI"
#prdLv   = "L2A12"
#prdVer  = "07"
#var     = "surfacePrecipitation"

prj     = "TRMM.PR"
prdLv   = "L2A25"
prdVer  = "07"
#var     = "e_SurfRain"
var     = "rain"

gpm     = GPyM.GPM(prj, prdLv, prdVer)
gpmobt  = gpm(var, iDTime, eDTime, BBox)
print "*"*50
print "Done"


