import  os
import calendar
import myfunc.util as util
from   datetime import datetime
from   numpy import *
import myfunc.IO.CloudSat as CloudSat

prdLev  = "2B"
#prdName = "GEOPROF"
prdName = "CLDCLASS"
prdVer  = "P_R04"
varName = {"GEOPROF" :"Radar_Reflectivity"
#          ,"CLDCLASS":"cloud_scenario"
          ,"CLDCLASS":"Height"
          }[prdName]

cs  = CloudSat.CloudSat(prdLev, prdName, prdVer)
BBox   = [[-90,0],[90,180]]
iDTime = datetime(2014,4,1,0)
eDTime = datetime(2014,4,1,2)
csobt = cs(varName, iDTime, eDTime, BBox=BBox)
data  = csobt.data
print data

