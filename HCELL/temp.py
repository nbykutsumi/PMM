import myfunc.IO.GPyML3 as GPyML3
from myfunc.fig import Fig
from numpy import *
import myfunc.IO.GSMaP as GSMaP



sensor = "TRMM.PR"
prdName= "L3A25"
version= "07"
Year= 1999
Mon = 1
#Var = "rainMean2"
Var = "rainMean1"
GRIDTYPE = 1

miss= -9999.
bnd = [0,2,4,6,8,10,15,20]
#
#BBox   = [[0,30],[30,160]]
BBox   = False
gpm = GPyML3.L3A25(version=version,GRIDTYPE=GRIDTYPE, crd="sa", BBox=BBox)

print gpm.Lat
print gpm.Lon

a   = gpm.load_var(Year,Mon,Var)
a2in= ma.masked_less(a[0],0)

figname = "./temp0.png"
cbarname= "./cbar0.png"

print "out",len(gpm.Lat), len(gpm.Lon)

Fig.DrawMapSimple(a2in=a[0], a1lat=gpm.Lat, a1lon=gpm.Lon, bnd=bnd, figname=figname, cbarname=cbarname)
#


#
#BBox   = [[0,30],[30,160]]
##BBox   = False
#gpm = GPyML3.L3A25(version=version,GRIDTYPE=2, crd="sa", BBox=BBox)
#
#a   = gpm.load_var(Year,Mon,Var)
#a2in= ma.masked_less(a[0],0)
#
#figname = "./temp1.png"
#Fig.DrawMapSimple(a2in=a[0], a1lat=gpm.Lat, a1lon=gpm.Lon, bnd=bnd, figname=figname)
##
