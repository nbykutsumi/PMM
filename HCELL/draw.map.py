import myfunc.IO.GPM as GPM
from myfunc.fig import Fig
from numpy import *


sensor = "TRMM.PR"
prdName= "L3A25"
version= "07"
Year= 1999
Mon = 7
#Var = "rainMean2"
Var = "stormHeightMean"
GRIDTYPE = 2

miss= -9999.
bnd = arange(0,9000,1000)
#
#BBox   = [[0,30],[30,160]]
BBox   = False
cmap= "Spectral_r"
gpm = GPM.L3A25(version=version,GRIDTYPE=GRIDTYPE, crd="sa", BBox=BBox)

print gpm.Lat
print gpm.Lon

a   = gpm.load_var(Year,Mon,Var)


figname = "/tank/utsumi/PMM/HCELL/pict/temp.M%02d.0.png"%(Mon)
cbarname= "/tank/utsumi/PMM/HCELL/pict/cbar0.png"
stitle  = "%04d-%02d layer:%d"%(Year,Mon,0)
Fig.DrawMapSimple(a2in=a[0], a1lat=gpm.Lat, a1lon=gpm.Lon, bnd=bnd, figname=figname, cbarname=cbarname, stitle=stitle, cmap=cmap)



figname = "/tank/utsumi/PMM/HCELL/pict/temp.M%02d.1.png"%(Mon)
cbarname= "/tank/utsumi/PMM/HCELL/pict/cbar1.png"
stitle  = "%04d-%02d layer:%d"%(Year,Mon,1)
Fig.DrawMapSimple(a2in=a[1], a1lat=gpm.Lat, a1lon=gpm.Lon, bnd=bnd, figname=figname, cbarname=cbarname, stitle=stitle, cmap=cmap)

