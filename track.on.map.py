from numpy import *
import pmm_func
from mpl_toolkits.basemap import Basemap
import matplotlib
import matplotlib.pyplot as plt
#---------------------
year,mon,day,hour  = 2001,4,1,12
hlen     = 24
miss     = -9999.


a2trmm   = pmm_func.timeave_backward_saone(year,mon,day,hour,hlen, prtype="PR", relaxflag=True)

a2trmm   = pmm_func.trmm2global_one(a2trmm,miss)[30:150]
a2obt    = ma.masked_where( a2trmm==miss, ones([120,360],float32)*9999. )

a2gsmap  = pmm_func.timeave_backward_saone(year,mon,day,hour,hlen, prtype="GSMaP", relaxflag=True)
a2gsmap  = ma.masked_equal(a2gsmap, 0.0)

a2shade  = ma.masked_where(a2gsmap !=miss, ones([120,360], float32)* -9999)

#----------
a2figdat = ma.masked_equal(a2gsmap, miss)*60*60. # mm/hour
lllat    = -60.
lllon    = 0.0
urlat    = 60.
urlon    = 360
cmap     = "jet"
cobt     = mycm = matplotlib.colors.ListedColormap(array([255,215,0])/255. )
figmap   = plt.figure()
axmap    = figmap.add_axes([0.1, 0.1, 0.8, 0.8])
M        = Basemap( resolution="l", llcrnrlat=lllat, llcrnrlon=lllon, urcrnrlat=urlat, urcrnrlon=urlon, ax=axmap)
im       = M.imshow(a2figdat, origin="lower", cmap=cmap, interpolation="nearest", vmax=10.)
#--- coastline ------
M.drawcoastlines()

#--- colorbar -------
plt.colorbar(im, orientation="horizontal")

#--- orbit ----------
im       = M.imshow(a2obt, origin="lower", cmap=cobt, interpolation="nearest")

#--- title ----------
stitle = "%04d-%02d-%02d-%02d"%(year,mon,day,hour)
plt.title(stitle)
#--- save -----------
figname  = "/media/disk2/out/temp/GSMaP.png"
plt.savefig(figname)
print figname

