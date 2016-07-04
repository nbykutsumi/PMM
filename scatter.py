from numpy import *
import matplotlib.pyplot as plt

Y  = 2014
M  = 3

idir1  = "/home/utsumi/members/azariah/out/RA"
idir2  = "/home/utsumi/members/azariah/out/GSMaP.JPN"

#ipath1 = idir1 + "/%04d" + "/RadarAmedas.201401.280x320"
#ipath2 = idir2 + "/%04d" + "/gsmap_mvk.20140300.0000.v6.0000.0.dat"

ipath1 = idir1 + "/%04d"%(Y) + "/RadarAmedas.%04d%02d.280x320"%(Y,M)
ipath2 = idir2 + "/%04d"%(Y) + "/gsmap_mvk.%04d%02d00.0000.v6.0000.0.dat"%(Y,M)

print ipath1
print ipath2


a1  = fromfile(ipath1, float32).reshape(280,320)
a2  = fromfile(ipath2, float32).reshape(280,320)

#--- restrict region ----
#  20N-35N
a1  = a1[130:]
a2  = a2[130:]    

#--- mask where RA==0 ----
a1  = ma.masked_equal( a1, 0)
a2  = ma.masked_where( a1.mask, a2)
#------------------------

print a1
print a2

#--- figure ----
figplot  = plt.figure(figsize=(5,5))
axplot   = figplot.add_axes([0.1, 0.1, 0.8, 0.8])
axplot.plot(a1, a2, ".", color="gray")

#-- axes ----
axplot.set_xlim(0, 600)
axplot.set_ylim(0, 600)

#-- 1x1 line ---
axplot.plot([0,600], [0,600],"-", color="k" )


#-- name -----
figname = "./scatter.RA.GS.%04d.%02d.png"%(Y,M)
plt.savefig(figname)


print figname
