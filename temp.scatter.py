from numpy import *
import matplotlib.pyplot as plt

Y  = 2014
M  = 7

#dattype1 = "GSMaP.gauge"
dattype1 = "RA"
#dattype2 = "PR.2A25"
dattype2 = "DPR.L2.DPR"
#ny, nx   = 1200, 3600  # GSMaP
ny, nx   = 2800, 3200   # RA

#idir  = "/tank/utsumi/PMM/COMP.MAP.OBT/%s.vs.%s/%04d"%(dattype1, dattype2, Y)
idir = "/home/azariah/AH/out/COMP.MAP.OBT/%s.vs.%s/%04d"%(dattype1, dattype2, Y)

sumpath1 = idir +  "/sum.%s.%02d.%d.%d"%(dattype1, M, ny, nx)
sumpath2 = idir +  "/sum.%s.%02d.%d.%d"%(dattype2, M, ny, nx)
numpath  = idir +  "/num.%02d.%d.%d"%(M, ny, nx)

print sumpath1
print sumpath2

a2sum1  = fromfile(sumpath1, float32).reshape(ny,nx)
a2sum2  = fromfile(sumpath2, float32).reshape(ny,nx)
a2num   = fromfile(numpath,  float32).reshape(ny,nx)

a2pr1   = ma.masked_where(a2num==0.0, a2sum1) / a2num
a2pr2   = ma.masked_where(a2num==0.0, a2sum2) / a2num

a2pr1   = a2pr1 * 60*60.  # mm/sec --> mm/hour
a2pr2   = a2pr2 * 60*60.  # mm/sec --> mm/hour

#--- figure ----
figplot  = plt.figure(figsize=(5,5))
axplot   = figplot.add_axes([0.1, 0.1, 0.8, 0.8])
axplot.plot(a2pr1, a2pr2, ".", color="gray")

##-- axes ----
axplot.set_xlim(0, 100)
axplot.set_ylim(0, 100)

#-- 1x1 line ---
axplot.plot([0,600], [0,600],"-", color="k" )


#-- name -----
figdir  = "."
figname = figdir + "/scatter.%s.%s.%04d.%02d.png"%(dattype1, dattype2, Y,M)
plt.savefig(figname)


print figname
