from numpy import *
import matplotlib.pyplot as plt

a=fromfile("/home/utsumi/temp/gsmap.2014.04.bn",float32)
b=fromfile("/home/utsumi/temp/ra.2014.04.bn"   ,float32)


lbin   = arange(2,50,2)
llbin  = [(0,2)] + zip( lbin[:-1], lbin[1:])
print llbin

lave = []
lmid = []
for Bin in llbin:
  binmin = Bin[0]
  binmax = Bin[1]
  lmid.append((binmin+binmax)*0.5)
  atmp   = ma.masked_where( b<binmin, a)
  atmp   = ma.masked_where( binmax<=b,atmp)
  lave.append(atmp.mean())

figplot = plt.figure()
axplot  = figplot.add_axes([0.1,0.1,0.8,0.8])
axplot.plot(b,a,"o")
axplot.plot(lmid, lave, "-", linewidth=10)
axplot.plot([0,100], [0,100], "-", color="k")
axplot.set_ylim([0,50])
axplot.set_xlim([0,50])

figname ="/home/utsumi/temp/plot.png"
plt.savefig(figname)
print figname
plt.close()
