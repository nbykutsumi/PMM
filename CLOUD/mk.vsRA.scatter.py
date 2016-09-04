from numpy import *
from collections import deque
import myfunc.util as util
import matplotlib.pylab as plt
import myfunc.IO.CLOUDTYPE as CLOUDTYPE

iYM    = [2014,4]
eYM    = [2015,6]
#eYM    = [2014,4]
lYM    = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [11,12,1,2,3]]
ldattype = ["KuPR","GMI","IMERG","IMERG.IR","IMERG.MW","GSMaP","GSMaP.IR","GSMaP.MW"]
#ldattype = ["GSMaP","GSMaP.IR","GSMaP.MW","IMERG"]
#ldattype = ["GSMaP.IR","GSMaP.MW","IMERG.IR","IMERG.MW"]
#ldattype = ["GSMaP.IR","IMERG.IR"]

#dattype= "RA"
#dattype= "GSMaP"
#dattype= "GSMaP.MW"
#dattype= "GSMaP.IR"
#dattype= "IMERG"
#dattype= "KuPR"
#dattype= "GMI"

#clVer = "JMA1"
#clVer = "MyWNP1"
clVer = "MyWNP2"

#rootDir = "/tank/utsumi"
rootDir = "/home/utsumi/mnt/well.share"
if clVer   == "JMA1":
  cl         = CLOUDTYPE.CloudWNP()
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.JMA"
  ibaseDirCL = "/tank/utsumi/CLOUDTYPE/WNPAC"

elif clVer[:5] == "MyWNP":
  ver        = int(clVer[5:])
  cl         = CLOUDTYPE.MyCloudWNP(ver=ver)
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.My%d"%(ver)
  ibaseDirCL = rootDir + "/CLOUDTYPE/MyWNP%d"%(ver)

lcltype = cl.licl
ncltype = len(lcltype)
dclName = cl.dclName
dclShortName = cl.dclShortName
dclid   = cl.dclid
LatUp   = cl.Lat
LonUp   = cl.Lon

#*******************************
def loadData(dattype, icl, lYM):
  lpr = deque([])
  for Year,Mon in lYM:
    baseDir  = ibaseDir
    sDir     = baseDir + "/VsRA.CL.%s"%(dattype)
    prPath   = sDir + "/%s.%04d.%02d.%s.bn"%(dattype,Year,Mon,dclShortName[icl])
    a1pr     = fromfile(prPath, float32)
    lpr.extend(a1pr)

  return array(lpr)

def loadRA(dattype, icl, lYM):
  lpr = deque([])
  for Year,Mon in lYM:
    baseDir  = ibaseDir
    sDir     = baseDir + "/VsRA.CL.%s"%(dattype)
    prPath   = sDir + "/RA.%04d.%02d.%s.bn"%(Year,Mon,dclShortName[icl])
    a1pr     = fromfile(prPath, float32)
    lpr.extend(a1pr)

  return array(lpr)

#*******************************

vlim    = 40  # mm/hour

for dattype in ldattype:
  #for icl in lcltype + [99]:
  for icl in [99]:
    if   icl != 99:
      lpr = loadData(dattype, icl, lYM)
      lra = loadRA(dattype  , icl, lYM)

    elif icl ==99:
      lpr = concatenate([loadData(dattype, icltmp, lYM) for icltmp in lcltype])
      lra = concatenate([loadRA(dattype,   icltmp, lYM) for icltmp in lcltype])

    #** Figure **********
    figplot = plt.figure(figsize=(4.5,4.5))

    # Scatter plot ----------
    axplot1  = figplot.add_axes([0.1,0.1,0.6,0.6])
    axplot1.scatter(lra, lpr, color="gray")
    axplot1.plot([0,100],[0,100],"--",color="k")

    axplot1.set_ylim(0.0, vlim)
    axplot1.set_xlim(0.0, vlim)

    # Average line
    bins   = arange(0, 100,1)
    BINS   = zip(bins[:-1],bins[1:])
    lmean  = [ma.masked_where(logical_or( lra<binmin,  binmax<=lra), lpr).mean() for (binmin,binmax) in BINS]
    lx     = [mean(BIN) for BIN in BINS]
    axplot1.plot(lx,lmean,"-",color="k", linewidth=2) 

    # Contour
    bins   = arange(0, 100,1)
    H, ybins, xbins = histogram2d(lpr,lra, [bins,bins], normed=True)  # histogram2d(y,x) <-- keep this order of y and x !!!

    mbins = (bins[:-1]+bins[1:])*0.5
    X,Y   = meshgrid(mbins,mbins)

    levels= arange(0.002,0.05,0.002)
    #axplot1.contour(X,Y,H,20,color="k")
    axplot1.contour(X,Y,H,levels,color="k")

    # PDF & CDF for precip --------
    axplot2  = figplot.add_axes([0.77,0.1,0.16,0.6])

    # PDF
    bins     = arange(0,70,1)
    y, bins  = histogram(lpr, bins=bins, density=True)
    x        = (bins[:-1]+bins[1:])*0.5
    axplot2.plot(y, x, "-", color="k")

    # options
    axplot2.set_ylim(0,vlim)
    for itick, tick in enumerate(axplot2.xaxis.get_major_ticks()):
       if itick not in [0, len(axplot2.xaxis.get_major_ticks())-2]:
         tick.label.set_visible(False)

    axplot2.tick_params(labelleft="off")

    # CDF
    axplot22 = axplot2.twiny()
    x        = sort(lpr)
    y        = cumsum(x)/sum(x)
    axplot22.plot(y, x, "--", color="red")

    # options
    axplot22.set_ylim(0,vlim)
    axplot22.set_xlim(0,1.0)
    labels   = list(axplot22.get_xticks())
    nlabels  = len(labels)
    for i in range(1,nlabels-1):
      labels[i] = ""
    axplot22.set_xticklabels(labels, color="red")
    axplot22.tick_params(labelleft="off")


    ## PDF & CDF for RA -----------
    axplot3  = figplot.add_axes([0.1,0.77,0.6,0.16])

    bins     = arange(0,70,1)
    y, bins  = histogram(lra, bins=bins, density=True)
    x        = (bins[:-1]+bins[1:])*0.5
    axplot3.plot(x, y, "-", color="k")
    # options
    axplot3.set_xlim(0,vlim)
    axplot3.tick_params(labelbottom="off")

    # CDF
    axplot32 = axplot3.twinx()
    x        = sort(lra)
    y        = cumsum(x)/sum(x)
    axplot32.plot(x,y, "--", color="red")

    # options
    axplot32.set_xlim(0,vlim)
    axplot32.set_ylim(0,1.0)
    labels   = list(axplot32.get_yticks())
    nlabels  = len(labels)
    for i in range(1,nlabels-1):
      labels[i] = ""
    axplot32.set_yticklabels(labels, color="red")
    axplot32.tick_params(labelbottom="off")





    # Add title
    stitle = "%04d/%02d-%04d/%02d CL=%s %s"%(iYM[0],iYM[1],eYM[0],eYM[1],dclShortName[icl],dattype)
    plt.title(stitle, fontsize=10)
    # Save
    figDir = ibaseDir + "/pict"
    figPath= figDir  + "/scatter.%s.vsRA.%s.png"%(dattype,dclShortName[icl])

    util.mk_dir(figDir)
    plt.savefig(figPath)
    print figPath
    #plt.show() 
    #plt.close()
