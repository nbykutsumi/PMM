import matplotlib
matplotlib.use("Agg")
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
#ldattype = ["KuPR","GMI","IMERG","IMERG.IR","IMERG.MW","GSMaP","GSMaP.IR","GSMaP.MW"]
#ldattype = ["GMI","IMERG","IMERG.IR","IMERG.MW","GSMaP.IR","GSMaP.MW","GSMaP"]
#ldattype = ["GMI","IMERG.IR","IMERG.MW","GSMaP.IR","GSMaP.MW"]
#ldattype = ["IMERG.IR","IMERG.MW","GSMaP.IR","GSMaP.MW"]
ldattype = ["IMERG.MW","GSMaP.IR","GSMaP.MW"]
#ldattype = ["GSMaP"]
#ldattype = ["IMERG.IR"]

#clVer = "JMA1"
clVer = "MyWNP.M.3"

Xvar  = "KuPR"
#Xvar  = "RA"

#pdfflag = True
pdfflag = False

#llndsea = ["any","lnd","sea","cst"]
#llndsea = ["sea","lnd"]
llndsea = ["sea"]
#llndsea = ["sea","lnd"]


#rootDir = "/tank/utsumi"
rootDir = "/home/utsumi/mnt/wellshare"
if clVer   == "JMA1":
  cl         = CLOUDTYPE.CloudWNP()
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.JMA"
  ibaseDirCL = "/tank/utsumi/CLOUDTYPE/WNPAC"

elif clVer[:5] == "MyWNP":
  ver        = clVer[5:]
  cl         = CLOUDTYPE.MyCloudWNP(ver=ver)
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.My%s"%(ver)
  ibaseDirCL = rootDir + "/CLOUDTYPE/MyWNP%s"%(ver)

if clVer[5:8]==".M.":
  MidFlag = True
else:
  MidFlag = False
print "*"*50+"\n"
print "MidFlag=",MidFlag
print "*"*50

lcltype = cl.licl
ncltype = len(lcltype)
dclName = cl.dclName
dclShortName = cl.dclShortName
dclid   = cl.dclid
LatUp   = cl.Lat
LonUp   = cl.Lon

#*******************************
def loadData(dattype, lndsea, icl, lYM):
  lpr = deque([])
  for Year,Mon in lYM:
    baseDir  = ibaseDir
    #sDir     = baseDir + "/VsKuPR.CL.%s/%04d"%(dattype,Year)
    sDir     = baseDir + "/Vs%s.CL.%s/%04d"%(Xvar,dattype,Year)

    prPath   = sDir + "/%s.%04d.%02d.%s.%s.bn"%(dattype,Year,Mon,lndsea,dclShortName[icl])
    a1pr     = fromfile(prPath, float32)
    lpr.extend(a1pr)

  return array(lpr)

def loadKu(dattype, lndsea, icl, lYM):
  lpr = deque([])
  for Year,Mon in lYM:
    baseDir  = ibaseDir
    sDir     = baseDir + "/VsKuPR.CL.%s/%04d"%(dattype,Year)

    prPath   = sDir + "/KuPR.%04d.%02d.%s.%s.bn"%(Year,Mon,lndsea,dclShortName[icl])
    a1pr     = fromfile(prPath, float32)
    lpr.extend(a1pr)

  return array(lpr)

def loadRA(dattype, lndsea, icl, lYM):
  lpr = deque([])
  for Year,Mon in lYM:
    baseDir  = ibaseDir
    sDir     = baseDir + "/VsRA.CL.%s/%04d"%(dattype,Year)

    prPath   = sDir + "/RA.%04d.%02d.%s.%s.bn"%(Year,Mon,lndsea,dclShortName[icl])
    a1pr     = fromfile(prPath, float32)
    lpr.extend(a1pr)

  return array(lpr)


if   Xvar == "KuPR": loadX = loadKu
elif Xvar == "RA":   loadX = loadRA

#*******************************
vmin    = 0.01
vmax    = 100  # mm/hour
#cols    = plt.matplotlib.cm.jet(linspace(0,1,len(ldattype)))
#coldarks= [col*array([1,1,1,0.5]) for col in cols]

for lndsea in llndsea:
  for dattype in ldattype:
    for i, icl in enumerate(lcltype[1:] + [99]):

      if   pdfflag != True:
        figplot = plt.figure(figsize=(4.5,4.5))
        axplot1  = figplot.add_axes([0.11,0.1,0.6,0.6])

      elif pdfflag == True:
        figplot = plt.figure(figsize=(4.5,4.5))
        axplot1  = figplot.add_axes([0.11,0.1,0.6,0.6])

      axplot1.set_xscale("log",nonposx="clip")
      axplot1.set_yscale("log",nonposy="clip")
      axplot1.set_ylim(vmin, vmax)
      axplot1.set_xlim(vmin, vmax)
    
      # Load data
      if   icl != 99:
        lpr = loadData(dattype, lndsea, icl, lYM)
        lref= loadX   (dattype, lndsea, icl, lYM)
      
      elif icl ==99:
        lpr = concatenate([loadData(dattype, lndsea, icltmp, lYM) for icltmp in lcltype])
        lref= concatenate([loadX   (dattype, lndsea, icltmp, lYM) for icltmp in lcltype])
   
      # Scatter plot
      axplot1.scatter(lref, lpr, color="gray")

      # tick labels for axplot1
      axplot1.tick_params(top="on",right="on",direction="in", labelsize=14)
    

   
      # Average line
      lk     = arange(-2,2.6+0.01,0.2)
      bins   = array([10**k for k in lk])
      BINS   = zip(bins[:-1],bins[1:])
      lmean  = [ma.masked_where(logical_or( lref<binmin,  binmax<=lref), lpr).mean() for (binmin,binmax) in BINS]
      lx     = [mean(BIN) for BIN in BINS]
      axplot1.plot(lx, lmean,"-", linewidth=2, color="k")
   
      ## Contour
      #bins   = arange(0.02,100,1)
      #H, ybins, xbins = histogram2d(lpr, lref, [bins, bins], normed=False)   # histogram2d(y,x) <-- keep this order of y and x !!!
  
      #mbins  = (bins[:-1] + bins[1:])*0.5
      #X,Y    = meshgrid(mbins, mbins)
      #levels = arange(5, 10000, 100)
      ##levels = [10**k for k in [-2,2.6+0.01,0.5]]
      ##axplot1.contour(X,Y,H,levels,color="k")
      ##axplot1.contour(X,Y,H,levels,cmap="jet")
      #axplot1.contour(X,Y,H,levels,cmap="gist_ncar")
   
      # Draw 1-1 line
      axplot1.plot([0.01,100],[0.01,100],"--",color="k")


      if pdfflag==True: 
        # PDF & CDF for precip --------
        #axplot2  = figplot.add_axes([0.77,0.1,0.16,0.6])
        axplot2  = figplot.add_axes([0.78,0.1,0.16,0.6])
        axplot2.set_xscale("log",nonposx="clip")

  
        # PDF
        bins     = arange(0,100,0.1)
        y, bins  = histogram(lpr, bins=bins, density=True)
        x        = (bins[:-1]+bins[1:])*0.5
        axplot2.plot(y, x, "-", color="k")
  
        # options
        axplot2.set_ylim(vmin,vmax)
        #for itick, tick in enumerate(axplot2.xaxis.get_major_ticks()):
        #   if itick not in [0, len(axplot2.xaxis.get_major_ticks())-2]:
        #     tick.label.set_visible(False)
 
        axplot2.tick_params(top="on",right="on",labelleft="off", direction="in", labelsize=14)
  
        # CDF
        axplot22 = axplot2.twiny()
        axplot22.set_yscale("log",nonposy="clip")
        x        = sort(lpr)
        y        = cumsum(x)/sum(x)
        axplot22.plot(y, x, "--", color="red")
 
        # options
        axplot22.set_ylim(vmin,vmax)
        axplot22.set_xlim(0.0,1.0)
        labels   = list(axplot22.get_xticks())
        nlabels  = len(labels)
        for i in range(1,nlabels-1):
          labels[i] = ""
        axplot22.set_xticklabels(labels, color="red")
        axplot22.tick_params(labelleft="off",direction="in", labelsize=14, labelcolor="r")
  
        ## PDF & CDF for Ref -----------
        #axplot3  = figplot.add_axes([0.1,0.77,0.6,0.16])
        axplot3  = figplot.add_axes([0.11,0.77,0.6,0.16])
  
        bins     = arange(0,100,0.1)
        y, bins  = histogram(lref, bins=bins, density=True)
        x        = (bins[:-1]+bins[1:])*0.5
        axplot3.plot(x, y, "-", color="k")
        # options
        axplot3.set_xlim(vmin,vmax)
        axplot3.tick_params(labelbottom="off",direction="in", labelsize=14)
  
        # CDF
        axplot32 = axplot3.twinx()
        x        = sort(lref)
        y        = cumsum(x)/sum(x)
        axplot32.plot(x,y, "--", color="red")
  
        # options
        axplot32.set_xlim(vmin,vmax)
        axplot32.set_ylim(0,1.0)
        #labels   = list(axplot32.get_yticks())
        #nlabels  = len(labels)
        #for i in range(1,nlabels-1):
        #  labels[i] = ""
        #axplot32.set_yticklabels(labels, color="red")
        axplot32.tick_params(labelbottom="off", direction="in", labelsize=14, labelcolor="r")



      # Add title
      stitle = "%04d/%02d-%04d/%02d %s CL=%s [%s]"%(iYM[0],iYM[1],eYM[0],eYM[1],dattype, dclShortName[icl], lndsea)
      plt.title(stitle, fontsize=10)
      # Save
      figDir = ibaseDir + "/pict"
      figPath= figDir  + "/scatter.log.%s.vs%s.%s.%s.png"%(dattype,Xvar,lndsea,dclShortName[icl])
    
      util.mk_dir(figDir)
      plt.savefig(figPath)
      print figPath
      #plt.show() 
      #plt.close()
 

