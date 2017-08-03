from numpy import *
from collections import deque
import myfunc.util as util
import matplotlib.pylab as plt
import myfunc.IO.CLOUDTYPE as CLOUDTYPE

#TrackFlag = True
TrackFlag = False

iYM    = [2014,4]
eYM    = [2015,6]
#eYM    = [2014,4]
lYM    = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [11,12,1,2,3]]
#ldattype = ["KuPR","GMI","IMERG","IMERG.IR","IMERG.MW","GSMaP","GSMaP.IR","GSMaP.MW"]
ldattype = ["KuPR","GMI","IMERG.IR","IMERG.MW","GSMaP.IR","GSMaP.MW"]

#clVer = "JMA1"
clVer = "MyWNP.M.3"

Xvar  = "RA"
#Xvar  = "KuPR"

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

#llndsea = ["any","lnd","sea","cst"]
llndsea = ["sea"]
#*******************************
def loadData(dattype, lndsea, icl, lYM):
  lpr = deque([])
  for Year,Mon in lYM:
    baseDir  = ibaseDir
    if   TrackFlag == False:
      sDir     = baseDir + "/VsRA.CL.%s/%04d"%(dattype,Year)
    elif TrackFlag == True:
      sDir     = baseDir + "/Tr.VsRA.CL.%s/%04d"%(dattype,Year)

    prPath   = sDir + "/%s.%04d.%02d.%s.%s.bn"%(dattype,Year,Mon,lndsea,dclShortName[icl])
    a1pr     = fromfile(prPath, float32)
    lpr.extend(a1pr)

  return array(lpr)

def loadRA(dattype, lndsea, icl, lYM):
  lpr = deque([])
  for Year,Mon in lYM:
    baseDir  = ibaseDir
    if   TrackFlag == False:
      sDir     = baseDir + "/VsRA.CL.%s/%04d"%(dattype,Year)
    elif TrackFlag == True:
      sDir     = baseDir + "/Tr.VsRA.CL.%s/%04d"%(dattype,Year)

    prPath   = sDir + "/RA.%04d.%02d.%s.%s.bn"%(Year,Mon,lndsea,dclShortName[icl])
    a1pr     = fromfile(prPath, float32)
    lpr.extend(a1pr)

  return array(lpr)

def loadKu(dattype, lndsea, icl, lYM):
  lpr = deque([])
  for Year,Mon in lYM:
    baseDir  = ibaseDir
    if   TrackFlag == False:
      sDir     = baseDir + "/VsRA.CL.%s/%04d"%(dattype,Year)
    elif TrackFlag == True:
      sDir     = baseDir + "/Tr.VsRA.CL.%s/%04d"%(dattype,Year)

    prPath   = sDir + "/KuPR.%04d.%02d.%s.%s.bn"%(Year,Mon,lndsea,dclShortName[icl])
    a1pr     = fromfile(prPath, float32)
    lpr.extend(a1pr)

  return array(lpr)

if   Xvar == "RA"  : loadX = loadRA
elif Xvar == "KuPR": loadX = loadKu

#*******************************
vlim    = 40  # mm/hour
for lndsea in llndsea:
  for i, icl in enumerate(lcltype[1:] + [99]):
  #for icl in [99]:
    #** Figure **********
    figplot = plt.figure(figsize=(3,3))
    axplot1  = figplot.add_axes([0.2,0.1,0.8,0.8])
    #axplot1.set_ylim(0.0, vlim)
    #axplot1.set_xlim(0.0, vlim)

    dmean = {}
    for dattype in ldattype:
      # Load data
      if   icl != 99:
        lpr = loadData(dattype, lndsea, icl, lYM)
        lra = loadX   (dattype, lndsea, icl, lYM)
  
      elif icl ==99:
        lpr = concatenate([loadData(dattype, lndsea, icltmp, lYM) for icltmp in lcltype])
        lra = concatenate([loadX   (dattype, lndsea, icltmp, lYM) for icltmp in lcltype])

      # Average line
      #bins   = arange(0, 100,2)
      #bins   = r_[arange(0, 10+0.1,2), arange(15,100,5)]
      #BINS   = zip(bins[:-1],bins[1:])
      BINS = [[0,999]]
      dmean[dattype] = [ma.masked_where(logical_or( lra<binmin,  binmax<=lra), lpr).mean() for (binmin,binmax) in BINS]
      lx     = [mean(BIN) for BIN in BINS]


    # Draw average lines
    lines = [axplot1.plot(lx,dmean[dattype],"o", linewidth=2) 
                for dattype in ldattype] 

    if icl==1:
      print dmean

#    # Draw 1-1 line
#    axplot1.plot([0,100],[0,100],"--",color="k")
#
#
    # Add title
    stitle = "%04d/%02d-%04d/%02d CL=%s [%s]"%(iYM[0],iYM[1],eYM[0],eYM[1],dclShortName[icl], lndsea)
    plt.title(stitle, fontsize=10)
    # Save
    figDir = ibaseDir + "/pict"
    if   TrackFlag == False:
      figPath= figDir  + "/temp.lines.mulProd.vs%s.%s.%s.png"%(Xvar,lndsea,dclShortName[icl])
    elif TrackFlag == True:
      figPath= figDir  + "/temp.Tr.lines.mulProd.vs%s.%s.%s.png"%(Xvar,lndsea,dclShortName[icl])
#
#    util.mk_dir(figDir)
    plt.savefig(figPath)
    print figPath
#    #plt.show() 
#    #plt.close()
#
    # Legend file
    if i==0:
      legPath = figDir + "/temp.legend.mulProd.vs%s.png"%(Xvar)
      figleg  = plt.figure(figsize=(2,3))
      lines   = [line[0] for line in lines]   # 2D list to 1D list
      figleg.legend(lines, ldattype)
      figleg.savefig(legPath)
      plt.close()

####
####
