from numpy import *
from collections import deque
import myfunc.util as util
import matplotlib.pylab as plt
import myfunc.IO.CLOUDTYPE as CLOUDTYPE

iYM    = [2014,4]
eYM    = [2015,6]
#eYM    = [2014,7]
lYM    = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [11,12,1,2,3]]
ldattype = ["KuPR","GMI","IMERG","IMERG.IR","IMERG.MW","GSMaP","GSMaP.IR","GSMaP.MW","IMERG"]
#ldattype = ["GSMaP","GSMaP.IR","GSMaP.MW","IMERG"]
#ldattype = ["GSMaP.IR","GSMaP.MW","IMERG.IR","IMERG.MW"]
#ldattype = ["IMERG.IR"]

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

lbin = [0,3,6,9,12,16,20]
lBin = map(list,zip(lbin[:-1],lbin[1:])) + [[lbin[-1],999]]

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


vlim    = 40  # mm/hour

for dattype in ldattype:
  #for icl in lcltype:
  for icl in [1,3,4]:
    lpr = deque([])
    lra = deque([])
    for Year,Mon in lYM:
      baseDir  = ibaseDir
      sDir     = baseDir + "/VsRA.CL.%s"%(dattype)
      prPath   = sDir + "/%s.%04d.%02d.%s.bn"%(dattype,Year,Mon,dclShortName[icl])
      raPath   = sDir + "/RA.%04d.%02d.%s.bn"%(Year,Mon,dclShortName[icl])

      a1pr     = fromfile(prPath, float32)
      a1ra     = fromfile(raPath, float32)
      lpr.extend(a1pr)
      lra.extend(a1ra)

    lpr = array(lpr)
    lra = array(lra)

    valLow = []
    valMid = []
    valUpp = []

    botMid = []
    botUpp = []

    for Bin in lBin:
      BinMin, BinMax = Bin
      lratmp = ma.masked_where        (logical_or((lpr< BinMin),(lpr>BinMax)), lra)
      lraLow = ma.masked_greater_equal( lratmp, BinMin)
      lraMid = ma.masked_outside      ( lratmp, BinMin, BinMax)
      lraUpp = ma.masked_less_equal   ( lratmp, BinMax)

      vTot   = lratmp.count()
      vLow   = lraLow.count()/ float(vTot)
      vMid   = lraMid.count()/ float(vTot)
      vUpp   = lraUpp.count()/ float(vTot)

      valLow.append(vLow)
      valMid.append(vMid)
      valUpp.append(vUpp)

      botMid.append(vLow)
      botUpp.append((vLow+vMid))

    #** Figure **********
    figplot = plt.figure(figsize=(4.1,1.0))
    axplot  = figplot.add_axes([0.11,0.20, 0.84, 0.58])

    # bar plot ----------
    X       = range(len(lBin))
    axplot.bar(X, valLow, color="blue")
    axplot.bar(X, valMid, color="gray",bottom=botMid)
    axplot.bar(X, valUpp, color="red" ,bottom=botUpp)

    # X-ticks
    plt.xticks(X, map(str, lbin[:-1]) + ["%d<"%(lbin[-1])])

    # Add title
    stitle = "%04d/%02d-%04d/%02d %s CL=%s"%(iYM[0],iYM[1],eYM[0],eYM[1],dattype, dclShortName[icl])
    plt.title(stitle, fontsize=10)
    # Save
    figDir = baseDir + "/pict"
    figPath= figDir  + "/bar.vsRA.%s.%s.png"%(dattype,dclShortName[icl])

    util.mk_dir(figDir)
    plt.savefig(figPath)
    print figPath
    #plt.close()
