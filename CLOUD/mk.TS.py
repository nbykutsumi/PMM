from numpy import *
import myfunc.util as util
import matplotlib.pyplot as plt
import sys

iYM = [2014,4]
eYM = [2015,3]
#eYM = [2014,5]
lYM = util.ret_lYM(iYM, eYM)
ldattype = ["KuPR","GMI","GSMaP","GSMaP.IR","GSMaP.MW","IMERG"]
#ldattype = ["GSMaP"]

ny,nx = 261, 265

lbinPr  = range(1,9+1) + range(10,48+1,2)
ncl      = 8
lcl     = range(8)
dclName ={0:"Clear Sky",   1:"Cumulonimbus(Cb)",  2:"High Cloud",3:"Mid Cloud"
         ,4:"Cumulus(Cu)", 5:"Stratocumulus(Sc)", 6:"Fog/St"    ,7:"Cloudy"}

dclShortName={0:"no", 1:"Cb",  2:"hi",3:"md"
             ,4:"Cu", 5:"Sc",  6:"St",7:"cw"}

llbinPr = [[0, lbinPr[0]]] + zip(lbinPr[0:-1], lbinPr[1:])
#----------------------
def loadSum(Year,Mon,binPr):
  baseDir = "/tank/utsumi/PMM/WNP.261x265"
  srcDir    = baseDir + "/CL.Pr.%s"%(dattype)
  iPath     = srcDir + "/sum.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncl,ny,nx)

  if binPr == 0.0:
    return zeros([ncl, ny, nx], int32)
  else:
    return fromfile(iPath, float32).reshape(ncl, ny, nx)

def loadNum(Year,Mon,binPr):
  baseDir = "/tank/utsumi/PMM/WNP.261x265"
  srcDir    = baseDir + "/CL.Pr.%s"%(dattype)
  iPath     = srcDir + "/num.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncl,ny,nx)

  if binPr == 0.0:
    return zeros([ncl, ny, nx], int32)
  else:
    return fromfile(iPath, int32).reshape(ncl, ny, nx)

#----------------------
dPR = {dattype:empty([ncl,ny,nx],float32) for dattype in ldattype}
for dattype in ldattype:
  dPR[dattype] = array([[loadSum(Year,Mon,999)[icl].sum()/(loadNum(Year,Mon,999)[icl].sum())
                  for (Year,Mon) in lYM]
                  for icl in lcl])   # mm/h, Average precip. rate

# Figure
for icl in lcl: 
  figplot = plt.figure()
  axplot  = figplot.add_axes([0.2,0.2,0.7,0.7])

  # Colors
  colors = plt.matplotlib.cm.jet(linspace(0,1,len(ldattype)))

  # Line style
  linestyle=["--" if dattype in ["KuPR","GMI"] else "-"
             for dattype in ldattype]
  # Plot
  lx = range(len(lYM))
  [axplot.plot(lx, dPR[dattype][icl],"-", linewidth=4
               ,c=colors[idattype]
               ,linestyle=linestyle[idattype])
              for idattype, dattype in enumerate(ldattype)]
  
  # X-tick label
  dMonName= util.ret_dMonName()
  lxlabel = ["%04d %s"%(Year,dMonName[Mon]) for (Year,Mon) in lYM]
  axplot.xaxis.set_ticks(lx[::3])
  axplot.xaxis.set_ticklabels(lxlabel[::3])
  
  # Legends
  plt.legend(["%s"%(dattype) for dattype in ldattype])
  
  # Title
  plt.title("precipitation rate [mm/h] %s"%(dclName[icl]))
  
  # Save figure
  sDir  = "/tank/utsumi/PMM/WNP.261x265/pict"
  sPath = sDir + "/TS.PRate.%s.png"%(dclShortName[icl]) 
  plt.savefig(sPath) 
  print sPath
  #plt.show()
