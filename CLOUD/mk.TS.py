from numpy import *
from bisect import bisect, bisect_left, bisect_right
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.util as util
import matplotlib.pyplot as plt
import sys

iYM = [2014,4]
#eYM = [2014,12]
eYM = [2015,6]
lYM = util.ret_lYM(iYM, eYM)
#ldattype = ["RA","KuPR","GMI","IMERG","GSMaP","GSMaP.IR","GSMaP.MW"]
ldattype = ["KuPR","GMI","IMERG","GSMaP","GSMaP.IR","GSMaP.MW"]
#ldattype = ["GSMaP"]

cl    = CLOUDTYPE.CloudWNP()
ny,nx = cl.ny, cl.nx   #261, 265
Lat   = cl.Lat
Lon   = cl.Lon

ny,nx = 261, 265
#dommask  = "RA"
dommask  = None

#BBox    = [[-0.1, 113.875],[52.1, 180.125]]  # WN.Pac
BBox    = [[20., 118.],[48., 150.]]    # RadarAMeDAS
miss  = -9999.

lbinPr  = range(1,9+1) + range(10,48+1,2)
ncl     = 8
lcl     = range(8)
dclName ={0:"Clear Sky",   1:"Cumulonimbus(Cb)",  2:"High Cloud",3:"Mid Cloud"
         ,4:"Cumulus(Cu)", 5:"Stratocumulus(Sc)", 6:"Fog/St"    ,7:"Cloudy", 99:"All"}

dclShortName={0:"no", 1:"Cb",  2:"hi",3:"md"
             ,4:"Cu", 5:"Sc",  6:"St",7:"cw", 99:"All"}

llbinPr = [[0, lbinPr[0]]] + zip(lbinPr[0:-1], lbinPr[1:])
#----------------------
iX      = bisect_left (Lon,BBox[0][1])
eX      = bisect_right(Lon,BBox[1][1])
iY      = bisect_left (Lat,BBox[0][0])
eY      = bisect_left (Lat,BBox[1][0])
a2BBox  = ones([ny,nx],float32)*miss
a2BBox[iY:eY+1, iX:eX+1]= 1.0

if dommask =="RA":
  maskPath  = "/tank/utsumi/data/RadarAMeDAS/mask/RAmask.kubota.0.20x0.25WNP.261x265"
  a2dommask = fromfile(maskPath,float32).reshape(ny,nx)
  a2dommask = ma.masked_where(a2BBox ==miss, a2dommask)
  a3dommask = array([a2dommask for i in range(len(lcl))])
else:
  a2dommask = a2BBox
  a3dommask = array([a2dommask for i in range(len(lcl))])

#----------------------
def loadSum(Year,Mon,binPr):
  baseDir = "/tank/utsumi/PMM/WNP.261x265"
  srcDir    = baseDir + "/CL.Pr.%s"%(dattype)
  iPath     = srcDir + "/sum.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncl,ny,nx)

  if binPr == 0.0:
    a3out = zeros([ncl, ny, nx], int32)
  else:
    a3out = fromfile(iPath, float32).reshape(ncl, ny, nx)

  return ma.masked_where(a3dommask==-9999., a3out)

def loadNum(Year,Mon,binPr):
  baseDir = "/tank/utsumi/PMM/WNP.261x265"
  srcDir    = baseDir + "/CL.Pr.%s"%(dattype)
  iPath     = srcDir + "/num.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncl,ny,nx)

  if binPr == 0.0:
    a3out = zeros([ncl, ny, nx], int32)
  else:
    a3out = fromfile(iPath, int32).reshape(ncl, ny, nx)

  return ma.masked_where(a3dommask==-9999., a3out)

#----------------------
dPR    = {dattype:empty([ncl,ny,nx],float32) for dattype in ldattype}
dPRall = {dattype:empty([ny,nx],float32) for dattype in ldattype}
for dattype in ldattype:
  dPR[dattype] = array([[loadSum(Year,Mon,999)[icl].sum()/(loadNum(Year,Mon,999)[icl].sum())
                  for (Year,Mon) in lYM]
                  for icl in lcl])   # mm/h, Average precip. rate

  dPRall[dattype] =  array([loadSum(Year,Mon,999).sum()/(loadNum(Year,Mon,999).sum())
                  for (Year,Mon) in lYM])  # mm/h, Average precip. rate

# Figure
for icl in lcl+[99]: 
  figplot = plt.figure(figsize=(4.1, 3.2))
  axplot  = figplot.add_axes([0.15,0.21,0.8,0.7])

  # Colors
  colors = plt.matplotlib.cm.jet(linspace(0,1,len(ldattype)))

  # Line style
  linestyle=["--" if dattype in ["KuPR","GMI"] else "-"
             for dattype in ldattype]
  # Plot
  lx = range(len(lYM))
  if icl != 99:
    lines = [axplot.plot(lx, dPR[dattype][icl],"-", linewidth=4
                 ,c=colors[idattype]
                 ,linestyle=linestyle[idattype])
                for idattype, dattype in enumerate(ldattype)]
  else:
    lines = [axplot.plot(lx, dPRall[dattype],"-", linewidth=4
                 ,c=colors[idattype]
                 ,linestyle=linestyle[idattype])
                for idattype, dattype in enumerate(ldattype)]
 
  # X-tick label
  dMonName= util.ret_dMonName()
  lxlabel = ["%04d %s"%(Year,dMonName[Mon]) for (Year,Mon) in lYM]
  axplot.xaxis.set_ticks(lx[::3])
  axplot.xaxis.set_ticklabels(lxlabel[::3], fontsize=18, rotation=22)
 
  # Y-tick label
  for tick in axplot.yaxis.get_major_ticks():
    tick.label.set_fontsize(18)
 
  ## Legends
  #plt.legend(["%s"%(dattype) for dattype in ldattype])
  
  # Title
  plt.title("precipitation rate [mm/h] %s"%(dclName[icl]))
  
  # Save figure
  sDir  = "/tank/utsumi/PMM/WNP.261x265/pict"
  if dommask == None:
    sPath = sDir + "/TS.PRate.%s.png"%(dclShortName[icl]) 
  else:
    sPath = sDir + "/%sdom.TS.PRate.%s.png"%(dommask, dclShortName[icl]) 

  plt.savefig(sPath) 
  print sPath

  # Legend file
  legPath= sDir + "/legend.TS.png"
  figleg = plt.figure(figsize=(2,5))
  lines  = [line[0] for line in lines]  # 2D list to 1D list
  figleg.legend(lines, ldattype)
  figleg.savefig(legPath)

  plt.close()
