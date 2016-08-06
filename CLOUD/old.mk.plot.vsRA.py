from numpy import *
from myfunc.regrid import Regrid
from bisect import bisect, bisect_left, bisect_right
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.util as util
import matplotlib.pyplot as plt
import sys


iYM = [2014,4]
eYM = [2014,11]
lYM = util.ret_lYM(iYM, eYM)
#ldattype = ["KuPR","GMI","IMERG","IMERG.IR","IMERG.MW","GSMaP","GSMaP.IR","GSMaP.MW"]
ldattype = ["GSMaP","GSMaP.IR","GSMaP.MW","IMERG"]

cl    = CLOUDTYPE.CloudWNP()
ny,nx = cl.ny, cl.nx   #261, 265
Lat   = cl.Lat
Lon   = cl.Lon

#dommask  = "RA"
dommask  = None

#BBox    = [[-0.1, 113.875],[52.1, 180.125]]  # WN.Pac
BBox    = [[20., 118.],[48., 150.]]    # RadarAMeDAS
miss  = -9999.


#lbinPr  = [0.0, 0.1, 0.3, 0.5, 0.7] + range(1,9+1) + range(10,48+1,2) + [999]
lbinPr  = [0.1, 0.3, 0.5, 0.7] + range(1,9+1) + range(10,48+1,2) + [999]
ncltype = 8
lcltype = range(8)
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
  a3dommask = array([a2dommask for i in range(len(lcltype))])
else:
  a2dommask = a2BBox
  a3dommask = array([a2dommask for i in range(len(lcltype))])

#----------------------
def aveCLnum(iYM, eYM, icl):
  clName = dclShortName[icl]
  sDir   = "/tank/utsumi/CLOUDTYPE/WNPAC/num"
  lYM    = util.ret_lYM(iYM, eYM)
  a2out  = zeros([ny,nx], int32) 
  for YM in lYM:
    Year, Mon = YM
    sPath = sDir + "/num.%04d%02d.%s.%dx%d"%(Year,Mon,clName,ny,nx)
    a2out = a2out + fromfile(sPath, int32).reshape(ny,nx)
  a2out  = a2out / len(lYM)

  return ma.masked_where(a2dommask==-9999., a2out)


def loadSum(Year,Mon,binPr):
  #baseDir = "/tank/utsumi/PMM/WNP.261x265"
  baseDir = "/home/utsumi/mnt/well.share/PMM/WNP.261x265"
  srcDir    = baseDir + "/ByRA.CL.Pr.%s"%(dattype)
  iPath     = srcDir + "/sum.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncltype,ny,nx)
  if binPr == 0.0:
    a3out = zeros([ncltype,ny,nx],float32)
  else:
    a3out = fromfile(iPath, float32).reshape(ncltype, ny, nx)

  return ma.masked_where(a3dommask==-9999., a3out)

def loadNum(Year,Mon,binPr):
  #baseDir = "/tank/utsumi/PMM/WNP.261x265"
  baseDir = "/home/utsumi/mnt/well.share/PMM/WNP.261x265"
  srcDir    = baseDir + "/ByRA.CL.Pr.%s"%(dattype)
  iPath     = srcDir + "/num.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncltype,ny,nx)
  if binPr == 0.0:
    a3out = zeros([ncltype,ny,nx],int32)
  else:
    a3out = fromfile(iPath, int32).reshape(ncltype, ny, nx)

  return ma.masked_where(a3dommask==-9999., a3out)

def accDat(iYM, eYM, binPr, sumnum="num"):
  dfunc = {"sum":loadSum
          ,"num":loadNum
          }
  ddtype= {"sum":"float32"
          ,"num":"int32"
          }

  accDat = zeros([ncltype, ny, nx], ddtype[sumnum])

  lYM = util.ret_lYM(iYM,eYM)
  for (Year,Mon) in lYM:
    accDat = accDat + dfunc[sumnum](Year,Mon,binPr)
  return accDat

def ret_rat(lBin,lNum,lSum):
  llBin  = [[0,lBin[0]]] + zip(lBin[:-1],lBin[1:])
  lMid   = [mean(Bin) for Bin in llBin]
  lRat   = ma.masked_where(lNum==0.0, lSum)/lNum
  return lMid, lRat

def mk_Fig(lMid, dY, icltype):
  figplot = plt.figure(figsize=(4.1, 3.2))
  axplot  = figplot.add_axes([0.14, 0.10, 0.82,0.8])

  # Colors
  colors = plt.matplotlib.cm.jet(linspace(0,1,len(ldattype)))

  # Line style
  linestyle=["--" if dattype in ["KuPR","GMI"] else "-"
             for dattype in ldattype]

  # Line width
#  linewidth = [0 if dattype in ["GSMaP.IR","GSMaP.MW"] else 4
#                        for dattype in ldattype]

  linewidth = [4 if dattype in ["GSMaP.IR","GSMaP.MW"] else 4
                        for dattype in ldattype]




  # Plot
  lines =[axplot.plot(lMid, dY[dattype,icltype]
              ,linewidth=linewidth[idattype], c=colors[idattype]
              ,linestyle=linestyle[idattype])
         for idattype, dattype in enumerate(ldattype)]

  # Axis
  axplot.set_ylim(0.0, 50.0)
  axplot.set_xlim(0.0, 50.0)

  # 1-to-1 line
  axplot.plot([0,200],[0,200], "-", color="k")

  # Axis-tick labels
  for tick in axplot.xaxis.get_major_ticks():
    tick.label.set_fontsize(20)
  for tick in axplot.yaxis.get_major_ticks():
    tick.label.set_fontsize(17)

  # Add title
  stitle = "%04d/%02d-%04d/%02d CL=%s"\
             %(iYM[0],iYM[1],eYM[0],eYM[1],dclName[icltype])

  if dommask !=None:
    stitle = stitle + " %sdom"%(dommask)

  plt.title(stitle, fontsize=10) 

  # Save figure
  #sDir  = "/tank/utsumi/PMM/WNP.261x265/pict"
  sDir  = "/home/utsumi/mnt/well.share/PMM/WNP.261x265/pict"
  if dommask==None:
    sPath = sDir + "/rat.vsRA.%s.png"%(dclShortName[icltype]) 
  else:
    sPath = sDir + "/%sdom.rat.vsRA.%s.png"%(dommask, dclShortName[icltype]) 
  plt.savefig(sPath)
  print sPath

  plt.close()
  # Legend file
  legPath = sDir + "/legend.rat.vsRA.png"
  figleg  = plt.figure(figsize=(2,3))
  lines = [line[0] for line in lines]  # 2D list to 1D list
  figleg.legend(lines, ldattype)
  figleg.savefig(legPath)
 
  plt.close()

#----------------------
dNum = {(dattype,icltype):[] 
       for dattype in ldattype
       for icltype in lcltype + [99]
       }

dSum = {(dattype,icltype):[] 
       for dattype in ldattype
       for icltype in lcltype + [99]
       }

dRat = {}

for dattype in ldattype:
  print dattype
  for (binmin, binmax) in llbinPr:
    a3num0 = accDat(iYM, eYM, binmin, sumnum="num")
    a3num1 = accDat(iYM, eYM, binmax, sumnum="num")

    a3sum0 = accDat(iYM, eYM, binmin, sumnum="sum")
    a3sum1 = accDat(iYM, eYM, binmax, sumnum="sum")

    a3num  = a3num1 - a3num0
    a3sum  = a3sum1 - a3sum0

    for icltype in lcltype:
      dNum[dattype,icltype].append(a3num[icltype].sum())
      dSum[dattype,icltype].append(a3sum[icltype].sum())

    
    dNum[dattype,99].append(a3num.sum())
    dSum[dattype,99].append(a3sum.sum())

  
  for icltype in lcltype + [99]:
    lMid, dRat[dattype,icltype] = ret_rat(lbinPr, dNum[dattype,icltype], dSum[dattype,icltype])

#-----------
# Figure 
for icltype in lcltype + [99]:
  mk_Fig(lMid, dRat, icltype)
