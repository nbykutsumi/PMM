from numpy import *
from myfunc.regrid import Regrid
from bisect import bisect, bisect_left, bisect_right
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.util as util
import matplotlib.pyplot as plt
import sys

#clVer = "JMA1"
#clVer = "MyWNP1"
#clVer = "MyWNP2"
#clVer = "MyWNP3"
clVer = "MyWNP.M.3"

iYM = [2014,4]
eYM = [2015,6]
#iYM = [2014,4]
#eYM = [2014,5]

lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [11,12,1,2,3]]

ldattype = ["RA","KuPR","GMI","IMERG","IMERG.IR","IMERG.MW","GSMaP","GSMaP.IR","GSMaP.MW"]
#ldattype = ["RA","KuPR"]
#ldattype = ["RA"]

# 
#rootDir = "/tank/utsumi"
rootDir = "/home/utsumi/mnt/well.share"
if clVer   == "JMA1":
  cl         = CLOUDTYPE.CloudWNP()
  ncltype = 8
  lcltype = range(ncltype)
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.JMA"
  ibaseDirCL = "/tank/utsumi/CLOUDTYPE/WNPAC"

elif clVer[:5] == "MyWNP":
  ver     = clVer[5:]
  cl      = CLOUDTYPE.MyCloudWNP(ver=ver)
  ncltype = cl.ncl
  lcltype = cl.licl
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.My%s"%(ver)
  ibaseDirCL = rootDir + "/CLOUDTYPE/MyWNP%s"%(ver)


ny,nx = cl.ny, cl.nx   #261, 265
Lat   = cl.Lat
Lon   = cl.Lon

#dommask  = "RA"
dommask  = None

#BBox    = [[-0.1, 113.875],[52.1, 180.125]]  # WN.Pac
BBox    = [[20., 118.],[48., 150.]]    # RadarAMeDAS
miss  = -9999.


#lbinPr  = [0.0, 0.1, 0.3, 0.5, 0.7] + range(1,9+1) + range(10,48+1,2) + [999]
#lbinPr  = [0.1, 0.3, 0.5, 0.7] + range(1,9+1) + range(10,48+1,2) + [999]
#lbinPr  = [0.1, 0.3, 0.5, 0.7] + range(1,9+1) + range(10,48+1,2) + [999]
lbinPr  = [0.0, 0.1, 0.3, 0.5, 0.7] + range(1,9+1,1) + range(10,18+1,2) + range(20, 40+1,4) +[999]  # Maxs of bin range
dclName = cl.dclName
dclShortName = cl.dclShortName

llbinPr = [[0, lbinPr[0]]] + zip(lbinPr[0:-1], lbinPr[1:])

llndsea = ["any","lnd","sea","cst"]
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
  #sDir   = "/tank/utsumi/CLOUDTYPE/WNPAC/num"
  sDir   = ibaseDirCL + "/num"
  lYM    = util.ret_lYM(iYM, eYM)
  a2out  = zeros([ny,nx], int32) 
  for YM in lYM:
    Year, Mon = YM
    sPath = sDir + "/num.%04d%02d.%s.%dx%d"%(Year,Mon,clName,ny,nx)
    a2out = a2out + fromfile(sPath, int32).reshape(ny,nx)
  a2out  = a2out / len(lYM)

  return ma.masked_where(a2dommask==-9999., a2out)


def loadSum(Year,Mon,binPr):
  srcDir    = ibaseDir + "/CL.Pr.%s/%04d"%(dattype,Year)
  iPath     = srcDir + "/sum.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncltype,ny,nx)
  if binPr == 0.0:
    a3out = zeros([ncltype,ny,nx],float32)
  else:
    a3out = fromfile(iPath, float32).reshape(ncltype, ny, nx)

  return ma.masked_where(a3dommask==-9999., a3out)

def loadNum(Year,Mon,binPr):
  srcDir    = ibaseDir + "/CL.Pr.%s/%04d"%(dattype,Year)
  iPath     = srcDir + "/num.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncltype,ny,nx)
  if binPr == 0.0:
    a3out = zeros([ncltype,ny,nx],int32)
  else:
    a3out = fromfile(iPath, int32).reshape(ncltype, ny, nx)

  return ma.masked_where(a3dommask==-9999., a3out)

def accDat(lYM, binPr, sumnum="num"):
  dfunc = {"sum":loadSum
          ,"num":loadNum
          }
  ddtype= {"sum":"float32"
          ,"num":"int32"
          }

  accDat = zeros([ncltype, ny, nx], ddtype[sumnum])

  for (Year,Mon) in lYM:
    accDat = accDat + dfunc[sumnum](Year,Mon,binPr)
  return accDat

def ret_pdf(lBin,lNum):
  llBin  = [[0,lBin[0]]] + zip(lBin[:-1],lBin[1:])
  lMid   = [mean(Bin) for Bin in llBin]
  lWidth = array([Bin[1]-Bin[0] for Bin in llBin], float32)
  lpdf   = array(lNum)/sum(lNum).astype(float32)/lWidth
  return lMid, lpdf

def ret_Cnt(lBin,lNum, lSum, icltype):
  if icltype != 99:
    a2CLnum = aveCLnum(iYM,eYM,icltype)
  elif icltype ==99:
    a2CLnum = array([aveCLnum(iYM,eYM,icltype) for icltype in lcltype]).sum(axis=0)

  llBin  = [[0,lBin[0]]] + zip(lBin[:-1],lBin[1:])
  lMid   = [mean(Bin) for Bin in llBin]
  lWidth = array([Bin[1]-Bin[0] for Bin in llBin], float32)
  lpdf   = array(lNum)/sum(lNum).astype(float32)/lWidth

  lMean  = (ma.masked_where(lNum==0.0, lSum)/lNum).filled(0.0)
  lCnt   = lpdf *  a2CLnum.mean() * lMean
  return lMid, lCnt 

def mk_Fig(lMid, dY, icltype, figtype="pdf"):
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

  linewidth = [2 if dattype in ["GSMaP.IR","GSMaP.MW"] else 2
                        for dattype in ldattype]

  # Plot
  lines =[axplot.plot(lMid, dY[dattype,icltype]
              ,linewidth=linewidth[idattype], c=colors[idattype]
              ,linestyle=linestyle[idattype])
         for idattype, dattype in enumerate(ldattype)]

  # Axis
  if clVer == "JMA1":
    if figtype=="pdf":
      if   icltype ==1: 
        axplot.set_ylim(0.0, 0.3)
        axplot.set_xlim(0.0, 20.0)
      elif icltype ==7:
        axplot.set_ylim(0.0, 0.05)
        axplot.set_xlim(0.0, 20.0)
  
      else:
        axplot.set_ylim(0.0, 0.05)
        axplot.set_xlim(0.0, 6.0)
  
    elif figtype=="cnt":
      if   icltype ==1: 
        axplot.set_xlim(0.0, 20.0)
      elif icltype ==7:
        axplot.set_xlim(0.0, 20.0)
      else:
        axplot.set_xlim(0.0, 6.0)
  
    elif figtype =="Weak.pdf":
      axplot.set_xlim(0.0, 3.0)
  
    elif figtype =="Heavy.pdf":
      if icltype==1:
        axplot.set_ylim(0.0, 0.05)
      else:
        axplot.set_ylim(0.0, 1e-3)
  
      axplot.set_xlim(10, lbinPr[-2])
  
    elif figtype =="log.Weak.pdf":
      plt.yscale("log")
      if icltype==1:
        axplot.set_ylim(1.e-3, 1.e+0)
        axplot.set_xlim(0.0, 3.0)
      else:
        axplot.set_ylim(1.e-6, 1.e+1)
        axplot.set_xlim(0.0, 3.0)
  
    elif figtype =="log.Heavy.pdf":
      plt.yscale("log")
      if icltype==1:
        axplot.set_ylim(1.e-5, 1.e-1)
      else:
        axplot.set_ylim(1.e-8, 1.e-3)
      axplot.set_xlim(10.0, lbinPr[-2])
    

  elif clVer[:5] == "MyWNP":
    if figtype=="pdf":
      axplot.set_ylim(0.0, 0.3)
      axplot.set_xlim(0.0, 20.0)
 
    elif figtype == "pdf.focus":
      axplot.set_ylim(0.0, 0.05)
      axplot.set_xlim(0.0, 6.0)
 
    elif figtype=="cnt":
        axplot.set_xlim(0.0, 20.0)

    elif figtype=="cnt.focus":
        axplot.set_xlim(0.0, 6.0)
  
    elif figtype =="Weak.pdf":
      axplot.set_xlim(0.0, 3.0)
  
    elif figtype =="Heavy.pdf":
      if icltype in [1,2]:
        axplot.set_ylim(0.0, 0.05)
      else:
        axplot.set_ylim(0.0, 1e-3)
  
      axplot.set_xlim(10, lbinPr[-2])
  
    elif figtype =="log.Weak.pdf":
      plt.yscale("log")
      if icltype in [1,2]:
        axplot.set_ylim(1.e-3, 1.e+0)
        axplot.set_xlim(0.0, 3.0)
      else:
        axplot.set_ylim(1.e-6, 1.e+1)
        axplot.set_xlim(0.0, 3.0)
  
    elif figtype =="log.Heavy.pdf":
      plt.yscale("log")
      if icltype in [1,2]:
        axplot.set_ylim(1.e-5, 1.e-1)
      else:
        axplot.set_ylim(1.e-8, 1.e-3)
      axplot.set_xlim(10.0, lbinPr[-2])
     
  else:
    print "check figtype", figtype
    sys.exit()

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
  #sDir  = "/home/utsumi/mnt/well.share/PMM/WNP.261x265/pict"
  #sDir  = "/home/utsumi/mnt/well.share/PMM/WNP.261x265/pict.MyCL"
  sDir  = ibaseDir + "/pict"
  util.mk_dir(sDir)
  if dommask==None:
    sPath = sDir + "/%s.%s.png"%(figtype, dclShortName[icltype]) 
  else:
    sPath = sDir + "/%sdom.%s.%s.png"%(dommask, figtype, dclShortName[icltype]) 
  plt.savefig(sPath)
  print sPath

  plt.close()
  # Legend file
  legPath = sDir + "/legend.pdfs.png"
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

dPDF = {}
dCnt = {}

for dattype in ldattype:
  print dattype
  for (binmin, binmax) in llbinPr:
    a3num0 = accDat(lYM, binmin, sumnum="num")
    a3num1 = accDat(lYM, binmax, sumnum="num")

    a3sum0 = accDat(lYM, binmin, sumnum="sum")
    a3sum1 = accDat(lYM, binmax, sumnum="sum")

    a3num  = a3num1 - a3num0
    a3sum  = a3sum1 - a3sum0

    for icltype in lcltype:
      dNum[dattype,icltype].append(a3num[icltype].sum())
      dSum[dattype,icltype].append(a3sum[icltype].sum())

    
    dNum[dattype,99].append(a3num.sum())
    dSum[dattype,99].append(a3sum.sum())

  
  for icltype in lcltype + [99]:
    lMid, dPDF[dattype,icltype] = ret_pdf(lbinPr, dNum[dattype,icltype])
    lMid, dCnt[dattype,icltype] = ret_Cnt(lbinPr, dNum[dattype,icltype], dSum[dattype,icltype], icltype)

#-----------
# Figure 
for icltype in lcltype + [99]:
  mk_Fig(lMid, dPDF, icltype, figtype="pdf")
  mk_Fig(lMid, dPDF, icltype, figtype="pdf.focus")
  #mk_Fig(lMid, dPDF, icltype, figtype="Weak.pdf")
  #mk_Fig(lMid, dPDF, icltype, figtype="Heavy.pdf")
  mk_Fig(lMid, dPDF, icltype, figtype="log.Weak.pdf")
  mk_Fig(lMid, dPDF, icltype, figtype="log.Heavy.pdf")
  mk_Fig(lMid, dCnt, icltype, figtype="cnt")
  mk_Fig(lMid, dCnt, icltype, figtype="cnt.focus")
  #mk_Fig(lMid, dCnt, icltype, figtype="cnt"
