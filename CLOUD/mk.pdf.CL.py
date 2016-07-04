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
ncltype = 8
lcltype = range(8)
dclName ={0:"Clear Sky",   1:"Cumulonimbus(Cb)",  2:"High Cloud",3:"Mid Cloud"
         ,4:"Cumulus(Cu)", 5:"Stratocumulus(Sc)", 6:"Fog/St"    ,7:"Cloudy"}

dclShortName={0:"no", 1:"Cb",  2:"hi",3:"md"
             ,4:"Cu", 5:"Sc",  6:"St",7:"cw"}

llbinPr = [[0, lbinPr[0]]] + zip(lbinPr[0:-1], lbinPr[1:])
#----------------------
def loadSum(Year,Mon,binPr):
  baseDir = "/tank/utsumi/PMM/WNP.261x265"
  srcDir    = baseDir + "/CL.Pr.%s"%(dattype)
  iPath     = srcDir + "/sum.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncltype,ny,nx)

  if binPr == 0.0:
    return zeros([ncltype, ny, nx], int32)
  else:
    return fromfile(iPath, float32).reshape(ncltype, ny, nx)

def loadNum(Year,Mon,binPr):
  baseDir = "/tank/utsumi/PMM/WNP.261x265"
  srcDir    = baseDir + "/CL.Pr.%s"%(dattype)
  iPath     = srcDir + "/num.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncltype,ny,nx)

  if binPr == 0.0:
    return zeros([ncltype, ny, nx], int32)
  else:
    return fromfile(iPath, int32).reshape(ncltype, ny, nx)

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

def ret_pdf(lBin,lNum):
  llBin  = [[0,lBin[0]]] + zip(lBin[:-1],lBin[1:])
  lMid   = [mean(Bin) for Bin in llBin]
  lWidth = array([Bin[1]-Bin[0] for Bin in llBin], float32)
  lpdf   = array(lNum)/sum(lNum).astype(float32)/lWidth
  return lMid, lpdf

def mk_Fig(lMid, dPDF, icltype, yscale=None):
  figplot = plt.figure()
  axplot  = figplot.add_axes([0.2,0.2,0.7,0.7])

  # Colors
  colors = plt.matplotlib.cm.jet(linspace(0,1,len(ldattype)))

  # Line style
  linestyle=["--" if dattype in ["KuPR","GMI"] else "-"
             for dattype in ldattype]

  # Plot
  [axplot.plot(lMid, dPDF[dattype,icltype]
              ,linewidth=4, c=colors[idattype]
              ,linestyle=linestyle[idattype])
   for idattype, dattype in enumerate(ldattype)]

  # Axis
  if yscale ==None:
    if   icltype ==1: 
      axplot.set_ylim(0.0, 0.3)
      axplot.set_xlim(0.0, 30.0)
    elif icltype ==7:
      axplot.set_ylim(0.0, 0.02)
      axplot.set_xlim(0.0, 30.0)
    else:
      axplot.set_ylim(0.0, 0.02)
      axplot.set_xlim(0.0, 15.0)

  elif yscale =="log":
    plt.yscale("log")
  else:
    print "check yscale",yscale
    sys.exit()


  # Add legend
  plt.legend(["%s"%(dattype) for dattype in ldattype])

  # Add title
  plt.title("%04d/%02d-%04d/%02d CL=%s"
             %(iYM[0],iYM[1],eYM[0],eYM[1],dclName[icltype]))

  # Save figure
  sDir  = "/tank/utsumi/PMM/WNP.261x265/pict"
  if   yscale ==None:
    sPath = sDir + "/pdf.%s.png"%(dclShortName[icltype]) 
  elif yscale =="log":
    sPath = sDir + "/log.pdf.%s.png"%(dclShortName[icltype]) 
  plt.savefig(sPath)
  print sPath
  #plt.show()


#----------------------
dNum = {(dattype,icltype):[] 
       for dattype in ldattype
       for icltype in lcltype
       }

dPDF = {}

for dattype in ldattype:
  print dattype
  for (binmin, binmax) in llbinPr:
    a3num0 = accDat(iYM, eYM, binmin, sumnum="num")
    a3num1 = accDat(iYM, eYM, binmax, sumnum="num")
    a3num  = a3num1 - a3num0
    for icltype in lcltype:
      dNum[dattype,icltype].append(a3num[icltype].sum())
  
  for icltype in lcltype:
    lMid, dPDF[dattype,icltype] = ret_pdf(lbinPr, dNum[dattype,icltype])

# Figure 
for icltype in lcltype:
  mk_Fig(lMid, dPDF, icltype, yscale=None)
  mk_Fig(lMid, dPDF, icltype, yscale="log")

