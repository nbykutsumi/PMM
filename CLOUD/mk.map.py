from numpy import *
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.fig.Fig as Fig
import myfunc.util as util
import matplotlib.pyplot as plt

iYM    = [2014,4]
eYM    = [2014,12]
lYM    = util.ret_lYM(iYM, eYM)
ldattype = ["RA","KuPR","GMI","GSMaP","GSMaP.IR","GSMaP.MW","IMERG"]
#ldattype = ["KuPR"]

ny,nx   = 261, 265
miss    = -9999.
lcltype = range(0,7+1)
ncltype = len(lcltype)
dclName ={0:"Clear Sky",   1:"Cumulonimbus(Cb)",  2:"High Cloud",3:"Mid Cloud"
         ,4:"Cumulus(Cu)", 5:"Stratocumulus(Sc)", 6:"Fog/St"    ,7:"Cloudy"}

dclShortName={0:"no", 1:"Cb",  2:"hi",3:"md"
             ,4:"Cu", 5:"Sc",  6:"St",7:"cw"}

cl   = CLOUDTYPE.CloudWNP()
Lat  = cl.Lat
Lon  = cl.Lon
ny   = cl.ny
nx   = cl.nx
BBox = cl.BBox
   
#******************************************
def loadSum(dattype,Year,Mon,binPr):
  baseDir = "/tank/utsumi/PMM/WNP.261x265"
  srcDir    = baseDir + "/CL.Pr.%s"%(dattype)
  iPath     = srcDir + "/sum.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncltype,ny,nx)

  if binPr == 0.0:
    return zeros([ncltype, ny, nx], int32)
  else:
    return fromfile(iPath, float32).reshape(ncltype, ny, nx)

def loadNum(dattype,Year,Mon,binPr):
  baseDir = "/tank/utsumi/PMM/WNP.261x265"
  srcDir    = baseDir + "/CL.Pr.%s"%(dattype)
  iPath     = srcDir + "/num.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncltype,ny,nx)

  if binPr == 0.0:
    return zeros([ncltype, ny, nx], int32)
  else:
    return fromfile(iPath, int32).reshape(ncltype, ny, nx)

def accDat(dattype,iYM, eYM, binPr, sumnum="num"):
  dfunc = {"sum":loadSum
          ,"num":loadNum
          }
  ddtype= {"sum":"float32"
          ,"num":"int32"
          }

  accDat = zeros([ncltype, ny, nx], ddtype[sumnum])

  lYM = util.ret_lYM(iYM,eYM)
  for (Year,Mon) in lYM:
    accDat = accDat + dfunc[sumnum](dattype,Year,Mon,binPr)
  return accDat

#******************************************
for dattype in ldattype:
  a3sum = accDat(dattype, iYM, eYM, 999.0, sumnum="sum")
  a3num = accDat(dattype, iYM, eYM, 999.0, sumnum="num")
  a3rat = (ma.masked_where(a3num==0.0, a3sum)/a3num).filled(miss)
  for icl in lcltype:
  #for icl in [7]:
    a2rat = a3rat[icl]

    # Names
    sDir     = "/tank/utsumi/PMM/WNP.261x265/pict"
    figname  = sDir + "/Map.prate.%s.%s.png"%(dattype, dclShortName[icl])

    # Colorbar
    if icl in [1,7]:
      bnd      = [1,3,5,7,9,11,13,15,17,19,21]
      cbarname = sDir + "/cbar.prate.cbcw.png"
    else: 
      bnd      = [0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
      cbarname = sDir + "/cbar.prate.png"
    mycm     = "jet"

    # Title
    stitle   = "%s %s"%(dattype,dclName[icl])

    # Draw
    Fig.DrawMap(a2rat, Lat, Lon, BBox, bnd=bnd, lowest_white=True,mycm=mycm,figname=figname, stitle=stitle, cbarname=cbarname )
    print figname
    plt.close()
    plt.clf()
