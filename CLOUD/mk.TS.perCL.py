from numpy import *
from datetime import datetime, timedelta
from bisect import bisect, bisect_left, bisect_right
from matplotlib import cm
import calendar
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.util as util
import matplotlib.pyplot as plt
import calendar
import sys
from scipy.stats import ttest_ind

#clVer  = "JMA1"
#clVer  = "MyWNP1"
#clVer  = "MyWNP2"
clVer  = "MyWNP.M.3"

iYM    = [2014,4]
eYM    = [2015,6]
#eYM    = [2014,8]
lYM    = util.ret_lYM(iYM, eYM)
lYM    = [YM for YM in lYM if YM[1] not in [11,12,1,2,3]]
print lYM
lMName = [calendar.month_abbr[M] for (Y,M) in lYM]
ldate  = ["%04d "%(Y) + "%s"%(calendar.month_abbr[M])
            for (Y,M) in lYM]

ldattype = ["KuPR","GMI","IMERG","IMERG.IR","IMERG.MW","GSMaP","GSMaP.IR","GSMaP.MW"]
#ldattype = ["KuPR","IMERG"]

rootDir = "/home/utsumi/mnt/well.share"
if clVer   == "JMA1":
  cl         = CLOUDTYPE.CloudWNP()
  ncltype = 8
  lcltype = range(ncltype)
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.JMA"

elif clVer[:5] == "MyWNP":
  ver        = clVer[5:]
  cl         = CLOUDTYPE.MyCloudWNP(ver=ver)
  ncltype = cl.ncl
  lcltype = cl.licl
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.My%s"%(ver)
  ibaseDirCL = rootDir + "/CLOUDTYPE/MyWNP%s"%(ver)

#
ny,nx = cl.ny, cl.nx   #261, 265
Lat   = cl.Lat
Lon   = cl.Lon

dommask  = "RA"  # Do NOT Change!!

#BBox    = [[-0.1, 113.875],[52.1, 180.125]]  # WN.Pac
BBox    = [[20., 118.],[40., 150.]]    # RadarAMeDAS
miss  = -9999.

dclName      = cl.dclName
dclShortName = cl.dclShortName
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
  a2dommask = ma.masked_equal(a2dommask, miss)
  a2dommask = ma.masked_where(a2BBox ==miss, a2dommask)
  a3dommask = array([a2dommask.mask for i in range(len(lcltype))])
else:
  a2dommask = a2BBox
  a3dommask = array([a2dommask.mask for i in range(len(lcltype))])



# land sea mask
landfracDir = "/home/utsumi/mnt/well.share/PMM/WNP.261x265/MASK"
landfracPath= landfracDir + "/landfrac.%dx%d"%(cl.ny, cl.nx)
a2lndfrc    = fromfile(landfracPath, float32).reshape(261,265)

a2lndmask   = ma.masked_less(a2lndfrc,1.0).mask  # mask sea&coast
a2seamask   = ma.masked_greater(a2lndfrc,0.0).mask # mask land&coast
a2cstmask   = ma.masked_equal(a2lndfrc, 1.0).mask \
             +ma.masked_equal(a2lndfrc, 0.0).mask  # mask land&sea
a2anymask   = array([False]*ny*nx).reshape(ny,nx)

a3lndmask   = array([a2lndmask for i in range(len(lcltype))])
a3seamask   = array([a2seamask for i in range(len(lcltype))])
a3cstmask   = array([a2cstmask for i in range(len(lcltype))])
a3anymask   = array([a2anymask for i in range(len(lcltype))])

da2lsmask   = {"lnd":a2lndmask, "sea":a2seamask, "cst":a2cstmask,"any":a2anymask}

da3lsmask   = {"lnd":a3lndmask, "sea":a3seamask, "cst":a3cstmask,"any":a3anymask}

llndsea     = ["lnd","sea","cst","any"]

#******************************************
def loadSum(dattype,Year,Mon,binPr):
  srcDir    = ibaseDir + "/CL.Pr.%s/%04d"%(dattype,Year)
  iPath     = srcDir + "/sum.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncltype,ny,nx)
  a3out = fromfile(iPath, float32).reshape(ncltype, ny, nx)

  return ma.masked_where(a3dommask==-9999., a3out)

def loadNum(dattype,Year,Mon,binPr):
  srcDir    = ibaseDir + "/CL.Pr.%s/%04d"%(dattype,Year)
  iPath     = srcDir + "/num.P%05.1f.%04d.%02d.%dx%dx%d"%(binPr, Year,Mon, ncltype,ny,nx)
  a3out = fromfile(iPath, int32).reshape(ncltype, ny, nx)
  return ma.masked_where(a3dommask==-9999., a3out)

def biasTS(lYM):
  binPr = 999
  dbias = {(icl,dattype,lndsea): [] 
            for icl     in lcltype 
            for dattype in ldattype
            for lndsea in ["any","lnd","sea"]
           }
  for dattype in ldattype:
    for Year,Mon in lYM:
      a3sum   = loadSum(dattype,Year,Mon,binPr)
      a3num   = loadNum(dattype,Year,Mon,binPr)
      a3sumRA = loadSum("RA"   ,Year,Mon,binPr)
      a3numRA = loadNum("RA"   ,Year,Mon,binPr)
  
      for lndsea in ["any","lnd","sea"]:
        a3mask= da3lsmask[lndsea] + a3dommask
        a1prat= ma.masked_where(a3mask, a3sum).sum(axis=(1,2))\
               /ma.masked_where(a3mask, a3num).sum(axis=(1,2))
  
        a1pratRA= ma.masked_where(a3mask, a3sumRA).sum(axis=(1,2))\
               /ma.masked_where(a3mask, a3numRA).sum(axis=(1,2))
  
        a1bias = ma.masked_invalid((a1prat - a1pratRA)/a1pratRA)*100.  # [%]
        [dbias[(icl,dattype,lndsea)].append(a1bias[icl])
             for icl in lcltype]

  return dbias


def ratTS(ldattype, lYM, binMin=False):
  dout = {(icl,dattype,lndsea): [] 
            for icl     in lcltype + [99]
            for dattype in ldattype
            for lndsea in ["any","lnd","sea"]
           }
  for dattype in ldattype:
    for Year,Mon in lYM:
      a3sum   = loadSum(dattype,Year,Mon,999)
      a3num   = loadNum(dattype,Year,Mon,999)

      if binMin not in [False, 0]:
        a3sumMin   = loadSum(dattype,Year,Mon,binMin)
        a3numMin   = loadNum(dattype,Year,Mon,binMin)

        a3sum      = a3sum - a3sumMin 
        a3num      = a3num - a3numMin 
 
      for lndsea in ["any","lnd","sea"]:
        a3mask= da3lsmask[lndsea] + a3dommask
        a1prat= ma.masked_where(a3mask, a3sum).sum(axis=(1,2))\
               /ma.masked_where(a3mask, a3num).sum(axis=(1,2))

        a1out = a1prat # [mm/hour] 
 
        [dout[(icl,dattype,lndsea)].append(a1out[icl])
             for icl in lcltype]

        dout[(99,dattype,lndsea)].append(\
          ma.masked_where(a3mask[0], a3sum.sum(axis=0)).sum()\
         /ma.masked_where(a3mask[0], a3num.sum(axis=0)).sum()
         ) 

  return dout

 
#******************************************

# Color
#lcol = cm.rainbow(linspace(0,1, 4))
lcol = ["r","cyan","purple","orange"]
dcol = {"RA": "k"
       ,"KuPR": lcol[0] 
       ,"GMI" : lcol[1]
       ,"IMERG":lcol[2]
       ,"IMERG.IR":lcol[2]
       ,"IMERG.MW":lcol[2]
       ,"GSMaP":lcol[3]
       ,"GSMaP.IR":lcol[3]
       ,"GSMaP.MW":lcol[3]
       }

# Line width
w0 = 1
w1 = 2.5
w2 = 1
dwidth = {"RA": w0
       ,"KuPR": w0
       ,"GMI" : w0
       ,"IMERG":w0
       ,"IMERG.IR":w1
       ,"IMERG.MW":w2
       ,"GSMaP":   w0
       ,"GSMaP.IR":w1
       ,"GSMaP.MW":w2
       }

# Line style
dstyle={"RA": "-"
       ,"KuPR": "-"
       ,"GMI" : "-"
       ,"IMERG":"-"
       ,"IMERG.IR":"-"
       ,"IMERG.MW":"--"
       ,"GSMaP":   "-"
       ,"GSMaP.IR":"-"
       ,"GSMaP.MW":"--"
       }

lcltype_tmp = [99] + lcltype[1:]
llndsea_tmp = ["lnd"]
ldattype_tmp= ["RA"] + ldattype
binMin = 0
#binMin= 0.1
#binMin= 0.5

drat = ratTS(ldattype_tmp, lYM, binMin=binMin)
print drat
for lndsea in llndsea_tmp:
  for icl in lcltype_tmp:
  #for icl in [3]:
    figplot = plt.figure(figsize=(5,3))
    axplot  = figplot.add_axes([0.1,0.24,0.8,0.66])

    X       = range(len(lYM))
    [plt.plot(X, drat[icl,dattype,lndsea]
         ,color     = dcol  [dattype]
         ,linewidth = dwidth[dattype]
         ,linestyle = dstyle[dattype]
         )
         for dattype in ldattype_tmp]

    # zero line
    plt.plot([X[0]-0.1,X[-1]+1.1],[0,0],"-",color="k")

    # X-limit
    axplot.set_xlim(X[0],X[-1])

    # X-ticks
    plt.xticks(X, ldate, fontsize=10, rotation=60)

    # Title
    plt.title("%s (%s) min=%.1f[mm/h]"%(dclName[icl],lndsea, binMin), fontsize=12)

    # Names & Save
    sDir  = ibaseDir + "/pict"
    sPath = sDir + "/TS.rat.min%.1f.%s.%s.png"%(binMin,lndsea,dclShortName[icl])
    figplot.savefig(sPath)

    print sPath 

# Legend file
figlegend  = plt.figure(figsize=(2,3))
lines      = axplot.get_lines()
figlegend.legend(lines, ldattype_tmp)
params     = {"legend.fontsize":30
            , "lines.linewidth":2}
plt.rcParams.update(params)
legendPath = sDir + "/legend.TS.rat.perCL.png"
figlegend.savefig(legendPath) 

