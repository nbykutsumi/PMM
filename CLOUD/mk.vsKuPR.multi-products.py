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
#ldattype = ["RA","GMI","IMERG.MW","IMERG.IR","GSMaP.MW","GSMaP.IR"]
ldattype = ["GMI","IMERG","IMERG.MW","IMERG.IR","GSMaP","GSMaP.MW","GSMaP.IR"]

#clVer = "JMA1"
clVer = "MyWNP.M.3"

Xvar  = "KuPR"

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
llndsea = ["sea","lnd"]
#*******************************
def ret_ldatname(ldattype):
    ldatname = []
    for dattype in ldattype:
        if dattype.split(".")[-1] == "MW":
            datname = dattype.split(".")[0]+".PMW"
        else:
            datname = dattype
        ldatname.append(datname)
    return ldatname

def loadData(dattype, lndsea, icl, lYM):
  lpr = deque([])
  for Year,Mon in lYM:
    baseDir  = ibaseDir
    sDir     = baseDir + "/VsKuPR.CL.%s/%04d"%(dattype,Year)

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

if Xvar == "KuPR": loadX = loadKu

#*******************************
vlim    = 40  # mm/hour
#cols    = plt.matplotlib.cm.jet(linspace(0,1,len(ldattype)))
cols    = plt.matplotlib.cm.jet(linspace(0,1,len(ldattype)+1))
#coldarks= [col*array([1,1,1,0.5]) for col in cols]

for lndsea in llndsea:
  for iicl, icl in enumerate(lcltype[1:] + [99]):
  #for icl in [99]:
    #** Figure **********
    figplot = plt.figure(figsize=(3,3))
    axplot1  = figplot.add_axes([0.1,0.1,0.8,0.8])
    axplot1.set_ylim(0.0, vlim)
    axplot1.set_xlim(0.0, vlim)

    # ticks
    axplot1.xaxis.set_ticks(arange(0,vlim+1,5)) 
    axplot1.yaxis.set_ticks(arange(0,vlim+1,5)) 

    dmean = {}
    dstd  = {}
    for idattype, dattype in enumerate(ldattype):
      # Load data
      if   icl != 99:
        lpr = loadData(dattype, lndsea, icl, lYM)
        lref= loadX   (dattype, lndsea, icl, lYM)
  
      elif icl ==99:
        lpr = concatenate([loadData(dattype, lndsea, icltmp, lYM) for icltmp in lcltype])
        lref= concatenate([loadX   (dattype, lndsea, icltmp, lYM) for icltmp in lcltype])

      # No plot for GSMaP and IMERG
      if dattype in ["GSMaP","IMERG"]:
        print "*"*50
        print "no plotting for GSMaP and IMERG"
        lpr = array([])
        lref= array([])
        print "*"*50

      # Average line
      bins   = arange(0, 50,2)
      #bins   = r_[arange(0, 10+0.1,2), arange(15,100,5)]
      BINS   = zip(bins[:-1],bins[1:])
      dmean[dattype] = [ma.masked_where(logical_or( lref<binmin,  binmax<=lref), lpr).mean() for (binmin,binmax) in BINS]
      dstd[dattype]  = [ma.masked_where(logical_or( lref<binmin,  binmax<=lref), lpr).std() for (binmin,binmax) in BINS]
      lx     = [mean(BIN) for BIN in BINS]

    # Draw average lines
    #lines = [axplot1.plot(lx,dmean[dattype],"-", linewidth=2, color=cols[idattype]) 
    lines = [axplot1.plot(lx,dmean[dattype],"-", linewidth=2, color=cols[idattype+1]) 
                for idattype,dattype in enumerate(ldattype)]

    ## Error Range
    #for idattype, dattype in enumerate(ldattype):
    #  lx_err = array(lx[2:3] + lx[5::2])+0.5*(i-2)
    #  ly_err = dmean[dattype][2:3] + dmean[dattype][5::2]
    #  lerr   = dstd[dattype][2:3] + dstd[dattype][5::2]
    #  axplot1.errorbar(lx_err, ly_err, lerr, ecolor=coldarks[idattype]) 
    #print "*"*50
    #print "lx=",lx
    #print "lx_err=",lx_err
    #print "*"*50

    # Draw 1-1 line
    axplot1.plot([0,100],[0,100],"--",color="k")


    # Add title
    stitle = "%04d/%02d-%04d/%02d CL=%s [%s]"%(iYM[0],iYM[1],eYM[0],eYM[1],dclShortName[icl], lndsea)
    plt.title(stitle, fontsize=10)
    # Save
    figDir = ibaseDir + "/pict"
    figPath= figDir  + "/lines.mulProd.vs%s.%s.%s.png"%(Xvar,lndsea,dclShortName[icl])

    util.mk_dir(figDir)
    plt.savefig(figPath)
    print figPath
    #plt.show() 
    #plt.close()

    # Legend file
    if iicl==0:
      legPath = figDir + "/legend.mulProd.vs%s.png"%(Xvar)
      figleg  = plt.figure(figsize=(2.2,2.4))
      lines   = [line[0] for line in lines]   # 2D list to 1D list
      figleg.legend(lines, ret_ldatname(ldattype))
      figleg.savefig(legPath)
      plt.close()
      print legPath
