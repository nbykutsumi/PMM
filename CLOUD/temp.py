import matplotlib
matplotlib.use("Agg")
from numpy import *
from collections import deque
import myfunc.util as util
import matplotlib.pylab as plt
import myfunc.IO.CLOUDTYPE as CLOUDTYPE


expr = "std"
#expr = "sht"
#expr = "old"

#logax = True
logax = False

if expr =="sht":
    iYM    = [2014,4]
    eYM    = [2014,10]
else:
    iYM    = [2014,4]
    eYM    = [2015,6]
    #eYM    = [2014,4]

lYM    = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [11,12,1,2,3]]
#ldattype = ["KuPR","GMI","IMERG","IMERG.IR","IMERG.MW","GSMaP","GSMaP.IR","GSMaP.MW"]
#ldattype = ["RA","GMI","IMERG.MW","IMERG.IR","GSMaP.MW","GSMaP.IR"]
#ldattype = ["GMI","IMERG","IMERG.MW","IMERG.IR","GSMaP","GSMaP.MW","GSMaP.IR"]
ldattype = ["GSMaP.MW"]

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
  if expr=="old":
    ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.My%s/old@20170925"%(ver)
  else:
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
#llndsea = ["sea","lnd"]
llndsea = ["sea"]

dwidth = {}
for dattype in ldattype:
    if dattype in ["GMI"]:
        dwidth[dattype] = 3
    elif dattype in ["GSMaP.MW"]:
        dwidth[dattype] = 1.2
    else:
        dwidth[dattype] = 2
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
if logax ==True:
    vmin    = 0.05
    vmax    = 100  # mm/hour
elif logax==False:
    vmin    = 0.
    vmax    = 3  # mm/hour


#cols    = plt.matplotlib.cm.jet(linspace(0,1,len(ldattype)))
#cols    = plt.matplotlib.cm.jet(linspace(0,1,len(ldattype)+1))
cols  = {"KuPR":"k"
        ,"GMI" :"0.2"
        ,"IMERG":"dodgerblue"
        ,"IMERG.MW":"royalblue"
        ,"IMERG.IR":"limegreen"
        ,"GSMaP":"red"
        ,"GSMaP.MW":"magenta"
        ,"GSMaP.IR":"darkorange"
        }



for lndsea in llndsea:
  #for iicl, icl in enumerate(lcltype[1:] + [99]):
  for icl in [4]:
    #** Figure **********
    figplot = plt.figure(figsize=(3,3))
    axplot1  = figplot.add_axes([0.15,0.1,0.8,0.8])
    axplot1.set_ylim(vmin, vmax)
    axplot1.set_xlim(vmin, vmax)
    if logax == True:
        axplot1.set_xscale("log",nonposx="clip")
        axplot1.set_yscale("log",nonposy="clip")

    # ticks
    if logax==False:
      axplot1.xaxis.set_ticks(arange(vmin,vmax+1,0.5))
      axplot1.yaxis.set_ticks(arange(vmin,vmax+1,0.5))


    dmean = {}
    dstd  = {}
    for idattype, dattype in enumerate(ldattype):
      # No plot for GSMaP and IMERG
      if dattype in ["GSMaP","IMERG"]:
        print "*"*50
        print "no plotting for GSMaP and IMERG"
        lpr = array([])
        lref= array([])
        print "*"*50
      # Load data
      elif   icl != 99:
        lpr = loadData(dattype, lndsea, icl, lYM)
        lref= loadX   (dattype, lndsea, icl, lYM)
  
      elif icl ==99:
        lpr = concatenate([loadData(dattype, lndsea, icltmp, lYM) for icltmp in lcltype])
        lref= concatenate([loadX   (dattype, lndsea, icltmp, lYM) for icltmp in lcltype])


      # Average line
      if logax==True:
        lk     = arange(-2,2.6+0.01,0.2)
        bins   = [10**k for k in lk]
      else:
        bins   = hstack([arange(-0.1,2,0.2),arange(2,50,1)])
      BINS   = zip(bins[:-1],bins[1:])
      dmean[dattype] = [ma.masked_where(logical_or( lref<binmin,  binmax<=lref), lpr).mean() for (binmin,binmax) in BINS]
      dstd[dattype]  = [ma.masked_where(logical_or( lref<binmin,  binmax<=lref), lpr).std() for (binmin,binmax) in BINS]
      lx     = [mean(BIN) for BIN in BINS]


    # Draw average lines
    #lines = [axplot1.plot(lx,dmean[dattype],"-", linewidth=dwidth[dattype], color=cols[idattype+1]) 
    lines = [axplot1.plot(lx,dmean[dattype],"-", linewidth=dwidth[dattype], color=cols[dattype]) 
                for idattype,dattype in enumerate(ldattype)]

    for i,x in enumerate(lx):
        print lx[i], dmean[dattype][i]
    print lx    
    # Draw 1-1 line
    axplot1.plot([vmin,vmax],[vmin,vmax],"--",color="k")

    # Add grid
    plt.grid(which="major",linewidth=0.6)
    plt.grid(which="minor",linewidth=0.6)

    # Add title
    stitle = "%04d/%02d-%04d/%02d CL=%s [%s]"%(iYM[0],iYM[1],eYM[0],eYM[1],dclShortName[icl], lndsea)

    if expr !="std":
        stitle = stitle + " %s"%(expr)
    plt.title(stitle, fontsize=10)


    # Save
    if expr =="std":
        figDir = ibaseDir + "/pict"
    elif expr=="sht":
        figDir = ibaseDir + "/pict.sht"
    elif expr=="old":
        figDir = ibaseDir + "/pict.old"
    else:
        print "check expr",expr
        sys.exit()


    if logax==True:
    #    figPath= figDir  + "/lines.mulProd.log.vs%s.%s.%s.png"%(Xvar,lndsea,dclShortName[icl])
        figPath = figDir + "/temp.log.%s.png"%(dclShortName[icl])
    else:
    #    figPath= figDir  + "/lines.mulProd.vs%s.%s.%s.png"%(Xvar,lndsea,dclShortName[icl])
        figPath = figDir + "/temp.%s.png"%(dclShortName[icl])

    util.mk_dir(figDir)
    plt.savefig(figPath)
    print figPath
    #plt.show() 
    plt.close()

    ## Legend file
    #if iicl==0:
    #  legPath = figDir + "/legend.mulProd.vs%s.png"%(Xvar)
    #  figleg  = plt.figure(figsize=(3,3))
    #  lines   = [line[0] for line in lines]   # 2D list to 1D list
    #  figleg.legend(lines, ret_ldatname(ldattype), fontsize=15)
    #  figleg.savefig(legPath)
    #  plt.close()
    #  print legPath
