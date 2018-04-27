import matplotlib
matplotlib.use("Agg")
from numpy import *
import os, sys
import myfunc.util         as util
import myfunc.IO           as IO
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.IO.CloudSat  as CloudSat
import myfunc.IO.CloudSat.util as CSutil
import matplotlib.pyplot   as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable


#-- CloudSat ---
prdLev  = "2B"
#prdName = "GEOPROF"
prdName = "CLDCLASS"
prdVer  = "P_R04"
varName = {"GEOPROF" :"Radar_Reflectivity"
          ,"CLDCLASS":"cloud_scenario"
          }[prdName]

cs  = CloudSat.CloudSat(prdLev, prdName, prdVer)
nbin= cs.nbin
liz = arange(nbin)
lz  = arange(-4810, 24950+1, 240) /1000. # [km]

#-- JMA-Cloud ---
clVer = "MyWNP.M.3"
rootDir = "/home/utsumi/mnt/wellshare"
if clVer[:5] == "MyWNP":
  ver        = clVer[5:]
  cl         = CLOUDTYPE.MyCloudWNP(ver=ver)
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.My%s"%(ver)
  ibaseDirCL = rootDir + "/CLOUDTYPE/%s"%(clVer)
  [[lllat,lllon],[urlat,urlon]] = cl.BBox
  dlat       = cl.dLat
  dlon       = cl.dLon
#-- CloudNames --
dclName      = cl.dclName
dclShortName = cl.dclShortName
lclid       = sort(dclName.keys())

dcsName      = CSutil.ret_dclName()
dcsShortName = CSutil.ret_dclShortName()
lcsid        = sort(dcsName.keys())

#iYM    = [2014,4]
#eYM    = [2015,7]
iYM    = [2014,4]
eYM    = [2014,4]



lYM    = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [11,12,1,2,3]]

da2Num   = {clid: zeros([len(lcsid),nbin], int)
                   for clid in lclid}

dEvents  = {clid: 0 for clid in lclid}

da2Num2  = {clid: zeros(2, int) for clid in lclid}

for Year,Mon in lYM:
  srcDir = os.path.join(ibaseDir, "Trc.%s"%(prdName), "%04d.%02d"%(Year,Mon))
  strYM    = "%04d.%02d"%(Year,Mon)
  pathMapcl= os.path.join(srcDir, "mapcl.%s.bn"%(strYM))
  pathProf = os.path.join(srcDir, "ctype.%s.%dlevs.bn"%(strYM,nbin))
  #pathLat  = os.path.join(srcDir, "lat.%s.bn"  %(strYM))
  #pathLon  = os.path.join(srcDir, "lon.%s.bn"  %(strYM))
  #pathTime = os.path.join(srcDir, "tstmp.%s.bn"%(strYM))

  aMapcl   = fromfile(pathMapcl,int16)
  aProf    = fromfile(pathProf, int16).reshape(-1,nbin)
  #aLat     = fromfile(pathLat,  float64)
  #aLon     = fromfile(pathLon,  float64)
  #aDTime   = util.tstmp2dtime(fromfile(pathTime,int32))


  ## test -------
  #itmp    = 10000
  #aMapcl  = aMapcl[itmp:itmp+300]
  #aProf   = aProf[itmp:itmp+300]
  #-------------
  ithres = 79  # 24950 - 240*(79-1) = 6230m
  aProf2 = vstack([aProf[:,ithres:].sum(axis=1)
                  ,aProf[:,:ithres].sum(axis=1)])  # mid-low, then high
  aProf2 = aProf2.T

  aProf2 = ma.masked_greater(aProf2, 0).filled(1)

  for clid in lclid:
  #for clid in [4]:
    print Year,Mon,clid
    clidx  = ma.masked_where(aMapcl !=clid, arange(len(aMapcl)))  # make index array for clid
    aprof  = aProf[clidx.compressed(),:]

    dEvents[clid] += sum(~clidx.mask)

    if len(aprof) == 0:
      a2num  = zeros([len(lcsid), nbin],int)
    elif aprof.ndim==2:
      a2num  = array([histogram(aprof[:,iz],list(lcsid)+[99999])[0]
                  for iz in liz], int).T
    elif aprof.ndim==1:      
      a2num =  array([ma.masked_where(aprof !=csid, ones([nbin],int)).filled(0)
                  for csid in lcsid])

    da2Num[clid] = da2Num[clid] + a2num

    #-----------------------------
    # count cloud existence probatility at two levels
    #-----------------------------
    aprof2  = aProf2[clidx.compressed(),:]
    if len(aprof2) ==0:
      a1num2 = zeros(2,int)
    else:
      a1num2 = aprof2.sum(axis=0)

    da2Num2[clid] = da2Num2[clid] + a1num2



## test --------------
#da2Num ={}
#dEvents={}
#for clid in lclid:
#    #da2Num[clid] = arange(nbin*(len(lcsid)+1)).reshape(len(lcsid)+1,nbin)
#    #dEvents[clid]= nbin*(len(lcsid)+1)
#    da2Num[clid], tmp = meshgrid(arange(nbin), ones(len(lcsid)))
#    dEvents[clid]= 124
##--------------------

# Draw
print "Figure"
for clid in lclid[:-1]:
  fig  = plt.figure(figsize=(3.4,1.4))
  ax   = fig.add_axes([0.1,0.19,0.62,0.67])
  x    = arange(len(lcsid)+1)[1:] -0.5   # don't show csid=0:no cloud
  #y    = arange(nbin+1) - 0.5
  y    = lz - (lz[1]-lz[0])*0.5      # height bin size=240 m
  X,Y  = meshgrid(x,y)
  Z    = flipud(da2Num[clid][1:].T)  # don't show csid=0: no cloud
  print "*"*50
  Z    = ma.masked_equal(Z,0) 
  Z    = Z / float(dEvents[clid]) * 100. # [%]
  print Z
  #cmap = "binary"
  cmap = "jet"
  if clid in [1,2]:
    im   = ax.pcolormesh(X, Y, Z, cmap=cmap, vmin=0, vmax=100)
  else:
    im   = ax.pcolormesh(X, Y, Z, cmap=cmap, vmin=0, vmax=50)
   
  # Axis limit
  #ax.set_ylim([0,y[-1]])
  ax.set_ylim([0,20])
  ax.set_xlim([x[0],x[-1]])

  # X-tick label (for meshgrid)
  lxlabel = [dcsShortName[i] for i in lcsid[1:]]
  ax.set_xticks(x[:-1]+0.5)
  ax.xaxis.set_ticklabels(lxlabel, fontsize=13)

  # Colorbar  -------------------------------
  divider = make_axes_locatable(ax)
  cbpos= fig.add_axes([0.9,0.19,0.01,0.67])
  cb   = plt.colorbar(im, cax=cbpos)
  cb.ax.tick_params(labelsize=9)

  # Total count profile  --------------------
  ax2  = fig.add_axes([0.75,0.19,0.04,0.67])
  x2   = Z.sum(axis=1)
  y2   = y
  ax2.plot(x2,y2,"-",color="k",linewidth=2)
  ax2.tick_params(labelleft="off")
  #ax2.set_ylim([0,y[-1]]) 
  ax2.set_ylim([0,20]) 
  if clid in [1,2]:
    xmax2 =100
  else:
    xmax2 =59
  ax2.set_xlim([0,xmax2]) 

  # Total count at 2-layers  ----------------
  ax3  = fig.add_axes([0.82,0.19,0.04,0.67])
  x3   = da2Num2[clid] / float(dEvents[clid])*100.
  y3   = array([0,6])+1
  lwidth= [4,12]
  ax3.barh(y3[0],x3[0],4,align="edge",color="gray")  # 3rd para=width
  ax3.barh(y3[1],x3[1],12,align="edge",color="gray") # 3rd para=width
  ax3.tick_params(labelleft="off")
  ax3.set_ylim([0,20]) 
  if clid in [1,2]:
    xmax3 =100
  else:
    #xmax3 =59
    xmax3 =100
  ax3.set_xlim([0,xmax3]) 
  print x3

 
  # Title
  plt.suptitle(dclName[clid])
  
  # Vertical lines
  [ax.plot([xx, xx],[y[0],y[-1]+1],"-",color="k",linewidth=0.8) for xx in x]
  
  # Save
  figDir  = os.path.join(ibaseDir,"pict")
  util.mk_dir(figDir)
  figPath = os.path.join(figDir, "prof.CL.2prob.%s.png"%(dclShortName[clid]))
  plt.savefig(figPath)
  print figPath
  #plt.show()
  
