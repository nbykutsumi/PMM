from numpy import *
import os, sys
import myfunc.util         as util
import myfunc.IO           as IO
import myfunc.IO.CLOUDTYPE as CLOUDTYPE
import myfunc.IO.CloudSat  as CloudSat
import myfunc.IO.CloudSat.util as CSutil
import matplotlib.pyplot   as plt
import matplotlib

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
rootDir = "/home/utsumi/mnt/well.share"
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

iYM    = [2014,4]
eYM    = [2015,7]
#iYM    = [2014,4]
#eYM    = [2014,4]



lYM    = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [11,12,1,2,3]]

da2Num   = {clid: zeros([len(lcsid),nbin], int)
                   for clid in lclid}

dEvents  = {clid: 0 for clid in lclid}

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

  for clid in lclid:
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

# Draw
for clid in lclid[:-1]:
  fig  = plt.figure(figsize=(3.5,1.8))
  ax   = fig.add_axes([0.1,0.15,0.78,0.64])
  x    = arange(len(lcsid)+1)[1:] -0.5   # don't show csid=0:no cloud
  #y    = arange(nbin+1) - 0.5
  y    = lz - (lz[1]-lz[0])*0.5      # height bin size=240 m
  X,Y  = meshgrid(x,y)
  Z    = flipud(da2Num[clid][1:].T)  # don't show csid=0: no cloud
  Z    = ma.masked_equal(Z,0) 
  Z    = Z / float(dEvents[clid]) * 100. # [%]
  if clid in [1,2]:
    im   = ax.pcolormesh(X, Y, Z, vmin=0, vmax=100)
  else:
    im   = ax.pcolormesh(X, Y, Z, vmin=0, vmax=50)
  
  # Colorbar
  cb   = plt.colorbar(im)
  cb.ax.tick_params(labelsize=9)

  # Axis limit
  plt.ylim([0,y[-1]])
  plt.xlim([x[0],x[-1]])
  
  # X-tick label
  lxlabel = [dcsShortName[i] for i in lcsid]
  ax.xaxis.set_ticklabels(lxlabel, fontsize=14)
  
  # Title
  plt.title(dclName[clid])
  
  # Vertical lines
  [ax.plot([xx, xx],[y[0],y[-1]+1],"-",color="k") for xx in x]
  
  # Save
  figDir  = os.path.join(ibaseDir,"pict")
  util.mk_dir(figDir)
  figPath = os.path.join(figDir, "prof.CL.%s.png"%(dclShortName[clid]))
  plt.savefig(figPath)
  print figPath
  #plt.show()
  
