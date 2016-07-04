from numpy import *
import matplotlib.pyplot as plt
import ctrack_func
import ctrack_para

iyear,eyear = 2014,2014
imon,emon   = 3,6
#imon,emon   = 3,5
#imon,emon   = 3,4
lyear       = range(iyear,eyear+1)
#lmon        = ctrack_para.ret_lmon(season)
lmon        = range(imon, emon+1)

prj = ["GPM.KaPR","L2"]
#prj = ["GPM.GMI","L2"]
#prj = ["TRMM.TMI","2A12"]
#prj = ["TRMM.PR","2A25"]
#lprmaptype = ["RA","GSMaP.gauge"]
#lprmaptype = ["RA"]
lprmaptype = ["GSMaP.gauge"]

for prmaptype in lprmaptype: 
  a1Obt  = array([],float32)
  a1Map  = array([],float32)
  for year in lyear:
    for mon in lmon:
      idir  = "/mnt/mizu.tank/utsumi/PMM/COMP.ORB.ORG/%04d/%02d"%(year,mon)
      iname = idir + "/%s.x.%s.bin"%(".".join(prj), prmaptype)
      print iname
      a2dat = fromfile(iname, float32).reshape(-1,2)
  
      a1Obt = r_[a1Obt, a2dat[:,0]]
      a1Map = r_[a1Map, a2dat[:,1]]
  
  a1Obt = a1Obt * 60.*60.  # mm/sec --> mm/hour
  a1Map = a1Map * 60.*60.  # mm/sec --> mm/hour
  
  #-- plot ----
  figplot  = plt.figure(figsize=(5.5, 5.5))
  axplot   = figplot.add_axes([0.2, 0.2, 0.7, 0.7])
  
  axplot.plot(a1Map, a1Obt, "o", color="k", markersize=2)
  
  #-- axis ----
  axmax  = 80.
  axplot.set_ylim(0.0, axmax)
  axplot.set_xlim(0.0, axmax)
  
  #-- axis label --
  axplot.set_ylabel("%s [mm/hour]"%(prj[0].split(".")[1]), fontsize=18)
  axplot.set_xlabel("%s [mm/hour]"%(prmaptype), fontsize=18)
  
  #-- ticks ----
  plt.xticks(fontsize=18)
  plt.yticks(fontsize=18)
  
  #-- title ----
  Corr   = corrcoef(a1Map,a1Obt)[0][1]
  stitle = "%s vs %s"%(prj[0].split(".")[1], prmaptype)
  stitle = stitle + "\n%04d-%04d M%02d-%02d R=%5.2f"%(iyear,eyear,imon,emon,Corr)
  axplot.set_title(stitle, fontsize=15)
  
  #-- save -----
  figdir = "/mnt/mizu.tank/utsumi/PMM/COMP.ORB.ORG/pict"
  figname= figdir + "/%s.x.%s.png"%(".".join(prj), prmaptype)
  ctrack_func.mk_dir(figdir)
  plt.savefig(figname)
  print "SAVE"
  print figname
