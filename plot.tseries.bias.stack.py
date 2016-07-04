from numpy import *
import ioPMM
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import myfunc.util as util
import sys

iYM = [2014,4]
eYM = [2014,12]
lYM = util.ret_lYM(iYM,eYM)
ltag= ["c","fbc","tc","ot"]
nclass  =4

prtype1 = "RA"
#lprtype2= ["GPM.KuPR","GPM.KaPR","GSMaP"]
lprtype2= ["GPM.KuPR","GSMaP"]
#lprtype2= ["GPM.KuPR"]

dTrc    = {"GPM.KuPR":"ALL"
          ,"GPM.KaPR":"ALL"
          ,"GSMaP"   :"GPM.KuPR"}

dMonName={1:"Jan",2:"Feb",3:"Mar",4 :"Apr",5 :"May",6 :"Jun"
         ,7:"Jul",8:"Aug",9:"Sep",10:"Oct",11:"Nov",12:"Dec"}

dTagName={"tc":"tropical cyclone", "c":"extratropical cyclone (center)"
         ,"fbc":"front"          ,"ot":"others"
         ,"plain": "total"
         }

dcolor = {"c":"royalblue","fbc":"lime","tc":"r","ot":"yellow"}
match  = ioPMM.matchRAdom()
ny, nx = match.ny, match.nx

hfig   = 6
wfig   = 4

iyfig0 = 0.05
wyfig0 = 0.3
iyfig1 = iyfig0+wyfig0+0.1
wyfig1 = 0.3
iyfig2 = iyfig1 + wyfig1 + 0.05
wyfig2 = (1.0-iyfig2)*0.8

#iyfig1 = 0.1
#wyfig1 = 0.45
#iyfig2 = iyfig1 + wyfig1 + 0.05
#wyfig2 = (1.0-iyfig2)*0.8

ixfig  = 0.15
wxfig  = 0.82
#-----------------------
def calVal(prtype1, prtype2, tag,Year,Mon):
  Trc   = dTrc[prtype2]
  match(prtype1=prtype1, prtype2=prtype2, Trc=Trc, nclass=nclass)
  a2pr1 = match.ret_monPrc(prtype1, tag, Year, Mon)
  a2pr2 = match.ret_monPrc(prtype2, tag, Year, Mon)

  if a2pr1.mean !=0.0:
    val = a2pr2.mean() - a2pr1.mean()
  else:
    val = None
  return val
#------------------------
def ret_bottom(lv, idat):
  if idat == 0:
    bottom = 0
  else:
    bottom_posi  = 0.0
    bottom_nega  = 0.0
    for i,v in enumerate(lv[:idat]):
      if v >=0:
        bottom_posi = bottom_posi + v
      elif v <0:
        bottom_nega = bottom_nega + v
    #-------
    if lv[idat] >=0:
      bottom = bottom_posi
    elif lv[idat] <0:
      bottom = bottom_nega
  return bottom

#-----------------------
for prtype2 in lprtype2:
  dlv = {}
  for tag in ltag:
    dlv[tag] = []

  for tag in ltag:
    for Year,Mon in lYM:
      v = calVal(prtype1, prtype2, tag, Year, Mon)
      dlv[tag].append(v)

  #- bottom ---
  a2v      = empty([len(ltag), len(lYM)])
  a2bottom = empty([len(ltag), len(lYM)])
  for itag, tag in enumerate(ltag):
    a2v[itag] = dlv[tag]

  for i, YM in enumerate(lYM):
    ltmp = a2v[:,i]
    for itag, tag in enumerate(ltag):
      a2bottom[itag, i] = ret_bottom(ltmp, itag)

  #-- plain bias --
  lBiasTot = []
  for i, YM in enumerate(lYM):
    Year,Mon = YM
    v = calVal(prtype1, prtype2, "plain", Year, Mon)
    lBiasTot.append(v)

  #*****************************
  # system-wise precipitation
  #-----------------------------
  a2prcp        = empty([len(ltag), len(lYM)])
  a2bottom_prcp = empty([len(ltag), len(lYM)])
  for itag, tag in enumerate(ltag):
    a2prcp[itag] = dlv[tag]


  for i, YM in enumerate(lYM):
    Year,Mon = YM
    for itag, tag in enumerate(ltag):
      prcp  = match.ret_monPrc(prtype2, tag, Year, Mon).mean()
      a2prcp[itag, i] = prcp
  #-- bottom --
  for i, YM in enumerate(lYM):
    ltmp  = a2prcp[:,i]
    for itag, tag in enumerate(ltag):
      a2bottom_prcp[itag,i] = ret_bottom(ltmp, itag)

  #*****************************
  # plain precipitation
  #-----------------------------
  lprcp_plain = []
  for i, YM in enumerate(lYM):
    Year,Mon = YM
    lprcp_plain.append(match.ret_monPrc("RA", "plain", Year, Mon).mean())
  #*****************************
  # plot
  #*****************************
  figplot = plt.figure(figsize=(wfig,hfig))

  #-----------------
  # (axplot0) System-wise precip (axplot0)
  #----------------- 
  axplot0 = figplot.add_axes([ixfig, iyfig0, wxfig, wyfig0])
  lx      = range(len(lYM))

  for itag, tag in enumerate(ltag):
    axplot0.bar(lx, a2prcp[itag], bottom=a2bottom_prcp[itag],  align="center", color=dcolor[tag])

  axplot0.plot(lx, lprcp_plain, "-", color="k", linewidth=3)
  #-----------------
  # (axplot1) Bars
  #----------------- 
  axplot1 = figplot.add_axes([ixfig, iyfig1, wxfig, wyfig1])
  lx      = range(len(lYM))
  for itag, tag in enumerate(ltag):
    axplot1.bar(lx, dlv[tag], bottom=a2bottom[itag],  align="center", color=dcolor[tag])

  #-----------------
  # (axplot2) Total bias
  #----------------- 
  axplot2 = figplot.add_axes([ixfig, iyfig2, wxfig, wyfig2])
  lx      = range(len(lYM))
  axplot2.plot(lx, lBiasTot, "-", color="k", linewidth=3)

  #-------------------------
  # common
  #-------------------------
  #- axis limit -
  xmin = lx[0] -0.5
  xmax = lx[-1]+0.5
  axplot0.set_xlim(xmin, xmax)
  axplot1.set_xlim(xmin, xmax)
  axplot2.set_xlim(xmin, xmax)

  axplot0.set_ylim(0.0,0.35)
  axplot1.set_ylim(-0.1,0.06)
  axplot2.set_ylim(-0.1,0.06)

  #- zero line  --
  axplot1.plot([xmin,xmax],[0,0],color="k",linestyle="-",linewidth=2)
  axplot2.plot([xmin,xmax],[0,0],color="k",linestyle="-",linewidth=2)

  #- grid lines --
  axplot0.yaxis.grid(b=True, which="major", color="k",linestyle="-")
  axplot1.yaxis.grid(b=True, which="major", color="k",linestyle="-")
  axplot2.grid(b=True, which="major", color="k",linestyle="-")

  #- axis ticks ---
  axplot0.xaxis.set_ticks(lx)
  axplot0.xaxis.set_ticklabels([dMonName[Mon] for Mon in zip(*lYM)[1]], fontsize=11)

  axplot1.xaxis.set_ticks(lx)
  axplot1.xaxis.set_ticklabels([dMonName[Mon] for Mon in zip(*lYM)[1]], fontsize=11)
  axplot2.xaxis.set_ticks(lx)
  axplot2.xaxis.set_ticklabels([dMonName[Mon] for Mon in zip(*lYM)[1]], fontsize=11)

  #- title ----
  stitle = "bias [mm/hour] (%s - %s)"%(prtype2,prtype1)
  axplot2.set_title(stitle, fontsize=15)

  stitle = "precipitation [mm/hour]"
  axplot0.set_title(stitle, fontsize=15)


  #- save -----
  figDir  = "/tank/utsumi/PMM/RAdom/pict/%dclass"%(nclass)
  figPath = figDir + "/TS.Stack.Bias.Prcp.%s.png"%(prtype2)
  util.mk_dir(figDir)
  plt.savefig(figPath)
  print figPath


