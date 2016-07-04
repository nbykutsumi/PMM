from numpy import *
import ioPMM
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import myfunc.util as util

iYM = [2014,4]
eYM = [2014,12]

lYM = util.ret_lYM(iYM, eYM)
ltag= ["plain","c","fbc","tc","ot"]
#ltag= ["ot"]
nclass=4

prtype1 = "RA"
#prtype1 = "GSMaP"
lprtype2 = ["RA","GPM.KuPR"]  # <-- spatial constrain
#lprtype2 = ["GPM.KuPR"]  # <-- spatial constrain

Trc     = "ALL"

dMonName={1:"Jan",2:"Feb",3:"Mar",4 :"Apr",5 :"May",6 :"Jun"
         ,7:"Jul",8:"Aug",9:"Sep",10:"Oct",11:"Nov",12:"Dec"}

dTagName={"tc":"tropical cyclone", "c":"extratropical cyclone (center)"
         ,"fbc":"front"          ,"ot":"others"
         ,"plain": "total"
         }

#* Functions ***********************
def calVal( prtype, tag,Year,Mon):
  #-- calc monthly mean precip. --
  a2sum   = match.ret_monSum(prtype, tag, Year, Mon)
  a2numPln= match.ret_monNum("plain", Year, Mon)
  a2pr    = ma.masked_where(a2numPln==0.0, a2sum)/a2numPln
  return a2pr.mean()

#***** main ************************   

match  = ioPMM.matchRAdom()
ny, nx = match.ny, match.nx
for tag in ltag:
  dv = {}
  for prtype2 in lprtype2:
    dv[prtype2]  = []
    match(prtype1=prtype1, prtype2=prtype2, Trc=Trc, nclass=nclass)

    for Year,Mon in lYM:
      print Year, Mon
      v  = calVal(prtype1,tag,Year,Mon)
      dv[prtype2].append(v)
  #- plot -------
  figplot = plt.figure(figsize=(4,2))
  axplot  = figplot.add_axes([0.2, 0.1, 0.75, 0.68])
  lx      = range(len(lYM))
  for iprtype2, prtype2 in enumerate(lprtype2):
    if   iprtype2 == 0:
      linestyle = "-"
    elif iprtype2 == 1:
      linestyle = "--"
    axplot.plot( dv[prtype2], linestyle=linestyle,color="k",linewidth=3)
  
  #- axis label --
  #axplot.set_ylabel()

  #- axis limit --
  axplot.set_ylim(0.0,0.4) 
  axplot.set_xlim(lx[0]-0.5, lx[-1]+0.5) 
  #- ticks -------
  axplot.xaxis.set_ticks(lx)
  axplot.xaxis.set_ticklabels([dMonName[Mon] for Mon in zip(*lYM)[1]])
  
  #- grid lines --
  axplot.grid(b=True, which="major", color="k",linestyle="-")
  #- zero line  --
  axplot.plot(zeros(len(lYM)),color="k",linestyle="-",linewidth=2)
  
  #-- axis label --
  #axplot.set_ylabel("[mm/hour]", fontsize=16)
  
  #-- title -------
  stitle  = "%s"%(dTagName[tag])
  stitle  = stitle + "\nmonthly mean %s [mm/hour]"%(prtype1)
  axplot.set_title(stitle)
  
  #- save --------
  figDir  = "/tank/utsumi/PMM/RAdom/pict/%dclass"%(nclass)
  figPath = figDir + "/TS.Prcp.%s.%s.png"%(prtype1, tag)
  util.mk_dir(figDir)
  plt.savefig(figPath)
  print figPath 
  
#- legend ------
lines     = axplot.get_lines()
figlegend = plt.figure(figsize=(3,3))
figlegend.legend(lines, lprtype2)
params = {'legend.fontsize': 50,
          'legend.linewidth':30}
plt.rcParams.update(params)
plt.title("spatial constraint")
figlegend.savefig(figDir+"/legend.TS.Prcp.png")

