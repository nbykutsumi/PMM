from numpy import *
import ioPMM
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import myfunc.util as util
import sys

iYM = [2014,4]
eYM = [2014,12]

lYM = util.ret_lYM(iYM, eYM)
ltag= ["plain","c","fbc","tc","ot"]
#ltag= ["fbc"]
nclass=4

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

#lvtype = ["Bias.Rte","Bias.Prc","R.Bias.Rte","TR.Bias.Prc","Corr"]
#lvtype = ["R.Bias.Rte","TR.Bias.Prc","Corr"]
#lvtype = ["Bias.Prc","Bias.Rte"]
lvtype = ["Bias.Prc"]

dUnit  = {"Bias.Prc":"mm/hour", "Bias.Rte":"mm/hour", "R.Bias.Rte":"-","TR.Bias.Prc":"-"
         ,"Corr":"-"
         }

def calVal(vtype,tag,Year,Mon):
  if vtype   == "Bias.Prc":
    a2pr1 = match.ret_monPrc(prtype1, tag, Year, Mon)
    a2pr2 = match.ret_monPrc(prtype2, tag, Year, Mon)

    if a2pr1.mean !=0.0:
      val = a2pr2.mean() - a2pr1.mean()
    else:
      val = None 

  elif vtype == "Bias.Rte":
    a2pr1 = match.ret_monRte(prtype1, tag, Year, Mon)
    a2pr2 = match.ret_monRte(prtype2, tag, Year, Mon)

    if a2pr1.mean !=0.0:
      val = a2pr2.mean() - a2pr1.mean()
    else:
      val = None 

  elif vtype == "R.Bias.Rte":
    a2pr1 = match.ret_monRte(prtype1, tag, Year, Mon)
    a2pr2 = match.ret_monRte(prtype2, tag, Year, Mon)
    if (a2pr1.mean() !=0.0):
      val   = (a2pr2.mean() - a2pr1.mean())/a2pr1.mean()
    else:
      val = None

  elif vtype == "TR.Bias.Prc":
    a2plain = match.ret_monPrc(prtype1, "plain", Year, Mon)
    a2sum1  = match.ret_monPrc(prtype1, tag, Year, Mon)
    a2sum2  = match.ret_monPrc(prtype2, tag, Year, Mon)
    if (a2plain.mean() !=0.0):
      val = (a2sum2.mean() - a2sum1.mean())/a2plain.mean()
    else:
      val = None

  elif vtype == "Corr":
    a2pr1 = match.ret_monPrc(prtype1, tag, Year, Mon)
    a2pr2 = match.ret_monPrc(prtype2, tag, Year, Mon)
    val = corrcoef(a2pr2.flatten(), a2pr1.flatten())[0][1]
  else:
    print "check vtype:",vtype
    sys.exit()
  return val
#***** main ************************   

match  = ioPMM.matchRAdom()
ny, nx = match.ny, match.nx
for vtype in lvtype:
  for tag in ltag:
    dv = {}
    for prtype2 in lprtype2:
      Trc         = dTrc[prtype2]
      match(prtype1=prtype1, prtype2=prtype2, Trc=Trc, nclass=nclass)
      dv[prtype2] = []
      for Year,Mon in lYM:
        print Year, Mon
        #a2pr1 = match.ret_monAve(prtype1, tag, Year, Mon)
        #a2pr2 = match.ret_monAve(prtype2, tag, Year, Mon)
        v  = calVal(vtype,tag,Year,Mon)
        dv[prtype2].append(v)

    #- plot -------
    figplot = plt.figure(figsize=(4,3))
    axplot  = figplot.add_axes([0.2, 0.1, 0.75, 0.75])
    for prtype2 in lprtype2:
      axplot.plot( dv[prtype2], "-")
  
    #- axis label --
    #axplot.set_ylabel()

    #- axis limit --
    if   vtype=="Bias.Prc":
      axplot.set_ylim(-0.1,0.1) 
    elif vtype=="TR.Bias.Prc":
      axplot.set_ylim(-0.5,0.5) 
    elif vtype=="R.Bias.Rte":
      axplot.set_ylim(-1.0,1.0) 
 
    #- ticks -------
    axplot.xaxis.set_ticklabels([dMonName[Mon] for Mon in zip(*lYM)[1]])
  
    #- grid lines --
    axplot.grid(b=True, which="major", color="k",linestyle="-")
    #- zero line  --
    axplot.plot(zeros(len(lYM)),color="k",linestyle="-",linewidth=2)
  
    #-- axis label --
    #axplot.set_ylabel("[mm/hour]", fontsize=16)
  
    #-- title -------
    stitle  = "%s"%(dTagName[tag])
    stitle  = stitle + "\n%s [%s]"%(vtype,dUnit[vtype])
    axplot.set_title(stitle)
  
    #- save --------
    figDir  = "/tank/utsumi/PMM/RAdom/pict/%dclass"%(nclass)
    figPath = figDir + "/TS.%s.%s.png"%(vtype, tag)
    util.mk_dir(figDir)
    plt.savefig(figPath)
    print figPath 
  
  #- legend ------
  lines     = axplot.get_lines()
  figlegend = plt.figure(figsize=(3,3))
  figlegend.legend(lines, lprtype2)
  params = {'legend.fontsize': 20,
            'legend.linewidth':5}
  plt.rcParams.update(params)
  figlegend.savefig(figDir+"/legend.TS.plot.png")

