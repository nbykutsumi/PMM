from numpy import *
import pmm_para
import ctrack_func, tag_func
import ctrack_para
from matplotlib import pyplot as plt


iY, eY  = 2002,2009
#iY, eY  = 2002,2003
lY      = range(iY,eY+1)
season  = "ALL"
lM      = ctrack_para.ret_lmon(season)
ltag    = ["plain","tc.comb","cf.comb","ot"]
var     = "rat"
lregion = ["EAS","SAS.SEA","JPN"]
#region  = "WN.PAC"
#region  = "JPN"
prtype1 = "GSMaP.v5"
prtype2 = "PR.2A25"
ny,nx     = pmm_para.ret_nynx(prtype1)
para_tag  = pmm_para.Para_tag()
para_tag.thpr = 0.0
para_tag.wnflag = False
para_tag.ny   = ny
para_tag.nx   = nx
para_tag.ext  = "%s.%s"%(ny,nx)
miss   = -9999.

#*****************************
def mk_plot_III(da1v, da1std, ltag, region):
  figplot  = plt.figure(figsize=(6,3))
  axplot   = figplot.add_axes([0.1, 0.15, 0.87, 0.62])
  #**** sdate **********
  #lsdate    = range(1,12+1)
  lsdate    = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
  lx  =  arange(len(lsdate))
  dcolor = ret_dcolor()
  #*********************
  #------------
  for i,tag in enumerate(ltag):
    lv   = da1v[tag]
    lstd = da1std[tag]
    if tag =="plain":
      continue
    axplot.bar(lx-0.5+0.2*i, lv, width=0.2, color=dcolor[tag]) 
    axplot.errorbar(lx-0.5 +0.2*i +0.1, lv, lstd, color=dcolor[tag], linewidth=0, elinewidth=1.0, ecolor="k")
  #**** plain ***********
  axplot.plot(lx, lv, "-o", color=dcolor["plain"])
  axplot.errorbar(lx, lv, lstd, linewidth=0., elinewidth=2, ecolor=dcolor["plain"])
  #**** zero-line *******
  axplot.plot([-1,len(lM)+2],[0,0],"-",color="k",linewidth=1.0)

  #**** x-ticks *********
  axplot.xaxis.set_ticks(lx[0::2], minor=True)
  axplot.xaxis.set_ticks(lx[0::2])
  axplot.xaxis.set_ticklabels( lsdate[0::2], minor=False, fontsize=22)
  #**** y-ticks *********
  for tick in axplot.yaxis.get_major_ticks():
    tick.label.set_fontsize(18)

  #**** axis limit ******
  #axplot.set_ylim(bottom=0)
  ltmp = []
  for tag in ltag:
    ltmp =  ltmp+da1v[tag]
  ylim_max     = max(ltmp) *1.2
  axplot.set_ylim(top=ylim_max)
  #-----------------------
  axplot.set_xlim([-0.5,11.5])

  #**** grid ************
  #plt.grid(b=True, which="major", linestyle="-",  linewidth=0.3)
  #plt.grid(b=True, which="minor", linestyle="--", linewidth=0.3)

  #**** title ***********
  #ssuptitle  = "  %s"%(regionname)
  #plt.suptitle(ssuptitle,fontsize=23)

  stitle = "d.rat [mm/h] %s %s minus %s %04d-%04d"%(region, prtype1, prtype2, iY, eY)
  plt.title(stitle, fontsize=10)

  #*********************
  season   = "ALL"

  sdir_root, sdir, amtname1, amtname2, ratname1, ratname2, sumname1, sumname2, numname = ret_name_clim(prtype1, prtype2, tag, para_tag, iY, eY, season)

  figdir  = sdir + "/pict.PMM"
  figname = figdir + "/drat.season.%s-%s.%s.png"%(prtype1,prtype2,region)

  ctrack_func.mk_dir(figdir)
  figplot.savefig(figname)
  print figname

#*****************************
def ret_dcolor():
  dcolor    = {"plain":"k","tc.comb":"r", "cf.comb":"royalblue", "ot":"orange"}

  #--------
  return dcolor

#------------------------------------------
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
    print "idat",idat
    if lv[idat] >=0:
      bottom = bottom_posi
    elif lv[idat] <0:
      bottom = bottom_nega
  return bottom


#------------------------------------------
#------------------------------------------
def ret_name_clim(prtype1,prtype2,tag,para_tag,iY,eY,season):
  lname = tag_func.ret_name_tagthpr_match_sumnum_clim\
          (\
           prtype1, prtype2, tag, para_tag, iY, eY, season\
          )

  sdir_root, sdir, amtname1, amtname2, ratname1, ratname2, sumname1, sumname2, numname\
        = lname
  return sdir_root, sdir, amtname1, amtname2, ratname1, ratname2, sumname1, sumname2, numname
#------------------------------------------
def ret_la2dat(var, Y,M,prtype1,prtype2, tag):
  sdir_root, sdir, varpath1, varpath2\
      =tag_func.ret_name_tagthpr_match_mon\
      (\
         var, prtype1, prtype2, tag, para_tag, Y, M\
      )
  a2var1 = fromfile(varpath1,  float32).reshape(ny,nx)
  a2var2 = fromfile(varpath2,  float32).reshape(ny,nx)
  return a2var1, a2var2

#***********************************
da1diff = {}
da1sig  = {}
for tag in ltag:
  for region in lregion:
    da1diff[region,tag] = []
    da1sig [region,tag] = []

da1v   = {}

for Y,M in [[Y,M] for Y in lY for M in lM]:
  for tag in ltag:
    a2var1, a2var2 = ret_la2dat(var,Y,M,prtype1,prtype2,tag)
    a2dvar  = a2var1 - a2var2
    a2mask  = ma.masked_where(a2var1==miss, zeros([ny,nx],float32)).filled(miss)
    a2mask  = ma.masked_where(a2var1==miss, a2mask).filled(miss)

    #--- extract region ---
    for region in lregion:
      lllon, urlon, lllat, urlat, iyfig, eyfig, ixfig, exfig\
               = pmm_para.ret_region_para(region)
  
      a2dvar_reg = a2dvar[iyfig:eyfig+1, ixfig:exfig+1]  # 60S-60N -> region
      a2mask_reg = a2mask[iyfig:eyfig+1, ixfig:exfig+1]  # 60S-60N -> region

      dvar = ma.masked_where(a2mask_reg==miss, a2dvar_reg).mean()
      if da1v.has_key((tag,region,M))==True:
        da1v[tag,region,M].append(dvar)
      else:
        da1v[tag,region,M] = [dvar]
      #----
for region in lregion:
  #----------------------
  dlv = {}
  dlstd = {}
  for tag in ltag:
    for M in lM:
      V = mean(da1v[tag,region,M])*60*60.
      S = std(da1v[tag,region,M])*60*60.
      if dlv.has_key(tag):
        dlv  [tag].append(V)
        dlstd[tag].append(S)
      else:
        dlv  [tag]=[V]
        dlstd[tag]=[S]

  mk_plot_III(dlv, dlstd, ltag, region) 
  
  print iY,"-",eY

