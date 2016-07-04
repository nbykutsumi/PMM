from numpy import *
import pmm_para
import ctrack_func, tag_func
from matplotlib import pyplot as plt


iY, eY  = 2002,2009
#iY, eY  = 2002,2002
#lseason = ["MAM","JJA","SON","DJF"]
lseason = [1,2,3,4,5,6,7,8,9,10,11,12]
#lseason = [1]
ltag    = ["plain","tc.comb","cf.comb","ot"]
sumnum  = "amt"
#region  = "WN.PAC"
region  = "JPN"
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
def mk_plot_III(da1sv, da1std, ltag, region, diftype):
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
    #-------
    if tag =="plain":
      continue
    #-------
    axplot.bar(lx-0.5+0.2*i, da1sv[tag], width=0.2, color=dcolor[tag]) 

  #**** plot plain ******
  axplot.plot(da1sv[tag], "-o",color="k")
  #**** zero-line *******
  axplot.plot([-1,len(lseason)+2],[0,0],"-",color="k",linewidth=1.0)

  #**** x-ticks *********
  axplot.xaxis.set_ticks(lx[0::2], minor=True)
  axplot.xaxis.set_ticks(lx[0::2])
  axplot.xaxis.set_ticklabels( lsdate[0::2], minor=False, fontsize=22)
  #**** y-ticks *********
  for tick in axplot.yaxis.get_major_ticks():
    tick.label.set_fontsize(18)

  #**** axis limit ******
  #axplot.set_ylim(bottom=0)
  ylim_max     = 100
  axplot.set_ylim(top=ylim_max)
  #-----------------------
  axplot.set_xlim([-0.5,11.5])

  #**** grid ************
  #plt.grid(b=True, which="major", linestyle="-",  linewidth=0.3)
  #plt.grid(b=True, which="minor", linestyle="--", linewidth=0.3)

  #**** title ***********
  #ssuptitle  = "  %s"%(regionname)
  #plt.suptitle(ssuptitle,fontsize=23)

  stitle = "%spr %s %s minus %s %04d-%04d"%(diftype, tag, prtype1, prtype2, iY, eY)
  plt.title(stitle, fontsize=10)

  #*********************
  season   = "ALL"

  figname  = "./test2.png" 
  #ctrack_func.mk_dir(figdir)
  figplot.savefig(figname)
  print figname

#*****************************
def ret_dcolor():
  dcolor    = {"plain":"gray","tc.comb":"r", "cf.comb":"royalblue", "ot":"orange"}

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
def ret_region_para(region):
  if region == "GLOB":
    lllon  = 0.05
    urlon  = 359.95
    lllat  = -39.95
    urlat  = 39.95
    iyfig  = 200    # 60S-60N -> 40S-40N
    eyfig  = 1000
    ixfig  = 0
    exfig  = None

  elif region == "WN.PAC":
    lllon  = 100
    urlon  = 180
    lllat  = 15
    urlat  = 40
    iyfig  = (lllat + 60)*10  # 15N-40N
    eyfig  = (urlat + 60)*10
    ixfig  = (lllon)*10
    exfig  = (urlon)*10
  elif region == "JPN":
    lllon  = 125
    urlon  = 135
    lllat  = 30
    urlat  = 40
    iyfig  = (lllat + 60)*10  # 15N-40N
    eyfig  = (urlat + 60)*10
    ixfig  = (lllon)*10
    exfig  = (urlon)*10


  else:
    print "check region!"
    sys.exit()
  return lllon, urlon, lllat, urlat, iyfig, eyfig, ixfig, exfig

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
def ret_la2clim(tag,sumnum):
  sdir_root, sdir, amtname1, amtname2, ratname1, ratname2, sumname1, sumname2, numname\
     =(\
       ret_name_clim(prtype1,prtype2,tag,para_tag,iY,eY,season)\
      )
  if sumnum=="amt":
    a2dat1  = fromfile(amtname1, float32).reshape(ny,nx)
    a2dat2  = fromfile(amtname2, float32).reshape(ny,nx)
    return a2dat1, a2dat2

  if sumnum=="rat":
    a2dat1  = fromfile(ratname1, float32).reshape(ny,nx)
    a2dat2  = fromfile(ratname2, float32).reshape(ny,nx)
    return a2dat1, a2dat2

  elif sumnum=="sum":
    a2dat1  = fromfile(sumname1, float32).reshape(ny,nx)
    a2dat2  = fromfile(sumname2, float32).reshape(ny,nx)
    return a2dat1, a2dat2

  elif sumnum=="num":
    a2dat  = fromfile(numname, float32).reshape(ny,nx)
    return a2dat

#***********************************
da1diff = {}
da1frac = {}
for tag in ltag:
  da1diff[tag] = []
  da1frac[tag] = []

da2pr1 = {}
da2pr2 = {}
for season in lseason:
  da2diff = {}
  da2frac = {}
  for tag in ltag:
    print season,tag
    #--- extract region ---
    lllon, urlon, lllat, urlat, iyfig, eyfig, ixfig, exfig\
             = ret_region_para(region)

    a2pr1, a2pr2 = ret_la2clim(tag,sumnum)
    a2pr1_reg = a2pr1[iyfig:eyfig+1, ixfig:exfig+1]  # 60S-60N -> region
    a2pr2_reg = a2pr2[iyfig:eyfig+1, ixfig:exfig+1]  # 60S-60N -> region

    dfrac = ma.masked_where(a2pr2_reg==miss, (a2pr1_reg - a2pr2_reg)).mean() / ma.masked_equal(a2pr2_reg, miss).mean()
    #----
    da1diff[tag].append(dfrac*100.)

#--- error bar --------
da1diff_std  = {}
for tag in ltag:
  da1diff_std[tag] = [0.]*len(lseason)
#----------------------
mk_plot_III(da1diff, da1diff_std, ltag, region,"d") 

print iY,"-",eY

#    sdir_root, sdir, amtname1, amtname2, ratname1, ratname2, sumname1, sumname2, numname = ret_name_clim(prtype1, prtype2, tag, para_tag, iY, eY, season)
#
#
#    figdir     = sdir + "/pict.%s"%(region)
#    sfracname  = sdir + "/frac.d%s.%s.th.%s.%s"%(sumnum, tag, para_tag.thpr, para_tag.ext)
#
