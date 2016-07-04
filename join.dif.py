import Image, shutil
from numpy import *
import tag_func
import pmm_para

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
iY,eY     = 2002,2004
ny,nx     = 1200,3600
region    = "WN.PAC"
#region    = "GLOB"
para_tag  = pmm_para.Para_tag()
para_tag.thpr = 0.0
para_tag.wnflag = False
para_tag.ny   = ny
para_tag.nx   = nx
para_tag.ext  = "%s.%s"%(ny,nx)
ltag  = ["plain","tc.comb","cf.comb","ot"]
prtype1  = "GSMaP.v5"
prtype2  = "PR.2A25"
lseason  = ["JJA","DJF","SON","MAM"]
sumnum   = "amt"
#ldiftype = ["d","frac.d"] 
ldiftype = ["d","frac.d"] 

def ret_figxy(region):
  if region == "GLOB":
    iyfig = 200  # top 
    eyfig = -223 # bottom
    ixfig = 10
    exfig = -10
  elif region =="WN.PAC":
    iyfig = 180
    eyfig = -197
    ixfig = 10
    exfig = -10
  return iyfig,eyfig,ixfig,exfig

for season in lseason:
  for diftype in ldiftype:
    iyfig,eyfig,ixfig,exfig  = ret_figxy(region)
    da2dat  ={}
    for i,tag in enumerate(ltag):
      sdir_root, sdir, amtname1, amtname2, ratname1, ratname2, sumname1, sumname2, numname = ret_name_clim(prtype1, prtype2, tag, para_tag, iY, eY, season)
      
      sname   = sdir + "/pict.%s"%(region) + "/%s%s.%s.th.%s.%s.png"%(diftype, sumnum, tag, para_tag.thpr, para_tag.ext)
      a2png   = Image.open(sname)
      a2array = asarray(a2png)[iyfig:eyfig, ixfig:exfig]
      da2dat[i] = a2array
    #-----------
    a2oarray  = vstack([da2dat[0], da2dat[1], da2dat[2], da2dat[3]])
    oimg      = Image.fromarray(uint8(a2oarray))
    #-----------

    oname     = sdir + "/pict.%s"%(region) + "/join.%s%s.th.%s.%s.png"%(diftype, sumnum, para_tag.thpr, para_tag.ext)
    oimg.save(oname)
    print oname

