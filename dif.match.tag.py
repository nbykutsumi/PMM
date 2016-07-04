from numpy import *
import pmm_func, ctrack_func, tag_func
import pmm_para, ctrack_para
import ctrack_fig
import sys

iY,eY   = 2002, 2004
lseason = ["DJF","JJA","SON","MAM"]
#lseason = ["SON"]
lregion = ["WN.PAC"]

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
    urlat  = 39
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
#------------------------------------------
prtype1   = "GSMaP.v5"
prtype2   = "PR.2A25"
ny,nx     = pmm_para.ret_nynx(prtype1)
para_tag  = pmm_para.Para_tag()
para_tag.thpr = 0.0
para_tag.wnflag = False
para_tag.ny   = ny
para_tag.nx   = nx
para_tag.ext  = "%s.%s"%(ny,nx)
ltag  = ["plain","tc.comb","cf.comb","ot"]
miss  = -9999.
sumnum= "amt"

for season in lseason:
  da2d     = {}
  da2fracd = {}
  for tag in ltag:
    a2pr1, a2pr2  = ret_la2clim(tag,sumnum)
    da2d[tag]     = a2pr1 - a2pr2 
    da2fracd[tag] = (ma.masked_where(a2pr2 <= 0.0, (a2pr1-a2pr2))/a2pr2).filled(miss)
    #-- write -----
    sdir_root, sdir, amtname1, amtname2, ratname1, ratname2, sumname1, sumname2, numname = ret_name_clim(prtype1, prtype2, tag, para_tag, iY, eY, season)
  
    
    sname      = sdir + "/d%s.%s.th.%s.%s"%(sumnum, tag, para_tag.thpr, para_tag.ext)
    sfracname  = sdir + "/frac.d%s.%s.th.%s.%s"%(sumnum, tag, para_tag.thpr, para_tag.ext)
    da2d[tag].tofile(sname)
    da2fracd[tag].tofile(sfracname)
    print sname
    #*************************************
    # Draw
    #*************************************
    for region in lregion:
      #--- fig : difference ---
      lllon, urlon, lllat, urlat, iyfig, eyfig, ixfig, exfig\
               = ret_region_para(region)
      a2dat    = da2d[tag][iyfig:eyfig+1, ixfig:exfig+1]  # 60S-60N -> 40S-40N
      a2shade  = a2pr1[iyfig:eyfig, ixfig:exfig]
    
      if sumnum == "amt":
        coef = 60.*60.*24.*30.
        bnd  = [-50,-30,-10,-5,0,5,10,30,50]
      elif sumnum == "rat":
        coef = 60.*60.
        bnd  = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1, 1.5, 2]   
    
      #---- title ----
      stitle = "%s minus %s %s %04d-%04d %s"%(prtype1, prtype2, tag, iY, eY, season)
      #---- names ----
      figdir = sdir   + "/pict.%s"%(region)
  
      figname= figdir + "/" +sname.split("/")[-1] + ".png"
      ctrack_func.mk_dir(figdir)
      cbarname = figdir + "/cbar.d%s.png"%(sumnum)
      #----- color ---
      mycm   = "RdYlBu_r"
      #----------------
      ctrack_fig.mk_map_symm(a2dat,bnd=bnd, mycm=mycm,soname=figname, stitle=stitle, cbarname=cbarname, miss=miss, lllat=lllat, lllon=lllon, urlat=urlat, urlon=urlon, a2shade=a2shade, coef=coef)
      print figname
    
      #*************************************
      #--- fig : fractional difference ---
      lllon, urlon, lllat, urlat, iyfig, eyfig, ixfig, exfig\
               = ret_region_para(region)
      a2dat    = da2fracd[tag][iyfig:eyfig+1,ixfig:exfig+1]  # 60S-60N -> 40S-40N
      a2shade  = a2dat

      print a2dat
    
      coef = 1.0
      bnd  = [-1.0, -0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8, 1.0] 
      #bnd=False
      #---- title ----
      stitle = "frac.d %s %s %s minus %s %04d-%04d"%(tag, season, prtype1, prtype2, iY, eY)
      #---- names ----
      figdir = sdir   + "/pict.%s"%(region)
      figname= figdir + "/" +sfracname.split("/")[-1] + ".png"
      ctrack_func.mk_dir(figdir)
      cbarname = figdir + "/cbar.frac.d%s.png"%(sumnum)
      #----- color ---
      mycm   = "RdYlBu_r"
      #----------------
      ctrack_fig.mk_map_symm(a2dat,bnd=bnd, mycm=mycm,soname=figname, stitle=stitle, cbarname=cbarname, miss=miss, lllat=lllat, lllon=lllon, urlat=urlat, urlon=urlon, a2shade=a2shade, coef=coef)
      print figname
  
  
