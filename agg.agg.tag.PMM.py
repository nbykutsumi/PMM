from numpy import *
from tag_func import *
import ctrack_func, tag_func, pmm_func
import ctrack_para, pmm_para
import ctrack_fig

#calcflag = True
calcflag = False
#figflag  = True
figflag  = False
#lseason = ["DJF","JJA"]
#lseason = ["SON"]
#lseason = ["MAM","DJF","JJA","SON","ALL"]
#lseason = [1,2,3,4,5,6,7,8,9,10,11,12]
lseason = [1]
#iY, eY = 2002, 2009
iY, eY = 2002, 2002
lY     = range(iY, eY+1)

prtype  = "GSMaP.v5"
miss    = -9999.

para_tag = pmm_para.Para_tag()
ny,nx    = pmm_para.ret_nynx(prtype)

para_tag.thpr = 0.0
para_tag.ny   = ny
para_tag.nx   = nx
para_tag.wnflag = False
if para_tag.wnflag == False:
  #----------
  ltag      = ["plain","tc","c","fbc","te","tf","ef","tef","ot"]
#  ltag      = ["plain"]
  ltag_comb = ["tc.comb","cf.comb"]

def ret_la2dat(year,mon,prtype, tag, para_tag):
  thpr   = para_tag.thpr
  thdura_tc = para_tag.thdura_tc
  thdura_c  = para_tag.thdura_c
  dist_tc   = para_tag.dist_tc
  dist_c    = para_tag.dist_c
  dist_f    = para_tag.dist_f
  ny        = para_tag.ny
  nx        = para_tag.nx
  ext       = "%sx%s"%(ny,nx)
  sresol    = para_tag.sresol
  bstflag_tc= para_tag.bstflag_tc
  bstflag_f = para_tag.bstflag_f
  loname = tag_func.ret_name_tagthpr_sumnum_mon\
           (\
            prtype, tag, thpr, sresol, thdura_tc, thdura_c, dist_tc, dist_c, dist_f, bstflag_tc, bstflag_f, year, mon, tstep="1hr",ext=ext\
           )
  odir    = loname[1]
  sumname = loname[2]
  numname = loname[3]
  print sumname
  a2sum = fromfile(sumname, float32).reshape(ny,nx) 
  a2num  = fromfile(numname,  float32).reshape(ny,nx) 
  return a2sum, a2num

def ret_name_clim(prtype,tag,para_tag,iyear,eyear,season):
  thpr   = para_tag.thpr
  thdura_tc = para_tag.thdura_tc
  thdura_c  = para_tag.thdura_c
  dist_tc   = para_tag.dist_tc
  dist_c    = para_tag.dist_c
  dist_f    = para_tag.dist_f
  ny        = para_tag.ny
  nx        = para_tag.nx
  ext       = "%sx%s"%(ny,nx)
  sresol    = para_tag.sresol
  bstflag_tc= para_tag.bstflag_tc
  bstflag_f = para_tag.bstflag_f

  lname = tag_func.ret_name_tagthpr_sumnum_clim\
          (\
           prtype, tag, thpr, sresol, thdura_tc, thdura_c, dist_tc, dist_c, dist_f, bstflag_tc, bstflag_f, iyear, eyear, season, tstep="1hr", ext=ext\
          )

  sdir_root, sdir, sumname, numname, prname\
        = lname

  ratname = sdir + "/rat."+ prname.split("/")[-1][2:]

  return sdir_root, sdir, sumname, numname, prname, ratname

def ret_a2clim(tag,sumnum):
  sdir_root, sdir, sumname, numname, prname, ratname = ret_name_clim(prtype,tag,para_tag,iY,eY,season)

  if sumnum=="pr":
    a2dat  = fromfile(prname, float32).reshape(ny,nx)

  elif sumnum == "rat":
    a2rat  = (ma.masked_where(a2num_tag==0., a2sum)/a2num_tag).filled(miss)

  elif sumnum=="sum":
    a2dat  = fromfile(sumname, float32).reshape(ny,nx)

  elif sumnum=="num":
    a2dat  = fromfile(numname, float32).reshape(ny,nx)

  return a2dat

def mk_la2comb(tag_comb):
  if tag_comb == "tc.comb":
    a2sum_tc  = ret_a2clim("tc","sum")
    a2sum_te  = ret_a2clim("te","sum")
    a2sum_tf  = ret_a2clim("tf","sum")
    a2sum_tef = ret_a2clim("tef","sum")

    a2num_tc  = ret_a2clim("tc","num")
    a2num_te  = ret_a2clim("te","num")
    a2num_tf  = ret_a2clim("tf","num")
    a2num_tef = ret_a2clim("tef","num")
    a2num_pln = ret_a2clim("plain","num")

    print "tc"
    print a2sum_tc
    print "te"
    print a2sum_te
    print "tf"
    print a2sum_tf
    print "tef"
    print a2sum_tef 
    a2sum = a2sum_tc + a2sum_te*0.5 + a2sum_tf*0.5 + a2sum_tef*(1./3.)

    a2num_tag  = a2num_tc + a2num_te + a2num_tf + a2num_tef


    a2pr    = (ma.masked_where(a2num_pln==0., a2sum)/a2num_pln).filled(miss)

    a2rat  = (ma.masked_where(a2num_tag==0., a2sum)/a2num_tag).filled(miss)

  elif tag_comb == "cf.comb":
    a2sum_c   = ret_a2clim("c","sum")
    a2sum_fbc = ret_a2clim("fbc","sum")
    a2sum_ef  = ret_a2clim("ef","sum")
    a2sum_te  = ret_a2clim("te","sum")
    a2sum_tf  = ret_a2clim("tf","sum")
    a2sum_tef = ret_a2clim("tef","sum")

    a2num_c   = ret_a2clim("c","num")
    a2num_fbc = ret_a2clim("fbc","num")
    a2num_ef  = ret_a2clim("ef","num")
    a2num_te  = ret_a2clim("te","num")
    a2num_tf  = ret_a2clim("tf","num")
    a2num_tef = ret_a2clim("tef","num")
    a2num_pln = ret_a2clim("plain","num")

    a2sum = a2sum_c + a2sum_fbc + a2sum_ef\
            +a2sum_te*0.5+ a2sum_tf*0.5 + a2sum_tef*(1./3.) 

    a2num_tag = a2num_c + a2num_fbc + a2num_ef\
               +a2num_te + a2num_tf + a2num_tef

    a2pr   = (ma.masked_where(a2num_pln==0., a2sum)/a2num_pln).filled(miss)

    a2rat  = (ma.masked_where(a2num_tag==0., a2sum)/a2num_tag).filled(miss)

  return a2pr, a2rat, a2sum, a2num_tag

#***************************************************
if calcflag == True:
  for season in lseason:
    lM = ctrack_para.ret_lmon(season)
    #---- num plain ---------
    tag = "plain"
    a2num_pln   = zeros([ny,nx],float32)
    for Y,M in [[Y,M] for Y in lY for M in lM]:
      a2num_pln_tmp = ret_la2dat(Y,M,prtype,"plain",para_tag)[1]
      a2num_pln     = a2num_pln + a2num_pln_tmp
    #---- calc precip (mean amount and rate) --
    for tag in ltag:
      a2sum  = zeros([ny,nx],float32)
      a2sum2  = zeros([ny,nx],float32)
      a2num_tag   = zeros([ny,nx],float32)
      for Y,M in [[Y,M] for Y in lY for M in lM]:
  #      print Y,M,tag
        a2sum_tmp, a2num_tmp = ret_la2dat(Y,M,prtype,tag,para_tag)
        a2sum = a2sum + a2sum_tmp
        a2num_tag  = a2num_tag  + a2num_tmp
  
        print Y,M,tag, a2sum_tmp.mean()
      #-------------
      a2pr  = (ma.masked_where(a2num_pln==0., a2sum)/a2num_pln).filled(miss)
      a2rat = (ma.masked_where(a2num_tag==0., a2sum)/a2num_tag).filled(miss)
  
      #-- write ----
      sdir_root, sdir, sumname, numname, prname, ratname\
            = ret_name_clim\
              (\
               prtype,tag,para_tag,iY,eY,season\
              )
  
      ctrack_func.mk_dir(sdir)
      a2pr.tofile(prname)
      a2rat.tofile(ratname)
      a2sum.tofile(sumname)
      a2num_tag.tofile(numname)
      print numname
  
    #-- tag.comb -----
    if len(ltag_comb) ==0:
      continue
    for tag_comb in ltag_comb:
      a2pr, a2rat, a2sum, a2num_tag = mk_la2comb(tag_comb)
  
      sdir_root, sdir, sumname, numname, prname, ratname\
      = ret_name_clim\
        (\
         prtype, tag_comb, para_tag, iY, eY, season\
        )

      a2pr.tofile(prname)
      a2rat.tofile(ratname)
      a2sum.tofile(sumname)
      a2num_tag.tofile(numname)
      print numname
 
#--------- Figures -----------
if figflag == True:
  for season in lseason:
    for tag in ["plain","ot"] + ltag_comb:
      #for sumnum in ["amt","rat"]:
      for sumnum in ["pr","rat"]:
        #-- names --
        sdir_root, sdir, sumname, numname, prname, ratname\
           = ret_name_clim\
             (\
              prtype,tag,para_tag,iY,eY,season\
             )
  
        if sumnum == "pr":
          iname = prname 
        elif sumnum == "rat":
          iname = ratname
        elif sumnum == "num":
          iname = numname
        elif sumnum == "sum":
          iname = sumname

 
        figdir  = sdir + "/pict"
        ctrack_func.mk_dir(figdir)
        figname = figdir + "/" +  iname.split("/")[-1] + ".png"
        cbarname= figdir + "/" +  "cbar.%s.png"%(sumnum)
        #-- load ---
        a2dat = fromfile(iname, float32).reshape(ny,nx)
    
        a2dat = a2dat[200:1000]  # 60S-60N -> 40S-40N
        #-- mask ---
        a2shade= a2dat
        #-- lat lon --
        lllon  = 0.05
        urlon  = 359.95
        lllat  = -39.95
        urlat  = 39.95

        if sumnum == "pr":
          coef = 60.*60.*24*30.
          bnd  = [0,10,30,50,70,100,200,300] 
        elif sumnum == "rat":
          coef = 60.*60.
          bnd  = [0.0,0.1,0.2,0.3,0.5,0.7,1,2] 
        elif sumnum == "num":
          coef = 1.
          bnd  = [1,10,30,50,70,100,200,300] 

        #---- mycm  ----
        mycm   = "rainbow"                  
        #---- title ----
        stitle = "%s %s %04d-%04d %s"%(prtype, tag, iY, eY, season) 
        ctrack_fig.mk_map(a2dat,bnd=bnd, mycm="jet",soname=figname, stitle=stitle, cbarname=cbarname, miss=miss, lllat=lllat, lllon=lllon, urlat=urlat, urlon=urlon, a2shade=a2shade, coef=coef)
        print figname
#  

##---- test: check data: clim ----------
#season = 1
#a2tot = zeros([ny,nx],float32)
#
#for tag in ["tc.comb","cf.comb","ot","plain"]:
##for tag in ["plain","tc","c","fbc","te","tf","ef","tef","ot"]:
##for tag in ltag:
#  a2dat = ret_a2clim(tag,"pr")
#  a2dat  = ma.masked_equal(a2dat, miss).filled(0.0)
#
#  print tag,a2dat.sum()
#  if tag !="plain":
#    a2tot = a2tot + a2dat
#  else:
#    a2plain = a2dat
#
#a2d  = a2tot - a2plain
#
#print "-- a2d.mean, a2d.min, a2d,max --"
#print a2d.min(), a2d.mean(), a2d.max()



#---- test: check data: mon -----------
Y = 2002
#lM   = [6,7,8]
lM   = [1]
for M in lM:
  da2sum   = {}
  da2num   = {}
  a2tot    = zeros([ny,nx],float32)
  a2totnum = zeros([ny,nx],float32)
  
  for tag in ltag:
    a2sum, a2num = ret_la2dat(Y,M,prtype, tag, para_tag)
    da2sum[tag]  = a2sum
    da2num[tag]  = a2num
    if tag != "plain":
      a2tot = a2tot + a2sum
      a2totnum  = a2totnum  + a2num
 
    elif tag == "plain":
      a2plain = a2sum
      a2plainnum = a2num
  
  a2d  = a2tot - a2plain
  a2dnum= a2totnum - a2plainnum
  print "---------------"
  print "Mon=",M
  print "d.mean(), d.min(), dmax()"
  print a2d.mean(), a2d.min(), a2d.max()
  print a2dnum.mean(), a2dnum.min(), a2dnum.max() 
 
