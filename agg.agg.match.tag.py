from numpy import *
from tag_func import *
import ctrack_func, tag_func, pmm_func
import ctrack_para, pmm_para
import ctrack_fig

calcflag = True
#calcflag = False
#figflag  = True
figflag  = False
#lseason = ["DJF","JJA"]
lseason = ["ALL"]
#lseason = ["MAM","DJF","JJA","SON","ALL"]
lseason = [1,2,3,4,5,6,7,8,9,10,11,12]
#lseason = [1]
iY, eY = 2002, 2009
#iY, eY = 2002, 2002
lY     = range(iY, eY+1)

prtype1 = "GSMaP.v5"
prtype2 = "PR.2A25"
miss    = -9999.

para_tag = pmm_para.Para_tag()
ny,nx    = pmm_para.ret_nynx(prtype1)

para_tag.thpr = 0.0
para_tag.ny   = ny
para_tag.nx   = nx
para_tag.wnflag = False
if para_tag.wnflag == False:
  #----------
  ltag      = ["plain","tc","c","fbc","te","tf","ef","tef","ot","tc.comb","cf.comb"]

def ret_la2dat(var, Y,M,prtype1,prtype2, tag):
  sdir_root, sdir, varpath1, varpath2\
      =tag_func.ret_name_tagthpr_match_mon\
      (\
         var, prtype1, prtype2, tag, para_tag, Y, M\
      )
  a2var1 = fromfile(varpath1,  float32).reshape(ny,nx) 
  a2var2 = fromfile(varpath2,  float32).reshape(ny,nx) 
  return a2var1, a2var2

def ret_name_clim(prtype1,prtype2,tag,para_tag,iY,eY,season):
  lname = tag_func.ret_name_tagthpr_match_sumnum_clim\
          (\
           prtype1, prtype2, tag, para_tag, iY, eY, season\
          )

  sdir_root, sdir, amtname1, amtname2, ratname1, ratname2, sumname1, sumname2, numname\
        = lname
  return sdir_root, sdir, amtname1, amtname2, ratname1, ratname2, sumname1, sumname2, numname

def ret_la2clim(var,tag):
  sdir_root, sdir, varname1, varname2\
    = tag_func.ret_name_tagthpr_match_clim(var,prtype1, prtype2, tag, para_tag, iY, eY, season)
  a2var1  = fromfile(varname1, float32).reshape(ny,nx)
  a2var2  = fromfile(varname2, float32).reshape(ny,nx)
  return a2var1, a2var2

#***************************************************
if calcflag == True:
  for season in lseason:
    lM = ctrack_para.ret_lmon(season)
    #---- num plain ---------
    tag = "plain"
    a2num_pln   = zeros([ny,nx],float32)
    for Y,M in [[Y,M] for Y in lY for M in lM]:
      a2num_pln_tmp = ret_la2dat("num",Y,M,prtype1,prtype2,"plain")[0]
      a2num_pln     = a2num_pln + a2num_pln_tmp
    #---- calc precip (mean amount and rate) --
    for tag in ltag:
      a2sum1  = zeros([ny,nx],float32)
      a2sum2  = zeros([ny,nx],float32)
      a2num_tag   = zeros([ny,nx],float32)

      for Y,M in [[Y,M] for Y in lY for M in lM]:
  #      print Y,M,tag
        a2sum_tmp1, a2sum_tmp2 = ret_la2dat("sum",Y,M,prtype1,prtype2,tag)
        a2num_tmp = ret_la2dat("num",Y,M,prtype1,prtype2,tag)[0]

        a2sum1 = a2sum1 + a2sum_tmp1
        a2sum2 = a2sum2 + a2sum_tmp2
        a2num_tag  = a2num_tag  + a2num_tmp

        print Y,M,tag, a2sum_tmp1.mean(), a2sum_tmp2.mean()
      #-------------
      a2amt1 = (ma.masked_where(a2num_pln==0., a2sum1)/a2num_pln).filled(miss)
      a2amt2 = (ma.masked_where(a2num_pln==0., a2sum2)/a2num_pln).filled(miss)

      a2rat1 = (ma.masked_where(a2num_tag==0., a2sum1)/a2num_tag).filled(miss)
      a2rat2 = (ma.masked_where(a2num_tag==0., a2sum2)/a2num_tag).filled(miss)

      #-- write ----
      sdir_root, sdir, amtname1, amtname2, ratname1, ratname2, sumname1, sumname2, numname\
            = ret_name_clim\
              (\
               prtype1,prtype2,tag,para_tag,iY,eY,season\
              )
  
      ctrack_func.mk_dir(sdir)
      a2amt1.tofile(amtname1)
      a2amt2.tofile(amtname2)
      a2rat1.tofile(ratname1)
      a2rat2.tofile(ratname2)
      #a2sum1.tofile(sumname1)
      #a2sum2.tofile(sumname2)
      #a2num_tag.tofile(numname)
      print numname

 
#--------- Figures -----------
if figflag == True:
  for season in lseason:
    for tag in ["plain","ot","tc.comb","cf.comb"]:
      #for sumnum in ["amt","rat"]:
      for sumnum in ["pr","rat"]:
        #-- names --
        sdir_root, sdir, iname1, iname2\
          = tag_func.ret_name_tagthpr_match_clim(sumnum, prtype1, prtype2, tag, para_tag, iY, eY, season)

        figdir  = sdir + "/pict"
        ctrack_func.mk_dir(figdir)
        figname1= figdir + "/" +  iname1.split("/")[-1] + ".png"
        figname2= figdir + "/" +  iname2.split("/")[-1] + ".png"
        cbarname= figdir + "/" +  "cbar.%s.png"%(sumnum)
        #-- load ---
        a2dat1 = fromfile(iname1, float32).reshape(ny,nx)
        a2dat2 = fromfile(iname2, float32).reshape(ny,nx)
    
        a2dat1 = a2dat1[200:1000]  # 60S-60N -> 40S-40N
        a2dat2 = a2dat2[200:1000]
        #-- mask ---
        a2shade= a2dat1
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
          bnd  = False
        elif sumnum == "num":
          coef = 1.
          bnd  = [1,10,30,50,70,100,200,300] 
        #---- mycm  ----
        #mycm = "gist_heat_r" 
        mycm = "jet" 
        #---- title ----
        stitle1 = "%s %s %04d-%04d %s"%(prtype1, tag, iY, eY, season) 
        stitle2 = "%s %s %04d-%04d %s"%(prtype2, tag, iY, eY, season) 
        ctrack_fig.mk_map(a2dat1,bnd=bnd, mycm=mycm,soname=figname1, stitle=stitle1, cbarname=cbarname, miss=miss, lllat=lllat, lllon=lllon, urlat=urlat, urlon=urlon, a2shade=a2shade, coef=coef)
        ctrack_fig.mk_map(a2dat2,bnd=bnd, mycm=mycm,soname=figname2, stitle=stitle2, cbarname=cbarname, miss=miss, lllat=lllat, lllon=lllon, urlat=urlat, urlon=urlon, a2shade=a2shade, coef=coef)
        print figname2
#  

#---- test: check data: clim ----------
season = 1
a2tot1 = zeros([ny,nx],float32)
a2tot2 = zeros([ny,nx],float32)

for tag in ["tc.comb","cf.comb","ot","plain"]:
#for tag in ["plain","tc","c","fbc","te","tf","ef","tef","ot"]:
#for tag in ltag:
  a2dat1, a2dat2 = ret_la2clim("pr",tag)
  

  a2dat1  = ma.masked_equal(a2dat1, miss).filled(0.0)
  a2dat2  = ma.masked_equal(a2dat2, miss).filled(0.0)

  print tag,a2dat1.sum(), a2dat2.sum()
  if tag !="plain":
    a2tot1 = a2tot1 + a2dat1
    a2tot2 = a2tot2 + a2dat2
  else:
    a2plain1 = a2dat1
    a2plain2 = a2dat2

a2d1  = a2tot1 - a2plain1
a2d2  = a2tot2 - a2plain2

print "-- a2d.mean, a2d.min, a2d,max --"
print a2d1.min(), a2d1.mean(), a2d1.max()
print a2d2.min(), a2d2.mean(), a2d2.max()



##---- test: check data: mon -----------
#Y = 2002
##lM   = [6,7,8]
#lM   = [6]
#for M in lM:
#  tot1 = 0.
#  tot2 = 0.
#  totnum  = 0.
#  a2tot1 = zeros([ny,nx],float32)
#  a2tot2 = zeros([ny,nx],float32)
#  a2totnum = zeros([ny,nx],float32)
#  
#  for tag in ltag:
#    a2sum1, a2sum2, a2num = ret_la2dat(Y,M,prtype1,prtype2, tag)
#    if tag != "plain":
#      tot1 = tot1 + a2sum1.sum()
#      tot2 = tot2 + a2sum2.sum()
#      totnum  = totnum + a2num.sum()
#  
#      a2tot1 = a2tot1 + a2sum1
#      a2tot2 = a2tot2 + a2sum2
#      a2totnum  = a2totnum  + a2num
# 
#    elif tag == "plain":
#      plain1 = a2sum1.sum()
#      plain2 = a2sum2.sum()
#  
#      a2plain1 = a2sum1
#      a2plain2 = a2sum2
#      a2plainnum = a2num
#  #print "tot, plain"
#  #print tot1, plain1
#  #print tot2, plain2
#  
#  a2d1  = a2tot1 - a2plain1
#  a2d2  = a2tot2 - a2plain2
#  a2dnum= a2totnum - a2plainnum
#  print "---------------"
#  print "Mon=",M
#  print "d.mean(), d.min(), dmax()"
#  print a2d1.mean(), a2d1.min(), a2d1.max()
#  print a2d2.mean(), a2d2.min(), a2d2.max()
#  print a2dnum.mean(), a2dnum.min(), a2dnum.max() 
# 
