from numpy import *
import datetime
import sys, os
#import tag_func
#import gsmap_func
#from gsmap_fsub import *
#***********************************
def projection_a2tag(a2tag_one,prtype="GSMaP.v5"):
  if prtype in ["GSMaP.v5","GSMaP.gauge"]:
    a2tag_one = gsmap_func.global2gsmap_one(a2tag_one)
    a2tag     = gsmap_fsub.saone_gsmap2dec_gsmap(a2tag_one.T).T

  return a2tag

#***************************************
def global2ra_one(a2glob):
  lllon  = 118
  urlon  = 150
  lllat  = 20
  urlat  = 48
  iyfig  = (lllat + 90) 
  eyfig  = (urlat + 90)
  ixfig  = (lllon)
  exfig  = (urlon)

  a2ra  = a2glob[iyfig:eyfig,ixfig:exfig]
  return a2ra
#***************************************
def ret_tagtime(year,mon,day,hour):
  now = datetime.datetime(year,mon,day,hour)
  if hour in [0,6,12,18]:
    tagtime = now
  elif hour in [1,7,13,19]:
    tagtime = now + datetime.timedelta(hours=-1)
  elif hour in [2,8,14,20]:
    tagtime = now + datetime.timedelta(hours=-2)
  elif hour in [3,9,15,21]:
    tagtime = now + datetime.timedelta(hours=-3)
  elif hour in [4,10,16,22]:
    tagtime = now + datetime.timedelta(hours=+2)
  elif hour in [5,11,17,23]:
    tagtime = now + datetime.timedelta(hours=+1)
  return tagtime.year, tagtime.month, tagtime.day, tagtime.hour

#***************************************
def ret_da2tag(year, mon, day, hour, para):
  #-----------------------------
  # parameters
  #-----------------------------
  dist_tc    = para.dist_tc
  dist_c     = para.dist_c
  dist_f     = para.dist_f
  thdura_tc  = para.thdura_tc
  thdura_c   = para.thdura_c
  bstflag_tc = para.bstflag_tc
  bstflag_f  = para.bstflag_f
  prtype     = para.prtype
  sresol     = para.sresol
  wnflag     = para.wnflag
  #-----------------------------
  if wnflag == False:
    ltag      = ["tc","c","fbc","te","tf","ef","tef","ot"]
    ltag_elem = ["tc","c","fbc","te","tf","ef","tef"]  
  else:
    print "check wnflag!!  wnflag=",wnflag
    print "by pmm_func"
    sys.exit()
  #-----------------------------
  # load corresponding tags
  #-----------------------------
  year_target, mon_target, day_target, hour_target\
     = ret_tagtime(year,mon,day,hour)
  #------------------
  #-- tag name ---
  tagname  = tag_func.ret_name_tagraw(sresol, thdura_tc, thdura_c, dist_tc, dist_c, dist_f, bstflag_tc, bstflag_f, year_target, mon_target, day_target, hour_target)

  if not os.access(tagname[0], os.F_OK):
    print "AAAA"
    print "nofile", tagname

  #-- load -------
  la2tagtmp    =   tag_func.ret_a2tagmask(sresol, thdura_tc, thdura_c, dist_tc, dist_c, dist_f, bstflag_tc, bstflag_f, year_target, mon_target, day_target, hour_target)


  a2tagtres_tc  =la2tagtmp[0]
  a2tagtres_c   =la2tagtmp[1]
  a2tagtres_fbc =la2tagtmp[2]
  a2tagtres_nbc =la2tagtmp[3]

  a2zero_tres = zeros(shape(la2tagtmp[0]))
  da2tag_tres      = {}
  if wnflag == False:
    da2tag_tres["tc" ] = ma.masked_where( reduce(logical_or, (a2tagtres_fbc ==1.0, a2tagtres_c   ==1.0)), a2tagtres_tc ).filled(0.0)

    da2tag_tres["c"  ] = ma.masked_where( reduce(logical_or, (a2tagtres_tc  ==1.0, a2tagtres_fbc ==1.0)), a2tagtres_c  ).filled(0.0)

    da2tag_tres["fbc"] = ma.masked_where( reduce(logical_or, (a2tagtres_tc  ==1.0, a2tagtres_c   ==1.0)), a2tagtres_fbc).filled(0.0)

    da2tag_tres["te" ] = ma.masked_where( (a2tagtres_tc ==1.0)&(a2tagtres_c ==1.0)&(a2tagtres_fbc ==0.0), a2zero_tres).filled(1.0)

    da2tag_tres["tf" ] = ma.masked_where( (a2tagtres_tc ==1.0)&(a2tagtres_c ==0.0)&(a2tagtres_fbc ==1.0), a2zero_tres).filled(1.0)

    da2tag_tres["ef" ] = ma.masked_where( (a2tagtres_tc ==0.0)&(a2tagtres_c ==1.0)&(a2tagtres_fbc ==1.0), a2zero_tres).filled(1.0)

    da2tag_tres["tef"] = ma.masked_where( (a2tagtres_tc ==1.0)&(a2tagtres_c ==1.0)&(a2tagtres_fbc ==1.0), a2zero_tres).filled(1.0)

  if wnflag == True:
    da2tag_tres["wn.tc" ] = ma.masked_where( reduce(logical_or, (a2tagtres_c  ==1.0, a2tagtres_fbc ==1.0, a2tagtres_nbc ==1.0)), a2tagtres_tc ).filled(0.0)

    da2tag_tres["wn.c"  ] = ma.masked_where( reduce(logical_or, (a2tagtres_tc ==1.0, a2tagtres_fbc ==1.0, a2tagtres_nbc ==1.0)), a2tagtres_c  ).filled(0.0)

    da2tag_tres["wn.fbc"] = ma.masked_where( reduce(logical_or, (a2tagtres_tc ==1.0, a2tagtres_c   ==1.0, a2tagtres_nbc ==1.0)), a2tagtres_fbc).filled(0.0)

    da2tag_tres["wn.nbc"] = ma.masked_where( reduce(logical_or, (a2tagtres_tc ==1.0, a2tagtres_c   ==1.0, a2tagtres_fbc ==1.0)), a2tagtres_nbc).filled(0.0)

    da2tag_tres["wn.te" ] = ma.masked_where( (a2tagtres_tc ==1.0)&(a2tagtres_c  ==1.0)&( a2tagtres_fbc ==0.0)&( a2tagtres_nbc ==0.0), a2zero_tres).filled(1.0)

    da2tag_tres["wn.tf" ] = ma.masked_where( (a2tagtres_tc ==1.0)&(a2tagtres_fbc==1.0)&( a2tagtres_c   ==0.0)&( a2tagtres_nbc ==0.0), a2zero_tres).filled(1.0)

    da2tag_tres["wn.tq" ] = ma.masked_where( (a2tagtres_tc ==1.0)&(a2tagtres_nbc==1.0)&( a2tagtres_c   ==0.0)&( a2tagtres_fbc ==0.0), a2zero_tres).filled(1.0)

    da2tag_tres["wn.ef" ] = ma.masked_where( (a2tagtres_c  ==1.0)&(a2tagtres_fbc==1.0)&( a2tagtres_tc  ==0.0)&( a2tagtres_nbc ==0.0), a2zero_tres).filled(1.0)

    da2tag_tres["wn.eq" ] = ma.masked_where( (a2tagtres_c  ==1.0)&(a2tagtres_nbc==1.0)&( a2tagtres_tc  ==0.0)&( a2tagtres_fbc ==0.0), a2zero_tres).filled(1.0)

    da2tag_tres["wn.fq" ] = ma.masked_where( (a2tagtres_fbc==1.0)&(a2tagtres_nbc==1.0)&( a2tagtres_tc  ==0.0)&( a2tagtres_c   ==0.0), a2zero_tres).filled(1.0)

    da2tag_tres["wn.tef"] = ma.masked_where( (a2tagtres_tc ==1.0)&(a2tagtres_c ==1.0)&(a2tagtres_fbc ==1.0)&(a2tagtres_nbc ==0.0), a2zero_tres).filled(1.0)

    da2tag_tres["wn.teq"] = ma.masked_where( (a2tagtres_tc ==1.0)&(a2tagtres_c ==1.0)&(a2tagtres_fbc ==0.0)&(a2tagtres_nbc ==1.0), a2zero_tres).filled(1.0)

    da2tag_tres["wn.tfq"] = ma.masked_where( (a2tagtres_tc ==1.0)&(a2tagtres_c ==0.0)&(a2tagtres_fbc ==1.0)&(a2tagtres_nbc ==1.0), a2zero_tres).filled(1.0)

    da2tag_tres["wn.efq"] = ma.masked_where( (a2tagtres_tc ==0.0)&(a2tagtres_c ==1.0)&(a2tagtres_fbc ==1.0)&(a2tagtres_nbc ==1.0), a2zero_tres).filled(1.0)

    da2tag_tres["wn.tefq"]= ma.masked_where( (a2tagtres_tc ==1.0)&(a2tagtres_c ==1.0)&(a2tagtres_fbc ==1.0)&(a2tagtres_nbc ==1.0), a2zero_tres).filled(1.0)

  #------- ot -------------
  a2tag_non_ot    = a2zero_tres
  for tag in ltag_elem:
    a2tag_non_ot  = a2tag_non_ot + da2tag_tres[tag]

  da2tag_tres["ot"]= ma.masked_where(a2tag_non_ot==0.0, a2zero_tres).filled(1.0)

  #------- projection -----
  da2tag = {}
  for tag in ltag:
    da2tag[tag] = projection_a2tag(da2tag_tres[tag], prtype)

  #-----------------------------------
  return da2tag

#***********************************
def trmm2global_one(a2org_one, miss):
  a2glob = ones([180,360], float32)*miss
  #a2glob[30:149+1,:] = a2org_one
  a2glob[50:130,:] = a2org_one
  return a2glob

#***********************************
def global2trmm_one(a2glob_one):
  a2org  = a2glob_one[50:130,:]
  return a2org


#------------------------------
def ret_a2backward_org(year,mon,day,hour,minute=0,prtype="RA",maskmiss=False):
  miss  = -9999.
  if prtype == "RA":
    #idir_root = "/mnt/iis.data1/hjkim/AMeDAS/ra_0.01"
    idir_root = "/data1/hjkim/AMeDAS/ra_0.01"
    idir      = idir_root + "/%04d%02d"%(year,mon)
    iname     = idir + "/RadarAmedas.%04d%02d%02d%02d%02d00.2800x3200"%(year,mon,day,hour,minute)
    a2pr      = ma.masked_equal(fromfile(iname, float32).reshape(2800,3200), -999.)
    a2pr      = flipud(a2pr) / (60.*60.)  # mm/hour --> mm/sec
  if prtype == "GSMaP.gauge":
    #idir_root = "/media/disk2/data/GSMaP/org.standard.gauge.v6/hourly"
    idir_root = "/data2/GSMaP/standard_gauge/v6/hourly"
    idir      = idir_root + "/%04d/%02d/%02d"%(year,mon,day)
    iname     = idir + "/gsmap_mvk.20140429.0500.v6.0000.0.dat"
    iname     = idir + "/gsmap_mvk.%04d%02d%02d.%02d00.v6.0000.0.dat"%(year,mon,day,hour)
    a2pr      = flipud(fromfile(iname, float32).reshape(1200,3600))
    a2pr      = ma.masked_less(a2pr, 0.0)/(60.*60.)  # mm/hour --> mm/sec

#  if prtype == "GSMaP.v5":
#    idir_root = "/media/disk2/data/GSMaP/org.v5/hourly"
#    idir      = idir_root + "/%04d/%02d/%02d"%(year,mon,day)
#    iname     = idir + "/gsmap_mvk.%04d%02d%02d.%02d00.v5.222.1.dat"%(year,mon,day,hour)
#    a2pr      = flipud(fromfile(iname, float32).reshape(1200,3600))
#    a2pr      = ma.masked_less(a2pr, 0.0)/(60.*60.)  # mm/hour --> mm/sec

  if maskmiss == False:
    a2pr = a2pr.filled(miss)
  return a2pr
#------------------------------

def timeave_backward_saone(year,mon,day,hour, hlen, prtype="PR",relaxflag=False):
  #***********************************
  # e.g. "backward" means:
  # Case for year,mon,day =2001,8,14, hour=2, hlen=2
  # average of 00UTC-02UTC
  # hour=0 ( mean precip of 00UTC-01UTC ) 
  # hour=1 ( mean precip of 01UTC-02UTC ) 
  #***********************************
  #lhlen = [1,2,3,6,12,24,72,168,672]
  #if not hlen in lhlen:
  #  print "'hlen' should be" ,lhlen
  #  print "by gpm_func.timeave_xxxx"
  #  sys.exit()

  def ret_iname(prtype,year,mon,day,hour):
    if prtype=="GSMaP":
      idir_root = "/media/disk2/data/GSMaP/sa.one.FWD/1hr/ptot"
      idir      = idir_root + "/%04d%02d"%(year,mon)
      iname     = idir + "/gsmap_mvk.1rh.%04d%02d%02d.%02d00.v5.222.1.sa.one"%(year_t,mon_t,day_t,hour_t)
    elif prtype=="PR": 
      idir_root = "/media/disk2/data/TRMM/PR.2A25.sa.one.FWD"
      idir      = idir_root + "/%04d/%02d"%(year,mon)
      iname     = idir + "/PR.%04d.%02d.%02d.%02d.sa.one"%(year,mon,day,hour)
    elif prtype=="TMI": 
      idir_root = "/media/disk2/data/TRMM/TMI.2A12.sa.one.FWD"
      idir      = idir_root + "/%04d/%02d"%(year,mon)
      iname     = idir + "/TMI.%04d.%02d.%02d.%02d.sa.one"%(year,mon,day,hour)

    return iname

  def ret_a2in(prtype,year,mon,day,hour):
    if prtype =="GSMaP.mtch.PR":
      iname   = ret_iname("GSMaP",year,mon,day,hour)
      mskname = ret_iname("PR",year,mon,day,hour)
      a2in    = fromfile(iname,   float32).reshape(120,360)
      a2mask  = fromfile(mskname, float32).reshape(80,360)
      a2mask  = trmm2global_one(a2mask, -9999.)[30:150]
      a2in    = ma.masked_where(a2mask==-9999., a2in).filled(-9999.)
    elif prtype =="TMI.mtch.PR":
      iname   = ret_iname("TMI",year,mon,day,hour)
      mskname = ret_iname("PR",year,mon,day,hour)
      a2in    = fromfile(iname,   float32).reshape(80,360)
      a2mask  = fromfile(mskname, float32).reshape(80,360)
      a2in    = ma.masked_where(a2mask==-9999., a2in).filled(-9999.)

    else:
      iname   = ret_iname(prtype,year,mon,day,hour)
      a2in    = fromfile(iname, float32).reshape(-1,360)

    #print iname
    return a2in

  #-------------
  def ret_FBinc(prtype):
    if prtype in ["PR","TMI","TMI.mtch.PR","GSMaP","GSMaP.mtch.PR"]:
      # input data= Forward average/accumulation
      FBinc = 1
    elif prtype in ["RA"]:
      # input data= Backward average/accumulation
      FBinc = 0
    return FBinc


  #-------------
  FBinc     = ret_FBinc(prtype)
  lh_inc    = range(hlen)
  now       = datetime.datetime(year,mon,day,hour)

  if prtype in ["PR","TMI","TMI.mtch.PR"]:
    ny,nx     = 80,360
  if prtype in ["GSMaP","GSMaP.mtch.PR"]:
    ny,nx     = 120,360

  a2ave     = zeros([ny,nx],float32)
  #---------------
  #---------------
  if   relaxflag == False:
    for h_inc in lh_inc:
      #dhour   = datetime.timedelta(hours = -h_inc)
      #dhour   = datetime.timedelta(hours = -h_inc-1)  # input data is "forward" average
      dhour   = datetime.timedelta(hours = -h_inc -FBinc)  # forward-input: FBinc= 1, backward-input: FBinc=0
      target  = now + dhour
      year_t  = target.year
      mon_t   = target.month
      day_t   = target.day
      hour_t  = target.hour
      a2in    = ret_a2in(prtype,year_t,mon_t,day_t,hour_t)
      a2ave     = a2ave + ma.masked_equal(a2in, -9999.0)
  elif relaxflag == True:
    a2mask    = ones([ny,nx],float32)* (-9999.0)
    for h_inc in lh_inc:
      #dhour   = datetime.timedelta(hours = -h_inc)
      #dhour   = datetime.timedelta(hours = -h_inc -1)  # input data is "forward" average
      dhour   = datetime.timedelta(hours = -h_inc -FBinc)  # input data is "forward" average
      target  = now + dhour
      year_t  = target.year
      mon_t   = target.month
      day_t   = target.day
      hour_t  = target.hour
      a2in    = ret_a2in(prtype,year_t,mon_t,day_t,hour_t)
      a2ave   = a2ave + ma.masked_equal(a2in, -9999.0).filled(0.0)
      a2mask  = ma.masked_where(a2in != -9999.0, a2mask).filled(1.0)
  #---------------
  a2ave       = a2ave / hlen
  if relaxflag == False:
    a2ave     = a2ave.filled(-9999.0)
  elif relaxflag == True:
    a2ave     = ma.masked_where(a2mask==-9999.0, a2ave).filled(-9999.0)
  #---------------
  return a2ave

