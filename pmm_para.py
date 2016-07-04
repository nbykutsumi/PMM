import sys

class Para_tag(object):
  dist_c  = 1000.      # km
  dist_tc = 1000.      # km
  dist_f  = 500.       # km
  thdura_tc  = 48      # hours
  thdura_c   = 48      # hours
  sresol     = "anl_p"
  bstflag_tc = "bst"
  bstflag_f  = ""

  def __init__(self):
    print "Para_tag"



def ret_nynx(prtype):
  if prtype in ["GSMaP.v5","GSMaP.gauge"]:
    ny,nx = 1200,3600
  elif prtype in ["RA"]:
    ny,nx = 2800,3200
  else:
    print "check prtype!!", prtype
    print "by pmm_para.ret_nynx"
    sys.exit()
  return ny,nx


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

  elif region == "EAS":
    lllon  = 115
    urlon  = 150
    lllat  = 30
    urlat  = 40
    iyfig  = (lllat + 60)*10  # 15N-40N
    eyfig  = (urlat + 60)*10
    ixfig  = (lllon)*10
    exfig  = (urlon)*10

  elif region == "SAS.SEA":
    lllon  = 65
    urlon  = 110
    lllat  = 0.0
    urlat  = 25
    iyfig  = (lllat + 60)*10  # 15N-40N
    eyfig  = (urlat + 60)*10
    ixfig  = (lllon)*10
    exfig  = (urlon)*10


  else:
    print "check region!"
    print "by pmm_para"
    sys.exit()
  return lllon, urlon, lllat, urlat, iyfig, eyfig, ixfig, exfig  
