from numpy import *
from pyhdf import SD
import h5py
import os, calendar, sys
from datetime import *
from collections import deque
from glob import glob
from pmm_fsub import *
import pmm_para
import pmm_func
from cf2.io.GPM import GPM
from cf2.io.TRMM import TRMM


#THIS CODE IS YET NOT READY TO RUN GSMaP

#------------------
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except:
    pass
#------------------
def ret_a2pr(prtype, year, mon, day, hour):
  # ret_a2backward_org returns saone like array (from south to north)
  if prtype in ["GSMaP.v5","GSMaP.gauge","RA"]:
    a2pr  = pmm_func.ret_a2backward_org(year,mon,day,hour,minute=0, prtype=prtype, maskmiss=False)
  #---------
  return a2pr 

#------------------
def global2prmap(a2global, maptype):
  if maptype in ["GSMaP.v5", "GSMaP.gauge"]:
    a2prmap  = a2global[300:1500] # 60S-60N (0-1800*0.1)
  elif maptype in ["RA"]:
    a2prmap  = a2global
  return a2prmap

#******************************
def ret_a2prmap_forward(maptype,Y,M,D,H):
  ttime   = datetime(Y,M,D,H) + timedelta(hours=1)
  tY,tM,tD,tH = ttime.year, ttime.month, ttime.day, ttime.hour
  a2prmap = pmm_func.ret_a2backward_org(tY,tM,tD,tH,minute=0,prtype=maptype,maskmiss=False)
  return a2prmap
#------------------------------
def ret_coord(maptype):
  if maptype == "RA":
    ny,nx     = 2800,3200
    dlat,dlon = 0.01, 0.01
    lat_first = 20.0  +0.005   # grid center
    lon_first = 118.0 +0.005  # grid center
  elif maptype in ["GSMaP.gauge", "GSMaP.v5"]:
    ny,nx     = 1200,3600
    dlat,dlon = 0.1, 0.1
    lat_first = -59.95   # grid center
    lon_first = 0.05  # grid center
  return lat_first, lon_first, dlat, dlon, ny, nx

#------------------------------
#  Main part
#------------------------------
iY,eY = 2014,2014
iM,eM = 4,4    # test
lY    = range(iY,eY+1)
lM    = range(iM,eM+1)
iD    = 1

#--- parameters ----
#obttype = "PR.2A25"
#obttype = "DPR.L2.DPR"
#obttype = "KaPR.L2"
#obttype = "KuPR.L2"
#obttype = "GMI.L2"

#maptype = "GSMaP.gauge"
maptype = "RA"

#obttype  = "GPM.KuPR"
#prdLv    = "L2"
#prdVer   = "02"
#varname  = "NS/SLV/precipRateESurface"

#*******************************************************
obttype_a  = "GPM.DPR"
prdLv_a    = "L2"
prdVer_a   = "03"
varname_a  = "NS/SLV/precipRateESurface"

#obttype_b  = "GPM.GMI"
#prdLv_b    = "L2"
#prdVer_b   = "03"  #IT WAS 02 BEFORE, I HAVE TO RE-RUN
#varname_b  = "S1/surfacePrecipitation"

obttype_c  = "TRMM.PR"
prdLv_c    = "L2A25"
prdVer_c   = "07"
varname_c  = "e_SurfRain"

#obttype_d  = "TRMM.TMI"
#prdLv_d    = "L2A12"
#prdVer_d   = "07"
#varname_d  = "surfacePrecipitation"
#********************************************************

#obt1 = 
#obt2

print obttype_a
print obttype_c

#ny,nx      = pmm_para.ret_nynx(maptype)#no maptype really, just same japan base area
#print "ny,nx",ny,nx

ny,nx       = 280, 320

#---- load orbit data ---------------
gpm_a = GPM(obttype_a, prdLv_a,prdVer_a)
#gpm_b = GPM(obttype_b, prdLv_b,prdVer_b)
#print gpm_a
#print gpm_b
gpm_c = TRMM(obttype_c, prdLv_c,prdVer_c)
#gpm_d = TRMM(obttype_d, prdLv_d,prdVer_d)

#print type(gpm_a)#:class 'cf2.io.GPM.gpm.GPM'

BBox= [[20.0, 118.0], [48.0, 150.0]]#Japan (type:list)

#----------------------------
a2one   = ones([ny,nx],float32)
a2zero  = zeros([ny,nx],float32)

#lat_first, lon_first, dlat, dlon, ny, nx = ret_coord(maptype)
gzflag  = False #type: boolean
miss    = -9999.
#hdftype, Grp, vname = ret_prjinfo(obttype)

print "ny,nx",ny,nx

for Y,M in [[Y,M] for Y in lY for M in lM]:
  eD = calendar.monthrange(Y,M)[1]
#  eD = 1 #test
  #-- Initialize ---------------
  a2sum_obt1  = zeros([ny,nx],float32)
  a2sum_obt2  = zeros([ny,nx],float32) #before map
  a2num       = zeros([ny,nx],float32)

  #**DICTIONARY**

  da2pr_tmp  = {}
  for dattype in ["obt1","obt2"]:
    da2pr_tmp[dattype]   = zeros([ny,nx],float32)

  a2num_tmp =  zeros([ny,nx],float32)
 
  print "ny,nx",[ny,nx]
  print "OBT1", shape(da2pr_tmp["obt1"])
  print "OBT2", shape(da2pr_tmp["obt2"])
 
  #-----------------------------
  for D in range(iD,eD+1):
#  for D in range(iD,iD+2):   # test
    #aDat,aLat,aLon,aDtime  = get_dat_days(Y,M,D,Y,M,D)

    idtime = datetime(Y,M,D,0)
    edtime = datetime(Y,M,D,0) + timedelta(days=1)
    ddtime = edtime - idtime
    dH     = int(ddtime.total_seconds() / 3600.)
    ldtime = [idtime + timedelta(hours=H) for H in range(dH)]
    for dtime in ldtime:
#    for dtime in ldtime[4:23+1]:   # test
      print "***********************"
      print dtime
      print "***********************"
      Y = dtime.year
      M = dtime.month
      D = dtime.day
      H = dtime.hour

 #----- cf module ----------
      sDTime = dtime
      eDTime = sDTime + timedelta(hours=1)

#******************************************************************************************************
#CHANGE a b c d
#******************************************************************************************************
      try:
        obt1   = gpm_c(varname_c, sDTime, eDTime, [[-89.99,-180.0],[89.99,179.99]])
        obt2   = gpm_a(varname_a, sDTime, eDTime, [[-89.99,-180.0],[89.99,179.99]])
#        obt    = gpm(varname, sDTime, eDTime, BBox)
      except ValueError:
        print ""
        print "SKIP!! ValueError", sDTime
        print ""
        continue
      except UnboundLocalError:
        print ""
        print "SKIP!! ValueError", sDTime
        print ""
        continue
      except IndexError :
        print ""
        print "SKIP!! IndexError", sDTime
        continue

#The following gets the data from the obt1 -> gpm_a -> GPM -> cf2.io.GPM
#It gets the global data of GPM and matches it to the area of Japan. At the same time the 2 retrieved areas (sum and num) are to make the average of precipitation from the high resolution of gpm data to a (0.1 resolution for radar amedas)
      #lout2_1   = pmm_fsub.obt2jp2800x3200(obt1.data.T, obt1.lon.T, obt1.lat.T)
      lout2_1   = pmm_fsub.obt2jp280x320(obt1.data.T, obt1.lon.T, obt1.lat.T)
      a2s2_1    = lout2_1[0].T
      a2n2_1    = lout2_1[1].T
      a2pr_1    = ma.masked_where(a2n2_1==0.0, a2s2_1)/a2n2_1
      a2pr_1    = a2pr_1 /60./60. # mm/h-->mm/s  # org
      a2pr_1    = a2pr_1.filled(miss)
 

      #lout2_2   = pmm_fsub.obt2jp2800x3200(obt2.data.T, obt2.lon.T, obt2.lat.T)
      lout2_2   = pmm_fsub.obt2jp280x320(obt2.data.T, obt2.lon.T, obt2.lat.T)
      a2s2_2    = lout2_2[0].T
      a2n2_2    = lout2_2[1].T
      a2pr_2    = ma.masked_where(a2n2_2==0.0, a2s2_2)/a2n2_2
      a2pr_2    = a2pr_2 /60./60. # mm/h-->mm/s  # org
      a2pr_2    = a2pr_2.filled(miss)

#  print "a2s2",a2s2.shape
#  print "a2n2",a2n2.shape
#  print "a2pr******",a2pr.shape
#  print "OBT1", shape(da2pr_tmp["obt1"])
#  print "OBT2", shape(da2pr_tmp["obt2"])
#  print "final"

      #------------------------------
      da2pr_tmp["obt1"] = a2pr_1 #loading the obt style pr
      da2pr_tmp["obt2"] = a2pr_2 #I ADDED ONE MORE
     
      #---------------------
      # load map-style precipitation NOT NECESSARY?
      #---------------------
      #da2pr_tmp["map"]  = ret_a2pr(maptype, Y, M, D, H)
      #da2pr_tmp["map"]  = ret_a2prmap_forward(maptype, Y, M, D, H)

      #--------------------
      # mask out obt1 where obt2==miss
      #--------------------
      da2pr_tmp["obt1"] = ma.masked_where(da2pr_tmp["obt2"]==miss, da2pr_tmp["obt1"]).filled(miss)

      #---------------------
      # mask out obt2  where obt1==miss
      #---------------------
      da2pr_tmp["obt2"] = ma.masked_where(da2pr_tmp["obt1"] ==miss, da2pr_tmp["obt2"]).filled(miss)

      #---------------------
      # mask out num where obt2==miss
      #---------------------
      #a2num_tmp        = ma.masked_where(da2pr_tmp["map"] ==miss, a2num_tmp).filled(0.0)
      a2num_tmp        = ma.masked_where(da2pr_tmp["obt2"] ==miss, a2one).filled(0.0)
      #---------------------
      # replace miss by zero
      #---------------------
      da2pr_tmp["obt1"] = ma.masked_equal(da2pr_tmp["obt1"], miss).filled(0.0)
      da2pr_tmp["obt2"] = ma.masked_equal(da2pr_tmp["obt2"], miss).filled(0.0)
      #--------------------
      a2sum_obt1  = a2sum_obt1 + da2pr_tmp["obt1"]
      a2sum_obt2  = a2sum_obt2 + da2pr_tmp["obt2"]
      a2num      = a2num     + a2num_tmp

  #******************
  # write data
  #******************
  odir_base = "/home/azariah/AH/out/COMP.OBT1.OBT2/%s.vs.%s"%(obttype_c, obttype_a)#in GMI DPR vv
  odir      = odir_base + "/%04d"%(Y)
  mk_dir(odir)
  print odir

  NY,NX     = shape(da2pr_tmp["obt2"])
  sumname_obt1= odir + "/sum.%s.%02d.%d.%d"%(obttype_c, M, NY, NX)
  sumname_obt2= odir + "/sum.%s.%02d.%d.%d"%(obttype_a, M, NY, NX)
  numname     = odir + "/num.%02d.%d.%d"%(M, NY, NX)

  a2sum_obt1  = array(a2sum_obt1, float32)
  a2sum_obt2  = array(a2sum_obt2, float32)
  a2num_obt2  = array(a2num, float32)#a2num_obt2 should have the same shape as obt1

  a2sum_obt1.tofile(sumname_obt1)
  a2sum_obt2.tofile(sumname_obt2)
  a2num.tofile(numname)
  print sumname_obt1
                          

