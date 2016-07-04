from numpy import *
from pyhdf import SD
import h5py
from datetime import datetime, timedelta
from pmm_fsub import *
from glob import glob
import pmm_func
import sys, os
import ctrack_func
#------------------------------
def ret_obtdir(prj,Y,M):
  if prj     == ["TRMM.PR","2A25"]:
    idir     = "/media/disk2/data/TRMM.PR/L2A25/%04d/%02d"%(Y,M)
    shead    = "T1PR"
    stail    = ".01"

  elif prj   == ["TRMM.TMI","2A12"]:
    idir     = "/media/disk2/data/TRMM.TMI/L2A12/%04d/%02d"%(Y,M)
    shead    = "T1TMI"
    stail    = ".01"

  elif prj   == ["GPM.DPR","L2.DPR"]:
    idir     = "/mnt/mizu.tank/hjkim/GPM/GPM.DPR/L2.DPR/%04d/%02d"%(Y,M)
    shead    = "GPMCOR_DPR"
    stail    = "B.h5"
  
  elif prj   == ["GPM.KaPR","L2"]:
    idir     = "/mnt/mizu.tank/hjkim/GPM/GPM.KaPR/L2/02/%04d/%02d"%(Y,M)
    shead    = "GPMCOR_KAR"
    stail    = "B.h5"

  elif prj   == ["GPM.KuPR","L2"]:
    idir     = "/mnt/mizu.tank/hjkim/GPM/GPM.KuPR/L2/%04d/%02d"%(Y,M)
    #obtpath   = "/mnt/mizu.tank/hjkim/GPM/GPM.KuPR/L2/2014/03/GPMCOR_KUR_1403121359_1531_000201_L2S_DU2_02B.h5"
    shead    = "GPMCOR_KUR"
    stail    = "B.h5"

  elif prj   == ["GPM.GMI","L2"]:
    idir     = "/mnt/mizu.tank/hjkim/GPM/GPM.GMI/L2/02/%04d/%02d"%(Y,M)
    shead    = "GPMCOR_GMI"
    stail    = "A.h5"


  return idir, shead, stail


#------------------------------
def ret_lobtpath(prj, Y,M, listflag=False):
  obtdir, shead, stail = ret_obtdir(prj, Y, M)
  if listflag==False:
    lobtpath = sorted(glob( os.path.join(obtdir, "%s*%s"%(shead,stail))))

  elif listflag == True:
    print "listflag==",listflag
    listdir_root = "/mnt/mizu.tank/utsumi/PMM/COMP.ORB.ORG/GranuleList"
    if prj   == ["GPM.DPR","L2.DPR"]:
      listdir   = listdir_root + "/GPM.DPR.L2"
    elif prj  == ["GPM.KaPR","L2"]:
      listdir   = listdir_root + "/GPM.KaPR.L2"
    elif prj  == ["GPM.KuPR","L2"]:
      listdir   = listdir_root + "/GPM.KuPR.L2"
    elif prj  == ["GPM.GMI","L2"]:
      listdir   = listdir_root + "/GPM.GMI.L2"

    sBD   = "LAT%05.2f-%05.2f.LON%05.2f-%05.2f"%(iBDlat, eBDlat, iBDlon, eBDlon)
    listname = listdir + "/Glist.%s.%04d.%02d.txt"%(sBD,Y,M)
    f = open(listname, "r"); lobtname = f.readlines(); f.close()   

    lobtpath = [os.path.join(obtdir, obtname.strip()) for obtname in lobtname]

  return lobtpath


#------------------------------
def get_dat(prj, spath):
  #------ HDF4 -----------------
  if prj in [["TRMM.PR","2A25"],["TRMM.TMI","2A12"]]:
    if prj == ["TRMM.PR","2A25"]:
      vname = "e_SurfRain"
    elif prj == ["TRMM.TMI","2A12"]:
      vname = "surfacePrecipitation"
    else:
      print "check prj",prj
      sys.exit()

    h4          = SD.SD(spath)
    a2dat       = h4.select(vname)[:]
    a2lat       = h4.select("Latitude")[:]
    a2lon       = h4.select("Longitude")[:]
    aY,aM,aD,aH,aMin = h4.select("Year")[:], h4.select("Month")[:], h4.select("DayOfMonth")[:], h4.select("Hour")[:], h4.select("Minute")[:]
    h4.end()

  #------ HDF5 -----------------
  elif prj in [ ["GPM.DPR","L2.DPR"]\
               ,["GPM.KaPR","L2"]   \
               ,["GPM.KuPR","L2"]   \
               ,["GPM.GMI","L2"]   \
              ]:

    if prj == ["GPM.DPR","L2.DPR"]:
      Grp = "NS"
      vname = "SLV/precipRateESurface"
    elif prj == ["GPM.KaPR","L2"]:
      Grp = "MS"
      vname = "SLV/precipRateESurface"
    elif prj == ["GPM.KuPR","L2"]:
      Grp = "NS"
      vname = "SLV/precipRateESurface"
    elif prj == ["GPM.GMI","L2"]:
      Grp = "S1"
      vname = "surfacePrecipitation"
    else:
      print "check prj",prj
      sys.exit()

    h5   = h5py.File(spath, "r")
    a2dat= h5["%s/%s"%(Grp,vname)][:]
    a2lat= h5["%s/Latitude"%(Grp)][:]
    a2lon= h5["%s/Longitude"%(Grp)][:]
    aY   = h5['%s/ScanTime/Year'%(Grp)][:].astype('int')
    aM   = h5['%s/ScanTime/Month'%(Grp)][:].astype('int')
    aD   = h5['%s/ScanTime/DayOfMonth'%(Grp)][:].astype('int')
    aH   = h5['%s/ScanTime/Hour'%(Grp)][:].astype('int')
    aMin = h5['%s/ScanTime/Minute'%(Grp)][:].astype('int')
    #aSec = h5['%s/ScanTime/Second'%(Grp)][:].astype('int')
    #aMicSec= h5['%s/ScanTime/MilliSecond'%(Grp)][:].astype('int')*1000
    h5.close()
  #-----
  adtime      = [datetime(Y,M,D,H,Min) for (Y,M,D,H,Min) in zip(aY,aM,aD,aH,aMin)]

  return a2dat, a2lat, a2lon, adtime
#------------------------------
def ret_a2prmap_forward(prmaptype,Y,M,D,H):
  ttime   = datetime(Y,M,D,H) + timedelta(hours=1)
  tY,tM,tD,tH = ttime.year, ttime.month, ttime.day, ttime.hour
  a2prmap = pmm_func.ret_a2backward_org(Y,M,D,H,minute=0,prtype=prmaptype,maskmiss=True).filled(miss)
  return a2prmap 
#------------------------------
def ret_coord(prmaptype):
  if prmaptype == "RA":
    ny,nx     = 2800,3200
    dlat,dlon = 0.01, 0.01
    lat_first = 20.0  +0.005   # grid center
    lon_first = 118.0 +0.005  # grid center
  elif prmaptype == "GSMaP.gauge":
    ny,nx     = 1200,3600
    dlat,dlon = 0.1, 0.1
    lat_first = -59.95   # grid center  
    lon_first = 0.05  # grid center

  return lat_first, lon_first, dlat, dlon, ny, nx

#------------------------------
def filter_domain(a1dat, a1lat, a1lon, prmaptype):
  lat_first, lon_first, dlat, dlon, ny, nx = ret_coord(prmaptype)
  lat_last = lat_first + dlat*(ny-1)
  lon_last = lon_first + dlon*(nx-1)

  #-- lon: negative to positive ----
  a1lon   = (ma.masked_greater_equal(a1lon, 0.0) + 360).data

  #-- filter outside target domain --
  #a1lat   = ma.masked_outside(a1lat, lat_first-0.5*dlat, lat_last+0.5*dlat) 
  #a1lon   = ma.masked_outside(a1lon, lon_first-0.5*dlon, lon_last+0.5*dlon) 
  a1lat   = ma.masked_outside(a1lat, iBDlat, eBDlat) 
  a1lon   = ma.masked_outside(a1lon, iBDlon, eBDlon) 


  a1lat   = ma.masked_where(a1lon.mask==True, a1lat)
  a1lon   = ma.masked_where(a1lat.mask==True, a1lon)

  a1dat   = ma.masked_where(a1lat.mask==True, a1dat)

  a1dat   = a1dat.compressed() 
  a1lat   = a1lat.compressed() 
  a1lon   = a1lon.compressed() 

 
  return a1dat, a1lat, a1lon
#------------------------------
def or_filter_miss(a1,a2,miss):
  a1 = ma.masked_equal(a1,miss)
  a1 = ma.masked_where(a2==miss, a1)
  a2 = ma.masked_where(a1.mask==True, a2)

  a1 = a1.compressed()
  a2 = a2.compressed()
  return a1, a2
#-------------------------------
def and_filter_miss(a1,a2,miss):
  am = ma.masked_where((a1==miss)&(a2==miss), zeros(a1.shape))
  a1 = ma.masked_where(am.mask==True, a1)
  a2 = ma.masked_where(am.mask==True, a2)

  a1 = a1.compressed()
  a2 = a2.compressed()
  return a1, a2
#------------------------------
def filter_negative(a1,a2,akey):
  a1 = ma.masked_less(a1, 0.0)
  a2 = ma.masked_where(a1.mask==True, a2)
  a1 = a1.compressed()
  a2 = a2.compressed()
  return a1, a2
#------------------------------
def obtXmap(prj,obtpath,prmaptype):
  a1ORB   = array([])
  a1MAP   = array([])
  #a1LAT   = array([])
  #a1LON   = array([])

  a2Dat, a2Lat, a2Lon, aDtime = get_dat(prj,obtpath)
  if len(a2Dat)==0.0:
    return a1ORB, a1MAP
    sys.exit()

  iY,iM,iD,iH = aDtime[0].year, aDtime[0].month, aDtime[0].day, aDtime[0].hour
  iDtime           = datetime(iY,iM,iD,iH)
  dMin = 60
  ii   = 0
  
  lat_first, lon_first, dlat, dlon, ny, nx = ret_coord(prmaptype)
  for dn in range(int(60*2/dMin)+2):
  #for dn in [0]:
    idtime = iDtime + timedelta(minutes=dMin*dn)
    edtime = iDtime + timedelta(minutes=dMin*(dn+1))
    ii = searchsorted(aDtime, idtime) 
    ei = searchsorted(aDtime, edtime) 
    if ii >= len(aDtime):
      break
  
    a2dat   = a2Dat [ii:ei]
    a2lat   = a2Lat [ii:ei]
    a2lon   = a2Lon [ii:ei]
    adtime  = aDtime[ii:ei]
   
    a1dat   = a2dat.flatten()
    a1lat   = a2lat.flatten()
    a1lon   = a2lon.flatten()
    a1dat, a1lat, a1lon = filter_domain(a1dat, a1lat, a1lon, prmaptype)
    if len(a1dat)==0:
      continue
    
    a2map = ret_a2prmap_forward(prmaptype, idtime.year, idtime.month, idtime.day, idtime.hour) 
    #-- For Validation -----
    #a2map = array([[lat]*3200 for lat in arange(20.0, 20.0+(2800)*0.01, 0.01)])
    #a2map = array([arange(118.0, 118.0+(3200)*0.01, 0.01)]*2800)
    #a1dat = a1lat
    #-----------------------
    a2map  = a2map.astype("float64")
    a2lat  = a2lat.astype("float64")
    a2lon  = a2lon.astype("float64")



    a1map = pmm_fsub.obt_match_map(a2map.T, lat_first, lon_first, dlat, dlon, a1lat, a1lon).astype("float32")
    #--- filter ---
    a1dat, a1map = filter_negative(a1dat, a1map, akey=a1dat) 
    a1dat, a1map = or_filter_miss(a1dat, a1map, miss)
    a1dat, a1map = and_filter_miss(a1dat, a1map, 0)
  
    #--- unit -----
    a1dat        = a1dat / (60.*60.)  # mm/hour --> mm/sec
    a1ORB = r_[a1ORB, a1dat]
    a1MAP = r_[a1MAP, a1map]
    #a1LAT = r_[a1LAT, a1lat]
    #a1LON = r_[a1LON, a1lon]
  return a1ORB, a1MAP 
#-------------------------------------

SetBD  = True
#listflag = False
listflag = True
#prj = ["TRMM.TMI","2A12"]
#prj = ["TRMM.PR","2A25"]
#prj = ["GPM.KaPR","L2"]
#prj = ["GPM.KuPR","L2"]
#prj = ["GPM.GMI","L2"]
#lprj = [ ["GPM.GMI","L2"], ["GPM.KaPR","L2"] ]
lprj = [ ["GPM.GMI","L2"] ]

#lprmaptype = ["RA", "GSMaP.gauge"]
#lprmaptype = ["GSMaP.gauge"]
lprmaptype = ["RA"]

iY,eY  = 2014,2014
lY     = range(iY, eY+1)
#iM,eM  = 3,6
iM,eM  = 5,6
lM     = range(iM, eM+1)
miss      = -9999.

for prj in lprj:
  for prmaptype in lprmaptype:
    lat_first, lon_first, dlat, dlon, ny, nx = ret_coord(prmaptype)
    if SetBD == True:
      iBDlat, eBDlat = 20.0, 48.0    # RA domain
      iBDlon, eBDlon = 118.0, 150.0  # RA domain
      #iBDlat, eBDlat = -59.95, 59.95    # RA domain
      #iBDlon, eBDlon = 0.05, 359.95  # RA domain
    else:
      iBDlat, eBDlat = lat_first, lat_first-0.5*dlat
      iBDlon, eBDlon = lon_first, lon_first-0.5*dlon
    
    for Y in lY:
      for M in lM:
        lobtpath = sorted(ret_lobtpath(prj, Y, M, listflag))
    
        a1Obt    = array([],float32)
        a1Map    = array([],float32)
        for obtpath in lobtpath:
          print obtpath
          a1obt, a1map = obtXmap(prj, obtpath, prmaptype)
          a1Obt = r_[a1Obt, a1obt]
          a1Map = r_[a1Map, a1map]
        
        a1Obt   = a1Obt.astype("float32")
        a1Map   = a1Map.astype("float32")
        a1Pair  = c_[a1Obt, a1Map]
        #---- write -----
        odir  = "/mnt/mizu.tank/utsumi/PMM/COMP.ORB.ORG/%04d/%02d"%(Y,M)
        oname = odir + "/%s.x.%s.bin"%(".".join(prj), prmaptype)
        
        ctrack_func.mk_dir(odir)
        a1Pair.tofile(oname)
        print a1Pair
        print oname
        
