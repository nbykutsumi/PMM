from ftplib import FTP
from datetime import datetime
import myfunc.util as util
import calendar
import os, sys
import socket
from datetime import datetime, timedelta
import numpy as np

hostname  = "arthurhou.pps.eosdis.nasa.gov"
irootDir = "/gpmallversions"
#irootDir = "/sm/730/gpmdata"

myhost = socket.gethostname()
if  myhost =="mizu":
    orootDir   = "/home/utsumi/mnt/wellshare/data/GPM"
elif myhost =="shui":
    orootDir   = "/work/hk02/PMM/NASA"
    #orootDir   = "/tank/utsumi/data/PMM/NASA"
elif myhost =="well":
    #orootDir   = "/media/disk2/share/data/GPM"
    #orootDir   = "/media/disk2/share/data/PMM/NASA"
    orootDir   = "/home/utsumi/mnt/lab_work/hk02/PMM/NASA"

#GPM/TRMM.TMI/L2A12/07/2014/


iDTime  = datetime(2018,1,1)
eDTime  = datetime(2018,12,31)
dDTime  = timedelta(days=1)
lDTime  = util.ret_lDTime(iDTime, eDTime, dDTime)

#spec      = ["TRMM","TMI","2A-CLIM","gprof","V05","A"]
#spec      = ["GPM","GMI","1C","1C","V05","A"]
#spec      = ["GPM","GMI","2A-CLIM","gprof","V05","A"] # input=ECMWF
#spec      = ["GPM","GMI","2A","gprof","V05","A"] # input=GANAL
#spec      = ["GPM","Ku","2A","radar","V06","A"]
#spec      = ["GPM","Ka","2A","radar","V06","A"]
#spec      = ["GPM","DPR","2A","radar","V06","A"]
#spec      = ["GPM","DPRGMI","2B","radar","V06","A"]

#gpr_gmi        = ["GPM","GMI","2A-CLIM","gprof","V05"] # input=GANAL
#gpr_amsr2      = ["GCOMW1","AMSR2","2A-CLIM","gprof","V05"]
#gpr_ssmis_f16  = ["F16","SSMIS","2A-CLIM","gprof","V05"]
#gpr_ssmis_f17  = ["F17","SSMIS","2A-CLIM","gprof","V05"]
#gpr_ssmis_f18  = ["F18","SSMIS","2A-CLIM","gprof","V05"]
#gpr_atms_npp   = ["NPP","ATMS","2A-CLIM","gprof","V05"]
#gpr_atms_noaa20= ["NOAA20","ATMS","2A-CLIM","gprof","V05"]
#gpr_mhs_metopa = ["METOPA","MHS","2A-CLIM","gprof","V05"]
#gpr_mhs_metopb = ["METOPB","MHS","2A-CLIM","gprof","V05"]
#gpr_mhs_noaa18 = ["NOAA18","MHS","2A-CLIM","gprof","V05"]
#gpr_mhs_noaa19 = ["NOAA19","MHS","2A-CLIM","gprof","V05"]


gmi       = ["GPM","GMI","1C","1C","V05"]
amsr2     = ["GCOMW1","AMSR2","1C","1C","V05"]
ssmis_f16 = ["F16","SSMIS","1C","1C","V05"]
ssmis_f17 = ["F17","SSMIS","1C","1C","V05"]
ssmis_f18 = ["F18","SSMIS","1C","1C","V05"]
atms_npp  = ["NPP","ATMS","1C","1C","V05"]
atms_noaa20= ["NOAA20","ATMS","1C","1C","V05"]

mhs_metopa= ["METOPA","MHS","1C","1C","V05"]
mhs_metopb= ["METOPB","MHS","1C","1C","V05"]
mhs_noaa18= ["NOAA18","MHS","1C","1C","V05"]
mhs_noaa19= ["NOAA19","MHS","1C","1C","V05"]

#lspec = [gmi, amsr2, ssmis_f16,ssmis_f17,ssmis_f18,atms_npp, atms_noaa20, mhs_metopa, mhs_metopb, mhs_noaa18, mhs_noaa19]
lspec = [gmi, amsr2, ssmis_f16,ssmis_f17,ssmis_f18,atms_npp, atms_noaa20, mhs_metopa, mhs_metopb, mhs_noaa18, mhs_noaa19]

#lspec = [gpr_gmi, gpr_amsr2, gpr_ssmis_f16, gpr_ssmis_f17, gpr_ssmis_f18, gpr_atms_npp, gpr_atms_noaa20, gpr_mhs_metopa, gpr_mhs_metopb, gpr_mhs_noaa18, gpr_mhs_noaa19]


myid      = "nbyk.utsumi@gmail.com"
mypass    = "nbyk.utsumi@gmail.com"

#----------------------------------
def mk_dir(sodir):
  try:
    os.makedirs(sodir)
  except OSError:
    pass


ftp = FTP(hostname)
ftp.login(myid, mypass)

#----------------------------------
for spec in lspec:
    sate      = spec[0]
    sensor    = spec[1]
    prdName   = spec[2]
    prj       = spec[3]
    ver       = spec[4]

    for DTime in lDTime:
        Year,Mon,Day = DTime.timetuple()[:3]


        #if (sate=='GPM')&(): continue # test
        #if (sate=='GCOMW1')&(DTime< datetime(2018,9,18)): continue # test
        ##--- Skip missing or problematic data -----
        #ldataoutage=[
        # [2014,5,20]  # 1274-1277
        #,[2014,10,22] # 3682-
        #,[2014,10,23] # all
        #,[2014,10,24] # -3718
        #,[2014,12,9]  # 004425, 004426(missing some lat&lon)
        #,[2014,12,10]  # 004447-004452
        #,[2017,9,26]
        #,[2017,9,27]
        #,[2017,9,28]
        #,[2017,9,29]
        #]


        #if [Year,Mon,Day] in ldataoutage:
        #    print 'Skip data outage day'
        #    continue

        #if Mon==7: continue  # test
        #--- path and directory: Remote -----------
        iDir = irootDir + "/%s/%04d/%02d/%02d/%s"%(ver,Year,Mon,Day,prj)
        #oDir = orootDir + "/%s.%s/%s/%s/%04d/%02d"%(sate,sensor,prdName,ver,Year,Mon) 
        oDir = orootDir + "/%s.%s/%s/%s/%04d/%02d/%02d"%(sate,sensor,prdName,ver,Year,Mon,Day) 
        #GPM/TRMM.TMI/L2A12/07/2014/
        mk_dir(oDir)
        #--- list --------------
        lPath = np.sort(ftp.nlst(iDir))
        print iDir
        print lPath
        for sPath in lPath:
            oid = int(sPath.split('.')[-3])


            fName = os.path.basename(sPath)
            prdNameTmp, sateTmp, sensorTmp, algFullTmp, dtime, gNum, verFullTmp, sfx = fName.split('.')
            ymd = dtime.split('-')[0]
            if fName.split('.')[3]=='GPM-SLH': continue
            if ymd        !='%04d%02d%02d'%(Year,Mon,Day): continue
            if prdNameTmp !=prdName: continue
            if sateTmp    !=sate   : continue
            if sensorTmp  !=sensor  : continue
            if verFullTmp[:-1] !=ver: continue
            oPath = oDir + "/" + sPath.split("/")[-1]
            ftp.retrbinary("RETR %s"%(sPath), open(oPath,"wb").write)
            print oPath

ftp.close()

