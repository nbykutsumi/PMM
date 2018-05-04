import os, sys
import glob
import calendar
from datetime import datetime, timedelta
import myfunc.util as util
from numpy import *

def ret_lYM(iYM, eYM):
  """
  iYM = [iYear, iMon], eYM = [eYear, eMon]
  """
  iYear, iMon = iYM
  eYear, eMon = eYM
  lYM = []
  for Year in range(iYear, eYear+1):
    if iYear == eYear:
      lMon = range(iMon,eMon+1)
    elif Year == iYear:
      lMon = range(iMon,12+1)
    elif Year == eYear:
      lMon = range(1,eMon+1)
    else:
      lMon = range(1,12+1)

    for Mon in lMon:
      lYM.append([Year,Mon])
  return lYM


class GPMGV(object):
    def __init__(self):
        self.rootDir = "/work/a01/utsumi/data/GPMGV"

    def load_sitelist(self):
        self.listDir = self.rootDir + "/sitelist"
        self.listPath= self.listDir + "/sitelist.csv"
        f=open(self.listPath, "r"); lines = f.readlines(); f.close()

        diYM = {}
        deYM = {}
        dlYM = {}
        dnwNames = {}
        dgNames  = {}
        dnwCode = {}
        dlatlon = {}
        DNoYM   = {}
        Lkey    = []
        Lregion = [] 

        for line in lines[1:]:
            Line   = line.strip().split(",")
            Region = line[0]
            NwName = line[1]
            nwCode = line[2]
            gName  = line[3]
            lat    = float(line[4])
            lon    = float(line[5])
            iYear  = int(line[6])
            iMon   = int(line[7])
            eYear  = int(line[8])
            eMon   = int(line[9])
            nNoYM   = int(line[10])
            if nNoYM !=0:
                lNoYM = [[int(YM[:4]), int(YM[4:])] for YM in line[11:]]
            else:
                lNoYM = []
            #-----
            lkey.append((region,nwName,gName))
            key = (region,nwName,gName)

            if region not in dnwNames.keys():
                dnwNames[region] = [nwName]
                lregion.append(region)
            else:
                if nwName not in dnwNames[region]:
                    dnwNames[region].append(nwName)

            if (region,nwName) not in dgNames.keys():
                dgNames[region, nwName] = [gName]
            else:
                if gName not in dgNames[region,nwName]:   
                     dgNames[region, nwName].append(gName)

            dnwCode[key] = nwCode
            dlatlon[key] = [lat,lon]
            diYM[key]    = [iYear,iMon]
            deYM[key]    = [eYear,eMon]

            # lYM             
            lYMtmp = ret_lYM([iYear,iMon],[eYear,eMon])
            lYM    = []
            for YM in lYMtmp:
                if YM not in lNoYM: 
                    lYM.append(YM)

            dlYM[key] = lYM
 
        #---------
        self.diYM = diYM
        self.diYM = deYM
        self.dlYM = dlYM
        self.dnwNames = dnwNames
        self.dgNames  = dgNames
        self.dnwCode  = dnwCode
        self.dlatlon  = dlatlon
        self.keys     = lkey
        self.regions  = lregion


    def load_sitelist_reclassified(self):
        self.listDir = self.rootDir + "/sitelist"
        self.listPath= self.listDir + "/sitelist_reclassified.csv"
        f=open(self.listPath, "r"); lines = f.readlines(); f.close()

        diYM = {}
        deYM = {}
        dlYM = {}
        dgNames  = {}
        dnwCode = {}
        dlatlon = {}
        dNoYM   = {}
        lkey    = []
        ldomain = [] 

        for line in lines[1:]:
            line   = line.strip().split(",")
            region = line[0]
            nwName = line[1]
            nwCode = line[2]
            domain = line[3]
            gName  = line[4]
            lat    = float(line[5])
            lon    = float(line[6])
            iYear  = int(line[7])
            iMon   = int(line[8])
            eYear  = int(line[9])
            eMon   = int(line[10])
            NoYM   = int(line[11])
            if NoYM !=0:
                lNoYM = [[int(YM[:4]), int(YM[4:])] for YM in line[12:]]
            else:
                lNoYM = []
            #-----
            lkey.append((domain,gName))
            key = (domain,gName)

            if domain not in ldomain:
                ldomain.append(domain)

            if domain not in dgNames.keys():
                dgNames[domain] = [gName]
            else:
                if gName not in dgNames[domain]:   
                     dgNames[domain].append(gName)

            dnwCode[key] = nwCode
            dlatlon[key] = [lat,lon]
            diYM[key]    = [iYear,iMon]
            deYM[key]    = [eYear,eMon]

            # lYM             
            lYMtmp = ret_lYM([iYear,iMon],[eYear,eMon])
            lYM    = []
            for YM in lYMtmp:
                if YM not in lNoYM: 
                    lYM.append(YM)

            dlYM[key] = lYM
 
        #---------
        self.diYM = diYM
        self.deYM = deYM
        self.dlYM = dlYM
        self.dgNames  = dgNames
        self.dnwCode  = dnwCode
        self.dlatlon  = dlatlon
        self.keys     = lkey
        self.domains  = ldomain


    def ret_ddomYM2gName(self):
        # load list
        self.load_sitelist_reclassified() 

        ldomain = self.domains
        dgNames = self.dgNames
        dlYM = self.dlYM
        
        dgOut = {}
        lYM   = ret_lYM([1995,1],[2018,12])
        for domain in ldomain:
            lgName = dgNames[domain]
            for YM in lYM:
                Year, Mon = YM
                lout = []
                for gName in lgName:
                    if YM in dlYM[domain, gName]:
                        lout.append(gName)
                if len(lout) !=0:
                    dgOut[domain, Year, Mon] = lout
        return dgOut


    def load_obsfile_asc(self, srcPath):
        f=open(srcPath,"r"); lines=f.readlines(); f.close()
        head   = lines[0].split()
        prdID  = head[0] 
        nwCode = head[1] + "_" + head[2]
        gID    = head[3]
        locName= head[4]
        gType  = head[5]
        timeres= head[6]
        lat    = head[7]
        lon    = head[8]
        radar  = head[9]
        radar_distance = head[10]
        azimuth= head[11]
        gX     = head[12]
        gY     = head[13] 
        radarElev=head[14]
        sDate  = head[15]
        sTime  = head[16]
        eDate  = head[17]
        eTime  = head[18]

        MonFile = int(sDate.split("/")[0])    
        YearFile= int(sDate.split("/")[-1])
        eDay   = calendar.monthrange(YearFile,MonFile)[1]

        sDTime = datetime(YearFile,MonFile,1,0,0)
        eDTime = datetime(YearFile,MonFile,eDay,23,59)
        dDTime = timedelta(seconds=60)
        aDTime = array(util.ret_lDTime(sDTime,eDTime,dDTime))

        ndat   = len(aDTime)
        miss   = -9999.
        aPrcp  = zeros(ndat)

        for line in lines[1:]:
            line = line.split()
            Year,Mon,Day,doy,Hour,Mnt,Sec = map(int,line[:7])
            prcp = float(line[7])

            DTime= datetime(Year,Mon,Day,Hour,Mnt,Sec)
            idtime= int((DTime-sDTime).total_seconds()/60.)

            aPrcp[idtime] = prcp

        return head, aDTime, aPrcp

    def load_obsfile_2a56(self, srcPath):
        f=open(srcPath,"r"); lines=f.readlines(); f.close()
        head   = lines[0].split()
        Year   = head[0]
        nwCode = head[1]
        gID    = head[2]
        gType  = head[3]
        gsize  = head[4]
        sDate  = head[6]
        sTime  = head[7]
        eDate  = head[8]
        eTime  = head[9]

        MonFile = int(sDate.split("/")[0])    
        YearFile= int(sDate.split("/")[-1])
        eDay   = calendar.monthrange(YearFile,MonFile)[1]

        sDTime = datetime(YearFile,MonFile,1,0,0)
        eDTime = datetime(YearFile,MonFile,eDay,23,59)
        dDTime = timedelta(seconds=60)
        aDTime = array(util.ret_lDTime(sDTime,eDTime,dDTime))

        ndat   = len(aDTime)
        miss   = -9999.
        aPrcp  = zeros(ndat)

        for line in lines[2:]:
            line = line.split()
            Year,Mon,Day,doy,Hour,Mnt = map(int,line[:6])
            prcp = float(line[6])

            DTime= datetime(Year,Mon,Day,Hour,Mnt)
            idtime= int((DTime-sDTime).total_seconds()/60.)

            aPrcp[idtime] = prcp

        return head, aDTime, aPrcp





if __name__ == "__main__":
    gv = GPMGV()

    #srcPath = "/work/a01/utsumi/data/GPMGV/2A56/FLORIDA/SFL/2014/2A56_MELB_SFL_0310_201412_4.asc"
    #head, aDTime, aPrcp = gv.load_obsfile_asc(srcPath)

    srcPath = "/work/a01/utsumi/data/GPMGV/2A56/MARYLAND/GSFC/2014/GSFC-Bldg33_Z-201412.2a56"
    head, aDTime, aPrcp = gv.load_obsfile_2a56(srcPath)

    sDTime = aDTime[0]
    DTime  = datetime(2014,12,29,9,0)
    i      = int((DTime - sDTime).total_seconds()/60)
    print aDTime[i],aPrcp[i]

