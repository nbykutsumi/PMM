import os, sys
import glob

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
        dNoYM   = {}
        lkey    = []
        lregion = [] 

        for line in lines[1:]:
            line   = line.strip().split(",")
            region = line[0]
            nwName = line[1]
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
            import globnNoYM   = int(line[11])
            if nNoYM !=0:
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
        self.diYM = deYM
        self.dlYM = dlYM
        self.dgNames  = dgNames
        self.dnwCode  = dnwCode
        self.dlatlon  = dlatlon
        self.keys     = lkey
        self.domains  = ldomain

    def load_obsfile(domain, nwCode, gCode, Year, Mon):
        region, nwName = domain.split("-")
        srcDir = self.rootDir + "/%s/%s/%04d"%(region,nwName,Year)
        ssearch= srcDir + "/2A56_%s_%s_%04d%02d_*.asc"%(
        srcPath= glob.glob(
 

if __name__ == "__main__":
    gv = GPMGV()
    gv.load_sitelist() 

    #print gv.diYM
    #print gv.deYM
    #print gv.dlYM
    #print gv.dnwNames 
    print gv.dgNames['N.Carolina', 'IPHEx_NASA'] 
