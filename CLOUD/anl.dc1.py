from numpy import *
import myfunc.util as util

#clVer = "JMA1"
#clVer = "MyWNP1"
#clVer = "MyWNP3"
clVer = "MyWNP.M.3"

#iYM    = [2014,4]
#eYM    = [2015,6]
iYM    = [2014,5]
eYM    = [2015,6]

lYM    = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [11,12,1,2,3]]

#ldattype = ["RA","KuPR","GMI","GSMaP","GSMaP.IR","GSMaP.MW","IMERG","IMERG.IR","IMERG.MW"]
ldattype = ["IMERG.IR"]
#ldattype = ["GSMaP","GSMaP.IR","GSMaP.MW","IMERG","IMERG.IR","IMERG.MW"]

miss = -9999.
rootDir = "/home/utsumi/mnt/wellshare"
if clVer   == "JMA1":
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.JMA"
  ibaseDirCL = "/tank/utsumi/CLOUDTYPE/WNPAC"

elif clVer[:5] == "MyWNP":
  ver        = clVer[5:]
  ibaseDir   = rootDir + "/PMM/WNP.261x265/CL.My%s"%(ver)

icl  = 1

for dattype in ldattype:
    apr   = array([], float32)
    aku   = array([], float32)
    aLocalTime = array([], float32)

    for [Year,Mon] in lYM:
        print Year,Mon
        baseDir = ibaseDir
        sDir    = baseDir + "/VsKuPR.LocTime/%04d"%(Year)

        PrPath   = sDir + "/%s.%04d.%02d.bn"%(dattype,Year,Mon)
        KuPath   = sDir + "/KuPR.%04d.%02d.bn"%(Year,Mon)
        CLPath   = sDir + "/CL.%04d.%02d.bn"%(Year,Mon)
        LSPath   = sDir + "/LS.%04d.%02d.bn"%(Year,Mon)
        LonPath  = sDir + "/Lon.%04d.%02d.bn"%(Year,Mon)
        LatPath  = sDir + "/Lat.%04d.%02d.bn"%(Year,Mon)
        DayPath  = sDir + "/Day.%04d.%02d.bn"%(Year,Mon)
        HourPath = sDir + "/Hour.%04d.%02d.bn"%(Year,Mon)
        MinPath  = sDir + "/Min.%04d.%02d.bn"%(Year,Mon)
        
        util.mk_dir(sDir)
        
        aprtemp  =fromfile(PrPath  , dtype="float32")  
        akutemp  =fromfile(KuPath  , dtype="float32")  
        acltemp  =fromfile(CLPath  , dtype="int32"  )  
        alstemp  =fromfile(LSPath  , dtype="float32") 
        alontemp =fromfile(LonPath , dtype="float32") 
        alattemp =fromfile(LatPath , dtype="float32") 
        adaytemp =fromfile(DayPath , dtype="int32"  ) 
        ahourtemp=fromfile(HourPath, dtype="int32"  ) 
        amintemp =fromfile(MinPath , dtype="int32"  ) 
    
        # Mask: Cloud  
        acltemp  = ma.masked_not_equal(acltemp, icl)
        #akutemp  = ma.masked_outside(akutemp, 0., 3.)

        #aprtempG  = ma.masked_less(aprtemp, akutemp*0.5)    # Good
        #aprtempG  = ma.masked_greater(aprtempG, akutemp*1.5)    # Good

        #aprtempB  = ma.masked_less(aprtemp, 10)   # Bad

        # Localtime
        aLocalTimetemp = ahourtemp.astype(float32) + amintemp.astype(float32)/60. + alontemp/360. * 24.
        aLocalTimetemp = aLocalTimetemp % 24.


        # Mask
        print "*"*50
        amask   = acltemp.mask
        amask   = amask + ma.masked_equal(akutemp,miss).mask
        amask   = amask + ma.masked_equal(aprtemp,miss).mask

        akutemp = ma.masked_where(amask, akutemp).compressed()
        aprtemp = ma.masked_where(amask, aprtemp).compressed()
        aLocalTimetemp= ma.masked_where(amask, aLocalTimetemp).compressed()


        # Concatenate
        aku   = r_[aku, akutemp]
        apr   = r_[apr, aprtemp]
        aLocalTime= r_[aLocalTime, aLocalTimetemp]


    # output
    sOut = ""
    for i, ku in enumerate(aku):
        sOut = sOut + "%s,"%(aku[i])
        sOut = sOut + "%s,"%(apr[i])
        sOut = sOut + "%s\n"%(aLocalTime[i])
    f = open("./temp.csv", "w") ; f.write(sOut); f.close()
