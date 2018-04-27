import os
import socket
import pysftp
import myfunc.util as util
from datetime import datetime, timedelta

hostname = socket.gethostname()
if  hostname =="mizu":
    #orootDir   = "/home/utsumi/mnt/wellshare/data/GPM"
    #orootDir   = "/data4/common/GPM"
    orootDir   = "/work/a01/utsumi/data/GPM"
elif hostname=="well":
    orootDir   = "/media/disk2/share/data/GPM"


iYM     = [2015,9]
eYM     = [2016,12]
#iYM     = [2017,12]
#eYM     = [2017,12]

lYM	= util.ret_lYM(iYM,eYM)

ow  = False # OverWrite
#ow  = True # OverWrite
#prj	= ["TRMM","TRMM.PR","L3A25","07"]
#prj	= ["TRMM","TRMM.PR","L2A23","07"]
#prj	= ["TRMM","TRMM.PR","L2A23","07"]
prj	= ["GPM","GPM.KuPR","L2","05"]
#prj	= ["GPM","GPM.GMI","L2","05"]
#prj	= ["GPM","GPM.DPR","L2.DPR","05"]
#prj	= ["GPM","GPM.GMI","L1B","05"]
sate	= prj[0]
sensor	= prj[1]
prdName = prj[2]
version = prj[3]

host        = 'sftp.gportal.jaxa.jp'
username    = 'nbyk.utsumi'
port        = 2051
ibaseDir    = os.path.join("/" 'standard', sate, sensor, prdName, version)
obaseDir    = os.path.join(orootDir, sensor, prdName, version)

sftp	= pysftp.Connection(host=host,username=username, port=port)
sftp.cwd(ibaseDir)
print ibaseDir
print sftp.listdir()

for Year, Mon in lYM:
    #if Mon in [11,12,1,2,3]: continue
    print Year,Mon
    srcDir = os.path.join(ibaseDir, "%04d"%(Year), "%02d"%(Mon))
    outDir = os.path.join(obaseDir, "%04d"%(Year), "%02d"%(Mon))

    try:
        sftp.cwd(srcDir)
    except IOError:
        print "No remote Directory: %s"%(srcDir)

    lfileName = sftp.listdir()
    print "srcDir=",srcDir
    util.mk_dir(outDir)

    for fileName in lfileName:
        srcPath = os.path.join(srcDir,fileName)
        outPath = os.path.join(outDir,fileName)
        stime = fileName.split("_")[2]
        etime = fileName.split("_")[3]

        YearFile  = 2000+int(stime[:2])
        MonFile   = int(stime[2:4])
        DayFile   = int(stime[4:6])
        sHourFile = int(stime[6:8])
        eHourFile = int(etime[0:2])
        # check day
        #if (DayFile<5)or(8<DayFile): continue



        '''
        #- check day --

        subDB   = 7263
        csvPath = "/home/utsumi/mnt/wellshare/ENSPR/JPLDB/csv/cmp.Simil.Emis.%05d.csv"%(subDB)
        f=open(csvPath,"r"); lines=f.readlines(); f.close()
        ltarget = []
        for line in lines[1:]:
            line = line.split(",")
            ltarget.append( map(int, line[4:7+1]))
        """
        ltarget = [
                     #[2015,5,22,23]
                     [2015,4,6,0]
                    ,[2015,6,25,14]
                    ,[2014,9,11,16]
                    ,[2015,8,23,5]
                    ,[2015,6,9,17]
                    ,[2014,12,23,6]
                    ,[2014,10,19,7]
                    ,[2015,2,20,8]
                    ,[2014,9,21,20]
                    ]
        """

        foundFlag=0
        for target in ltarget:
            tYear, tMon, tDay, tHour = target
            tDTime = datetime(tYear,tMon,tDay,tHour)
            sDTimeFile = datetime(YearFile,MonFile,DayFile,sHourFile)
            if ((sDTimeFile -timedelta(seconds=60*60) <=tDTime)&(tDTime <=sDTimeFile+timedelta(seconds=60*60*2))):
                foundFlag = foundFlag+1
            
            else:
                continue

        if foundFlag ==0: continue

#        if (YearFile==2015)&(MonFile ==5)&(DayFile==22)&(sHourFile>=23-2)&(eHourFile<=23+2): pass
#        else: continue

        '''

        #--------------
        if os.path.exists(outPath) and ow==False:
            print "Skip (ow:%5s) :: %s"%(ow, outPath)

            continue


   
        with pysftp.cd(outDir):
            try:
                sftp.get(srcPath)
            except IOError:
                print "No remote File: %s"%(srcPath)
                continue
            print outPath



