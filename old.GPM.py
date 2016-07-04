#! /usr/bin/python
#--------------------------------------------------------------------
# PROGRAM    : GPM.py
# CREATED BY : hjkim @IIS.2014-07-23 15:37:53.314375
# MODIFED BY :
#
# USAGE      : $ ./GPM.py
#
# DESCRIPTION:
#------------------------------------------------------cf0.2@20120401


import  os,sys,re
from    optparse        import OptionParser
from    cf.util.LOGGER  import *

import  h5py
from    pyhdf           import SD

from    numpy           import array,concatenate,searchsorted,sqrt

import  pysftp
from    datetime        import datetime, timedelta


class GPortal(object):

    def __init__(self, prjName):
        '''
        prjName     in ['']
        '''

        host        = 'sftp.gportal.jaxa.jp'
        username    = 'AP05708'
        port        = 2051

        baseDir     = os.path.join( 'standard', prjName.split('.')[0], prjName )

        self.prjName= prjName
        self.version= {'TRMM':'07',
                       'GPM':'03',
#                       'GPM':'02',
                       }[ prjName.split('.')[0] ]

        self.dtFunc = {'TRMM': lambda s: re.findall(r'\d{8}', s.split('/')[-1])[0],
                       'GPM' : lambda s: ('20' + s.split('/')[-1].split('_')[2] + '01')[:8]
                       }[ prjName.split('.')[0] ]

        self.sftp   = pysftp.Connection(host=host,username=username, port=port)
        self.sftp.cwd(baseDir)

        self.dirs, self.files   = self.sep_files_dirs()


    def __getitem__(self,k):

        self.sftp.cwd(k)
        self.sftp.cwd(self.version)
        self.dirs, self.files   = self.sep_files_dirs()

        self.subPrj     = k

        return self


    def sep_files_dirs(self):
        FILE    = []
        DIR     = []

        for fd in self.sftp.listdir():
            if self.sftp.isdir(fd):
                DIR.append(fd)

            else:
                FILE.append(fd)

        return DIR,FILE


    def get_dtime(self, sDTime, eDTime, outDir='./', ow=False):
        srcDIR  = ['%i/%02d/'%(y,m) for y in range(sDTime.year, eDTime.year+1) for m in range(1,13)]
        srcDIR  = srcDIR[sDTime.month-1 : None if eDTime.month==12 else eDTime.month-12]      # tentaive patch to resolve skipping December
#        srcDIR  = srcDIR[sDTime.month-1 : eDTime.month-12]

        srcPATH = []

        for srcDir in srcDIR:
            if not self.sftp.exists(srcDir):
                print 'Not Exists!', srcDir
                continue

            with self.sftp.cd(srcDir):
                srcPATH.extend([os.path.join(srcDir,srcFName) for srcFName in self.sftp.listdir()])

        DTIME   = [datetime.strptime( self.dtFunc(srcPath),'%Y%m%d') for srcPath in srcPATH]

        for dtime, srcPath in map(None, DTIME, srcPATH):
            if sDTime <= dtime and eDTime >= dtime: pass
            else:                                   continue

            subDir  = dtime.strftime('%Y/%m')

            saveDir = os.path.join(outDir, self.prjName, self.subPrj, self.version, subDir)

            if not os.path.exists( saveDir ): os.makedirs( saveDir )

            outPath     = os.path.join(outDir, self.prjName, self.subPrj, self.version, srcPath)

            if os.path.exists(outPath) and ow == False:
                print 'Skip (ow:%5s) ::'%ow, outPath
                continue

            with pysftp.cd(saveDir):
                self.sftp.get(srcPath)
                print 'Save (ow:%5s) ->'%ow, outPath



class GPM(object):
    def __init__(self, prjName):
        self.baseDir    = '/tank/hjkim/GPM'
        self.prjName    = prjName

        self.readFunc   = {'GPM':self.read_hdf5,
                           'TRMM':self.read_hdf4,
                           }[ prjName.split('.')[0] ]

        self.dtFunc = {'TRMM': lambda s: re.findall(r'\d{8}', s.split('/')[-1])[0],
                       'GPM' : lambda s: ('20' + s.split('/')[-1].split('_')[2] + '01')[:8]
                       }[ prjName.split('.')[0] ]


    def __getitem__(self,k):

        self.srcDir     = os.path.join(self.baseDir,
                                       self.prjName,
                                       k)
        return self

    def read_hdf5(self, dtime, srcPath, varName, outMeta=True):

        print dtime,srcPath

        h5  = h5py.File(srcPath)

        aOut    = h5[varName][:]

        if outMeta:

            Lon     = h5['NS/Longitude'][:]
            Lat     = h5['NS/Latitude'][:]

            Year    = h5['NS/ScanTime/Year'][:].astype('int')
            Month   = h5['NS/ScanTime/Month'][:].astype('int')
            Day     = h5['NS/ScanTime/DayOfMonth'][:].astype('int')
            Hour    = h5['NS/ScanTime/Hour'][:].astype('int')
            Minute  = h5['NS/ScanTime/Minute'][:].astype('int')
            Second  = h5['NS/ScanTime/Second'][:].astype('int')
            MicSec  = h5['NS/ScanTime/MilliSecond'][:].astype('int')*1000
            DTime   = [datetime(y,m,d,H,M,S)+timedelta(microseconds=int(uS))
                    for y,m,d,H,M,S,uS in map(None,Year,Month,Day,Hour,Minute,Second,MicSec)]
            DTime   = array(DTime)

            h5.close()
            return aOut, [DTime, Lat, Lon]

        else:
            h5.close()
            return aOut, None

    def read_hdf4(self, dtime, srcPath, varName, outMeta=True):
        print dtime,srcPath

        h4  = SD.SD(srcPath)

        aOut    = h4.select(varName)[:]

        if outMeta:

            Lon     = h4.select('Longitude')[:]
            Lat     = h4.select('Latitude')[:]

            Year    = h4.select('Year')[:]
            Month   = h4.select('Month')[:]
            Day     = h4.select('DayOfMonth')[:]
            Hour    = h4.select('Hour')[:]
            Minute  = h4.select('Minute')[:]

            Sec     = h4.select('Second')[:]
            uSec    = h4.select('MilliSecond')[:]*1000

            '''
            scanTime= h4.select('scanTime_sec')[:]
            uSec    = ( (scanTime-scanTime.astype('int'))*1.E6 ).astype('int')
            '''

            DTime   = [datetime(y,m,d,H,M,S)+timedelta(microseconds=int(uS))
                    for y,m,d,H,M,S,uS in map(None,Year,Month,Day,Hour,Minute,Sec,uSec)]
            DTime   = array(DTime)

            h4.end()
            return aOut, [DTime, Lat, Lon]

        else:
            h4.end()
            return aOut, None


    def get_dtime(self, varName, sDTime, eDTime):
        srcDIR  = [os.path.join(self.srcDir, '%i/%02d'%(y,m))
                                for y in range(sDTime.year,eDTime.year+1)
                                    for m in range(1,13)]

        srcDIR  = srcDIR[sDTime.month-1 : eDTime.month-12]

        srcPATH = []

        for srcDir in srcDIR:
            if not os.path.exists(srcDir):
                print 'Not Exists!', srcDir
                continue

            srcPATH.extend([os.path.join(srcDir, srcFName)
                                for srcFName in sorted(os.listdir(srcDir))])

        DTIME   = [datetime.strptime( self.dtFunc(srcPath),'%Y%m%d') for srcPath in srcPATH]

        aOut    = []

        for dtime, srcPath in map(None, DTIME, srcPATH):
            if sDTime <= dtime and eDTime >= dtime: pass
            else:                                   continue

            aOut.append( self.readFunc(dtime, srcPath, varName) )

        aOut, dim   = zip(*aOut)

        if dim[0] != None:
            DTime, Lat, Lon = zip(*dim)

            dim     = [concatenate(DTime),
                       concatenate(Lat),
                       concatenate(Lon)]

        return concatenate( aOut ), dim


class TRMM(GPM):
    pass



@ETA
def main(args,opts):
    print args
    print opts

    year    = int(args[0])
    month   = int(args[1])

    outDir  = args[2]

#    gpm     = GPortal('GPM.COMB')['L3']
#    gpm     = GPortal('GPM.DPR')['L3.Daily']
#    gpm.get_dtime( datetime(2014,3,1), datetime(2014,3,31), outDir='../' )

    for d in range(1,32):
#        gpm     = GPortal('GPM.KaPR')['L2']
#        gpm     = GPortal('GPM.KuPR')['L2']
#        gpm     = GPortal('GPM.DPR')['L2.DPR']
        gpm     = GPortal('GPM.GMI')['L2']
#        gpm     = GPortal('GPM.COMB')['L2']
#        gpm     = GPortal('TRMM.COMB')['L2B31']
#        gpm     = GPortal('TRMM.PR')['L2A21']
#        gpm     = GPortal('TRMM.PR')['L2A25']
#        gpm     = GPortal('TRMM.TMI')['L2A12']
        gpm.get_dtime( datetime(year,month,d), datetime(year,month,d), outDir=outDir )
        gpm.sftp.close()

#    gpm     = GPortal('TRMM.COMB')['L3B42']
#    gpm     = GPortal('TRMM.COMB')['L3B43']
#    gpm     = GPortal('TRMM.COMB')['L3B31']
#    gpm.get_dtime( datetime(year,7,1), datetime(year,7,31), outDir='../' )

    return


if __name__=='__main__':
    usage   = 'usage: %prog [options] arg'
    version = '%prog 1.0'

    parser  = OptionParser(usage=usage,version=version)

#    parser.add_option('-r','--rescan',action='store_true',dest='rescan',
#                      help='rescan all directory to find missing file')

    (options,args)  = parser.parse_args()

#    if len(args) == 0:
#        parser.print_help()
#    else:
#        main(args,options)

    LOG     = LOGGER()
    main(args,options)


