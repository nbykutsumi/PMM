#! /usr/bin/python
#--------------------------------------------------------------------
# PROGRAM    : get.GPM.py
# CREATED BY : hjkim @IIS.2014-07-23 15:37:53.314375
# MODIFED BY : n.utsumi @IIS.2014-11-22
#
# USAGE      : $ python ./get.GPM.py
#
# DESCRIPTION:
#------------------------------------------------------cf0.2@20120401

import  os,sys,re
from    optparse        import OptionParser
#from    cf.util.LOGGER  import *

import  pysftp
from    datetime        import datetime, timedelta

#--------------------------------------------------------
#outDir   = "/media/disk2/data"
outDir   = "/home/utsumi/mnt/wellshare/data/TRMM"
#--------------------------------------------------------
class GPortal(object):

    def __init__(self, prjName):
        '''
        prjName     in ['']
        '''

        host        = 'sftp.gportal.jaxa.jp'
        username    = 'nbyk.utsumi'
        port        = 2051

        baseDir     = os.path.join( 'standard', prjName.split('.')[0], prjName )

        self.prjName= prjName
        self.version= {'TRMM':'07',
                       'GPM':'02',
                       }[ prjName.split('.')[0] ]

        self.dtFunc = {'TRMM': lambda s: re.findall(r'\d{8}', s.split('/')[-1])[0],
                       'GPM' : lambda s: ('20' + s.split('/')[-1].split('_')[2] + '01')[:8]
                       }[ prjName.split('.')[0] ]

        self.sftp   = pysftp.Connection(host=host,username=username, port=port)
        self.sftp.cwd(baseDir)

        print "in Gportal, 1"
        self.dirs, self.files   = self.sep_files_dirs()
        print "in Gportal, 2"

        #print self.sftp.listdir()

    def __getitem__(self,k):

        print "__getitem__"
        self.sftp.cwd(k)
        self.sftp.cwd(self.version)
        self.dirs, self.files   = self.sep_files_dirs()

        self.subPrj     = k

        return self


    def sep_files_dirs(self):
        FILE    = []
        DIR     = []

        print "sep_files_dirs"
        print "AA", self.sftp.listdir()
        for fd in self.sftp.listdir():
            if self.sftp.isdir(fd):
                DIR.append(fd)

            else:
                FILE.append(fd)

        return DIR,FILE


    #def get_dtime(self, sDTime, eDTime, outDir='./', ow=False):
    def get_dtime(self, sDTime, eDTime, outDir=outDir, ow=False):
        srcDIR  = ['%i/%02d/'%(y,m) for y in range(sDTime.year, eDTime.year+1) for m in range(1,13)]
        srcDIR  = srcDIR[sDTime.month-1 : None if eDTime.month==12 else eDTime.month-12]      # tentaive patch to resolve skipping December
        #srcDIR  = srcDIR[sDTime.month-1 : eDTime.month-12]

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

            saveDir = os.path.join(outDir, self.prjName, self.subPrj, subDir)

            if not os.path.exists( saveDir ): os.makedirs( saveDir )

            outPath     = os.path.join(outDir, self.prjName, self.subPrj, srcPath)

            if os.path.exists(outPath) and ow == False:
                print 'Skip (ow:%5s) ::'%ow, outPath
                continue

            with pysftp.cd(saveDir):
                self.sftp.get(srcPath)
                print 'Save (ow:%5s) ->'%ow, outPath




#@ETA
def main(args,opts):
    print args
    print opts

    print "main"
    year    = int(args[0])
    month   = int(args[1])

#    gpm     = GPortal('GPM.COMB')['L3']
#    gpm     = GPortal('GPM.DPR')['L3.Daily']
#    gpm.get_dtime( datetime(2014,3,1), datetime(2014,3,31), outDir='../' )

    for d in range(1,32):
        print "day ==",d
        print "call GPortal"
#        gpm     = GPortal('GPM.KaPR')['L2']
#        gpm     = GPortal('GPM.KuPR')['L2']
#        gpm     = GPortal('GPM.DPR')['L2.DPR']
#        gpm     = GPortal('GPM.COMB')['L2']
#        gpm     = GPortal('TRMM.COMB')['L2B31']
#        gpm     = GPortal('TRMM.PR')['L2A21']
#        gpm     = GPortal('TRMM.PR')['L2A25']
        gpm     = GPortal('TRMM.PR')['L3A25']
#        gpm     = GPortal('TRMM.TMI')['L2A12']
        #gpm.get_dtime( datetime(year,month,d), datetime(year,month,d), outDir='../' )
        print "call gpm.get_dtime"
        gpm.get_dtime( datetime(year,month,d), datetime(year,month,d), outDir=outDir)
        print "call gpm.sftp.close()"
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

    #LOG     = LOGGER()
    print "0000"
    main(args,options)


