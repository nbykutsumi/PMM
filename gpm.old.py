#! /usr/bin/python
#--------------------------------------------------------------------
# PROGRAM    : gpm.py
# CREATED BY : hjkim @IIS.2015-01-14 11:52:17.992599
# MODIFED BY :
#
# USAGE      : $ ./gpm.py
#
# DESCRIPTION:
#------------------------------------------------------cf0.2@20120401


import  os,sys
from    datetime        import datetime, timedelta
from    optparse        import OptionParser
from    cf.util.LOGGER  import *

from    numpy           import array, ma

import  h5py
from    pyhdf           import SD

from    cf.util         import nearest_idx


class GPM(object):
    def __init__(self,prjName,prdVer='02'):
        '''
        prjName     : e.g.) 'GPM.KuPR.L2'
        prdVer      : e.g.) '02'
        '''

        self.baseDir        = '/tank/hjkim/GPM/'

        self.prjMain,   \
        self.prjSub,    \
        self.prdLv          = prjName.split('.')
        self.prdVer         = prdVer

        self.prdDir         = '%s/%s.%s/%s/%s'%(self.baseDir,
                                            self.prjMain, self.prjSub,
                                            self.prdLv,
                                            self.prdVer)


    def get_dtime(self,h5,Idx=None):
        '''
        h5      : h5 file object
        Idx     : indices in int/list/slice  /* default: slice(None,None,None) */
        '''

        Idx     = [Idx]                 if type(Idx) == int else Idx
        Idx     = slice(None,None,None) if Idx == None      else Idx

        Year    = h5['NS/ScanTime/Year'][Idx].astype('int')
        Month   = h5['NS/ScanTime/Month'][Idx].astype('int')
        Day     = h5['NS/ScanTime/DayOfMonth'][Idx].astype('int')
        Hour    = h5['NS/ScanTime/Hour'][Idx].astype('int')
        Minute  = h5['NS/ScanTime/Minute'][Idx].astype('int')
        Second  = h5['NS/ScanTime/Second'][Idx].astype('int')
        MicSec  = h5['NS/ScanTime/MilliSecond'][Idx].astype('int')*1000
        DTime   = [datetime(y,m,d,H,M,S)+timedelta(microseconds=int(uS))
                for y,m,d,H,M,S,uS in map(None,Year,Month,Day,Hour,Minute,Second,MicSec)]
        DTime   = array(DTime)

        return DTime


    def parse_fname(self,fName,ATTR):
        '''
        fname   : GPM HDF filename
        ATTR    : list of attributes
        '''

        fName   = fName.split('_')

        dictFunc= {'sDTime': datetime.strptime(fName[2], '%y%m%d%H%M'),
                   'eDTime': datetime.strptime(fName[2][:6]+fName[3], '%y%m%d%H%M')
                   }

        return [dictFunc[attr] for attr in ATTR]


    def search_granules(self, BBox=None, sDTime=None, eDTime=None, thresh=0.001):

#        sDTime  = self.
#        eDTime  = self. ...
#        BBox    = [[],[]]


        import time

        dictGrp = {'GPM.GMI':'S1',
                   'GPM.DPR':'NS',      # HS, MS, NS
                   'GPM.KaPR':'MS',     # HS, MS
                   'GPM.KuPR':'NS',}

        grpCode = dictGrp['.'.join( [self.prjMain, self.prjSub] )]

        srcDIR  = [ os.path.join( self.prdDir, '%i/%02d'%(y,m) )
                                                        for y in range(sDTime.year, eDTime.year+1)
                                                        for m in range(1,13)]
        srcDIR  = filter(os.path.exists, srcDIR)

        targetPATH  = []
        for srcDir in srcDIR:
            for fName in os.listdir( srcDir ):

                sDTime_orb, eDTime_orb  = self.parse_fname( fName, ['sDTime','eDTime'])

                # screening by filename
                if (sDTime <= sDTime_orb <= eDTime) or (sDTime <= eDTime_orb <= eDTime):
                    srcPath     = os.path.join(srcDir,fName)
                    print srcPath

                    try:
                      h5      = h5py.File( srcPath, 'r' )

                      Lat     = h5['%s/Latitude'%grpCode]
                      Lon     = h5['%s/Longitude'%grpCode]
                    except IOError:
                      print "Fail!! in reading file="
                      print srcPath
                      continue
  
                    #---- By NU ------
                    if len(Lat)<24:
                      print "Blank!!"
                      print srcPath
                      continue
                    #-----------------

                    try:
                      mask    = ma.masked_outside( Lat[:,24], BBox[0][0], BBox[1][0] ).mask
                      sect    = ma.array( Lon[:,24], mask=mask).compressed()
                    except:
                      continue

                    for s in sect:
                        if BBox[0][1] <= s <= BBox[1][1]:
                            targetPATH.append( srcPath )
                            break

                    h5.close()

        return targetPATH


@ETA
def main(args,opts):
    print args
    print opts

    prjName = 'GPM.KuPR.L2'
    prdVer  = '02'

    BBox    = [[20.0, 118.0], [48.0, 150.0]]
    sDTime  = datetime( 2014,4,15,12 )
    sDTime  = datetime( 2014,4,30,3 )
    eDTime  = datetime( 2014,5,5,8 )

    gpm     = GPM(prjName, prdVer)
    print gpm.search_granules(BBox, sDTime, eDTime)


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

#    LOG     = LOGGER()
    main(args,options)


