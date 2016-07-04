#! /usr/bin/python
#--------------------------------------------------------------------
# PROGRAM    : comp.GPM_TRMM.py
# CREATED BY : hjkim @IIS.2014-07-28 11:44:32.391196
# MODIFED BY :
#
# USAGE      : $ ./comp.GPM_TRMM.py
#
# DESCRIPTION:
#------------------------------------------------------cf0.2@20120401


import  os,sys,re
from    datetime        import *
from    itertools       import imap
from    optparse        import OptionParser
#from    cf.util.LOGGER  import *

from    numpy           import array,concatenate,searchsorted,sqrt

import  h5py
from    pyhdf           import SD

from    cf.TimeSeries   import TSData


class GPM(object):
    def __init__(self, prjName):
        #self.baseDir    = '/tank/hjkim/GPM'
        self.baseDir    = '/mnt/mizu.tank/hjkim/GPM'
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


def get_cross(aTRMM, dimTR, aGPM, dimGPM):

    '''
    for ts in TSData(dimGPM[0],
                     dimGPM[0][0],
                     delT=timedelta(microseconds=70000),
                     iterT=timedelta(seconds=70)):
        print ts.data
    '''

    sDTime  = max(dimTR[0][0],  dimGPM[0][0])
    eDTime  = min(dimTR[0][-1], dimGPM[0][-1])

    sIdxTR  = searchsorted(dimTR[0],sDTime)
    eIdxTR  = searchsorted(dimTR[0],eDTime)
    sIdxGPM = searchsorted(dimGPM[0],sDTime)
    eIdxGPM = searchsorted(dimGPM[0],eDTime)

    trmm    = aTRMM[sIdxTR:eIdxTR]
    gpm     = aGPM[sIdxGPM:eIdxGPM]

    dtimeT  = dimTR[0][sIdxTR:eIdxTR]
    latT    = dimTR[1][sIdxTR:eIdxTR,25]
    lonT    = dimTR[2][sIdxTR:eIdxTR,25]

    dtimeG  = dimGPM[0][sIdxGPM:eIdxGPM]
    latG    = dimGPM[1][sIdxGPM:eIdxGPM,25]
    lonG    = dimGPM[2][sIdxGPM:eIdxGPM,25]



    '''
    stT     = 1000
    stG     = 856
    for t,g,lat0,lat1,lon0,lon1 in map(None,
                                       dtimeT[::stT], dtimeG[::stG],
                                       latT[::stT], latG[::stT],
                                       lonT[::stT], lonG[::stT]):
        print t,g, t-g, lat0, lon0, lat1,lon1#sqrt((lon1-lon0)**2+(lat1-lat0)**2)

    '''
    tsGPM   = TSData(array([latG,lonG]).T,
                     dtimeG[0],
                     delT=timedelta(microseconds=700000),
                     iterT=timedelta(seconds=42*8))

    tsTRMM  = TSData(array([latT,lonT]).T,
                     dtimeT[0],
                     delT=timedelta(microseconds=600000),
                     iterT=timedelta(seconds=42*8))

    outTRMM = []
    outGPM  = []
    for tsT, tsG  in imap(None, tsTRMM, tsGPM):
        distDeg     = sqrt(sum((tsT.data[0]-tsG.data[0])**2))

        #print tsT.data[0], tsG.data[0], '%5.1f'%distDeg

        if distDeg < 5:
            print tsT.DTIME[0], tsT.data[0]
            print tsG.DTIME[0], tsG.data[0]
            print '***',distDeg

#            outTRMM.append(tsT)
#            outGPM.append(tsG)
            crdT, crdG  = get_superimposed_crd(tsT, tsG)

            if crdT == None: continue

            idxT    = tsTRMM.data.tolist().index(crdT.tolist())
            idxG    = tsGPM.data.tolist().index(crdG.tolist())

            outTRMM.append( [dtimeT[idxT-100:idxT+100],
                             dimTR[1][sIdxTR:eIdxTR][idxT-100:idxT+100],
                             dimTR[2][sIdxTR:eIdxTR][idxT-100:idxT+100],
                             trmm[idxT-100:idxT+100]] )

            outGPM.append( [dtimeG[idxG-100:idxG+100],
                             dimGPM[1][sIdxGPM:eIdxGPM][idxG-100:idxG+100],
                             dimGPM[2][sIdxGPM:eIdxGPM][idxG-100:idxG+100],
                             gpm[idxG-100:idxG+100]] )

            print '+'*100
            print outTRMM[-1][0][100], outGPM[-1][0][100]
            print '+'*100

    return outTRMM, outGPM


def get_superimposed_crd(tsTRMM, tsGPM):

    threshDis   = 1.0

    get_closest_idx = lambda crd,CRD: ((CRD-crd)**2).sum(1).argmin()


    IdxGPM  = [get_closest_idx(a, tsGPM.data)
                        for a in tsTRMM.data]

    DIS     = [sqrt(((tsGPM.data[i]-a)**2).sum())
                        for a,i in map(None, tsTRMM.data, IdxGPM)]

    idxTRMM = array(DIS).argmin()
    idxGPM  = IdxGPM[idxTRMM]

    dtTRMM, crdTRMM = tsTRMM.DTIME[idxTRMM], tsTRMM.data[idxTRMM]
    dtGPM,  crdGPM  = tsGPM.DTIME[idxGPM], tsGPM.data[idxGPM]

    dtDiff          = dtTRMM-dtGPM
    crdDiff         = sqrt( ((crdTRMM-crdGPM)**2).sum() )

    if crdDiff > threshDis:
        return None, None

    print '='*80
    print '%22s | %22s | %20s'%('TRMM','GPM','diff')
    print '%22s | %22s | %16s'%(dtTRMM, dtGPM, dtDiff)
    print '%22s | %22s | %6.3f'%(crdTRMM, crdGPM, crdDiff)
    print '='*80


    return crdTRMM, crdGPM


def mapping(CrossT,CrossG):

    from pylab import *
    from mpl_toolkits.basemap import Basemap
    '''
    fig=figure(figsize=(12,6))
    ax=fig.add_subplot(111)
    M = Basemap(resolution='c')
    '''

    for crossT, crossG in map(None, CrossT, CrossG):

        dtT, latT, lonT, aT = crossT
        dtG, latG, lonG, aG = crossG

        #M.pcolor(lonT,latT,aT)
        #M.pcolor(lonG,latG,aG)
        fig=figure(figsize=(3,9))
        subplot(311);pcolor(lonT,latT,aT);pcolor(lonG,latG,aG);colorbar()
        subplot(312);pcolor(lonT,latT,aT);colorbar()
        subplot(313);pcolor(lonG,latG,aG);colorbar()
    '''
    M.scatter(lonT,latT,2,color='dodgerblue')
    M.scatter(lonG,latG,2,color='lime')
    M.scatter(tsT.data[:,1],tsT.data[:,0],30,color='navy')
    M.scatter(tsG.data[:,1],tsG.data[:,0],30,color='darkgreen')
    '''

    #M.drawcoastlines()
    show()

    return


#@ETA
def main(args,opts):
    print args
    print opts

    #year        = int(args[0])
    #month       = int(args[1])

    year,month  = 2014, 5

    # 2014/05/14
    # 2014/07/07 #3, 6, 7
    sDTime      = datetime(year,month,14)
    eDTime      = datetime(year,month,14)
    delT        = timedelta(days=1)

    trmm        = TRMM('TRMM.PR')['L2A25']
    gpm         = GPM('GPM.DPR')['L2.DPR']

    while sDTime <= eDTime:
        aTRMM, dimTR= trmm.get_dtime('e_SurfRain',
                                    sDTime, sDTime)

        #print aTRMM.shape, dimTR[0].shape

        aGPM, dimGPM= gpm.get_dtime('NS/SLV/precipRateESurface',
                                    sDTime, sDTime)

        #print aGPM.shape, dimGPM[0].shape

        CrossT, CrossG  = get_cross(aTRMM, dimTR, aGPM, dimGPM)

        mapping(CrossT[0:], CrossG[0:])

        sDTime += delT


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
    main(args,options)


