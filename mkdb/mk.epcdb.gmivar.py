from numpy import *
import calendar
import myfunc.util as util
import numpy as np
from collections import deque
import sys, os, glob
import gc
import h5py
import socket


prod     = '1C'
verGMI   = '05'
subverGMI= 'A'
fullverGMI='%s%s'%(verGMI, subverGMI)

iYM = [2017,1]
eYM = [2017,12]
lYM = util.ret_lYM(iYM,eYM)
#lepcid_range = [[0,2500],[2500,5000],[5000,7500],[7500,10000],[10000,12500],[12500,15624]]  # 25*25*25 = 15625
#lepcid_range = [[0,5000],[5000,10000],[10000,15000],[15000,20000],[20000,25000]]  # 29*29*29 = 24389
lepcid_range = [[0,2500],[2500,5000],[5000,7500],[7500,10000],[10000,12500],[12500,15000],[15000,17500],[17500,20000],[20000,22500],[22500,25000]]  # 29*29*29 = 24389


cx  = 110  # GMI center angle bin (py-idx)
cw  = 15    # extract this width around center
w   = int(cw/2)
worg= 221  # GMI total angle bins

#lvar = [['gmi','S1/Latitude'],['gmi','S1/Longitude'],['gprof','S1/surfaceTypeIndex'],['gprof','S1/surfacePrecipitation'], ['gmi','S1/Tc']]
#lvar = [['gmi','S1/Latitude'],['gmi','S1/Longitude']]
#lvar = [['gmi','S1/Latitude']]
#lvar = [['gmi','S1/Longitude']]
#lvar = [['gmi','epc'],['gmi','S1/Tc'],['gmi','S1/pYXpmw'],['gprof','S1/surfaceTypeIndex']]
#lvar = [['gmi','epc'],['gmi','Tc'],['gmi','S1/pYXpmw']]
#lvar = [['gprof','S1/surfacePrecipitation']]
#lvar = [['gmi','S1/Longitude']]
#lvar = [['gmi','S1/Latitude']]
#lvar = [['gmi','S1/pYXpmw'],['gprof','S1/surfaceTypeIndex']]
#lvar = [['gprof','S1/surfaceTypeIndex']]
#lvar = [['gprof','S1/precipTotWaterCont']]
lvar = [['gprof','S1/qualityFlag']]
#lvar = [['gmi','epc']]
#lvar = [['gmi','S1/pYXpmw'],['gmi','S1/gNum'],['gmi','S1/Tc'],['gmi','S1/ScanTime/Year'],['gmi','S1/ScanTime/mdhms']]
#lvar = [['gmi','S1/pYXpmw'],['gmi','S1/gNum'],['gmi','S1/ScanTime/Year'],['gmi','S1/ScanTime/mdhms']]
#lvar = [['gmi','S1/ScanTime/Year'],['gmi','S1/ScanTime/mdhms']]
#lvar = [['gmi','S1/gNum']]
#lvar = [['gmi','S1/ScanTime/Year']]
#lvar = [['gmi','S1/ScanTime/mdhms']]
#lvar = [['gmi','S1/LatLon']]

'''
int8 : -128 ~ +127
int16: -32768 ~ +32767
int32: -2147483648 ~ +2147483647
'''

dattype={
 'S1/Latitude' :'float32'
,'S1/Longitude':'float32'
,'S1/ScanTime/Year':'int16'
,'S1/ScanTime/mdhms':'int8'
,'S1/qualityFlag':'int8'
,'S1/surfaceTypeIndex':'int32'
,'S1/surfacePrecipitation':'float32'
,'S1/precipTotWaterCont':'int16'
,'S1/pYXpmw': 'int16'
,'S1/gNum': 'int16'
,'epc': 'float32'
,'Tc' : 'float32'
}

dnvect ={
 'S1/Latitude' :1
,'S1/Longitude':1
,'S1/ScanTime/Year':1
,'S1/ScanTime/mdhms':5
,'S1/qualityFlag':1
,'S1/surfaceTypeIndex':1
,'S1/surfacePrecipitation':1
,'S1/precipTotWaterCont':28
,'S1/pYXpmw': 2
,'S1/gNum': 1
,'epc': 13
,'Tc' : 13
}

#-- if vector size=11 --
maxmem = 3*1000*1000*1000 # N GB
dmaxrec= {}
for (prodName, var) in lvar:
    nvect = dnvect[var]
    dmaxrec[var] = int(maxmem / (4*nvect*10000))

#-----------------------
myhost = socket.gethostname()
if myhost =='shui':
    listDir  = '/work/hk01/utsumi/PMM/EPCDB/list'
    extractidDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V%s/%s.ABp%03d-%03d.GMI.%s'%(fullverGMI, 'S1', cx-w, cx+w, 'epcid-s1')
    gmibaseDir   = '/work/hk01/PMM/NASA/GPM.GMI/%s/V%s'%(prod,verGMI)
    gprofbaseDir = '/work/hk01/PMM/NASA/GPM.GMI/2A/V05'
    outbaseDir   = '/work/hk01/utsumi/PMM/EPCDB/GMI.V%s.%s.ABp%03d-%03d'%(fullverGMI, 'S1', cx-w, cx+w)
elif myhost =='well':
    listDir      = '/home/utsumi/mnt/lab_tank/utsumi/PMM/EPCDB/list'
    extractidDir = '/home/utsumi/mnt/lab_tank/utsumi/PMM/MATCH.GMI.V%s/%s.ABp%03d-%03d.GMI.%s'%(fullverGMI, 'S1', cx-w, cx+w, 'epcid-s1')
    gmibaseDir   = '/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.GMI/%s/V%s'%(prod,verGMI)
    gprofbaseDir = '/home/utsumi/mnt/lab_work/hk01/PMM/NASA/GPM.GMI/2A/V05'
    outbaseDir   = '/home/utsumi/mnt/lab_tank/utsumi/PMM/EPCDB/GMI.V%s.%s.ABp%03d-%03d'%(fullverGMI, 'S1', cx-w, cx+w)


#-- Functions --
def csv2list(srcPath):
    f=open(srcPath,'r'); lines=f.readlines(); f.close()
    lout = []
    for line in lines:
        line = line.strip().split(',')
        lout.append(line)
    return lout

def ret_aprof(a4clusterProf, a2tIndex, a3profNum, a3profScale, lspecies=[0,2,3,4]):
    nh = 28
    ny,nx = a2tIndex.shape
    a4out = empty([len(lspecies),ny, nx, nh], dtype='float32')

    for i,species in enumerate(lspecies):
        a1profScale= a3profScale[:,:,species].flatten()
        a1profNum  = a3profNum[:,:,species].flatten()
        a1tIndex   = a2tIndex.flatten()

        #-- Handle non-precipitation pixels --
        a1flag = ma.masked_equal(a1profNum, 0).mask
        a1profNum[a1flag] = 1
        a1tIndex[a1flag] = 1

        a2prof = a1profScale.reshape(-1,1) * a4clusterProf[a1profNum-1,:, a1tIndex-1, species]
        a2prof[a1flag,:] = 0.0
        a4out[i] = a2prof.reshape(ny,nx,nh)

    return a4out


#---------------
for Year,Mon in lYM:
    for (prodName, var) in lvar:
        for epcid_range in lepcid_range:

            #if ((Mon==4)&(epcid_range==[0,2500])): continue  # test

            epcid_min, epcid_max = epcid_range
            maxrec   = dmaxrec[var]
            grp      = '/'.join(var.split('/')[:-1])
            varName  = var.split('/')[-1]
            listPath = listDir + '/list.shuffle.%s.V%s.%04d%02d.csv'%(prod,verGMI,Year,Mon)
            lobt = csv2list(listPath)
        
            dstack = {}
            dnum   = {}
            for epcid in range(epcid_min,epcid_max): 
                dstack[epcid] = deque([])
                dnum  [epcid] = 0
    
            ##-- test --
            #lobt = lobt[:3]
            ##----------    
            for (obtnum, Year,Mon,Day,time0, time1) in lobt:   
                Year,Mon,Day = map(int, [Year,Mon,Day])
                print var, obtnum,Year,Mon,Day
                #-- read EPC-id --
                extractDir = extractidDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
                ssearch = extractDir + '/epcid-s1.%s.npy'%(obtnum)
                epcidPath = glob.glob(ssearch)[0]
                a2epcid   = np.load(epcidPath)
                a1epcid   = a2epcid.flatten()
    
                #-- Search GMI granules --
                if  (prodName=='gmi')&(varName !='epc'):
                    srcDir = gmibaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
                    ssearch= srcDir + '/*.%s.V%s.HDF5'%(obtnum,fullverGMI)
                    srcPath= glob.glob(ssearch)[0]
                    print srcPath
                elif (prodName=='gprof')&(varName !='epc'):
                    srcDir = gprofbaseDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
                    ssearch= srcDir + '/*.%s.V05?.HDF5'%(obtnum)
                    srcPath= glob.glob(ssearch)[0]
                    print srcPath
    
                #-- Read GMI variables --
                if varName  =='epc':
                    NEM     = 12
                    extractepcDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp103-117.GMI.epc-s1'
                    extractDir = extractepcDir + '/%04d/%02d/%02d'%(Year,Mon,Day)
                    ssearch = extractDir + '/epc-s1.%s.npy'%(obtnum)
                    epcPath = glob.glob(ssearch)[0]
                    Dat     = np.load(epcPath).reshape(-1,NEM)
 

                elif varName=='Tc':
                    NTB1 = 9
                    NTB2 = 4
                    extractDir1 = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp103-117.GMI.Tc/%04d/%02d/%02d'%(Year,Mon,Day)
                    extractDir2 = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A/S1.ABp103-117.GMI.TcS2/%04d/%02d/%02d'%(Year,Mon,Day)
                    datPath1= extractDir1 + '/Tc.%s.npy'%(obtnum)
                    datPath2= extractDir2 + '/TcS2.1.%s.npy'%(obtnum)
                    Dat1    = np.load(datPath1).reshape(-1,NTB1)
                    Dat2    = np.load(datPath2).reshape(-1,NTB2)
                    Dat = np.concatenate([Dat1,Dat2],axis=1) 
    
                elif varName=='mdhms':
                    with h5py.File(srcPath) as h:
                        mm   = h[grp+'/Month'][:]
                        dd   = h[grp+'/DayOfMonth'][:]
                        hh   = h[grp+'/Hour'][:]
                        mn   = h[grp+'/Minute'][:]
                        ss   = h[grp+'/Second'][:]
    
                    nytmp,nxtmp = a2epcid.shape 
                    Dat = np.zeros([nytmp,nxtmp,5],int8)
                    for xtmp in range(nxtmp):
                        Dat[:,xtmp,0] = mm
                        Dat[:,xtmp,1] = dd
                        Dat[:,xtmp,2] = hh
                        Dat[:,xtmp,3] = mn
                        Dat[:,xtmp,4] = ss
                    Dat = Dat.reshape(-1,5)
    
                elif varName=='pYXpmw':
                    nytmp,nxtmp = a2epcid.shape
                    a1y = range(nytmp)
                    a1x = range(worg)
                    X,Y = np.meshgrid(a1x,a1y)
                    X   = X[:,cx-w:cx+w+1]
                    Y   = Y[:,cx-w:cx+w+1]

                    Dat = zeros([nytmp,nxtmp,2],int16)
                    Dat[:,:,0] = Y
                    Dat[:,:,1] = X
                    Dat = Dat.reshape(-1,2).astype(int16)
    
                elif varName=='gNum':
                    nytmp,nxtmp = a2epcid.shape
                    Dat = empty(nytmp*nxtmp,int16)
                    Dat[:] = obtnum

                elif varName=='Year':
                    with h5py.File(srcPath) as h:
                        a1yyyy= h[grp+'/Year'][:]
   
                    nytmp,nxtmp = a2epcid.shape
                    Dat = empty([nytmp,nxtmp],int16)
                    for xtmp in range(nxtmp):
                        Dat[:,xtmp] = a1yyyy
                    Dat = Dat.flatten()

                elif varName in ['qualityFlag','surfaceTypeIndex','surfacePrecipitation','Latitude','Longitude']:
                    with h5py.File(srcPath) as h:
                        Dat = h[var][:]

                    if len(Dat.shape)==2:
                        Dat = Dat[:,cx-w:cx+w+1].flatten()
                    elif len(Dat.shape)==3:
                        nyTmp,nxTmp,nzTmp = Dat.shape
                        Dat = Dat[:,cx-w:cx+w+1,:].reshape(-1,nzTmp)
                    else:
                        print 'check Dat.shape',Dat.shape
                        sys.exit()

                elif varName in ['precipTotWaterCont']:
                    with h5py.File(srcPath) as h:
                        a4clusterProf= h['GprofDHeadr/clusterProfiles'][:]  # (profNumber, nlev, nT, nspecies) = (80, 28, 12, 5)
                        a2tIndex   = h['S1/profileTemp2mIndex']     [:,cx-w:cx+w+1] # (Y,X)  zero=missing value?
                        a3profNum  = h['S1/profileNumber']          [:,cx-w:cx+w+1,:] # (Y,X, nspecies)
                        a3profScale= h['S1/profileScale']           [:,cx-w:cx+w+1,:]  # (Y,X, nspecies)


                    a4profg = ret_aprof(a4clusterProf, a2tIndex, a3profNum, a3profScale)
                    a3profg = ma.masked_less(a4profg,0).sum(axis=0)[:,:,::-1].filled(-9.999) # Bottom to top --> Top to bottom
                    a3profg = ma.masked_greater(a3profg, 32.7).filled(32.7)  # int16: -32767 ~ +32767
                    a3profg = (ma.masked_less(a3profg, 0) *1000).filled(-9999).astype('int16')  # Scaled by 1000, Convert to int16
                    nztmp = a3profg.shape[2]
                    Dat = a3profg.reshape(-1,nztmp)

                else:
                    print 'check varName', varName
                    print 'exit' 
                    sys.exit()
                     
                
                #-- stack to epc-classifiled database ---
                lepcidset = sort(list(set(a1epcid)))
                lepcidset = ma.masked_less(lepcidset, epcid_min)
                lepcidset = ma.masked_greater_equal(lepcidset, epcid_max)
                lepcidset = lepcidset.compressed()

                print 'epc_min,epc_max,setlen=',epcid_min,epcid_max,len(lepcidset)
                for epcid in lepcidset:
                    if epcid ==-9999: continue
                    a1bool = ma.masked_equal(a1epcid, epcid).mask
                    DatTmp  = Dat[a1bool]
                    dstack[epcid].extend(DatTmp.tolist())
                    #-- save large arrays ---
                    if (len(dstack[epcid]))>maxrec:
                        dnum[epcid] = dnum[epcid] +1 
        
                        outDir  = outbaseDir + '/%s/%04d%02d'%(varName,Year,Mon)
                        outPath = outDir + '/%s.%05d.%d.npy'%(varName,epcid, dnum[epcid])
                        util.mk_dir(outDir)
                        np.save(outPath, array(dstack[epcid]).astype(dattype[var]))
        
                        del dstack[epcid]
                        dstack[epcid] = deque([])
                        gc.collect()
                         
                        print outPath
        
            #--- save ---
            for epcid in range(epcid_min,epcid_max):
                dnum[epcid] = dnum[epcid] + 1 
                if len(dstack[epcid])>0:
                    outDir  = outbaseDir + '/%s/%04d%02d'%(varName,Year,Mon)
                    outPath = outDir + '/%s.%05d.%d.npy'%(varName,epcid,dnum[epcid])
                    util.mk_dir(outDir)
    
                    np.save(outPath, array(dstack[epcid]).astype(dattype[var]))
    
                    del dstack[epcid]
                    gc.collect()
    
    
            #--- Joint segments ---
            for epcid in range(epcid_min,epcid_max):
                outDir  = outbaseDir + '/%s/%04d%02d'%(varName,Year,Mon)
                ssearch = outDir + '/%s.%05d.*.npy'%(varName,epcid)
                lsrcPath= sort(glob.glob(ssearch))
    
                if len(lsrcPath)==0: continue
    
                aout = deque([])
                for srcPath in lsrcPath:
                    atmp = np.load(srcPath)
                    aout.extend(atmp)
    
                outPath =  outDir + '/%s.%05d.npy'%(varName,epcid)
    
                np.save(outPath, aout)
                print outPath
    
                #-- delete temporary files --
                for srcPath in lsrcPath:
                    os.remove(srcPath)
    


