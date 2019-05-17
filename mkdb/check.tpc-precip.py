import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from numpy import *
import numpy as np
import myfunc.util as util
import sys, os

#calcFlag= True
calcFlag= False
figFlag = True
iYM = [2017,1]
eYM = [2017,12]
lYM = util.ret_lYM(iYM,eYM)
#NPC_USE = 3
#NBINS   = 22
#pctype  = 'tpc'

NPC_USE = 3
NBINS   = 25
pctype  = 'epc'

lbnd = arange(0,NBINS**NPC_USE+1) -0.5
#-- Read orbit list  ---
listDir = '/work/hk01/utsumi/PMM/TPCDB/list'
lorbit  = []
for (Year,Mon) in lYM:
    listPath = listDir + '/list.1C.V05.%04d%02d.csv'%(Year,Mon)
    f=open(listPath,'r'); lines = f.readlines(); f.close()
    for line in lines:
        line = map(int, line.split(','))
        #print line
        lorbit.append(line)

#-----------------------
a1countAll = zeros(len(lbnd)-1, int32)
a1countPr  = zeros(len(lbnd)-1, int32)
a1sumAll   = zeros(len(lbnd)-1, float32)
a1sumPr    = zeros(len(lbnd)-1, float32)

for orbinfo in lorbit:
    if calcFlag !=True: continue
    oid,Year,Mon,Day,itime,etime = orbinfo
    print oid,Year,Mon,Day
    matchBaseDir = '/work/hk01/utsumi/PMM/MATCH.GMI.V05A'
    dprDir = matchBaseDir + '/S1.ABp083-137.Ku.V06A.9ave.precipRateESurface/%04d/%02d/%02d'%(Year,Mon,Day)
    dprPath = dprDir + '/precipRateESurface.%06d.npy'%(oid)
    a2dpr   = np.load(dprPath)[:,103-83:103-83+15]

    if pctype == 'epc':
        pcidDir = matchBaseDir + '/S1.ABp103-117.GMI.epcid-s1/%04d/%02d/%02d'%(Year,Mon,Day)
        pcidPath= pcidDir + '/epcid-s1.%06d.npy'%(oid)
    elif pctype=='tpc':
        pcidDir = matchBaseDir + '/S1.ABp103-117.GMI.tpcid-s1-%dpc-%dbin/%04d/%02d/%02d'%(NPC_USE,NBINS,Year,Mon,Day)
        pcidPath= pcidDir + '/tpcid-s1.%06d.npy'%(oid)


    a2pcid = np.load(pcidPath)

    a1dprAll = ma.masked_less(a2dpr,0).compressed()
    a1pcAll  = ma.masked_where(a2dpr<0, a2pcid).compressed()

    a1dprPr  = ma.masked_less_equal(a2dpr,0).compressed()
    a1pcPr   = ma.masked_where(a2dpr<=0, a2pcid).compressed()

    a1countAllTmp = np.histogram(a1pcAll, bins=lbnd)[0]
    a1sumAllTmp   = np.histogram(a1pcAll, bins=lbnd, weights=a1dprAll)[0]

    a1countPrTmp  = np.histogram(a1pcPr, bins=lbnd)[0]
    a1sumPrTmp    = np.histogram(a1pcPr, bins=lbnd, weights=a1dprPr)[0]

    a1countAll= a1countAll + a1countAllTmp
    a1sumAll  = a1sumAll   + a1sumAllTmp
    a1countPr = a1countPr + a1countPrTmp
    a1sumPr   = a1sumPr   + a1sumPrTmp


#-- Save ---
outDir = '/home/utsumi/temp/TPC'
countAllPath = outDir + '/countAll.%s.%dpc-%dbin.npy'%(pctype,NPC_USE,NBINS)
countPrPath  = outDir + '/countPr.%s.%dpc-%dbin.npy'%(pctype,NPC_USE,NBINS)
sumAllPath   = outDir + '/sumAll.%s.%dpc-%dbin.npy'%(pctype,NPC_USE,NBINS)
sumPrPath    = outDir + '/sumPr.%s.%dpc-%dbin.npy'%(pctype,NPC_USE,NBINS)

if calcFlag==True:
    np.save(countAllPath, a1countAll)
    np.save(countPrPath, a1countPr)
    np.save(sumAllPath, a1sumAll)
    np.save(sumPrPath, a1sumPr)
    print countAllPath
   
#*************************
# Figure 
#*************************
a1countAll = np.load(countAllPath)
a1countPr  = np.load(countPrPath)
a1sumAll   = np.load(sumPrPath)
a1sumPr    = np.load(sumPrPath)
a1precAll  = ma.masked_invalid(a1sumAll / a1countAll)
a1precPr   = ma.masked_invalid(a1sumPr  / a1countPr )


a1id = arange(0,NBINS**NPC_USE).astype(int32)
if NPC_USE==3:
    a1id1= (a1id/(NBINS**2)).astype(int32)%NBINS
    a1id2= (a1id/(NBINS   )).astype(int32)%NBINS
    a1id3= a1id %NBINS

elif NPC_USE==4:
    a1id1= (a1id/(NBINS**3)).astype(int32)%NBINS
    a1id2= (a1id/(NBINS**2)).astype(int32)%NBINS
    a1id3= (a1id/(NBINS   )).astype(int32)%NBINS
    a1id4= a1id %NBINS

else:
    print 'check NPC_USE',NPC_USE
    sys.exit()

#-- Make mask ---
dmask1 = {}
dmask2 = {}
dmask3 = {}
dmask4 = {}
for i in range(NBINS):
    dmask1[i] = ma.masked_equal(a1id1,i).mask
    dmask2[i] = ma.masked_equal(a1id2,i).mask
    dmask3[i] = ma.masked_equal(a1id3,i).mask
if NPC_USE==4:
    dmask4[i] = ma.masked_equal(a1id4,i).mask

#****************************
# make 2D map ( PC1 vs PC2) --
a2countAll12 = zeros([NBINS,NBINS],int32)
a2countPr12  = zeros([NBINS,NBINS],int32)
a2countRat12 = zeros([NBINS,NBINS],float32)
a2precAll12  = zeros([NBINS,NBINS],float32)
a2precPr12   = zeros([NBINS,NBINS],float32)

for i in range(NBINS):
    for j in range(NBINS):
        a1mask = dmask1[i] * dmask2[j]
        countAll= ma.masked_less(a1countAll[a1mask],0).sum()
        countPr = ma.masked_less(a1countPr [a1mask],0).sum()
        precAll = ma.masked_less(a1precAll [a1mask],0).mean()
        precPr  = ma.masked_less(a1precPr  [a1mask],0).mean()

        if countAll==0:
            countRat = 0
        else:

            countRat= float(countPr) / countAll 
 
        a2countAll12[j,i] = countAll
        a2countPr12[j,i]  = countPr
        a2countRat12[j,i] = countRat
        a2precAll12 [j,i] = precAll
        a2precPr12  [j,i] = precPr

#-- Draw PC1 vs PC2 -------
x = 1
y = 2
for varName in ['countAll','countPr','countRat','precAll','precPr']:
    vmin, vmax = None, None
    if varName == 'countAll':
        a2fig = a2countAll12
    elif varName == 'countPr':
        a2fig = a2countPr12
    elif varName == 'countRat':
        a2fig = a2countRat12
    elif varName == 'precAll':
        vmin, vmax = 0,2
        a2fig = a2precAll12
    elif varName == 'precPr':
        vmin, vmax = 0,5
        a2fig = a2precPr12

    a2fig = ma.masked_less_equal(a2fig,0)
    a2fig = ma.masked_invalid(a2fig)
    print a2fig.min(), a2fig.max()
    fig = plt.figure(figsize=(4,4))
    ax  = fig.add_axes([0.1,0.1,0.8,0.8])
    im  = ax.imshow(a2fig, origin='lower',interpolation='nearest',vmin=vmin,vmax=vmax)
    plt.colorbar(im, orientation='horizontal')
    plt.title('%s x:%d %d'%(varName,x,y))
    plt.xlabel('PC%d'%(x))
    plt.ylabel('PC%d'%(y))

    figPath = outDir + '/twod.%s.%dpc-%dbin.pc%d%d.%s.png'%(pctype, NPC_USE,NBINS,x,y,varName)
    plt.savefig(figPath)
    plt.clf()
    print figPath

#****************************
#-- Draw PC2 vs PC3 -------
if NBINS ==25:
    lpcrange1 = [[0,4],[5,9],[10,14],[15,19],[20,24]]
elif NBINS ==22:
    lpcrange1 = [[0,0],[1,5],[6,10],[11,15],[16,20],[21,21]]
elif NBINS ==10:
    lpcrange1 = [[0,1],[2,3],[4,5],[6,7],[8,9]]


dcountAll = {}
dcountPr  = {}
dcountRat = {}
dprecAll  = {}
dprecPr   = {}
for ipcrange1,pcrange1 in enumerate(lpcrange1):
    ipc1,epc1 = pcrange1

    #--------------------------
    a2countAll = zeros([NBINS,NBINS],int32)
    a2countPr  = zeros([NBINS,NBINS],int32)
    a2countRat = zeros([NBINS,NBINS],float32)
    a2precAll  = zeros([NBINS,NBINS],float32)
    a2precPr   = zeros([NBINS,NBINS],float32)

    a1mask1 = array([False]*NBINS**NPC_USE)
    for pc1 in range(ipc1,epc1+1): 
        a1mask1 = a1mask1 + dmask1[pc1]

    for i in range(NBINS):
        for j in range(NBINS):
            #a1mask = a1mask1 * dmask2[i] * dmask3[j]   # X=PC2, Y=PC3
            a1mask = a1mask1 * dmask3[i] * dmask2[j]   # X=PC3, Y=PC2
            countAll= ma.masked_less(a1countAll[a1mask],0).sum()
            countPr = ma.masked_less(a1countPr [a1mask],0).sum()
            precAll = ma.masked_less(a1precAll [a1mask],0).mean()
            precPr  = ma.masked_less(a1precPr  [a1mask],0).mean()
    
            if countAll==0:
                countRat = 0
            else:
                countRat= float(countPr) / countAll 
     
            a2countAll[j,i] = countAll
            a2countPr[j,i]  = countPr
            a2countRat[j,i] = countRat
            a2precAll [j,i] = precAll
            a2precPr  [j,i] = precPr

    dcountAll[ipcrange1] = a2countAll
    dcountPr [ipcrange1] = a2countPr
    dcountRat[ipcrange1] = a2countRat
    dprecAll [ipcrange1] = a2precAll
    dprecPr  [ipcrange1] = a2precPr


x,y = 3,2
for varName in ['countAll','countPr','countRat','precAll','precPr']:
#for varName in ['countPr']:

    fig = plt.figure(figsize=(20,4))
    for ipcrange1,pcrange1 in enumerate(lpcrange1):
        vmin, vmax = None, None
        if varName == 'countAll':
            a2fig = dcountAll[ipcrange1]
        elif varName == 'countPr':
            a2fig = dcountPr[ipcrange1]
        elif varName == 'countRat':
            a2fig = dcountRat[ipcrange1]
        elif varName == 'precAll':
            a2fig = dprecAll[ipcrange1]
            vmin, vmax = 0,2
        elif varName == 'precPr':
            a2fig = dprecPr[ipcrange1]
            vmin, vmax = 0,5

        a2fig = ma.masked_less_equal(a2fig,0)

        w  = 0.8/len(lpcrange1)
        ax = fig.add_axes([0.1+w*ipcrange1,0.1,w, 0.8])
        im = ax.imshow(a2fig, origin='lower', interpolation='nearest',vmin=vmin,vmax=vmax)
        plt.colorbar(im,orientation='horizontal') 
        plt.xlabel('PC%d'%(x))
        if ipcrange1 ==0:
            plt.ylabel('PC%d'%(y))

    #-- Save figure ----
    figPath = outDir + '/twod.%s.%dpc-%dbin.pc%d%d.%s.png'%(pctype,NPC_USE,NBINS,x,y,varName)
    plt.savefig(figPath) 
    print figPath
    plt.clf()
