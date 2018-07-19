from numpy import *
import myfunc.util as util
import GPMGV
import numpy as np

calc = True
#calc = False
iYM = [2005,4]
#iYM = [2014,4]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
print lYM
figtype = 'bar'
#figtype = 'cont'
#corrFlag= 'CR'
corrFlag= 'NO'

#cls = 'RH'
cls = 'RainType'

#thdist = 2.5
thdist = 5.0
minNum = 3
prdName = 'L2A25'
gv = GPMGV.GPMGV()
gv.load_sitelist_reclassified()

dgName = gv.ret_ddomYM2gName()

ldomain = gv.domains
#ldomain = ['FLORIDA-STJ','FLORIDA-SFL-N','N.Carolina-IPHEx_Duke','N.Carolina-IPHEx_NASA','KWAJALEIN-KWA']
#ldomain = ['FLORIDA-STJ']

offset_bef = 15  # 'bef' should be identical to that in mk_match.py
#offset_aft = 30

ldt  = [1, 5,10,15,20,25,30]
ldh  = [0] + range(2,34+1,2)

nt = len(ldt)
nh = len(ldh)
#nh = 5

lcls = ['RH','RainType']
dlclstype = {}
dlclstype['RH']       = ['all','dry','hum']
dlclstype['RainType'] = ['strat','conv']

ldattype = ['cc','brat','rmse']
#-------------------------------
def load_data(srcPath, nlev=1,ntime=1):
    f=open(srcPath,'r'); lines= f.readlines(); f.close()
    eSurf= float( lines[1].strip().split(',')[ntime+1])
    satelev= float(lines[nlev+1].strip().split(',')[ntime+1])
    return eSurf, satelev
     
#-------------------------------
def load_bestlev(srcPath,dattype='cc', ntime=1):
    f=open(srcPath,'r'); lines= f.readlines(); f.close()

    esurf = float(lines[1].strip().split(',')[ntime+1])
    lsatelev = []

    for line in lines[1:]:
        satelev = float(line.strip().split(',')[ntime+1])
        lsatelev.append(satelev)


    if dattype in ['cc']:
        ibest = np.argmax(lsatelev)
    else:
        ibest = np.argmin( map(abs, lsatelev) )

    satelev = lsatelev[ibest]
    rdif    = (abs(satelev)-abs(esurf))/abs(esurf)

    return ibest, rdif
 

#-------------------------------
lntime = [3,6]
lnlev  = [6,7,8,9,10,11,12,13,14]

sout = ''
for ntime in lntime:
    time= ldt[ntime]

    drdif  = {}
    dibest = {}
    for cls in lcls:
        lclstype = dlclstype[cls]
        for clstype in lclstype:
            for dattype in ldattype:
                figDir  = '/work/a01/utsumi/GPMGV/fig'
                csvPath = figDir + '/table.nt-nlev.%s.%.1fkm.minNum.%d.%s.%s.%s.%s.csv'%(prdName, thdist,minNum, dattype, corrFlag, cls, clstype)
    
                ibest, rdif = load_bestlev(csvPath, dattype=dattype, ntime=ntime)
   
                print cls, clstype, dattype, ibest 
                drdif [dattype,clstype] = rdif
                dibest[dattype,clstype] = ibest
 
    #-- sout ---
    sout  = sout + '\n' 
    sout  = sout + 'time=%d'%(time) + '\n'
    llabel= ['']+['all','hum','dry','conv','strat']
    slabel= ','.join(llabel)
    sout  = sout + slabel + '\n'
    for dattype in ldattype:
        line = []

        for clstype in ['all','hum','dry','conv','strat']:
            ibest = dibest[dattype,clstype]
            lev = 0.25*ldh[ibest]
            line.append(lev)
    
        sline = dattype + ',' +  ','.join(map(str, line))
        sout  = sout + sline + '\n'
 

        line = []
        for clstype in ['all','hum','dry','conv','strat']:
            rdif = drdif[dattype,clstype]*100
            line.append(rdif)
    
        sline = '%' + ',' +  ','.join(map(str, line))
        sout  = sout + sline + '\n'

   
outPath = figDir + '/summary.%s.bestlev.csv'%(corrFlag)
f=open(outPath,'w'); f.write(sout); f.close()
print outPath
    
    

