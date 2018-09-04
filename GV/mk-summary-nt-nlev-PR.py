from numpy import *
import myfunc.util as util
import GPMGV

calc = True
#calc = False
iYM = [2005,4]
#iYM = [2014,4]
eYM = [2014,10]
lYM = util.ret_lYM(iYM, eYM)
lYM = [YM for YM in lYM if YM[1] not in [1,2,3,11,12]]
print lYM
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
ldh  = arange(0, 8.5+0.01, 0.5)

nt = len(ldt)
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
    h    = lines[nlev+1].strip().split(',')[0]
    t    = int(lines[0].strip().split(',')[ntime+1])
    return t, eSurf, h, satelev
     

#-------------------------------
lntime = [3,6]
lnlev  = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]

sout = ''
for (ntime,nlev) in [[ntime,nlev] for ntime in lntime for nlev in lnlev]:
    time= ldt[ntime]
    lev = ldh[nlev]
    print time, lev

    drdif  = {}
    desurf = {}
    dprof  = {}
    for cls in lcls:
        lclstype = dlclstype[cls]
        for clstype in lclstype:
            for dattype in ldattype:
                figDir  = '/work/a01/utsumi/GPMGV/fig'
                csvPath = figDir + '/table.nt-nlev.%s.%.1fkm.minNum.%d.%s.%s.%s.%s.csv'%(prdName, thdist,minNum, dattype, corrFlag, cls, clstype)
    
                t, esurf, h, satelev= load_data(csvPath, nlev=nlev, ntime=ntime)
                print h, cls, clstype, dattype, esurf, satelev
    
                rdif = (abs(satelev) - abs(esurf))/abs(esurf) *100
    
                drdif[dattype,clstype] = rdif
                desurf[dattype,clstype]= esurf
                dprof [dattype,clstype]= satelev

    #-- sout ---
    sout  = sout + '\n' 
    sout  = sout + 'time=%d ,lev=%.1fkm'%(time, lev) + '\n'
    #llabel= ['']+['All(eSurf)', 'All(VeticalAverage)','Humid(eSurf)','Humid(VerticalAverage)','Dry(eSurf)','Dry(VerticalAverage)','Convective(eSurf)','Convective(VerticalAverage)','Stratiform(eSurf)','Stratiform(VerticalAverage)']

    llabel= [''] + ['All','Humid','Dry','Convective','Stratiform']
    slabel= ','.join(llabel)
    sout  = sout + slabel + '\n'
    for dattype in ldattype:
        line = []
        for clstype in ['all','hum','dry','conv','strat']:
                line.append(desurf[dattype,clstype])
        sline = dattype +'(eSurf)' + ',' +  ','.join(map(str, line))
        sout  = sout + sline + '\n'

        line = []
        for clstype in ['all','hum','dry','conv','strat']:
                line.append(dprof[dattype,clstype])
        sline = dattype +'(VerticalAve)' + ',' +  ','.join(map(str, line))
        sout  = sout + sline + '\n'

        line = []
        for clstype in ['all','hum','dry','conv','strat']:
                line.append(drdif[dattype,clstype])
        sline = dattype +'(Change)' + ',' +  ','.join(map(str, line))
        sout  = sout + sline + '\n'



    ##-- sout ---
    #sout  = sout + '\n' 
    #sout  = sout + 'time=%d ,lev=%.1fkm'%(time, lev) + '\n'
    #llabel= ['']+['all','hum','dry','conv','strat']
    #slabel= ','.join(llabel)
    #sout  = sout + slabel + '\n'
    #for dattype in ldattype:
    #    line = []
    #    for clstype in ['all','hum','dry','conv','strat']:
    #        line.append(drdif[dattype,clstype])
    #
    #    sline = dattype + ',' +  ','.join(map(str, line))
    #    sout  = sout + sline + '\n'
    
outPath = figDir + '/summary.%s.gain.csv'%(corrFlag)
f=open(outPath,'w'); f.write(sout); f.close()
print outPath
    
    

