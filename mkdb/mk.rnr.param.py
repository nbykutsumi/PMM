#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
from numpy import *
from scipy.stats import gaussian_kde
import socket
import random
import os, sys
import numpy as np
import epcfunc
import EPCDB
import myfunc.util as util

DB_MAXREC  = 10000
DB_MINREC  = 1000

sensor = 'GMI'
myhost = socket.gethostname()
if myhost =="shui":
    dbDir   = '/tank/utsumi/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
    outDir  = '/home/utsumi/temp/ret'
elif myhost =="well":
    coefDir = '/media/disk2/share/PMM/EPCDB/EPC_COEF/%s'%(sensor)
    dbDir   = '/media/disk2/share/PMM/EPCDB/samp.%d.GMI.V05A.S1.ABp103-117.01-12'%(DB_MAXREC)
    outDir  = '/home/utsumi/temp/ret'


def read_table(srcPath, type=float):
    f=open(srcPath,'r'); lines=f.readlines(); f.close()
    lout = []
    for line in lines:
        line = line.strip().split()
        line = map(type, line)
        lout.append(line)
    return array(lout)

def read_nrain(idx_db):
    '''
     /* first six= Ntotal then Nrain exceeding 1 5 10 20 50 mm/hr when T2m < 278K */
     /* second six= same for when T2m > 278K */
    '''
    #srcPath = dbDir + '/db_%05d.bin.nrain.txt'%(idx_db)
    nrainDir= dbDir + '/nrain'
    srcPath = nrainDir + '/db_%05d.bin.nrain.txt'%(idx_db)
    a1nrain = read_table(srcPath,type=int32)[0]
    a1nrain_cold = a1nrain[:6]
    a1nrain_warm = a1nrain[6:]
    return a1nrain_warm, a1nrain_cold


#* Chose database *****
dbtype= 'my'
if dbtype =='JPL':
    db    = JPLDB.JPLDB()
elif dbtype=='my':
    db    = EPCDB.EPCDB()
else:
    print 'check dbtype',dbtype
#**********************
NEM = 12

#lthnall = [[1,10000],[10001,20000]][::-1]
lthnall = [[1,10000]]
lthpr = [0, 0.1, 1, 5, 10]
thrma0 = 0.05
thrma5 = 0
NEM_USE    = 3
NPCHIST    = 29
NDB_EXPAND = 20

for [nallmin, nallmax] in lthnall:
    ##** Make idx_db list ***
    ndb = 29*29*29  # 24388
    lidx_db = []
    for idx_db in range(ndb):
        nrainDir= dbDir + '/nrain'
        srcPath = nrainDir + '/db_%05d.bin.nrain.txt'%(idx_db)
        if not os.path.exists(srcPath):
            continue
        _,a1nrain = read_nrain(idx_db)
        nall = a1nrain[0]
        if (nallmin <= nall) and (nall <=nallmax):
            lidx_db.append(idx_db)

    #lidx_db = np.sort(lidx_db)
    #***********************
    #lidx_db = np.sort(random.sample(lidx_db, 10))
    #lidx_db = [3312,3948,4750,8859,9503,9865,13668,13752,14425,19618,20666,22617]
    #lidx_db = [22617]
    #lidx_db = [20666]
    #lidx_db = [3312]
    #lidx_db = [4750]
    #lidx_db = [9865,4750]
    #lidx_db = [9865]
    #lidx_db = [4120]

    icount = 0
    #******************************
    # Start pribary idx loop
    #------------------------------
    for idx_db in lidx_db:
        #******************************
        # Make expanded idx_db list
        #------------------------------
        lidx_db_expand_tmp = [idx_db] + [idx_db + i*sign for i in range(1,NDB_EXPAND) for sign in [-1,1]]

        nevent_all = 0

        lidx_db_expand = []
        for idx_db_expand in lidx_db_expand_tmp:

            #-- If idx_db == -9999 --
            if ((idx_db_expand <0)or(pow(NPCHIST, NEM_USE)-1<idx_db_expand)):
                print 'No matching database'
                print 'idx_db=',idx_db_expand
                continue
    
            #-- Read nrain file --
            try:
                a1nrain_warm, a1nrain_cold = read_nrain(idx_db_expand)
            except IOError:
                print 'No DB file for idx_db=',idx_db_expand
                print 'SKIP'
                continue
    
            nevent_all  = nevent_all  + a1nrain_cold[0] + a1nrain_warm[0]
            lidx_db_expand.append(idx_db_expand)
            if nevent_all >DB_MINREC: break
        #******************************
        # Stack data
        #------------------------------
        a2epcdb = []
        a1prec  = []
        a1t2m   = []
        for idx_db_expand in lidx_db_expand:
            dbPath = dbDir + '/db_%05d.bin'%(idx_db_expand)
            if   dbtype == 'JPL':
                db.set_file(dbPath)
            elif dbtype == 'my':
                db.set_idx_db(dbDir, idx_db_expand)
        
            a2epcdb.append(db.get_var('pc_emis')[:,:NEM])  # (nrec, 12)
            a1prec .append(db.get_var('precip_NS_cmb'))
            a1t2m  .append(db.get_var('t2m'))
 
        a2epcdb = np.concatenate(a2epcdb, axis=0)
        a1prec  = np.concatenate(a1prec, axis=0)
        a1t2m   = np.concatenate(a1t2m,  axis=0)
        #------------------------------

        print 'idx_db=',idx_db
        irregular = 0
        #-- Screening ----------
        a1flagepc= ma.masked_not_equal(a2epcdb,-9999.).all(axis=1).mask
        a1flagT  = ma.masked_greater(a1t2m,0).mask   # screen T2m=0
        a1flagP1 = ma.masked_greater_equal(a1prec,0).mask
        a1flagP2 = ~ma.masked_invalid(a1prec).mask
        a1flag_good = a1flagepc * a1flagT * a1flagP1 * a1flagP2
   
 
        #-- Only valid data ---
        a2epcdb      = a2epcdb[a1flag_good,:]
        a1prec = a1prec[a1flag_good]
        a1t2m        = a1t2m[a1flag_good]
        #----------------------
        a1flag_rain = ma.masked_greater(a1prec, 0).mask
        a1flag_no   = ma.masked_equal(a1prec,0).mask
    
        #-------------
        if a1flag_rain.sum() + a1flag_no.sum() ==0:
            print 'No valid data in idx_db',idx_db

        if a1flag_good.sum() <=1: irregular = 1
        if type(a1flag_no) is np.bool_: irregular=1
        if type(a1flag_rain) is np.bool_: irregular=2
    
   
        if irregular ==0: 
            icount = icount + 1
    
            a1d_all = None
            a1d_rain= None
            a1d_no  = None

            lm_org   = []
            lm_rain  = []
            lm_no    = []
            lstd_org = []
            lstd_rain= []
            lstd_no  = []

            for iepc in range(NEM+1):
                if iepc==NEM:
                    a1epc_all  = a1t2m
                    a1epc_rain = a1t2m[a1flag_rain]
                    a1epc_no   = a1t2m[a1flag_no  ]
                else:
                    a1epc_all  = a2epcdb[:, iepc]
                    a1epc_rain = a2epcdb[a1flag_rain,iepc]
                    a1epc_no   = a2epcdb[a1flag_no  ,iepc]
    
                #-- Normalize ----
                m_org   = a1epc_all.mean()
                std_org = a1epc_all.std()
                a1epc_all = (a1epc_all - m_org)/std_org
                a1epc_rain= (a1epc_rain- m_org)/std_org
                a1epc_no  = (a1epc_no  - m_org)/std_org
    
                m_rain  = a1epc_rain.mean()
                m_no    = a1epc_no.mean()
                std_rain= a1epc_rain.std()
                std_no  = a1epc_no.std()
    
                a1d_allTmp  = (m_no - m_rain)/(std_no + std_rain)*(m_no - a1epc_all)
                a1d_rainTmp = (m_no - m_rain)/(std_no + std_rain)*(m_no - a1epc_rain)
                a1d_noTmp   = (m_no - m_rain)/(std_no + std_rain)*(m_no - a1epc_no)
    
                if a1d_all is None:
                    a1d_all  = a1d_allTmp
                    a1d_rain = a1d_rainTmp
                    a1d_no   = a1d_noTmp
                else: 
                    a1d_all  = a1d_all  + a1d_allTmp
                    a1d_rain = a1d_rain + a1d_rainTmp
                    a1d_no   = a1d_no   + a1d_noTmp


                lm_org  .append(m_org)
                lm_rain .append(m_rain)
                lm_no   .append(m_no) 
                lstd_org  .append(std_org)
                lstd_rain .append(std_rain)
                lstd_no   .append(std_no) 
 
            #** Descrimination *****
            vmin = np.percentile(a1d_all,0.05)
            m_rain = a1d_rain.mean()
            lthd = np.linspace(vmin, m_rain, 30)[::-1] # larger to smaller
            lrma0= []
            lrma5= []
            lsr  = []
            lwer = []
            for thd in lthd:
                a1flag_est_no   = ma.masked_less(a1d_all, thd).mask
                a1flag_est_rain = ~a1flag_est_no
                a1flag_rain_5mm = ma.masked_greater_equal(a1prec, 5).mask
                sr  = a1flag_est_no.sum() / float(a1flag_est_no.sum() + a1flag_est_rain.sum())  # screening ratio
                wer = a1flag_rain.sum() / float(len(a1flag_rain))  # wet event ratio
    
                rma0= a1prec[a1flag_est_no].sum() / a1prec.sum()  # Ratio of missing amount
                rma5= a1prec[a1flag_est_no * a1flag_rain_5mm].sum() / a1prec[a1flag_rain_5mm].sum()

                lsr.append(sr)
                lwer.append(wer)
                lrma0.append(rma0)
                lrma5.append(rma5)
    
            #-- Find best threshold --
            #print ''
            for i in range(len(lthd)):  # loop from large thd to small thd
                thd = lthd[i]
                sr  = lsr[i]
                wer = lwer[i]
                rma0= lrma0[i]
                rma5= lrma5[i]
                #print idx_db, thd, sr, rma
                if (rma0<=thrma0)and(rma5<=thrma5):
                    break
            #-------------
            a1flag_est_no   = ma.masked_less(a1d_all, thd).mask
            a1flag_est_rain = ~a1flag_est_no
    
    
            a1flag_est_no   = ma.masked_less(a1d_all, thd).mask
            a1flag_est_rain = ~a1flag_est_no
    
            dpod = {}
            dfar = {}
            dmr  = {}
            drma = {}
            sline = '%d,%f,%f,%f'%(idx_db,thd,sr,wer)
    
            #print a1prec.shape, a1flag_rain.shape, a1flag_est_rain.shape
            for thpr in lthpr:
                a1flagP  = ma.masked_greater_equal(a1prec, thpr).mask
                
                OE = (a1flagP * a1flag_rain * a1flag_est_rain).sum()
                Oe = (a1flagP * a1flag_rain * a1flag_est_no).sum()
                oE = (a1flagP * a1flag_no   * a1flag_est_rain).sum()
                oe = (a1flagP * a1flag_no   * a1flag_est_no).sum()
                pod = OE/float(OE+Oe)
                far = oE/float(OE+oE)
                mr  = Oe/float(OE+Oe)
    
                rma = a1prec[a1flag_est_no*a1flagP].sum() / a1prec[a1flagP].sum()  # Ratio of missing amount
                nall= OE+Oe+oe+oE
                #print '----------------------'
                #print 'thrp=',thpr
                #print 'nall=',nall
                #print 'OE  =',OE
                #print 'Oe  =',Oe
                #print 'oe  =',oe
                #print 'oE  =',oE
                #print ''
                #print 'POD =',pod
                #print 'FAR =',far 
                #print 'MR  =',mr
                #print 'RMA =',rma
                #print 'random, OE/(OE+oE)', p_rain, OE/float(OE+oE)
                #print ''
                dpod[thpr] = pod
                dfar[thpr] = far
                dmr [thpr] = mr
                drma[thpr] = rma
    
                lout  = [nall,OE,Oe,oe,oE]
                sline = sline + ',' + ','.join(map(str,lout))
            
            lout = [d[thpr] for d in [dpod,dfar,dmr,drma] for thpr in lthpr]
            sline = sline + ',' + ','.join(map(str,lout))
             
    
            if (icount==1):
                llabel = ['db_idx','D','SR','WER'] + ['%s %.1f'%(var,thpr) for thpr in lthpr for var in ['nall','OE','Oe','oe','oE']]
                llabel = llabel + ['%s %.1f'%(var,thpr) for var in ['POD','FAR','MR','RMA'] for thpr in lthpr]
                slabel = ','.join(map(str,llabel))
                sout   = slabel + '\n'
                sout   = sout + sline + '\n'
            else:
                sout = sout + sline + '\n'
         
        #--- Save RNR paremeter files --
        
        print idx_db,irregular
        if irregular==1:
            calcflag = 2  # calc all pixels
        elif irregular==2:
            calcflag = 0  # skip all pixels
        elif np.isinf(thd):
            calcflag = 2  # calc all pixels
        else:
            calcflag = 1 


        if calcflag in [0,2]:
            thd  = -9999.
            sr   = -9999.
            wer  = -9999.
            rma  = -9999.
            drma = {}
            for thpr in lthpr: drma[thpr] = -9999.
            lm_org  = [-9999.]*(NEM+1)
            lm_rain = [-9999.]*(NEM+1)
            lm_no   = [-9999.]*(NEM+1)
            lstd_org  = [-9999.]*(NEM+1)
            lstd_rain = [-9999.]*(NEM+1)
            lstd_no   = [-9999.]*(NEM+1)


        for thpr in lthpr:
            if np.isnan(drma[thpr]):
                drma[thpr] = 0.0 


        sout2 = 'calc\t%d\n'%(calcflag)  # 0: Skip all  1:Screening  2:Calc all
        sout2 = sout2 + 'D\t%.5f\n'%(thd)
        sout2 = sout2 + 'SR\t%.5f\n'%(sr)
        sout2 = sout2 + 'WER\t%.5f\n'%(wer)
        sout2 = sout2 + 'RMA\t' + '\t'.join(map(str, ['%.5f'%(drma[thpr]) for thpr in lthpr])) + '\n'
        sout2 = sout2 + 'm_org\t'  + '\t'.join(map(str, lm_org)) + '\n'
        sout2 = sout2 + 'm_rain\t' + '\t'.join(map(str, lm_rain)) + '\n'
        sout2 = sout2 + 'm_no\t'   + '\t'.join(map(str, lm_no)) + '\n'
        sout2 = sout2 + 's_org\t'  + '\t'.join(map(str, lstd_org)) + '\n'
        sout2 = sout2 + 's_rain\t' + '\t'.join(map(str, lstd_rain)) + '\n'
        sout2 = sout2 + 's_no\t'   + '\t'.join(map(str, lstd_no)) + '\n'

        
        #rnrDir = dbDir + '/rnr'
        rnrDir = dbDir + '/rnr.minrec%d'%(DB_MINREC)
        util.mk_dir(rnrDir)
        rnrPath = rnrDir + '/rnr.%05d.txt'%(idx_db)
        f=open(rnrPath,'w'); f.write(sout2); f.close()
        #-------------------------------    
    #-- Save csv table -------------
    outPath = outDir + '/RNR.linear.maxrec%d.minrec%d.%05d-%05d.csv'%(DB_MAXREC,DB_MINREC,nallmin,nallmax)
    f=open(outPath,'w'); f.write(sout); f.close()
    print outPath



