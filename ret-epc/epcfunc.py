import numpy as np
import numpy.ma as ma
#-------------------
def mk_epc_12pc(a3tb, a2coef):
    NEM = 12
    NTBREG = a3tb.shape[2]
    print 'NTBREG=',NTBREG
    ny,nx,ntmp = a3tb.shape
    a3epc = np.zeros([ny,nx,NEM])

    # +a0
    kt  = 0
    for iem in range(NEM):
        a3epc[:,:,iem] = a2coef[kt,iem]

        ##-- test ----
        #print iem,a3epc[0,0,iem]
        ##------------

    for i in range(NTBREG):   # NTBRERG=number of TBs
        #print 'iTBREG=',i,'out of',NTBREG
        kt =kt+1

        # + bi*Tbi
        for iem in range(NEM): # NEM = number of output emissivity
            b = a2coef[kt, iem]
            a3epc[:,:,iem] = a3epc[:,:,iem] + b*a3tb[:,:,i]

            #-- test --
            #print 'i,iem',i,iem,b,a3tb[0,0,i]
            #print 'i,iem',i,iem,a3epc[0,0,iem]
            #print 'i,iem',i,iem,b
            #----------

        # + cij * Tbi * Tbj
        for j in range(i,NTBREG):
            kt = kt+1
            #print i,j
            a2tmp = a3tb[:,:,i] * a3tb[:,:,j]

            for iem in range(NEM):
                c = a2coef[kt, iem]
                a3epc[:,:,iem] = a3epc[:,:,iem] + c*a3tb[:,:,i]*a3tb[:,:,j]

                ##-- test --
                #if iem==emtemp:
                #    print 'i,j,c, Tbi*Tbj',i,j,c,a3tb[iy,ix,i]*a3tb[iy,ix,j]
                ##----------


    # + d * (Tbi - Tbi-1)(Tbi + Tbi-1)
    for i in range(1,NTBREG):
        kt = kt+1
        for iem in range(NEM):
            d = a2coef[kt,iem]
            a3epc[:,:,iem] = a3epc[:,:,iem] \
                 +d*(a3tb[:,:,i]-a3tb[:,:,i-1])/(a3tb[:,:,i]+a3tb[:,:,i-1])

    return a3epc

#------------------------------------------
def mk_epc_id_25bins(a3epc, a2pc_edge):
    NEM_USE   = 3
    ny,nx,nz  = a3epc.shape
    a2id_db  = np.zeros([ny,nx],np.int32)
    for iem in range(NEM_USE):
        a1bin = a2pc_edge[iem]
        a2idTmp = np.digitize(a3epc[:,:,iem], a1bin, right=False) - 1

        a2idTmp = ma.masked_outside(a2idTmp,0,24)
        a2id_db = a2id_db + a2idTmp*pow(25, NEM_USE-1-iem)

    a2id_db = a2id_db.filled(-9999)
    return a2id_db



