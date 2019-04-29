import numpy as np
import numpy.ma as ma
from bisect import bisect_left
import sys
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


#------------------------------------------
def extract_domain_2D(a2dat, a2lat, a2lon, clat, clon, dlatlon, dscan, returnidx=False):

    nyTmp, nxTmp = a2lat.shape
    a1lat = a2lat[:,nxTmp/2]
    a1lon = a2lon[:,nxTmp/2]

    idx_latmax = np.argmax(a1lat)
    a1lat0 = a1lat[:idx_latmax+1]
    a1lat1 = a1lat[idx_latmax+1:]
    a1lon0 = a1lon[:idx_latmax+1]
    a1lon1 = a1lon[idx_latmax+1:]

    if (-180<=clat)and(clat <=180):
        #-- search first half: ascending --
        found_domain = 0
        if len(a1lat0) !=1:
            idx_c  = bisect_left(a1lat0, clat)
            #print 'clat=',clat
            #print 'a1lat0=',a1lat0
            #print 'a1lat1=',a1lat1
            latTmp = a1lat0[idx_c]
            lonTmp = a1lon0[idx_c]
            if (clat-dlatlon<=latTmp)&(latTmp <=clat+dlatlon)&(clon-dlatlon<=lonTmp)&(lonTmp<=clon+dlatlon):
                found_domain = 1

        if found_domain !=1:
            #-- search second half: descending --
            idx_c  = bisect_left(a1lat1[::-1], clat)
            idx_c  = len(a1lat) - idx_c -1
            latTmp = a1lat[idx_c]
            lonTmp = a1lon[idx_c]
            #print 'A',idx_c, len(a1lat1)
            if (clat-dlatlon<=latTmp)&(latTmp <=clat+dlatlon)&(clon-dlatlon<=lonTmp)&(lonTmp<=clon+dlatlon):
                found_domain =1

                #print 'B',idx_c, len(a1lat1)

        if found_domain==1:
            idx_first = idx_c - dscan
            idx_last  = idx_c + dscan
            if idx_first<0: idx_first=0
            if idx_first>=len(a1lat): idx_first=len(a1lat)-1
            if idx_last<0: idx_last=0
            if idx_last>=len(a1lat): idx_last=len(a1lat)-1


            a2odat  = a2dat[idx_first:idx_last+1,:]
            a2olat  = a2lat[idx_first:idx_last+1,:]
            a2olon  = a2lon[idx_first:idx_last+1,:]

        else:
            print 'No matching scans in the target domain are found.'
            print 'Exit'
            sys.exit()

        print 'Extract target domain'
        print 'Extracted array size=', a2odat.shape
        print 'idx_c='         ,idx_c
        print 'idx_first/last=',idx_first, idx_last
    else:
        pass
    if returnidx==False:
        return a2odat, a2olat, a2olon
    else:
        return a2odat, a2olat, a2olon, idx_first, idx_last


