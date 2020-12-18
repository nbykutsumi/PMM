# %%
import numpy as np
from collections import OrderedDict as odict
import struct
import os, sys
import gzip
from datetime import datetime, timedelta

class level2(object):
    """
    Read GSMaP level 2 data
    """

    def __init__(self):
        self.npix = 90  # fixed. MHS


    def get_var(self, srcPath=None, lvname=[], nrec=None, origin=0, compressed=True):
        npix    = self.npix

        lout = []
        if compressed is True:
            f= gzip.open(srcPath, 'rb')
        else:
            f= open(srcPath, 'rb')

        shead  = f.read(8)   # Only 8?
        ncount = struct.unpack_from('<i',shead, offset=4)[0]  # number of scans


        dictvars    = odict((
                            ('year',                '%di'%(ncount)),       #
                            ('day',                 '%di'%(ncount)),       #
                            ('utc',                 '%di'%(ncount)),       #
                            ('startdate',           '%di'%(3)),       #
                            ('enddate',             '%di'%(3)),       #
                            ('lat',                 '%df'%(ncount*npix)),       #
                            ('lon',                 '%df'%(ncount*npix)),       #
                            ('sfc',                 '%df'%(ncount*npix)),       #
                            ('si',                  '%df'%(ncount*npix)),       #
                            ('lz',                  '%df'%(ncount*npix)),       #
                            ('qc',                  '%di'%(ncount*npix)),       #
                            ('snowm',               '%di'%(ncount*npix)),       #
                            ('rr',                  '%df'%(ncount*npix)),       #
        ))
        
        vars, fmts  = list(zip(*list(dictvars.items())))
        self.dictvars = dictvars.copy()
        #self.fmtsize  = struct.calcsize( '<'+ ''.join(self.fmts))   # bytes

        if 'dtime' in lvname:
            lvname_tmp = lvname.copy()
            lvname_tmp.remove('dtime')
        else:
            lvname_tmp = lvname

        for vname in lvname_tmp:
            vidx = vars.index(vname)
            fmtsize_pre= struct.calcsize( '<'+ ''.join(fmts[:vidx]))

            f.seek(8+fmtsize_pre, 0)

            varsize = struct.calcsize( '<'+ ''.join(fmts[vidx]))
            sdat = f.read(varsize)

            if vidx <= 2:  # year, day, utc
                sfmt_all = '4s'*ncount
                adat = np.array(struct.unpack(sfmt_all, sdat))
            elif (vidx >2)&(vidx<=4):  # startdate, enddate
                sfmt_all = '4s'*3
                adat = np.array(struct.unpack(sfmt_all, sdat))
            else:
                sfmt_all = '4s'*ncount*npix
                adat = np.array(struct.unpack(sfmt_all, sdat)).reshape(ncount,-1)
            adat.dtype = fmts[vidx][-1]

            if (vidx >2)&(vidx<=4):  # startdate, enddate
                lout.append(adat)
            else: 
                if nrec is None:
                    lout.append(adat[origin:])
                else: 
                    lout.append(adat[origin:origin+nrec])
        f.close()

#        if 'dtime' in lvname:
#            dtime 
    

        print(lvname)
        print(lout)

        #lout = []
        #for vname in lvname:
        #    vidx = vars.index(vname)
        #    atmp = sdat[:,vidx]
        #    atmp.dtype = fmts[vidx]
        #    if vname =='dtime':
        #        adtime = np.array([datetime(1801,1,1,0) + timedelta(minutes=i) for i in atmp.tolist()])
        #        lout.append(adtime.reshape(-1,npix)[:,0])  # 1-dimensional array
        #    else:
        #        lout.append(atmp.reshape(-1,npix))

        return lout
# %%
