# %%
import numpy as np
from collections import OrderedDict as odict
import struct
import os
import gzip
from datetime import datetime, timedelta

class level2(object):
    """
    Read GSMaP level 2 data
    """

    def __init__(self):
        self.npix = 180  # fixed. SSMIS

        dictvars    = odict((
                            ('pad1',                '4x'),
                            ('lon',                 'f'),       #
                            ('lat',                 'f'),       #
                            ('dtime',               'i'),       #
                            ('itoil',               'f'),       #
                            ('irflg',               'f'),       #
                            ('rainfg',              'f'),       #
                            ('snowprb',             'f'),       #
                            ('pad2',                '4x'),
        ))

        self.vars, self.fmts  = list(zip(*list(dictvars.items())))
        self.dictvars = dictvars.copy()
        self.fmtsize  = struct.calcsize( '<'+ ''.join(self.fmts))   # bytes

    def get_var(self, srcPath=None, lvname=[], nrec=None, origin=0, compressed=True):
        fmts    = self.fmts
        fmtsize = self.fmtsize
        vars    = self.vars
        npix    = self.npix


        if compressed is True:
            with gzip.open(srcPath, 'rb') as f:
                shead  = f.read(12)
                ncount = struct.unpack_from('<i',shead, offset=4)[0]  # number of scans

                if nrec is None:
                    nrec = ncount - origin
                f.seek(fmtsize*origin*npix+12, 0)
                sdat = f.read(fmtsize*nrec*npix)

        else:
            with open(srcPath, 'rb') as f:
                shead  = f.read(12)
                ncount = struct.unpack_from('<i',shead, offset=4)[0]  # number of scans

                if nrec is None:
                    nrec = ncount - origin
                f.seek(fmtsize*origin*npix+12, 0)
                sdat = f.read(fmtsize*nrec*npix)

        sfmt_all = '4s'*len(vars)*nrec*npix
        sdat = np.array(struct.unpack(sfmt_all, sdat)).reshape(nrec*npix,-1)

        lout = []
        for vname in lvname:
            vidx = vars.index(vname)
            atmp = sdat[:,vidx]
            atmp.dtype = fmts[vidx]
            if vname =='dtime':
                adtime = np.array([datetime(1801,1,1,0) + timedelta(minutes=i) for i in atmp.tolist()])
                lout.append(adtime.reshape(-1,npix)[:,0])  # 1-dimensional array
            else:
                lout.append(atmp.reshape(-1,npix))

        return lout
# %%
