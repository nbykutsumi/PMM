import os, sys
from numpy import *
from collections import OrderedDict as odict
import struct
import numpy as np

class EPCDB(object):
    def __init__(self):

        NCHAN   = 13        # 13 channels of GMI Tb
        #NEM     = 11       # 9 emissivity [10V 10H 19V 19H 23V 37V 37H 89V 89H], Ts, Column vapor
        NEM     = 16

        NLEV_MODEL= 42      # 42 levels of MERRA2
        NLEV_DPR= 50        # 50 lowest bins of DPR (250m each)
        NLEV_PRECIP= 60     # 60 lowest levels for preicp (250m each)



        self.dictvars    = odict((
            ('satid',            ['',   'h']),                 # PMW
            ('satid2',           ['',   'h']),                 # Radar:GPM
            ('rev',              ['gNum',   'i']),             # PMW
            ('rev2',             ['',   'i']),                 # Radar:GPM
            ('SC_orientation',   ['',   'h']),                 #
            ('SC_orientation2',  ['',   'h']),                 #
            ('i_S1',             ['',   'h']),                 #
            ('j_S1',             ['',   'h']),                 #
            ('i_NS',             ['',   'h']),                 #
            ('j_NS',             ['',   'h']),                 #
            ('yyyy',             ['',   'h']),                 #
            ('mm',               ['',   'h']),                 #
            ('dd',               ['',   'h']),                 #
            ('hh',               ['',   'h']),                 #
            ('mn',               ['',   'h']),                 #
            ('ss',               ['',   'h']),                 #
            ('timediff',         ['',   'i']),                 # Time difference at each pixel
            ('kmdiff',           ['',   'f']),                 # Distance between radar and PMW pixel [km]
            ('glat',             ['',   'f']),                 # Lat for PMW pixel
            ('glon',             ['',   'f']),                 #
            ('slat',             ['',   'f']),                 # Spacecraft lat for PMW
            ('slon',             ['',   'f']),                 #
            ('salt',             ['',   'f']),                 # Spacecraft altitude
            ('slat2',            ['',   'f']),                 # Spacecraft lat for radar
            ('slon2',            ['',   'f']),                 # 
            ('salt2',            ['',   'f']),                 #

            ('inc_S1',           ['',   'f']),                 # PMW incidence angle (1)
            ('inc_S2',           ['',   'f']),                 # PMW incidence angle (2)
            ('zen_NS',           ['',   'f']),                 # Radar zenith angle
            ('tb',               ['Tc', '%if'%NCHAN]),         # Always 13, currently

            ('pc_emis',          ['epc', '%if'%NEM]),          #
            ('emis',             ['', '%if'%NEM]),             #
            ('emis_NS_cmb',      ['', '%if'%NCHAN]),           # surfEmissivity from CMB product
            ('s0_NS',            ['',   'f']),                 # Sigma0 for NS (Ku)
            ('s0_MS',            ['',   'f']),                 # Sigma0 for MS (Ka)
            ('sfc_class',        ['',   'h']),                 # Surface class from GPROF
            ('sfc_min',          ['',   'h']),                 # Max of landSurfaceType in a PMW footprint from DRP product
            ('sfc_max',          ['',   'h']),                 #
            ('elev',             ['gtopo',  'h']),             # Elevation of the middle pixel (from GTOPO30)

            ('ndpr_NS',          ['',   'h']),                 # Number of radar pixels in a PMW pixel (Ku)
            ('ndpr_MS',          ['',   'h']),                 # Number of radar pixels in a PMW pixel (Ka)

            ('nku10',            ['',   'h']),                 # Number of observation in 3-D radar (10<=Z<15, no attenuattion corrected)
            ('nka10',            ['',   'h']),                 # For screening cloudy scenes
            ('nku15',            ['',   'h']),                 # 15<=Z<20
            ('nka15',            ['',   'h']),                 #
            ('nku20',            ['',   'h']),                 # 20<=Z<25
            ('nka20',            ['',   'h']),                 #
            ('nku25',            ['',   'h']),                 # 25<=Z
            ('nka25',            ['',   'h']),                 #

            ('pia_NS',           ['',   'f']),                 # PIA from DPR product
            ('pia_MS',           ['',   'f']),                 #
            ('pia_NS_cmb',       ['',   'f']),                 #
            ('pia_MS_cmb',       ['',  '2f']),                 #


            ('precip_nsfc_NS',     ['Ku_NS_precipRateNearSurface', 'f']),   # NS: Ku
            ('precip_nsfc_max_NS', ['', 'f']),                 #
            ('precip_esfc_NS',     ['', 'f']),                 #
            ('precip_esfc_max_NS', ['', 'f']),                 #
            ('precip_nsfc_MS',     ['Ka_MS_precipRateNearSurface', 'f']),   # MS: Ka
            ('precip_nsfc_max_MS', ['', 'f']),                 #
            ('precip_esfc_MS',     ['', 'f']),                 #
            ('precip_esfc_max_MS', ['', 'f']),                 #

            ('precip_NS_cmb',      ['DPRGMI_NS_surfPrecipTotRate', 'f']),      # Precip from Combined product: near surface, not esurf. (surfPrecipTotRate)
            ('precip_max_NS_cmb',  ['', 'f']),                 #
            ('precip_MS_cmb',      ['DPRGMI_MS_surfPrecipTotRate', 'f']),                 #
            ('precip_max_MS_cmb',  ['', 'f']),                 #

            ('type_precip_NS',     ['','3h']),                 # Number of radar pixels for each typePrecip (from DPR product) in PMW pixel.
            ('shallow_rain_NS',    ['','5h']),                 # Number of radar pixels for each flagShallowRain (from DPR product) in PMW pixel.
            ('type_precip_MS',     ['','3h']),                 #
            ('shallow_rain_MS',    ['','5h']),                 #

            ('precip_GPROF',       ['', 'f']),                 #
            ('prob_precip_GPROF',  ['', 'f']),                 #
            ('frozen_precip_GPROF',['', 'f']),                 #

            ('ts',                 ['', 'f']),                 #
            ('t2m',                ['t2m', 'f']),              # MERRA2
            ('t2m_dew',            ['', 'f']),                 #
            ('t2m_wet',            ['', 'f']),                 #
            ('tqv',                ['tqv', 'f']),              # MERRA2
            ('hs',                 ['', 'f']),                 #
            ('ps',                 ['', 'f']),                 #
            ('u850',               ['', 'f']),                 #
            ('u500',               ['', 'f']),                 #
            ('u250',               ['', 'f']),                 #
            ('v850',               ['', 'f']),                 #
            ('v500',               ['', 'f']),                 #
            ('v250',               ['', 'f']),                 #
            ('p_prof',            ['', '%ih'%NLEV_MODEL]),     #
            ('h_prof',            ['', '%if'%NLEV_MODEL]),     #
            ('t_prof',            ['', '%if'%NLEV_MODEL]),     #
            ('qv_prof',           ['', '%if'%NLEV_MODEL]),     #

            ('bin_sfc_ku',        ['',  'h']),                 # First surface bin above clutter (center pixel of PMW footprint)
            ('bin_sfc_ka',        ['',  'h']),                 #

            ('precip_prof_NS',    ['', '%ih'%NLEV_PRECIP]),    # Lowest 22 precipitation 250m-bins (i.e., two 125m bins are vertically averaged for DPR) (averaged over PMW pixel)
            ('precip_prof_MS',    ['', '%ih'%NLEV_PRECIP]),    # Scaled by 100
            ('precip_prof_NS_cmb',['', '%ih'%NLEV_PRECIP]),    # Orderd from higher to lower bins
            ('precip_prof_MS_cmb',['', '%ih'%NLEV_PRECIP]),    #
            ('nclutter_ku',       ['', '%ih'%NLEV_PRECIP]),    # Number of 3D pixels with clutter
            ('nclutter_ka',       ['', '%ih'%NLEV_PRECIP]),    #
            ('z_ku',              ['Ku_NS_zFactorMeasured', '%ih'%NLEV_DPR]), # Reflectivity Zm (zFactorMeasured) from DPR product, scaled by 100
            ('z_ka',              ['Ka_MS_zFactorMeasured', '%ih'%NLEV_DPR]), #

            ('precip_water_prof_NS',    ['DPRGMI_NS_precipTotWaterCont', '%if'%NLEV_PRECIP]),  # precipitation water content profile from Comb product

            ))


        self.ddattype = {
                        'h':'int16'
                       ,'i':'int32'
                       ,'f':'float32'
                        }

        self.dmiss    = {
                        'h':-9999
                       ,'i':-9999
                       ,'f':-9999.9
                        }


    def set_idx_db(self, baseDir, idx_db):
        self.baseDir = baseDir
        self.idx_db  = idx_db

    #def set_file(self, srcPath):
    #    self.filesize   = os.stat( srcPath ).st_size
    #    self.nchunks    = self.filesize / self.fmtsize

    #    self.dbmmap     = np.memmap( srcPath, dtype='S1', mode='r',
    #                                 shape=(self.nchunks, self.fmtsize) )

    #    #self.curr       = 0



    def get_var(self, vname, nrec=None, origin=0):
        filevname, fmt = self.dictvars[vname]
        srcPath = self.baseDir + '/%s/%05d.npy'%(filevname, self.idx_db)
        if os.path.exists(srcPath):
      
            if nrec is None: 
                data = np.load(srcPath)
            else:
                data = np.load(srcPath,mmap_mode='r')
                data = np.array(data[origin:origin+nrec])

        else:
            ' read nrain file and return'
            nrainPath = self.baseDir + '/nrain/db_%05d.bin.nrain.txt'%(self.idx_db)
            f=open(nrainPath,'r'); lines=f.readlines(); f.close()
            line = lines[0].split()
            ndat = int(line[0])

            dattype = self.ddattype[fmt[-1]]
            miss    = self.dmiss[fmt[-1]]
            if nrec is not None:
                if ndat > nrec:
                    ndat = nrec

            if len(fmt)==1:
                data = np.ones(ndat, dtype=dattype) * miss
                nlev = 1
            else:
                nlev = int(fmt[:-1])
                data = np.ones([ndat,nlev], dtype=dattype) * miss

        return data 


    #def get_var(self, vname, nrec=None, origin=0):
    #    filevname, fmt = self.dictvars[vname]
    #    srcPath = self.baseDir + '/%s/%05d.npy'%(filevname, self.idx_db)
    #    if os.path.exists(srcPath):
    #   
    #        data = np.load(srcPath)
    #        if nrec is not None:
    #            data = data[origin:origin+nrec]

    #    else:
    #        ' read nrain file and return'
    #        nrainPath = self.baseDir + '/nrain/db_%05d.bin.nrain.txt'%(self.idx_db)
    #        f=open(nrainPath,'r'); lines=f.readlines(); f.close()
    #        line = lines[0].split()
    #        ndat = int(line[0])

    #        dattype = self.ddattype[fmt[-1]]
    #        miss    = self.dmiss[fmt[-1]]
    #        if nrec is not None:
    #            if ndat > nrec:
    #                ndat = nrec

    #        if len(fmt)==1:
    #            data = np.ones(ndat, dtype=dattype) * miss
    #            nlev = 1
    #        else:
    #            nlev = int(fmt[:-1])
    #            data = np.ones([ndat,nlev], dtype=dattype) * miss

    #    return data 


