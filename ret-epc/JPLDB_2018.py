import os, sys
from numpy import *
from collections import OrderedDict as odict
import struct
import numpy as np

'''
Reads JPL observational database format
February 2019

***************************************************************
** The notes below are for older version of the database. ***
** They should be updatad ***
***************************************************************

Example input file syntax   20150725_204633_GPM_007985.strm  is provided
In this example, the first record is at 2015/07/25 20:46:33 UTC    7985= GPM orbit revolution

The procedure and of the DPR-GMI matching and calculation of coeffficients to estimate the emissivity principal components (EPC)
from the GMI TB is documented in:

Turk, F.J., Haddad, Z.S. & You, Y., 2016, Estimating Non-Raining Surface Parameters to Assist
GPM Constellation Radiometer Precipitation Algorithms, Journal of Atmospheric and Oceanic Technology, 33(2016), pp. 1333-53

Use "-fpack-struct" gcc compiler option to pack all structure members together without holes
gcc -fpack-struct -o FILE FILE.c  (where FILE.c is the filename of this program)

To test, run as:  FILE 20150725_204633_GPM_007985.strm

sfc_class value:
1= ocean
2= sea ice
3-7 = decreasing vegetation
8-11 = decreasing snow cover
12 = inland water
13 = coast
14 = ocean/sea-ice boundary

NS refers to the retrievals done across the full Ku-band radar swath (49 beam positions), no consideration of Ka-band
MS refers to the retrievals done within NS beam positions 12-36, where both Ku and Ka band radars are available
HS refers to the retrievals done with the Ka-band High Sensivity radar (interleaved with Ka-band)

DESCRIPTION OF EACH RECORD

satid= identifier for the satellite with the radiometer  0=GPM 16-F16 17=F17 18=F18 19=F19 100=ATMS
rev= orbit revolution number for the satellite with the radar (always GPM)
rev2= orbit rev number for satellite with the radiometer (obviously rev2=rev for GMI, but rev2 != rev for SSMIS and ATMS)
SC_orientation= GPM spacecraft orientation, 0 or 180 degrees
i_S1,j_S1= scan and pixel number from radiometer.  For GMI, i_S1 ranges from 0-3000 or so. j_S1 ranges from 0-220.
i_NS,j_NS= scan and pixel number from radar.  For DPR, i_NS ranges from 0-9500 or so. j_NS ranges from 0-48.
sfc_class= TELSEM surface class index, listed above
yyyy,mm,dd,hh,mn,ss= date
timediff = time offset between radar and radiometer, in seconds
slat, slon= spacecraft coordinates of the satellite with the radar, here GPM
glat1, glon1= GMI coordinates
slat2, slon2= spacecraft coordinates of the satellite with the radiometer (slat=slat2 and slon=slon2 for GMI)
tb= GMI TB, from 1B.GPM.GMI, in Kelvin
pc_emis= emissivity principal components (empty, computed below), length 11 array
emis = emissivity that is reconstructed from pc_emis (empty, computed below), length 11 array
sfc_min, sfc_max= land surface type min and max values encountered from the DPR 3x3 profiles
elev = elevation in meters

nku10, nka10, nka10_HS= Number of bins where Z > 10 dBZ within 3x3 DPR profile, for Ku, Ka and Ka-HighSens radars
nku15, nka15, nka15_HS= Number of bins where Z > 15 dBZ within 3x3 DPR profile, for Ku, Ka and Ka-HighSens radars
nku20, nka20, nka20_HS= Number of bins where Z > 20 dBZ within 3x3 DPR profile, for Ku, Ka and Ka-HighSens radars
nku25, nka25, nka25_HS= Number of bins where Z > 25 dBZ within 3x3 DPR profile, for Ku, Ka and Ka-HighSens radars
zen_NS= DPR zenith angle, from 2A.GPM.DPR, in degrees
pia_NS, pia_MS= path integrated attenuation in dB for DPR Ku-band only (NS) and Ku+Ka (MS) retrievals
pia_NS_cmb= path integrated attenuation in dB for Ku-band only (NS) combined (DPR+GMI) retrievals
pia_MS_cmb[2]= path integrated attenuation in dB for Ku+Ka band (MS) combined (DPR+GMI) retrievals

precip_NS= Estimated surface precip rate from DPR-only NS retrieval, precipRateESurface in 2A.GPM.DPR
precip_NS_max= max rain rate encountered within the 3x3 profiles

precip2_NS= Near surface precip rate from DPR-only NS retrieval, precipRateNearSurface in 2A.GPM.DPR
precip2_NS_max= max precip rate encountered within the 3x3 profiles

precip_MS= Estimated surface precip rate from DPR-only MS retrieval, precipRateESurface in 2A.GPM.DPR
precip_MS_max= max rain rate encountered within the 3x3 profiles

precip2_MS= Near surface precip rate from DPR-only MS retrieval, precipRateNearSurface in 2A.GPM.DPR
precip2_MS_max= max precip rate encountered within the 3x3 profiles

precip_NS_cmb= Surface precip rate from DPR+GMI combined NS retrieval, surfPrecipTotRate in 2B.GPM.DPRGMI.CORRA
precip_NS_cmb_max= max precip rate encountered within the 3x3 profiles

precip_MS_cmb= Surface precip rate from DPR+GMI combined MS retrieval, surfPrecipTotRate in 2B.GPM.DPRGMI.CORRA
precip_MS_cmb_max= max precip rate encountered within the 3x3 profiles

precip_GPROF=  GPROF precip, surfacePrecipitation, copied over from 2A.GPM.GMI.GPROF
frozen_precip_GPROF= same as above but frozenPrecipitation, copied over from 2A.GPM.GMI.GPROF
prob_precip_GPROF= probability of precipitation, probabilityOfPrecip copied over from 2A.GPM.GMI.GPROF

ts, t2m= surface temperature (Kelvin), 2-m air temp (Kelvin), interpolated from MERRA2 1-hour reanalysis
tqv= total precip column water (mm), interpolated from MERRA2 3-hourly reanalysis
hs,ps= surface geopotential height (meters) and surface pressure (hPa), interpolated from MERRA2 3-hourly reanalysis
p_prof= pressure levels of MERRA2 (42 levels), interpolated from MERRA2 3-hourly reanalysis, scaled by 10
h_prof= geopotential height profile, in meters
t_prof= temperature profile, in Kelvin
qv_prof= specific humidity profile, in g/g

z_ku= Ku-band uncorrected radar reflectivity profile, zFactorMeasured in 2A.GPM.DPR, aggregated to 88 bins, scaled by 100
z_ka= Ka-band uncorrected radar reflectivity profile, zFactorMeasured in 2A.GPM.DPR, aggregated to 88 bins, scaled by 100
'''


class JPLDB(object):
    """
    Based on jpldb.py by Hyungjun Kim @ IIS, U-Tokyo, 2018-01-18
    Updated for 2018 version @ JPL, USA, 2019-02-19, by N.UTSUMI
    Updated for 2020-March varsion @ JPL, USA, 2020-03-01 by N.UTSUMI
    """
    def __init__(self, sensor): 
        #self.srcDir = "/home/utsumi/mnt/wellshare/data/JPLDB"

        NCHAN   = 13        # 13 channels of GMI Tb
        #NEM     = 11       # 9 emissivity [10V 10H 19V 19H 23V 37V 37H 89V 89H], Ts, Column vapor
        NEM     = 16  

        NLEV_MODEL= 42      # 42 levels of MERRA2
        NLEV_DPR= 88        # 88 vertical bins of DPR (250m each)
        NLEV_PRECIP= 22     # 22 lowest levels for preicp (250m each)
        NLEV_GPROF = 28

        dictvars    = odict((
                            ('satid',               'h'),           # PMW
                            ('satid2',              'h'),           # Radar:GPM
                            ('rev',                 'i'),           # PMW
                            ('rev2',                'i'),           # Radar:GPM
                            ('SC_orientation',      'h'),           #
                            ('SC_orientation2',     'h'),           #
                            ('i_S1',                'h'),           #
                            ('j_S1',                'h'),           #
                            ('i_NS',                'h'),           #
                            ('j_NS',                'h'),           #
                            ('yyyy',                'h'),           #
                            ('mm',                  'h'),           #
                            ('dd',                  'h'),           #
                            ('hh',                  'h'),           #
                            ('mn',                  'h'),           #
                            ('ss',                  'h'),           #
                            ('timediff',            'i'),           # Time difference at each pixel
                            ('kmdiff',              'f'),           # Distance between radar and PMW pixel [km]
                            ('glat',                'f'),           # Lat for PMW pixel
                            ('glon',                'f'),           #
                            ('slat',                'f'),           # Spacecraft lat for PMW
                            ('slon',                'f'),           #
                            ('salt',                'f'),           # Spacecraft altitude
                            ('slat2',               'f'),           # Spacecraft lat for radar
                            ('slon2',               'f'),           #
                            ('salt2',               'f'),           #

                            ('inc_S1',              'f'),           # PMW incidence angle (1)
                            ('inc_S2',              'f'),           # PMW incidence angle (2)
                            ('zen_NS',              'f'),           # Radar zenith angle
                            ('tb',                '%if'%NCHAN),     # Always 13, currently
                            ('tb_min',            '%if'%NCHAN),     # Always 13, currently
                            ('tb_max',            '%if'%NCHAN),     # Always 13, currently

                            ('pc_emis',           '%if'%NEM),       #
                            ('emis',              '%if'%NEM),       #
                            ('emis_NS_cmb',       '%if'%NCHAN),     # surfEmissivity from CMB product
                            ('s0_NS',               'f'),           # Sigma0 for NS (Ku)
                            ('s0_MS',               'f'),           # Sigma0 for MS (Ka)
                            ('sfc_class',           'h'),           # Surface class from GPROF
                            ('sfc_min',             'h'),           # Max of landSurfaceType in a PMW footprint from DRP product
                            ('sfc_max',             'h'),           #
                            ('elev',                'h'),           # Elevation of the middle pixel (from DPR product)

                            ('ndpr_NS',             'h'),           # Number of radar pixels in a PMW pixel (Ku)
                            ('ndpr_MS',             'h'),           # Number of radar pixels in a PMW pixel (Ka)

                            ('nku10',               'h'),           # Number of observation in 3-D radar (10<=Z<15, no attenuattion corrected)
                            ('nka10',               'h'),           # For screening cloudy scenes
                            ('nku15',               'h'),           # 15<=Z<20
                            ('nka15',               'h'),           #
                            ('nku20',               'h'),           # 20<=Z<25
                            ('nka20',               'h'),           #
                            ('nku25',               'h'),           # 25<=Z
                            ('nka25',               'h'),           #

                            ('pia_NS',              'f'),           # PIA from DPR product
                            ('pia_MS',              'f'),           #
                            ('pia_NS_cmb',          'f'),           #
                            ('pia_MS_cmb',         '2f'),           #


                            ('precip_nsfc_NS',      'f'),           # NS: Ku
                            ('precip_nsfc_max_NS',  'f'),           #
                            ('precip_nsfc_MS',      'f'),           # MS: Dural freq.(Ku & Ka)
                            ('precip_nsfc_max_MS',  'f'),           #
                            ('vfrac_conv_NS',       'f'),           #
                            ('vfrac_conv_MS',       'f'),           #

                            ('precip_NS_cmb',       'f'),           # Precip from Combined product: near surface, not esurf. (surfPrecipTotRate)
                            ('precip_max_NS_cmb',   'f'),           #
                            ('precip_MS_cmb',       'f'),           #
                            ('precip_max_MS_cmb',   'f'),           #
                            ('vfrac_conv_NS_cmb',   'f'),           #
                            ('vfrac_conv_MS_cmb',   'f'),           #

                            ('type_precip_NS',     '3h'),           # Number of radar pixels for each typePrecip (from DPR product) in PMW pixel.
                            ('shallow_rain_NS',    '5h'),           # Number of radar pixels for each flagShallowRain (from DPR product) in PMW pixel.
                            ('type_precip_MS',     '3h'),           #
                            ('shallow_rain_MS',    '5h'),           #
                            ('storm_height_ku',    'f'),            #
                            ('storm_height_max_ku','f'),            #
                            ('storm_height_ka',    'f'),            #
                            ('storm_height_max_ka','f'),            #

                            ('precip_GPROF',        'f'),           #
                            ('prob_precip_GPROF',   'f'),           #
                            ('frozen_precip_GPROF', 'f'),           #
                            ('conv_precip_GPROF',   'f'),           #
                            ('pixel_status_GPROF',  'h'),           #
                            ('qual_flag_GPROF',     'h'),           #
                            ('prof_GPROF',        '%ih'%NLEV_GPROF),   # Top to bottom (different from original GPROF), Scaled by 1000
                    

                            ('ts',                  'f'),           # MERRA2 (Interpolated to space and time)
                            ('t2m',                 'f'),           #
                            ('t2m_dew',             'f'),           #
                            ('t2m_wet',             'f'),           #
                            ('tqv',                 'f'),           #
                            ('hs',                  'f'),           #
                            ('ps',                  'f'),           #
                            ('u850',                'f'),           #
                            ('u500',                'f'),           #
                            ('u250',                'f'),           #
                            ('v850',                'f'),           #
                            ('v500',                'f'),           #
                            ('v250',                'f'),           #
                            ('p_prof',            '%ih'%NLEV_MODEL),#
                            ('h_prof',            '%if'%NLEV_MODEL),#
                            ('t_prof',            '%if'%NLEV_MODEL),#
                            ('qv_prof',           '%if'%NLEV_MODEL),#

                            ('bin_sfc_ku',          'h'),           # First surface bin above clutter (center pixel of PMW footprint)
                            ('bin_sfc_ka',          'h'),           #

                            ('precip_prof_NS',    '%ih'%NLEV_PRECIP),# Lowest 22 precipitation 250m-bins (i.e., two 125m bins are vertically averaged for DPR) (averaged over PMW pixel)
                            ('precip_prof_MS',    '%ih'%NLEV_PRECIP),# Scaled by 100
                            ('precip_prof_NS_cmb','%ih'%NLEV_PRECIP),# Orderd from higher to lower bins
                            ('precip_prof_MS_cmb','%ih'%NLEV_PRECIP),#
                            ('nclutter_ku',       '%ih'%NLEV_PRECIP),# Number of 3D pixels with clutter
                            ('nclutter_ka',       '%ih'%NLEV_PRECIP),#
                            ('z_ku',              '%ih'%NLEV_DPR),# Reflectivity Zm (zFactorMeasured) from DPR product, scaled by 100
                            ('z_ka',              '%ih'%NLEV_DPR),#
                            ('precip_water_prof_NS', '%ih'%NLEV_DPR), # Combined, Top to bottom, Scaled by 1000
                            ('precip_water_prof_MS', '%ih'%NLEV_DPR), # Combined, Top to bottom, Scaled by 1000


                           ))



        self.vars,      \
        self.fmts       = zip(*dictvars.items())
        self.dictvars   = dictvars.copy()

        #self.fmtsizes           = [struct.calcsize( fmt ) for fmt in self.fmts]
        self.fmtsize    = struct.calcsize( '<'+''.join(self.fmts) )



    #def __call__(self, dbName="GMI_dbase_PC_prof2_V5_NPC4_BIN10_OVERLAP1"):
    #    self.srcDir   = self.rootDir + "/" + dbName


    def __getattr__(self, k):

        if k in self.dictvars:
            return self.get_var( k )

        else:
            raise KeyError

    def set_file(self, srcPath):
        self.filesize   = os.stat( srcPath ).st_size
        self.nchunks    = self.filesize / self.fmtsize

        self.dbmmap     = np.memmap( srcPath, dtype='S1', mode='r',
                                     shape=(self.nchunks, self.fmtsize) )

        #self.curr       = 0



    def get_var( self, vname, nrec=None, origin=0 ):

        nrec    = self.nchunks if nrec == None       \
             else nrec

        vidx    = self.vars.index( vname )
        fmt     = self.fmts[ vidx ]
        size    = struct.calcsize( fmt )
        seek    = struct.calcsize( '<'+''.join(self.fmts[:vidx]) )

        data    = np.array( self.dbmmap[origin : origin+nrec,
                                        seek   : seek+size] )

        data.dtype  = fmt[-1]

        if data.shape[-1] == 1:
            data.shape  = data.shape[:-1]

        return data



