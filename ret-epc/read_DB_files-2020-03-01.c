/*
Reads JPL observational database format
September 2017

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

*/
/*------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <strings.h>
#include <time.h>

#define NLEV_NS 88  /* DPR 250-m bins */
#define NLEV 42  /* from MERRA2  42 levels */
#define NCHAN 13  /* GMI TB */

#define NEM 16  /* 13= 13 emis then Ts then total column vapor T2m */
#define NTBREG 9
#define NREG 57

#define NEM_USE 4
#define NPCHIST 10


   typedef struct {
     short satid, satid2;
     int rev, rev2;
     short SC_orientation, SC_orientation2;
     short i_S1,j_S1,i_NS,j_NS;
     short yyyy,mm,dd,hh,mn,ss;
     int timediff;
     float kmdiff;
     float glat,glon;
     float slat,slon,salt,slat2,slon2,salt2;

     float inc_S1, inc_S2, zen_NS;
     float tb[13], tb_min[13], tb_max[13];

     float pc_emis[16], emis[16];
     float emis_NS_cmb[13];
     float s0_NS, s0_MS;
     short sfc_class,sfc_min,sfc_max,elev;

     short ndpr_NS, ndpr_MS;

     short nku10, nka10;
     short nku15, nka15;
     short nku20, nka20;
     short nku25, nka25;

     float pia_NS, pia_MS;
     float pia_NS_cmb, pia_MS_cmb[2];

     float precip_nsfc_NS, precip_nsfc_max_NS;
     float precip_nsfc_MS, precip_nsfc_max_MS;
     float vfrac_conv_NS, vfrac_conv_MS;

     float precip_NS_cmb, precip_max_NS_cmb;
     float precip_MS_cmb, precip_max_MS_cmb;
     float vfrac_conv_NS_cmb, vfrac_conv_MS_cmb;

     short type_precip_NS[3], shallow_rain_NS[5];
     short type_precip_MS[3], shallow_rain_MS[5];
     float stormtop_NS, stormtop_max_NS;
     float stormtop_MS, stormtop_max_MS;

     float precip_GPROF, prob_precip_GPROF, frozen_precip_GPROF, conv_precip_GPROF;
     short pixel_status_GPROF, qual_flag_GPROF;
     short prof_GPROF[28];

     float ts, t2m, t2m_dew, t2m_wet, tqv, hs, ps;
     float u850, u500, u250;
     float v850, v500, v250;
     short p_prof[42];
     float h_prof[42], t_prof[42], qv_prof[42];

     short bin_sfc_ku, bin_sfc_ka;

     short precip_prof_NS[22];
     short precip_prof_MS[22];
     short precip_prof_NS_cmb[22];
     short precip_prof_MS_cmb[22];
     short nclutter_ku[22];
     short nclutter_ka[22];
     short z_ku[88];
     short z_ka[88];
     short watercont_prof_NS_cmb[88];
     short watercont_prof_MS_cmb[88];
   } strm_record;



strm_record prof;


main(argc,argv)
int argc;
char *argv[];
{

  FILE *fdb;
  int i, j, j2, k, m, nbad, bytes, irec, nrec=0, nrec_all=0, ndbf=0, nf=0, nrain=0, i1, iclutter;
  float thick, qv, delp, kPa, abstemp, densair, rmsd, ht, z_ku, z_ka, hs1, ps1, qv1;
  char pdate[32];

  FILE *fpc;
  float ave_emis[NEM], std_emis[NEM];
  float ave_pc[NEM], std_pc[NEM];
  float b_all[NREG+1][NEM];  /* regression coeffs for all NTBREG */
  float tmp[NREG+1], psum;
  float pc_e[NEM], e[NEM], u[NEM][NEM];
  char svar[512];
  int id, k1, k2, kt;
  float coeff_disc[NEM], disc, disc_mid;
  float x1,x2,y1,y2,z1,z2, watercont;

  float pc_e_min[NEM], pc_e_max[NEM];
  float pcmin, pcmax, pc_range[NEM][NPCHIST][2], pc_lo, pc_hi;
  int idx, idxn[NEM_USE], found, idx2, i2;

  int n, ndb_files;
  int p, pow2[NEM_USE], ndb_files_used;

  int *nrec_dbfile;
  float gmi_freq[13]= {10.7,10.7,18.6,18.6,21.65,36.5,36.5,89.0,89.0,166.0,166.0,186.31,190.31};

  /*---------------------------------------------------------------*/

  bytes= sizeof(strm_record);
  printf("record length= %d\n", bytes);

  if (argc != 2 ) exit(-1);

  if ((fdb = fopen(argv[1],"r")) == 0) {
    printf("Unable to open %s\n", argv[1]);
    exit(-1);
  }

  nrec= nrain= 0;
  while ( fread(&prof,bytes,1,fdb) == 1 ) {
    nrec++;

    /*if ( prof.elev < 500 ) continue;*/
    /*if ( prof.precip_NS_cmb < 2.0 ) continue;*/

    /*n= 0;
    for (j= 0; j< 28; j++) {
      if ( prof.prof_GPROF[j] > 0 ) n++;
    }
    if ( n != 0 ) continue;*/

    /*---------------------------------------------------*/
    /* quality control on each DB record */


    hs1= ps1= qv1= -99.99;
    for (i= 1; i< NLEV; i++) {
      if ( prof.h_prof[i] < 0.0 ) {
        hs1= prof.h_prof[i-1];
        ps1= 0.1*prof.p_prof[i-1];
        qv1= prof.qv_prof[i-1];
        break;
      }
    }
    if ( hs1 < 0.0 ) { hs1= prof.h_prof[NLEV-1]; ps1= 0.1*prof.p_prof[NLEV-1]; qv1= prof.qv_prof[NLEV-1]; }
    /*if ( prof.hs > hs1 || prof.ps < ps1 || qv1 <= 0.0 ) {*/
    if ( prof.hs > hs1 || qv1 <= 0.0 ) {
      /*printf("hs ps= %f %f  Lowest geopotential= %f %f", prof.hs, prof.ps, hs1, ps1);*/
      for (j= 0; j< NLEV; j++) {
        printf("%2d %8.2f %7.2f %6.2f %11.8f\n",j,prof.h_prof[j],0.1*prof.p_prof[j],prof.t_prof[j], prof.qv_prof[j]);
      }
      printf("sfc  %8.3f %8.3f %8.3f %8.3f\n",prof.hs,prof.ps,prof.ts, prof.t2m);
      exit(1);
    }

    /*---------------------------------------------------*/

    sprintf(pdate,"%04d/%02d/%02d %02d:%02d:%02d", prof.yyyy, prof.mm, prof.dd, prof.hh, prof.mn, prof.ss);

/*
    printf("%6d ", nrec);
    printf("%6d ", prof.satid);
    printf("%s ", pdate);
    printf("%8.3f %8.3f ", prof.glat, prof.glon);
    printf("%6.2f %6.2f %6.2f ", prof.ts, prof.t2m, prof.tqv);
    printf("%8.3f ", prof.precip_nsfc_NS);
    for (i= 0; i< 3; i++) printf("%7.3f ", prof.pc_emis[i]);
    printf("\n");
    printf("MIN ");
    for (i= 0; i< 13; i++) printf("%7.2f ", prof.tb_min[i]);
    printf("\n");
    printf("    ");
    for (i= 0; i< 13; i++) printf("%7.2f ", prof.tb[i]);
    printf("\n");
    printf("MAX ");
    for (i= 0; i< 13; i++) printf("%7.2f ", prof.tb_max[i]);
    printf("\n");
    printf("INC=%8.3f %8.3f %8.3f\n", prof.inc_S1, prof.inc_S2, prof.zen_NS);
*/

    printf("\n\n\n");
    printf("%06d %06d %4d %4d %4d %4d ", prof.rev, prof.rev2, prof.i_S1, prof.j_S1, prof.i_NS, prof.j_NS);
    printf("%6d %6d ", prof.satid, prof.satid2);
    printf("%8.3f %8.3f ", prof.glat, prof.glon);
    printf("NDPR=%3d %3d %6.2f ", prof.ndpr_NS, prof.ndpr_MS, prof.inc_S1);
    printf("%5d %6.2f ", prof.timediff, prof.kmdiff);
    printf("%8.3f %8.3f %8.3f %6d ", prof.slat, prof.slon, prof.salt, prof.SC_orientation);
    printf("\n");
    printf("%8.3f %8.3f %8.3f %6d ", prof.slat2, prof.slon2, prof.salt2, prof.SC_orientation2);
    printf("\n");
    printf("%3d %4d %4d %5d ", prof.sfc_class, prof.sfc_min, prof.sfc_max, prof.elev);
    printf("NDPR=%3d %3d ", prof.ndpr_NS, prof.ndpr_MS);
    printf("INC=%6.2f %6.2f %6.2f ", prof.inc_S1, prof.inc_S2, prof.zen_NS);
    printf("SFCBIN=%3d %3d ", prof.bin_sfc_ku, prof.bin_sfc_ka);
    printf("\n");
    printf("%6.2f %6.2f %6.2f ", prof.ts, prof.t2m, prof.tqv);
    printf("%6.2f %6.2f ", prof.t2m_wet, prof.t2m_dew);
    printf("\n");
    printf("%6.2f %6.2f ", prof.u850, prof.v850);
    printf("%6.2f %6.2f ", prof.u500, prof.v500);
    printf("%6.2f %6.2f ", prof.u250, prof.v250);
    printf("\n");
    printf("NS=%8.3f %8.3f ", prof.precip_nsfc_NS, prof.precip_nsfc_max_NS);
    printf("%8.3f %8.3f ", prof.precip_NS_cmb, prof.precip_max_NS_cmb);
    printf("\n");
    printf("MS=%8.3f %8.3f ", prof.precip_nsfc_MS, prof.precip_nsfc_max_MS);
    printf("%8.3f %8.3f ", prof.precip_MS_cmb, prof.precip_max_MS_cmb);
    printf("\n");
    for (i= 0; i< 13; i++) printf("%6.2f ", prof.tb[i]);
    printf("\n");
    printf("MIN ");
    for (i= 0; i< 13; i++) printf("%7.2f ", prof.tb_min[i]);
    printf("\n");
    printf("MAX ");
    for (i= 0; i< 13; i++) printf("%7.2f ", prof.tb_max[i]);
    printf("\n");
    for (i= 0; i< 13; i++) printf("%6.3f ", prof.emis[i]);
    printf("\n");
    for (i= 0; i< 12; i++) printf("%7.3f ", prof.pc_emis[i]);
    printf("\n\n");


    /*continue;*/

    for (j= 0; j< NLEV; j++) {
      printf("%2d %8.2f %7.2f %6.2f %11.8f\n",j,prof.h_prof[j],0.1*prof.p_prof[j],prof.t_prof[j], prof.qv_prof[j]);
    }
    printf("\n");

    /* DPR radar reflectivity profiles */
    j2= 0;
    for (j= 0; j< NLEV_NS; j++) {
      if ( j >= 3*NLEV_NS/4 ) j2++;
      ht= 0.25*(NLEV_NS-j-1);
      printf("%3d %6.3f %6d  ", j, ht, prof.elev); 

      z_ku= 0.01*prof.z_ku[j]; 
      watercont= 0.001*prof.watercont_prof_NS_cmb[j]; 
      printf("%6.2f %6.3f ", z_ku, watercont);
      if ( prof.j_NS >= 12 && prof.j_NS <= 36 ) {
        z_ka= 0.01*prof.z_ka[j]; 
        watercont= 0.001*prof.watercont_prof_MS_cmb[j]; 
        printf("%6.2f %6.3f  ", z_ka, watercont);
      }

      if ( j >= 3*NLEV_NS/4 ) {
        printf("%3d ", prof.nclutter_ku[j2-1]);
        printf("%8.3f %8.3f ", 0.01*prof.precip_prof_NS[j2-1], 0.01*prof.precip_prof_NS_cmb[j2-1]);
        if ( prof.j_NS >= 12 && prof.j_NS <= 36 ) {
          printf("%3d ", prof.nclutter_ka[j2-1]);
          printf("%8.3f ", 0.01*prof.precip_prof_MS_cmb[j2-1]);
        }
      }
      printf("\n");
    }

    /* GPROF profiles */
    for (j= 0; j< 28; j++) {
      printf("%3d %6.3f \n", j, 0.001*prof.prof_GPROF[j]); 
    }

    printf("GPROF PixelStatus=%3d  QualityFlag=%3d ", prof.pixel_status_GPROF, prof.qual_flag_GPROF); 
    printf("GPROF %8.3f  F=%8.3f C=%8.3f\n", prof.precip_GPROF, prof.frozen_precip_GPROF, prof.conv_precip_GPROF);
    printf("SFC  DPR-NS=%8.3f %8.3f ", prof.precip_nsfc_NS, prof.precip_nsfc_max_NS);
    printf("DPR-MS=%8.3f %8.3f ", prof.precip_nsfc_MS, prof.precip_nsfc_max_MS);
    printf("CMB-NS %8.3f %8.3f ", prof.precip_NS_cmb, prof.precip_max_NS_cmb);
    printf("CMB-MS %8.3f %8.3f ", prof.precip_MS_cmb, prof.precip_max_MS_cmb);
    printf("\n");

    printf("S1=%4d NS=%4d  Ndpr=%3d %3d  STORMTOP  %8.2f %8.2f VolFracConv=%8.3f %8.3f ",
      prof.j_S1, prof.j_NS, prof.ndpr_NS, prof.ndpr_MS, prof.stormtop_NS, prof.stormtop_max_NS, prof.vfrac_conv_NS, prof.vfrac_conv_NS_cmb);
    printf("Type= ");
    for (j= 0; j< 3; j++) printf("%3d ", prof.type_precip_NS[j]);
    /*if ( prof.j_NS >= 12 && prof.j_NS <= 36 ) {*/
    if ( prof.ndpr_MS > 0 ) {
      printf(" %8.2f %8.2f VolFracConv=%8.3f %8.3f ", prof.stormtop_MS, prof.stormtop_max_MS, prof.vfrac_conv_MS, prof.vfrac_conv_MS_cmb);
      printf(" Type= ");
      for (j= 0; j< 3; j++) printf("%3d ", prof.type_precip_MS[j]);
    }
    printf("\n");

  }  /* next record */

  fclose(fdb);
  printf("Nrec=%d Nrain=%d\n", nrec, nrain);


  exit(0);

}

