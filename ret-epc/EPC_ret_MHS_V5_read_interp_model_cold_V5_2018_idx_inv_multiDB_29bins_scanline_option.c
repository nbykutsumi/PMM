/*
Reads PPS 1C production files

Run interp_model_EPC on the 1C file to get the model file

Performs EPC-based precip retrievals using 1C TB channels

2018/03/25
Added optional DPR and DPRGMI (combined) files to save DPR precip, Ku/Ka profiles and 20-dB height.
Added CCFADS output file for NS (Ku-band database) and MS (Ku/Ka-band database) retrievals

2018/09/12
Added 22-level estimated precip profiles 

Three run options:
1) Set isat_start and isat_end to the scanline range desired
2) Otherwise set tlat and tlon to be the target coordinates and process nlines on either side of the scanline closest to the target
3) Otherwise set tlat > 90 to process entire orbit

Compile as:
gcc -L/opt/local/lib -l netcdf -I /opt/local/include amsr2.o atms.o gmi.o mhs.o saphir.o ssmis.o -o  file file.c

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <netcdf.h>


#define RADEARTH 6371
#define RTD 57.29578
#define DTR 0.017453

#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(-1);}

#define NCLASS 15	/* surface index classification used in GPROF */
#define NLEV_DPR 88
#define NLEV_PRECIP 22	/* lowest 22 levels */

#define NEM_USE 3	/* 25 bins for the first 3 EPC */
#define NPCHIST 29

#define PRECIP1 0  /* for precip PDF */
#define PRECIP2 100
#define DPRECIP 0.5
#define NPRECIP 200

#define ZMIN 18   /* For CCFADS */
#define ZMAX 58
#define DZ 0.5
#define NZ 80

#define HTMIN 0   /* For CCFADS */
#define HTMAX 20
#define DHT 0.5
#define NHT 40

/*-------------------------------------------*/
#define DB_MAXREC 10000   	/* max records to read from DB files */
#define DB_MINREC 100      	/* If DB_MAXREC not reached for initial DB index, min records needed from expanded search */
#define MAX_DB_EXPAND 20	/* expand database index search this far in + and - directions, in order to get at least DB_MINREC entries */
#define DB_RAINFRAC 0.10	/* If fraction(DB) > 1 mm/hr is less than this, set all pixels to R=0 */
#define N2MAX 200		/* Max records to weight after sorting */
#define WT2MIN 0.01		/* Exit when weight reaches this level */
#define MAX_T2M_DIFF 20		/* If have_model, then max T2m abs difference between observed and DB pixel */
#define MAX_TB_RMSD 50		/* max TB RMS difference between observed and DB pixel */
#define DB_MAX_INC_DIFF 20	/* max incidence angle difference for crosstrack sounders */
/*-------------------------------------------*/

/* record structure (1728 bytes) */

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
     float tb[13];

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
     float precip_esfc_NS, precip_esfc_max_NS;
     float precip_nsfc_MS, precip_nsfc_max_MS;
     float precip_esfc_MS, precip_esfc_max_MS;

     float precip_NS_cmb, precip_max_NS_cmb;
     float precip_MS_cmb, precip_max_MS_cmb;

     short type_precip_NS[3], shallow_rain_NS[5];
     short type_precip_MS[3], shallow_rain_MS[5];

     float precip_GPROF, prob_precip_GPROF, frozen_precip_GPROF;

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
   } strm_record;

strm_record prof, prof_idx[DB_MAXREC];

int *tmp1;
int *tmp2;
float *tmp3;
float *tmp7;

int max (int *a, int n, int i, int j, int k) {
    int m = i;
    if (j < n && a[j] > a[m]) {
        m = j;
    }
    if (k < n && a[k] > a[m]) {
        m = k;
    }
    return m;
}
 
void downheap (int *a, int *b, int n, int i) {
    while (1) {
        int j = max(a, n, i, 2 * i + 1, 2 * i + 2);
        if (j == i) {
            break;
        }
        int t = a[i];
        int t2 = b[i];
        a[i] = a[j];
        a[j] = t;
        b[i] = b[j];
        b[j] = t2;
        i = j;
    }
}
 
void heapsort2 (int *a, int *b, int n) {
    int i;
    for (i = (n - 2) / 2; i >= 0; i--) {
        downheap(a, b, n, i);
    }
    for (i = 0; i < n; i++) {
        int t = a[n - i - 1];
        int t2 = b[n - i - 1];
        a[n - i - 1] = a[0];
        a[0] = t;
        b[n - i - 1] = b[0];
        b[0] = t2;
        downheap(a, b, n - i - 1, 0);
    }
}
 

/*------------------ sensors ----------------------*/
char *em_names_atms[16]= {"e_23.8QV", "e_31.4QV", "e_88.2QV", "TS", "TQV", "T2M", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-"};
char *tb_names_atms[13]= {"23.8QV", "31.4QV", "88.2QV", "165.5QH", "183.31+/-7QH", "183.31+/-4.5QH", "183.31+/-3QH", "183.31+/-1.8QH", "183.31+/-1QH", "-", "-", "-", "-"};
int em_compare_atms[16]= {1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0};
int tb_compare_atms[13]= {1,1,1,1,1,1,1,1,1,0,0,0,0};

char *em_names_ssmis[16]= {"e_19.35V", "e_19.35H", "e_22.235V", "e_37.1V", "e_37.1H", "e_91.7V", "e_91.7H", "TS", "TQV", "T2M", "-", "-", "-", "-", "-", "-"};
char *tb_names_ssmis[13]= {"19.35V", "19.35H", "22.235V", "37.1V", "37.1H", "91.0V", "91.0H", "150H", "183.31+/-1", "183.31+/-3", "183.31+/-7", "-", "-"}; 
int em_compare_ssmis[16]= {1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0};
int tb_compare_ssmis[13]= {1,1,1,1,1,1,1,0,0,0,0,0,0};	

char *em_names_amsr2[16]= {"e_10.7V", "e_10.7H", "e_18.7V", "e_18.7H", "e_23.8V", "e_23.8H", "e_36.5V", "e_36.5H", "e_89.0V", "e_89.0H", "TS", "TQV", "T2M", "-", "-", "-"};
char *tb_names_amsr2[13]= {"10.7V", "10.7H", "18.7V", "18.7H", "23.8V", "23.8H", "36.5V", "36.5H", "89.0V", "89.0H", "-", "-", "-"};
int em_compare_amsr2[16]= {1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0};
int tb_compare_amsr2[13]= {1,1,1,1,1,1,1,1,1,1,0,0,0};

char *em_names_gmi[16]= {"e_10.7V", "e_10.7H", "e_18.7V", "e_18.7H", "e_23.8V", "e_36.5V", "e_36.5H", "e_89.0V", "e_89.0H", "TS", "TQV", "T2M", "-", "-", "-", "-"};
char *tb_names_gmi[13]= {"10.7V", "10.7H", "18.7V", "18.7H", "23.8V", "36.5V", "36.5H", "89.0V", "89.0H", "166V", "166H", "183.31+/-3V", "183.31+/-7V"};
int em_compare_gmi[16]= {1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0};
int tb_compare_gmi[13]= {1,1,1,1,1,1,1,1,1,1,1,1,1};

char *em_names_mhs[16]= {"e_89V", "TS", "TQV", "T2M", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-"};
char *tb_names_mhs[13]= {"89V", "157V", "183.31+/-1H", "183.31+/-3H", "190.31V", "-", "-", "-", "-", "-", "-", "-", "-"};
int em_compare_mhs[16]= {1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0};
int tb_compare_mhs[13]= {1,1,1,1,1,0,0,0,0,0,0,0,0};

char *em_names_saphir[16]= {"e_183.31+/-11.0", "TS", "TQV", "T2M", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-"};
char *tb_names_saphir[13]= {"183.31+/-0.2", "183.31+/-1.1", "183.31+/-2.8", "183.31+/-4.2", "183.31+/-6.8", "183.31+/-11.0", "-", "-", "-", "-", "-", "-", "-"};
int em_compare_saphir[16]= {1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0};
int tb_compare_saphir[13]= {1,1,1,1,1,1,0,0,0,0,0,0,0};

void tc_atms();
void tc_ssmis();
void tc_amsr2();
void tc_gmi();
void tc_mhs();
void tc_saphir();
/*----------------------------------------------*/




/*-------------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------*/

int main (argc, argv)
int argc;
char *argv[];

{
   int ncid;
   int grpid1;
   int grpid1_scanstatus, grpid1_scantime, grpid1_nav;
   int grpid2;
   int grpid2_scanstatus;
   int grpid_NS, grpid_MS;
   int grpid3;
   int grpid4;
   int grpid5;
   int grpid6;

   FILE *flist;
   char tbfile[256];

   int NEM, NCHAN, NTBREG, NTB, NREG;
   short em_compare[16]= {0}, tb_compare[13]= {0};
   float epc_min[16]= {1.0E6}, epc_max[16]= {-1.0E6};
   char dbdir[256], dbfile[256];
   char *tb_names[13], *em_names[16];

   int i, j, k, m, n, i1,i2,j1, j2, k1, status, isat, jsat, ip;
   int isat2, isat2_min, jsat2, jsat2_min, itest;
   int isat_start=-1, isat_end=-1, nsat=0;
   size_t nl, ns;
   time_t t0=-1.0, t1, t2;
   struct tm tm, tm0, tm1, tm2;
   char tdate[128], pdate[128], cdate[128], outdate[32], mdate[32];
   char cdate1[128], cdate2[128];
   char cdate12[128];
   long secs, secs1, secs2;


   /*------------------------------ 1C-------------------------------------------*/


   int sclat_S1_id, sclon_S1_id, FractionalGranuleNumber_S1_id, SCorientation_S1_id;
   float *sclat_S1_buf, *sclon_S1_buf;
   double *FractionalGranuleNumber_S1_buf;
   short *SCorientation_S1_buf;

   /*--- geolocation for each channel set ---*/
   int lat_S1_id, lat_S2_id, lat_S3_id, lat_S4_id, lat_S5_id, lat_S6_id;
   float *lat_S1_buf, *lat_S2_buf, *lat_S3_buf, *lat_S4_buf, *lat_S5_buf, *lat_S6_buf;
   int lon_S1_id, lon_S2_id, lon_S3_id, lon_S4_id, lon_S5_id, lon_S6_id;
   float *lon_S1_buf, *lon_S2_buf, *lon_S3_buf, *lon_S4_buf, *lon_S5_buf, *lon_S6_buf;

   /*--- incidence angle for each channel set ---*/
   int incidenceAngle_S1_id, incidenceAngle_S2_id, incidenceAngle_S3_id;
   int incidenceAngle_S4_id, incidenceAngle_S5_id, incidenceAngle_S6_id;
   float *incidenceAngle_S1_buf, *incidenceAngle_S2_buf, *incidenceAngle_S3_buf;
   float *incidenceAngle_S4_buf, *incidenceAngle_S5_buf, *incidenceAngle_S6_buf;

   /*--- Tc for the channel sets ---*/
   int Tc_S1_id, Tc_S2_id, Tc_S3_id, Tc_S4_id, Tc_S5_id, Tc_S6_id;
   float *Tc_S1_buf, *Tc_S2_buf, *Tc_S3_buf, *Tc_S4_buf, *Tc_S5_buf, *Tc_S6_buf;
   nc_type Tc_S1_type, Tc_S2_type, Tc_S3_type, Tc_S4_type, Tc_S5_type, Tc_S6_type;
   int Tc_S1_ndims, Tc_S2_ndims, Tc_S3_ndims, Tc_S4_ndims, Tc_S5_ndims, Tc_S6_ndims;
   int Tc_S1_dimids[NC_MAX_VAR_DIMS], Tc_S2_dimids[NC_MAX_VAR_DIMS], Tc_S3_dimids[NC_MAX_VAR_DIMS];
   int Tc_S4_dimids[NC_MAX_VAR_DIMS], Tc_S5_dimids[NC_MAX_VAR_DIMS], Tc_S6_dimids[NC_MAX_VAR_DIMS];
   int Tc_S1_natts, Tc_S2_natts, Tc_S3_natts, Tc_S4_natts, Tc_S5_natts, Tc_S6_natts;
   float *Tc_buf;

   int qual_S1_id, qual_S2_id, qual_S3_id, qual_S4_id, qual_S5_id, qual_S6_id;
   char *qual_S1_buf, *qual_S2_buf, *qual_S3_buf, *qual_S4_buf, *qual_S5_buf, *qual_S6_buf;

   int yyyy_S1_id, mm_S1_id, dd_S1_id, hh_S1_id, mn_S1_id, ss_S1_id;
   short *yyyy_S1_buf;
   short *mm_S1_buf;
   short *dd_S1_buf;
   short *hh_S1_buf;
   short *mn_S1_buf;
   short *ss_S1_buf;
   double *secs_S1_buf;

   size_t nl_S1, ns_S1, nc_S1;
   size_t nl_S2, ns_S2, nc_S2;
   size_t nl_S3, ns_S3, nc_S3;
   size_t nl_S4, ns_S4, nc_S4;
   size_t nl_S5, ns_S5, nc_S5;
   size_t nl_S6, ns_S6, nc_S6;

   long secs1_S1, secs2_S1;
   float tb_atms[13], tb_reg[13];
   float tb_atms_min[13], tb_atms_max[13];
   char gdate[32], sdate[32], ddate[32];

   float incidenceAngle;
   int irev, isatname;
   char satname[32], sensor[32];
   char *header;
   size_t att_len;
   unsigned short *att_val;

   int nbad;
   float percent;

   int *pixel_index, *n_index, *ncold_index;
   int npix=0;
   float ts_MERRA2_min= 1.0E6, ts_MERRA2_max= -1.0E6;
   float tqv_MERRA2_min= 1.0E6, tqv_MERRA2_max= -1.0E6;
   float t2m_MERRA2_min= 1.0E6, t2m_MERRA2_max= -1.0E6;

   float pixres, pixres_nadir= -1.0;
   int ndpr, ndpr_total, npixres;
  char svar[512], epc_coef_dir[512], outdir[512];




   /*--------------------------------------------------------------------------------*/

   if ( argc < 7 ) {
     printf("Usage: %s <XCAL-1C-FILE> target_lat target_lon scan_start scan_end OUTDIR <MERRA2-FILE> <GPROF-FILE> <DPR-FILE> <CMB-FILE> \n", argv[0]);
     exit(1);
   }


   /*--- Open 1C Tc file and read the header to determine the sensor  ---*/

   sprintf(tbfile,"%s", argv[1]);
   status = nc_open(tbfile, NC_NETCDF4, &ncid);
   if ( status != 0 ) {
     printf("File not found %s\n", tbfile);
     exit(1);
   }
   printf("Opened %s\n", tbfile);

   /* dumb way to get the satellite name */
   if ((status= nc_inq_attlen(ncid, NC_GLOBAL, "FileHeader", &att_len))) ERR(status);
   header= (char *) malloc(att_len*sizeof(char));
   if ((status= nc_get_att_text(ncid, NC_GLOBAL, "FileHeader", header))) ERR(status);

   isatname= -1;
   if ( isatname < 0 && strstr(header, "SatelliteName=TRMM") != NULL ) {isatname= 0; sprintf(satname,"%s","TRMM"); sprintf(sensor,"%s","TMI");}
   if ( isatname < 0 && strstr(header, "SatelliteName=GPM") != NULL ) {isatname= 1; sprintf(satname,"%s","GPM"); sprintf(sensor,"%s","GMI");}
   if ( isatname < 0 && strstr(header, "SatelliteName=F16") != NULL ) {isatname= 16; sprintf(satname,"%s","F16"); sprintf(sensor,"%s","SSMIS");}
   if ( isatname < 0 && strstr(header, "SatelliteName=F17") != NULL ) {isatname= 17; sprintf(satname,"%s","F17"); sprintf(sensor,"%s","SSMIS");}
   if ( isatname < 0 && strstr(header, "SatelliteName=F18") != NULL ) {isatname= 18; sprintf(satname,"%s","F18"); sprintf(sensor,"%s","SSMIS");}
   if ( isatname < 0 && strstr(header, "SatelliteName=F19") != NULL ) {isatname= 19; sprintf(satname,"%s","F19"); sprintf(sensor,"%s","SSMIS");}
   if ( isatname < 0 && strstr(header, "SatelliteName=GCOMW1") != NULL ) {isatname= 20; sprintf(satname,"%s","GCOMW1"); sprintf(sensor,"%s","AMSR2");}
   if ( isatname < 0 && strstr(header, "SatelliteName=NPP") != NULL ) {isatname= 100; sprintf(satname,"%s","NPP"); sprintf(sensor,"%s","ATMS");}
   if ( isatname < 0 && strstr(header, "SatelliteName=NOAA20") != NULL ) {isatname= 101; sprintf(satname,"%s","NOAA20"); sprintf(sensor,"%s","ATMS");}
   if ( isatname < 0 && strstr(header, "SatelliteName=METOPA") != NULL ) {isatname= 201; sprintf(satname,"%s","METOPA"); sprintf(sensor,"%s","MHS");}
   if ( isatname < 0 && strstr(header, "SatelliteName=METOPB") != NULL ) {isatname= 202; sprintf(satname,"%s","METOPB"); sprintf(sensor,"%s","MHS");}
   if ( isatname < 0 && strstr(header, "SatelliteName=NOAA18") != NULL ) {isatname= 318; sprintf(satname,"%s","NOAA18"); sprintf(sensor,"%s","MHS");}
   if ( isatname < 0 && strstr(header, "SatelliteName=NOAA19") != NULL ) {isatname= 319; sprintf(satname,"%s","NOAA19"); sprintf(sensor,"%s","MHS");}
   if ( isatname < 0 && strstr(header, "SatelliteName=MT1") != NULL ) {isatname= 400; sprintf(satname,"%s","MT1"); sprintf(sensor,"%s","SAPHIR");}

   if ((status= nc_close(ncid))) ERR(status);
   printf("Closed %s\n", tbfile);

   if ( isatname < 0 ) {
     printf("SatelliteName not recognized anywhere in \n\n%s\n", header);
     exit(1);
   }
   printf("%s  satname=%s sensor=%s\n", tbfile, satname, sensor);

   if ( strcmp(sensor, "ATMS") == 0 ) {
     NEM= 6;
     NCHAN= 9;
     NTBREG= 9;
     /*sprintf(dbdir,"%s", "/Volumes/GDRIVE4/projects/EPC_DB/ATMS_EPC_DATABASE");*/
     /*sprintf(dbdir,"%s", "/Volumes/15TB-RAID/projects/EPC_DB/ATMS_EPC_DATABASE");*/
     sprintf(dbdir,"%s", "/Volumes/6TB-1/projects/EPC_DB/ATMS_EPC_DATABASE_TEST29");
     sprintf(epc_coef_dir,"%s", "/Users/jturk/code/GPM/EPC_COEF/ATMS");
     for (i= 0; i< NCHAN; i++) tb_compare[i]= tb_compare_atms[i];
     for (i= 0; i< NCHAN; i++) tb_names[i]= tb_names_atms[i];
     for (i= 0; i< NEM; i++) em_names[i]= em_names_atms[i];
     for (i= 0; i< NEM; i++) em_compare[i]= em_compare_atms[i];
     pixres_nadir= 30.0;
   }
   else if ( strcmp(sensor, "SSMIS") == 0 ) {
     NEM= 10;
     NCHAN= 11;
     NTBREG= 7;
     /*sprintf(dbdir,"%s", "/Volumes/15TB-RAID/projects/EPC_DB/SSMIS_EPC_DATABASE");*/
     sprintf(dbdir,"%s", "/Volumes/GDRIVE4/projects/EPC_DB/SSMIS_EPC_DATABASE");
     /*sprintf(dbdir,"%s", "/Users/jturk/projects/EPC_DB/SSMIS_EPC_DATABASE");*/
     sprintf(epc_coef_dir,"%s", "/Users/jturk/code/GPM/EPC_COEF/SSMIS");
     for (i= 0; i< NCHAN; i++) tb_compare[i]= tb_compare_ssmis[i];
     for (i= 0; i< NCHAN; i++) tb_names[i]= tb_names_ssmis[i];
     for (i= 0; i< NEM; i++) em_names[i]= em_names_ssmis[i];
     for (i= 0; i< NEM; i++) em_compare[i]= em_compare_ssmis[i];
     pixres= 25.0;
   }
  else if ( strcmp(sensor, "AMSR2") == 0 ) {
     NEM= 13;
     NCHAN= 10;
     NTBREG= 10;
     sprintf(dbdir,"%s", "/Volumes/6TB-1/projects/EPC_DB/AMSR2_EPC_DATABASE");
     sprintf(epc_coef_dir,"%s", "/Users/jturk/code/GPM/EPC_COEF/AMSR2");
     for (i= 0; i< NCHAN; i++) tb_compare[i]= tb_compare_amsr2[i];
     for (i= 0; i< NCHAN; i++) tb_names[i]= tb_names_amsr2[i];
     for (i= 0; i< NEM; i++) em_names[i]= em_names_amsr2[i];
     for (i= 0; i< NEM; i++) em_compare[i]= em_compare_amsr2[i];
     pixres= 15.0;
   }
  else if ( strcmp(sensor, "GMI") == 0 ) {
     NEM= 12;
     NCHAN= 13;
     NTBREG= 13;
     /*sprintf(dbdir,"%s", "/Volumes/6TB-1/projects/EPC_DB/GMI_EPC_DATABASE");*/
     sprintf(dbdir,"%s", "/Volumes/6TB-1/projects/EPC_DB/GMI_EPC_DATABASE_TEST29");
     sprintf(epc_coef_dir,"%s", "/Users/jturk/code/GPM/EPC_COEF/GMI");
     for (i= 0; i< NCHAN; i++) tb_compare[i]= tb_compare_gmi[i];
     for (i= 0; i< NCHAN; i++) tb_names[i]= tb_names_gmi[i];
     for (i= 0; i< NEM; i++) em_names[i]= em_names_gmi[i];
     for (i= 0; i< NEM; i++) em_compare[i]= em_compare_gmi[i];
     pixres= 15.0;
   }
  else if ( strcmp(sensor, "MHS") == 0 ) {
     NEM= 4;
     NCHAN= 5;
     NTBREG= 5;
     sprintf(dbdir,"%s", "/Volumes/GDRIVE4/projects/EPC_DB/MHS_EPC_DATABASE");
     sprintf(epc_coef_dir,"%s", "/Users/jturk/code/GPM/EPC_COEF/MHS");
     for (i= 0; i< NCHAN; i++) tb_compare[i]= tb_compare_mhs[i];
     for (i= 0; i< NCHAN; i++) tb_names[i]= tb_names_mhs[i];
     for (i= 0; i< NEM; i++) em_names[i]= em_names_mhs[i];
     for (i= 0; i< NEM; i++) em_compare[i]= em_compare_mhs[i];
     pixres_nadir= 15.0;
   }
  else if ( strcmp(sensor, "SAPHIR") == 0 ) {
     NEM= 4;
     NCHAN= 6;
     NTBREG= 6;
     sprintf(dbdir,"%s", "/Volumes/GDRIVE4/projects/EPC_DB/SAPHIR_EPC_DATABASE");
     sprintf(epc_coef_dir,"%s", "/Users/jturk/code/GPM/EPC_COEF/SAPHIR");
     for (i= 0; i< NCHAN; i++) tb_compare[i]= tb_compare_saphir[i];
     for (i= 0; i< NCHAN; i++) tb_names[i]= tb_names_saphir[i];
     for (i= 0; i< NEM; i++) em_names[i]= em_names_saphir[i];
     for (i= 0; i< NEM; i++) em_compare[i]= em_compare_saphir[i];
     pixres_nadir= 15.0;
   }
   else {
     printf("%s not supported\n", sensor);
     exit(1);
   }

   n= 0;
   for (i= 0; i< NCHAN; i++) n+= tb_compare[i];
   printf("Sensor=%s  TB channels used=%d  NCHAN=%d\n", sensor, n, NCHAN);
   if ( n != NTBREG ) exit(1);

   NREG= NTBREG + (NTBREG*NTBREG-NTBREG)/2 + NTBREG + (NTBREG-1);
   printf("Sensor=%s  NEM=%d  NCHAN=%d  NTBREG=%d  NREG=%d\n", sensor, NEM, NCHAN, NTBREG, NREG);
   for (i= 0; i< NCHAN; i++) printf("%s ", tb_names[i]); printf("\n");
   for (i= 0; i< NEM; i++) printf("%s ", em_names[i]); printf("\n");
   printf("EPC DB directory=%s\n", dbdir);

   /*------------------------------ GPROF -------------------------------------------*/

   char gprof_file[256];
   int have_gprof= -1;
   int qualityFlag_SG_id, surfacePrecipitation_SG_id, probabilityOfPrecip_SG_id;
   int sfc_SG_id, totalColumnWaterVaporIndex_SG_id, temp2mIndex_SG_id;
   int frozenPrecipitation_SG_id;
   signed char *qualityFlag_GPROF_buf;
   signed char *sfc_GPROF_buf;
   float *surfacePrecipitation_GPROF_buf, *probabilityOfPrecip_GPROF_buf;
   float *lat_GPROF_buf, *lon_GPROF_buf;
   signed char *totalColumnWaterVaporIndex_GPROF_buf;
   float *frozenPrecipitation_GPROF_buf;
   short *temp2mIndex_GPROF_buf;

   signed char *sfc_GPROF_S1_buf;
   float *surfacePrecipitation_GPROF_S1_buf, *probabilityOfPrecip_GPROF_S1_buf;
   signed char *totalColumnWaterVaporIndex_GPROF_S1_buf;
   float *frozenPrecipitation_GPROF_S1_buf;
   short *temp2mIndex_GPROF_S1_buf;

   nc_type lat_SG_type, lon_SG_type;
   size_t nl_SG, ns_SG;
   int lat_SG_id, lon_SG_id, lat_SG_ndims, lat_SG_dimids[NC_MAX_VAR_DIMS], lat_SG_natts;

   float precip_GPROF, prob_precip_GPROF, frozen_precip_GPROF;
   int tqv_GPROF, t2m_GPROF;
   int t2m_GPROF_varid, tqv_GPROF_varid;
   int sfc_GPROF_varid;

   /*------------------------------ 2A.GPM.DPR -------------------------------------------*/

   char dpr_file[256];
   int have_dpr= -1;
   size_t nl_NS, ns_NS, nlev_NS;
   size_t nl_MS, ns_MS, nlev_MS;
   int grpid1_pre, grpid1_slv, grpid1_csf;
   int grpid2_pre, grpid2_slv;


   /*--- NS grid ---*/
   int lat_NS_id, lon_NS_id;
   float *lat_NS_buf, *lon_NS_buf;
   nc_type lat_NS_type, lon_NS_type;
   int lat_NS_ndims, lon_NS_ndims;
   int lat_NS_dimids[NC_MAX_VAR_DIMS], lon_NS_dimids[NC_MAX_VAR_DIMS];
   int lat_NS_natts, lon_NS_natts;

   int elev_NS_id;  /* elevation */
   float *elev_NS_buf;
   nc_type elev_NS_type;
   int elev_NS_ndims;
   int elev_NS_dimids[NC_MAX_VAR_DIMS];
   int elev_NS_natts;

   int zFactorMeasured_NS_id;  /* Z uncorrected */
   float *zFactorMeasured_NS_buf;
   nc_type zFactorMeasured_NS_type;
   int zFactorMeasured_NS_ndims;
   int zFactorMeasured_NS_dimids[NC_MAX_VAR_DIMS];
   int zFactorMeasured_NS_natts;

   int binRealSurface_NS_id;
   short *binRealSurface_NS_buf;
   int binClutterFreeBottom_NS_id;
   short *binClutterFreeBottom_NS_buf;

   int precipRateESurface_NS_id;  
   float *precipRateESurface_NS_buf;

   int flagShallowRain_NS_id, typePrecip_NS_id;
   int *flagShallowRain_NS_buf, *typePrecip_NS_buf;

   int yyyy_NS_id, mm_NS_id, dd_NS_id, hh_NS_id, mn_NS_id, ss_NS_id;
   short *yyyy_NS_buf;
   short *mm_NS_buf;
   short *dd_NS_buf;
   short *hh_NS_buf;
   short *mn_NS_buf;
   short *ss_NS_buf;
   double *secs_NS_buf;

   int qual_NS_id; 
   signed char *qual_NS_buf;


   /*--- MS grid ---*/
   int lat_MS_id, lon_MS_id;
   float *lat_MS_buf, *lon_MS_buf;
   nc_type lat_MS_type, lon_MS_type;
   int lat_MS_ndims, lon_MS_ndims;
   int lat_MS_dimids[NC_MAX_VAR_DIMS], lon_MS_dimids[NC_MAX_VAR_DIMS];
   int lat_MS_natts, lon_MS_natts;

   int zFactorMeasured_MS_id;  /* Z uncorrected */
   float *zFactorMeasured_MS_buf;
   nc_type zFactorMeasured_MS_type;
   int zFactorMeasured_MS_ndims;
   int zFactorMeasured_MS_dimids[NC_MAX_VAR_DIMS];
   int zFactorMeasured_MS_natts;

   int binRealSurface_MS_id;
   short *binRealSurface_MS_buf;
   int binClutterFreeBottom_MS_id;
   short *binClutterFreeBottom_MS_buf;

   int precipRateESurface_MS_id;  
   float *precipRateESurface_MS_buf;

   int qual_MS_id;  /* scan data quality flags */
   signed char *qual_MS_buf;

   /*------------------------------ COMBINED  -------------------------------------------*/

   char cmb_file[256];
   int have_cmb= -1;
   size_t nl_NS_cmb, ns_NS_cmb, nlev_NS_cmb;
   size_t nl_MS_cmb, ns_MS_cmb, nlev_MS_cmb;

   int PrecipTotRate_NS_cmb_id, PrecipTotRate_MS_cmb_id;
   float *PrecipTotRate_NS_cmb_buf, *PrecipTotRate_MS_cmb_buf;
   nc_type PrecipTotRate_NS_cmb_type, PrecipTotRate_MS_cmb_type;
   int PrecipTotRate_NS_cmb_ndims, PrecipTotRate_MS_cmb_ndims;
   int PrecipTotRate_NS_cmb_dimids[NC_MAX_VAR_DIMS], PrecipTotRate_MS_cmb_dimids[NC_MAX_VAR_DIMS];
   int PrecipTotRate_NS_cmb_natts, PrecipTotRate_MS_cmb_natts;

   int lat_NS_cmb_id, lon_NS_cmb_id;
   float *lat_NS_cmb_buf, *lon_NS_cmb_buf;
   nc_type lat_NS_cmb_type;
   int lat_NS_cmb_ndims;
   int lat_NS_cmb_dimids[NC_MAX_VAR_DIMS];
   int lat_NS_cmb_natts;

   int lat_MS_cmb_id, lon_MS_cmb_id;
   float *lat_MS_cmb_buf, *lon_MS_cmb_buf;
   nc_type lat_MS_cmb_type;
   int lat_MS_cmb_ndims;
   int lat_MS_cmb_dimids[NC_MAX_VAR_DIMS];
   int lat_MS_cmb_natts;

   int surfPrecipTotRate_NS_cmb_id, surfPrecipTotRate_MS_cmb_id;
   float *surfPrecipTotRate_NS_cmb_buf, *surfPrecipTotRate_MS_cmb_buf;

   int dataQuality_NS_cmb_id, dataQuality_MS_cmb_id;
   signed char *dataQuality_NS_cmb_buf, *dataQuality_MS_cmb_buf;

   float precip_NS_cmb, precip_MS_cmb;

   /*---------------------------------------------------------------------------------*/

   int nlines= 50;		/* Process this many lines on either side of the center scan line */
   int save_top_ranked= 1;	/* Set < zero to NOT save top-ranked database profiles and TB in output file and make files smaller */
   int save_dpr_profiles= 1;	/* Set < zero to NOT save DPR profiles in output file and make files smaller */
   int save_model= 1;		/* Set < zero to NOT save MERRA2 variables in output file and make files smaller */
   int save_pc= 1;		/* Set < zero to NOT save emis PC variables in output file and make files smaller */

   /*---------------------------------------------------------------------------------*/

   double km, km_min, km_max, secs_min;
   int qual;
   int iyyyy, imm, idd, ihh, imn, iss;
   float pcent, pcent1, tb, tb_gmi[20];
   float glat1, glon1, glat2, glon2, slat, slon, hrs, clat, clon;
   float clat1, clon1, clat2, clon2;
   int ilat, ilon, idx, idx1, idx2, idx3, s1s2, sfc_class, idx_db, idx_db1, idx_db2, ndb_expand, idx_db_expand, sign, ndb, idx_out, idx3_out;
   float pcent_db;

   /*--- output file ---*/
   int ncid2;
   char outfile[512];
   char outfile2[512];
   char outfilename[256];
   int dimids6[NC_MAX_VAR_DIMS];
   int dimids5[NC_MAX_VAR_DIMS];
   int dimids4[NC_MAX_VAR_DIMS];
   int dimids3[NC_MAX_VAR_DIMS];
   int dimids2[NC_MAX_VAR_DIMS];
   int scan_dimid, pix_dimid, tb_dimid, em_dimid, lev_dimid, lev_precip_dimid;
   char ctmp[256];
   size_t start3[3], count3[3];
   size_t start2[2], count2[2];
   int retval;
   char units[256];
   int pixel_index_varid;
   int lat_varid, lon_varid, inc_varid;
   int secs_varid, tb_varid, emis_varid, pc_emis_varid;
   int precip1_NS_varid, precip1_prob_NS_varid;
   int precip2_NS_varid, precip2_prob_NS_varid;
   int precip1_MS_varid, precip1_prob_MS_varid;
   int precip2_MS_varid, precip2_prob_MS_varid;
   int precip1_prof_NS_varid;
   int precip2_prof_NS_varid;
   int precip1_prof_MS_varid;
   int precip2_prof_MS_varid;
   int precip_GPROF_varid, precip_prob_GPROF_varid;
   int frozen_precip_GPROF_varid;
   int ts_NS_varid, t2m_NS_varid, tqv_NS_varid, h273_NS_varid, ku_ztop_NS_varid;
   int ts_MS_varid, t2m_MS_varid, tqv_MS_varid, h273_MS_varid, ku_ztop_MS_varid, ka_ztop_MS_varid;
   int ts_MERRA2_varid, t2m_MERRA2_varid, tqv_MERRA2_varid;
   int t2m_wet_MERRA2_varid, t2m_dew_MERRA2_varid;
   int z_ku_NS_varid;
   int z_ku_MS_varid, z_ka_MS_varid;
   int tb1_NS_varid;
   int tb1_MS_varid;
   int h273_MERRA2_varid;
   int quality_NS_varid;
   int quality_MS_varid;

   int ishift, jshift;
   int ishift_min=1000, jshift_min=1000;
   int ishift_max=-1000, jshift_max=-1000;

  FILE *fpc;
  float ave_pc[NEM], std_pc[NEM];
  float ave_tb[NCHAN], std_tb[NCHAN];
  float ave_emis[NEM], std_emis[NEM];
  float b_all[NREG+1][NEM];  /* regression coeffs for all NTBREG */
  float tmp5[NREG+1], psum;
  float pc_emis[NEM], emis[NEM], u[NEM+1][NEM+1];
  char emis_name[32], tb_name[32];
  int k2, kt, nt, kc, id;

   float *emis_buf;
   float *pc_emis_buf;

   float *precip1_NS_buf;
   float *precip2_NS_buf;
   short *precip1_prof_NS_buf;
   short *precip2_prof_NS_buf;
   float *precip1_prob_NS_buf;
   float *precip2_prob_NS_buf;
   int *epc_qual_NS_buf;
   float *precip1_MS_buf;
   float *precip2_MS_buf;
   short *precip1_prof_MS_buf;
   short *precip2_prof_MS_buf;
   float *precip1_prob_MS_buf;
   float *precip2_prob_MS_buf;
   int *epc_qual_MS_buf;

   FILE *fdb2, *fdb_nrain;
   int ndb_files= pow(NPCHIST, NEM_USE);
   int p, pow2[NEM_USE], ndb_files_used, bytes;
   float pcmin, pcmax, pc_range[NEM][NPCHIST+1], pc_lo, pc_hi;
   int idxn[NEM_USE], found, kuka, idb;
   int nem_compare, nem_compare_varid, ntb_compare_varid;
   int ntb_compare;

   int ix, iht;
   float x, rmsd_pc, rmsd_tb, range=0.0, wt, wt0, wt2, std;
   float ht, ku_ztop, ka_ztop;
   int maxrec, nrain_db_total[6][2], nrain_db[6][2], nrain, nrec, nrec1, nrec2, n1, n2, class[NCLASS], iclass, nmodel;
   float z_ku[NLEV_DPR], z_ka[NLEV_DPR], zlog, zlog1, zlog2, zlog3, zlin, zmax;
   float prob_precip1, prob_precip2, precip1, precip2;
   int nz_ku[NLEV_DPR], nz_ka[NLEV_DPR], hist_precip[NPRECIP], ku_clutter, ka_clutter, iz;
   int nrmsd_tb, nrmsd_pc, iprecip;
   int nfound, nfound_all=0;
   float rmsd_pc_min, rmsd_pc_max, rmsd_tb_min, rmsd_tb_max;
   short *z_ku_1_NS_buf;
   short *z_ku_1_MS_buf, *z_ka_1_MS_buf;
   float *tb_1_NS_buf, *tb_1_MS_buf;
   short *h273_buf;

   float tlat= -999.0, tlon= -999.0;
   int mrms= 0;  /* 1= look for MRMS data near the target time  */
   int n_mrms=0;

   int imax, i1max, i2max, nmax;
   float fmax;
  
   FILE *fccfads, *fccfads2;		/* CCFADS at the end */
   int sum1, sum2, ccfads_ku[NHT][NZ][2], ccfads_ka[NHT][NZ][2]; 
   float cdf_ccfads[10], ccfads1[NZ], ccfads2[NZ]; 
   int cdf_total, cdf_sum;

   float ts, t2m, tqv;
   float ts_sum, t2m_sum, tqv_sum, h273_sum, precip1_sum, precip2_sum;
   float wt_ts_sum, wt_t2m_sum, wt_tqv_sum, wt_h273_sum, wt_precip_sum;
   int n_precip1_sum, n_precip2_sum;
   int n_precip_all_sum;
   float ku_ztop_sum, ka_ztop_sum;
   float wt_ku_ztop_sum, wt_ka_ztop_sum;
   float precip1_prof_sum[NLEV_PRECIP], precip2_prof_sum[NLEV_PRECIP];
   int n_precip1_prof_sum[NLEV_PRECIP], n_precip2_prof_sum[NLEV_PRECIP];
   short precip1_prof[NLEV_PRECIP], precip2_prof[NLEV_PRECIP];

   int typePrecip[3], flagShallowRain[5];
   float type_precip[3], shallow_rain[5];

   int *ibuf;
   short *sbuf;
   char *bbuf;
   float *fbuf, *fbuf3, *fbuf4;
   double *dbuf;
   int nl_out, nc_out;

   FILE *fmrms;
   char mrms_name[512];


   /*---- variables for MERRA2 ----*/

   int have_model= -1;
   char model_file[256];
   char mfile[256], m1dir[256], m2dir[256], m3dir[256];
   float lat, lon, temp, t, qv, h, h1, ps, hs, hs1, p_lowest, h_lowest, slp, ps2;

   float *p_buf;
   float *lat_buf, *lon_buf;
   float *h_buf, *t_buf, *qv_buf;
   float *ps_buf, *hs_buf;
   float *ts_buf, *t2m_buf, *tqv_buf;
   float *t2m_wet_buf, *t2m_dew_buf;
   float *ts_epc_NS_buf, *t2m_epc_NS_buf, *tqv_epc_NS_buf;
   short *h273_epc_NS_buf;
   short *ku_ztop_epc_NS_buf;
   float *ts_epc_MS_buf, *t2m_epc_MS_buf, *tqv_epc_MS_buf;
   short *h273_epc_MS_buf;
   short *ku_ztop_epc_MS_buf, *ka_ztop_epc_MS_buf;
   double *time_buf;

   short badshort= -32768;
   float badfloat=-1000.0, weight, total_weight, lat_ratio, lon_ratio, time_ratio;
   int ntotal_weight, ilev_lowest;

   size_t nlat, nlon, nlevel;

   float lat1, lon1, lat2, lon2;
   int icube, iday, ilat1, ilon1, ilat2, ilon2, itime1, itime2;
   int lat_id, lon_id, time_id, height_id;
   int p_id, h_id, t_id, qv_id, lev_id, ps_id, hs_id, ts_id, t2m_id, slp_id, tqv_id;
   int t2m_wet_id, t2m_dew_id;

   nc_type t_type, ts_type, tqv_type;
   int t_ndims, ts_ndims, tqv_ndims;
   int t_dimids[NC_MAX_VAR_DIMS], ts_dimids[NC_MAX_VAR_DIMS], tqv_dimids[NC_MAX_VAR_DIMS];
   int t_natts, ts_natts, tqv_natts;

   float ts_MERRA2, t2m_MERRA2, tqv_MERRA2;
   short h273_MERRA2;

   int plot_files= 0, elev, nku20, nku25;
   float e0, e1, frac, frac0, frac1, frac5, h273;

   int ibin_sfc, ibin_clutter_free;
   int i1_start=-1, i1_end=-1;
   int imin_NS, jmin_NS;
   int imin_MS, jmin_MS;
   float precip_NS, precip_MS;
   float nslat1, nslon1, nslon1_360, glon1_360, lon_diff, z, diff;
   float nslat2, nslon2;
   int nclutter_ku[NLEV_DPR];
   int nclutter_ka[NLEV_DPR];

   short *z_ku_DPR_buf, *z_ka_DPR_buf;
   int z_ku_DPR_varid, z_ka_DPR_varid;
   short *elev_DPR_buf;
   int elev_DPR_varid;
   float *lat_NS_DPR_buf, *lon_NS_DPR_buf;
   int lat_NS_DPR_varid, lon_NS_DPR_varid;
   float *precip_NS_DPR_buf, *precip_MS_DPR_buf;
   int precip_NS_DPR_varid, precip_MS_DPR_varid;
   float *precip_NS_CMB_buf, *precip_MS_CMB_buf;
   int precip_NS_CMB_varid, precip_MS_CMB_varid;
   short *ztop_ku_DPR_buf, *ztop_ka_DPR_buf;
   int ztop_ku_DPR_varid, ztop_ka_DPR_varid;

   int precip_prof_NS_CMB_varid, precip_prof_MS_CMB_varid;
   short *precip_prof_NS_CMB_buf, *precip_prof_MS_CMB_buf;

   /* estimated - for output file */
   int type_precip_NS_varid, shallow_rain_NS_varid;
   short *type_precip_NS_buf, *shallow_rain_NS_buf;
   int type_precip_MS_varid, shallow_rain_MS_varid;
   short *type_precip_MS_buf, *shallow_rain_MS_buf;

   /* actual DPR - for output file */
   int type_precip_DPR_varid, shallow_rain_DPR_varid;
   short *type_precip_DPR_buf, *shallow_rain_DPR_buf;

   double now;
   char processing_date_start[64], processing_date_end[64];

   int inside_NS=-1, inside_MS=-1;
   int *dpr_NS_i, *dpr_NS_j;

   /*--------------------------------------------*/

   now= time(NULL);
   t0= now;
   tm= *gmtime(&t0);
   sprintf(processing_date_start,"%4d/%02d/%02d %02d:%02d:%02d", 1900+tm.tm_year, 1+tm.tm_mon, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
   printf("Processing started at %s\n", processing_date_start);

   sprintf(model_file,"%s", "none");
   sprintf(gprof_file,"%s", "none");
   sprintf(dpr_file,"%s", "none");
   sprintf(cmb_file,"%s", "none");

   tlat= atof(argv[2]);
   tlon= atof(argv[3]);
   isat_start= atoi(argv[4]);
   isat_end= atoi(argv[5]);

   if ( isat_start >= 0 && isat_end > isat_start )
     printf("Processing orbit between scan lines %d %d\n", isat_start, isat_end);
   else if ( fabs(tlat) < 90 && fabs(tlon) < 180 )
     printf("Processing orbit for %d scan lines nearest to %f %f\n", nlines, tlat, tlon);
   else if ( fabs(tlat) > 90 || fabs(tlon) > 180 || isat_start < 0 || isat_end < 0 )
     printf("Processing entire orbit\n");
   else {
     printf("Unrecognized options %s %s %s %s\n", argv[2], argv[3], argv[4], argv[5]);
     exit(-1);
   }

   sprintf(outdir,"%s", argv[6]);
   /*sprintf(outdir,"%s", "/Users/jturk/code/GPM/GPM_DB_V5/EPC_DB_ENTRIES_V5_MRMS");*/
   printf("outdir= %s\n", outdir);

   if ( argc >= 8 ) {
     sprintf(model_file,"%s", argv[7]);
     printf("MERRA2= %s\n", model_file);
   }
   if ( argc >= 9 ) {
     sprintf(gprof_file,"%s", argv[8]);
     printf("GPROF= %s\n", gprof_file);
   }
   if ( argc >= 10 ) {
     sprintf(dpr_file,"%s", argv[9]);
     printf("DPR= %s\n", dpr_file);
   }
   if ( argc == 11 ) {
     sprintf(cmb_file,"%s", argv[10]);
     printf("CMB= %s\n", cmb_file);
   }

   bytes= sizeof(strm_record);
   printf("record length= %d\n", bytes);

   tmp1= (int *) malloc(DB_MAXREC*sizeof(int));
   tmp2= (int *) malloc(DB_MAXREC*sizeof(int));
   tmp3= (float *) malloc(DB_MAXREC*sizeof(float));
   tmp7= (float *) malloc(DB_MAXREC*sizeof(float));

   for (k= 0; k< NHT; k++) {
     for (k1= 0; k1< NZ; k1++) {
       ccfads_ku[k][k1][0]= 0;
       ccfads_ka[k][k1][0]= 0;
       ccfads_ku[k][k1][1]= 0;
       ccfads_ka[k][k1][1]= 0;
     }
   }

   for (k= 0; k< NEM_USE; k++)
     pow2[k]= pow(NPCHIST,k);


   /*--- Read in EPC database index values ---*/


   sprintf(svar,"%s/%s", epc_coef_dir, "PC_MIN_MAX_29.txt");
   if ((fpc= fopen(svar,"r")) == 0) {
     printf("Unable to open PC min-max file %s\n", svar);
     exit(-1);
   }
   printf("Opened PC min-max file %s\n", svar);
   for (i=0; i< NEM; i++) {
     fscanf(fpc,"%d ", &k);
     printf("%2d ", k);
     for (j=0; j< NPCHIST+1; j++) {
       fscanf(fpc,"%f ", &pc_range[i][j]) ;
       if ( j == 0 ) pc_range[i][j]-= 1.0E6;
       if ( j == NPCHIST ) pc_range[i][j]+= 1.0E6;
       printf("%12.6f ", pc_range[i][j]);
     }
     fscanf(fpc,"\n");
     printf("\n");
   }
   fclose(fpc);


   /*--- Read EPC coefficients file ---*/

   sprintf(svar,"%s/%s", epc_coef_dir, "coef_pc.txt");
   if ((fpc= fopen(svar,"r")) == 0) {
     printf("Unable to open %s\n", svar);
     exit(-1);
   }
   printf("Opened %s\n", svar);
   for (i=0; i<= NREG; i++) {
     fscanf(fpc,"%d ", &k);
     for (j=0; j< NEM; j++) {
       fscanf(fpc,"%f ", &b_all[i][j]);
       printf("%f ", b_all[i][j]);
     }
     fscanf(fpc,"\n");
     printf("\n");
   }
   fclose(fpc);

   /*--- Read eigenvalues wmatrix file ---*/

   sprintf(svar,"%s/%s", epc_coef_dir, "wmatrix.txt");
   if ((fpc= fopen(svar,"r")) == 0) {
     printf("Unable to open %s\n", svar);
     exit(-1);
   }
   printf("Opened %s\n", svar);
   for (i = 0; i < NEM; i++) {
     for (j = 0; j < NEM; j++) fscanf(fpc,"%f ",&u[i+1][j+1]);
     fscanf(fpc,"\n");
     for (j = 0; j < NEM; j++) printf("%f ",u[i+1][j+1]);
     printf("\n");
   }
   fclose(fpc);

   /*--- Read ave emissivity file ---*/

   sprintf(svar,"%s/%s", epc_coef_dir, "ave_emis.txt");
   if ((fpc= fopen(svar,"r")) == 0) {
     printf("Unable to open %s\n", svar);
     exit(-1);
   }
   printf("Opened %s\n", svar);
   for (i = 0; i < NEM; i++) {
     fscanf(fpc,"%d %s %f %f\n", &j, &emis_name, &ave_emis[i], &std_emis[i]);
     printf("%3d %16s %12.6f %12.6f\n", i, emis_name, ave_emis[i], std_emis[i]);
   }
   fclose(fpc);

   /*--- Read ave EPC file ---*/

   sprintf(svar,"%s/%s", epc_coef_dir, "ave_pc.txt");
   if ((fpc= fopen(svar,"r")) == 0) {
     printf("Unable to open %s\n", svar);
     exit(-1);
   }
   printf("Opened %s\n", svar);
   for (i = 0; i < NEM; i++) {
     fscanf(fpc,"%d %f %f\n", &j, &ave_pc[i], &std_pc[i]);
     printf("%2d %12.6f %12.6f\n", i, ave_pc[i], std_pc[i]);
   }
   fclose(fpc);

   /*--- Read ave TB file ---*/

   sprintf(svar,"%s/%s", epc_coef_dir, "ave_tb.txt");
   if ((fpc= fopen(svar,"r")) == 0) {
     printf("Unable to open %s\n", svar);
     exit(-1);
   }
   printf("Opened %s\n", svar);
   for (i = 0; i < NCHAN; i++) {
     fscanf(fpc,"%d %s %f %f\n", &j, &tb_name, &ave_tb[i], &std_tb[i]);
     printf("%3d %16s %12.6f %12.6f\n", i, tb_name, ave_tb[i], std_tb[i]);
   }
   fclose(fpc);


  /*---------------------------------------------------*/
  /* Read external 2-minute elevation database */

  int nx_id, ny_id;
  int x_id, y_id, z_id, gelev, gelev1;
  float *elat_buf, *elon_buf, *elev_buf;
  size_t nx_elev, ny_elev;
  char elevfile[64];

  sprintf(elevfile, "%s", "ETOPO2v2g_f4.nc");

  status = nc_open(elevfile, NC_NETCDF4, &ncid);
  if ( status != 0 ) {
    printf("File not found %s\n", elevfile);
    exit(1);
  }

  status = nc_inq_dimid(ncid, "x", &nx_id);
  status = nc_inq_dimlen(ncid, nx_id, &nx_elev);
  status = nc_inq_dimid(ncid, "y", &ny_id);
  status = nc_inq_dimlen(ncid, ny_id, &ny_elev);

  status= nc_inq_varid(ncid, "x", &x_id);
  status= nc_inq_varid(ncid, "y", &y_id);
  status= nc_inq_varid(ncid, "z", &z_id);

  elon_buf= (float *) malloc(nx_elev*sizeof(float));
  elat_buf= (float *) malloc(ny_elev*sizeof(float));
  elev_buf= (float *) malloc(nx_elev*ny_elev*sizeof(float));

  status = nc_get_var_float(ncid, x_id, elon_buf);
  status = nc_get_var_float(ncid, y_id, elat_buf);
  status = nc_get_var_float(ncid, z_id, elev_buf);

  status= nc_close(ncid);

  /*for (i= 0; i< ny_elev; i++) {
    for (j= 0; j< nx_elev; j++) {
      elev= elev_buf[i*nx_elev+j];
      if ( elev > 0 ) printf("%4d %4d %8.3f %8.3f %d\n", i, j, elat_buf[i], elon_buf[j], elev);
    }
  }*/

  /*---------------------------------------------------*/



   /*--------------------------Read 1C.NPP.ATMS -------------------------------------------*/

   status = nc_open(tbfile, NC_NETCDF4, &ncid);
   if ( status != 0 ) {
     printf("File not found %s\n", tbfile);
     exit(1);
   }
   printf("Opened %s\n", tbfile);

   if ((status= nc_inq_ncid(ncid, "S1", &grpid1))) ERR(status);
   if ( strcmp(sensor, "GMI") == 0 || strcmp(sensor, "ATMS") == 0 || strcmp(sensor, "SSMIS") == 0 || strcmp(sensor, "AMSR2") == 0 ) {
     if ((status= nc_inq_ncid(ncid, "S2", &grpid2))) ERR(status);
     if ( strcmp(sensor, "ATMS") == 0 || strcmp(sensor, "SSMIS") == 0 || strcmp(sensor, "AMSR2") == 0 ) {
       if ((status= nc_inq_ncid(ncid, "S3", &grpid3))) ERR(status);
       if ((status= nc_inq_ncid(ncid, "S4", &grpid4))) ERR(status);
       if ( strcmp(sensor, "AMSR2") == 0 ) {
         if ((status= nc_inq_ncid(ncid, "S5", &grpid5))) ERR(status);
         if ((status= nc_inq_ncid(ncid, "S6", &grpid6))) ERR(status);
       }
     }
   }

   if ((status= nc_inq_ncid(grpid1, "ScanTime", &grpid1_scantime))) ERR(status);
   if ((status= nc_inq_ncid(grpid1, "SCstatus", &grpid1_scanstatus))) ERR(status);

   /*--- S1 Tc lat lon incidence quality date --*/

   if ((status= nc_inq_varid(grpid1, "Tc", &Tc_S1_id))) ERR(status);
   if ((status = nc_inq_var(grpid1, Tc_S1_id, 0, &Tc_S1_type, &Tc_S1_ndims, Tc_S1_dimids, &Tc_S1_natts))) ERR(status);
   status = nc_inq_dimlen(grpid1, Tc_S1_dimids[0], &nl_S1);
   status = nc_inq_dimlen(grpid1, Tc_S1_dimids[1], &ns_S1);
   status = nc_inq_dimlen(grpid1, Tc_S1_dimids[2], &nc_S1);
   printf("Tc S1 dims= %d %d %d\n", (int) nl_S1, (int) ns_S1, (int) nc_S1);
   Tc_S1_buf= (float *) malloc(nl_S1*ns_S1*nc_S1*sizeof(float));
   if ((status = nc_get_var_float(grpid1, Tc_S1_id, Tc_S1_buf))) ERR(status);

   if ((status= nc_inq_varid(grpid1, "Latitude", &lat_S1_id))) ERR(status);
   if ((status= nc_inq_varid(grpid1, "Longitude", &lon_S1_id))) ERR(status);
   lat_S1_buf= (float *) malloc(nl_S1*ns_S1*sizeof(float));
   lon_S1_buf= (float *) malloc(nl_S1*ns_S1*sizeof(float));
   if ((status = nc_get_var_float(grpid1, lat_S1_id, lat_S1_buf))) ERR(status);
   if ((status = nc_get_var_float(grpid1, lon_S1_id, lon_S1_buf))) ERR(status);

   if ((status= nc_inq_varid(grpid1, "Quality", &qual_S1_id))) ERR(status);
   qual_S1_buf= (char *) malloc(nl_S1*ns_S1*sizeof(char));
   if ((status = nc_get_var_schar(grpid1, qual_S1_id, qual_S1_buf))) ERR(status);

   if ((status= nc_inq_varid(grpid1, "incidenceAngle", &incidenceAngle_S1_id))) ERR(status);
   incidenceAngle_S1_buf= (float *) malloc(nl_S1*ns_S1*sizeof(float));
   if ((status = nc_get_var_float(grpid1, incidenceAngle_S1_id, incidenceAngle_S1_buf))) ERR(status);

   /*--- S1 date ---*/

   status= nc_inq_varid(grpid1_scantime, "Year", &yyyy_S1_id);
   if ( status != 0 )  ERR(status);
   yyyy_S1_buf= (short *) malloc(nl_S1*sizeof(short));
   status = nc_get_var_short(grpid1_scantime, yyyy_S1_id, yyyy_S1_buf);
   if ( status != 0 ) ERR(status);

   status= nc_inq_varid(grpid1_scantime, "Month", &mm_S1_id);
   if ( status != 0 )  ERR(status);
   mm_S1_buf= (short *) malloc(nl_S1*sizeof(short));
   status = nc_get_var_short(grpid1_scantime, mm_S1_id, mm_S1_buf);
   if ( status != 0 ) ERR(status);

   status= nc_inq_varid(grpid1_scantime, "DayOfMonth", &dd_S1_id);
   if ( status != 0 )  ERR(status);
   dd_S1_buf= (short *) malloc(nl_S1*sizeof(short));
   status = nc_get_var_short(grpid1_scantime, dd_S1_id, dd_S1_buf);
   if ( status != 0 ) ERR(status);

   status= nc_inq_varid(grpid1_scantime, "Hour", &hh_S1_id);
   if ( status != 0 )  ERR(status);
   hh_S1_buf= (short *) malloc(nl_S1*sizeof(short));
   status = nc_get_var_short(grpid1_scantime, hh_S1_id, hh_S1_buf);
   if ( status != 0 ) ERR(status);

   status= nc_inq_varid(grpid1_scantime, "Minute", &mn_S1_id);
   if ( status != 0 )  ERR(status);
   mn_S1_buf= (short *) malloc(nl_S1*sizeof(short));
   status = nc_get_var_short(grpid1_scantime, mn_S1_id, mn_S1_buf);
   if ( status != 0 ) ERR(status);

   status= nc_inq_varid(grpid1_scantime, "Second", &ss_S1_id);
   if ( status != 0 )  ERR(status);
   ss_S1_buf= (short *) malloc(nl_S1*sizeof(short));
   status = nc_get_var_short(grpid1_scantime, ss_S1_id, ss_S1_buf);
   if ( status != 0 ) ERR(status);

   /*--- S1 spacecraft position ---*/

   status= nc_inq_varid(grpid1_scanstatus, "SClatitude", &sclat_S1_id);
   if ( status != 0 )  ERR(status);
   sclat_S1_buf= (float *) malloc(nl_S1*sizeof(float));
   status = nc_get_var_float(grpid1_scanstatus, sclat_S1_id, sclat_S1_buf);
   if ( status != 0 ) ERR(status);

   status= nc_inq_varid(grpid1_scanstatus, "SClongitude", &sclon_S1_id);
   if ( status != 0 )  ERR(status);
   sclon_S1_buf= (float *) malloc(nl_S1*sizeof(float));
   status = nc_get_var_float(grpid1_scanstatus, sclon_S1_id, sclon_S1_buf);
   if ( status != 0 ) ERR(status);

   status= nc_inq_varid(grpid1_scanstatus, "FractionalGranuleNumber", &FractionalGranuleNumber_S1_id);
   if ( status != 0 )  ERR(status);
   FractionalGranuleNumber_S1_buf= (double *) malloc(nl_S1*sizeof(double));
   status = nc_get_var_double(grpid1_scanstatus, FractionalGranuleNumber_S1_id, FractionalGranuleNumber_S1_buf);
   if ( status != 0 ) ERR(status);

   status= nc_inq_varid(grpid1_scanstatus, "SCorientation", &SCorientation_S1_id);
   if ( status != 0 )  ERR(status);
   SCorientation_S1_buf= (short *) malloc(nl_S1*sizeof(short));
   status = nc_get_var_short(grpid1_scanstatus, SCorientation_S1_id, SCorientation_S1_buf);
   if ( status != 0 ) ERR(status);

   if ( strcmp(sensor, "GMI") == 0 || strcmp(sensor, "ATMS") == 0 || strcmp(sensor, "SSMIS") == 0 || strcmp(sensor, "AMSR2") == 0 ) {
     if ((status= nc_inq_ncid(ncid, "S2", &grpid2))) ERR(status);

     /*--- S2 Tc --*/
     if ((status= nc_inq_varid(grpid2, "Tc", &Tc_S2_id))) ERR(status);
     if ((status = nc_inq_var(grpid2, Tc_S2_id, 0, &Tc_S2_type, &Tc_S2_ndims, Tc_S2_dimids, &Tc_S2_natts))) ERR(status);
     status = nc_inq_dimlen(grpid2, Tc_S2_dimids[0], &nl_S2);
     status = nc_inq_dimlen(grpid2, Tc_S2_dimids[1], &ns_S2);
     status = nc_inq_dimlen(grpid2, Tc_S2_dimids[2], &nc_S2);
     printf("Tc S2 dims= %d %d %d\n", (int) nl_S2, (int) ns_S2, (int) nc_S2);
     Tc_S2_buf= (float *) malloc(nl_S2*ns_S2*nc_S2*sizeof(float));
     if ((status = nc_get_var_float(grpid2, Tc_S2_id, Tc_S2_buf))) ERR(status);

     if ((status= nc_inq_varid(grpid2, "Latitude", &lat_S2_id))) ERR(status);
     if ((status= nc_inq_varid(grpid2, "Longitude", &lon_S2_id))) ERR(status);
     lat_S2_buf= (float *) malloc(nl_S2*ns_S2*sizeof(float));
     lon_S2_buf= (float *) malloc(nl_S2*ns_S2*sizeof(float));
     if ((status = nc_get_var_float(grpid2, lat_S2_id, lat_S2_buf))) ERR(status);
     if ((status = nc_get_var_float(grpid2, lon_S2_id, lon_S2_buf))) ERR(status);

     if ((status= nc_inq_varid(grpid2, "Quality", &qual_S2_id))) ERR(status);
     qual_S2_buf= (char *) malloc(nl_S2*ns_S2*sizeof(char));
     if ((status = nc_get_var_schar(grpid2, qual_S2_id, qual_S2_buf))) ERR(status);

     if ((status= nc_inq_varid(grpid2, "incidenceAngle", &incidenceAngle_S2_id))) ERR(status);
     incidenceAngle_S2_buf= (float *) malloc(nl_S2*ns_S2*sizeof(float));
     if ((status = nc_get_var_float(grpid2, incidenceAngle_S2_id, incidenceAngle_S2_buf))) ERR(status);

     if ( strcmp(sensor, "ATMS") == 0 || strcmp(sensor, "SSMIS") == 0 || strcmp(sensor, "AMSR2") == 0 ) {
       if ((status= nc_inq_ncid(ncid, "S3", &grpid3))) ERR(status);
       if ((status= nc_inq_ncid(ncid, "S4", &grpid4))) ERR(status);

       /*--- S3 Tc --*/
       if ((status= nc_inq_varid(grpid3, "Tc", &Tc_S3_id))) ERR(status);
       if ((status = nc_inq_var(grpid3, Tc_S3_id, 0, &Tc_S3_type, &Tc_S3_ndims, Tc_S3_dimids, &Tc_S3_natts))) ERR(status);
       status = nc_inq_dimlen(grpid3, Tc_S3_dimids[0], &nl_S3);
       status = nc_inq_dimlen(grpid3, Tc_S3_dimids[1], &ns_S3);
       status = nc_inq_dimlen(grpid3, Tc_S3_dimids[2], &nc_S3);
       printf("Tc S3 dims= %d %d %d\n", (int) nl_S3, (int) ns_S3, (int) nc_S3);
       Tc_S3_buf= (float *) malloc(nl_S3*ns_S3*nc_S3*sizeof(float));
       if ((status = nc_get_var_float(grpid3, Tc_S3_id, Tc_S3_buf))) ERR(status);

       if ((status= nc_inq_varid(grpid3, "Latitude", &lat_S3_id))) ERR(status);
       if ((status= nc_inq_varid(grpid3, "Longitude", &lon_S3_id))) ERR(status);
       lat_S3_buf= (float *) malloc(nl_S3*ns_S3*sizeof(float));
       lon_S3_buf= (float *) malloc(nl_S3*ns_S3*sizeof(float));
       if ((status = nc_get_var_float(grpid3, lat_S3_id, lat_S3_buf))) ERR(status);
       if ((status = nc_get_var_float(grpid3, lon_S3_id, lon_S3_buf))) ERR(status);

       if ((status= nc_inq_varid(grpid3, "Quality", &qual_S3_id))) ERR(status);
       qual_S3_buf= (char *) malloc(nl_S3*ns_S3*sizeof(char));
       if ((status = nc_get_var_schar(grpid3, qual_S3_id, qual_S3_buf))) ERR(status);

       /*--- S4 Tc --*/
       if ((status= nc_inq_varid(grpid4, "Tc", &Tc_S4_id))) ERR(status);
       if ((status = nc_inq_var(grpid4, Tc_S4_id, 0, &Tc_S4_type, &Tc_S4_ndims, Tc_S4_dimids, &Tc_S4_natts))) ERR(status);
       status = nc_inq_dimlen(grpid4, Tc_S4_dimids[0], &nl_S4);
       status = nc_inq_dimlen(grpid4, Tc_S4_dimids[1], &ns_S4);
       status = nc_inq_dimlen(grpid4, Tc_S4_dimids[2], &nc_S4);
       printf("Tc S4 dims= %d %d %d\n", (int) nl_S4, (int) ns_S4, (int) nc_S4);
       Tc_S4_buf= (float *) malloc(nl_S4*ns_S4*nc_S4*sizeof(float));
       if ((status = nc_get_var_float(grpid4, Tc_S4_id, Tc_S4_buf))) ERR(status);

       if ((status= nc_inq_varid(grpid4, "Latitude", &lat_S4_id))) ERR(status);
       if ((status= nc_inq_varid(grpid4, "Longitude", &lon_S4_id))) ERR(status);
       lat_S4_buf= (float *) malloc(nl_S4*ns_S4*sizeof(float));
       lon_S4_buf= (float *) malloc(nl_S4*ns_S4*sizeof(float));
       if ((status = nc_get_var_float(grpid4, lat_S4_id, lat_S4_buf))) ERR(status);
       if ((status = nc_get_var_float(grpid4, lon_S4_id, lon_S4_buf))) ERR(status);

       if ((status= nc_inq_varid(grpid4, "Quality", &qual_S4_id))) ERR(status);
       qual_S4_buf= (char *) malloc(nl_S4*ns_S4*sizeof(char));
       if ((status = nc_get_var_schar(grpid4, qual_S4_id, qual_S4_buf))) ERR(status);

       if ( strcmp(sensor, "AMSR2") == 0 ) {
         if ((status= nc_inq_ncid(ncid, "S5", &grpid5))) ERR(status);
         if ((status= nc_inq_ncid(ncid, "S6", &grpid6))) ERR(status);

         /*--- S5 Tc and coordinates --*/
         if ((status= nc_inq_varid(grpid5, "Tc", &Tc_S5_id))) ERR(status);
         if ((status = nc_inq_var(grpid5, Tc_S5_id, 0, &Tc_S5_type, &Tc_S5_ndims, Tc_S5_dimids, &Tc_S5_natts))) ERR(status);
         status = nc_inq_dimlen(grpid5, Tc_S5_dimids[0], &nl_S5);
         status = nc_inq_dimlen(grpid5, Tc_S5_dimids[1], &ns_S5);
         status = nc_inq_dimlen(grpid5, Tc_S5_dimids[2], &nc_S5);
         printf("Tc S5 dims= %d %d %d\n", (int) nl_S5, (int) ns_S5, (int) nc_S5);
         Tc_S5_buf= (float *) malloc(nl_S5*ns_S5*nc_S5*sizeof(float));
         if ((status = nc_get_var_float(grpid5, Tc_S5_id, Tc_S5_buf))) ERR(status);

         if ((status= nc_inq_varid(grpid5, "Latitude", &lat_S5_id))) ERR(status);
         if ((status= nc_inq_varid(grpid5, "Longitude", &lon_S5_id))) ERR(status);
         lat_S5_buf= (float *) malloc(nl_S5*ns_S5*sizeof(float));
         lon_S5_buf= (float *) malloc(nl_S5*ns_S5*sizeof(float));
         if ((status = nc_get_var_float(grpid5, lat_S5_id, lat_S5_buf))) ERR(status);
         if ((status = nc_get_var_float(grpid5, lon_S5_id, lon_S5_buf))) ERR(status);

         if ((status= nc_inq_varid(grpid5, "Quality", &qual_S5_id))) ERR(status);
         qual_S5_buf= (char *) malloc(nl_S5*ns_S5*sizeof(char));
         if ((status = nc_get_var_schar(grpid5, qual_S5_id, qual_S5_buf))) ERR(status);

         /*--- S6 Tc and coordinates --*/
         if ((status= nc_inq_varid(grpid6, "Tc", &Tc_S6_id))) ERR(status);
         if ((status = nc_inq_var(grpid6, Tc_S6_id, 0, &Tc_S6_type, &Tc_S6_ndims, Tc_S6_dimids, &Tc_S6_natts))) ERR(status);
         status = nc_inq_dimlen(grpid6, Tc_S6_dimids[0], &nl_S6);
         status = nc_inq_dimlen(grpid6, Tc_S6_dimids[1], &ns_S6);
         status = nc_inq_dimlen(grpid6, Tc_S6_dimids[2], &nc_S6);
         printf("Tc S6 dims= %d %d %d\n", (int) nl_S6, (int) ns_S6, (int) nc_S6);
         Tc_S6_buf= (float *) malloc(nl_S6*ns_S6*nc_S6*sizeof(float));
         if ((status = nc_get_var_float(grpid6, Tc_S6_id, Tc_S6_buf))) ERR(status);

         if ((status= nc_inq_varid(grpid6, "Latitude", &lat_S6_id))) ERR(status);
         if ((status= nc_inq_varid(grpid6, "Longitude", &lon_S6_id))) ERR(status);
         lat_S6_buf= (float *) malloc(nl_S6*ns_S6*sizeof(float));
         lon_S6_buf= (float *) malloc(nl_S6*ns_S6*sizeof(float));
         if ((status = nc_get_var_float(grpid6, lat_S6_id, lat_S6_buf))) ERR(status);
         if ((status = nc_get_var_float(grpid6, lon_S6_id, lon_S6_buf))) ERR(status);

         if ((status= nc_inq_varid(grpid6, "Quality", &qual_S6_id))) ERR(status);
         qual_S6_buf= (char *) malloc(nl_S6*ns_S6*sizeof(char));
         if ((status = nc_get_var_schar(grpid6, qual_S6_id, qual_S6_buf))) ERR(status);
       }
     }
   }

   if ((status= nc_close(ncid))) ERR(status);
   printf("Closed %s\n", tbfile);

   secs_S1_buf= (double *) malloc(nl_S1*sizeof(double));

   if ( strcmp(sensor, "ATMS") == 0 || strcmp(sensor, "SSMIS") == 0 ) nc_out= nc_S1 + nc_S2 + nc_S3 + nc_S4;
   if ( strcmp(sensor, "AMSR2") == 0 ) nc_out= nc_S1 + nc_S2 + nc_S3 + nc_S4 + nc_S5;
   if ( strcmp(sensor, "GMI") == 0 ) nc_out= nc_S1 + nc_S2;
   if ( strcmp(sensor, "MHS") == 0 || strcmp(sensor, "SAPHIR") == 0 ) nc_out= nc_S1;

   Tc_buf= (float *) malloc(nl_S1*ns_S1*nc_out*sizeof(float));
   for (i=0; i < nl_S1*ns_S1*nc_out; i++) Tc_buf[k]= -999.99;

   if ( strcmp(sensor, "ATMS") == 0 ) {
     printf("%d %d %d %d \n",  (int) nl_S1, (int) ns_S1, (int) nc_S1, nc_out);
     tc_atms(Tc_S1_buf, Tc_S2_buf, Tc_S3_buf, Tc_S4_buf,
          (int) nl_S1, (int) ns_S1, (int) nc_S1, 
          (int) nl_S2, (int) ns_S2, (int) nc_S2, 
          (int) nl_S3, (int) ns_S3, (int) nc_S3, 
          (int) nl_S4, (int) ns_S4, (int) nc_S4, 
          lat_S1_buf, lon_S1_buf,
          qual_S1_buf, qual_S2_buf, qual_S3_buf, qual_S4_buf, 
          incidenceAngle_S1_buf,
          Tc_buf);
   }
   else if ( strcmp(sensor, "SSMIS") == 0 ) {
     printf("%d %d %d %d \n",  (int) nl_S1, (int) ns_S1, (int) nc_S1, nc_out);
     tc_ssmis(Tc_S1_buf, Tc_S2_buf, Tc_S3_buf, Tc_S4_buf,
          (int) nl_S1, (int) ns_S1, (int) nc_S1, 
          (int) nl_S2, (int) ns_S2, (int) nc_S2, 
          (int) nl_S3, (int) ns_S3, (int) nc_S3, 
          (int) nl_S4, (int) ns_S4, (int) nc_S4, 
          lat_S1_buf, lon_S1_buf, lat_S4_buf, lon_S4_buf,
          qual_S1_buf, qual_S2_buf, qual_S3_buf, qual_S4_buf, 
          incidenceAngle_S1_buf,
          Tc_buf);
   }
   else if ( strcmp(sensor, "AMSR2") == 0 ) {
     printf("%d %d %d %d \n",  (int) nl_S1, (int) ns_S1, (int) nc_S1, nc_out);
     tc_amsr2(Tc_S1_buf, Tc_S2_buf, Tc_S3_buf, Tc_S4_buf, Tc_S5_buf, Tc_S6_buf,
          (int) nl_S1, (int) ns_S1, (int) nc_S1,
          (int) nl_S2, (int) ns_S2, (int) nc_S2,
          (int) nl_S3, (int) ns_S3, (int) nc_S3,
          (int) nl_S4, (int) ns_S4, (int) nc_S4,
          (int) nl_S5, (int) ns_S5, (int) nc_S5,
          (int) nl_S6, (int) ns_S6, (int) nc_S6,
          lat_S1_buf, lon_S1_buf, lat_S5_buf, lon_S5_buf, lat_S6_buf, lon_S6_buf,
          qual_S1_buf, qual_S2_buf, qual_S3_buf, qual_S4_buf, qual_S5_buf, qual_S6_buf,
          incidenceAngle_S1_buf,
          Tc_buf);
   }
   else if ( strcmp(sensor, "GMI") == 0 ) {
     printf("%d %d %d %d \n",  (int) nl_S1, (int) ns_S1, (int) nc_S1, nc_out);
     tc_gmi(Tc_S1_buf, Tc_S2_buf, 
          (int) nl_S1, (int) ns_S1, (int) nc_S1,
          (int) nl_S2, (int) ns_S2, (int) nc_S2,
          lat_S1_buf, lon_S1_buf, lat_S2_buf, lon_S2_buf, 
          qual_S1_buf, qual_S2_buf, 
          incidenceAngle_S1_buf, incidenceAngle_S2_buf, 
          Tc_buf);
   }
   else if ( strcmp(sensor, "MHS") == 0 || strcmp(sensor, "SAPHIR") == 0 ) {
     printf("%d %d %d %d \n",  (int) nl_S1, (int) ns_S1, (int) nc_S1, nc_out);
     tc_mhs(Tc_S1_buf,
          (int) nl_S1, (int) ns_S1, (int) nc_S1,
          lat_S1_buf, lon_S1_buf, 
          qual_S1_buf, 
          incidenceAngle_S1_buf,
          Tc_buf);
   }
   else {
     printf("%s not supported\n", sensor);
     exit(1);
   }

   /*  -400 < TB < -20 implies that 1C quality flag is > 0  useable */
   /*  TB < -400 implies that 1C quality flag is < 0  not useable */
   for (i=0; i < nl_S1; i++) {
     for (j=0; j< ns_S1; j++) {
       for (k=0; k< nc_out; k++) {
         idx1= i*ns_S1*nc_out + j*nc_out + k;
         if ( Tc_buf[idx1] < -20 && Tc_buf[idx1] > -400 ) Tc_buf[idx1]= -1.0*Tc_buf[idx1];
       }
     }
   }

   irev= FractionalGranuleNumber_S1_buf[nl_S1/2];

   for (i=0; i < nl_S1; i++) {
     iyyyy= yyyy_S1_buf[i];
     imm= mm_S1_buf[i];
     idd= dd_S1_buf[i];
     ihh= hh_S1_buf[i];
     imn= mn_S1_buf[i];
     iss= ss_S1_buf[i];
     tm.tm_year= iyyyy-1900;
     tm.tm_mon= imm-1;
     tm.tm_mday= idd;
     tm.tm_hour= ihh;
     tm.tm_min= imn;
     tm.tm_sec= iss;
     t0= timegm(&tm);
     secs_S1_buf[i]= ( double ) t0;

     tm= *gmtime(&t0);
     sprintf(tdate,"%04d/%02d/%02d %02d:%02d:%02d",1900+tm.tm_year,1+tm.tm_mon,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);
     idx= i*ns_S1 + (ns_S1/2);
     lat1= lat_S1_buf[idx];
     lon1= lon_S1_buf[idx];
     printf("%s scan=%4d rev=%6d orient=%6d date=%s %8.3f %8.3f\n", sensor, i, irev, (int) SCorientation_S1_buf[i], tdate, lat1, lon1);
   }

   free(Tc_S1_buf);
   free(qual_S1_buf);
   if ( strcmp(sensor, "GMI") == 0 || strcmp(sensor, "ATMS") == 0 || strcmp(sensor, "SSMIS") == 0 || strcmp(sensor, "AMSR2") == 0 ) {
     free(Tc_S2_buf);
     free(qual_S2_buf);
     free(lat_S2_buf); free(lon_S2_buf);
     if ( strcmp(sensor, "ATMS") == 0 || strcmp(sensor, "SSMIS") == 0 || strcmp(sensor, "AMSR2") == 0 ) {
       free(Tc_S3_buf);
       free(qual_S3_buf);
       free(lat_S3_buf); free(lon_S3_buf);
       free(Tc_S4_buf);
       free(qual_S4_buf);
       free(lat_S4_buf); free(lon_S4_buf);
       if ( strcmp(sensor, "AMSR2") == 0 ) {
         free(Tc_S5_buf);
         free(qual_S5_buf);
         free(lat_S5_buf); free(lon_S5_buf);
         free(Tc_S6_buf);
         free(qual_S6_buf);
         free(lat_S6_buf); free(lon_S6_buf);
       }
     }
   }

   /*---------------- GPROF for ATMS -------------------------------*/

   /* GPROF for ATMS is produced at the 1C coordinate system 96 spots across scan */

   status = nc_open(gprof_file, NC_NETCDF4, &ncid);
   if ( status != 0 ) {
     printf("File not found %s\n", gprof_file);
   }
   else {
     have_gprof= 1;
     printf("Opened GPROF %s\n", gprof_file);
     if ((status= nc_inq_ncid(ncid, "S1", &grpid1)))
     ERR(status);

     /*--- GPROF latitude and longitude ---*/

     status= nc_inq_varid(grpid1, "Latitude", &lat_SG_id);
     if ( status != 0 )  ERR(status);
     status= nc_inq_varid(grpid1, "Longitude", &lon_SG_id);
     if ( status != 0 )  ERR(status);
     status= nc_inq_varid(grpid1, "qualityFlag", &qualityFlag_SG_id);
     if ( status != 0 )  ERR(status);
     status= nc_inq_varid(grpid1, "surfacePrecipitation", &surfacePrecipitation_SG_id);
     if ( status != 0 )  ERR(status);

     if ( strcmp(sensor, "SAPHIR") != 0 ) {
       status= nc_inq_varid(grpid1, "surfaceTypeIndex", &sfc_SG_id);
       if ( status != 0 )  ERR(status);
       status= nc_inq_varid(grpid1, "probabilityOfPrecip", &probabilityOfPrecip_SG_id);
       if ( status != 0 )  ERR(status);
       status= nc_inq_varid(grpid1, "temp2mIndex", &temp2mIndex_SG_id);
       if ( status != 0 )  ERR(status);
       status= nc_inq_varid(grpid1, "totalColumnWaterVaporIndex", &totalColumnWaterVaporIndex_SG_id);
       if ( status != 0 )  ERR(status);
       status= nc_inq_varid(grpid1, "frozenPrecipitation", &frozenPrecipitation_SG_id);
       if ( status != 0 )  ERR(status);
     }

     status = nc_inq_var(grpid1, lat_SG_id, 0, &lat_SG_type, &lat_SG_ndims, lat_SG_dimids, &lat_SG_natts);
     if ( status != 0 )  ERR(status);
     printf("latitude ndims= %d\n", lat_SG_ndims);
     status = nc_inq_dimlen(grpid1, lat_SG_dimids[0], &nl_SG);
     status = nc_inq_dimlen(grpid1, lat_SG_dimids[1], &ns_SG);
     printf("GPROF latitude dims= %d %d\n", (int) nl_SG, (int) ns_SG);

     lat_GPROF_buf= (float *) malloc(nl_SG*ns_SG*sizeof(float));
     lon_GPROF_buf= (float *) malloc(nl_SG*ns_SG*sizeof(float));
     qualityFlag_GPROF_buf= (signed char *) malloc(nl_SG*ns_SG*sizeof(signed char));
     surfacePrecipitation_GPROF_buf= (float *) malloc(nl_SG*ns_SG*sizeof(float));
     sfc_GPROF_buf= (signed char *) malloc(nl_SG*ns_SG*sizeof(signed char));
     probabilityOfPrecip_GPROF_buf= (float *) malloc(nl_SG*ns_SG*sizeof(float));
     totalColumnWaterVaporIndex_GPROF_buf= (signed char *) malloc(nl_SG*ns_SG*sizeof(signed char));
     temp2mIndex_GPROF_buf= (short *) malloc(nl_SG*ns_SG*sizeof(short));
     frozenPrecipitation_GPROF_buf= (float *) malloc(nl_SG*ns_SG*sizeof(float));

     status = nc_get_var_float(grpid1, lat_SG_id, lat_GPROF_buf);
     if ( status != 0 ) ERR(status);
     status = nc_get_var_float(grpid1, lon_SG_id, lon_GPROF_buf);
     if ( status != 0 ) ERR(status);
     status = nc_get_var_schar(grpid1, qualityFlag_SG_id, qualityFlag_GPROF_buf);
     if ( status != 0 ) ERR(status);
     status = nc_get_var_float(grpid1, surfacePrecipitation_SG_id, surfacePrecipitation_GPROF_buf);
     if ( status != 0 ) ERR(status);

     if ( strcmp(sensor, "SAPHIR") != 0 ) {
       status = nc_get_var_schar(grpid1, sfc_SG_id, sfc_GPROF_buf);
       if ( status != 0 ) ERR(status);
       status = nc_get_var_float(grpid1, probabilityOfPrecip_SG_id, probabilityOfPrecip_GPROF_buf);
       if ( status != 0 ) ERR(status);
       status = nc_get_var_short(grpid1, temp2mIndex_SG_id, temp2mIndex_GPROF_buf);
       if ( status != 0 ) ERR(status);
       status = nc_get_var_schar(grpid1, totalColumnWaterVaporIndex_SG_id, totalColumnWaterVaporIndex_GPROF_buf);
       if ( status != 0 ) ERR(status);
       status = nc_get_var_float(grpid1, frozenPrecipitation_SG_id, frozenPrecipitation_GPROF_buf);
       if ( status != 0 ) ERR(status);
     }
     else {
       for (i=0; i < nl_SG*ns_SG; i++) {
         sfc_GPROF_buf[i]= -99;
         probabilityOfPrecip_GPROF_buf[i]= -999.0;
         temp2mIndex_GPROF_buf[i]= -999.0;
         totalColumnWaterVaporIndex_GPROF_buf[i]= -99;
         frozenPrecipitation_GPROF_buf[i]= -999.0;
       }
     }

     if ((retval = nc_close(ncid))) ERR(retval);
     printf("Closed %s\n", gprof_file);

     sfc_GPROF_S1_buf= (signed char *) malloc(nl_S1*ns_S1*sizeof(signed char));
     surfacePrecipitation_GPROF_S1_buf= (float *) malloc(nl_S1*ns_S1*sizeof(float));
     probabilityOfPrecip_GPROF_S1_buf= (float *) malloc(nl_S1*ns_S1*sizeof(float));
     totalColumnWaterVaporIndex_GPROF_S1_buf= (signed char *) malloc(nl_S1*ns_S1*sizeof(signed char));
     temp2mIndex_GPROF_S1_buf= (short *) malloc(nl_S1*ns_S1*sizeof(short));
     frozenPrecipitation_GPROF_S1_buf= (float *) malloc(nl_S1*ns_S1*sizeof(float));

     for (i=0; i < nl_S1; i++) {
       for (j=0; j < ns_S1; j++) {
         idx= i*ns_S1 + j;
         if ( strcmp(sensor, "SSMIS") == 0 ) idx2= i*ns_SG + (2*j+1);
         else if ( strcmp(sensor, "AMSR2") == 0 ) idx2= i*ns_SG + (2*j+1);
         else idx2= idx;

         sfc_GPROF_S1_buf[idx]= sfc_GPROF_buf[idx2];
         surfacePrecipitation_GPROF_S1_buf[idx]= surfacePrecipitation_GPROF_buf[idx2];
         probabilityOfPrecip_GPROF_S1_buf[idx]= probabilityOfPrecip_GPROF_buf[idx2];
         totalColumnWaterVaporIndex_GPROF_S1_buf[idx]= totalColumnWaterVaporIndex_GPROF_buf[idx2];
         temp2mIndex_GPROF_S1_buf[idx]= temp2mIndex_GPROF_buf[idx2];
         frozenPrecipitation_GPROF_S1_buf[idx]= frozenPrecipitation_GPROF_buf[idx2];

         lat1= lat_S1_buf[idx];
         lon1= lon_S1_buf[idx];
         lat2= lat_GPROF_buf[idx2];
         lon2= lon_GPROF_buf[idx2];
         km= RADEARTH*acos(cos(DTR*lon1-DTR*lon2)*cos(DTR*lat1)*cos(DTR*lat2) + sin(DTR*lat1)*sin(DTR*lat2));

/*
         printf("%4d %4d ", i,j);
         printf("%8.3f %8.3f ", lat_S1_buf[idx], lon_S1_buf[idx]);
         printf("%8.3f %8.3f ", lat_GPROF_buf[idx2], lon_GPROF_buf[idx2]);
         printf("%8.3f ", km);
         printf("%6d ", sfc_GPROF_S1_buf[idx]);
         printf("%9.3f ", surfacePrecipitation_GPROF_S1_buf[idx]);
         printf("%9.3f ", probabilityOfPrecip_GPROF_S1_buf[idx]);
         printf("%9.3f ", frozenPrecipitation_GPROF_S1_buf[idx]);
         printf("%6d ", totalColumnWaterVaporIndex_GPROF_S1_buf[idx]);
         printf("%6d ", temp2mIndex_GPROF_S1_buf[idx]);
         printf("\n");
*/
       }
     }

     free(lat_GPROF_buf);
     free(lon_GPROF_buf);
     free(sfc_GPROF_buf);
     free(qualityFlag_GPROF_buf);
     free(surfacePrecipitation_GPROF_buf);
     free(probabilityOfPrecip_GPROF_buf);
     free(temp2mIndex_GPROF_buf);
     free(totalColumnWaterVaporIndex_GPROF_buf);
     free(frozenPrecipitation_GPROF_buf);
   }

   /*---------------- Read DPR -----------------------------------*/

   status = nc_open(dpr_file, NC_NETCDF4, &ncid);
   if ( status != 0 ) {
     printf("File not found %s\n", dpr_file);
   }
   else {
     have_dpr= 1;
     printf("Opened %s\n", dpr_file);

     if ((status= nc_inq_ncid(ncid, "NS", &grpid1))) ERR(status);
     if ((status= nc_inq_ncid(grpid1, "scanStatus", &grpid1_scanstatus))) ERR(status);
     if ((status= nc_inq_ncid(grpid1, "PRE", &grpid1_pre))) ERR(status);
     if ((status= nc_inq_ncid(grpid1, "SLV", &grpid1_slv))) ERR(status);
     if ((status= nc_inq_ncid(grpid1, "CSF", &grpid1_csf))) ERR(status);
     if ((status= nc_inq_ncid(grpid1, "ScanTime", &grpid1_scantime))) ERR(status);

     if ((status= nc_inq_ncid(ncid, "MS", &grpid2))) ERR(status);
     if ((status= nc_inq_ncid(grpid2, "scanStatus", &grpid2_scanstatus))) ERR(status);
     if ((status= nc_inq_ncid(grpid2, "PRE", &grpid2_pre))) ERR(status);
     if ((status= nc_inq_ncid(grpid2, "SLV", &grpid2_slv))) ERR(status);

     /*--- NS zFactorMeasured ---*/

     status= nc_inq_varid(grpid1_pre, "zFactorMeasured", &zFactorMeasured_NS_id);
     if ( status != 0 )  ERR(status);

     status = nc_inq_var(grpid1_pre, zFactorMeasured_NS_id, 0, &zFactorMeasured_NS_type, &zFactorMeasured_NS_ndims, zFactorMeasured_NS_dimids, &zFactorMeasured_NS_natts);
     if ( status != 0 )  ERR(status);
     printf("NS zFactorMeasured ndims= %d\n", zFactorMeasured_NS_ndims);
     status = nc_inq_dimlen(grpid1_pre, zFactorMeasured_NS_dimids[0], &nl_NS);
     status = nc_inq_dimlen(grpid1_pre, zFactorMeasured_NS_dimids[1], &ns_NS);
     status = nc_inq_dimlen(grpid1_pre, zFactorMeasured_NS_dimids[2], &nlev_NS);
     printf("NS zFactorMeasured dims= %d %d %d\n", (int) nl_NS, (int) ns_NS, (int) nlev_NS);
     if ( nlev_NS != 2*NLEV_DPR ) {
       printf("ERROR  nlevels in DPR NS=%d   nlevels in DB=%d\n", nlev_NS, NLEV_DPR);
       exit(1);
     }

     zFactorMeasured_NS_buf= (float *) malloc(nl_NS*ns_NS*nlev_NS*sizeof(float));
     status = nc_get_var_float(grpid1_pre, zFactorMeasured_NS_id, zFactorMeasured_NS_buf);
     if ( status != 0 ) ERR(status);

     /*--- NS surface bin index ---*/
     status= nc_inq_varid(grpid1_pre, "binRealSurface", &binRealSurface_NS_id);
     if ( status != 0 )  ERR(status);
     binRealSurface_NS_buf= (short *) malloc(nl_NS*ns_NS*sizeof(short));
     status = nc_get_var_short(grpid1_pre, binRealSurface_NS_id, binRealSurface_NS_buf);
     if ( status != 0 ) ERR(status);

     /*--- NS no-clutter bin index ---*/
     status= nc_inq_varid(grpid1_pre, "binClutterFreeBottom", &binClutterFreeBottom_NS_id);
     if ( status != 0 )  ERR(status);
     binClutterFreeBottom_NS_buf= (short *) malloc(nl_NS*ns_NS*sizeof(short));
     status = nc_get_var_short(grpid1_pre, binClutterFreeBottom_NS_id, binClutterFreeBottom_NS_buf);
     if ( status != 0 ) ERR(status);

     /*--- NS data quality ---*/

     status= nc_inq_varid(grpid1_scanstatus, "dataQuality", &qual_NS_id);
     if ( status != 0 )  ERR(status);
     qual_NS_buf= (signed char *) malloc(nl_NS*sizeof(unsigned char));
     status = nc_get_var_schar(grpid1_scanstatus, qual_NS_id, qual_NS_buf);
     if ( status != 0 ) ERR(status);

     /*--- NS surface rain ---*/

     status= nc_inq_varid(grpid1_slv, "precipRateESurface", &precipRateESurface_NS_id);
     if ( status != 0 )  ERR(status);
     precipRateESurface_NS_buf= (float *) malloc(nl_NS*ns_NS*sizeof(float));
     status = nc_get_var_float(grpid1_slv, precipRateESurface_NS_id, precipRateESurface_NS_buf);
     if ( status != 0 ) ERR(status);

     /*--- NS latitude and longitude ---*/

     status= nc_inq_varid(grpid1, "Latitude", &lat_NS_id);
     if ( status != 0 )  ERR(status);
     status= nc_inq_varid(grpid1, "Longitude", &lon_NS_id);
     if ( status != 0 )  ERR(status);
     lat_NS_buf= (float *) malloc(nl_NS*ns_NS*sizeof(float));
     lon_NS_buf= (float *) malloc(nl_NS*ns_NS*sizeof(float));
     status = nc_get_var_float(grpid1, lat_NS_id, lat_NS_buf);
     if ( status != 0 ) ERR(status);
     status = nc_get_var_float(grpid1, lon_NS_id, lon_NS_buf);
     if ( status != 0 ) ERR(status);

     /*--- NS elevation ---*/

     status= nc_inq_varid(grpid1_pre, "elevation", &elev_NS_id);
     if ( status != 0 )  ERR(status);
     elev_NS_buf= (float *) malloc(nl_NS*ns_NS*sizeof(float));
     status = nc_get_var_float(grpid1_pre, elev_NS_id, elev_NS_buf);
     if ( status != 0 ) ERR(status);

     /*--- NS flag shallow rain ---*/

     status= nc_inq_varid(grpid1_csf, "flagShallowRain", &flagShallowRain_NS_id);
     if ( status != 0 )  ERR(status);
     flagShallowRain_NS_buf= (int *) malloc(nl_NS*ns_NS*sizeof(int));
     status = nc_get_var_int(grpid1_csf, flagShallowRain_NS_id, flagShallowRain_NS_buf);
     if ( status != 0 ) ERR(status);

     /*--- NS type precip ---*/

     status= nc_inq_varid(grpid1_csf, "typePrecip", &typePrecip_NS_id);
     if ( status != 0 )  ERR(status);
     typePrecip_NS_buf= (int *) malloc(nl_NS*ns_NS*sizeof(int));
     status = nc_get_var_int(grpid1_csf, typePrecip_NS_id, typePrecip_NS_buf);
     if ( status != 0 ) ERR(status);

  /*--- NS date ---*/

   status= nc_inq_varid(grpid1_scantime, "Year", &yyyy_NS_id);
   if ( status != 0 )  ERR(status);
   yyyy_NS_buf= (short *) malloc(nl_NS*sizeof(short));
   status = nc_get_var_short(grpid1_scantime, yyyy_NS_id, yyyy_NS_buf);
   if ( status != 0 ) ERR(status);

   status= nc_inq_varid(grpid1_scantime, "Month", &mm_NS_id);
   if ( status != 0 )  ERR(status);
   mm_NS_buf= (short *) malloc(nl_NS*sizeof(short));
   status = nc_get_var_short(grpid1_scantime, mm_NS_id, mm_NS_buf);
   if ( status != 0 ) ERR(status);

   status= nc_inq_varid(grpid1_scantime, "DayOfMonth", &dd_NS_id);
   if ( status != 0 )  ERR(status);
   dd_NS_buf= (short *) malloc(nl_NS*sizeof(short));
   status = nc_get_var_short(grpid1_scantime, dd_NS_id, dd_NS_buf);
   if ( status != 0 ) ERR(status);

   status= nc_inq_varid(grpid1_scantime, "Hour", &hh_NS_id);
   if ( status != 0 )  ERR(status);
   hh_NS_buf= (short *) malloc(nl_NS*sizeof(short));
   status = nc_get_var_short(grpid1_scantime, hh_NS_id, hh_NS_buf);
   if ( status != 0 ) ERR(status);

   status= nc_inq_varid(grpid1_scantime, "Minute", &mn_NS_id);
   if ( status != 0 )  ERR(status);
   mn_NS_buf= (short *) malloc(nl_NS*sizeof(short));
   status = nc_get_var_short(grpid1_scantime, mn_NS_id, mn_NS_buf);
   if ( status != 0 ) ERR(status);

   status= nc_inq_varid(grpid1_scantime, "Second", &ss_NS_id);
   if ( status != 0 )  ERR(status);
   ss_NS_buf= (short *) malloc(nl_NS*sizeof(short));
   status = nc_get_var_short(grpid1_scantime, ss_NS_id, ss_NS_buf);
   if ( status != 0 ) ERR(status);


     /*--- MS Z measured ---*/

     status= nc_inq_varid(grpid2_pre, "zFactorMeasured", &zFactorMeasured_MS_id);
     if ( status != 0 )  ERR(status);

     status = nc_inq_var(grpid2_pre, zFactorMeasured_MS_id, 0, &zFactorMeasured_MS_type, &zFactorMeasured_MS_ndims, zFactorMeasured_MS_dimids, &zFactorMeasured_MS_natts);
     if ( status != 0 )  ERR(status);
     printf("MS zFactorMeasured ndims= %d\n", zFactorMeasured_MS_ndims);
     status = nc_inq_dimlen(grpid2_pre, zFactorMeasured_MS_dimids[0], &nl_MS);
     status = nc_inq_dimlen(grpid2_pre, zFactorMeasured_MS_dimids[1], &ns_MS);
     status = nc_inq_dimlen(grpid2_pre, zFactorMeasured_MS_dimids[2], &nlev_MS);
     printf("MS zFactorMeasured dims= %d %d %d\n", (int) nl_MS, (int) ns_MS, (int) nlev_MS);
     if ( nlev_MS != 2*NLEV_DPR ) {
       printf("ERROR  nlevels in DPR MS=%d   nlevels in DB=%d\n", nlev_MS, NLEV_DPR);
       exit(1);
     }

     zFactorMeasured_MS_buf= (float *) malloc(nl_MS*ns_MS*nlev_MS*sizeof(float));
     status = nc_get_var_float(grpid2_pre, zFactorMeasured_MS_id, zFactorMeasured_MS_buf);
     if ( status != 0 ) ERR(status);

     /*--- MS surface bin index ---*/
     status= nc_inq_varid(grpid2_pre, "binRealSurface", &binRealSurface_MS_id);
     if ( status != 0 )  ERR(status);
     binRealSurface_MS_buf= (short *) malloc(nl_MS*ns_MS*sizeof(short));
     status = nc_get_var_short(grpid2_pre, binRealSurface_MS_id, binRealSurface_MS_buf);
     if ( status != 0 ) ERR(status);

     /*--- MS no-clutter bin index ---*/
     status= nc_inq_varid(grpid2_pre, "binClutterFreeBottom", &binClutterFreeBottom_MS_id);
     if ( status != 0 )  ERR(status);
     binClutterFreeBottom_MS_buf= (short *) malloc(nl_MS*ns_MS*sizeof(short));
     status = nc_get_var_short(grpid2_pre, binClutterFreeBottom_MS_id, binClutterFreeBottom_MS_buf);
     if ( status != 0 ) ERR(status);

     /*--- MS surface rain ---*/

     status= nc_inq_varid(grpid2_slv, "precipRateESurface", &precipRateESurface_MS_id);
     if ( status != 0 )  ERR(status);
     precipRateESurface_MS_buf= (float *) malloc(nl_MS*ns_MS*sizeof(float));
     status = nc_get_var_float(grpid2_slv, precipRateESurface_MS_id, precipRateESurface_MS_buf);
     if ( status != 0 ) ERR(status);

     /*--- MS data quality ---*/

     status= nc_inq_varid(grpid2_scanstatus, "dataQuality", &qual_MS_id);
     if ( status != 0 )  ERR(status);
     qual_MS_buf= (signed char *) malloc(nl_MS*sizeof(unsigned char));
     status = nc_get_var_schar(grpid2_scanstatus, qual_MS_id, qual_MS_buf);
     if ( status != 0 ) ERR(status);


     if ((retval = nc_close(ncid))) ERR(retval);
     printf("Closed %s\n", dpr_file);

     secs_NS_buf= (double *) malloc(nl_NS*sizeof(double));
     for (i=0; i < nl_NS; i++) {
       iyyyy= yyyy_NS_buf[i];
       imm= mm_NS_buf[i];
       idd= dd_NS_buf[i];
       ihh= hh_NS_buf[i];
       imn= mn_NS_buf[i];
       iss= ss_NS_buf[i];
       tm.tm_year= iyyyy-1900;
       tm.tm_mon= imm-1;
       tm.tm_mday= idd;
       tm.tm_hour= ihh;
       tm.tm_min= imn;
       tm.tm_sec= iss;
       t0= timegm(&tm);
       secs_NS_buf[i]= ( double ) t0;
     }

   }

   /*---------------- Read COMBINED -----------------------------------*/

   status = nc_open(cmb_file, NC_NETCDF4, &ncid);
   if ( status != 0 ) {
     printf("File not found %s\n", cmb_file);
   }
   else {
     have_cmb= 1;
     printf("Opened %s\n", cmb_file);

     if ((status= nc_inq_ncid(ncid, "NS", &grpid1))) ERR(status);
     if ((status= nc_inq_ncid(ncid, "MS", &grpid2))) ERR(status);
     if ((status= nc_inq_ncid(grpid1, "scanStatus", &grpid1_scanstatus))) ERR(status);
     if ((status= nc_inq_ncid(grpid2, "scanStatus", &grpid2_scanstatus))) ERR(status);

     /*--- latitude and longitude ---*/

     status= nc_inq_varid(grpid1, "Latitude", &lat_NS_cmb_id);
     if ( status != 0 )  ERR(status);
     status= nc_inq_varid(grpid1, "Longitude", &lon_NS_cmb_id);
     if ( status != 0 )  ERR(status);
     status= nc_inq_varid(grpid1, "surfPrecipTotRate", &surfPrecipTotRate_NS_cmb_id);
     if ( status != 0 )  ERR(status);
     status= nc_inq_varid(grpid1, "precipTotRate", &PrecipTotRate_NS_cmb_id);
     if ( status != 0 )  ERR(status);
     status= nc_inq_varid(grpid1_scanstatus, "dataQuality", &dataQuality_NS_cmb_id);
     if ( status != 0 )  ERR(status);

     status= nc_inq_varid(grpid2, "Latitude", &lat_MS_cmb_id);
     if ( status != 0 )  ERR(status);
     status= nc_inq_varid(grpid2, "Longitude", &lon_MS_cmb_id);
     if ( status != 0 )  ERR(status);
     status= nc_inq_varid(grpid2, "surfPrecipTotRate", &surfPrecipTotRate_MS_cmb_id);
     if ( status != 0 )  ERR(status);
     status= nc_inq_varid(grpid2, "precipTotRate", &PrecipTotRate_MS_cmb_id);
     if ( status != 0 )  ERR(status);
     status= nc_inq_varid(grpid2_scanstatus, "dataQuality", &dataQuality_MS_cmb_id);
     if ( status != 0 )  ERR(status);

     status = nc_inq_var(grpid1, PrecipTotRate_NS_cmb_id, 0, &PrecipTotRate_NS_cmb_type, &PrecipTotRate_NS_cmb_ndims, PrecipTotRate_NS_cmb_dimids, &PrecipTotRate_NS_cmb_natts);
     if ( status != 0 )  ERR(status);
     printf("PrecipTotRate_NS_cmb ndims= %d\n", PrecipTotRate_NS_cmb_ndims);
     status = nc_inq_dimlen(grpid1, PrecipTotRate_NS_cmb_dimids[0], &nl_NS_cmb);
     status = nc_inq_dimlen(grpid1, PrecipTotRate_NS_cmb_dimids[1], &ns_NS_cmb);
     status = nc_inq_dimlen(grpid1, PrecipTotRate_NS_cmb_dimids[2], &nlev_NS_cmb);
     printf("COMBINED NS precip profile dims= %d %d %d\n", (int) nl_NS_cmb, (int) ns_NS_cmb, (int) nlev_NS_cmb);
     if ( nlev_NS_cmb != NLEV_DPR ) {
       printf("ERROR  nlevels in PrecipTotRate_NS_cmb=%d   nlevels in DB=%d\n", nlev_NS_cmb, NLEV_DPR);
       exit(1);
     }

     status = nc_inq_var(grpid1, PrecipTotRate_MS_cmb_id, 0, &PrecipTotRate_MS_cmb_type, &PrecipTotRate_MS_cmb_ndims, PrecipTotRate_MS_cmb_dimids, &PrecipTotRate_MS_cmb_natts);
     if ( status != 0 )  ERR(status);
     printf("PrecipTotRate_MS_cmb ndims= %d\n", PrecipTotRate_MS_cmb_ndims);
     status = nc_inq_dimlen(grpid1, PrecipTotRate_MS_cmb_dimids[0], &nl_MS_cmb);
     status = nc_inq_dimlen(grpid1, PrecipTotRate_MS_cmb_dimids[1], &ns_MS_cmb);
     status = nc_inq_dimlen(grpid1, PrecipTotRate_MS_cmb_dimids[2], &nlev_MS_cmb);
     printf("COMBINED MS precip profile dims= %d %d %d\n", (int) nl_MS_cmb, (int) ns_MS_cmb, (int) nlev_MS_cmb);
     if ( nlev_MS_cmb != NLEV_DPR ) {
       printf("ERROR  nlevels in PrecipTotRate_MS_cmb=%d   nlevels in DB=%d\n", nlev_MS_cmb, NLEV_DPR);
       exit(1);
     }

     status = nc_inq_var(grpid1, lat_NS_cmb_id, 0, &lat_NS_cmb_type, &lat_NS_cmb_ndims, lat_NS_cmb_dimids, &lat_NS_cmb_natts);
     if ( status != 0 )  ERR(status);
     printf("lat_NS_cmb ndims= %d\n", lat_NS_cmb_ndims);
     status = nc_inq_dimlen(grpid1, lat_NS_cmb_dimids[0], &nl_NS_cmb);
     status = nc_inq_dimlen(grpid1, lat_NS_cmb_dimids[1], &ns_NS_cmb);
     printf("COMBINED NS dims= %d %d\n", (int) nl_NS_cmb, (int) ns_NS_cmb);

     status = nc_inq_var(grpid2, lat_MS_cmb_id, 0, &lat_MS_cmb_type, &lat_MS_cmb_ndims, lat_MS_cmb_dimids, &lat_MS_cmb_natts);
     if ( status != 0 )  ERR(status);
     printf("lat_MS_cmb ndims= %d\n", lat_MS_cmb_ndims);
     status = nc_inq_dimlen(grpid2, lat_MS_cmb_dimids[0], &nl_MS_cmb);
     status = nc_inq_dimlen(grpid2, lat_MS_cmb_dimids[1], &ns_MS_cmb);
     printf("COMBINED MS dims= %d %d\n", (int) nl_MS_cmb, (int) ns_MS_cmb);

     lat_NS_cmb_buf= (float *) malloc(nl_NS_cmb*ns_NS_cmb*sizeof(float));
     lon_NS_cmb_buf= (float *) malloc(nl_NS_cmb*ns_NS_cmb*sizeof(float));
     surfPrecipTotRate_NS_cmb_buf= (float *) malloc(nl_NS_cmb*ns_NS_cmb*sizeof(float));
     PrecipTotRate_NS_cmb_buf= (float *) malloc(nl_NS_cmb*ns_NS_cmb*nlev_NS_cmb*sizeof(float));
     dataQuality_NS_cmb_buf= (signed char *) malloc(nl_NS_cmb*sizeof(signed char));

     lat_MS_cmb_buf= (float *) malloc(nl_MS_cmb*ns_MS_cmb*sizeof(float));
     lon_MS_cmb_buf= (float *) malloc(nl_MS_cmb*ns_MS_cmb*sizeof(float));
     surfPrecipTotRate_MS_cmb_buf= (float *) malloc(nl_MS_cmb*ns_MS_cmb*sizeof(float));
     PrecipTotRate_MS_cmb_buf= (float *) malloc(nl_MS_cmb*ns_MS_cmb*nlev_MS_cmb*sizeof(float));
     dataQuality_MS_cmb_buf= (signed char *) malloc(nl_MS_cmb*sizeof(signed char));

     status = nc_get_var_float(grpid1, lat_NS_cmb_id, lat_NS_cmb_buf);
     if ( status != 0 ) ERR(status);
     status = nc_get_var_float(grpid1, lon_NS_cmb_id, lon_NS_cmb_buf);
     if ( status != 0 ) ERR(status);
     status = nc_get_var_float(grpid1, surfPrecipTotRate_NS_cmb_id, surfPrecipTotRate_NS_cmb_buf);
     if ( status != 0 ) ERR(status);
     status = nc_get_var_float(grpid1, PrecipTotRate_NS_cmb_id, PrecipTotRate_NS_cmb_buf);
     if ( status != 0 ) ERR(status);
     status = nc_get_var_schar(grpid1_scanstatus, dataQuality_NS_cmb_id, dataQuality_NS_cmb_buf);
     if ( status != 0 ) ERR(status);

     status = nc_get_var_float(grpid2, lat_MS_cmb_id, lat_MS_cmb_buf);
     if ( status != 0 ) ERR(status);
     status = nc_get_var_float(grpid2, lon_MS_cmb_id, lon_MS_cmb_buf);
     if ( status != 0 ) ERR(status);
     status = nc_get_var_float(grpid2, surfPrecipTotRate_MS_cmb_id, surfPrecipTotRate_MS_cmb_buf);
     if ( status != 0 ) ERR(status);
     status = nc_get_var_float(grpid2, PrecipTotRate_MS_cmb_id, PrecipTotRate_MS_cmb_buf);
     if ( status != 0 ) ERR(status);
     status = nc_get_var_schar(grpid2_scanstatus, dataQuality_MS_cmb_id, dataQuality_MS_cmb_buf);
     if ( status != 0 ) ERR(status);

     if ((retval = nc_close(ncid))) ERR(retval);
     printf("Closed %s\n", cmb_file);

     if ( nl_NS != nl_NS_cmb || ns_NS != ns_NS_cmb ||  nl_MS != nl_MS_cmb || ns_MS != ns_MS_cmb ) {
       printf("CMB dimensions are not parallel to DPR\n");
       have_cmb= -1;
     }
     else {
       for (i=0; i < nl_NS; i++) {
         idx= i*ns_NS + (ns_NS/2);
         glat1= lat_NS_buf[idx];
         glon1= lon_NS_buf[idx];
         glat2= lat_NS_cmb_buf[idx];
         glon2= lon_NS_cmb_buf[idx];
         if ( fabs(glat1-glat2) > 0.001 ) {
           km= RADEARTH*acos(cos(DTR*glon1-DTR*glon2)*cos(DTR*glat1)*cos(DTR*glat2) + sin(DTR*glat1)*sin(DTR*glat2));
           printf("%4d %8.3f %8.3f %8.3f %8.3f %f\n",i,glat1,glon1,glat2,glon2, km);
           printf("CMB coordinates are not the same as DPR coordinates\n");
           have_cmb= -1;
           break;
         }
       }
     }

   }

   /*--------------------------------------------*/

   irev= FractionalGranuleNumber_S1_buf[nl_S1/2];

   if ( isat_start >= 0 && isat_end > isat_start ) {
     printf("Processing orbit between scan lines %d %d\n", isat_start, isat_end);
   }
   else if ( fabs(tlat) < 90 && fabs(tlon) < 180 ) {
     km_min= 1.0E8;
     i= -1;
     for (isat=0; isat < nl_S1; isat++) {
       idx= isat*ns_S1 + (ns_S1/2);
       glat1= lat_S1_buf[idx];
       glon1= lon_S1_buf[idx];
       km= RADEARTH*acos(cos(DTR*glon1-DTR*tlon)*cos(DTR*glat1)*cos(DTR*tlat) + sin(DTR*glat1)*sin(DTR*tlat));
       if ( km < km_min ) {
         km_min= km;
         i= isat;
         clat= glat1;
         clon= glon1;
       }
     }
     if ( i < 0 || km_min > 800 ) {
       printf("Starting line not found  km_min=%f\n", km_min);
       exit(1);
     }
     printf("Starting line %d   km_min=%f\n", i, km_min);
     isat_start= i-nlines;
     if ( isat_start < 0 ) isat_start= 0;
     isat_end= i+nlines;
     if ( isat_end >= nl_S1 ) isat_end= nl_S1-1;
   }
   else {
     isat_start= 0;
     isat_end= nl_S1-1;
   }

   if ( isat_start < 0 || isat_end < 0 ) {
     printf("Starting line not found\n");
     exit(1);
   }

   nl_out= isat_end - isat_start + 1;
   printf("Lines=%d  Range=%d %d\n", nl_out, isat_start, isat_end);

   i= (isat_end + isat_start)/2;
   idx= i*ns_S1 + (ns_S1/2);
   clat= lat_S1_buf[idx];
   clon= lon_S1_buf[idx];
   iyyyy= yyyy_S1_buf[i];
   imm= mm_S1_buf[i];
   idd= dd_S1_buf[i];
   ihh= hh_S1_buf[i];
   imn= mn_S1_buf[i];
   iss= ss_S1_buf[i];
   sprintf(cdate,"%04d/%02d/%02d %02d:%02d:%02d", iyyyy, imm, idd, ihh, imn, iss);
   sprintf(outdate,"%04d%02d%02d_%02d%02d", iyyyy, imm, idd, ihh, imn);
   sprintf(outfilename, "%s/%s_EPC_%06d_%s_%04d_%04d", outdir, satname, irev, outdate, isat_start, isat_end);

   i= isat_start;
   idx= i*ns_S1 + (ns_S1/2);
   clat1= lat_S1_buf[idx];
   clon1= lon_S1_buf[idx];
   iyyyy= yyyy_S1_buf[i];
   imm= mm_S1_buf[i];
   idd= dd_S1_buf[i];
   ihh= hh_S1_buf[i];
   imn= mn_S1_buf[i];
   iss= ss_S1_buf[i];
   sprintf(cdate1,"%04d/%02d/%02d %02d:%02d:%02d", iyyyy, imm, idd, ihh, imn, iss);
   sprintf(cdate12,"%04d%02d%02d_%02d%02d", iyyyy, imm, idd, ihh, imn, iss);

   i= isat_end;
   idx= i*ns_S1 + (ns_S1/2);
   clat2= lat_S1_buf[idx];
   clon2= lon_S1_buf[idx];
   iyyyy= yyyy_S1_buf[i];
   imm= mm_S1_buf[i];
   idd= dd_S1_buf[i];
   ihh= hh_S1_buf[i];
   imn= mn_S1_buf[i];
   iss= ss_S1_buf[i];
   sprintf(cdate2,"%04d/%02d/%02d %02d:%02d:%02d", iyyyy, imm, idd, ihh, imn, iss);
   sprintf(cdate12,"%s_%04d%02d%02d_%02d%02d", cdate12, iyyyy, imm, idd, ihh, imn, iss);

   /* For full orbit, name the file with beginning and end times, rather than center time */
   if ( fabs(tlat) > 90 ) {
     sprintf(outdate,"%s", cdate12);
     sprintf(outfilename, "%s/%s_EPC_%06d_%s", outdir, satname, irev, outdate);
   }

   printf("First line coords=%8.3f %8.3f  date=%s\n", clat1, clon1, cdate1);
   printf("Last line coords= %8.3f %8.3f  date=%s\n", clat2, clon2, cdate2);
   printf("Center coords=    %8.3f %8.3f  date=%s\n", clat, clon, cdate);
   printf("Output filenames prefix= %s\n", outfilename);

   /*--------------------------------------------*/

   /*-- map DPR NS onto Radiometer S1 grid --*/

   if ( have_dpr > 0 ) {
     printf("Finding DPR NS locations on %s S1 grid\n", sensor);

     dpr_NS_i= (int *) malloc(nl_S1*ns_S1*sizeof(int));
     dpr_NS_j= (int *) malloc(nl_S1*ns_S1*sizeof(int));
     for (i= 0; i< nl_S1*ns_S1; i++) {
       dpr_NS_i[i]= -1;
       dpr_NS_j[i]= -1;
     }

     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         glat1= lat_S1_buf[idx];
         glon1= lon_S1_buf[idx];
         glon1_360= ( glon1 < 0.0 ) ? glon1+360.0 : glon1;

         /* Find nearest NS radar */
         /* speed up search - initialize */
         if ( i1_start < 0 || i1_end < 0 ) {
           i1_start= 0;
           i1_end= nl_NS;
         }

         km_min= 1000.0;
         imin_NS= jmin_NS= -1;
         for (i1=i1_start; i1 < i1_end; i1++) {
           if ( i1 < 0 || i1 >= nl_NS ) continue;
           if ( km_min < 2.0 ) break;
           for (j1=0; j1 < ns_NS; j1++) {
             if ( km_min < 2.0 ) break;
             idx1= i1*ns_NS + j1;
             nslat1= lat_NS_buf[idx1];
             nslon1= lon_NS_buf[idx1];
             if ( abs(nslat1-glat1) > 1.0 ) continue;
             nslon1_360= ( nslon1 < 0.0 ) ? nslon1+360.0 : nslon1;
             lon_diff= 180.0 - fabs(180.0 - abs(nslon1_360 - glon1_360));
             if ( lon_diff > 1.0 ) continue;
             km= RADEARTH*acos(cos(DTR*glon1-DTR*nslon1)*cos(DTR*glat1)*cos(DTR*nslat1) + sin(DTR*glat1)*sin(DTR*nslat1));
             /*printf("LOOP  S1=%4d %4d %8.3f %8.3f  NS=%4d %4d %8.3f %8.3f  km=%8.2f\n", isat,jsat,glat1,glon1, i1,j1,nslat1,nslon1, km);*/
             if ( km < km_min ) {
               km_min= km;
               imin_NS= i1; jmin_NS=j1;
             }
           }
         }
         if ( km_min > 5.0 ) continue;

         /* speed up search - update starting index */
         i1_start= imin_NS - 30;
         if ( i1_start < 0 ) {
           i1_start= 0;
           i1_end= nl_NS;
         }
         else
           i1_end= i1_start + 60;

         idx1= imin_NS*ns_NS + jmin_NS;
         nslat1= lat_NS_buf[idx1];
         nslon1= lon_NS_buf[idx1];
         dpr_NS_i[idx]= imin_NS;
         dpr_NS_j[idx]= jmin_NS;
         secs_min= secs_S1_buf[isat] - secs_NS_buf[imin_NS];

         n++;

         if ( jmin_NS == ns_NS/2 ) 
           printf("S1=%4d %4d %8.3f %8.3f  NS=%4d %4d %8.3f %8.3f   km=%8.2f secs=%6d\n",
           isat,jsat,glat1,glon1, imin_NS,jmin_NS,nslat1,nslon1, km_min, (int) secs_min);
       }
     }
     printf("Total S1 pixels=%d    DPR NS pixels inside swath=%d\n", nl_out*ns_S1, n);
     if ( n < 10 ) have_dpr= -1;
   }

   /*--------------------------------------------*/

   /* emis and pc_emis computed from input TB */
   emis_buf= (float *) malloc(nl_out*ns_S1*NEM*sizeof(float));
   pc_emis_buf= (float *) malloc(nl_out*ns_S1*NEM*sizeof(float));

   /* estimated precipitation   NS and MS database */
   precip1_NS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   precip2_NS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   precip1_prof_NS_buf= (short *) malloc(nl_out*ns_S1*NLEV_PRECIP*sizeof(short));
   precip2_prof_NS_buf= (short *) malloc(nl_out*ns_S1*NLEV_PRECIP*sizeof(short));
   precip1_prob_NS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   precip2_prob_NS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   epc_qual_NS_buf= (int *) malloc(nl_out*ns_S1*sizeof(int));
   type_precip_NS_buf= (short *) malloc(nl_out*ns_S1*sizeof(short));
   shallow_rain_NS_buf= (short *) malloc(nl_out*ns_S1*sizeof(short));

   precip1_MS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   precip2_MS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   precip1_prof_MS_buf= (short *) malloc(nl_out*ns_S1*NLEV_PRECIP*sizeof(short));
   precip2_prof_MS_buf= (short *) malloc(nl_out*ns_S1*NLEV_PRECIP*sizeof(short));
   precip1_prob_MS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   precip2_prob_MS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   epc_qual_MS_buf= (int *) malloc(nl_out*ns_S1*sizeof(int));
   type_precip_MS_buf= (short *) malloc(nl_out*ns_S1*sizeof(short));
   shallow_rain_MS_buf= (short *) malloc(nl_out*ns_S1*sizeof(short));

   /* additional estimated model parameters */
   ts_epc_NS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   t2m_epc_NS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   tqv_epc_NS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   h273_epc_NS_buf= (short *) malloc(nl_out*ns_S1*sizeof(short));
   ku_ztop_epc_NS_buf= (short *) malloc(nl_out*ns_S1*sizeof(short));

   ts_epc_MS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   t2m_epc_MS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   tqv_epc_MS_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   h273_epc_MS_buf= (short *) malloc(nl_out*ns_S1*sizeof(short));
   ku_ztop_epc_MS_buf= (short *) malloc(nl_out*ns_S1*sizeof(short));
   ka_ztop_epc_MS_buf= (short *) malloc(nl_out*ns_S1*sizeof(short));

   /* top-ranked Z and TB candidates for the NS and MS estimates */
   z_ku_1_NS_buf= (short *) malloc(nl_out*ns_S1*NLEV_DPR*sizeof(short));
   z_ku_1_MS_buf= (short *) malloc(nl_out*ns_S1*NLEV_DPR*sizeof(short));
   z_ka_1_MS_buf= (short *) malloc(nl_out*ns_S1*NLEV_DPR*sizeof(short));
   tb_1_NS_buf= (float *) malloc(nl_out*ns_S1*nc_out*sizeof(float));
   tb_1_MS_buf= (float *) malloc(nl_out*ns_S1*nc_out*sizeof(float));

   for (i=0; i< nl_out*ns_S1*NEM; i++) {
     emis_buf[i]= -999.0;
     pc_emis_buf[i]= -999.0;
   }
   for (i=0; i< nl_out*ns_S1; i++) {
     precip1_NS_buf[i]= -999.0;
     precip2_NS_buf[i]= -999.0;
     precip1_prob_NS_buf[i]= -999.0;
     precip2_prob_NS_buf[i]= -999.0;
     ts_epc_NS_buf[i]= -999.0;
     t2m_epc_NS_buf[i]= -999.0;
     h273_epc_NS_buf[i]= -9999;
     tqv_epc_NS_buf[i]= -999.0;
     epc_qual_NS_buf[i]= -999;  /* -2= no DB index found, set this as default */
     ku_ztop_epc_NS_buf[i]= -9999;
     type_precip_NS_buf[i]= -9999;
     shallow_rain_NS_buf[i]= -9999;

     precip1_MS_buf[i]= -999.0;
     precip2_MS_buf[i]= -999.0;
     precip1_prob_MS_buf[i]= -999.0;
     precip2_prob_MS_buf[i]= -999.0;
     ts_epc_MS_buf[i]= -999.0;
     t2m_epc_MS_buf[i]= -999.0;
     h273_epc_MS_buf[i]= -9999;
     tqv_epc_MS_buf[i]= -999.0;
     epc_qual_MS_buf[i]= -999;
     ku_ztop_epc_MS_buf[i]= -9999;
     ka_ztop_epc_MS_buf[i]= -9999;
     type_precip_MS_buf[i]= -9999;
     shallow_rain_MS_buf[i]= -9999;
   }
   for (i=0; i < nl_out*ns_S1*NLEV_DPR; i++) {
     z_ku_1_NS_buf[i]= -9999;
     z_ku_1_MS_buf[i]= z_ka_1_MS_buf[i]= -9999;
   }
   for (i=0; i< nl_out*ns_S1*nc_out; i++) {
     tb_1_NS_buf[i]= -999.0;
     tb_1_MS_buf[i]= -999.0;
   }
   for (i=0; i< nl_out*ns_S1*NLEV_PRECIP; i++) {
     precip1_prof_NS_buf[i]= -9999;
     precip2_prof_MS_buf[i]= -9999;
     precip1_prof_NS_buf[i]= -9999;
     precip2_prof_MS_buf[i]= -9999;
   }


   pixel_index = (int *) malloc(nl_S1*ns_S1*sizeof(int));
   for (i=0; i< nl_S1*ns_S1; i++) {
     pixel_index[i]= -1;
   }
   n_index = (int *) malloc(ndb_files*sizeof(int));
   for (i=0; i< ndb_files; i++) n_index[i]= 0;

   ncold_index = (int *) malloc(ndb_files*sizeof(int));
   for (i=0; i< ndb_files; i++) ncold_index[i]= 0;

   /*---------------------------------------------------------------*/

   /* Read MERRA2 interpolated model fields */

   status = nc_open(model_file, NC_NETCDF4, &ncid);
   if ( status != 0 ) {
     printf("Model file not found %s\n", model_file);
   }
   else {
     have_model= 1;
     printf("Opened %s\n", model_file);

     if ((status= nc_inq_varid(ncid, "t", &t_id))) ERR(status);
     if ((status = nc_inq_var(ncid, t_id, 0, &t_type, &t_ndims, t_dimids, &t_natts))) ERR(status);
     status = nc_inq_dimlen(ncid, t_dimids[0], &nlat);
     status = nc_inq_dimlen(ncid, t_dimids[1], &nlon);
     status = nc_inq_dimlen(ncid, t_dimids[2], &nlevel);
     t_buf= (float *) malloc(nlat*nlon*nlevel*sizeof(float));
     if ((status = nc_get_var_float(ncid, t_id, t_buf))) ERR(status);
     printf("T dims= %d %d %d\n", (int) nlat, (int) nlon, (int) nlevel);

     /* number of levels of the interpolated model fields is the same as the number of levels in the database entries */
     /* For MERRA2, there are 42 levels which are indexed from top to surface */
     if ( nlevel != 42 ) {
       printf("MERRA2 levels=%d are different than the 42-level database entries\n", (int) nlevel);
       exit(1);
     }
     /* interpolated model fields are on the identical grid as the S1 channels in the radiometer 1C file */
     if ( nlat != nl_S1 || nlon != ns_S1 ) {
       printf("Dimensions are not parallel  MERRA2=%d %d  SAT=%d %d\n", (int) nlat, (int) nlon, (int) nl_S1, (int) ns_S1);
       exit(1);
     }

     if ((status= nc_inq_varid(ncid, "time", &time_id))) ERR(status);
     time_buf= (double *) malloc(nlat*sizeof(double));
     if ((status = nc_get_var_double(ncid, time_id, time_buf))) ERR(status);
     printf("TIME dims= %d\n", (int) nlat);

     if ((status= nc_inq_varid(ncid, "h", &h_id))) ERR(status);
     h_buf= (float *) malloc(nlat*nlon*nlevel*sizeof(float));
     if ((status = nc_get_var_float(ncid, h_id, h_buf))) ERR(status);
     printf("H dims= %d %d %d\n", (int) nlat, (int) nlon, (int) nlevel);

     if ((status= nc_inq_varid(ncid, "qv", &qv_id))) ERR(status);
     qv_buf= (float *) malloc(nlat*nlon*nlevel*sizeof(float));
     if ((status = nc_get_var_float(ncid, qv_id, qv_buf))) ERR(status);
     printf("QV dims= %d %d %d\n", (int) nlat, (int) nlon, (int) nlevel);

     if ((status= nc_inq_varid(ncid, "ts", &ts_id))) ERR(status);
     ts_buf= (float *) malloc(nlat*nlon*sizeof(float));
     if ((status = nc_get_var_float(ncid, ts_id, ts_buf))) ERR(status);
     printf("TS dims= %d %d\n", (int) nlat, (int) nlon);

     if ((status= nc_inq_varid(ncid, "t2m", &t2m_id))) ERR(status);
     t2m_buf= (float *) malloc(nlat*nlon*sizeof(float));
     if ((status = nc_get_var_float(ncid, t2m_id, t2m_buf))) ERR(status);
     printf("T2M dims= %d %d\n", (int) nlat, (int) nlon);

     if ((status= nc_inq_varid(ncid, "t2m_dew", &t2m_dew_id))) ERR(status);
     t2m_dew_buf= (float *) malloc(nlat*nlon*sizeof(float));
     if ((status = nc_get_var_float(ncid, t2m_dew_id, t2m_dew_buf))) ERR(status);
     printf("T2M_DEW dims= %d %d\n", (int) nlat, (int) nlon);

     if ((status= nc_inq_varid(ncid, "t2m_wet", &t2m_wet_id))) ERR(status);
     t2m_wet_buf= (float *) malloc(nlat*nlon*sizeof(float));
     if ((status = nc_get_var_float(ncid, t2m_wet_id, t2m_wet_buf))) ERR(status);
     printf("T2M_WET dims= %d %d\n", (int) nlat, (int) nlon);

     if ((status= nc_inq_varid(ncid, "tqv", &tqv_id))) ERR(status);
     tqv_buf= (float *) malloc(nlat*nlon*sizeof(float));
     if ((status = nc_get_var_float(ncid, tqv_id, tqv_buf))) ERR(status);
     printf("TQV dims= %d %d\n", (int) nlat, (int) nlon);

     if ((status= nc_inq_varid(ncid, "hs", &hs_id))) ERR(status);
     hs_buf= (float *) malloc(nlat*nlon*sizeof(float));
     if ((status = nc_get_var_float(ncid, hs_id, hs_buf))) ERR(status);
     printf("HS dims= %d %d\n", (int) nlat, (int) nlon);

     if ((status= nc_inq_varid(ncid, "ps", &ps_id))) ERR(status);
     ps_buf= (float *) malloc(nlat*nlon*sizeof(float));
     if ((status = nc_get_var_float(ncid, ps_id, ps_buf))) ERR(status);
     printf("PS dims= %d %d\n", (int) nlat, (int) nlon);

     if ((status= nc_inq_varid(ncid, "latitude", &lat_id))) ERR(status);
     lat_buf= (float *) malloc(nlat*nlon*sizeof(float));
     if ((status = nc_get_var_float(ncid, lat_id, lat_buf))) ERR(status);
     printf("LAT dims= %d %d\n", (int) nlat, (int) nlon);

     if ((status= nc_inq_varid(ncid, "longitude", &lon_id))) ERR(status);
     lon_buf= (float *) malloc(nlat*nlon*sizeof(float));
     if ((status = nc_get_var_float(ncid, lon_id, lon_buf))) ERR(status);
     printf("LON dims= %d %d\n", (int) nlat, (int) nlon);

     if ((status= nc_inq_varid(ncid, "levels", &p_id))) ERR(status);
     p_buf= (float *) malloc(nlevel*sizeof(float));
     if ((status = nc_get_var_float(ncid, p_id, p_buf))) ERR(status);
     printf("P dims= %d\n", (int) nlevel);

     status = nc_close(ncid);
     if ( status != 0 ) ERR(status);

     nbad= 0;
     for (i= 0; i< nlat; i++) {
       for (j= 0; j< nlon; j++) {
         idx= i*nlon + j;
         lat1= lat_buf[idx];
         lat2= lat_S1_buf[idx];
         lon1= lon_buf[idx];
         lon2= lon_S1_buf[idx];
         if ( ( fabs(lat1-lat2) > 0.2 ) || ( fabs(lon) < 178 && fabs(lon1-lon2) > 0.2 )) {
           printf("SAT coords=%8.3f %8.3f  model=%8.3f %8.3f\n", lat_S1_buf[idx], lon_S1_buf[idx], lat_buf[idx], lon_buf[idx]);
           nbad++;
         }
       }
     }
     if ( nbad > 0 ) {
       printf("SAT and model coordinates are not aligned  N=%d\n", nbad);
       exit(1);
     }

     /* Freezing level height */

     h273_buf= (short *) malloc(nlat*nlon*sizeof(short));
     for (i= 0; i< nlat; i++) {
       for (j= 0; j< nlon; j++) {
         idx= i*nlon + j;
         h273_buf[idx]= -999.0;
         h273= -1.0;
         for (j1= 0; j1 < nlevel; j1++) {
           idx3= i*ns_S1*nlevel + j*nlevel + j1;
           /*printf("%2d ModelH= %9.4f  ModelP=%9.4f  ModelT=%9.4f\n", j1, h_buf[idx3], p_buf[j1], t_buf[idx3]);*/
           if ( j1 == 0 ) continue;
           idx2= i*ns_S1*nlevel + j*nlevel + (j1-1);
           if ( h_buf[idx3] < -2000 || h_buf[idx3] > 10000 ) continue;
           if ( h_buf[idx2] < -2000 || h_buf[idx2] > 10000 ) continue;
           if ( t_buf[idx2] < 273 && t_buf[idx3] >= 273 ) {
             frac= ( 273.0 - t_buf[idx2])/(t_buf[idx3] - t_buf[idx2]);
             h273= h_buf[idx2] + frac*(h_buf[idx3] - h_buf[idx2]);
             /*printf("FREEZLEV  %7.2f %7.2f ModelH=%9.4f %9.4f  ModelT=%9.4f %9.4f Frac=%9.4f H273=%9.4f\n",
               lat_buf[idx], lon_buf[idx], h_buf[idx3], h_buf[idx2], t_buf[idx3], t_buf[idx2], frac, h273);*/
             h273_buf[idx]= h273;
             break;
           }
         }
       }
     }

   }  

/*
   if ( have_model > 0 ) {
     for (i= 0; i< nlat; i++) {
       t1= time_buf[i];
       tm= *gmtime(&t1);
       sprintf(tdate,"%04d/%02d/%02d %02d:%02d:%02d",1900+tm.tm_year,1+tm.tm_mon,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);
       for (j= 0; j< nlon; j++) {
         idx= i*nlon + j;
         printf("%4d %4d %8.3f %8.3f %s  %8.3f %8.3f %8.3f  %8.3f %8.3f %6d\n", i, j, lat_buf[idx], lon_buf[idx], tdate, ts_buf[idx], t2m_buf[idx], tqv_buf[idx], hs_buf[idx], ps_buf[idx], h273_buf[idx]);
         for (k= 0; k< nlevel; k++) {
           idx1= i*nlon*nlevel + j*nlevel + k;
           printf(" %4d %8.2f %8.2f %8.2f %12.8f\n", k, p_buf[k], h_buf[idx1], t_buf[idx1], qv_buf[idx1]);
         }
       }
     } 
   }
*/


   /*---------------------------------------------------------------*/

   /* Assign DB file index for each pixel */

   for (k= 0; k< nc_out; k++) {
     tb_atms_min[k]= 1.0E6;
     tb_atms_max[k]= -1.0E6;
   }

   n1= n2= 0;
   printf("Locating number of DB files required \n");

   nsat= 0;
   for (isat= isat_start; isat <= isat_end; isat++) {

     iyyyy= yyyy_S1_buf[isat];
     imm= mm_S1_buf[isat];
     idd= dd_S1_buf[isat];
     ihh= hh_S1_buf[isat];
     imn= mn_S1_buf[isat];
     iss= ss_S1_buf[isat];
     tm.tm_year= iyyyy-1900;
     tm.tm_mon= imm-1;
     tm.tm_mday= idd;
     tm.tm_hour= ihh;
     tm.tm_min= imn;
     tm.tm_sec= iss;
     t1= timegm(&tm);
     if ( t1 <= 0 ) {
       for (jsat=0; jsat < ns_S1; jsat++) {
         idx_out= nsat*ns_S1 + jsat;
         epc_qual_NS_buf[idx_out]= -3;
         epc_qual_MS_buf[idx_out]= -3;
       }
       continue;  /* bad time */
     }
     if ( t1 > t0 ) t0= t1;  /* oldest time in the files */
     tm= *gmtime(&t1);
     sprintf(tdate,"%04d/%02d/%02d %02d:%02d:%02d",1900+tm.tm_year,1+tm.tm_mon,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);

     for (jsat=0; jsat < ns_S1; jsat++) {
       idx= isat*ns_S1 + jsat;

       /*printf("%4d %4d %s\n", isat, jsat, tdate);*/

       glat1= lat_S1_buf[idx];
       glon1= lon_S1_buf[idx];
       if ( fabs(glat1) > 90 || fabs(glon1) > 180 ) {
         idx_out= nsat*ns_S1 + jsat;
         epc_qual_NS_buf[idx_out]= -4;
         epc_qual_MS_buf[idx_out]= -4;
         continue;
       }

       slat= sclat_S1_buf[isat];
       slon= sclon_S1_buf[isat];
       if ( fabs(slat) > 90 || fabs(slon) > 180 ) {
         for (j=0; j< ns_S1; j++) {
           idx_out= nsat*ns_S1 + j;
           epc_qual_NS_buf[idx_out]= -5;
           epc_qual_MS_buf[idx_out]= -5;
         }
         continue;
       }
       /*printf("%4d %3d %8d %8.3f %8.3f %8.3f %8.3f \n", isat,jsat, idx, glat1, glon1, slat, slon);*/

       nt= nbad= 0;
       for (k= 0; k< nc_out; k++) {
         idx1= isat*ns_S1*nc_out + jsat*nc_out + k;
         if ( Tc_buf[idx1] > 20 && Tc_buf[idx1] < 400 ) tb_atms[k]= Tc_buf[idx1];
         else continue;
         if ( tb_compare[k] == 1 ) {
           tb_reg[nt]= tb_atms[k];
           nt++;
         }
         if ( tb_atms[k] < tb_atms_min[k] ) tb_atms_min[k]= tb_atms[k];
         if ( tb_atms[k] > tb_atms_max[k] ) tb_atms_max[k]= tb_atms[k];
       }
       if ( nt != NTBREG ) {
         for (j=0; j< ns_S1; j++) {
           idx_out= nsat*ns_S1 + j;
           epc_qual_NS_buf[idx_out]= -1;
           epc_qual_MS_buf[idx_out]= -1;
         }
         continue;
       }

       /*if ( glat1 < -32 || glat1 > -31 || glon1 < -53 || glon1 > -52 ) continue;*/
       /*if ( tb_gmi[6] > 200 && tb_gmi[8] > 180 ) continue;*/

       /* If available, use interpolated MERRA2 fields */

       if ( have_model > 0 ) {
         ts_MERRA2= ts_buf[idx];
         t2m_MERRA2= t2m_buf[idx];
         tqv_MERRA2= tqv_buf[idx];
         h273_MERRA2= h273_buf[idx];
       }


       /*----------------------------------------------------*/

       tmp5[0]= 1.0;
       kt= 0;
       for (k= 0; k< NTBREG; k++) {
         tmp5[kt+1]= tb_reg[k];
         kt++;
         for (k1= k; k1< NTBREG; k1++) {
           tmp5[kt+1]= tb_reg[k]*tb_reg[k1];
           kt++;
         }
       }
       for (k= 1; k< NTBREG; k++) {
         tmp5[kt+1]= (tb_reg[k]-tb_reg[k-1])/(tb_reg[k]+tb_reg[k-1]);
         kt++;
       }
       if ( kt != NREG ) {
         printf("nterms=%d  NREG=%d\n", kt, NREG);
         exit(1);
       }

       /* emissivity PC */
       for (k= 0; k< NEM; k++) {
         psum= 0.0;
         for (k2= 0; k2< NREG+1; k2++) psum+= b_all[k2][k]*tmp5[k2];
         pc_emis[k]= psum;
         /*printf("%10.6f ", pc_emis[k]);*/
       }
       /*printf("\n");*/

       /* emissivites from PC */
       for (k= 0; k< NEM; k++) {
         psum= 0.0;
         for (k2= 0; k2< NEM; k2++) psum+= u[k+1][k2+1] * pc_emis[k2];
         emis[k]= psum + ave_emis[k];
         /*printf("%10.6f ", emis[k]);*/
       }
       /*printf("\n");*/

       /* Save for output file */
       for (k=0; k < NEM; k++) {
         idx_out= nsat*ns_S1*NEM + jsat*NEM + k;
         emis_buf[idx_out]= emis[k];
         pc_emis_buf[idx_out]= pc_emis[k];
       }


       for (k= 0; k< NEM_USE; k++) {
         found= 0;
         for (j= 0; j < NPCHIST; j++) {
           pc_lo= pc_range[k][j];
           pc_hi= pc_range[k][j+1];
           /*printf("%d   pc=%f  range=%f %f   j=%d \n", k, pc_emis[k], pc_lo, pc_hi, j);*/
           if ( pc_emis[k] > pc_lo && pc_emis[k] <= pc_hi ) {found=1; break;}
         }
         if ( found != 1 ) {
           printf("DB index %d not located  pc=%f  range=%f %f\n", k, pc_emis[k], pc_lo, pc_hi);
           exit(1);
         }
         idxn[k]= j;
       }
       idx_db= 0;
       for (k= 0; k< NEM_USE; k++)
         idx_db+= idxn[NEM_USE-1-k]*pow2[k];

       /*for (k= 0; k< NEM_USE; k++) {
         j= idxn[k];
         printf("%1d ", j);
         printf("%12.6f %12.6f  %12.6f\n", pc_range[k][j], pc_range[k][j+1], pc_emis[k]);
       }
       for (k= 0; k< NEM_USE; k++) printf("%d ", idxn[k]);
       printf("  index=%d\n", idx_db);*/

       if ( idx_db < 0 || idx_db >= ndb_files ) {
         printf("Exceeded DB index=%d\n", idx_db);
         exit(1);
       }

       pixel_index[idx]= idx_db;
       if ( n_index[idx_db] == 0 ) n1++;
       n_index[idx_db]++;
       if ( have_model > 0 ) {
         if ( ts_MERRA2 < ts_MERRA2_min ) ts_MERRA2_min= ts_MERRA2;
         if ( ts_MERRA2 > ts_MERRA2_max ) ts_MERRA2_max= ts_MERRA2;
         if ( t2m_MERRA2 < t2m_MERRA2_min ) t2m_MERRA2_min= t2m_MERRA2;
         if ( t2m_MERRA2 > t2m_MERRA2_max ) t2m_MERRA2_max= t2m_MERRA2;
         if ( tqv_MERRA2 < tqv_MERRA2_min ) tqv_MERRA2_min= tqv_MERRA2;
         if ( tqv_MERRA2 > tqv_MERRA2_max ) tqv_MERRA2_max= tqv_MERRA2;
         if ( t2m_MERRA2 < 273 ) {
           if ( ncold_index[idx_db] == 0 ) n2++;
           ncold_index[idx_db]++;
         }
       }
       for (k= 0; k< NEM; k++) {
         if ( pc_emis[k] < epc_min[k] ) epc_min[k]= pc_emis[k];
         if ( pc_emis[k] > epc_max[k] ) epc_max[k]= pc_emis[k];
       }



     }  /* next radiometer 1C sample */

     nsat++;

   }  /* next radiometer 1C line */

   printf("Number of DB files required=%d\n", n1);
   if ( n1 == 0 ) exit(1);

/*
   nsat= 0;
   for (isat= isat_start; isat <= isat_end; isat++) {
     for (jsat=0; jsat < ns_S1; jsat++) {
       idx1= isat*ns_S1 + jsat;
       idx_out= nsat*ns_S1 + jsat;
       if ( pixel_index[idx1] < 0 ) {
         printf("DB_INDEX_NOT_DETERMINED  Occured at %6d %6d %10.3f %10.3f  TB=", isat, jsat, lat_S1_buf[idx1], lon_S1_buf[idx1]);
         for (k=0; k< nc_out; k++) {
           idx3= isat*ns_S1*nc_out + jsat*nc_out + k;
           printf("%7.2f ", Tc_buf[idx3]);
         }
         printf("\n");
       }
     }
     nsat++;
   }
*/

   if ( have_model > 0 ) {
     printf("Number of DB files with T2m < 273K= %d\n", n2);
     printf("Ts  range= %8.3f %8.3f\n", ts_MERRA2_min, ts_MERRA2_max);
     printf("T2m range= %8.3f %8.3f\n", t2m_MERRA2_min, t2m_MERRA2_max);
     printf("TQV range= %8.3f %8.3f\n", tqv_MERRA2_min, tqv_MERRA2_max);
   }
   for (k= 0; k< NEM; k++) {
     printf("EPC %2d range= %12.6f %12.6f\n", k, epc_min[k], epc_max[k]);
   }



   /*---------------------------------------------------*/
   /*--- Process all data for each DB file only once ---*/

   for (idx_db= 0; idx_db< ndb_files; idx_db++) {

     if ( n_index[idx_db] == 0 ) continue;
     ndb= 0;

     for (i= 0; i< 6; i++) {
       for (j= 0; j< 2; j++) nrain_db[i][j]= 0;
     }
     idx_db1= ndb_files+1;
     idx_db2= -1;
     nrec= ndb_expand= 0;

     for (n= 0; n< 2*MAX_DB_EXPAND; n++) {
       sign= ( n % 2 == 0 ) ? 1 : -1;
       idx_db_expand= idx_db + sign*((n+1)/2) ;
       if ( idx_db_expand < 0 || idx_db_expand >= ndb_files ) continue;
     
       /*--- Read the number of rain events for cold and warm T2m conditions ---*/
       /* first six=  Ntotal then Nrain exceeding 1 5 10 20 50 mm/hr   when T2m < 273K */
       /* second six= Ntotal then Nrain exceeding 1 5 10 20 50 mm/hr   when T2m >= 273K */

       sprintf(svar,"%s/db_%05d.bin.nrain.txt", dbdir, idx_db_expand);
       if ((fdb_nrain = fopen(svar,"r")) == 0) {
         /*printf("Unable to locate database precip CDF file %s\n", svar);*/
         continue;
       }
       if ( idx_db_expand < idx_db1 ) idx_db1= idx_db_expand;
       if ( idx_db_expand > idx_db2 ) idx_db2= idx_db_expand;
       ndb_expand++;
       
       /*printf("%s  COLD= ", dbfile);*/
       for (i= 0; i< 6; i++) {
         fscanf(fdb_nrain, "%d ", &j);
         nrain_db[i][0]+= j;
         /*printf("%8d ", j);*/
       }
       /*printf("  WARM= ");*/
       for (i= 0; i< 6; i++) {
         fscanf(fdb_nrain, "%d ", &j);
         nrain_db[i][1]+= j;
         /*printf("%8d ", j);*/
       }
       fclose(fdb_nrain);

       nrec= nrain_db[0][0] + nrain_db[0][1];;
       if ( nrec >= DB_MAXREC ) break;
       if ( n > 0 && nrec >= DB_MINREC ) break;
     }
     if ( ndb_expand == 0 ) {
       printf("DB ERROR idx_db=%6d   N_db_expanded=%6d\n", idx_db, ndb_expand);
       continue;
     }
     else {
       printf("idx_db=%6d   N_db_expanded=%6d  idx_db range=%6d %6d  ", idx_db, idx_db1, idx_db2);

       frac0= ( nrain_db[0][0] > 0 ) ? 1.0*nrain_db[1][0]/nrain_db[0][0] : 0.0;  /* cold T2m */
       frac1= ( nrain_db[0][1] > 0 ) ? 1.0*nrain_db[1][1]/nrain_db[0][1] : 0.0;  /* warm T2m */
       for (i= 0; i< 6; i++) {
         for (j= 0; j< 2; j++) printf("%8d ", nrain_db[i][j]);
       }
       printf("  Frac1=%9.5f %9.5f", frac0, frac1);
       printf("\n");
     }

     /*--- If the DB file has insufficient records, set epc_qual flags to -2 */

     if ( nrec < DB_MINREC ) {
       printf("DB_ERROR idx_db=%6d  N insufficient=%d\n", idx_db, nrec);
       nsat= 0;
       for (isat= isat_start; isat <= isat_end; isat++) {
         for (jsat=0; jsat < ns_S1; jsat++) {
           idx1= isat*ns_S1 + jsat;
           idx_out= nsat*ns_S1 + jsat;
           if ( pixel_index[idx1] == idx_db ) {
             printf("  DB_ERROR  Occured for %4d %4d %8.3f %8.3f  TB=", isat, jsat, lat_S1_buf[idx1], lon_S1_buf[idx1]);
             for (k=0; k< nc_out; k++) {
               idx3= isat*ns_S1*nc_out + jsat*nc_out + k;
               printf("%7.2f ", Tc_buf[idx3]);
             }
             printf("\n");
             epc_qual_NS_buf[idx_out]= -2;
             epc_qual_MS_buf[idx_out]= -2;
           }
         }
         nsat++;
       }
       npix+= n_index[idx_db];
       continue;
     }

     /*--- If the DB file contains very little rain for cold and warm T2m, set to zero-rain */

     else if ( nrain_db[2][0] == 0 && nrain_db[2][1] == 0 && frac0 < DB_RAINFRAC && frac1 < DB_RAINFRAC ) {
       printf("ZERO_RAIN idx_db=%6d  Insufficient fraction of DB entries exceeding 1 mm/hr= %9.5f %9.5f\n", idx_db, frac0, frac1);

       nsat= 0;
       for (isat= isat_start; isat <= isat_end; isat++) {
         for (jsat=0; jsat < ns_S1; jsat++) {
           idx1= isat*ns_S1 + jsat;
           idx_out= nsat*ns_S1 + jsat;
           if ( pixel_index[idx1] != idx_db ) continue;
           printf("  ZERO RAIN  Occured for %4d %4d %8.3f %8.3f  TB=", isat, jsat, lat_S1_buf[idx1], lon_S1_buf[idx1]);
           for (k=0; k< nc_out; k++) {
             idx3= isat*ns_S1*nc_out + jsat*nc_out + k;
             printf("%7.2f ", Tc_buf[idx3]);
           }
           printf("\n");
           precip1_NS_buf[idx_out]= 0.0;
           precip2_NS_buf[idx_out]= 0.0;
           precip1_MS_buf[idx_out]= 0.0;
           precip2_MS_buf[idx_out]= 0.0;
           precip1_prob_NS_buf[idx_out]= 0.0;
           precip2_prob_NS_buf[idx_out]= 0.0;
           precip1_prob_MS_buf[idx_out]= 0.0;
           precip2_prob_MS_buf[idx_out]= 0.0;
           epc_qual_NS_buf[idx_out]= nrec;
           epc_qual_MS_buf[idx_out]= nrec;
         }
         nsat++;
       }
       npix+= n_index[idx_db];
       continue;
     }

     maxrec= ( nrec < DB_MAXREC ) ? nrec : DB_MAXREC;

     /*-----------------------------------------------------------------------*/
     /* Read the index database file one time */

     nrec= nrec1= nrec2= nrain= 0;

     /*--- Open DB files ---*/
     for (idx_db_expand= idx_db1; idx_db_expand<= idx_db2; idx_db_expand++) {

       sprintf(dbfile,"%s/db_%05d.bin", dbdir, idx_db_expand);
       if ((fdb2 = fopen(dbfile,"r")) == 0) {
         /*printf("DB_ERROR  Unable to locate %s\n", dbfile);*/
         continue;
       }

       while ( fread(&prof,bytes,1,fdb2) == 1 ) {
         nrec++;

         /*if ( have_model > 0 ) {
           if ( prof.t2m < t2m_MERRA2_min - MAX_T2M_DIFF ) continue;
           if ( prof.t2m > t2m_MERRA2_max + MAX_T2M_DIFF ) continue;
         }*/

         if ( prof.j_NS < 0 || prof.j_NS >= 49 ) continue;  /* NS positions 0-48 */
         /* NS positions 12-36 correspond to MS positions 0-24 */
         /*kuka= ( prof.j_NS < 12 || prof.j_NS > 36 ) ? 0 : 1;*/
         kuka= ( prof.j_NS <= 12 || prof.j_NS >= 36 ) ? 0 : 1;

         /*--- Reject bad entries ----*/
         nbad= 0;
         for (i= 0; i< NCHAN; i++) {
           if ( prof.tb[i] < -20.0 && prof.tb[i] > -400 ) prof.tb[i]*= -1.0; 
           if ( tb_compare[i] == 1 && ( prof.tb[i] < 20.0 || prof.tb[i] > 400 )) nbad++;
         }
         /*if ( prof.sfc_class < 0 ) nbad++;*/
         if ( prof.yyyy < 1900 || prof.yyyy > 2100 ) nbad++;
         if ( prof.mm < 1 || prof.mm > 12 ) nbad++;
         if ( prof.dd < 1 || prof.dd > 31 ) nbad++;
         if ( prof.hh < 0 || prof.hh > 23 ) nbad++;
         if ( prof.mn < 0 || prof.mn > 59 ) nbad++;
         if ( prof.ss < 0 || prof.ss > 59 ) nbad++;
         if ( fabs(prof.glat) > 90 || fabs(prof.glon) > 180 ) nbad++;
         if ( prof.nku10 == 0 && prof.nku15 == 0 && prof.nku20 == 0 && prof.nku25 == 0 ) nbad++;
         if ( prof.precip_nsfc_NS < 0.0 || prof.precip_NS_cmb < 0.0 ) nbad++;
         if ( kuka == 1 ) {
           if ( prof.nka10 == 0 && prof.nka15 == 0 && prof.nka20 == 0 && prof.nka25 == 0 ) nbad++;
           if ( prof.precip_nsfc_MS < 0.0 || prof.precip_MS_cmb < 0.0 ) nbad++;
         }
         if ( prof.ts < 0 || prof.ts > 400 ) nbad++;
         if ( prof.t2m < 0 || prof.t2m > 400 ) nbad++;
         if ( prof.tqv < 0 || prof.tqv > 100 ) nbad++;

         if ( nbad > 0 ) {
           printf("BAD_DPR %2d %6d %5d %3d %8.3f %8.3f ", kuka, prof.rev, prof.i_NS, prof.j_NS, prof.glat, prof.glon);
           printf("%04d/%02d/%02d %02d:%02d:%02d ", prof.yyyy, prof.mm, prof.dd, prof.hh, prof.mn, prof.ss);
           for (k= 0; k< 3; k++) printf("%8.3f ", prof.pc_emis[k]);
           for (k= 0; k< NCHAN; k++) printf("%7.2f ", prof.tb[k]);
           printf("%3d %6.2f %6.2f %6.2f ", prof.sfc_class, prof.ts, prof.t2m, prof.tqv);
           printf("%4d %4d %4d %4d ", prof.nku10, prof.nku15, prof.nku20, prof.nku25);
           if ( kuka == 1 ) printf("%4d %4d %4d %4d ", prof.nka10, prof.nka15, prof.nka20, prof.nka25);
           printf("%8.2f %8.2f ", prof.precip_nsfc_NS, prof.precip_NS_cmb);
           if ( kuka == 1 ) printf("%8.2f %8.2f ", prof.precip_nsfc_MS, prof.precip_MS_cmb);
           printf("\n");
           continue;
         }
         /*---------------------------*/

         if ( prof.precip_NS_cmb > 1 ) nrain++;

         prof_idx[nrec1]= prof;
         nrec1++;
         if ( kuka == 1 ) nrec2++;

         if ( nrec1 == maxrec ) {
           /*printf("Reached max records=%d   Ku=%d KuKa=%d\n", nrec, nrec1, nrec2);*/
           break;
         }
       }
       fclose(fdb2);
       printf("%s  NRead=%d  Nrec=%d %d Nrain=%d    Npix using this DB file=%d\n", dbfile, nrec, nrec1, nrec2, nrain, n_index[idx_db_expand]);
       if ( nrec1 == maxrec ) {
         printf("Reached max records=%d   Ku=%d KuKa=%d\n", nrec, nrec1, nrec2);
         break;
       }
     }

     if ( nrec1 == 0 || nrec2 == 0 ) continue;

     /*-----------------------------------------------------------------------*/

     nsat= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {

       t1= secs_S1_buf[isat];
       tm= *gmtime(&t1);
       sprintf(tdate,"%04d/%02d/%02d %02d:%02d:%02d",1900+tm.tm_year,1+tm.tm_mon,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);

       for (jsat=0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx_out= nsat*ns_S1 + jsat;
           
         if ( pixel_index[idx] != idx_db ) continue;
         ndb++;
         pcent_db= 100.0*ndb/n_index[idx_db];

         printf("%d %d\n", idx, idx_out);

         /* Check if this radiometer 1C pixel lies within the DPR NS swath */
         /* Accumulate CCFADS only when within the interesecting Ku and/or Ka-band swath area */
         if ( have_dpr > 0 ) {
           inside_NS= ( dpr_NS_i[idx] >= 0 && dpr_NS_j[idx] >= 0 ) ? 1 : -1;
           inside_MS= ( dpr_NS_i[idx] >= 0 && dpr_NS_j[idx] >= 12 && dpr_NS_j[idx] <= 36 ) ? 1 : -1;
         }

         glat1= lat_S1_buf[idx];
         glon1= lon_S1_buf[idx];
         if ( fabs(glat1) > 90 || fabs(glon1) > 180 ) continue;
         incidenceAngle= fabs(incidenceAngle_S1_buf[idx]);

         /*--- elevation from external database */
         ilat= ilon= -1;
         gelev= -9999;
         for (k= 1; k< ny_elev; k++) {
           if ( glat1 >= elat_buf[k-1] && glat1 <= elat_buf[k] ) {ilat= k; break; }
         }
         for (k= 1; k< nx_elev; k++) {
           if ( glon1 >= elon_buf[k-1] && glon1 <= elon_buf[k] ) {ilon= k; break; }
         }
         if ( ilat >= 0 && ilat < ny_elev && ilon >= 0 && ilon < nx_elev)
           gelev= elev_buf[ilat*nx_elev + ilon];
         /*printf("ELEV %8.3f %8.3f %8.3f %8.3f %6d %6d %d\n", glat1, glon1, elat_buf[ilat], elon_buf[ilon], ilat, ilon, gelev);*/


         if ( have_model > 0 ) {
           ts_MERRA2= ts_buf[idx];
           t2m_MERRA2= t2m_buf[idx];
           tqv_MERRA2= tqv_buf[idx];
           h273_MERRA2= h273_buf[idx];
         }
         else {
           ts_MERRA2= t2m_MERRA2= tqv_MERRA2= -1.0;
           h273_MERRA2= -9999;
         }

         if ( have_gprof > 0 ) {
           idx1= isat*ns_SG + jsat;
           sfc_class= sfc_GPROF_S1_buf[idx];
           precip_GPROF= surfacePrecipitation_GPROF_S1_buf[idx];
           prob_precip_GPROF= probabilityOfPrecip_GPROF_S1_buf[idx];
           tqv_GPROF= totalColumnWaterVaporIndex_GPROF_S1_buf[idx];
           t2m_GPROF= temp2mIndex_GPROF_S1_buf[idx];
           frozen_precip_GPROF= frozenPrecipitation_GPROF_S1_buf[idx];
           if ( precip_GPROF < 0.0 ) precip_GPROF= -1.0;
           if ( frozen_precip_GPROF < 0.0 ) frozen_precip_GPROF= -1.0;
           if ( prob_precip_GPROF < 0.0 || prob_precip_GPROF > 100.0 ) prob_precip_GPROF= -1.0;
         }
         else {
           sfc_class= tqv_GPROF= t2m_GPROF= -1;
           precip_GPROF= prob_precip_GPROF= frozen_precip_GPROF= -1.0;
         }


         nbad= 0;
         for (k= 0; k< nc_out; k++) {
           idx1= isat*ns_S1*nc_out + jsat*nc_out + k;
           /*tb_atms[k]= Tc_buf[idx1];*/
           tb_atms[k]= fabs(Tc_buf[idx1]);
           if ( tb_compare[k] == 1 && ( tb_atms[k] < 20 || tb_atms[k] > 400 )) nbad++;
         }
         if ( nbad != 0 ) {
           printf("%4d %4d  bad TB\n", isat, jsat);
           npix++;
           continue;
         }

         for (k=0; k < NEM; k++) {
           idx3= nsat*ns_S1*NEM + jsat*NEM + k;
           emis[k]= emis_buf[idx3];
           pc_emis[k]= pc_emis_buf[idx3];
         }

         /* START LOOKUP */

         printf("\n\n");
         x= -1.0;
         ix= -1;
         printf("OBS %6d %3d %6d %s %4d %3d %7.2f %7.2f ", idx_db, isatname, irev, tdate, isat, jsat, glat1, glon1);
         for (k= 0; k< 3; k++) printf("%8.3f ", pc_emis[k]);
         for (k= 0; k< 9; k++) printf("%6.2f ", tb_atms[k]);
         printf("%3d ", sfc_class);
         printf("%7.2f %7.2f %7.2f %6d ", ts_MERRA2, t2m_MERRA2, tqv_MERRA2, h273_MERRA2);
         printf("%4d %4d ", ix, ix);
         printf("%7.2f %7.2f %7.2f ", x,x,x);
         printf("%6d ", gelev);
         e0= ( emis[0] > 0.0 && emis[0] < 2.0 ) ? emis[0] : -1.0;
         e1= ( emis[1] > 0.0 && emis[1] < 2.0 ) ? emis[1] : -1.0;
         printf("%6.3f %6.3f ", e0, e1);
         printf("%6.3f %6.3f ", x, x );
         printf("\n");

         rmsd_tb_min= rmsd_pc_min= 1.0E9;
         rmsd_tb_max= rmsd_pc_max= 0.0;
         nrec2= 0;
         nmodel= 0;
         for (i= 0; i< nrec1; i++) {

           /*-------------------------------------------------------------------------------*/

           /* Over high terrain discard elevation using external database */
           if ( gelev > 500 ) {

             /*ilat= ilon= -1;
             gelev1= -9999;
             for (k= 1; k< ny_elev; k++) {
               if ( prof_idx[i].glat >= elat_buf[k-1] && prof_idx[i].glat <= elat_buf[k] ) {ilat= k; break; }
             }
             for (k= 1; k< nx_elev; k++) {
               if ( prof_idx[i].glon >= elon_buf[k-1] && prof_idx[i].glon <= elon_buf[k] ) {ilon= k; break; }
             }
             if ( ilat >= 0 && ilat < ny_elev && ilon >= 0 && ilon < nx_elev)
               gelev1= elev_buf[ilat*nx_elev + ilon];
             printf("ELEV  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f   %6d %6d %6d\n", glat1, glon1, prof_idx[i].glat, prof_idx[i].glon, elat_buf[ilat], elon_buf[ilon], gelev, prof_idx[i].elev, gelev1);*/
             
             if ( prof_idx[i].elev < gelev ) continue;
           }

           /* If TB are too far apart disregard */
           if ( MAX_TB_RMSD > 0 ) {
             rmsd_tb= 0.0;
             ntb_compare= 0;
             for (k= 0; k< NCHAN; k++) {
               if  (tb_atms[k] < 20 || tb_atms[k] > 400 ) continue;
               if ( prof_idx[i].tb[k] < 20 || prof_idx[i].tb[k] > 400 ) continue;
               x= (tb_atms[k] - prof_idx[i].tb[k]);
               rmsd_tb+= x*x;
               ntb_compare++;
             }
             if ( ntb_compare > 0 ) {
               rmsd_tb= sqrt(rmsd_tb/ntb_compare);
               if ( rmsd_tb > MAX_TB_RMSD ) {
                 /*for (k= 0; k< 7; k++) printf("%6.2f ", tb_atms[k]);
                 for (k= 0; k< 7; k++) printf("%6.2f ", prof_idx[i].tb[k]);
                 printf(" N=%3d %8.3f \n", ntb_compare, rmsd_tb);*/
                 continue;
               }
             }
           }

           /* If T2m is warm, discard entries where T2m is cold, and vice-versa */
           if ( have_model > 0 ) {
             if ( fabs(prof_idx[i].t2m - t2m_MERRA2) > MAX_T2M_DIFF ) {

/*
               printf("OBS MODEL=%8.3f %8.3f %8.3f TB=", ts_MERRA2, t2m_MERRA2, tqv_MERRA2);
               for (k= 0; k< 7; k++) printf("%6.2f ", tb_atms[k]);
               printf("  DB MODEL=%8.3f %8.3f %8.3f TB=", prof_idx[i].ts, prof_idx[i].t2m, prof_idx[i].tqv);
               for (k= 0; k< 7; k++) printf("%6.2f ", prof_idx[i].tb[k]);

               for (k= 0; k< 4; k++) printf("%8.3f ", pc_emis[k]);
               for (k= 0; k< NTBREG; k++) tb_reg[k]= prof_idx[i].tb[k];
               tmp5[0]= 1.0;
               kt= 0;
               for (k= 0; k< NTBREG; k++) {
                 tmp5[kt+1]= tb_reg[k];
                 kt++;
                 for (k1= k; k1< NTBREG; k1++) {
                   tmp5[kt+1]= tb_reg[k]*tb_reg[k1];
                   kt++;
                 }
               }
               for (k= 1; k< NTBREG; k++) {
                 tmp5[kt+1]= (tb_reg[k]-tb_reg[k-1])/(tb_reg[k]+tb_reg[k-1]);
                 kt++;
               }
               if ( kt != NREG ) {
                 printf("nterms=%d  NREG=%d\n", kt, NREG);
                 exit(1);
               }

               for (k= 0; k< 4; k++) {
                 psum= 0.0;
                 for (k2= 0; k2< NREG+1; k2++) psum+= b_all[k2][k]*tmp5[k2];
                 printf("%8.3f ", psum);
               }
               printf("\n");
*/

               nmodel++;
               continue;
             }
           }

           /*----- discard candidates where the incidence angle difference is too large ----*/
           diff= incidenceAngle - fabs(prof_idx[i].inc_S1);
           if ( fabs(diff) > DB_MAX_INC_DIFF ) {
             /*printf("SAT=%3d %8.3f   DB=%3d %8.3f  inc angle diff=%8.3f\n", jsat, incidenceAngle, prof_idx[i].j_S1, prof_idx[i].inc_S1, diff);*/
             /*exit(1);*/
             continue;
           }

           /*----- discard candidates from the same orbit rev ----*/
           if ( irev == prof_idx[i].rev ) continue;
           /*-------------------------------------------------------------------------------*/

           rmsd_tb= 0.0;
           ntb_compare= 0;
           for (k= 0; k< NCHAN; k++) {
             /*if  (tb_atms[k] < 20 || tb_atms[k] > 400 ) continue;
             if ( prof_idx[i].tb[k] < 20 || prof_idx[i].tb[k] > 400 ) continue;
             if ( std_tb[k] <= 0.0 ) continue;*/
             x= (tb_atms[k] - prof_idx[i].tb[k])/std_tb[k];
             rmsd_tb+= x*x;
             ntb_compare++;
           }
           if ( ntb_compare > 0 ) {
             rmsd_tb= sqrt(rmsd_tb/ntb_compare);
           }

           rmsd_pc= 0.0;
           nem_compare= 0;
           for (k= 0; k< NEM; k++) {
             if ( em_compare[k] != 1 ) continue;
             if ( std_pc[k] <= 0.0 ) continue;
             x= (pc_emis[k] - prof_idx[i].pc_emis[k])/std_pc[k];
             rmsd_pc+= x*x;
             nem_compare++;
           }
           if ( nem_compare == 0 ) continue;
           /*if ( rmsd_pc < 1.0E-12 ) continue;*/
           rmsd_pc= sqrt(rmsd_pc/nem_compare);

           wt= rmsd_pc;
           /*wt= sqrt(0.7*rmsd_tb*rmsd_tb + 0.3*rmsd_pc*rmsd_pc);*/
           /*wt= sqrt(0.3*rmsd_tb*rmsd_tb + 0.7*rmsd_pc*rmsd_pc);*/

           if ( wt >= INT_MAX/1E6 ) {
             printf("Exceeded INT_MAX=%d  rmsd_pc=%f rmsd_tb=%f wt=%f\n", INT_MAX, rmsd_pc, rmsd_tb, wt);
             exit(1);
           }
           tmp1[nrec2]= 1E6*wt;
           tmp2[nrec2]= i;
           tmp3[i]= wt;
           tmp7[nrec2]= wt;
           if ( rmsd_pc > rmsd_pc_max ) rmsd_pc_max= rmsd_pc;
           if ( rmsd_pc < rmsd_pc_min ) rmsd_pc_min= rmsd_pc;
           if ( rmsd_tb > rmsd_tb_max ) rmsd_tb_max= rmsd_tb;
           if ( rmsd_tb < rmsd_tb_min ) rmsd_tb_min= rmsd_tb;
           nrec2++;

         }
         /*printf("PC rmsd min max= %10.5f %10.5f ", rmsd_pc_min, rmsd_pc_max);
         printf("TB rmsd min max= %10.5f %10.5f\n", rmsd_tb_min, rmsd_tb_max);*/

         heapsort2(tmp1, tmp2, nrec2);
         /*for (i= 0; i< nrec2; i++) {
           k2= tmp2[i];
           if ( i < 10 ) printf("%6d %6d   %d %f   %d %f\n", i, k2, tmp1[i], tmp7[i], tmp1[k2], tmp7[k2]);
         }*/


         /*--- First the NS database then the MS database ---*/

         for (idb= 0; idb< 2; idb++) {

           if ( idb == 0 ) printf("NS database  N=%d %8.3f\n", nrec2, pcent_db);
           else printf("MS database  N=%d %8.3f\n", nrec2, pcent_db);

           zmax= 0.0;
           n1= n2= 0;
           precip1_sum= precip2_sum= ts_sum= t2m_sum= tqv_sum= h273_sum= ku_ztop_sum= ka_ztop_sum= 0.0;
           wt_precip_sum= wt_ts_sum= wt_t2m_sum= wt_tqv_sum= wt_h273_sum= wt_ku_ztop_sum= wt_ka_ztop_sum= 0.0;
           n_precip1_sum= n_precip2_sum= n_precip_all_sum= 0;
           for (iclass= 0; iclass< NCLASS; iclass++) class[iclass]= 0;
           for (iprecip= 0; iprecip< NPRECIP; iprecip++) hist_precip[iprecip]= 0;
           for (k= 0; k< NLEV_DPR; k++) {
             z_ku[k]= z_ka[k]= 0.0;
             nz_ku[k]= nz_ka[k]= 0;
           }
           for (k= 0; k< NLEV_PRECIP; k++) {
             precip1_prof_sum[k]= precip2_prof_sum[k]= 0.0;
             n_precip1_prof_sum[k]= n_precip2_prof_sum[k]= 0;
           }
           for (k= 0; k< 3; k++) type_precip[k]= 0.0;
           for (k= 0; k< 5; k++) shallow_rain[k]= 0.0;

           wt0= -1.0;
           nfound= 0;
           nku20= nku25= 0;

           for (i= 0; i< nrec2; i++) {
             k2= tmp2[i];
             /*wt= tmp7[k2];*/
             wt= tmp3[k2];

             /*kuka= ( prof_idx[k2].j_NS < 12 || prof_idx[k2].j_NS > 36 ) ? 0 : 1;*/
             kuka= ( prof_idx[k2].j_NS <= 12 || prof_idx[k2].j_NS >= 36 ) ? 0 : 1;
             if ( idb == 1 && kuka == 0 ) continue;

             if ( wt0 < 0.0 ) {
               wt0= wt;
               wt2= 1.0;
             }
             else {
               wt2= exp(-0.5*(wt/wt0)*(wt/wt0));
             }

             /* Freezing level height */
             h273= -9999;
             for (j1= 0; j1 < 42; j1++) {
               /*printf("%2d ModelH= %9.4f  ModelP=%9.4f  ModelT=%9.4f\n", j1, prof_idx[k2].h_prof[j1], 0.1*prof_idx[k2].p_prof[j1], prof_idx[k2].t_prof[j1]);*/
               if ( j1 == 0 ) continue;
               if ( prof_idx[k2].h_prof[j1] < -2000 || prof_idx[k2].h_prof[j1] > 10000 ) continue;
               if ( prof_idx[k2].h_prof[j1-1] < -2000 || prof_idx[k2].h_prof[j1-1] > 10000 ) continue;
               if ( prof_idx[k2].t_prof[j1-1] < 273 && prof_idx[k2].t_prof[j1] >= 273 ) {
                 frac= ( 273.0 - prof_idx[k2].t_prof[j1-1])/(prof_idx[k2].t_prof[j1] - prof_idx[k2].t_prof[j1-1]);
                 h273= prof_idx[k2].h_prof[j1-1] + frac*(prof_idx[k2].h_prof[j1] - prof_idx[k2].h_prof[j1-1]);
                 /*printf("  FREEZLEV  ModelH=%9.4f %9.4f  ModelT=%9.4f %9.4f Frac=%9.4f H273=%9.4f\n",
                   prof_idx[k2].h_prof[j1], prof_idx[k2].h_prof[j1-1], prof_idx[k2].t_prof[j1], prof_idx[k2].t_prof[j1-1], frac, h273);*/
                 break;
               }
             }
             /*if ( h273 < -2000 || h273 > 9000 ) {
               printf("Freezing level height out of bounds %f\n", h273);
               exit(1);
             }*/
             elev= prof_idx[k2].elev;

             if ( idb == 1 ) {
               precip1= prof_idx[k2].precip_MS_cmb;
               precip2= prof_idx[k2].precip_nsfc_MS;
               for (k= 0; k< NLEV_PRECIP; k++) {
                 precip1_prof[k]= prof_idx[k2].precip_prof_MS_cmb[k];
                 precip2_prof[k]= prof_idx[k2].precip_prof_MS[k];
               }
             }
             else {
               precip1= prof_idx[k2].precip_NS_cmb;
               precip2= prof_idx[k2].precip_nsfc_NS;
               for (k= 0; k< NLEV_PRECIP; k++) {
                 precip1_prof[k]= prof_idx[k2].precip_prof_NS_cmb[k];
                 precip2_prof[k]= prof_idx[k2].precip_prof_NS[k];
               }
             }
             if ( precip1 < 0.0 || precip2 < 0.0 ) {
               printf("Rainrates out of bounds %f %f\n", precip1, precip2);
               exit(1);
             }

             precip1_sum+= wt2*precip1;
             precip2_sum+= wt2*precip2;
             wt_precip_sum+= wt2;
             if ( precip1 > 0.0 ) n_precip1_sum++;
             if ( precip2 > 0.0 ) n_precip2_sum++;
             n_precip_all_sum++;
             for (k= 0; k< NLEV_PRECIP; k++) {
               if ( precip1_prof[k] >= 0 ) precip1_prof_sum[k]+= wt2*(0.01*precip1_prof[k]);
               if ( precip2_prof[k] >= 0 ) precip2_prof_sum[k]+= wt2*(0.01*precip2_prof[k]);
             }

             if ( idb == 1 ) {
               if ( prof_idx[k2].ndpr_MS > 0 ) {
                 for (k= 0; k< 3; k++) {
                   if ( prof_idx[k2].type_precip_MS[k] > 0 ) type_precip[k]+= wt2*prof_idx[k2].type_precip_MS[k]/prof_idx[k2].ndpr_MS;
                 }
                 for (k= 0; k< 5; k++) {
                   if ( prof_idx[k2].shallow_rain_MS[k] > 0 ) shallow_rain[k]+= wt2*prof_idx[k2].shallow_rain_MS[k]/prof_idx[k2].ndpr_MS;
                 }
               }
             } 
             else {
               if ( prof_idx[k2].ndpr_NS > 0 ) {
                 for (k= 0; k< 3; k++) {
                   if ( prof_idx[k2].type_precip_NS[k] > 0 ) type_precip[k]+= wt2*prof_idx[k2].type_precip_NS[k]/prof_idx[k2].ndpr_NS;
                 }
                 for (k= 0; k< 5; k++) {
                   if ( prof_idx[k2].shallow_rain_NS[k] > 0 ) shallow_rain[k]+= wt2*prof_idx[k2].shallow_rain_NS[k]/prof_idx[k2].ndpr_NS;
                 }
               }
             }

             ts= prof_idx[k2].ts;
             ts_sum+= wt2*ts;

             t2m= prof_idx[k2].t2m;
             t2m_sum+= wt2*t2m;

             tqv= prof_idx[k2].tqv;
             tqv_sum+= wt2*tqv;

             if ( h273 > -2000 && h273 < 9000 ) {
               h273_sum+= wt2*h273;
               wt_h273_sum+= wt2;
             }

             iclass= prof_idx[k2].sfc_class;
             if ( iclass >= 0 && iclass < NCLASS ) class[iclass]++;

             /* Store top-ranked Ku/Ka profiles and TB */
             if ( wt2 == 1.0 ) {
               for (k= 0; k< NLEV_DPR; k++) {
                 idx1= nsat*ns_S1*NLEV_DPR + jsat*NLEV_DPR + k;
                 if ( idb == 0 ) {
                   z_ku_1_NS_buf[idx1]= prof_idx[k2].z_ku[k];
                 }
                 else {
                   z_ku_1_MS_buf[idx1]= prof_idx[k2].z_ku[k];
                   z_ka_1_MS_buf[idx1]= prof_idx[k2].z_ka[k];
                 }
               }
               for (k=0; k< nc_out; k++) {
                 idx1= nsat*ns_S1*nc_out + jsat*nc_out + k;
                 if ( idb == 0 ) 
                   tb_1_NS_buf[idx1]= prof_idx[k2].tb[k];
                 else
                   tb_1_MS_buf[idx1]= prof_idx[k2].tb[k];
               }
             }

             /* CCFADS of database Ku and Ka Z profiles */
             /* NOTE below only accumulating CCFADS when wt2 > 0.1 */

             ku_clutter= ka_clutter= -1;
             for (k= 0; k< NLEV_DPR; k++) {
               kc= k - (NLEV_DPR - NLEV_PRECIP);
               ht= 0.25*(NLEV_DPR-k-1);
               iht= NHT*(ht-HTMIN)/(HTMAX-HTMIN);

               zlog= 0.01*prof_idx[k2].z_ku[k];
               if ( kc >= 0 && kc < NLEV_PRECIP ) {
                 if ( ku_clutter < 0 && prof_idx[k2].nclutter_ku[kc] > 0 ) ku_clutter= k;
               }
               if ( zlog > 0 ) z_ku[k]+= pow(10.0, 0.1*zlog);
               nz_ku[k]++;
               if ( ku_clutter < 0 ) {
                 if ( zlog > zmax ) zmax= zlog;
                 iz= NZ*(zlog-ZMIN)/(ZMAX-ZMIN);
                 if ( inside_NS > 0 && wt2 > 0.5 && iht >= 0 && iht < NHT && iz >= 0 && iz < NZ ) ccfads_ku[iht][iz][idb]++;
               }

               zlog= 0.01*prof_idx[k2].z_ka[k];
               if ( kc >= 0 && kc < NLEV_PRECIP ) {
                 if ( ka_clutter < 0 && prof_idx[k2].nclutter_ka[kc] > 0 ) ka_clutter= k;
               }
               if ( zlog > 0 ) z_ka[k]+= pow(10.0, 0.1*zlog);
               nz_ka[k]++;
               if ( ka_clutter < 0 ) {
                 if ( zlog > zmax ) zmax= zlog;
                 iz= NZ*(zlog-ZMIN)/(ZMAX-ZMIN);
                 if ( inside_MS > 0 && wt2 > 0.5 && iht >= 0 && iht < NHT && iz >= 0 && iz < NZ ) ccfads_ka[iht][iz][idb]++;
               }
             }


             /* 20-dB top of Ku and Ka Z profiles */
             if ( ku_clutter > 0 ) {
               ku_ztop= -1.0;
               for (k= 1; k< ku_clutter-3; k++) {
                 ht= 0.25*(NLEV_DPR-k-1);
                 zlog1= 0.01*prof_idx[k2].z_ku[k-1];
                 zlog2= 0.01*prof_idx[k2].z_ku[k];
                 zlog3= 0.01*prof_idx[k2].z_ku[k+1];
                 if ( zlog1 > 20 && zlog2 > 20 && zlog3 > 20 ) {
                   ku_ztop= ht+0.25;
                   break;
                 }
               }
               if ( ku_ztop > 0 ) ku_ztop_sum+= wt2*ku_ztop;
               wt_ku_ztop_sum+= wt2;
             }
             if ( ka_clutter > 0 ) {
               ka_ztop= -1.0;
               for (k= 1; k< ka_clutter-3; k++) {
                 ht= 0.25*(NLEV_DPR-k-1);
                 zlog1= 0.01*prof_idx[k2].z_ka[k-1];
                 zlog2= 0.01*prof_idx[k2].z_ka[k];
                 zlog3= 0.01*prof_idx[k2].z_ka[k+1];
                 if ( zlog1 > 20 && zlog2 > 20 && zlog3 > 20 ) {
                   ka_ztop= ht+0.25;
                   break;
                 }
               }
               if ( ka_ztop > 0 ) ka_ztop_sum+= wt2*ka_ztop;
               wt_ka_ztop_sum+= wt2;
             }

             if ( nfound < 10 ) {
               nku20+= prof_idx[k2].nku20;
               nku25+= prof_idx[k2].nku25;
             }

             /*  PDF precip */
             if ( precip1 < 0.1 ) iprecip= 0;
             else iprecip= NPRECIP*(precip1-PRECIP1)/(PRECIP2-PRECIP1) + 1;
             /*printf("precip=%8.3f iprecip=%d\n", precip1, iprecip);*/
             if ( iprecip >= 0 && iprecip < NPRECIP ) hist_precip[iprecip]++;

             sprintf(pdate,"%04d/%02d/%02d %02d:%02d:%02d", prof_idx[k2].yyyy, prof_idx[k2].mm, prof_idx[k2].dd, prof_idx[k2].hh, prof_idx[k2].mn, prof_idx[k2].ss);

             /* Also write to output file for scatter plots */

             if ( nfound < 10 ) {
               printf("DB  %6d %3d %6d %s %4d %3d %7.2f %7.2f ", idx_db, prof_idx[k2].satid, prof_idx[k2].rev, pdate, prof_idx[k2].j_S1, prof_idx[k2].j_NS, prof_idx[k2].glat, prof_idx[k2].glon);
               for (k= 0; k< 3; k++) printf("%8.3f ", prof_idx[k2].pc_emis[k]);
               for (k= 0; k< 9; k++) printf("%6.2f ", prof_idx[k2].tb[k]);
               printf("%3d ", prof_idx[k2].sfc_class);
               printf("%7.2f %7.2f %7.2f %6d ", prof_idx[k2].ts, prof_idx[k2].t2m, prof_idx[k2].tqv, (int) h273);
               printf("%4d %4d ", prof_idx[k2].nku20, prof_idx[k2].nku25);
               printf("%7.2f ", prof_idx[k2].precip_GPROF);
               printf("%7.2f %7.2f %6d ", precip1, precip2, elev);
               e0= e1= -1.0;
               if ( prof_idx[k2].nku20 == 0 && prof_idx[k2].nku25 == 0 ) {
                 e0= ( prof_idx[k2].emis[0] > 0.0 ) ? prof_idx[k2].emis[0] : -1.0;
                 e1= ( prof_idx[k2].emis[1] > 0.0 ) ? prof_idx[k2].emis[1] : -1.0;
               }
               printf("%6.3f %6.3f ", e0, e1);
               e0= ( prof_idx[k2].emis_NS_cmb[0] > 0.0 && prof_idx[k2].emis_NS_cmb[0] < 2.0 ) ? prof_idx[k2].emis_NS_cmb[0] : -1.0;
               e1= ( prof_idx[k2].emis_NS_cmb[1] > 0.0 && prof_idx[k2].emis_NS_cmb[1] < 2.0 ) ? prof_idx[k2].emis_NS_cmb[1] : -1.0;
               printf("%6.3f %6.3f ", e0, e1);
               printf("%7.4f ", wt2);
               for (k= 0; k< 3; k++) printf("%3d ", prof_idx[k2].type_precip_NS[k]);
               printf("%5.2f ", ku_ztop);
               if ( idb == 1 ) printf("%5.2f ", ka_ztop);
               printf("\n");
             }

             nfound++;
             nfound_all++;

             if ( wt2 < WT2MIN || n_precip_all_sum >= N2MAX ) {
               /*if ( idb == 0 ) epc_qual_NS_buf[idx_out]= nfound;
               else epc_qual_MS_buf[idx_out]= nfound;*/
               break;
             }

           }  /* nrec2 loop */

           /*--- All done with lookup ---*/

           if ( n_precip_all_sum == 0 ) {
             printf("EXPAND No database entries located  idx=%d idb=%d idx_db=%d nrec2=%d maxrec=%d wt2=%f nfound=%d nmodel=%d  TB=", idx, idb, idx_db, nrec2, maxrec, wt2, nfound, nmodel);
             for (k= 0; k< NCHAN; k++) printf("%7.2f ", tb_atms[k]);
             printf("%3d ", sfc_class);
             printf("%7.2f %7.2f %7.2f %6d ", ts_MERRA2, t2m_MERRA2, tqv_MERRA2, h273_MERRA2);
             printf("\n");

             /*exit(1);*/
   
             if ( idb == 0 ) epc_qual_NS_buf[idx_out]= -2;
             else epc_qual_MS_buf[idx_out]= -2;
             continue;
           }

           if ( idb == 0 ) epc_qual_NS_buf[idx_out]= nfound;
           else epc_qual_MS_buf[idx_out]= nfound;

           precip1= precip1_sum/wt_precip_sum;
           if ( idb == 0 ) {
             precip1_NS_buf[idx_out]= precip1;
             for (k= 0; k< NLEV_PRECIP; k++) {
               idx3_out= nsat*ns_S1*NLEV_PRECIP + jsat*NLEV_PRECIP + k;
               precip1_prof_NS_buf[idx3_out]= 100*(precip1_prof_sum[k]/wt_precip_sum);
             }
           }
           else {
             precip1_MS_buf[idx_out]= precip1;
             for (k= 0; k< NLEV_PRECIP; k++) {
               idx3_out= nsat*ns_S1*NLEV_PRECIP + jsat*NLEV_PRECIP + k;
               precip1_prof_MS_buf[idx3_out]= 100*(precip1_prof_sum[k]/wt_precip_sum);
             }
           }

           precip2= precip2_sum/wt_precip_sum;
           if ( idb == 0 ) {
             precip2_NS_buf[idx_out]= precip2;
             for (k= 0; k< NLEV_PRECIP; k++) {
               idx3_out= nsat*ns_S1*NLEV_PRECIP + jsat*NLEV_PRECIP + k;
               precip2_prof_NS_buf[idx3_out]= 100*(precip2_prof_sum[k]/wt_precip_sum);
             }
           }
           else {
             precip2_MS_buf[idx_out]= precip2;
             for (k= 0; k< NLEV_PRECIP; k++) {
               idx3_out= nsat*ns_S1*NLEV_PRECIP + jsat*NLEV_PRECIP + k;
               precip2_prof_MS_buf[idx3_out]= 100*(precip2_prof_sum[k]/wt_precip_sum);
             }
           }

           /* determine type precip by majority vote */
           i1max= -1; fmax= 0;
           for (k= 0; k< 3; k++) {
             if ( type_precip[k] > fmax ) { i1max= k; fmax= type_precip[k]; }
           }
           if ( i1max >= 0 ) {
             if ( idb == 0 ) type_precip_NS_buf[idx_out]= i1max;
             else type_precip_MS_buf[idx_out]= i1max;
           }

           /* determine shallow rain by majority vote, clump various cases together */
           i2max= -1; fmax= 0;
           if ( shallow_rain[0] > fmax ) { i2max= 0; fmax= shallow_rain[0]; }
           if ( shallow_rain[1] + shallow_rain[2] + shallow_rain[3] + shallow_rain[4] > fmax ) i2max= 1;
           if ( i2max >= 0 ) {
             if ( idb == 0 ) shallow_rain_NS_buf[idx_out]= i2max;
             else shallow_rain_MS_buf[idx_out]= i2max;
           }

           printf("EST N_type= ");
           for (k= 0; k< 3; k++) printf("%8.4f ", type_precip[k]/wt_precip_sum);
           printf("    N_shallow_rain= ");
           printf("%8.4f ", shallow_rain[0]/wt_precip_sum);
           printf("%8.4f ", (shallow_rain[1] + shallow_rain[2] + shallow_rain[3] + shallow_rain[4])/wt_precip_sum);
           printf("    R= %7.3f %6d %6d", precip2, i1max, i2max);
           printf("\n");



           prob_precip1= 100.0*n_precip1_sum/n_precip_all_sum;
           if ( idb == 0 ) precip1_prob_NS_buf[idx_out]= prob_precip1;
           else precip1_prob_MS_buf[idx_out]= prob_precip1;

           prob_precip2= 100.0*n_precip2_sum/n_precip_all_sum;
           if ( idb == 0 ) precip2_prob_NS_buf[idx_out]= prob_precip2;
           else precip2_prob_MS_buf[idx_out]= prob_precip2;
 
           ts= ts_sum/wt_precip_sum;
           if ( idb == 0 ) ts_epc_NS_buf[idx_out]= ts;
           else ts_epc_MS_buf[idx_out]= ts;

           t2m= t2m_sum/wt_precip_sum;
           if ( idb == 0 ) t2m_epc_NS_buf[idx_out]= t2m;
           else t2m_epc_MS_buf[idx_out]= t2m;

           tqv= tqv_sum/wt_precip_sum;
           if ( idb == 0 ) tqv_epc_NS_buf[idx_out]= tqv;
           else tqv_epc_MS_buf[idx_out]= tqv;

           if ( wt_h273_sum > 0.0 ) {
             h273= h273_sum/wt_h273_sum;
             if ( idb == 0 ) h273_epc_NS_buf[idx_out]= h273;
             else h273_epc_MS_buf[idx_out]= h273;
           }
           else
             h273= -1;

           if ( wt_ku_ztop_sum > 0.0 ) {
             ku_ztop= ku_ztop_sum/wt_ku_ztop_sum;
             if ( idb == 0 ) ku_ztop_epc_NS_buf[idx_out]= 1000*ku_ztop;
             else ku_ztop_epc_MS_buf[idx_out]= 1000*ku_ztop;
           }
           else
             ku_ztop= -1.0;

           if ( idb == 1 ) {
             if ( wt_ka_ztop_sum > 0.0 ) {
               ka_ztop= ka_ztop_sum/wt_ka_ztop_sum;
               ka_ztop_epc_MS_buf[idx_out]= 1000*ka_ztop;
             }
           }
           else
             ka_ztop= -1.0;


           printf("EST %6d %3d %6d %s %4d %3d %7.2f %7.2f ", idx_db, isatname, irev, tdate, isat, jsat, glat1, glon1);
           for (k= 0; k< 3; k++) printf("%8.3f ", pc_emis[k]);
           for (k= 0; k< 9; k++) printf("%6.2f ", tb_atms[k]);
           printf("%3d ", sfc_class);
           printf("%7.2f %7.2f %7.2f ", ts, t2m, tqv);
           printf("%6d ", (int) h273);
           printf("%4d %4d ", nku20, nku25);
           printf("%7.2f ", precip_GPROF);
           printf("%7.2f %7.2f ", precip1, precip2);
           printf("%7.2f %7.2f %7.2f ", prob_precip_GPROF, prob_precip1, prob_precip2);
           printf("%2d %2d ", i1max, i2max);
           printf("%5.2f ", ku_ztop);
           if ( idb == 1 ) printf("%5.2f ", ka_ztop);
           printf("\n");

           if ( zmax > 40 && precip1 > 25 ) {
             for (k= 0; k< NLEV_DPR; k++) {
               zlog1= ( z_ku[k] > 0 && nz_ku[k] > 0 ) ? 10.0*log10(z_ku[k]/nz_ku[k]) : -99.99;
               zlog2= ( z_ka[k] > 0 && nz_ka[k] > 0 ) ? 10.0*log10(z_ka[k]/nz_ka[k]) : -99.99;
               printf("%3d %6d %6.2f %6d %6.2f\n", k, nz_ku[k], zlog1, nz_ka[k], zlog2);
             }
             for (iprecip= 0; iprecip< NPRECIP; iprecip++) {
               x= ( iprecip == 0 ) ? 0.0 : PRECIP1 + iprecip*DPRECIP;
               if ( hist_precip[iprecip] > 0 ) printf("%3d %6.2f N=%d\n", iprecip, x, hist_precip[iprecip]);
             }
           }

         }   /* idb loop for NS or MS retrievals */

         npix++;

       }  /* next radiometer 1C sample */

       nsat++;  /* increment output buffer line index */

     }  /* next radiomter 1C line */

     printf("\n%s KuKa=%2d DBidx=%5d Frac=%6.4f %6.4f Nread=%7d Nrain=%7d Npix=%6d  Percent complete=%7.3f\n", tbfile, idb, idx_db, frac0, frac1, nrec, nrain, ndb, 100.0*npix/((isat_end-isat_start+1)*ns_S1) );

   }  /* next DB file */

   /*-------------------------------------------------------------------------*/
   /* Save CCFADS files */

   if ( have_dpr > 0 ) {

   for (idb= 0; idb< 2; idb++) {
     if ( idb == 0 ) {
       sprintf(outfile, "%s.ccfads_NS.txt", outfilename);
       sprintf(outfile2, "%s.ccfads_contours_NS.txt", outfilename);
     }
     if ( idb == 1 ) {
       sprintf(outfile, "%s.ccfads_MS.txt", outfilename);
       sprintf(outfile2, "%s.ccfads_contours_MS.txt", outfilename);
     }
     if ((fccfads= fopen(outfile,"w")) == 0) {
       printf("Unable to open %s\n", outfile);
       exit(-1);
     }
     printf("Opened %s\n", outfile);
     if ((fccfads2= fopen(outfile2,"w")) == 0) {
       printf("Unable to open %s\n", outfile2);
       exit(-1);
     }
     printf("Opened %s\n", outfile2);

     /* print from top-down */
     for (k= NHT-1; k>= 0; k--) {
       ht= DHT*(k+1);

       sum1= sum2= 0;
       for (k1= 0; k1< NZ; k1++) {
         sum1+= ccfads_ku[k][k1][idb];
         sum2+= ccfads_ka[k][k1][idb];
       }

       for (k1= 0; k1< NZ; k1++) ccfads1[k1]= -99.99;
       if ( sum1 > 2 ) {
         for (k1= 0; k1< NZ; k1++) ccfads1[k1]= 100.0*ccfads_ku[k][k1][idb]/sum1;
       }

       for (k1= 0; k1< NZ; k1++) ccfads2[k1]= -99.99;
       if ( sum2 > 2 ) {
         for (k1= 0; k1< NZ; k1++) ccfads2[k1]= 100.0*ccfads_ka[k][k1][idb]/sum2;
       }

       for (k1= 0; k1< NZ; k1++) {
         zlog= ZMIN + k1*DZ;
         fprintf(fccfads, "%6.2f %6.2f %8d %9.5f ", ht, zlog, ccfads_ku[k][k1][idb], ccfads1[k1]);
         if ( idb == 1 ) fprintf(fccfads, "%8d %9.5f", ccfads_ka[k][k1][idb], ccfads2[k1]);
         fprintf(fccfads,"\n");
       }

       /*-----------------------------------------------*/
       /*  CDF contour lines every 10 percent */

       for (j1= 0; j1< 10; j1++) cdf_ccfads[j1]= -1.0;
       if ( sum1 > 100 ) {
         cdf_sum= 0;
         for (k1= 0; k1< NZ; k1++) {
           zlog= ZMIN + k1*DZ;
           cdf_sum+= ccfads_ku[k][k1][idb];
           for (j1= 0; j1< 10; j1++) {
             pcent= 10.0 + 10.0*j1;
             if ( j1 == 9 ) pcent= 95.0;
             pcent1= 100.0*cdf_sum/sum1;
             if ( cdf_ccfads[j1] < 0.0 && pcent1 >= pcent ) {
               cdf_ccfads[j1]= zlog;
             }
           }
         }
       }
       fprintf(fccfads2,"%3d %7.2f %10d ", k, ht, sum1);
       for (j1= 0; j1< 10; j1++) fprintf(fccfads2, "%6.2f ", cdf_ccfads[j1]);

       if ( idb == 1 ) {
         for (j1= 0; j1< 10; j1++) cdf_ccfads[j1]= -1.0;
         if ( sum2 > 100 ) {
           cdf_sum= 0;
           for (k1= 0; k1< NZ; k1++) {
             zlog= ZMIN + k1*DZ;
             cdf_sum+= ccfads_ka[k][k1][idb];
             for (j1= 0; j1< 10; j1++) {
               pcent= 10.0 + 10.0*j1;
               if ( j1 == 9 ) pcent= 95.0;
               pcent1= 100.0*cdf_sum/sum2;
               if ( cdf_ccfads[j1] < 0.0 && pcent1 >= pcent ) {
                 cdf_ccfads[j1]= zlog;
               }
             }
           }
         }
         fprintf(fccfads2,"%10d ", sum2);
         for (j1= 0; j1< 10; j1++) fprintf(fccfads2, "%6.2f ", cdf_ccfads[j1]);
       }
       fprintf(fccfads2, "\n");
       /*-----------------------------------------------*/

     }
     fclose(fccfads);
     fclose(fccfads2);
     printf("Closed %s\n", outfile);
     printf("Closed %s\n", outfile2);
   }

   }

   /*-------------------------------------------------------------------------*/
   /*--- Add in actual DPR profiles if provided ---*/

   if ( have_dpr > 0 ) {
     printf("Adding DPR measured Ku and Ka profiles\n");

     z_ku_DPR_buf= (short *) malloc(nl_out*ns_S1*NLEV_DPR*sizeof(short));
     z_ka_DPR_buf= (short *) malloc(nl_out*ns_S1*NLEV_DPR*sizeof(short));
     elev_DPR_buf= (short *) malloc(nl_out*ns_S1*sizeof(short));
     precip_NS_DPR_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
     precip_MS_DPR_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
     lat_NS_DPR_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
     lon_NS_DPR_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
     ztop_ku_DPR_buf= (short *) malloc(nl_out*ns_S1*sizeof(short));
     ztop_ka_DPR_buf= (short *) malloc(nl_out*ns_S1*sizeof(short));
     if ( have_cmb > 0 ) {
       precip_NS_CMB_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
       precip_MS_CMB_buf= (float *) malloc(nl_out*ns_S1*sizeof(float));
       precip_prof_NS_CMB_buf= (short *) malloc(nl_out*ns_S1*NLEV_PRECIP*sizeof(short));
       precip_prof_MS_CMB_buf= (short *) malloc(nl_out*ns_S1*NLEV_PRECIP*sizeof(short));
     }
     type_precip_DPR_buf= (short *) malloc(nl_out*ns_S1*sizeof(short)); /* 0=convective 1=stratiform 2=other */
     shallow_rain_DPR_buf= (short *) malloc(nl_out*ns_S1*sizeof(short)); /* 0=no shallow 1=shallow */

     for (i=0; i < nl_out*ns_S1*NLEV_DPR; i++) {
       z_ku_DPR_buf[i]= z_ka_DPR_buf[i]= -9999;
     }
     for (i=0; i < nl_out*ns_S1; i++) {
       precip_NS_DPR_buf[i]= -9999.0;
       precip_MS_DPR_buf[i]= -9999.0;
       elev_DPR_buf[i]= -9999;
       if ( have_cmb > 0 ) {
         precip_NS_CMB_buf[i]= -9999.0;
         precip_MS_CMB_buf[i]= -9999.0;
       }
       lat_NS_DPR_buf[i]= -9999.0;
       lon_NS_DPR_buf[i]= -9999.0;
       ztop_ku_DPR_buf[i]= -9999;
       ztop_ka_DPR_buf[i]= -9999;
       type_precip_DPR_buf[i]= -9999;
       shallow_rain_DPR_buf[i]= -9999;
     }

     for (k= 0; k< NHT; k++) {
       for (k1= 0; k1< NZ; k1++) {
         ccfads_ku[k][k1][0]= 0;
         ccfads_ka[k][k1][0]= 0;
         ccfads_ku[k][k1][1]= 0;
         ccfads_ka[k][k1][1]= 0;
       }
     }

     nsat= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         glat1= lat_S1_buf[idx];
         glon1= lon_S1_buf[idx];

         incidenceAngle= fabs(incidenceAngle_S1_buf[idx]);
         if ( pixres_nadir >= 0.0 ) pixres= pixres_nadir/cos(DTR*incidenceAngle);
         npixres= ((pixres*pixres)/25) + 0.5;
         ndpr= (pixres/10.0) + 3;
         /*printf("%4d %4d pixres=%8.3f Npixres=%4d Ndpr=%4d N=%d\n", isat, jsat, pixres, npixres, ndpr, nrec);*/

         /* Use nearest NS radar pixel found earlier */
         imin_NS= dpr_NS_i[idx];
         jmin_NS= dpr_NS_j[idx];
         if ( imin_NS < 0 || jmin_NS < 0 ) continue;

         idx1= imin_NS*ns_NS + jmin_NS;
         nslat1= lat_NS_buf[idx1];
         nslon1= lon_NS_buf[idx1];
         km= RADEARTH*acos(cos(DTR*glon1-DTR*nslon1)*cos(DTR*glat1)*cos(DTR*nslat1) + sin(DTR*glat1)*sin(DTR*nslat1));
         /*  IS IT NEEDED   precip2= precipRateESurface_NS_buf[idx1];*/
         elev= elev_NS_buf[idx1];

         idx2= nsat*ns_S1 + jsat;
         elev_DPR_buf[idx2]= elev;
         lat_NS_DPR_buf[idx2]= nslat1;
         lon_NS_DPR_buf[idx2]= nslon1;

         if ( jmin_NS == ns_NS/2) printf("S1=%4d %4d %8.3f %8.3f  NS=%4d %4d %8.3f %8.3f   %8.2f %6d\n",
           isat,jsat,glat1,glon1, imin_NS,jmin_NS,nslat1,nslon1, km,elev);

         /* Initialize */

         for (k=0; k < nlev_NS; k++) {
           z_ku[k]= z_ka[k]= 0.0;
           nz_ku[k]= nz_ka[k]= 0;
           nclutter_ku[k]= nclutter_ka[k]= 0;
         }
         for (k=0; k < NLEV_PRECIP; k++) {
           n_precip2_prof_sum[k]= n_precip1_prof_sum[k]= 0;
           precip2_prof_sum[k]= precip1_prof_sum[k]= 0.0;
         }
         n_precip2_sum= n_precip1_sum= 0;
         precip2_sum= precip1_sum= 0.0;
         ku_ztop_sum= wt_ku_ztop_sum= 0.0;
         for (k=0; k< 3; k++) typePrecip[k]= 0;
         for (k=0; k< 5; k++) flagShallowRain[k]= 0;

         /* Average around nearest NS */

         ndpr_total= 0;
         for (i1=imin_NS-ndpr; i1 <= imin_NS+ndpr; i1++) {
           if ( i1 < 0 || i1 >= nl_NS ) continue;
           qual= qual_NS_buf[i1];
           if ( qual != 0 ) continue;
           for (j1=jmin_NS-ndpr; j1 <= jmin_NS+ndpr; j1++) {
             if ( j1 < 0 || j1 >= ns_NS ) continue;
             idx1= i1*ns_NS + j1;
             ibin_sfc= binRealSurface_NS_buf[idx1];
             ibin_clutter_free= binClutterFreeBottom_NS_buf[idx1];

             if ( ibin_sfc < 0 || ibin_sfc > nlev_NS ) continue;
             if ( ibin_clutter_free < 0 || ibin_clutter_free > nlev_NS ) continue;
 
             nslat2= lat_NS_buf[idx1];
             nslon2= lon_NS_buf[idx1];
             km= RADEARTH*acos(cos(DTR*nslon1-DTR*nslon2)*cos(DTR*nslat1)*cos(DTR*nslat2) + sin(DTR*nslat1)*sin(DTR*nslat2));
             /*printf(" NS inc=%6.2f sub=%5d %5d  %8.3f %8.3f %8.3f %8.3f  km=%8.3f N=%4d\n", incidenceAngle, imin_NS-i1,jmin_NS-j1,nslat1,nslon1,nslat2,nslon2, km, ndpr_total);*/
             if ( km > pixres/2.0 ) continue;
             ndpr_total++;

             precip2= precipRateESurface_NS_buf[idx1];
             if ( precip2 >= 0.0 ) {
               n_precip2_sum++;
               precip2_sum+= precip2;
             }
             if ( have_cmb > 0 ) {
               precip1= surfPrecipTotRate_NS_cmb_buf[idx1];
               if ( dataQuality_NS_cmb_buf[i1] == 0 && precip1 >= 0.0 ) {
                 n_precip1_sum++;
                 precip1_sum+= precip1;
               }
               /*  lowest NLEV_PRECIP bins of the nlev_NS_cmb total */
               for (k=0; k < NLEV_PRECIP; k++) {
                 idx3= i1*ns_NS*nlev_NS_cmb + j1*nlev_NS_cmb + (nlev_NS_cmb-NLEV_PRECIP + k);
                 precip1= PrecipTotRate_NS_cmb_buf[idx3];
                 if ( dataQuality_NS_cmb_buf[i1] == 0 && precip1 >= 0.0 ) {
                   precip1_prof_sum[k]+= precip1;
                   n_precip1_prof_sum[k]++;
                 }
               }
             }

             if ( typePrecip_NS_buf[idx1]/10000000 == 1 ) typePrecip[0]++;
             if ( typePrecip_NS_buf[idx1]/10000000 == 2 ) typePrecip[1]++;
             if ( typePrecip_NS_buf[idx1]/10000000 == 3 ) typePrecip[2]++;
             if ( flagShallowRain_NS_buf[idx1] == 0 ) flagShallowRain[0]++;
             if ( flagShallowRain_NS_buf[idx1] == 10 ) flagShallowRain[1]++;
             if ( flagShallowRain_NS_buf[idx1] == 11 ) flagShallowRain[2]++;
             if ( flagShallowRain_NS_buf[idx1] == 20 ) flagShallowRain[3]++;
             if ( flagShallowRain_NS_buf[idx1] == 21 ) flagShallowRain[4]++;

             ku_ztop= -1.0;
             for (k=0; k < ibin_sfc; k++) {
               if ( k >= nlev_NS ) continue;
               ht= 0.125*(nlev_NS-k-1);
               k2= k/2;  /* save z_ku at coarser resolution */
               idx1= i1*ns_NS*nlev_NS + j1*nlev_NS + k;
               z= zFactorMeasured_NS_buf[idx1];
               /*printf("NS bin=%5d sfc=%5d noClutter=%5d ZKu=%8.3f\n",k,ibin_sfc, ibin_clutter_free, z);*/
               if ( z < 10.0 ) continue;
               z_ku[k2]+= pow(10.0, 0.1*z);
               nz_ku[k2]++;
               if ( k >= ibin_clutter_free ) nclutter_ku[k2]++;
               if ( k > 0 && k < ibin_clutter_free-3 ) {
                 zlog1= zFactorMeasured_NS_buf[idx1-1];
                 zlog2= zFactorMeasured_NS_buf[idx1];
                 zlog3= zFactorMeasured_NS_buf[idx1+1];
                 if ( ku_ztop < 0 && zlog1 > 20 && zlog2 > 20 && zlog3 > 20 ) ku_ztop= ht+0.125;
               }
             }
             if ( ku_ztop > 0 ) ku_ztop_sum+= ku_ztop;
             wt_ku_ztop_sum+= 1.0;
           }
         }

         pcent= 100.0*ndpr_total/npixres;
         /*printf("NS PIXFILL %4d %4d  %4d %4d pixres=%8.3f Npixres=%4d Ndpr=%4d Ndpr_total=%d  percent=%8.3f\n", isat, jsat, imin_NS, jmin_NS, pixres, npixres, ndpr, ndpr_total, pcent);*/
         if ( pcent < 20 ) {
           printf("No good Ku-band  Percent= %8.3f\n", pcent);
           continue;
         }

         for (k=0; k < NLEV_DPR; k++) {
           ht= 0.25*(NLEV_DPR-k-1);
           iht= NHT*(ht-HTMIN)/(HTMAX-HTMIN);

           if ( nz_ku[k] == 0 ) continue;
           idx1= nsat*ns_S1*NLEV_DPR + jsat*NLEV_DPR + k;
           z= 10.0*log10(z_ku[k]/nz_ku[k]);
           if ( nclutter_ku[k] > 0 ) z*= -1.0;
           z_ku_DPR_buf[idx1]= 100*z;

           iz= NZ*(z-ZMIN)/(ZMAX-ZMIN);
           if ( iht >= 0 && iht < NHT && iz >= 0 && iz < NZ ) ccfads_ku[iht][iz][0]++;
         }

         if ( n_precip2_sum > 0 ) {
           precip2= precip2_sum/n_precip2_sum;
           precip_NS_DPR_buf[idx2]= precip2;
         }
         if ( n_precip1_sum > 0 ) {
           precip1= precip1_sum/n_precip1_sum;
           precip_NS_CMB_buf[idx2]= precip1;
         }
         for (k=0; k < NLEV_PRECIP; k++) {
           if ( n_precip1_prof_sum[k] > 0 ) {
             idx3= nsat*ns_S1*NLEV_PRECIP + jsat*NLEV_PRECIP + k;
             precip1= precip1_prof_sum[k]/n_precip1_prof_sum[k];
             precip_prof_NS_CMB_buf[idx3]= 100*precip1;
           }
         }

         if ( wt_ku_ztop_sum > 0.0 ) {
           ku_ztop= ku_ztop_sum/wt_ku_ztop_sum;
           ztop_ku_DPR_buf[idx2]= 1000*ku_ztop;
         }

         /* determine type precip by majority vote */ 
         imax= -1; nmax= 0; 
         for (k= 0; k< 3; k++) {
           if ( typePrecip[k] > nmax ) { imax= k; nmax= typePrecip[k]; }
         }
         if ( imax >= 0 ) type_precip_DPR_buf[idx2]= imax;

         /* determine shallow rain by majority vote, clump various cases together */ 
         imax= -1; nmax= 0; 
         if ( flagShallowRain[0] > nmax ) { imax= 0; nmax= flagShallowRain[0]; }
         if ( flagShallowRain[1] + flagShallowRain[2] + flagShallowRain[3] + flagShallowRain[4] > nmax ) imax= 1; 
         if ( imax >= 0 ) shallow_rain_DPR_buf[idx2]= imax;

         /*printf("OBS N_type= ");
         for (k= 0; k< 3; k++) printf("%8.4f ", 1.0*typePrecip[k]/ndpr_total);
         printf("    N_shallow_rain= ");
         printf("%8.4f ", 1.0*flagShallowRain[0]/ndpr_total);
         printf("%8.4f ", 1.0*(flagShallowRain[1] + flagShallowRain[2] + flagShallowRain[3] + flagShallowRain[4])/ndpr_total);
         printf("    R= %7.3f ", precip2);
         printf("\n");*/

         /* Next find nearest MS */
         /* NS positions 12-36 correspond to MS positions 0-24 */

         imin_MS= imin_NS;
         jmin_MS= jmin_NS - 12;
         if ( jmin_MS < 0 || jmin_MS >= ns_MS ) continue;

         n_precip2_sum= n_precip1_sum= 0;
         precip2_sum= precip1_sum= 0.0;
         for (k=0; k < NLEV_PRECIP; k++) {
           n_precip2_prof_sum[k]= n_precip1_prof_sum[k]= 0;
           precip2_prof_sum[k]= precip1_prof_sum[k]= 0.0;
         }
         ka_ztop_sum= wt_ka_ztop_sum= 0.0;

         for (i1=imin_MS-ndpr; i1 <= imin_MS+ndpr; i1++) {
           if ( i1 < 0 || i1 >= nl_MS ) continue;
           qual= qual_MS_buf[i1];
           if ( qual != 0 ) continue;
           for (j1=jmin_MS-ndpr; j1 <= jmin_MS+ndpr; j1++) {
             if ( j1 < 0 || j1 >= ns_MS ) continue;
             idx1= i1*ns_MS + j1;
             ibin_sfc= binRealSurface_MS_buf[idx1];
             ibin_clutter_free= binClutterFreeBottom_MS_buf[idx1];

             if ( ibin_sfc < 0 || ibin_sfc > nlev_MS ) continue;
             if ( ibin_clutter_free < 0 || ibin_clutter_free > nlev_MS ) continue;

             precip2= precipRateESurface_MS_buf[idx1];
             if ( precip2 >= 0.0 ) {
               n_precip2_sum++;
               precip2_sum+= precip2;
             }
             if ( have_cmb > 0 ) {
               precip1= surfPrecipTotRate_MS_cmb_buf[idx1];
               if ( dataQuality_MS_cmb_buf[i1] == 0 && precip1 >= 0.0 ) {
                 n_precip1_sum++;
                 precip1_sum+= precip1;
               }
               /*  lowest NLEV_PRECIP bins of the nlev_MS_cmb total */
               for (k=0; k < NLEV_PRECIP; k++) {
                 idx3= i1*ns_MS*nlev_MS_cmb + j1*nlev_MS_cmb + (nlev_MS_cmb-NLEV_PRECIP + k);
                 precip1= PrecipTotRate_MS_cmb_buf[idx3];
                 if ( dataQuality_MS_cmb_buf[i1] == 0 && precip1 >= 0.0 ) {
                   precip1_prof_sum[k]+= precip1;
                   n_precip1_prof_sum[k]++;
                 }
               }
             }

             ka_ztop= -1.0;
             for (k=0; k < ibin_sfc; k++) {
               if ( k >= nlev_NS ) continue;
               ht= 0.125*(nlev_MS-k-1);
               k2= k/2;  /* save z_ka at coarser resolution */
               idx1= i1*ns_MS*nlev_MS + j1*nlev_MS + k;
               z= zFactorMeasured_MS_buf[idx1];
               if ( z < 10.0 ) continue;
               z_ka[k2]+= pow(10.0, 0.1*z);
               nz_ka[k2]++;
               if ( k >= ibin_clutter_free ) nclutter_ka[k2]++;
               /*printf("MS bin=%5d sfc=%5d noClutter=%5d ZKu=%8.3f\n",k,ibin_sfc, ibin_clutter_free, z_ka[k]);*/
               if ( k > 0 && k < ibin_clutter_free-3 ) {
                 zlog1= zFactorMeasured_MS_buf[idx1-1];
                 zlog2= zFactorMeasured_MS_buf[idx1];
                 zlog3= zFactorMeasured_MS_buf[idx1+1];
                 if ( ka_ztop < 0 && zlog1 > 20 && zlog2 > 20 && zlog3 > 20 ) ka_ztop= ht+0.125;
               }
             }
             if ( ka_ztop > 0 ) ka_ztop_sum+= ka_ztop;
             wt_ka_ztop_sum+= 1.0;
           }
         }

         for (k=0; k < NLEV_DPR; k++) {
           ht= 0.25*(NLEV_DPR-k-1);
           iht= NHT*(ht-HTMIN)/(HTMAX-HTMIN);

           if ( nz_ka[k] == 0 ) continue;
           idx1= nsat*ns_S1*NLEV_DPR + jsat*NLEV_DPR + k;
           z= 10.0*log10(z_ka[k]/nz_ka[k]);
           if ( nclutter_ka[k] > 0 ) z*= -1.0;
           z_ka_DPR_buf[idx1]= 100*z;

           iz= NZ*(z-ZMIN)/(ZMAX-ZMIN);
           if ( iht >= 0 && iht < NHT && iz >= 0 && iz < NZ ) ccfads_ka[iht][iz][0]++;
         }

         if ( n_precip2_sum > 0 ) {
           precip2= precip2_sum/n_precip2_sum;
           precip_MS_DPR_buf[idx2]= precip2;
         }
         if ( n_precip1_sum > 0 ) {
           precip1= precip1_sum/n_precip1_sum;
           precip_MS_CMB_buf[idx2]= precip1;
         }
         for (k=0; k < NLEV_PRECIP; k++) {
           if ( n_precip1_prof_sum[k] > 0 ) {
             idx3= nsat*ns_S1*NLEV_PRECIP + jsat*NLEV_PRECIP + k;
             precip1= precip1_prof_sum[k]/n_precip1_prof_sum[k];
             precip_prof_MS_CMB_buf[idx3]= 100*precip1;
           }
         }
         if ( wt_ka_ztop_sum > 0.0 ) {
           ka_ztop= ka_ztop_sum/wt_ka_ztop_sum;
           ztop_ka_DPR_buf[idx2]= 1000*ka_ztop;
         }

       } /* next S1 location */
       nsat++;
     }

     /* Save DPR CCFADS */

     sprintf(outfile, "%s.ccfads_DPR.txt", outfilename);
     sprintf(outfile2, "%s.ccfads_contours_DPR.txt", outfilename);
     if ((fccfads= fopen(outfile,"w")) == 0) {
       printf("Unable to open %s\n", outfile);
       exit(-1);
     }
     printf("Opened %s\n", outfile);
     if ((fccfads2= fopen(outfile2,"w")) == 0) {
       printf("Unable to open %s\n", outfile2);
       exit(-1);
     }
     printf("Opened %s\n", outfile2);

     /* print from top-down */
     for (k= NHT-1; k>= 0; k--) {
       ht= DHT*(k+1);

       sum1= sum2= 0;
       for (k1= 0; k1< NZ; k1++) {
         sum1+= ccfads_ku[k][k1][0];
         sum2+= ccfads_ka[k][k1][0];
       }

       for (k1= 0; k1< NZ; k1++) ccfads1[k1]= -99.99;
       if ( sum1 > 2 ) {
         for (k1= 0; k1< NZ; k1++) ccfads1[k1]= 100.0*ccfads_ku[k][k1][idb]/sum1;
       }

       for (k1= 0; k1< NZ; k1++) ccfads2[k1]= -99.99;
       if ( sum2 > 2 ) {
         for (k1= 0; k1< NZ; k1++) ccfads2[k1]= 100.0*ccfads_ka[k][k1][idb]/sum2;
       }

       for (k1= 0; k1< NZ; k1++) {
         zlog= ZMIN + k1*DZ;
         fprintf(fccfads, "%6.2f %6.2f %8d %9.5f ", ht, zlog, ccfads_ku[k][k1][0], ccfads1[k1]);
         fprintf(fccfads, "%8d %9.5f", ccfads_ka[k][k1][0], ccfads2[k1]);
         fprintf(fccfads,"\n");
       }

       /*-----------------------------------------------*/
       /*  CDF contour lines every 10 percent */

       for (j1= 0; j1< 10; j1++) cdf_ccfads[j1]= -1.0;
       if ( sum1 > 100 ) {
         cdf_sum= 0;
         for (k1= 0; k1< NZ; k1++) {
           zlog= ZMIN + k1*DZ;
           cdf_sum+= ccfads_ku[k][k1][0];
           for (j1= 0; j1< 10; j1++) {
             pcent= 10.0 + 10.0*j1;
             if ( j1 == 9 ) pcent= 95.0;
             pcent1= 100.0*cdf_sum/sum1;
             if ( cdf_ccfads[j1] < 0.0 && pcent1 >= pcent ) {
               cdf_ccfads[j1]= zlog;
             }
           }
         }
       }
       fprintf(fccfads2,"%3d %7.2f %10d ", k, ht, sum1);
       for (j1= 0; j1< 10; j1++) fprintf(fccfads2, "%6.2f ", cdf_ccfads[j1]);

       for (j1= 0; j1< 10; j1++) cdf_ccfads[j1]= -1.0;
       if ( sum2 > 100 ) {
         cdf_sum= 0;
         for (k1= 0; k1< NZ; k1++) {
           zlog= ZMIN + k1*DZ;
           cdf_sum+= ccfads_ka[k][k1][0];
           for (j1= 0; j1< 10; j1++) {
             pcent= 10.0 + 10.0*j1;
             if ( j1 == 9 ) pcent= 95.0;
             pcent1= 100.0*cdf_sum/sum2;
             if ( cdf_ccfads[j1] < 0.0 && pcent1 >= pcent ) {
               cdf_ccfads[j1]= zlog;
             }
           }
         }
       }
       fprintf(fccfads2,"%10d ", sum2);
       for (j1= 0; j1< 10; j1++) fprintf(fccfads2, "%6.2f ", cdf_ccfads[j1]);
       fprintf(fccfads2, "\n");
       /*-----------------------------------------------*/

     }
     fclose(fccfads);
     fclose(fccfads2);
     printf("Closed %s\n", outfile);
     printf("Closed %s\n", outfile2);


   }  /* have_dpr */

   /*-------------------------------------------------------------------------*/


   /*-------------------------------------------------------------------------*/
   /*  output in netCDF */

   if ( have_model > 0 && save_model < 0 ) have_model= -1;

   now= time(NULL);
   t0= now;
   tm= *gmtime(&t0);
   sprintf(processing_date_end,"%4d/%02d/%02d %02d:%02d:%02d", 1900+tm.tm_year, 1+tm.tm_mon, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
   printf("End date %s\n", processing_date_end);

   fbuf= (float *) malloc(nl_out*ns_S1*sizeof(float));
   fbuf3= (float *) malloc(nl_out*ns_S1*nc_out*sizeof(float));
   fbuf4= (float *) malloc(nl_out*ns_S1*NEM*sizeof(float));
   ibuf= (int *) malloc(nl_out*ns_S1*sizeof(int));
   sbuf= (short *) malloc(nl_out*ns_S1*sizeof(short));
   bbuf= (char *) malloc(nl_out*ns_S1*sizeof(char));
   dbuf= (double *) malloc(nl_out*sizeof(double));

   sprintf(outfile, "%s.nc", outfilename);
   if ((retval = nc_create(outfile, NC_NETCDF4, &ncid))) ERR(retval);
   printf("Opened %s\n", outfile);

   sprintf(ctmp,"%s", argv[1]);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "XCAL_1C", strlen(ctmp), ctmp))) ERR(retval);
   if ( have_gprof > 0 ) {
     sprintf(ctmp,"%s", gprof_file);
     if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "GPROF", strlen(ctmp), ctmp))) ERR(retval);
   }
   if ( have_model > 0 ) {
     sprintf(ctmp,"%s", model_file);
     if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "MERRA2", strlen(ctmp), ctmp))) ERR(retval);
   }
   if ( have_dpr > 0 ) {
     sprintf(ctmp,"%s", dpr_file);
     if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "DPR", strlen(ctmp), ctmp))) ERR(retval);
   }
   if ( have_cmb > 0 ) {
     sprintf(ctmp,"%s", cmb_file);
     if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "CMB", strlen(ctmp), ctmp))) ERR(retval);
   }

   sprintf(ctmp,"%f", clat);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "center_lat", strlen(ctmp), ctmp))) ERR(retval);
   sprintf(ctmp,"%f", clon);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "center_lon", strlen(ctmp), ctmp))) ERR(retval);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "center_date", strlen(cdate), cdate))) ERR(retval);

   sprintf(ctmp,"%f", clat1);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "beginning_lat", strlen(ctmp), ctmp))) ERR(retval);
   sprintf(ctmp,"%f", clon1);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "beginning_lon", strlen(ctmp), ctmp))) ERR(retval);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "beginning_date", strlen(cdate1), cdate1))) ERR(retval);

   sprintf(ctmp,"%f", clat2);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "ending_lat", strlen(ctmp), ctmp))) ERR(retval);
   sprintf(ctmp,"%f", clon2);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "ending_lon", strlen(ctmp), ctmp))) ERR(retval);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "ending_date", strlen(cdate2), cdate2))) ERR(retval);

   sprintf(ctmp,"%06d", irev);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "orbit_rev", strlen(ctmp), ctmp))) ERR(retval);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "satname", strlen(satname), satname))) ERR(retval);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "sensor", strlen(sensor), sensor))) ERR(retval);

   sprintf(ctmp,"%d", DB_MAX_INC_DIFF );
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "DB_MAX_INC_DIFF", strlen(ctmp), ctmp))) ERR(retval);

   sprintf(ctmp,"%d", DB_MAXREC);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "DB_MAXREC", strlen(ctmp), ctmp))) ERR(retval);

   sprintf(ctmp,"%d", DB_MINREC);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "DB_MINREC", strlen(ctmp), ctmp))) ERR(retval);

   sprintf(ctmp,"%d", MAX_DB_EXPAND);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "MAX_DB_EXPAND", strlen(ctmp), ctmp))) ERR(retval);

   sprintf(ctmp,"%f", DB_RAINFRAC);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "DB_RAINFRAC", strlen(ctmp), ctmp))) ERR(retval);

   sprintf(ctmp,"%d", N2MAX);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "N2MAX", strlen(ctmp), ctmp))) ERR(retval);

   sprintf(ctmp,"%f", WT2MIN);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "WT2MIN", strlen(ctmp), ctmp))) ERR(retval);

   sprintf(ctmp,"%d", MAX_TB_RMSD);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "MAX_TB_RMSD", strlen(ctmp), ctmp))) ERR(retval);

   sprintf(ctmp,"%d", NEM_USE);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "NEM_USE", strlen(ctmp), ctmp))) ERR(retval);

   sprintf(ctmp,"%d", NPCHIST);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "NPCHIST", strlen(ctmp), ctmp))) ERR(retval);

   if ( have_model > 0 ) {
     sprintf(ctmp,"%d", MAX_T2M_DIFF);
     if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "MAX_T2M_DIFF", strlen(ctmp), ctmp))) ERR(retval);
   }

   sprintf(ctmp,"%s", dbdir);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "EPC_DB_DIR", strlen(ctmp), ctmp))) ERR(retval);
   sprintf(ctmp,"%s", epc_coef_dir);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "EPC_COEF_DIR", strlen(ctmp), ctmp))) ERR(retval);
   sprintf(ctmp,"%d", NREG);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "NREG", strlen(ctmp), ctmp))) ERR(retval);
   sprintf(ctmp,"%d", isat_start);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "scan_start", strlen(ctmp), ctmp))) ERR(retval);
   sprintf(ctmp,"%d", isat_end);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "scan_end", strlen(ctmp), ctmp))) ERR(retval);

   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "processing_start", strlen(processing_date_start), processing_date_start))) ERR(retval);
   if ((retval = nc_put_att_text(ncid, NC_GLOBAL, "processing_end", strlen(processing_date_end), processing_date_end))) ERR(retval);
   if ((retval = nc_put_att(ncid, NC_GLOBAL, "tb_order", NC_STRING, NCHAN, tb_names))) ERR(status);
   if ((retval = nc_put_att(ncid, NC_GLOBAL, "epc_order", NC_STRING, NEM, em_names))) ERR(status);

   if ((retval = nc_def_dim(ncid, "nscan", nl_out, &scan_dimid))) ERR(retval);
   if ((retval = nc_def_dim(ncid, "npix", ns_S1, &pix_dimid))) ERR(retval);
   if ((retval = nc_def_dim(ncid, "ntb", nc_out, &tb_dimid))) ERR(retval);
   if ((retval = nc_def_dim(ncid, "nem", NEM, &em_dimid))) ERR(retval);
   if ((retval = nc_def_dim(ncid, "nlev", NLEV_DPR, &lev_dimid))) ERR(retval);
   if ((retval = nc_def_dim(ncid, "nlev_precip", NLEV_PRECIP, &lev_precip_dimid))) ERR(retval);

   if ((retval = nc_def_grp(ncid, "NS", &grpid_NS))) ERR (retval);
   if ((retval = nc_def_grp(ncid, "MS", &grpid_MS))) ERR (retval);
   if ( have_gprof > 0 ) if ((retval = nc_def_grp(ncid, "GPROF", &grpid1))) ERR (retval);
   if ( have_model > 0 ) if ((retval = nc_def_grp(ncid, "MERRA2", &grpid2))) ERR (retval);
   if ( have_dpr > 0 ) if ((retval = nc_def_grp(ncid, "DPR", &grpid3))) ERR (retval);
   if ( have_cmb > 0 ) if ((retval = nc_def_grp(ncid, "CMB", &grpid4))) ERR (retval);

   if ((retval = nc_def_var(ncid, "time", NC_DOUBLE, 1, &scan_dimid, &secs_varid))) ERR(retval);
   sprintf(units,"seconds since 1 Jan 1970");
   if ((retval = nc_put_att_text(ncid, secs_varid, "units", strlen(units), units))) ERR(retval);

   if ((retval = nc_def_var(ncid, "epc_compare", NC_SHORT, 1, &em_dimid, &nem_compare_varid))) ERR(retval);
   sprintf(units,"Set to 1 for EPC elements used in EPC-distance calculation");
   if ((retval = nc_put_att_text(ncid, nem_compare_varid, "units", strlen(units), units))) ERR(retval);

   if ((retval = nc_def_var(ncid, "tb_compare", NC_SHORT, 1, &tb_dimid, &ntb_compare_varid))) ERR(retval);
   sprintf(units,"Set to 1 for TB used in EPC calculation");
   if ((retval = nc_put_att_text(ncid, ntb_compare_varid, "units", strlen(units), units))) ERR(retval);

   dimids2[0] = scan_dimid;
   dimids2[1] = pix_dimid;

   /* latitude longitude */
   if ((retval = nc_def_var(ncid, "latitude", NC_FLOAT, 2, dimids2, &lat_varid))) ERR(retval);
   sprintf(units,"degrees");
   if ((retval = nc_put_att_text(ncid, lat_varid, "units", strlen(units), units))) ERR(retval);

   if ((retval = nc_def_var(ncid, "longitude", NC_FLOAT, 2, dimids2, &lon_varid))) ERR(retval);
   sprintf(units,"degrees");
   if ((retval = nc_put_att_text(ncid, lon_varid, "units", strlen(units), units))) ERR(retval);

   /* radiometer incidence angle   0=nadir */
   if ((retval = nc_def_var(ncid, "incidence_angle", NC_FLOAT, 2, dimids2, &inc_varid))) ERR(retval);
   sprintf(units,"S1 incidence angle, in degrees");
   if ((retval = nc_put_att_text(ncid, inc_varid, "units", strlen(units), units))) ERR(retval);

   /* a-priori database index */
   if ((retval = nc_def_var(ncid, "db_index", NC_SHORT, 2, dimids2, &pixel_index_varid))) ERR(retval);
   sprintf(units,"DB index, ranges between 0 and (NPCHIST^NEM_USE)-1");
   if ((retval = nc_put_att_text(ncid, pixel_index_varid, "units", strlen(units), units))) ERR(retval);

   /*---- estimates with the NS database ----*/

   /* mean precip based on combined */
   if ((retval = nc_def_var(grpid_NS, "precip", NC_FLOAT, 2, dimids2, &precip1_NS_varid))) ERR(retval);
   sprintf(units,"estimate using DPR+GMI combined estimate, in mm/hr");
   if ((retval = nc_put_att_text(grpid_NS, precip1_NS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean precip based on radar-only */
   if ((retval = nc_def_var(grpid_NS, "precip2", NC_FLOAT, 2, dimids2, &precip2_NS_varid))) ERR(retval);
   sprintf(units,"estimate using DPR estimate, in mm/hr");
   if ((retval = nc_put_att_text(grpid_NS, precip2_NS_varid, "units", strlen(units), units))) ERR(retval);

   /* prob of precip */
   if ((retval = nc_def_var(grpid_NS, "precip_prob", NC_FLOAT, 2, dimids2, &precip1_prob_NS_varid))) ERR(retval);
   sprintf(units,"probability of precip for the DPR+GMI combined estimate, in percent");
   if ((retval = nc_put_att_text(grpid_NS, precip1_prob_NS_varid, "units", strlen(units), units))) ERR(retval);

   /* prob of precip2 */
   if ((retval = nc_def_var(grpid_NS, "precip2_prob", NC_FLOAT, 2, dimids2, &precip2_prob_NS_varid))) ERR(retval);
   sprintf(units,"probability of precip for the DPR estimate, in percent");
   if ((retval = nc_put_att_text(grpid_NS, precip2_prob_NS_varid, "units", strlen(units), units))) ERR(retval);

   /* quality flag */
   if ((retval = nc_def_var(grpid_NS, "quality", NC_INT, 2, dimids2, &quality_NS_varid))) ERR(retval);
   sprintf(units,"0=good 1=bad cooords 2=no DBfile located");
   if ((retval = nc_put_att_text(grpid_NS, quality_NS_varid, "units", strlen(units), units))) ERR(retval);

   /* type precip */
   if ((retval = nc_def_var(grpid_NS, "type_precip", NC_SHORT, 2, dimids2, &type_precip_NS_varid))) ERR(retval);
   sprintf(units,"type precip, 0=convective 1=stratiform 2=other");
   if ((retval = nc_put_att_text(grpid_NS, type_precip_NS_varid, "units", strlen(units), units))) ERR(retval);

   /* shallow rain */
   if ((retval = nc_def_var(grpid_NS, "shallow_rain", NC_SHORT, 2, dimids2, &shallow_rain_NS_varid))) ERR(retval);
   sprintf(units,"shallow rain flag, 0=no shallow 1=shallow");
   if ((retval = nc_put_att_text(grpid_NS, shallow_rain_NS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean Ts */
   if ((retval = nc_def_var(grpid_NS, "ts", NC_FLOAT, 2, dimids2, &ts_NS_varid))) ERR(retval);
   sprintf(units,"estimated surface temperature, in Kelvin");
   if ((retval = nc_put_att_text(grpid_NS, ts_NS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean T2m */
   if ((retval = nc_def_var(grpid_NS, "t2m", NC_FLOAT, 2, dimids2, &t2m_NS_varid))) ERR(retval);
   sprintf(units,"estimated 2-m air temperature, in Kelvin");
   if ((retval = nc_put_att_text(grpid_NS, t2m_NS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean column vapor  */
   if ((retval = nc_def_var(grpid_NS, "tqv", NC_FLOAT, 2, dimids2, &tqv_NS_varid))) ERR(retval);
   sprintf(units,"estimated total column water vapor, in mm");
   if ((retval = nc_put_att_text(grpid_NS, tqv_NS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean height 273K */
   if ((retval = nc_def_var(grpid_NS, "h273", NC_SHORT, 2, dimids2, &h273_NS_varid))) ERR(retval);
   sprintf(units,"estimated freezing level height, in meters");
   if ((retval = nc_put_att_text(grpid_NS, h273_NS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean height Ku-band 20 db level */
   if ((retval = nc_def_var(grpid_NS, "h20dB_ku", NC_SHORT, 2, dimids2, &ku_ztop_NS_varid))) ERR(retval);
   sprintf(units,"estimated 20 dB height from Ku-band, in meters");
   if ((retval = nc_put_att_text(grpid_NS, ku_ztop_NS_varid, "units", strlen(units), units))) ERR(retval);

   /*---- estimates with the MS database ----*/

   /* mean precip based on combined */
   if ((retval = nc_def_var(grpid_MS, "precip", NC_FLOAT, 2, dimids2, &precip1_MS_varid))) ERR(retval);
   sprintf(units,"estimate using DPR+GMI combined estimate, in mm/hr");
   if ((retval = nc_put_att_text(grpid_MS, precip1_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean precip based on radar-only */
   if ((retval = nc_def_var(grpid_MS, "precip2", NC_FLOAT, 2, dimids2, &precip2_MS_varid))) ERR(retval);
   sprintf(units,"estimate using DPR estimate, in mm/hr");
   if ((retval = nc_put_att_text(grpid_MS, precip2_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* prob of precip */
   if ((retval = nc_def_var(grpid_MS, "precip_prob", NC_FLOAT, 2, dimids2, &precip1_prob_MS_varid))) ERR(retval);
   sprintf(units,"probability of precip for the DPR+GMI combined estimate, in percent");
   if ((retval = nc_put_att_text(grpid_MS, precip1_prob_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* prob of precip2 */
   if ((retval = nc_def_var(grpid_MS, "precip2_prob", NC_FLOAT, 2, dimids2, &precip2_prob_MS_varid))) ERR(retval);
   sprintf(units,"probability of precip for the DPR estimate, in percent");
   if ((retval = nc_put_att_text(grpid_MS, precip2_prob_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* quality flag */
   if ((retval = nc_def_var(grpid_MS, "quality", NC_INT, 2, dimids2, &quality_MS_varid))) ERR(retval);
   sprintf(units,"0=good 1=bad cooords 2=no DBfile located");
   if ((retval = nc_put_att_text(grpid_MS, quality_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* type precip */
   if ((retval = nc_def_var(grpid_MS, "type_precip", NC_SHORT, 2, dimids2, &type_precip_MS_varid))) ERR(retval);
   sprintf(units,"type precip, 0=convective 1=stratiform 2=other");
   if ((retval = nc_put_att_text(grpid_MS, type_precip_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* shallow rain */
   if ((retval = nc_def_var(grpid_MS, "shallow_rain", NC_SHORT, 2, dimids2, &shallow_rain_MS_varid))) ERR(retval);
   sprintf(units,"shallow rain flag, 0=no shallow 1=shallow");
   if ((retval = nc_put_att_text(grpid_MS, shallow_rain_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean Ts */
   if ((retval = nc_def_var(grpid_MS, "ts", NC_FLOAT, 2, dimids2, &ts_MS_varid))) ERR(retval);
   sprintf(units,"estimated surface temperature, in Kelvin");
   if ((retval = nc_put_att_text(grpid_MS, ts_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean T2m */
   if ((retval = nc_def_var(grpid_MS, "t2m", NC_FLOAT, 2, dimids2, &t2m_MS_varid))) ERR(retval);
   sprintf(units,"estimated 2-m air temperature, in Kelvin");
   if ((retval = nc_put_att_text(grpid_MS, t2m_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean column vapor  */
   if ((retval = nc_def_var(grpid_MS, "tqv", NC_FLOAT, 2, dimids2, &tqv_MS_varid))) ERR(retval);
   sprintf(units,"estimated total column water vapor, in mm");
   if ((retval = nc_put_att_text(grpid_MS, tqv_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean height 273K */
   if ((retval = nc_def_var(grpid_MS, "h273", NC_SHORT, 2, dimids2, &h273_MS_varid))) ERR(retval);
   sprintf(units,"estimated freezing level height, in meters");
   if ((retval = nc_put_att_text(grpid_MS, h273_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean height Ku-band 20 db level */
   if ((retval = nc_def_var(grpid_MS, "h20dB_ku", NC_SHORT, 2, dimids2, &ku_ztop_MS_varid))) ERR(retval);
   sprintf(units,"estimated 20 dB height from Ku-band, in meters");
   if ((retval = nc_put_att_text(grpid_MS, ku_ztop_MS_varid, "units", strlen(units), units))) ERR(retval);

   /* mean height Ka-band 20 db level */
   if ((retval = nc_def_var(grpid_MS, "h20dB_ka", NC_SHORT, 2, dimids2, &ka_ztop_MS_varid))) ERR(retval);
   sprintf(units,"estimated 20 dB height from Ka-band, in meters");
   if ((retval = nc_put_att_text(grpid_MS, ka_ztop_MS_varid, "units", strlen(units), units))) ERR(retval);


   if ( have_gprof > 0 ) {
     /* GPROF precip */
     if ((retval = nc_def_var(grpid1, "precip", NC_FLOAT, 2, dimids2, &precip_GPROF_varid))) ERR(retval);
     sprintf(units,"precip from GPROF, in mm/hr");
     if ((retval = nc_put_att_text(grpid1, precip_GPROF_varid, "units", strlen(units), units))) ERR(retval);

     /* GPROF prob of precip */
     if ((retval = nc_def_var(grpid1, "precip_prob", NC_FLOAT, 2, dimids2, &precip_prob_GPROF_varid))) ERR(retval);
     sprintf(units,"probability of precip from GPROF, in percent");
     if ((retval = nc_put_att_text(grpid1, precip_prob_GPROF_varid, "units", strlen(units), units))) ERR(retval);

     /* GPROF T2m index */
     if ((retval = nc_def_var(grpid1, "t2m_index", NC_SHORT, 2, dimids2, &t2m_GPROF_varid))) ERR(retval);
     sprintf(units,"2-m air temperature index from GPROF, in Kelvin");
     if ((retval = nc_put_att_text(grpid1, t2m_GPROF_varid, "units", strlen(units), units))) ERR(retval);

     /* GPROF column vapor index */
     if ((retval = nc_def_var(grpid1, "tqv_index", NC_BYTE, 2, dimids2, &tqv_GPROF_varid))) ERR(retval);
     sprintf(units,"total column water vapor index from GPROF, in mm");
     if ((retval = nc_put_att_text(grpid1, tqv_GPROF_varid, "units", strlen(units), units))) ERR(retval);

     /* GPROF frozen precip */
     if ((retval = nc_def_var(grpid1, "frozen_precip", NC_FLOAT, 2, dimids2, &frozen_precip_GPROF_varid))) ERR(retval);
     sprintf(units,"frozen precip from GPROF, in mm/hr");
     if ((retval = nc_put_att_text(grpid1, frozen_precip_GPROF_varid, "units", strlen(units), units))) ERR(retval);

     /* GPROF surface class */
     if ((retval = nc_def_var(grpid1, "sfc_class", NC_BYTE, 2, dimids2, &sfc_GPROF_varid))) ERR(retval);
     sprintf(units,"TELSEM surface class index");
     if ((retval = nc_put_att_text(grpid1, sfc_GPROF_varid, "units", strlen(units), units))) ERR(retval);
   }

   if ( have_model > 0 ) {
     /* actual MERRA2 Ts */
     if ((retval = nc_def_var(grpid2, "ts", NC_FLOAT, 2, dimids2, &ts_MERRA2_varid))) ERR(retval);
     sprintf(units,"MERRA2 surface temperature, in Kelvin");
     if ((retval = nc_put_att_text(grpid2, ts_MERRA2_varid, "units", strlen(units), units))) ERR(retval);

     /* actual MERRA2 T2m */
     if ((retval = nc_def_var(grpid2, "t2m", NC_FLOAT, 2, dimids2, &t2m_MERRA2_varid))) ERR(retval);
     sprintf(units,"MERRA2 2-m air temperature, in Kelvin");
     if ((retval = nc_put_att_text(grpid2, t2m_MERRA2_varid, "units", strlen(units), units))) ERR(retval);

     /* actual MERRA2 T2m dewpoint */
     if ((retval = nc_def_var(grpid2, "t2m_dew", NC_FLOAT, 2, dimids2, &t2m_dew_MERRA2_varid))) ERR(retval);
     sprintf(units,"MERRA2 2-m dewpoint temperature, in Kelvin");
     if ((retval = nc_put_att_text(grpid2, t2m_dew_MERRA2_varid, "units", strlen(units), units))) ERR(retval);

     /* actual MERRA2 T2m wet bulb*/
     if ((retval = nc_def_var(grpid2, "t2m_wet", NC_FLOAT, 2, dimids2, &t2m_wet_MERRA2_varid))) ERR(retval);
     sprintf(units,"MERRA2 2-m wet bulb temperature, in Kelvin");
     if ((retval = nc_put_att_text(grpid2, t2m_wet_MERRA2_varid, "units", strlen(units), units))) ERR(retval);

     /* actual MERRA2 column vapor  */
     if ((retval = nc_def_var(grpid2, "tqv", NC_FLOAT, 2, dimids2, &tqv_MERRA2_varid))) ERR(retval);
     sprintf(units,"MERRA2 total column water vapor, in mm");
     if ((retval = nc_put_att_text(grpid2, tqv_MERRA2_varid, "units", strlen(units), units))) ERR(retval);

     /* actual MERRA2 273K level */
     if ((retval = nc_def_var(grpid2, "h273", NC_SHORT, 2, dimids2, &h273_MERRA2_varid))) ERR(retval);
     sprintf(units,"MERRA2 freezing level, in meters");
     if ((retval = nc_put_att_text(grpid2, h273_MERRA2_varid, "units", strlen(units), units))) ERR(retval);

   }

   if ( have_dpr > 0 ) {
     if ((retval = nc_def_var(grpid3, "latitude", NC_FLOAT, 2, dimids2, &lat_NS_DPR_varid))) ERR(retval);
     sprintf(units,"latitude of nearest DPR NS pixel, in degrees");
     if ((retval = nc_put_att_text(grpid3, lat_NS_DPR_varid, "units", strlen(units), units))) ERR(retval);

     if ((retval = nc_def_var(grpid3, "longitude", NC_FLOAT, 2, dimids2, &lon_NS_DPR_varid))) ERR(retval);
     sprintf(units,"longitude of nearest DPR NS pixel, in degrees");
     if ((retval = nc_put_att_text(grpid3, lon_NS_DPR_varid, "units", strlen(units), units))) ERR(retval);

     /* actual NS precip */
     if ((retval = nc_def_var(grpid3, "precip_NS", NC_FLOAT, 2, dimids2, &precip_NS_DPR_varid))) ERR(retval);
     sprintf(units,"DPR NS precip, in mm/hr");
     if ((retval = nc_put_att_text(grpid3, precip_NS_DPR_varid, "units", strlen(units), units))) ERR(retval);

     /* actual MS precip */
     if ((retval = nc_def_var(grpid3, "precip_MS", NC_FLOAT, 2, dimids2, &precip_MS_DPR_varid))) ERR(retval);
     sprintf(units,"DPR MS precip, in mm/hr");
     if ((retval = nc_put_att_text(grpid3, precip_MS_DPR_varid, "units", strlen(units), units))) ERR(retval);

     /* elevation from DPR */
     if ((retval = nc_def_var(grpid3, "elevation", NC_SHORT, 2, dimids2, &elev_DPR_varid))) ERR(retval);
     sprintf(units,"elevation from DPR, in meters");
     if ((retval = nc_put_att_text(grpid3, elev_DPR_varid, "units", strlen(units), units))) ERR(retval);

     /* height DPR Ku-band 20 db level */
     if ((retval = nc_def_var(grpid3, "h20dB_ku", NC_SHORT, 2, dimids2, &ztop_ku_DPR_varid))) ERR(retval);
     sprintf(units,"20 dB height from DPR Ku-band, in meters");
     if ((retval = nc_put_att_text(grpid3, ztop_ku_DPR_varid, "units", strlen(units), units))) ERR(retval);

     /* height DPR Ka-band 20 db level */
     if ((retval = nc_def_var(grpid3, "h20dB_ka", NC_SHORT, 2, dimids2, &ztop_ka_DPR_varid))) ERR(retval);
     sprintf(units,"20 dB height from DPR Ka-band, in meters");
     if ((retval = nc_put_att_text(grpid3, ztop_ka_DPR_varid, "units", strlen(units), units))) ERR(retval);

     /* NS type precip */
     if ((retval = nc_def_var(grpid3, "type_precip_NS", NC_SHORT, 2, dimids2, &type_precip_DPR_varid))) ERR(retval);
     sprintf(units,"DPR NS typePrecip, 0=convective 1=stratiform 2=other");
     if ((retval = nc_put_att_text(grpid3, type_precip_DPR_varid, "units", strlen(units), units))) ERR(retval);

     /* NS shallow rain */
     if ((retval = nc_def_var(grpid3, "shallow_rain_NS", NC_SHORT, 2, dimids2, &shallow_rain_DPR_varid))) ERR(retval);
     sprintf(units,"DPR NS flagShallowRain, 0=no shallow 1=shallow");
     if ((retval = nc_put_att_text(grpid3, shallow_rain_DPR_varid, "units", strlen(units), units))) ERR(retval);
   }

   if ( have_cmb > 0 ) {
     /* actual NS precip */
     if ((retval = nc_def_var(grpid4, "precip_NS", NC_FLOAT, 2, dimids2, &precip_NS_CMB_varid))) ERR(retval);
     sprintf(units,"DPR+GMI combined NS precip, in mm/hr");
     if ((retval = nc_put_att_text(grpid4, precip_NS_CMB_varid, "units", strlen(units), units))) ERR(retval);

     /* actual MS precip */
     if ((retval = nc_def_var(grpid4, "precip_MS", NC_FLOAT, 2, dimids2, &precip_MS_CMB_varid))) ERR(retval);
     sprintf(units,"DPR+GMI combined MS precip, in mm/hr");
     if ((retval = nc_put_att_text(grpid4, precip_MS_CMB_varid, "units", strlen(units), units))) ERR(retval);
   }



   dimids3[0] = scan_dimid;
   dimids3[1] = pix_dimid;
   dimids3[2] = tb_dimid;

   dimids4[0] = scan_dimid;
   dimids4[1] = pix_dimid;
   dimids4[2] = em_dimid;

   dimids5[0] = scan_dimid;
   dimids5[1] = pix_dimid;
   dimids5[2] = lev_dimid;

   dimids6[0] = scan_dimid;
   dimids6[1] = pix_dimid;
   dimids6[2] = lev_precip_dimid;

   /* TB */
   if ((retval = nc_def_var(ncid, "Tb", NC_FLOAT, 3, dimids3, &tb_varid))) ERR(retval);
   sprintf(units,"Kelvin");
   if ((retval = nc_put_att_text(ncid, tb_varid, "units", strlen(units), units))) ERR(retval);


   /* NS precip profile estimates */
   if ((retval = nc_def_var(grpid_NS, "precip_prof", NC_SHORT, 3, dimids6, &precip1_prof_NS_varid))) ERR(retval);
   sprintf(units,"precip profile estimate using DPR+GMI combined estimate, in mm/hr scaled by 100");
   if ((retval = nc_put_att_text(grpid_NS, precip1_prof_NS_varid, "units", strlen(units), units))) ERR(retval);

   if ((retval = nc_def_var(grpid_NS, "precip2_prof", NC_SHORT, 3, dimids6, &precip2_prof_NS_varid))) ERR(retval);
   sprintf(units,"precip profile estimate using DPR estimate, in mm/hr scaled by 100");
   if ((retval = nc_put_att_text(grpid_NS, precip2_prof_NS_varid, "units", strlen(units), units))) ERR(retval);

   /* MS precip profile estimates */
   if ((retval = nc_def_var(grpid_MS, "precip_prof", NC_SHORT, 3, dimids6, &precip1_prof_MS_varid))) ERR(retval);
   sprintf(units,"precip profile estimate using DPR+GMI combined estimate, in mm/hr scaled by 100");
   if ((retval = nc_put_att_text(grpid_MS, precip1_prof_MS_varid, "units", strlen(units), units))) ERR(retval);

   if ((retval = nc_def_var(grpid_MS, "precip2_prof", NC_SHORT, 3, dimids6, &precip2_prof_MS_varid))) ERR(retval);
   sprintf(units,"precip profile estimate using DPR estimate, in mm/hr scaled by 100");
   if ((retval = nc_put_att_text(grpid_MS, precip2_prof_MS_varid, "units", strlen(units), units))) ERR(retval);


   if ( save_pc > 0 ) {
     /* reconstructed EMIS vector everywhere */
     if ((retval = nc_def_var(ncid, "emis", NC_FLOAT, 3, dimids4, &emis_varid))) ERR(retval);
     sprintf(units,"- Kelvin mm Kelvin");
     if ((retval = nc_put_att_text(ncid, emis_varid, "units", strlen(units), units))) ERR(retval);

     /* EPC vector everywhere */
     if ((retval = nc_def_var(ncid, "pc_emis", NC_FLOAT, 3, dimids4, &pc_emis_varid))) ERR(retval);
     sprintf(units,"principal components of the emis vector");
     if ((retval = nc_put_att_text(ncid, pc_emis_varid, "units", strlen(units), units))) ERR(retval);
   }

   if ( save_top_ranked > 0 ) {
     /* top candidate Ku profile from NS estimate */
     if ((retval = nc_def_var(grpid_NS, "z_ku_1", NC_SHORT, 3, dimids5, &z_ku_NS_varid))) ERR(retval);
     sprintf(units,"ku-band DPR profile from top-ranked candidate, dB scaled by 100");
     if ((retval = nc_put_att_text(grpid_NS, z_ku_NS_varid, "units", strlen(units), units))) ERR(retval);

     /* top candidate TB from NS estimate */
     if ((retval = nc_def_var(grpid_NS, "Tb_1", NC_FLOAT, 3, dimids3, &tb1_NS_varid))) ERR(retval);
     sprintf(units,"TB from top-ranked profile, Kelvin");
     if ((retval = nc_put_att_text(grpid_NS, tb1_NS_varid, "units", strlen(units), units))) ERR(retval);

     /* top candidate Ku profile from MS estimate */
     if ((retval = nc_def_var(grpid_MS, "z_ku_1", NC_SHORT, 3, dimids5, &z_ku_MS_varid))) ERR(retval);
     sprintf(units,"DPR ku-band profile from top-ranked candidate, dB scaled by 100");
     if ((retval = nc_put_att_text(grpid_MS, z_ku_MS_varid, "units", strlen(units), units))) ERR(retval);

     /* top candidate Ka profile from MS estimate */
     if ((retval = nc_def_var(grpid_MS, "z_ka_1", NC_SHORT, 3, dimids5, &z_ka_MS_varid))) ERR(retval);
     sprintf(units,"DPR ka-band profile from top-ranked candidate, dB scaled by 100");
     if ((retval = nc_put_att_text(grpid_MS, z_ka_MS_varid, "units", strlen(units), units))) ERR(retval);

     /* top candidate TB from MS estimate */
     if ((retval = nc_def_var(grpid_MS, "Tb_1", NC_FLOAT, 3, dimids3, &tb1_MS_varid))) ERR(retval);
     sprintf(units,"TB from top-ranked profile, Kelvin");
     if ((retval = nc_put_att_text(grpid_MS, tb1_MS_varid, "units", strlen(units), units))) ERR(retval);
   }

   /* Actual DPR Ku and Ka profiles */
   if ( save_dpr_profiles > 0 && have_dpr > 0 ) {
     if ((retval = nc_def_var(grpid3, "z_ku", NC_SHORT, 3, dimids5, &z_ku_DPR_varid))) ERR(retval);
     sprintf(units,"ku-band DPR profile, dB scaled by 100");
     if ((retval = nc_put_att_text(grpid3, z_ku_DPR_varid, "units", strlen(units), units))) ERR(retval);

     if ((retval = nc_def_var(grpid3, "z_ka", NC_SHORT, 3, dimids5, &z_ka_DPR_varid))) ERR(retval);
     sprintf(units,"ka-band DPR profile, dB scaled by 100");
     if ((retval = nc_put_att_text(grpid3, z_ka_DPR_varid, "units", strlen(units), units))) ERR(retval);
   }

   /* Actual CMB precip profile */
   if ( save_dpr_profiles > 0 && have_cmb > 0 ) {
     if ((retval = nc_def_var(grpid4, "precip_prof_NS", NC_SHORT, 3, dimids6, &precip_prof_NS_CMB_varid))) ERR(retval);
     sprintf(units,"DPR+GMI combined NS precip profile, in mm/hr, scaled by 100");
     if ((retval = nc_put_att_text(grpid4, precip_prof_NS_CMB_varid, "units", strlen(units), units))) ERR(retval);

     if ((retval = nc_def_var(grpid4, "precip_prof_MS", NC_SHORT, 3, dimids6, &precip_prof_MS_CMB_varid))) ERR(retval);
     sprintf(units,"DPR+GMI combined MS precip profile, in mm/hr, scaled by 100");
     if ((retval = nc_put_att_text(grpid4, precip_prof_MS_CMB_varid, "units", strlen(units), units))) ERR(retval);
   }

   /* on the fly compression for the 3-D vars */
   if ((status= nc_def_var_deflate(ncid, tb_varid, 0, 1, 4))) ERR(status);
   if ( save_pc > 0 ) {
     if ((status= nc_def_var_deflate(ncid, emis_varid, 0, 1, 4))) ERR(status);
     if ((status= nc_def_var_deflate(ncid, pc_emis_varid, 0, 1, 4))) ERR(status);
   }
   if ( save_top_ranked > 0 ) {
     if ((status= nc_def_var_deflate(grpid_NS, z_ku_NS_varid, 0, 1, 4))) ERR(status);
     if ((status= nc_def_var_deflate(grpid_NS, tb1_NS_varid, 0, 1, 4))) ERR(status);
     if ((status= nc_def_var_deflate(grpid_MS, z_ku_MS_varid, 0, 1, 4))) ERR(status);
     if ((status= nc_def_var_deflate(grpid_MS, z_ka_MS_varid, 0, 1, 4))) ERR(status);
     if ((status= nc_def_var_deflate(grpid_MS, tb1_MS_varid, 0, 1, 4))) ERR(status);
   }
   if ( save_dpr_profiles > 0 && have_dpr > 0 ) {
     if ((status= nc_def_var_deflate(grpid3, z_ku_DPR_varid, 0, 1, 4))) ERR(status);
     if ((status= nc_def_var_deflate(grpid3, z_ka_DPR_varid, 0, 1, 4))) ERR(status);
   }


   if ((retval = nc_enddef(ncid))) ERR(retval);


   n= 0;
   for (isat= isat_start; isat <= isat_end; isat++) {
     dbuf[n]= secs_S1_buf[isat];
     n++;
   }
   if ((retval = nc_put_var_double(ncid, secs_varid, dbuf))) ERR(retval);

   if ((retval = nc_put_var_short(ncid, nem_compare_varid, em_compare))) ERR(retval);
   if ((retval = nc_put_var_short(ncid, ntb_compare_varid, tb_compare))) ERR(retval);

   count2[0] = nl_out;
   count2[1] = ns_S1;
   start2[0] = 0;
   start2[1] = 0;

   n= 0;
   for (isat= isat_start; isat <= isat_end; isat++) {
     for (jsat= 0; jsat < ns_S1; jsat++) {
       idx= isat*ns_S1 + jsat;
       idx1= n*ns_S1 + jsat;
       fbuf[idx1]= lat_S1_buf[idx];
     }
     n++;
   }
   if ((retval = nc_put_vara_float(ncid, lat_varid, start2, count2, fbuf))) ERR(retval);

   n= 0;
   for (isat= isat_start; isat <= isat_end; isat++) {
     for (jsat= 0; jsat < ns_S1; jsat++) {
       idx= isat*ns_S1 + jsat;
       idx1= n*ns_S1 + jsat;
       fbuf[idx1]= lon_S1_buf[idx];
     }
     n++;
   }
   if ((retval = nc_put_vara_float(ncid, lon_varid, start2, count2, fbuf))) ERR(retval);

   n= 0;
   for (isat= isat_start; isat <= isat_end; isat++) {
     for (jsat= 0; jsat < ns_S1; jsat++) {
       idx= isat*ns_S1 + jsat;
       idx1= n*ns_S1 + jsat;
       fbuf[idx1]= incidenceAngle_S1_buf[idx];
     }
     n++;
   }
   if ((retval = nc_put_vara_float(ncid, inc_varid, start2, count2, fbuf))) ERR(retval);

   n= 0;
   for (isat= isat_start; isat <= isat_end; isat++) {
     for (jsat= 0; jsat < ns_S1; jsat++) {
       idx= isat*ns_S1 + jsat;
       idx1= n*ns_S1 + jsat;
       sbuf[idx1]= pixel_index[idx];
     }
     n++;
   }
   if ((retval = nc_put_vara_short(ncid, pixel_index_varid, start2, count2, sbuf))) ERR(retval);

   /* estimated precip and precip2 for the NS and MS databases */

   if ((retval = nc_put_vara_float(grpid_NS, precip1_NS_varid, start2, count2, precip1_NS_buf))) ERR(retval);
   if ((retval = nc_put_vara_float(grpid_NS, precip2_NS_varid, start2, count2, precip2_NS_buf))) ERR(retval);
   if ((retval = nc_put_vara_float(grpid_NS, precip1_prob_NS_varid, start2, count2, precip1_prob_NS_buf))) ERR(retval);
   if ((retval = nc_put_vara_float(grpid_NS, precip2_prob_NS_varid, start2, count2, precip2_prob_NS_buf))) ERR(retval);
   if ((retval = nc_put_vara_int(grpid_NS, quality_NS_varid, start2, count2, epc_qual_NS_buf))) ERR(retval);
   if ((retval = nc_put_vara_short(grpid_NS, type_precip_NS_varid, start2, count2, type_precip_NS_buf))) ERR(retval);
   if ((retval = nc_put_vara_short(grpid_NS, shallow_rain_NS_varid, start2, count2, shallow_rain_NS_buf))) ERR(retval);

   if ((retval = nc_put_vara_float(grpid_MS, precip1_MS_varid, start2, count2, precip1_MS_buf))) ERR(retval);
   if ((retval = nc_put_vara_float(grpid_MS, precip2_MS_varid, start2, count2, precip2_MS_buf))) ERR(retval);
   if ((retval = nc_put_vara_float(grpid_MS, precip1_prob_MS_varid, start2, count2, precip1_prob_MS_buf))) ERR(retval);
   if ((retval = nc_put_vara_float(grpid_MS, precip2_prob_MS_varid, start2, count2, precip2_prob_MS_buf))) ERR(retval);
   if ((retval = nc_put_vara_int(grpid_MS, quality_MS_varid, start2, count2, epc_qual_MS_buf))) ERR(retval);
   if ((retval = nc_put_vara_short(grpid_MS, type_precip_MS_varid, start2, count2, type_precip_MS_buf))) ERR(retval);
   if ((retval = nc_put_vara_short(grpid_MS, shallow_rain_MS_varid, start2, count2, shallow_rain_NS_buf))) ERR(retval);

   /* estimated ts, t2m, tqv, h273, Ztop */

   if ((retval = nc_put_vara_float(grpid_NS, ts_NS_varid, start2, count2, ts_epc_NS_buf))) ERR(retval);
   if ((retval = nc_put_vara_float(grpid_NS, t2m_NS_varid, start2, count2, t2m_epc_NS_buf))) ERR(retval);
   if ((retval = nc_put_vara_float(grpid_NS, tqv_NS_varid, start2, count2, tqv_epc_NS_buf))) ERR(retval);
   if ((retval = nc_put_vara_short(grpid_NS, h273_NS_varid, start2, count2, h273_epc_NS_buf))) ERR(retval);
   if ((retval = nc_put_vara_short(grpid_NS, ku_ztop_NS_varid, start2, count2, ku_ztop_epc_NS_buf))) ERR(retval);

   if ((retval = nc_put_vara_float(grpid_MS, ts_MS_varid, start2, count2, ts_epc_MS_buf))) ERR(retval);
   if ((retval = nc_put_vara_float(grpid_MS, t2m_MS_varid, start2, count2, t2m_epc_MS_buf))) ERR(retval);
   if ((retval = nc_put_vara_float(grpid_MS, tqv_MS_varid, start2, count2, tqv_epc_MS_buf))) ERR(retval);
   if ((retval = nc_put_vara_short(grpid_MS, h273_MS_varid, start2, count2, h273_epc_MS_buf))) ERR(retval);
   if ((retval = nc_put_vara_short(grpid_MS, ku_ztop_MS_varid, start2, count2, ku_ztop_epc_MS_buf))) ERR(retval);
   if ((retval = nc_put_vara_short(grpid_MS, ka_ztop_MS_varid, start2, count2, ka_ztop_epc_MS_buf))) ERR(retval);

   /* GPROF variables if provided */

   if ( have_gprof > 0 ) {
     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx1= n*ns_S1 + jsat;
         fbuf[idx1]= surfacePrecipitation_GPROF_S1_buf[idx];
       }
       n++;
     }
     if ((retval = nc_put_vara_float(grpid1, precip_GPROF_varid, start2, count2, fbuf))) ERR(retval);

     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx1= n*ns_S1 + jsat;
         fbuf[idx1]= probabilityOfPrecip_GPROF_S1_buf[idx];
       }
       n++;
     }
     if ((retval = nc_put_vara_float(grpid1, precip_prob_GPROF_varid, start2, count2, fbuf))) ERR(retval);

     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx1= n*ns_S1 + jsat;
         sbuf[idx1]= temp2mIndex_GPROF_S1_buf[idx];
       }
       n++;
     }
     if ((retval = nc_put_vara_short(grpid1, t2m_GPROF_varid, start2, count2, sbuf))) ERR(retval);

     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx1= n*ns_S1 + jsat;
         bbuf[idx1]= totalColumnWaterVaporIndex_GPROF_S1_buf[idx];
       }
       n++;
     }
     if ((retval = nc_put_vara_schar(grpid1, tqv_GPROF_varid, start2, count2, bbuf))) ERR(retval);

     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx1= n*ns_S1 + jsat;
         fbuf[idx1]= frozenPrecipitation_GPROF_S1_buf[idx];
       }
       n++;
     }
     if ((retval = nc_put_vara_float(grpid1, frozen_precip_GPROF_varid, start2, count2, fbuf))) ERR(retval);

     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx1= n*ns_S1 + jsat;
         bbuf[idx1]= sfc_GPROF_S1_buf[idx];
       }
       n++;
     }
     if ((retval = nc_put_vara_schar(grpid1, sfc_GPROF_varid, start2, count2, bbuf))) ERR(retval);

   }

   /* Actual DPR variables, if provided */

   if ( have_dpr > 0 ) {
     if ((retval = nc_put_vara_float(grpid3, lat_NS_DPR_varid, start2, count2, lat_NS_DPR_buf))) ERR(retval);
     if ((retval = nc_put_vara_float(grpid3, lon_NS_DPR_varid, start2, count2, lon_NS_DPR_buf))) ERR(retval);
     if ((retval = nc_put_vara_float(grpid3, precip_NS_DPR_varid, start2, count2, precip_NS_DPR_buf))) ERR(retval);
     if ((retval = nc_put_vara_float(grpid3, precip_MS_DPR_varid, start2, count2, precip_MS_DPR_buf))) ERR(retval);
     if ((retval = nc_put_vara_short(grpid3, elev_DPR_varid, start2, count2, elev_DPR_buf))) ERR(retval);
     if ((retval = nc_put_vara_short(grpid3, ztop_ku_DPR_varid, start2, count2, ztop_ku_DPR_buf))) ERR(retval);
     if ((retval = nc_put_vara_short(grpid3, ztop_ka_DPR_varid, start2, count2, ztop_ka_DPR_buf))) ERR(retval);
     if ((retval = nc_put_vara_short(grpid3, type_precip_DPR_varid, start2, count2, type_precip_DPR_buf))) ERR(retval);
     if ((retval = nc_put_vara_short(grpid3, shallow_rain_DPR_varid, start2, count2, shallow_rain_DPR_buf))) ERR(retval);
   }
   if ( have_cmb > 0 ) {
     if ((retval = nc_put_vara_float(grpid4, precip_NS_CMB_varid, start2, count2, precip_NS_CMB_buf))) ERR(retval);
     if ((retval = nc_put_vara_float(grpid4, precip_MS_CMB_varid, start2, count2, precip_MS_CMB_buf))) ERR(retval);
   }

   /* Actual MERRA2 variables, if provided */

   if ( have_model > 0 ) {
     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx1= n*ns_S1 + jsat;
         fbuf[idx1]= ts_buf[idx];
       }
       n++;
     }
     if ((retval = nc_put_vara_float(grpid2, ts_MERRA2_varid, start2, count2, fbuf))) ERR(retval);

     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx1= n*ns_S1 + jsat;
         fbuf[idx1]= t2m_buf[idx];
       }
       n++;
     }
     if ((retval = nc_put_vara_float(grpid2, t2m_MERRA2_varid, start2, count2, fbuf))) ERR(retval);

     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx1= n*ns_S1 + jsat;
         fbuf[idx1]= t2m_wet_buf[idx];
       }
       n++;
     }
     if ((retval = nc_put_vara_float(grpid2, t2m_wet_MERRA2_varid, start2, count2, fbuf))) ERR(retval);

     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx1= n*ns_S1 + jsat;
         fbuf[idx1]= t2m_dew_buf[idx];
       }
       n++;
     }
     if ((retval = nc_put_vara_float(grpid2, t2m_dew_MERRA2_varid, start2, count2, fbuf))) ERR(retval);

     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx1= n*ns_S1 + jsat;
         fbuf[idx1]= tqv_buf[idx];
       }
       n++;
     }
     if ((retval = nc_put_vara_float(grpid2, tqv_MERRA2_varid, start2, count2, fbuf))) ERR(retval);

     n= 0;
     for (isat= isat_start; isat <= isat_end; isat++) {
       for (jsat= 0; jsat < ns_S1; jsat++) {
         idx= isat*ns_S1 + jsat;
         idx1= n*ns_S1 + jsat;
         sbuf[idx1]= h273_buf[idx];
       }
       n++;
     }
     if ((retval = nc_put_vara_short(grpid2, h273_MERRA2_varid, start2, count2, sbuf))) ERR(retval);

   }


   /*--- Observed radiometer TB and the top-ranked TB database entries for the NS and MS searches */

   count3[0] = nl_out;
   count3[1] = ns_S1;
   count3[2] = nc_out;
   start3[0] = 0;
   start3[1] = 0;
   start3[2] = 0;

   n= 0;
   for (isat= isat_start; isat <= isat_end; isat++) {
     for (jsat= 0; jsat < ns_S1; jsat++) {
       for (k= 0; k < nc_out; k++) {
         idx= isat*ns_S1*nc_out + jsat*nc_out + k;
         idx1= n*ns_S1*nc_out + jsat*nc_out + k;
         fbuf3[idx1]= Tc_buf[idx];
       }
     }
     n++;
   }
   if ((retval = nc_put_vara_float(ncid, tb_varid, start3, count3, fbuf3))) ERR(retval);

   if ( save_top_ranked > 0 ) {
     if ((retval = nc_put_vara_float(grpid_NS, tb1_NS_varid, start3, count3, tb_1_NS_buf))) ERR(retval);
     if ((retval = nc_put_vara_float(grpid_MS, tb1_MS_varid, start3, count3, tb_1_MS_buf))) ERR(retval);
   }

   /*--- emissivity vector and its EPC ---*/

   count3[0] = nl_out;
   count3[1] = ns_S1;
   count3[2] = NEM;
   start3[0] = 0;
   start3[1] = 0;
   start3[2] = 0;

   if ( save_pc > 0 ) {
     if ((retval = nc_put_vara_float(ncid, emis_varid, start3, count3, emis_buf))) ERR(retval);
     if ((retval = nc_put_vara_float(ncid, pc_emis_varid, start3, count3, pc_emis_buf))) ERR(retval);
   }

   /*--- top-ranked Ku profile database entry, for the NS search */
   /*--- top-ranked Ku and Ka profile database entries, for the MS search */

   count3[0] = nl_out;
   count3[1] = ns_S1;
   count3[2] = NLEV_DPR;
   start3[0] = 0;
   start3[1] = 0;
   start3[2] = 0;

   if ( save_top_ranked > 0 ) {
     if ((retval = nc_put_vara_short(grpid_NS, z_ku_NS_varid, start3, count3, z_ku_1_NS_buf))) ERR(retval);
     if ((retval = nc_put_vara_short(grpid_MS, z_ku_MS_varid, start3, count3, z_ku_1_MS_buf))) ERR(retval);
     if ((retval = nc_put_vara_short(grpid_MS, z_ka_MS_varid, start3, count3, z_ka_1_MS_buf))) ERR(retval);
   }

   /*--- Observed DPR Ku and Ka profiles */

   if ( save_dpr_profiles > 0 && have_dpr > 0 ) {
     if ((retval = nc_put_vara_short(grpid3, z_ku_DPR_varid, start3, count3, z_ku_DPR_buf))) ERR(retval);
     if ((retval = nc_put_vara_short(grpid3, z_ka_DPR_varid, start3, count3, z_ka_DPR_buf))) ERR(retval);
   }

   /*--- precip profiles based on CMB and DPR, for the NS and MS search */

   count3[0] = nl_out;
   count3[1] = ns_S1;
   count3[2] = NLEV_PRECIP;
   start3[0] = 0;
   start3[1] = 0;
   start3[2] = 0;

   if ((retval = nc_put_vara_short(grpid_NS, precip1_prof_NS_varid, start3, count3, precip1_prof_NS_buf))) ERR(retval);
   if ((retval = nc_put_vara_short(grpid_NS, precip2_prof_NS_varid, start3, count3, precip2_prof_NS_buf))) ERR(retval);
   if ((retval = nc_put_vara_short(grpid_MS, precip1_prof_MS_varid, start3, count3, precip1_prof_MS_buf))) ERR(retval);
   if ((retval = nc_put_vara_short(grpid_MS, precip2_prof_MS_varid, start3, count3, precip2_prof_MS_buf))) ERR(retval);

   /*--- Actual CMB precip profile */

   if ( have_cmb > 0 ) {
    if ((retval = nc_put_vara_short(grpid4, precip_prof_NS_CMB_varid, start3, count3, precip_prof_NS_CMB_buf))) ERR(retval);
    if ((retval = nc_put_vara_short(grpid4, precip_prof_MS_CMB_varid, start3, count3, precip_prof_MS_CMB_buf))) ERR(retval);
   }


   if ((retval = nc_close(ncid))) ERR(retval);
   printf("Closed %s\n", outfile);

   /*-------------------------------------------------------------------------*/

   free(yyyy_S1_buf);
   free(mm_S1_buf);
   free(dd_S1_buf);
   free(hh_S1_buf);
   free(mn_S1_buf);
   free(ss_S1_buf);

   free(lat_S1_buf);
   free(lon_S1_buf);

   free(sclat_S1_buf);
   free(sclon_S1_buf);

   free(qual_S1_buf);
   free(qual_S2_buf);

   free(Tc_buf);
   free(pc_emis_buf);
   free(emis_buf);

   free(z_ku_1_NS_buf);
   free(z_ku_1_MS_buf);
   free(z_ka_1_MS_buf);
   free(tb_1_NS_buf);
   free(tb_1_MS_buf);

   if ( have_model > 0 ) {
     free(h_buf);
     free(t_buf);
     free(qv_buf);
   }

   if ( have_dpr > 0 ) {
     free(lat_NS_buf);
     free(lon_NS_buf);
     free(zFactorMeasured_NS_buf);
     free(zFactorMeasured_MS_buf);
   }

   if ( have_cmb > 0 ) {
     free(lat_NS_cmb_buf);
     free(lon_NS_cmb_buf);
   }

   free(fbuf3);
   free(fbuf4);
   free(fbuf);
   free(sbuf);
   free(bbuf);
   free(dbuf);
   free(ibuf);


   exit(0);
}


