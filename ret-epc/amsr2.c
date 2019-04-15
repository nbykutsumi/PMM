/* tc array for AMSR2 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define RADEARTH 6371
#define RTD 57.29578
#define DTR 0.017453

void tc_amsr2(float *Tc_S1_buf, float *Tc_S2_buf, float *Tc_S3_buf, float *Tc_S4_buf, float *Tc_S5_buf, float *Tc_S6_buf,
          int nl_S1, int ns_S1, int nc_S1, 
          int nl_S2, int ns_S2, int nc_S2, 
          int nl_S3, int ns_S3, int nc_S3, 
          int nl_S4, int ns_S4, int nc_S4, 
          int nl_S5, int ns_S5, int nc_S5, 
          int nl_S6, int ns_S6, int nc_S6, 
          float *lat_S1_buf, float *lon_S1_buf, float *lat_S5_buf, float *lon_S5_buf, float *lat_S6_buf, float *lon_S6_buf,
          char *qual_S1_buf, char *qual_S2_buf, char *qual_S3_buf, char *qual_S4_buf, char *qual_S5_buf, char *qual_S6_buf,
          float *incidenceAngle_S1_buf, 
          float *tc)
{

  int nc_out= nc_S1 + nc_S2 + nc_S3 + nc_S4 + nc_S5;
  int i,j,k,idx,idx1, idx2, qual;
  float tb, tb2[nc_out];
  float lat1, lon1, lat2, lon2, km, km_max= -1.0, inc;

  for (i=0; i< nl_S1*ns_S1*nc_out; i++) tc[i]= -999.0;
  for (i=0; i< nc_out; i++) tb2[i]= -999.0;

  for (i=0; i< nl_S1; i++) {
    for (j=0; j<ns_S1; j++) {
      idx= i*ns_S1+j;
      idx2= i*ns_S5+(2*j+1);
      lat1= lat_S1_buf[idx];
      lon1= lon_S1_buf[idx];
      inc= incidenceAngle_S1_buf[idx];
      if ( fabs(lat1) > 90 || fabs(lon1) > 180 ) continue;
      lat2= lat_S5_buf[idx2];
      lon2= lon_S5_buf[idx2];
      if ( fabs(lat2) > 90 || fabs(lon2) > 180 ) continue;

      km= RADEARTH*acos(cos(DTR*lon1-DTR*lon2)*cos(DTR*lat1)*cos(DTR*lat2) + sin(DTR*lat1)*sin(DTR*lat2));
      if ( km > km_max ) km_max= km;
      printf("%4d %3d %8.3f %8.3f %8.3f %8.3f %8.3f km=%8.3f  ", i,j, inc, lat1, lon1, lat2, lon2, km);


      qual= qual_S1_buf[idx];
      printf("%4d ", qual);
      if ( qual >= 0 ) {
        for (k= 0; k< nc_S1; k++) {
          idx1= i*ns_S1*nc_S1 + j*nc_S1 + k;
          tb= Tc_S1_buf[idx1];
          if ( tb > 20 && tb < 400 ) {
            if ( qual > 0 ) tb*= -1.0;
            tb2[k]= tb;
          }
        }
      }

      qual= qual_S2_buf[idx];
      printf("%4d ", qual);
      if ( qual >= 0 ) {
        for (k= 0; k< nc_S2; k++) {
          idx1= i*ns_S2*nc_S2 + j*nc_S2 + k;
          tb= Tc_S2_buf[idx1];
          if ( tb > 20 && tb < 400 ) {
            if ( qual > 0 ) tb*= -1.0;
            tb2[k+nc_S1]= tb;
          }
        }
      }

      qual= qual_S3_buf[idx];
      printf("%4d ", qual);
      if ( qual >= 0 ) {
        for (k= 0; k< nc_S3; k++) {
          idx1= i*ns_S3*nc_S3 + j*nc_S3 + k;
          tb= Tc_S3_buf[idx1];
          if ( tb > 20 && tb < 400 ) {
            if ( qual > 0 ) tb*= -1.0;
            tb2[k+nc_S1+nc_S2]= tb;
          }
        }
      }

      qual= qual_S4_buf[idx];
      printf("%4d ", qual);
      if ( qual >= 0 ) {
        for (k= 0; k< nc_S4; k++) {
          idx1= i*ns_S4*nc_S4 + j*nc_S4 + k;
          tb= Tc_S4_buf[idx1];
          if ( tb > 20 && tb < 400 ) {
            if ( qual > 0 ) tb*= -1.0;
            tb2[k+nc_S1+nc_S2+nc_S3]= tb;
          }
        }
      }

      qual= qual_S5_buf[idx2];
      printf("%4d ", qual);
      if ( qual >= 0 ) {
        for (k= 0; k< nc_S5; k++) {
          idx1= i*ns_S5*nc_S5 + (2*j+1)*nc_S5 + k;
          tb= Tc_S5_buf[idx1];
          if ( tb > 20 && tb < 400 ) {
            if ( qual > 0 ) tb*= -1.0;
            tb2[k+nc_S1+nc_S2+nc_S3+nc_S4]= tb;
          }
        }
      }

      for (k= 0; k< nc_out; k++) {
        idx1= i*ns_S1*nc_out + j*nc_out + k;
        tc[idx1]= tb2[k];
        printf("%7.2f ", tb2[k]);
      }
      printf("\n");

    }
  }

}

